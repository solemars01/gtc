!######################################################################
! Module for reading plasma configuration file prepared by M3D-C1 code, 
! specifically with magnetic islands created by RMP. 
! Written by Lei Shi, Oct 11, 2016

! Subroutines: read_m3dc1_dims, read_m3dc1_profs, check_m3dc1, check_nc
!######################################################################

module m3dc1

use precision, only: lk
use global_parameters, only: r0, b0, pi, fielddir, gtcout
use equilibrium, only: lst, lsp, psiw, ped, spdpsi, spdtheta, qpsi, gpsi, &
     rpsi, bsp, xsp, zsp, gsp, jsp, rd, cpsi
use spline_function, only: construct_spline1d, construct_spline2d
use smooth_function, only: smooth_periodic
 
implicit none
private
public read_m3dc1_dims, read_m3dc1_eq, check_Jacobian

! The physical unit conversion coefficients from m3dc1 output file
! quantities to CGS quantities
  real(lk), parameter :: m3dc1_length_unit=100.0_lk, m3dc1_time_unit=1.0_lk, &
       m3dc1_magnetic_unit=10000.0_lk

contains

!#####################################################################

! Subroutine read_m3dc1_dims
! Purpose: read in the dimensions of psi and theta grids
! Input: m3dfile: A string contains the path to the M3DC1 file
! Output: lsp: number of grid points in psi direction
!         lst: number of grid points in theta direction 
subroutine read_m3dc1_dims(m3dfile)

include 'netcdf.inc'

! netcdf file name should be passing in
character(len=*), intent(in):: m3dfile

! some preset dimension names
character(len=*), parameter :: psi_dimname = 'npsi', theta_dimname = 'mpol' 

! netcdf ids
integer :: ncid, psi_dimid, theta_dimid


! read in lsp and lst from the netcdf file

! open the netcdf file
call check_nf( nf_open(m3dfile, nf_nowrite, ncid))

! get the dimension ids
call check_nf( nf_inq_dimid(ncid, psi_dimname, psi_dimid))
call check_nf( nf_inq_dimid(ncid, theta_dimname, theta_dimid))

! get the dimension lengths
call check_nf( nf_inq_dimlen(ncid, psi_dimid, lsp))
call check_nf( nf_inq_dimlen(ncid, theta_dimid, lst))

! close the netcdf file
call check_nf( nf_close(ncid))

! assign the global variables lsp and lst
! Note that in m3dc1 output files, inner most psi is at the first stepsize 
! delta_psi, and the out most psi is the one right inside the separatrix, psi_ed-delta_psi.
! So, to include the magnetic axis psi=0 and the Last Closed Flux Surface (LCFS) 
! psi=psi_ed, we need to add two value points, one at axis, one at LCFS. 
lsp = lsp+2

! For theta, m3dc1 output file range from 0 to 2pi-delta_theta, to include the periodic 
! data point at 2pi, we add one data point.
lst = lst+1

! write out the read in dimensions
write(gtcout, *), "*** Read in M3DC1 dimensions: lsp = ", lsp, ", lst = ", lst
end subroutine read_m3dc1_dims
!##############################################################################################


! Subroutine read_m3dc1_eq
! Purpose: readin the axisymmetric equilibrium quantities
! Input: m3dfile: A string contains the path to the M3DC1 file
! Fills:
!   Spline constants:
!     spdpsi
!   New radial coordinate spline coefficient
!     rpsi (outer mid-plane R-R0 as function of psi_p)
!   2D spline coefficients:
!     Total B field : bsp, 
!     R, Z coordinates : xsp, zsp, 
!     Jacobian : gsp,
!     Toroidal current: cpsi (1D), rd (essentially 1D, coeffs in dtheta are all 0), duplicated with cpsi
!   1D spline coefficients: 
!     Safety factor : qpsi, 
!     Poloidal current : gpsi, 
!     Minor radius : rpsi,
!     Toroidal current: cpsi 

subroutine read_m3dc1_eq(m3dfile)
include 'netcdf.inc'

! netcdf file name should be passing in
character(len=*), intent(in):: m3dfile

! netcdf variable names
character(len=*), parameter :: q_name='q', g_name='F', I_name='current', &
     psi_name='psi', psin_name='psi_norm', r_name='rpath', z_name='zpath', &
     bp_name='Bp'
! netcdf ids
integer :: ncid, q_id, g_id, I_id, psi_id, psin_id, r_id, z_id, bp_id

! Note that M3DC1 data is stored as single precision arrays
! temporary container for 1D variables
real(4), dimension(:), allocatable :: temp1d
! temporary container for 2D variables
real(4), dimension(:,:), allocatable :: temp2d

! loop index
integer :: i, j
! temporary psi information storage
! we need to use these values to calculate the correct psiw since poloidal flux in M3DC1 may have different definition on axis and on wall
real(lk) :: psi_0, psi_1, psin_0, psin_1

! The vacuum permeability constant in SI
real(lk), parameter:: mu0 = 1.2566370614D-6

! open m3dc1 netcdf file
call check_nf( nf_open(m3dfile, nf_nowrite, ncid))

! get the variable dimensions
call check_nf( nf_inq_varid(ncid, q_name, q_id))
call check_nf( nf_inq_varid(ncid, g_name, g_id))
call check_nf( nf_inq_varid(ncid, I_name, I_id))
call check_nf( nf_inq_varid(ncid, psi_name, psi_id))
call check_nf( nf_inq_varid(ncid, psin_name, psin_id))
call check_nf( nf_inq_varid(ncid, r_name, r_id))
call check_nf( nf_inq_varid(ncid, z_name, z_id))
call check_nf( nf_inq_varid(ncid, bp_name, bp_id))

! read in the variables

! 1D variables
allocate (temp1d(lsp-2))

! Poloidal flux in physical unit (weber * meter**2): psi0, psi1
call check_nf( nf_get_var(ncid, psi_id, temp1d))
! store the first and last values of psi
psi_0 = temp1d(1)
psi_1 = temp1d(lsp-2)

! Normalized poloidal flux in range [0, 1]: psin0, psin1
call check_nf( nf_get_var(ncid, psin_id, temp1d))
! store the first and last values of psi
psin_0 = temp1d(1)
psin_1 = temp1d(lsp-2)

! The last closed flux surface ped
! Calculate ped in terms of psi_0, psi_1, psin_0, psin_1
! (psi_0-psi_1)/ped = (psin_0 - psin_1)
ped = (psi_0-psi_1)/(psin_0-psin_1)
! convert unit
ped = ped * m3dc1_length_unit**2 * m3dc1_magnetic_unit
! set psiw to be the same as ped
psiw = ped

! Safety factor: qpsi
call check_nf( nf_get_var(ncid, q_id, temp1d))

! fill in the values of spline coefficients
qpsi(1, 2:lsp-1)=temp1d(:)

! on-axis and at psiw values need to be extrapolated
! both sides use linear extrapolation with derivative calculated by '431' formula 
call linear_extrapolation(qpsi(1,:), 1, -1)
call linear_extrapolation(qpsi(1,:), lsp, 1)

! Covariant toroidal magnetic field: gpsi
call check_nf( nf_get_var(ncid, g_id, temp1d))

! fill in the values of spline coefficients
! on-axis value set to be the same as the first point
gpsi(1, 2:lsp-1)=temp1d(:)
! extrapolated on both sides
call linear_extrapolation(gpsi(1,:), 1, -1)
call linear_extrapolation(gpsi(1,:), lsp, 1)
! convert unit
gpsi(1, :) = gpsi(1, :) * m3dc1_magnetic_unit * m3dc1_length_unit

! Covariant poloidal magnetic field: cpsi
call check_nf( nf_get_var(ncid, I_id, temp1d))
! fill in the values of spline coefficients
! on-axis I value should be zero since no toroidal current enclosed in the "zero width" flux surface
cpsi(1, 1) = 0
cpsi(1, 2:lsp-1)=temp1d(:)
! extrapolate on outside
call linear_extrapolation(cpsi(1,:), lsp, 1)

! The read in quantity is toroidal current, times 2pi and magnetic permeability mu0 to obtain the
! correct poloidal magnetic field quantity
cpsi(1, :) = cpsi(1, :)*mu0/(2*pi)
! convert the unit
cpsi(1, :) = cpsi(1, :)* m3dc1_length_unit * m3dc1_magnetic_unit  

! 2D variables
! Note that in M3DC1 output file, 2D quantities are stored so that the fastest changing index is in theta
! So, in Fortran, the first dimension should be theta, the second dimension is psi
allocate (temp2d(lst-1, lsp-2))

! Major Radius : R
call check_nf( nf_get_var(ncid, r_id, temp2d))
! Be careful that in GTC spline coefficients, psi is in front of theta

do i=2, lsp-1
   do j=1, lst-1
      xsp(1, i, j) = temp2d(j, i-1)
   enddo
   ! periodic in theta
   xsp(1, i, lst) = temp2d(1, i-1)
enddo
! on-axis value take the average of the inner most flux surface
xsp(1, 1, 1)=0
do j=1, lst-1
   xsp(1, 1, 1) = xsp(1,1,1)+ temp2d(j, 1)
enddo
xsp(1, 1, 1) = xsp(1, 1, 1)/(lst-1)
xsp(1, 1, :) = xsp(1, 1, 1)

! Extrapolate to obtain the values at LCFS.
! Estimate the first order derivative using the "431" formula. achieve the second order accuracy.
do j=1, lst
   call linear_extrapolation(xsp(1,:,j), lsp, 1)
enddo
! convert from M3DC1 unit to CGS unit
xsp(1,:,:) = xsp(1,:,:) * m3dc1_length_unit
! reset the magnetic axis R0 value
r0 = xsp(1, 1, 1)
! reset the magnetic axis B0 value
b0 = abs(gpsi(1, 1)/r0)

! Vertical Coordinate: Z
call check_nf( nf_get_var(ncid, z_id, temp2d))
! Be careful that in GTC spline coefficients, psi is in front of theta
do i=2, lsp-1
   do j=1, lst-1
      zsp(1, i, j) = temp2d(j, i-1)
   enddo
   ! periodic in theta
   zsp(1, i, lst) = temp2d(1, i-1)
enddo
! Since theta=0 line is defined as horizontal, on-axis Z value should equal Z values at theta=0 
zsp(1, 1, :)=zsp(1, 2, 1)

! Z value at separatrix needs to be extrapolated from inside. 
! Estimate the first order derivative using the "431" formula. achieve the second order accuracy.
do j=1, lst
   call linear_extrapolation(zsp(1,:,j), lsp, 1)
enddo
! convert the unit
   zsp(1,:,:) = zsp(1,:,:) * m3dc1_length_unit

! Poloidal magnetic field strength: Bp
! We'll temporarily store this quantity in bsp
call check_nf( nf_get_var(ncid, bp_id, temp2d))
! Be careful that in GTC spline coefficients, psi is in front of theta
do i=2, lsp-1
   do j=1, lst-1
      bsp(1, i, j) = temp2d(j, i-1)
   enddo
   ! periodic in theta
   bsp(1, i, lst) = temp2d(1, i-1)
enddo
! on-axis Bp should be zero
bsp(1, 1, :)=0

! Bp value at psiw needs to be extrapolated from inside. 
! Estimate the first order derivative using the "431" formula. achieve the second order accuracy.
do j=1, lst
   call linear_extrapolation(bsp(1,:,j), lsp, 1)
enddo

! convert the unit
bsp(1, :, :) = bsp(1, :, :)* m3dc1_magnetic_unit

! Calculate the total B based on B_p, gpsi and xsp
do i=1, lsp
   do j=1, lst
      bsp(1,i,j) = sqrt(bsp(1,i,j)**2 + (gpsi(1,i)/xsp(1,i,j))**2)
   enddo
enddo

! close netcdf file
call check_nf (nf_close(ncid))

! convert the coordinates
call coordinate_conversion

! Normalize the quantities into GTC unit
! GTC units:
! Length : R0
! Magnetic field : B0
! Time : 1/omega_cp
call normalize_profs 

! Smooth the periodic direction for noise control
do i=1, lsp
   call smooth_periodic(xsp(1,i,:), lst)
   call smooth_periodic(zsp(1,i,:), lst)
   call smooth_periodic(bsp(1,i,:), lst)
enddo
! generate the spline coefficients
call create_splines

! populate rd array using the calculated cpsi
! cpsi is actually used in GTC calculation, rd is just for equilibrium output
rd(1:3,:,1) = cpsi(:,:)
do i=1,3
   do j=1,lsp
      rd(i,j,:) = rd(i,j,1)
   enddo
enddo
rd(4:9,:,:)=0

! check the self-consistency of the equilibrium, write out warning when error is too large
! Jacobian splines will be created here.
! Jacobians are checked at the end of eqdata subroutine. No need to check here.
! call check_m3dc1

! write out a message
write(gtcout, *) '***M3DC1 equilibrium profiles successfully read.'
write(gtcout, *) 'q_axis = ', qpsi(1,2),'q_LCFS = ', qpsi(1,lsp-1)
write(gtcout, *) 'R0=', r0, 'R_min/R0=', xsp(1, lsp, lst/2), 'R_max/R0=', xsp(1, lsp, 1)
write(gtcout, *) 'Z0/R0=', zsp(1,1,1), 'Z_min(max)/R0=', zsp(1, lsp, lst/4), 'Z_max(min)/R0=', zsp(1, lsp, lst*3/4)
write(gtcout, *) 'B0=', b0, 'B_out/B0=', bsp(1, lsp, 1), 'B_in/B0=', bsp(1, lsp, lst/2)
if (fielddir == 0) then
   write(gtcout, *) 'fielddir=', fielddir, 'I:out, B:out'
else if (fielddir ==1) then
   write(gtcout, *) 'fielddir=', fielddir, 'I:out, B:in'
else if (fielddir ==2) then
   write(gtcout, *) 'fielddir=', fielddir, 'I:in, B:in'
else
   write(gtcout, *) 'fielddir=', fielddir, 'I:in, B:out'
endif
write(gtcout, *) '**********************************************'
end subroutine read_m3dc1_eq
!####################################################################################

! Subroutine normalize_profs
! Purpose: Normalize the equilibrium quantities into GTC units
! GTC units:
! Length : R0 (major radius of magnetic axis)
! Magnetic field : B0 (on-axis magnetic field)
! Time : 1/omega_cp (invert cyclotron frequency of proton)
subroutine normalize_profs
! some useful derived units
real(lk) :: flux_unit

! normalize lengths
xsp(:,:,:) = xsp(:,:,:)/r0
zsp(:,:,:) = zsp(:,:,:)/r0

! normalize magnetic field
bsp(:,:,:) = bsp(:,:,:)/b0

! normalize covariant field quantities
cpsi(1,:) = cpsi(1,:)/(r0*b0) 
gpsi(1,:) = gpsi(1,:)/(r0*b0)

! normalize flux 
flux_unit = r0*r0*b0
psiw = psiw/flux_unit
ped = ped/flux_unit

end subroutine normalize_profs
!#####################################################################################

! Subroutine coordinate_conversion
! Purpose: convert all the quantities in M3DC1 coordinates into GTC coordinates
! See the documentation for detailed convention information

subroutine coordinate_conversion

integer:: I_sign, B_sign, j
real(lk):: temp(lsp)

! The signs of two quantities are required to determine the Ip and B_T directions
! For simplicity, we just use the signs of toroidal current and the g(psi)
if (cpsi(1, lsp)>0) then
   I_sign = 1
else 
   I_sign = -1
endif

if (gpsi(1, 1)>0) then
   B_sign = 1
else
   B_sign = -1
endif
! Check if sign of qpsi is compatible with I_sign and B_sign
! sign of q should always equal I_sign*B_sign
if (I_sign+B_sign == 0 .and. qpsi(1,1)>0) then
   write(gtcout, *) 'M3DC1 EQUILIBRIUM ERROR: Signs of toroidal current and toroidal B &
        field don''t agree with sign of q. Check the equilibrium file.'
   write(gtcout, *) 'cpsi(1,lsp) = ', cpsi(1,lsp), ', I_sign = ', I_sign, ', gpsi(1,1)= ', gpsi(1,1), &
        ', B_sign = ', B_sign
   stop
elseif (I_sign+B_sign /= 0 .and. qpsi(1,1)<0) then
   write(gtcout, *) 'M3DC1 EQUILIBRIUM ERROR: Signs of toroidal current and toroidal B &
             field don''t agree with sign of q. Check the equilibrium file.'
    write(gtcout, *) 'cpsi(1,lsp) = ', cpsi(1,lsp), ', I_sign = ', I_sign, ', gpsi(1,1)= ', gpsi(1,1), &
        ', B_sign = ', B_sign
   stop
endif

gpsi(1,:) = I_sign * gpsi(1,:)
cpsi(1,:) = abs(cpsi(1,:))
psiw = abs(psiw)
ped = abs(ped)

! if I_sign is negative, then M3DC1 theta and GTC theta are opposite
! All 2D quantities need to be flipped in theta, 
! i.e. f_gtc(theta_gtc) = f_m3d(2*pi - theta_gtc)
! Need a buffer for switching the two
if (I_sign<0) then
   do j=1,lst/2  ! Switch the first half with the second half 
      temp(:) = xsp(1,:,j)
      xsp(1,:,j) = xsp(1,:,lst-j+1)
      xsp(1,:,lst-j+1) = temp(:)

      temp(:) = zsp(1,:,j)
      zsp(1,:,j) = zsp(1,:,lst-j+1)
      zsp(1,:,lst-j+1) = temp(:)

      temp(:) = bsp(1,:,j)
      bsp(1,:,j) = bsp(1,:,lst-j+1)
      bsp(1,:,lst-j+1) = temp(:)
   enddo
endif

! set fielddir value to corresponding cases, for constructing field line
! aligned GTC mesh
! 4 cases: fielddir =
! 0 : I out, B out
! 1 : I out, B in
! 2 : I in, B in
! 3 : I in, B out
if (I_sign < 0) then
   if (B_sign < 0) then
      fielddir = 0
   else
      fielddir = 1
   endif
else
   if (B_sign < 0) then
      fielddir = 3
   else
      fielddir = 2
   endif
endif

end subroutine coordinate_conversion
!############################################################################

! Subroutine create_splines
! Purpose: generate the spline coefficients for all the read in quantities

subroutine create_splines
integer:: i,j

! global parameter spdpsi will be used in the future 
spdpsi = psiw/(lsp-1)

! Create the new minor radius coordinate rpsi
! rpsi is simply the radial distance along the outer mid-plane
! note that xsp is already normalized to R0, so rpsi is automatically normalized
do i=1, lsp
   rpsi(1, i) = xsp(1,i,1)-1
enddo

call construct_spline1d(lsp, spdpsi, rpsi,1)
call construct_spline1d(lsp,spdpsi,qpsi,0)
call construct_spline1d(lsp,spdpsi,gpsi,0)
call construct_spline1d(lsp,spdpsi,cpsi,0)


! global paramter spdtheta will be used in the future
spdtheta = 2*pi/(lst-1)

call construct_spline2d(lsp, lst, spdpsi, spdtheta, xsp, 1, 2)
call construct_spline2d(lsp, lst, spdpsi, spdtheta, zsp, 1, 2)
call construct_spline2d(lsp, lst, spdpsi, spdtheta, bsp, 1, 2)

end subroutine create_splines
!##########################################################################################

! Subroutine check_m3dc1
! Purpose: check the consistence of the read in quantities
! Method:
! 1. Check Jacobian
!    1/Jacobian = (Grad psi) X (Grad theta) . (Grad zeta)
!    1/Jacobian = B.B /(gq+I)
!    The first formula can be evaluated from R(psi,theta), Z(psi,theta).
! 2. Check poloidal flux
!    Poloidal flux is a read-in quantity. On the other hand, it can be calculated by integrating 
!    B_pol along major radius on mid-plane 
subroutine check_m3dc1
call check_Jacobian
!call check_poloidal_flux
end subroutine check_m3dc1

subroutine check_Jacobian

integer:: j
real(lk):: error(lsp, lst), maxerr
! calculate Jacobian using J^-1 = (Grad psi) X (Grad theta) . Grad(zeta) formula
! Note that J = R*|(dR_dpsi*dZ_dtheta - dR_dtheta*dZ_dpsi)|

jsp(1,:,:) = abs((xsp(2,:,:)*zsp(4,:,:)-xsp(4,:,:)*zsp(2,:,:))*xsp(1,:,:))
do j=1,lst
   gsp(1,:,j) = (gpsi(1,:)*qpsi(1,:) + cpsi(1,:))/(bsp(1,:,j)**2)
enddo

if (mod(lst,2)==0) then
   call construct_spline2d(lsp, lst, spdpsi, spdtheta, gsp, 1, 2)
   call construct_spline2d(lsp, lst, spdpsi, spdtheta, jsp, 1, 2)
else
   call construct_spline2d(lsp, lst, spdpsi, spdtheta, gsp, 1, 0)
   call construct_spline2d(lsp, lst, spdpsi, spdtheta, jsp, 1, 0)
endif

!error = abs(jacobian_metric-jacobian_boozer)/jacobian_boozer
!maxerr = maxval(error)
!if ( maxval(error) > 1e-4 ) then
!   print *, 'Max relative Jacobian error exceeds 1e-4, maxerr= ', maxerr
!endif

end subroutine check_Jacobian

subroutine check_poloidal_flux
! Checking the self-consistency of the poloidal field and the poloidal flux
! NOT IMPLEMENTED YET
end subroutine check_poloidal_flux
!############################################################

!##########################################################################################
! subroutine to check the status of netcdf calls
subroutine check_nf(status)
  include 'netcdf.inc'
  integer, intent ( in) :: status
  
  if(status /= nf_noerr) then 
     write(gtcout,*) 'NetCDF ERROR:', nf_strerror(status)
     stop "Stopped"
  end if
end subroutine check_nf

!#######################################################################
! Some used short subroutines
!####################################################################

! Linear extrapolation using the '341' formula to estimate the derivative at boundary

! f: the array contains the values that needs to be extrapolated at point n
! n: the index on which the value is needed
! dir: extrapolation direction 
!      dir==1 means forward, using n-1, n-2, n-3 to extrapolate n
!      dir==-1 means backward, using n+1, n+2, n+3 to extrapolate n                        

! This formula is only valid for equally space function
! n must be greater than 3
subroutine linear_extrapolation(f, n, dir)
  real(lk), intent(inout) :: f(:)
  integer, intent(in) :: n, dir
  if (dir==1) then
     f(n) = 2.5_lk*f(n-1) - 2.0_lk*f(n-2) + 0.5_lk*f(n-3)
  else
     f(n) = 2.5_lk*f(n+1) - 2.0_lk*f(n+2) + 0.5_lk*f(n+3)
  endif
end subroutine linear_extrapolation


!#### ALPHA LOADING PART ##########################
! NOT FINISHED
!###################################################

! 3D data creation from m and n harmonics
! y = sum_nm (ycn*cos(n*zeta - m*theta) + ysn*sin(n*zeta-m*theta))

! first sum over m
subroutine poloidal_harmonic_sum(alphamn_re, alphamn_im, alphan_re, alphan_im, ndim, m)
  use precision, only: lk
  use global_parameters, only: pi
  use equilibrium, only: lsp, lst
  implicit none
  integer, intent(in):: ndim, m(lst)
  real(lk), intent(in):: alphamn_re(ndim,lsp,lst), alphamn_im(ndim,lsp,lst)
  real(lk), intent(out):: alphan_re(ndim,lsp,lst), alphan_im(ndim,lsp,lst)
  integer i,j,k,l
  real(lk) dtheta, theta, thetam

  dtheta = 2.0_lk*pi/real(lst-1, lk)
  alphan_re = 0.0_lk
  alphan_im = 0.0_lk
  do i=1, lst
     theta = real(i-1,lk)*dtheta
     do j=1, lst
        thetam = theta*real(m(j),lk)
        do k=1, lsp
           do l=1, ndim
              alphan_re(l, k, i) = alphan_re(l,k,i)+ alphamn_re(l,k,j)*cos(thetam) + &
                   alphamn_im(l,k,j)*sin(thetam)
              alphan_im(l, k, i) = alphan_im(l,k,i)- alphamn_re(l,k,j)*sin(thetam) + &
                   alphamn_im(l,k,j)*cos(thetam)
           enddo
        enddo
     enddo
  enddo
end subroutine poloidal_harmonic_sum

subroutine toroidal_harmonic_sum(alpha, alphan_re, alphan_im, ndim, n)
  use precision, only: lk
  use global_parameters, only: pi, mtoroidal
  use equilibrium, only: lsp, lst
  implicit none
  integer, intent(in):: ndim, n(ndim)
  real(lk), intent(in):: alphan_re(ndim, lsp, lst), alphan_im(ndim, lsp, lst)
  real(lk), intent(out):: alpha(lsp, lst, mtoroidal+1)

  integer i,j,k,l
  real(lk) dzeta, zeta, zetan
  dzeta = 2.0_lk*pi/real(mtoroidal,lk)
  alpha = 0.0_lk
  do i=1, mtoroidal
     zeta=dzeta*real(i-1,lk)
     do j=1, ndim
        zetan = zeta*real(n(j),lk)
        do k=1, lsp
           do l=1, ndim        
              alpha(l,k,i) = alpha(l,k,i) + alphan_re(l,k,i)*cos(zetan) - alphan_im(l,k,i)*sin(zetan)
           enddo
        enddo
     enddo
  enddo
end subroutine

end module m3dc1
