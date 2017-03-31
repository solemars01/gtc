!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODIFIED BY Lei Shi to have unified interface for 1D, 2D, and 3D quadratic spline
! Date: Nov. 19, 2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!###### MODULE SMOOTH_FUNCTION ###########################
! Contains generic functions for smoothing

module smooth_function
use precision, only: lk
implicit none

contains

! subroutine to smooth the periodic data
subroutine smooth_periodic(array, n)
  real(lk), intent(inout) :: array(:)
  integer, intent (in) :: n
  integer :: i
  real(lk) :: temp(n)
  temp(:) = array(:)
  do i=2, n-1
     array(i) = 0.25_lk *(temp(i-1) + 2.0_lk*temp(i) + temp(i+1))
  enddo
  array(1) = 0.25_lk *(temp(n-1) + 2.0_lk *temp(1) + temp(2))
  array(n) = array(1)
end subroutine smooth_periodic

end module smooth_function


!##############################################################################
!####  MODULE SPLINE FUNCTIONS
!##############################################################################

module spline_function
use precision, only: lk
use global_parameters, only: gtcout
use equilibrium, only: lsp, lst, spdpsi, spdtheta, zepp, ropp, erpp, spcos, spsin,&
     qpsi, gpsi, cpsi, ppsi, torpsi, psitor, bsp, xsp, zsp, gsp, jsp , spdtor, psitor,&
     rpsi, rgpsi, spdrg, psirg, spdim, fsp, rd, nu, dl
use smooth_function, only: smooth_periodic
implicit none

public
private:: spline0, dspline0, spline1, dspline1, spline_periodic, dspline_periodic, &
     spline2d_dep, construct_spline2d_dep, construct_spline0, construct_spline1, &
     construct_spline_periodic, construct_spline2d_periodic_dep
contains

! function of plasma profile
! electron density
real(lk) function spden(pdum,denpp)
  real(lk),intent(in):: pdum, denpp(3,lsp)
  spden = spline1d(pdum,0,lsp,spdpsi,denpp,0)
end function spden

! d_densitye/d_psi
real(lk) function den_dp(pdum,denpp)
  real(lk),intent(in):: pdum, denpp(3,lsp)
  den_dp = spline1d(pdum,1,lsp,spdpsi,denpp,0)
end function den_dp

! electron temperature
real(lk) function sptem(pdum,tpp)
  real(lk),intent(in):: pdum, tpp(3,lsp)
  sptem = spline1d(pdum,0,lsp,spdpsi,tpp,0)
end function sptem

! d_etemperature/d_psi
real(lk) function tem_dp(pdum,tpp)
  real(lk),intent(in):: pdum, tpp(3,lsp)
  tem_dp = spline1d(pdum,1,lsp,spdpsi,tpp,0)
end function tem_dp

! Z_eff
real(lk) function spzeff(pdum)
  real(lk), intent(in):: pdum
  spzeff = spline1d(pdum,0,lsp,spdpsi,zepp,0)
end function spzeff

! toroidal rotation
real(lk) function sptrot(pdum)
  real(lk),intent(in):: pdum
  sptrot = spline1d(pdum,0,lsp,spdpsi,ropp,0)
end function sptrot

! E_r
real(lk) function sper(pdum)
  real(lk),intent(in):: pdum
  sper = spline1d(pdum,0,lsp,spdpsi,erpp,0)
end function sper

! spline cos function
real(lk) function splcos(tdum)
  real(lk),intent(in):: tdum
  splcos = spline1d(tdum,0,lst,spdtheta,spcos,2)
end function splcos

! spline sin function
real(lk) function splsin(tdum)
  real(lk),intent(in):: tdum
  splsin=spline1d(tdum,0,lst,spdtheta,spsin,2)
end function splsin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! function from equilibrium data
! safety factor
real(lk) function spq(pdum)
  real(lk),intent(in):: pdum
  spq = spline1d(pdum,0,lsp,spdpsi,qpsi,0)
end function spq

! dq_dpsi
real(lk) function dq_dp(pdum)
  real(lk),intent(in):: pdum
  dq_dp = spline1d(pdum,1,lsp,spdpsi,qpsi,0)
end function dq_dp

! g and dg/dp
real(lk) function spgpsi(pdum)
  real(lk),intent(in):: pdum
  spgpsi = spline1d(pdum,0,lsp,spdpsi,gpsi,0)
end function spgpsi

real(lk) function spgpsi_dp(pdum)
  real(lk),intent(in):: pdum
  spgpsi_dp = spline1d(pdum,1,lsp,spdpsi,gpsi,0)
end function spgpsi_dp

! I and dI/dp
real(lk) function spcpsi(pdum)
  real(lk),intent(in):: pdum
  spcpsi = spline1d(pdum,0,lsp,spdpsi,cpsi,0)
end function spcpsi

real(lk) function spicpsi_dp(pdum)
  real(lk) pdum
  spicpsi_dp = spline1d(pdum,1,lsp,spdpsi,cpsi,0)
end function spicpsi_dp

! pressure
real(lk) function sppressure(pdum)
  real(lk) pdum
  sppressure = spline1d(pdum,0,lsp,spdpsi,ppsi,0)
end function sppressure

! poloidal to toroidal flux
real(lk) function sptorpsi(pdum)
  real(lk) pdum
  sptorpsi = spline1d(pdum,0,lsp,spdpsi,torpsi,0)
end function sptorpsi

! d(toroidal flux)/d(poloidal flux)
real(lk) function dtorpsi(pdum)
  real(lk) pdum
  dtorpsi = spline1d(pdum,1,lsp,spdpsi,torpsi,0)
end function dtorpsi

! toroidal to poloidal flux
real(lk) function sppsitor(pdum)
  real(lk) pdum
  sppsitor = spline1d(pdum,0,lsp,spdtor,psitor,0)
end function sppsitor

! minor radius
real(lk) function sprpsi(pdum)
  real(lk) pdum
  sprpsi=spline1d(pdum,0,lsp,spdpsi,rpsi,1)
end function sprpsi

! dr/dpsi
real(lk) function dr_dp(pdum)
  real(lk) pdum
  dr_dp=spline1d(pdum,1,lsp,spdpsi,rpsi,1)
end function dr_dp

! psi to radial grid
real(lk) function sprgpsi(pdum)
  real(lk) pdum
  sprgpsi=spline1d(pdum,0,lsp,spdpsi,rgpsi,1)
end function sprgpsi

! radial grid to psi
real(lk) function sppsirg(pdum)
  real(lk) pdum
  sppsirg=spline1d(pdum,0,lsp,spdrg,psirg,1)
end function sppsirg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! magnetic field amplitude
real(lk) function spb(pdum,tdum)
  real(lk) pdum,tdum
  spb = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,bsp,1,2)
end function spb

! db_field/dpsi
real(lk) function dbdp(pdum,tdum)
  real(lk) pdum,tdum
  dbdp = spline2d(pdum,tdum,1,lsp,lst,spdpsi,spdtheta,bsp,1,2)
end function dbdp

! db_field/dtheta
real(lk) function dbdt(pdum,tdum)
  real(lk) pdum,tdum
  dbdt = spline2d(pdum,tdum,2,lsp,lst,spdpsi,spdtheta,bsp,1,2)
end function dbdt

! transform from (psi,theta) to (X,Z)
real(lk) function spx(pdum,tdum)
  real(lk) pdum,tdum
  spx = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,xsp,1,2)
end function spx

! transform from (psi,theta) to (X,Z)
real(lk) function spz(pdum,tdum)
  real(lk) pdum,tdum
  spz = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,zsp,1,2)
end function spz

!transformation dR/dpsi
real(lk) function dxdp(pdum,tdum)
  real(lk) pdum,tdum
  dxdp = spline2d(pdum,tdum,1,lsp,lst,spdpsi,spdtheta,xsp,1,2)
end function dxdp

!dz/dpsi
real(lk) function dzdp(pdum,tdum)
  real(lk) pdum,tdum
  dzdp = spline2d(pdum,tdum,1,lsp,lst,spdpsi,spdtheta,zsp,1,2)
end function dzdp

!dx/dtheta
real(lk) function dxdt(pdum,tdum)
  real(lk) pdum,tdum
  dxdt = spline2d(pdum,tdum,2,lsp,lst,spdpsi,spdtheta,xsp,1,2)
end function dxdt

!dz/dtheta
real(lk) function dzdt(pdum,tdum)
  real(lk) pdum,tdum
  dzdt = spline2d(pdum,tdum,2,lsp,lst,spdpsi,spdtheta,zsp,1,2)
end function dzdt
!!!!!!!!!!!!!!!!!!!--3D-case--!!!!!!!!!!!!!!!!!!!!!!!!!!!
!transformation d\Phi/dpsi
real(lk) function dfdp(pdum,tdum)
  real(lk) pdum,tdum
  dfdp=spline2d_dep(1,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,fsp)
end function dfdp

!transformation d\Phi/dtheta
real(lk) function dfdt(pdum,tdum)
  real(lk) pdum,tdum
  dfdt=spline2d_dep(2,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,fsp)
end function dfdt

!transformation d\Phi/dzeta
real(lk) function dfdz(pdum,tdum)
  real(lk) pdum,tdum
  dfdz=1.0_lk + spline2d_dep(3,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,fsp) !!! 1.0 since \Phi is determined as delta\Phi
end function dfdz

!dx/dzeta
real(lk) function dxdz(pdum,tdum)
  real(lk) pdum,tdum
  dxdz=spline2d_dep(3,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,xsp)
end function dxdz

!dz/dzeta
real(lk) function dzdz(pdum,tdum)
  real(lk) pdum,tdum
  dzdz=spline2d_dep(3,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,zsp)
end function dzdz

! Jacobian calculated by Boozer formula
real(lk) function spjac(pdum,tdum)
  real(lk) pdum,tdum
  spjac = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,gsp,1,2)
end function spjac

! Jacobian calculated by metric tensor
real(lk) function spjacm(pdum,tdum)
  real(lk) pdum,tdum
  spjacm = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,jsp,1,2)
end function spjacm

! toroidal current I
real(lk) function currenti(pdum,tdum)
  real(lk) pdum,tdum
  currenti = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,rd,1,2)
end function currenti

real(lk) function didp(pdum,tdum)
  real(lk) pdum,tdum
  didp = spline2d(pdum,tdum,1,lsp,lst,spdpsi,spdtheta,rd,1,2)
end function didp

! difference between magnetic angle zeta and cylindrical angle phi
real(lk) function zeta2phi(pdum,tdum)
  real(lk) pdum,tdum
  zeta2phi = spline2d_dep(0,spdim,pdum,tdum,lsp,lst,spdpsi,spdtheta,nu)
end function zeta2phi

! delta in B-field contravariant representation
real(lk) function delb(pdum,tdum)
  real(lk) pdum,tdum
! WARNING: dl spline boundary condition not determined.
  delb = spline2d(pdum,tdum,0,lsp,lst,spdpsi,spdtheta,dl,1,2)
end function delb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Generic 1D spline evaluations
! Evaluate f(x) at given x 
! Boundary condition needs to be provided by programmer
! BC 0: use spline0
! BC 1: use spline1
! BC 2: use spline_periodic

real(lk) function spline1d(x, deriv, nx, delx, y, bcx)
  integer, intent(in) :: nx, deriv, bcx
  real(lk), intent(in) :: x, delx, y(3,nx)
  integer i
  real(lk) dx
  
  if (bcx == 0) then
     if (deriv==0) spline1d = spline0(x, nx, delx, y)
     if (deriv==1) spline1d = dspline0(x, nx, delx, y)
  else if (bcx == 1) then
     if (deriv==0) spline1d = spline1(x, nx, delx, y)
     if (deriv==1) spline1d = dspline1(x, nx, delx, y)
  else if (bcx == 2) then
     if (deriv==0) spline1d = spline_periodic(x, nx, delx, y)
     if (deriv==1) spline1d = dspline_periodic(x, nx, delx, y)
  else
     write(gtcout,*) 'Spline1d Error: Invalid boundary condition bcx=', bcx
     stop
  endif

end function spline1d


! 1D spline for radial profiles
real(lk) function spline0(pdum,nsp,delx,y)
  integer nsp,i
  real(lk) pdum,y(3,nsp),delx,dpx

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  spline0=y(1,i)+dpx*y(2,i)+dpx*dpx*y(3,i)

end function spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of 1D spline function
real(lk) function dspline0(pdum,nsp,delx,y)
  integer nsp,i
  real(lk) pdum,y(3,nsp),delx,dpx

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  dspline0=y(2,i)+2.0*dpx*y(3,i)

end function dspline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 1D spline with first point being linear function y=sqrt(x)
real(lk) function spline1(pdum,nsp,delx,y)
  integer nsp,i
  real(lk) pdum,y(3,nsp),delx,dpx,dp2

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
! expand y(x) using sprt(x) near x=0
  if(i==1)dpx=sqrt(dpx)
  dp2=dpx*dpx
  spline1=y(1,i)+dpx*y(2,i)+dp2*y(3,i)

end function spline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! derivative of 1D spline with first point being linear function y=sqrt(x)
real(lk) function dspline1(pdum,nsp,delx,y)
  integer nsp, i
  real(lk) pdum,y(3,nsp),delx,dpx

  i=max(1,min(nsp-1,ceiling(pdum/delx)))
  dpx=pdum-delx*real(i-1)
  if(i==1)dpx=sqrt(dpx)

  if(i==1)then
     dspline1=0.5*y(2,i)/dpx+y(3,i) ! y(x)=y1+y2*sprt(x)+y3*x near x=0
  else
     dspline1=y(2,i)+2.0*dpx*y(3,i) ! y(x)=y1+y2*x+y3*x*x otherwise
  endif

end function dspline1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! periodic spline evaluation 
! Same as spline0 but out of range values are wrapped back in to its periodic
! image location.
real(lk) function spline_periodic(x, nx, delx, y)
  integer, intent(in) :: nx
  real(lk), intent(in) :: x, delx, y(3,nx)
  integer i
  real(lk) dx, x_ev
  x_ev = modulo(x, (nx-1)*delx)
  i=max(1,min(nx-1,ceiling(x/delx)))
  dx=x-delx*real(i-1, lk)
  spline_periodic = y(1,i) + dx*y(2,i) + dx*dx*y(3,i)
end function spline_periodic

real(lk) function dspline_periodic(x, nx, delx, y)
  integer, intent(in) :: nx
  real(lk), intent(in) :: x, delx, y(3,nx)
  integer i
  real(lk) dx, x_ev
  x_ev = modulo(x, (nx-1)*delx)
  i=max(1,min(nx-1,ceiling(x/delx)))
  dx=x-delx*real(i-1, lk)
  dspline_periodic = y(2,i) + 2*dx*y(3,i)
end function dspline_periodic

!#############################################################################


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generic 2D spline function
! Purpose: Evaluate f(x,y) or its partial derivatives using the created spline 
!          coefficients for f

! Written by Lei Shi

! Arguments:
!   x,y : the coordinates where the value of f is evaluated
!   deriv: derivative flag. 0 function value, 1 for partial derivative respect
!          to x, 2 for partial derivative respect to y
!   nx, ny: number of spline nodes in x and y. 
!   delx, dely: stepsize for each spline cell
!               By default, spline starts from x=0,y=0, and get up to x=(nx-1)*delx,
!               and y=(ny-1)*dely. Out of range (x,y) are not checked
!   f: 2D spline coefficient array constructed by construct_spline2d subroutine
!   bcx, bcy: Boundary conditions used when f is created. Consistency is not checked. 

! WARNING: Boundary Condition consistency is not checked! It is the programer's 
!          responsibility to make sure bcx,bcy passed in is the same as those
!          used for construct_spline2d.

real(lk) function spline2d(x, y, deriv, nx, ny, delx, dely, f, bcx, bcy)
  integer, intent(in):: deriv, nx, ny, bcx, bcy
  real(lk), intent(in):: x, y, delx, dely
  real(lk), dimension(:,:,:), intent(in):: f
  integer i,j
  real(lk) dx, dx2, dy, dy2
  
! locate the indexes for the spline cell, out of range values are not checked!
  i = max(1, min(nx-1, ceiling(x/delx)))
  j = max(1, min(ny-1, ceiling(y/dely)))
  dx = x - delx*real(i-1, lk)
  dy = y - dely*real(j-1, lk)
! Boundary Condition 1 requires first cell in form f = f1+f2*sqrt(x)+f3*x
  if (i==1 .and. bcx==1) dx=sqrt(dx)
  if (j==1 .and. bcy==1) dy=sqrt(dy)
  dx2 = dx*dx
  dy2 = dy*dy
  
  spline2d=0.0_lk
  if (deriv==0) then
! evaluate f(x,y)
     spline2d = f(1,i,j)+f(2,i,j)*dx + f(3,i,j)*dx2 +&
          ( f(4,i,j)+f(5,i,j)*dx + f(6,i,j)*dx2)*dy +&
          (f(7,i,j) + f(8,i,j)*dx + f(9,i,j)*dx2)*dy2
  else if (deriv == 1) then
! evaluate df/dx
     if (i==1 .and. bcx==1) then
        spline2d = (f(2,i,j) + f(5,i,j)*dy + f(8,i,j)*dy2)/dx +&
             (f(3,i,j) + f(6,i,j)*dy + f(9,i,j)*dy2)
     else
        spline2d = (f(2,i,j) + f(5,i,j)*dy + f(8,i,j)*dy2) +&
             2*(f(3,i,j) + f(6,i,j)*dy + f(9,i,j)*dy2)*dx
     endif
  else if (deriv == 2) then
     if (j==1 .and. bcy==1) then
        spline2d = (f(4,i,j) + f(5,i,j)*dx + f(6,i,j)*dx2)/dy +&
             (f(7,i,j) + f(8,i,j)*dy + f(9,i,j))
     else
        spline2d = (f(4,i,j) + f(5,i,j)*dx + f(6,i,j)*dx2) +&
             2*(f(7,i,j) + f(8,i,j)*dx + f(9,i,j)*dx2)*dy
     endif
  else
     write(gtcout, *) "Spline2d Error: wrong derivative flag deriv=",deriv
     stop
  endif
end function spline2d

! Deprecated. should be removed when all changes have been made to other parts of the code.  
real(lk) function spline2d_dep(iflag,sd,pdum,tdum,nsp,nst,delp,delt,f)
  integer iflag,nsp,nst,i,j,is,sd
  real(lk) pdum,tdum,f(sd,nsp,nst),dpx,dp2,dtx,dt2,delp,delt,dx(sd),tempres

  i=max(1,min(nsp-1,ceiling(pdum/delp)))
  dpx=pdum-delp*real(i-1)
  ! Boundary Condition 1 corresponds to construct_spline1, first cell uses sqrt(x)
  if(i==1)dpx=sqrt(dpx)
  dp2=dpx*dpx

  j=max(1,min(nst-1,ceiling(tdum/delt)))
  dtx=tdum-delt*real(j-1)
  dt2=dtx*dtx

  dx=0.0_lk
  dx(1)=1.0_lk
  dx(2)=dpx
  dx(3)=dp2
  dx(4:6)=dx(1:3)*dtx
  dx(7:9)=dx(1:3)*dt2

  if(iflag==0)then !2D spline value
     spline2d_dep=f(1,i,j)    +f(2,i,j)*dpx    +f(3,i,j)*dp2 &
          +f(4,i,j)*dtx+f(5,i,j)*dtx*dpx+f(6,i,j)*dtx*dp2 &
          +f(7,i,j)*dt2+f(8,i,j)*dt2*dpx+f(9,i,j)*dt2*dp2
  elseif(iflag==1)then !derivative with respect to x
     if(i==1)then
        spline2d_dep=0.5_lk*(f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2)/dpx+f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2
     else
        spline2d_dep=f(2,i,j)+f(5,i,j)*dtx+f(8,i,j)*dt2+2.0_lk*dpx*(f(3,i,j)+f(6,i,j)*dtx+f(9,i,j)*dt2)
     endif
  elseif(iflag==2) then !derivative with respect to y
     spline2d_dep=f(4,i,j)+f(5,i,j)*dpx+f(6,i,j)*dp2+2.0_lk*dtx*(f(7,i,j)+f(8,i,j)*dpx+f(9,i,j)*dp2)
  elseif(iflag==3) then !derivative with respect to zeta
     if(sd==27)then
       dx(10:18)=dx(1:9)
       tempres=0
       do is = 10, 18
          tempres=tempres+f(is,i,j)*dx(is)
       enddo
       spline2d_dep = tempres
     else
       spline2d_dep=0.0_lk
     endif
  endif
end function spline2d_dep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generic 3D spline function
! Purpose: Evaluate f(x,y,z) or its partial derivatives using the created spline 
!          coefficients for f

! Written by Lei Shi

! Arguments:
!   x,y,z : the coordinates where the value of f is evaluated
!   deriv: derivative flag. 0 function value, 1 for partial derivative respect
!          to x, 2 for partial derivative respect to y, 3 for partial derivative 
!          respect to z
!   nx, ny, nz: number of spline nodes in x, y and z. 
!   delx, dely, delz: stepsize for each spline cell
!               By default, spline starts from x=0,y=0,z=0 and get up to x=(nx-1)*delx,
!               and y=(ny-1)*dely, z=(nz-1)*delz. Out of range (x,y,z) are not checked
!   f: 3D spline coefficient array constructed by construct_spline3d subroutine
!   bcx, bcy, bcz: Boundary conditions used when f is created. 
!                  Consistency is not checked. 

! WARNING: Boundary Condition consistency is not checked! It is the programer's 
!          responsibility to make sure bcx,bcy,bcz passed in is the same as those
!          used for construct_spline3d.

real(lk) function spline3d(x, y, z, deriv, nx, ny, nz, delx, dely, delz,  f, bcx, bcy, bcz)
  integer, intent(in):: deriv, nx, ny, nz, bcx, bcy, bcz
  real(lk), intent(in):: x, y, z, delx, dely, delz, f(27, nx, ny, nz)
  integer i, j, k
  real(lk) dx, dxinv, dy, dyinv, dz, dzinv, dvec(27)
  
! locate the indexes for the spline cell, out of range values are not checked!
  i = max(1, min(nx-1, ceiling(x/delx)))
  j = max(1, min(ny-1, ceiling(y/dely)))
  k = max(1, min(nz-1, ceiling(z/delz)))
  dx = x - delx*real(i-1, lk)
  dy = y - dely*real(j-1, lk)
  dz = z - delz*real(k-1, lk)
  dxinv = 1.0_lk/dx
  dyinv = 1.0_lk/dy
  dzinv = 1.0_lk/dz

! Boundary Condition 1 requires first cell in form f = f1+f2*sqrt(x)+f3*x
  if (i==1 .and. bcx==1) dx=sqrt(dx)
  if (j==1 .and. bcy==1) dy=sqrt(dy)
  if (k==1 .and. bcz==1) dz=sqrt(dz)
  
! Construct the spline vector
! take care of the boundary condition
  if (deriv==0) then
     dvec(1) = 1
     dvec(2) = dx
     dvec(3) = dx*dx
     dvec(4:6) = dvec(1:3)*dy
     dvec(7:9) = dvec(4:6)*dy
     dvec(10:18) = dvec(1:9)*dz
     dvec(19:27) = dvec(10:18)*dz
  else if (deriv==1) then
! partial derivative in x
     if (bcx==1) then
        dvec(1) = 0
        dvec(2) = 0.5_lk*dxinv
        dvec(3) = 1
        dvec(4:6) = dvec(1:3)*dy
        dvec(7:9) = dvec(4:6)*dy
        dvec(10:18) = dvec(1:9)*dz
        dvec(19:27) = dvec(10:18)*dz      
     else
        dvec(1) = 0
        dvec(2) = dx
        dvec(3) = 2*dx
        dvec(4:6) = dvec(1:3)*dy
        dvec(7:9) = dvec(4:6)*dy
        dvec(10:18) = dvec(1:9)*dz
        dvec(19:27) = dvec(10:18)*dz        
     endif
  else if (deriv == 2) then
! partial derivative in y
     if (bcy==1) then
        dvec(1:3) = 0
        dvec(7) = 1
        dvec(8) = dx
        dvec(9) = dx*dx
        dvec(4:6) = dvec(7:9)*0.5_lk*dyinv
        dvec(10:18) = dvec(1:9)*dz
        dvec(19:27) = dvec(10:18)*dz        
     else
        dvec(1:3) = 0
        dvec(4) = 1
        dvec(5) = dx
        dvec(6) = dx*dx
        dvec(7:9) = dvec(4:6)*2*dy
        dvec(10:18) = dvec(1:9)*dz
        dvec(19:27) = dvec(10:18)*dz  
     endif    
  else if (deriv == 3) then
     if (bcz==1) then
        dvec(1:9)=0
        dvec(19) = 1
        dvec(20) = dx
        dvec(21) = dx*dx
        dvec(22:24) = dvec(19:21)*dy
        dvec(25:27) = dvec(22:24)*dy
        dvec(10:18) = dvec(19:27)*dzinv*0.5_lk
     else
        dvec(1:9)=0
        dvec(10) = 1
        dvec(11) = dx
        dvec(12) = dx*dx
        dvec(13:15) = dvec(10:12)*dy
        dvec(16:18) = dvec(13:15)*dy
        dvec(19:27) = dvec(10:18)*dz*2
     endif
  else
     write(gtcout, *) "Spline2d Error: wrong derivative flag deriv=",deriv
     stop
  endif

  spline3d = sum(f(:,i,j,k)*dvec(:))

end function spline3d





!##############################################################################
! Construct spline subroutines
! Generic subroutines are provided for 1D, 2D, and 3D quadratic splines
!    construct_spline1d
!    construct_spline2d
!    construct_spline3d
! Three kinds of boundary conditions are available for each dimension
!    BC 0: Use 341 formula to obtain the estimated derivative at x=0
!    BC 1: Use 341 formula to obtain the estimated derivative at x=xmax, and use
!          f(x) = a + b*sqrt(x) + c*x formula for the first cell at x=0.
!    BC 2: Periodic condition, df/dx(0)=df/dx(xmax)
!###############################################################################

! Generic construct 1D spline 
subroutine construct_spline1d(nx, delx, f, bcx)
  integer, intent(in):: nx, bcx
  real(lk), intent(in):: delx
  real(lk), intent(inout):: f(3,nx)

  if (bcx == 0) then
     call construct_spline0(0,nx,delx,f)
  else if (bcx == 1) then
     call construct_spline1(nx,delx,f)
  else if (bcx == 2) then
     ! smooth the data before construct the spline
     ! periodic spline is very sensitive to the data smoothness
     call smooth_periodic(f(1,:), nx)
     call construct_spline_periodic(nx, delx, f)
  else
     print *,'Wrong spline boundary condition choice:', bcx
     stop 
  endif
end subroutine construct_spline1d

! Generic construct 2D spline
subroutine construct_spline2d(nx,ny,delx,dely,f, bcx, bcy)
  integer, intent(in):: nx, ny, bcx, bcy
  real(lk), intent(in):: delx, dely
  real(lk), intent(inout):: f(9,nx,ny)
  integer i,j,s
  real(lk) ddum1(3,nx),ddum2(3,ny)

  ddum1=0.0
  ddum2=0.0
  
! periodic condition needs to be enforced if bc is set to 2
  if (bcx==2) then
     f(1,1,:) = 0.5_lk*(f(1,1,:)+f(1,nx,:))
     f(1,nx,:) = f(1,1,:)
  endif
  if (bcy==2) then
     do j = 1, ny-1
        ddum1(1,:)=f(1,:,j)
        call construct_spline1d(nx, delx, ddum1, bcx)
        f(1,:,j)=ddum1(1,:) ! the function value may be altered due to smoothing
        f(2,:,j)=ddum1(2,:)
        f(3,:,j)=ddum1(3,:)
     enddo
     f(1,:,ny) = f(1,:,1)
     f(2,:,ny) = f(2,:,1)
     f(3,:,ny) = f(3,:,1)
  else
     do j = 1, ny
        ddum1(1,:)=f(1,:,j)
        call construct_spline1d(nx, delx, ddum1, bcx)
        f(1,:,j)=ddum1(1,:) ! the function value may be altered due to smoothing
        f(2,:,j)=ddum1(2,:)
        f(3,:,j)=ddum1(3,:)
     enddo
  endif

  do i = 1, nx
     do s = 1, 3
        ddum2(1,:)=f(s,i,:)
        call construct_spline1d(ny, dely, ddum2, bcy)
        f(s,i,:)=ddum2(1,:)
        f(s+3,i,:)=ddum2(2,:)
        f(s+6,i,:)=ddum2(3,:)
     enddo
  enddo
end subroutine construct_spline2d

! generic 3d quadratic spline
subroutine construct_spline3d(nx,ny,nz,delx,dely,delz,f,bcx,bcy,bcz)
  integer, intent(in):: nx, ny, nz, bcx, bcy, bcz
  real(lk), intent(in):: delx, dely, delz
  real(lk), intent(inout):: f(27,nx,ny,nz)

  integer i,j,k,s
  real(lk) temp1d(3,nz), temp2d(9,nx,ny)

! construct 2D spline on each x-y plane
! periodic in z will be enforced if bcz==2

  if (bcz==2) then
     do i=1,nz-1
        temp2d(1,:,:) = f(1,:,:,i)
        call construct_spline2d(nx,ny,delx,dely,temp2d,bcx,bcy)
        f(1:9,:,:,i) = temp2d(:,:,:)
     enddo
     f(1:9,:,:,nz)=f(1:9,:,:,1)
  else
     do i=1,nz
        temp2d(1,:,:) = f(1,:,:,i)
        call construct_spline2d(nx,ny,delx,dely,temp2d,bcx,bcy)
        f(1:9,:,:,i) = temp2d(:,:,:)
     enddo
  endif

! generate spline coefficient on z direction
  do i=1,nx
     do j=1,ny
        do s=1,9
           temp1d(1,:) = f(s,i,j,:)
           call construct_spline1d(nz, delz, temp1d, bcz)
           f(s,i,j,:) = temp1d(1,:)
           f(s+9,i,j,:) = temp1d(2,:)
           f(s+18,i,j,:) = temp1d(3,:)
        enddo
     enddo
  enddo
  
end subroutine construct_spline3d

! construct 1D spline with on x=[0,xmax], grid x_i = (i-1)*delx, delx=xmax/(nsp-1)
! in domain i, y = y(1,i) + y(2,i)*delx + y(3,i)*delx**2
subroutine construct_spline0(iflag,nsp,delx,y)
  integer i,nsp,ipp,iflag
  real(lk) delx,y(3,nsp)

! first point
  if(iflag==0)then
! iflag=0: first point being y=y1+y2*x+y3*x*x
     y(2,1)=(4.0_lk*y(1,2)-y(1,3)-3.0_lk*y(1,1))/(2.0_lk*delx)
     y(3,1)=(y(1,2)-y(1,1)-y(2,1)*delx)/(delx*delx)

  elseif(iflag==1)then
! iflag=1: first point being linear function y=y_1+y2*x
     y(2,1)=(y(1,2)-y(1,1))/delx
     y(3,1)=0.0_lk

  elseif(iflag==2)then
! iflag=2: first point being quadratic function y=y1+y3*x*x
     y(2,1)=0.0_lk
     y(3,1)=(y(1,2)-y(1,1))/(delx*delx)
  endif

  do i=2,nsp-2
     ipp=min(i+2,nsp)
     y(2,i)=-y(2,i-1)+2.0_lk*(y(1,i)-y(1,i-1))/delx

! smooth f1
     y(1,i+1)=0.5_lk*delx*y(2,i)+0.25_lk*y(1,ipp)+0.75_lk*y(1,i)
  enddo

  y(2,nsp-1)=-y(2,nsp-2)+2.0_lk*(y(1,nsp-1)-y(1,nsp-2))/delx
  y(2,nsp)=-y(2,nsp-1)+2.0_lk*(y(1,nsp)-y(1,nsp-1))/delx

  do i=2,nsp-1
     y(3,i)=(y(2,i+1)-y(2,i))/(2.0_lk*delx)
  enddo

! last point is not used;
  y(3,nsp)=0.0_lk

end subroutine construct_spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Construct 1D Quadratic Spline using periodic boundary condition
! Written by Lei Shi, Oct 16, 2016

! Check the documentation for a detailed discussion
subroutine construct_spline_periodic(nsp,delx,y)
  integer, intent(in) :: nsp
  real(lk), intent(in) :: delx
  real(lk), intent(inout):: y(3, nsp)
! local variables
  integer :: i
  real(lk) :: dxinv
! parameter to check if two real number are equal
  real(lk), parameter:: tol = 1D-5

! Periodic boundary condition requires the first order derivative to be continuous at the end points
! Check the documentation for a detailed discussion about the spline algorithm

dxinv = 1.0_lk/delx

! First, let's check if the input data satisfies the requirements
! Check grid point number
if (mod(nsp, 2) /= 0) then
   write(gtcout,*) 'Periodic Spline Error: Total number of data must be even. nsp=',nsp,' doesn''t work. Use normal spline or request a new input data set.'
   stop
endif

! check periodicity in raw data
if (abs(y(1,1)) >tol ) then
   ! use relative error
   if (abs((y(1, 1)-y(1, nsp))/y(1,1)) > tol) then
      write(gtcout,*) 'Periodic Spline Error: function values at end points must be the same. y(1)=',y(1,1),', y(n)=',y(1,nsp),' doesn''t work. Use normal spline or request a new input data set.'
      stop
   endif
else
   ! close to zero values, use absolute error
   if (abs(y(1, 1)-y(1, nsp)) > tol) then
      write(gtcout,*) 'Periodic Spline Error: function values at end points must be the same. y(1)=',y(1,1),', y(n)=',y(1,nsp),' doesn''t work. Use normal spline or request a new input data set.'
      stop
   endif
endif
! reset all the higher order spline coefficients
y(2:3,:) = 0.0_lk
! Now we let y = a + b*x + c*x**2 on each section, first step is to calculate b1

do i=1,nsp/2-1
   y(2,1) = y(2,1) + y(1,2*i) - y(1, 2*i+1)
enddo
y(2,1) = y(2,1)*2.0_lk*dxinv
! we use y(2,nsp) as a buffer to save the periodic b1
y(2, nsp)=y(2,1)
! Then, we can obtain bi one by one
do i=2, nsp-1
   y(2, i) = 2.0_lk*(y(1,i) - y(1,i-1))*dxinv - y(2,i-1)
enddo
! Then the ci's
do i=1, nsp-1
   y(3, i) = (y(2, i+1)-y(2, i))*0.5_lk*dxinv
enddo
end subroutine construct_spline_periodic


! inversion of spline0 y(x) to x(y)
subroutine invert_spline0(iflag,nsp,delx,dely,y,x)
  integer i,nsp,j,iflag
  real(lk) delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

! first point given by input
  x(1,1)=0.0_lk
! other points
  do i=2,nsp-1
! y grid
     ydum=dely*real(i-1)
! search x grid for ydum
     j=1
     do while (ydum>y(1,j+1))
        j=j+1
     enddo

! x(1,i)=j grid location + distance from j grid
     y0=y(1,j)
     y1=y(2,j)
     y2=y(3,j)
     if (abs(y2)>0.000001_lk) then
       x(1,i)=delx*real(j-1)+(sqrt(y1*y1+4.0_lk*y2*(ydum-y0))-y1)/(2.0_lk*y2)
     else
       x(1,i)=delx*real(j-1)+(ydum-y0)/y1
     endif

  enddo

! last point
  x(1,nsp)=x(1,1)+delx*real(nsp-1)

! spline fit x function
  if(iflag==0)call construct_spline0(0,nsp,dely,x)

end subroutine invert_spline0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! spline1 for cases with first point being function y=y1+y2*sqrt(x)+y3*x
subroutine construct_spline1(nsp,delx,f)
  integer i,nsp,ipp
  real(lk) delx,f(3,nsp)

! first point
  f(2,1)=(2.0_lk*f(1,2)-f(1,3)-f(1,1))/((2.0_lk-sqrt(2.0_lk))*sqrt(delx))
  f(3,1)=(f(1,2)-f(1,1)-f(2,1)*sqrt(delx))/delx

! second point
  f(2,2)=0.5_lk*f(2,1)/sqrt(delx)+f(3,1)
  f(3,2)=(f(1,3)-f(1,2)-delx*f(2,2))/(delx*delx)

  do i=3,nsp-2
     ipp=min(i+2,nsp)
     f(2,i)=-f(2,i-1)+2.0_lk*(f(1,i)-f(1,i-1))/delx

! smooth f1
     f(1,i+1)=0.5_lk*delx*f(2,i)+0.25_lk*f(1,ipp)+0.75_lk*f(1,i)
  enddo

  f(2,nsp-1)=-f(2,nsp-2)+2.0_lk*(f(1,nsp-1)-f(1,nsp-2))/delx
  f(2,nsp)=-f(2,nsp-1)+2.0_lk*(f(1,nsp)-f(1,nsp-1))/delx

  do i=3,nsp-1  !!change start from 1 to 3
     f(3,i)=(f(2,i+1)-f(2,i))/(2.0_lk*delx)
  enddo

! last point is not used;
  f(3,nsp)=0.0_lk

end subroutine construct_spline1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! inversion of spline1 y(x) to x(y)
subroutine invert_spline1(nsp,delx,dely,y,x)
  integer i,nsp,j
  real(lk) delx,dely,y(3,nsp),x(3,nsp),ydum,y0,y1,y2

! first point is given by inputs
  x(1,1)=0.0_lk
! other points
  do i=2,nsp-1
! y grid
     ydum=dely*real(i-1,lk)
! search x grid for ydum
     j=1
     do while (ydum>y(1,j+1))
        j=j+1
     enddo

! x(1,i)=j grid location + distance from j grid
     y0=y(1,j)
     y1=y(2,j)
     y2=y(3,j)

     if (abs(y2)>0.000001_lk) then
       x(1,i)=delx*real(j-1,lk)+(sqrt(y1*y1+4.0_lk*y2*(ydum-y0))-y1)/(2.0_lk*y2)
     else
       x(1,i)=delx*real(j-1,lk)+(ydum-y0)/y1
     endif

     if (j==1) x(1,i)=x(1,i)*x(1,i)

  enddo
! last point
  x(1,nsp)=delx*real(nsp-1,lk)

! spline fit x function
! call spline with first point being ~ sqrt(r)
  call construct_spline1(nsp,dely,x)

end subroutine invert_spline1

! Deprecated, should be removed after other parts are changed.
subroutine construct_spline2d_dep(nx,ny,delx,dely,y)
  integer i,j,s,nx,ny,ipp
  real(lk) delx,dely,y(9,nx,ny),ddum1(3,nx),ddum2(3,ny)

  ddum1=0.0_lk
  ddum2=0.0_lk

  do j = 1, ny
     ddum1(1,:)=y(1,:,j)     
     call construct_spline1(nx,delx,ddum1)
     y(1,:,j)=ddum1(1,:)
     y(2,:,j)=ddum1(2,:)
     y(3,:,j)=ddum1(3,:)
  enddo
  do i = 2, nx
     do s = 1, 3
        ddum2(1,:)=y(s,i,:)
        call construct_spline0(0,ny,dely,ddum2)
        y(s,i,:)=ddum2(1,:)
        y(s+3,i,:)=ddum2(2,:)
        y(s+6,i,:)=ddum2(3,:)
     enddo
  enddo
end subroutine construct_spline2d_dep
!!!!!!!!!!!!!

! Deprecated, should be removed after other parts are changed.
! Construct 2D spline with the second dimension periodic
! Written by Lei Shi, Oct 16, 2016
! See subroutine construct_spline_periodic for more details
subroutine construct_spline2d_periodic_dep(nx,ny,delx,dely,z)
  integer i,j,s,nx,ny,ipp
  real(lk) delx,dely,z(9,nx,ny),ddum1(3,nx),ddum2(3,ny)

  ddum1=0.0_lk
  ddum2=0.0_lk

  do j = 1, ny-1
     ddum1(1,:)=z(1,:,j)
     call construct_spline1(nx,delx,ddum1)
     z(1,:,j)=ddum1(1,:)
     z(2,:,j)=ddum1(2,:)
     z(3,:,j)=ddum1(3,:)
  enddo
  ! Values at periodic boundary are forced to be equal
  z(1,:,ny)=z(1,:,1)
  z(2,:,ny)=z(2,:,1)
  z(3,:,ny)=z(3,:,1)

  do i = 2, nx
     do s = 1, 3
        ddum2(1,:)=z(s,i,:)
        call construct_spline_periodic(ny,dely,ddum2)
        z(s,i,:)=ddum2(1,:)
        z(s+3,i,:)=ddum2(2,:)
        z(s+6,i,:)=ddum2(3,:)
     enddo
  enddo
end subroutine construct_spline2d_periodic_dep

!====================================================================
! 3D Data creation and spline subroutines for VMEC
! Written by Lei Shi, Nov 19, 2016
!==================================================================== 
! NOT IMPLEMENTED YET

! Deprecated, will be replaced by a 3D field generator, and apply
! construct_spline3d to created the coefficients

! TEMPORARY UPGRADED TO USE PERIODIC SPLINE IN ZETA
! By Lei Shi, Dec 13, 2016

! Issue #
! Slow data structure, try to change the data in (ndim, ny, nx) shape
subroutine spline3d_dep(nx,ny,nz,delx,dely,ndim,ntor,ycn,ysn,btemp,sgn)
  use precision
  use global_parameters
  implicit none

  integer i,j,k,n,s, bcx, bcy, bcz

  integer,intent(in) :: nx,ny,nz,ndim,sgn
  integer,intent(in) :: ntor(ndim)
  real(lk),intent(in) :: delx,dely,ycn(nx,ny,ndim),ysn(nx,ny,ndim)
  real(lk),intent(out) :: btemp(27,nx,ny,nz)
  real(lk) :: delz,zdum,dum,ddum(9,nx,ny),ddum3(3,nz+1)

  delz=2.0_lk*pi/real(nz,lk)
  ddum=0.0_lk
  ddum3=0.0_lk
  btemp=0.0_lk

  ! Determine the boundary condition based on the grid numbers
  ! x direction is always psi, needs BC No. 1
  bcx = 1
  ! ny is the total data points in poloidal direction
  if (mod(ny,2)==0) then
     bcy=2
  else
     bcy=0
     write(gtcout,*) "SPLINE WARNING: total poloidal spline grid number is odd, periodic spline not available. Non-periodic spline is used.", "lst = ", ny 
  endif
  ! nz is the total poloidal planes, thus is one less than the toroidal data points
  if (mod(nz,2)==0) then
     bcz=0
     write(gtcout,*) "SPLINE WARNING: total poloidal cross-section number is even, periodic spline in toroidal direction is not available. Non-periodic spline is used.", "mtoroidal = ", nz
  else
     bcz=2
  endif
  
! Now generate the 3D data and create spline coefficients
! for each toroidal grid the value is calculated as sum over cos and sin harmonics
! negative sign in zeta due to VMEC convention, sgn flipps directions of theta and zeta
! to ensure poloidal flux is positive
  if(sgn>0)then
    do k = 1, nz
      zdum= -real(k-1,lk)*delz 
      btemp(1,:,:,k)=0.0_lk
      do i = 1, nx
          do j = 1, ny 
             dum=0.0_lk
             do n = 1, ndim
                dum=dum+ycn(i,j,n)*cos(real(ntor(n),lk)*zdum)+ysn(i,j,n)*sin(real(ntor(n),lk)*zdum)
             enddo
             btemp(1,i,j,k)=dum
          enddo
      enddo
      ddum(1,:,:)=btemp(1,:,:,k)
      call construct_spline2d(nx,ny,delx,dely,ddum,bcx,bcy)
      btemp(1:9,:,:,k)=ddum(1:9,:,:)
    enddo
  else !VMEC original z-direction
    do k = 1, nz
      zdum= real(k-1,lk)*delz ! changed sign of zeta
      btemp(1,:,:,k)=0.0_lk
      do i = 1, nx
          do j = 1, ny ! change order for sgn<0
             dum=0.0_lk
             do n = 1, ndim
                dum=dum+ycn(i,ny-j+1,n)*cos(real(ntor(n),lk)*zdum)+ysn(i,ny-j+1,n)*sin(real(ntor(n),lk)*zdum)
             enddo
             btemp(1,i,j,k)=dum 
          enddo
      enddo
      ddum(1,:,:)=btemp(1,:,:,k)
      call construct_spline2d(nx,ny,delx,dely,ddum,bcx,bcy)
      btemp(1:9,:,:,k)=ddum(1:9,:,:)
    enddo
  endif

  do i = 1, nx
    do j = 1, ny
        do s = 1, 9
            ddum3(1,1:nz)=btemp(s,i,j,1:nz)
            ddum3(1,nz+1)=btemp(s,i,j,1) ! extra point to ensure periodicity in spline construction
            call construct_spline1d(nz+1,delz,ddum3,bcz)
            btemp(s,i,j,:)=ddum3(1,1:nz)
            btemp(s+9,i,j,:)=ddum3(2,1:nz)
            btemp(s+18,i,j,:)=ddum3(3,1:nz)
        enddo
    enddo
  enddo
end subroutine spline3d_dep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine cross(vec1,vec2,vec3)
  real(lk):: vec1(3),vec2(3),vec3(3)
  
  vec3(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
  vec3(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
  vec3(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
  
end subroutine cross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(lk) function dots(vec1,vec2,n)
  integer i,n
  real(lk) ::temp,vec1(n),vec2(n)
  temp=0.0_lk
  do i=1,n
    temp=temp+vec1(i)*vec2(i)
  enddo
  dots=temp
end function dots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine qdspline_r
! Purpose: For FRC field line mesh creation, given Z and psi values, 
!          calculate the corresponding R coordinate based on the psi(R,Z) 
!          spline function.
! Arguments: x, y, z: spline z(x, y) mesh and coefficients
!            yi, zi: the y, z value required for calculate xi
!            xi: output, the calculated x satisfying z(xi, yi)=zi to certain 
!                accuracy
!            zidum: output, the actual z at xi, yi, zidum=z(xi,yi)
!            iflag: 0 is to find root for Z function itself; 1 is to find root for derivative of Z function. (only 0 is available)

subroutine sppsiR1d(iflag,nx,ny,dx,dy,x,y,z,yi,zi,xi,zidum)
  use precision
  use global_parameters,only:gtcout
  implicit none

  integer,intent(in) :: nx,ny,iflag
  real(lk),intent(in) :: dx,dy,x(nx),y(ny),z(9,nx,ny),yi,zi
  real(lk),intent(out) :: xi,zidum

  integer :: i,j,ii
  real(lk) :: dpx,dpy,dx2,dy2,a,b,c,sign_r,zdum,phigh,plow,del_zdum,sperr,c1,c2
!Find the nearest j point (y direction)
  j=max(1,min(ny-1,ceiling(yi/dy)))
  dpy=yi-dy*real(j-1)
  dy2=dpy*dpy

!Find the nearest ii point(x direction),only valid for the monotonos Z funtion
  do i=1,nx-1
    c1=z(1,i,j)+z(4,i,j)*dpy+z(7,i,j)*dy2
    c2=z(1,i+1,j)+z(4,i+1,j)*dpy+z(7,i+1,j)*dy2
    sign_r=(zi-c1)*(zi-c2)
    if(sign_r .lt. 0)then
     ii=i
    exit
    elseif(sign_r==0 .and. c2>0)then
     ii=i
    exit
    endif
  enddo

!Find the x root of Z(x,y)=0 with known Z and y

if(iflag==0)then  !2D spline value
  a=z(3,ii,j)+z(6,ii,j)*dpy+z(9,ii,j)*dy2
  b=z(2,ii,j)+z(5,ii,j)*dpy+z(8,ii,j)*dy2
  c=z(1,ii,j)+z(4,ii,j)*dpy+z(7,ii,j)*dy2
  phigh=dx
  plow=0.0_lk
  dpx=0.5_lk*(phigh+plow)
  dx2=dpx*dpx
  zdum=a*dx2+b*dpx+c
  sperr=0.000001_lk*zi
  del_zdum=1.0_lk-zdum/zi
  i=1
  do while (abs(del_zdum)>sperr .and. i<1000)
    i=i+1
    if(del_zdum>0)then
    plow=dpx
    else
    phigh=dpx
    endif
    dpx=0.5_lk*(plow+phigh)
    dx2=dpx*dpx
    zdum=a*dx2+b*dpx+c
    del_zdum=1.0_lk-zdum/zi
  enddo

!Recalculate the Z value by using the x and y position
  xi=dpx+dx*real(ii-1)
  zidum= z(1,ii,j)    + z(2,ii,j)*dpx     + z(3,ii,j)*dx2     &
    + z(4,ii,j)*dpy  + z(5,ii,j)*dpx*dpy + z(6,ii,j)*dx2*dpy &
    + z(7,ii,j)*dy2 + z(8,ii,j)*dpx*dy2 + z(9,ii,j)*dx2*dy2
 
  if(abs(1.0_lk-zidum/zi)>1.0e3*sperr)then
  write(gtcout,*)"The x root of Z(x,y)=0 is not converged","error =",abs(1.0_lk-zidum/zi)*100,'%'
  endif

else
 
  write(gtcout,*)"IFLAG in sppsiR1d is only available for iflag==0"
  stop

endif

end subroutine sppsiR1d



end module spline_function
!###### MODULE SPLINE_FUNCTION END ######################################
