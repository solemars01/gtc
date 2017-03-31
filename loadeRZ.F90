subroutine loadeRZ
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use particle_tracking
  use LAMYRIDGEeq
  use cylindricalRZ
  use spline_function, only: spline2d,sppsiR1d
  use precision
  implicit none

  integer i,j,ij,ij0,ij1,ii,m,mtest,ierror,icount,mmin,mmax,nrsp,nzsp,marknum(numr*numz),nz,np
  real(lk) :: c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308
  real(lk) w_initial,rhoe,b,umax,pdum,tdum,rmin1,rmax1,rmin2,rmax2,zzmin,zzmax,rdum,zdum,dpr,dpz,del_r,del_z,rdum0,zdum0,vth,te,nonuni(numr*numz),anonuni,rg,psi_cy0,psi_cy1,deltapsi_cy,psidum,dr1,dr2,dr3,del_S,sdum(numr),zdum2,dz2,psifrc,del_p,jaco_frc,psidum1,psidum2,dr11,dr22,dr33

! allocate field memory
  allocate(pmarke(numr),markere(numr*numz),densitye(0:1,numr*numz),flowe(0:1,numr*numz),zonale(numr),zonalce(numr),STAT=mtest)

! allocate particle memory
  memax=me+100*ceiling(sqrt(real(me))) !ions array upper bound
  allocate(zelectron(1:nparame,memax),zelectron0(nparame,memax),zelectron1(nparame,memax),wzelectron(memax),wpelectron(ngyroe,memax),&
       jtelectron0(ngyroe,memax),jtelectron1(ngyroe,memax),wtelectron0(ngyroe,memax),&
       wtelectron1(ngyroe,memax),STAT=mtest)

  if(mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate ion: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
  pmarke=0.0
  markere=0.0
  densitye=0.0
  flowe=0.0
  zonale=0.0
  zonalce=0.0
  zelectron=0.0
  zelectron0=0.0
  zelectron1=0.0
  wzelectron=0.0
  wpelectron=0.0
  jtelectron0=0.0
  jtelectron1=0.0
  wtelectron0=0.0
  wtelectron1=0.0

  if(irun>0)return

  markere=0.0
  nonuni=0.0
  marknum=0
! number of marker per grid, Jacobian=(gq+I)/B^2
!$omp parallel do private(i,j,ij,ij1,jaco_frc,del_p,del_z)
  do i=2,numr-1
     do j=1,nz_frc(i)
        ij=igrid_frc(i)+j
        ij0=igrid_frc(i-1)+j
        ij1=igrid_frc(i+1)+j
        jaco_frc=RZ_reg(ij,1)/dpsifrcdr(ij)
        del_p=0.5*(psiEq0(ij1)-psiEq0(ij0))
        if(j==1 .or. j==nz_frc(i))then
        del_z=0.5*delta_Z
        else
        del_z=1.0*delta_Z
        endif
        markere(ij)=meshne_cy(ij)*jaco_frc*del_p*del_z   ! physical particle distribution profile
        nonuni(ij)=jaco_frc*del_p*del_z
     enddo
  enddo
!$omp parallel do private(j,ij,ij0,jaco_frc,del_p,del_z)
  do j=1,nz_frc(numr)
        ij=igrid_frc(numr)+j
        ij0=igrid_frc(numr-1)+j
        jaco_frc=RZ_reg(ij,1)/dpsifrcdr(ij)
        del_p=0.5*(psiEq0(ij)-psiEq0(ij0))
        if(j==1 .or. j==nz_frc(numr))then
        del_z=0.5*delta_Z
        else
        del_z=1.0*delta_Z
        endif
        markere(ij)=meshne_cy(ij)*jaco_frc*del_p*del_z   ! physical particle distribution profile
        nonuni(ij)=jaco_frc*del_p*del_z
   enddo

if(psi0_frc/=0.0)then
!$omp parallel do private(j,ij,ij1,jaco_frc,del_p,del_z)
  do j=1,nz_frc(1)
        ij=igrid_frc(1)+j
        ij1=igrid_frc(2)+j
        jaco_frc=RZ_reg(ij,1)/dpsifrcdr(ij)
        del_p=0.5*(psiEq0(ij1)-psiEq0(ij))
        if(j==1 .or. j==nz_frc(1))then
          del_z=0.5*delta_Z
        else
          del_z=1.0*delta_Z
        endif    
        markere(ij)=meshne_cy(ij)*jaco_frc*del_p*del_z   ! physical particle distribution profile
        nonuni(ij)=jaco_frc*del_p*del_z
     enddo
else
!$omp parallel do private(j,ij,psidum1,psidum2,del_z,dr1,dr2,dr3,dr11,dr22,dr33,del_S)
    do j=1,nz_frc(1)
        ij=igrid_frc(1)+j
        psidum1=psi_frc(1)+0.5*deltapsi_frc
        psidum2=psi_frc(1)
        del_z=0.5*delta_Z
       if(j==1)then
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psidum1,dr1,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax+del_z,psidum1,dr2,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psidum2,dr11,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax+del_z,psidum2,dr22,psidum)
         dr1=dr1-dr11
         dr2=dr2-dr22
         del_S=0.5*(dr1+dr2)*del_z*(dr11+0.5*dr1) !Jacobian volume
       elseif(j==nz_frc(1))then
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psidum1,dr1,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax-del_z,psidum1,dr2,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psidum2,dr11,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax-del_z,psidum2,dr22,psidum)
         dr1=dr1-dr11
         dr2=dr2-dr22
         del_S=0.5*(dr1+dr2)*del_z*(dr11+0.5*dr1) !Jacobian volume
       else
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax-del_z,psidum1,dr1,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax+del_z,psidum1,dr2,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax-del_z,psidum2,dr11,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax+del_z,psidum2,dr22,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psidum1,dr3,psidum)
         call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psidum2,dr33,psidum)
         dr1=dr1-dr11
         dr2=dr2-dr22
         dr3=dr3-dr33

         del_S=(dr1+dr2)*del_z*(dr33+0.5*dr3) !Jacobian volume
       endif
        markere(ij)=meshne_cy(ij)*del_S   ! physical particle distribution profile
        nonuni(ij)=del_S  !del_S is Jacobian volume here!!!! 
     enddo
endif

nonuni=1.0
!write(gtcout,*)"nonuni=",nonuni
!write(gtcout,*)"marknum=",marknum


! set uniform marker number per cell
  pdum=sum(nonuni)
!$omp parallel do private(i,j,ij)
  do i=1,numr
     do j=1,nz_frc(i)
        ij=igrid_frc(i)+j
        marknum(ij)=int(nonuni(ij)*real(me)/pdum) ! # of marker per cell
     enddo
  enddo

  pmarke=0.0
  do i=1,numr
     do j=1,nz_frc(i)
        ij=igrid_frc(i)+j
        pmarke(i)=pmarke(i)+real(marknum(ij))
     enddo
  enddo
  pmarke=pmarke*real(numberpe)

! marker number which is used for nonmalization, proportional to physcial density
  pdum=sum(markere)
!$omp parallel do private(i,j,ij)
  do i=1,numr
     do j=1,nz_frc(i)
        ij=igrid_frc(i)+j
        markere(ij)=(markere(ij)*real(me)/pdum) ! # of marker per cell
     enddo
  enddo
!$omp parallel do private(i,j,ij)
  do i=1,numr
     do j=1,nz_frc(i)
        ij=igrid_frc(i)+j
        markere(ij)=real(npartdom)*markere(ij)
       if(markere(ij)<1.0e-6)then
          markere(ij)=0.0
       else
          markere(ij)=1.0/markere(ij) !to avoid divide operation
       endif
     enddo
  enddo


!$omp parallel do private(i,j,ij) 
  do i=1,numr
     do j=1,nz_frc(i)
        ij=igrid_frc(i)+j
        nonuni(ij)=1.0/(real(marknum(ij))*markere(ij)*real(npartdom))
     enddo
  enddo


write(gtcout,*)"average nonuni=",sum(nonuni)/(real(numr)*real(numz))



!write(gtcout,*)"markere",markere
  me=sum(marknum)
  w_initial=1.0e-3_lk
  umax=4.0

! random # uniformly distributed between 0 and 1
  call random_number(zelectron(1,1:me))
  call random_number(zelectron(2,1:me))
  call random_number(zelectron(3,1:me))
  call random_number(zelectron(4,1:me))
  call random_number(zelectron(5,1:me))
  call random_number(zelectron(6,1:me))

  mmin=0
  mmax=0
  do i=1,numr
   do j=1,nz_frc(i)
     ij=igrid_frc(i)+j

     zzmin=RZ_reg(ij,2)
     zzmax=RZ_reg(ij+1,2)
!! range of particle #
     mmin=mmax+1
     mmax=mmin-1+marknum(ij)
     if(j==1 .or. j==nz_frc(numr))then
       del_z=0.5*delta_Z
     else
       del_z=1.0*delta_Z
     endif
     if(i==1)then
       del_p=0.5*(psiEq0(igrid_frc(i+1)+j)-psiEq0(ij))
     elseif(i==numr)then
       del_p=-0.5*(psiEq0(ij)-psiEq0(igrid_frc(i-1)+j))
     else
       del_p=0.5*(psiEq0(igrid_frc(i+1)+j)-psiEq0(igrid_frc(i-1)+j))
     endif

! load ion in radial direction: linear in psi within a grid cell
!$omp parallel do private(m,rdum,zdum)
      do m=mmin,mmax
       if(i==1 .or. i==numr)then
        rdum=zelectron(1,m)
       else
        rdum=zelectron(1,m)-0.5
       endif
       if(j==1)then
        zdum=zelectron(3,m)
       elseif(j==nz_frc(i))then
        zdum=-zelectron(3,m)
       else
        zdum=zelectron(3,m)-0.5
       endif
        zelectron(3,m)=zzmin+zdum*del_z

        psidum=psi_frc(i)+rdum*del_p
        call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,zelectron(3,m)+zmax,psidum,zelectron(1,m),psidum)

        zelectron(5,m)=w_initial
      enddo
 if(.false.)then
         if(mype==0 .and. ij==20)then
            open(9991,file='time_ext.out',status='replace')
         endif
         if(mype==0 .and. ij==20)then
            write(9991,*)zelectron(1,mmin:mmax)
         endif
         if(mype==0 .and. ij==20)then
            close(9991)
         endif

         if(mype==0 .and. ij==20)then
            open(9992,file='phi_ext.out',status='replace')
         endif
         if(mype==0 .and. ij==20)then
            write(9992,*)zelectron(3,mmin:mmax)
         endif
         if(mype==0 .and. ij==20)then
            close(9992)
         endif
 endif
   enddo
  enddo
!write(gtcout,*)"deltas_frc=",deltas_frc(20)
!write(gtcout,*)"sin_frc=",sin_frc(20)
!write(gtcout,*)"zelectron(8,100)=",zelectron(8,100)
!! load electron uniformly in zeta
!$omp parallel do private(m,tdum)
  do m=1,me
     tdum=zelectron(2,m)*deltaz
     zelectron(2,m)=zeta0+tdum
  enddo

write(gtcout,*)"deltaz =",deltaz
if(.false.)then
         if(mype==0)then
            open(9991,file='time_ext.out',status='replace')
         endif
         if(mype==0)then
            write(9991,*)zelectron(1,1:me)
         endif
         if(mype==0)then
            close(9991)
         endif

         if(mype==0)then
            open(9992,file='phi_ext.out',status='replace')
         endif
         if(mype==0)then
            write(9992,*)zelectron(3,1:me)
         endif
         if(mype==0)then
            close(9992)
         endif
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(.true.)then
!$omp parallel do private(m,zdum)
      do m=1,me
        zdum=zelectron(4,m)
        zelectron(4,m)=zelectron(4,m)-0.5
        zelectron0(4,m)=sign(1.0_lk,zelectron(4,m))
        zelectron(4,m)=sqrt(max(1.0e-20_lk,log(1.0_lk/max(1.0e-20_lk,zelectron(4,m)**2))))
        zelectron(4,m)=zelectron(4,m)-(c0+c1*zelectron(4,m)+c2*zelectron(4,m)**2)/&
          (1.0+d1*zelectron(4,m)+d2*zelectron(4,m)**2+d3*zelectron(4,m)**3)
        if(zelectron(4,m)>umax)zelectron(4,m)=zdum
        zelectron(4,m)=zelectron0(4,m)*zelectron(4,m)
        zelectron(6,m)=max(1.0e-20_lk,min(umax*umax,-log(max(1.0e-20_lk,zelectron(6,m)))))
      enddo

  !       if(mype==0)then
  !          open(9991,file='time_ext.out',status='replace')
  !       endif
  !       if(mype==0)then
  !          write(9991,*)zelectron(4,1:me)
  !       endif
  !       if(mype==0)then
  !          close(9991)
  !       endif
 
  !       if(mype==0)then
  !          open(9992,file='phi_ext.out',status='replace')
  !       endif
  !       if(mype==0)then
  !          write(9992,*)zelectron(2,1:me)
  !       endif
  !       if(mype==0)then
  !          close(9992)
  !       endif

!!$omp parallel do private(m,rdum,zdum,b,te,psidum,rhoe,vth)
  do m=1,me
     rdum=zelectron(1,m)
     zdum=zelectron(3,m)
     b = spline2d(rdum,zdum+zmax,0,lsr,lsz,spdR,spdZ,B0mag_Eq,0,0)
     te = spline2d(rdum,zdum+zmax,0,lsr,lsz,spdR,spdZ,lrte_Eq,0,0)
     psidum = spline2d(rdum,zdum+zmax,0,lsr,lsz,spdR,spdZ,psirz_Eq,0,0)
 !    b=1
 !    te=meshte_cy(1)
     rhoe=sqrt(te)*sqrt(aelectron)
     vth=sqrt(te)/sqrt(aelectron)
     zelectron(4,m)=vth*zelectron(4,m)
     zelectron(6,m)=vth*sqrt(aelectron*zelectron(6,m)/b)
    
     zelectron(7,m)=psidum
     zelectron(8,m)=0.5*aelectron*zelectron(4,m)*zelectron(4,m)+zelectron(6,m)*zelectron(6,m)*b

  enddo


endif


  if(.true.)then
      if(mype==0)then
          zelectron(1,1)=1.8
          zelectron(2,1)=zeta0+0.5*deltaz
          zelectron(3,1)=0.0
          zelectron(4,1)=1.662  !Te=29.9eV
          zelectron(6,1)=0.0205   !Te=29.9eV at R=1.8,Z=0.0
          b = spline2d(zelectron(1,1),zelectron(3,1)+zmax,0,lsr,lsz,spdR,spdZ,B0mag_Eq,0,0)
          psidum = spline2d(zelectron(1,1),zelectron(3,1)+zmax,0,lsr,lsz,spdR,spdZ,psirz_Eq,0,0)
          zelectron(7,1)=psidum
          zelectron(8,1)=0.5*aelectron*zelectron(4,1)*zelectron(4,1)+zelectron(6,1)*zelectron(6,1)*b

       write(gtcout,*)"zelectron1","R=",zelectron(1,1),"\zeta=",zelectron(2,1),"Z=",zelectron(3,1),"vpara=",zelectron(4,1),"\mu",zelectron(6,1)
      endif
  endif
       if(.false.)then
         if(mype==0)then
            open(9991,file='time_ext.out',status='replace')
         endif
         if(mype==0)then
            write(9991,*)zelectron(4,1:me)
         endif
         if(mype==0)then
            close(9991)
         endif

         if(mype==0)then
            open(9992,file='phi_ext.out',status='replace')
         endif
         if(mype==0)then
            write(9992,*)zelectron(6,1:me)
         endif
         if(mype==0)then
            close(9992)
         endif
       endif
  if(mype==0)then
     write(gtcout,*)'me=',me
     call FLUSH(gtcout)
  endif

end subroutine loadeRZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

