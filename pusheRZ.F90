subroutine pusheRZ
  use precision
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use LAMYRIDGEeq
  use cylindricalRZ
  use spline_function, only: spline2d,sppsiR1d
  implicit none

  integer i,j,m,igyro,ij
  real(lk) dtime,wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,&
       e1,e2,e3,cmratio,cinv,wpgc(3,me),rdot,zdot,zetadot,vparadot,wdot,energy0,gyrodum,perturb,&
       rdum,rdum_inv,zdum,thdum,vparadum,wdum,mudum,nrsp,dpr,dr2,nzsp,dpz,dz2,psidum,&
       psifrc,gradpsi_r,gradpsi_z,grad2psi_rr,grad2psi_rz,grad2psi_zr,grad2psi_zz,&
       br,bz,b,b_inv,gradbr_r,gradbr_z,gradbz_r,gradbz_z,gradB0R,gradB0Z,omegaP,spdR_inv,spdZ_inv

  cmratio=qelectron/aelectron
  perturb=real(nonlinear)
  spdR_inv=1.0/spdR
  spdZ_inv=1.0/spdZ
 
  if(irk==1)then
! 1st step of Runge-Kutta method
     dtime=0.5*tstep
     
!$omp parallel do private(m)
       do m=1,me
          zelectron0(1:nparame,m)=zelectron(1:nparame,m)
       enddo
     
! 2nd step of Runge-Kutta method
  else
     dtime=tstep
  endif

! gather e_field using ngyroi-point gyro-averaging
  gyrodum=1.0/real(ngyroe)

! update GC position. 
!$omp parallel do private(m,rdum,rdum_inv,zdum,thdum,vparadum,wdum,mudum,nrsp,dpr,dr2,nzsp,dpz,dz2,&
!$omp& psifrc,gradpsi_r,gradpsi_z,grad2psi_rr,grad2psi_rz,grad2psi_zr,grad2psi_zz,br,bz,b,b_inv,gradbr_r,gradbr_z,gradbz_r,gradbz_z,gradB0R,gradB0Z,energy0,&
!$omp& omegaP,rdot,zetadot,zdot,vparadot,wdot)
  do m=1,me
!~~~~~~~~~~initial particle quantities~~~~~~~~~~~~~~~~~~~~~~~
    rdum = max(zelectron(1,m),min_r)
    rdum_inv=1.0/rdum
    zdum = zelectron(3,m)
    thdum = zelectron(2,m)
    vparadum = zelectron(4,m)
    wdum = zelectron(5,m)
    mudum = zelectron(6,m)*zelectron(6,m)

    nrsp=max(1,min(lsr-1,ceiling(rdum*spdR_inv)))
    dpr=rdum-spdR*real(nrsp-1)
    dr2=dpr*dpr

    nzsp=max(1,min(lsz-1,ceiling((zdum+zmax)*spdZ_inv)))
    dpz=zdum+zmax-spdZ*real(nzsp-1)
    dz2=dpz*dpz

!!!!!!!!!!!!!!!!!!!!!!!!! psi !!!!!!!!!!!!!!!!!!!!!!!!!!
    psifrc=psirz_Eq(1,nrsp,nzsp)     + psirz_Eq(2,nrsp,nzsp)*dpr     + psirz_Eq(3,nrsp,nzsp)*dr2 &
        +psirz_Eq(4,nrsp,nzsp)*dpz + psirz_Eq(5,nrsp,nzsp)*dpr*dpz + psirz_Eq(6,nrsp,nzsp)*dr2*dpz &
        +psirz_Eq(7,nrsp,nzsp)*dz2 + psirz_Eq(8,nrsp,nzsp)*dpr*dz2 + psirz_Eq(9,nrsp,nzsp)*dr2*dz2

    gradpsi_r= psirz_Eq(2,nrsp,nzsp)+ 2.0*psirz_Eq(3,nrsp,nzsp)*dpr &
       + psirz_Eq(5,nrsp,nzsp)*dpz + 2.0*psirz_Eq(6,nrsp,nzsp)*dpr*dpz &
       + psirz_Eq(8,nrsp,nzsp)*dz2 + 2.0*psirz_Eq(9,nrsp,nzsp)*dpr*dz2

    gradpsi_z=psirz_Eq(4,nrsp,nzsp) + psirz_Eq(5,nrsp,nzsp)*dpr + psirz_Eq(6,nrsp,nzsp)*dr2 &
      +2.0*psirz_Eq(7,nrsp,nzsp)*dpz + 2.0*psirz_Eq(8,nrsp,nzsp)*dpr*dpz + 2.0*psirz_Eq(9,nrsp,nzsp)*dr2*dpz

    grad2psi_rr= 2.0*psirz_Eq(3,nrsp,nzsp)+ 2.0*psirz_Eq(6,nrsp,nzsp)*dpz &
                 + 2.0*psirz_Eq(9,nrsp,nzsp)*dz2
    grad2psi_rz= psirz_Eq(5,nrsp,nzsp) + 2.0*psirz_Eq(6,nrsp,nzsp)*dpr &
       + 2.0*psirz_Eq(8,nrsp,nzsp)*dpz + 4.0*psirz_Eq(9,nrsp,nzsp)*dpr*dpz

    grad2psi_zr= psirz_Eq(5,nrsp,nzsp) + 2.0*psirz_Eq(6,nrsp,nzsp)*dpr &
       + 2.0*psirz_Eq(8,nrsp,nzsp)*dpz + 4.0*psirz_Eq(9,nrsp,nzsp)*dpr*dpz
    grad2psi_zz= 2.0*psirz_Eq(7,nrsp,nzsp) + 2.0*psirz_Eq(8,nrsp,nzsp)*dpr + 2.0*psirz_Eq(9,nrsp,nzsp)*dr2
!!!!!!!!!!!!!!!!!!!!!!!!! psi !!!!!!!!!!!!!!!!!!!!!!!!!!
    br=-gradpsi_z*rdum_inv
    bz=gradpsi_r*rdum_inv
    b=sqrt(br**2+bz**2)
    b_inv=1.0/b
!unit magnetic field vector
    br=br*b_inv
    bz=bz*b_inv
!space derivative of Br and Bz
    gradbr_r=-grad2psi_zr*rdum_inv+gradpsi_z*rdum_inv*rdum_inv
    gradbr_z=-grad2psi_zz*rdum_inv
    gradbz_r=grad2psi_rr*rdum_inv-gradpsi_r*rdum_inv*rdum_inv
    gradbz_z=grad2psi_rz*rdum_inv
!space derivative of B0
    gradB0R=br*gradbr_r+bz*gradbz_r
    gradB0Z=br*gradbr_z+bz*gradbz_z
!space derivative of unit br and bz
    gradbr_z=gradbr_z*b_inv-br*gradB0Z*b_inv
    gradbz_r=gradbz_r*b_inv-bz*gradB0R*b_inv
    gradbr_r=gradbr_r*b_inv-br*gradB0R*b_inv
    gradbz_z=gradbz_z*b_inv-bz*gradB0Z*b_inv

!cyclotron frequency
    omegaP=qelectron*b/aelectron
    energy0=0.5*aelectron*vparadum*vparadum+mudum*b


    zelectron(7,m)=psifrc
    zelectron(8,m)=energy0

    wdot=0.0



    rdot=vparadum*br
    zetadot=((bz*gradB0R-br*gradB0Z)*mudum/(omegaP*aelectron)& !gradB drift
            +(gradbr_Z-gradbz_R)*vparadum*vparadum/omegaP)*rdum_inv  !curvature B drift
    zdot=vparadum*bz

    vparadot=-(br*gradB0R+bz*gradB0Z)*mudum/aelectron



    zelectron(1,m) = zelectron0(1,m)+dtime*rdot
    zelectron(2,m) = zelectron0(2,m)+dtime*zetadot
    zelectron(3,m) = zelectron0(3,m)+dtime*zdot
    zelectron(4,m) = zelectron0(4,m)+dtime*vparadot
    zelectron(5,m)=zelectron0(5,m)+dtime*wdot


    zelectron(2,m)=modulo(zelectron(2,m),torbound)

  enddo

if(irk==2)then

! out of boundary particle
!$omp parallel do private(m) 
    do m=1,me
    
      if(zelectron(7,m) > psi1_frc .or. zelectron(7,m)<psi0_frc)then    
        psidum = spline2d(zelectron0(1,m),zelectron0(3,m)+zmax,0,lsr,lsz,spdR,spdZ,psirz_Eq,0,0)
        call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,zelectron(3,m)+zmax,psidum,zelectron(1,m),psidum)
      endif

      if(zelectron(3,m) > Zsimulation1)then
        zelectron(1,m)=zelectron0(1,m)
        zelectron(3,m)=zelectron0(3,m)
        zelectron(4,m)=-zelectron0(4,m)
      elseif(zelectron(3,m) < Zsimulation0)then
        zelectron(1,m)=zelectron0(1,m)
        zelectron(3,m)=zelectron0(3,m)
        zelectron(4,m)=-zelectron0(4,m)
      endif
    enddo

endif

end subroutine pusheRZ
