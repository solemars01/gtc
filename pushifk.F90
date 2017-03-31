subroutine pushifk
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use magnetic_island
  implicit none

  integer i,j,k,m,ii,ierror,igyro,ij,irsmooth,icount,iir(mi),isp,jst
  real(lk) rdum,tdum,dtime,q,rg,b,g,gp,ri,rip,dbdp,dbdt,dedb,deni,gqi,&
       deltaf,upara,energy,kappa,dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,wdot,wdrift,&
       delr,dmark(0:mpsi),xdum,ydum,xdot,ydot,tem,&
       wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,ti_inv,epara,wpara,&
       perturb,paxis,pdum,wdrive,psitmp,thetatmp,zetatmp,&
       e1,e2,e3,cmratio,cinv,wpgc(3,mi),dtem(0:mpsi),dden(0:mpsi),ddum(0:mpsi),&
       vdrtmp(0:mpsi),diagtmp(mpdiag),data1dtmp(0:mpsi,mpdata1d),dpx,dp2,dtx,dt2,&
       b1,b2,b3,b4,b5,wbgc(5,mi),i1,i2,i3,i4,wigc(4,mi),dapdp,dapdt,dapdz,dapidt,dapidz,vap,gyrodum,&
       lam,dlamdp,dlamdt,dlamdz,plampt,delp(mpsi),dp1,kappan,kappat,fullf,&
       upara0,dupara0dp,kappav,dedb0,energy0,b_inv,majorr,dmajorrdp,r,angmom,er,cost,sint,&
       te,ne,kne,uperp,uperp0,angle,tandalpha,sindalpha,delvx,delvy,vxtemp,vytemp,utemp,drgdp,wisland,&
       Ex,Ey,vx1,vy1,vx2,vy2

  paxis=0.0
  if(mpsilow>1)paxis=psimesh(mpsilow)
  delr=1.0/deltar
  cmratio=qion/aion
  cinv=1.0/qion
  perturb=real(nonlinear)
  vdrtmp=0.0
  delp=1.0/deltap

  if(irk==1)then
! 1st step of Runge-Kutta method
     dtime=0.5*tstep
     
!$omp parallel do private(m)
       do m=1,mi
          zion0(1:nparami,m)=zion(1:nparami,m)
       enddo
     
! 2nd step of Runge-Kutta method
  else
     dtime=tstep
!by zhs
!     if(nonlinear==1)vdrtmp=pfluxi
  endif

! gather e_field using ngyroi-point gyro-averaging
  gyrodum=1.0/real(ngyroi)

! electrostatic fluctuations
  if(magnetic==0)then
!$omp parallel do private(m,igyro,e1,e2,e3,wz1,wz0,wp1,wp0,wt00,wt10,wt01,wt11,ij)
     do m=1,mi
        e1=0.0
        e2=0.0
        e3=0.0
        wz1=wzion(m)                !weight for upper toroidal grid
        wz0=1.0-wz1                 !weight for lower toroidal grid

        do igyro=1,ngyroi
           wp1=wpion(igyro,m)       !outer flux surface
           wp0=1.0-wp1              !inner flux surface

           wt10=wp0*wtion0(igyro,m) !upper poloidal grid on inner flux surface
           wt00=wp0-wt10            !lower poloidal grid on inner flux surface

           wt11=wp1*wtion1(igyro,m) !upper poloidal grid on outer flux surface
           wt01=wp1-wt11            !lower poloidal grid on outer flux surface

           ij=jtion0(igyro,m)       !lower poloidal grid on inner flux surface
           e1=e1+wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

           ij=ij+1                  !upper poloidal grid on inner flux surface
           e1=e1+wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

           ij=jtion1(igyro,m)       !lower poloidal grid on outer flux surface
           e1=e1+wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))

           ij=ij+1                  !upper poloidal grid on outer flux surface
           e1=e1+wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
        enddo

        wpgc(1,m)=gyrodum*e1
        wpgc(2,m)=gyrodum*e2
        wpgc(3,m)=gyrodum*e3
        wbgc(1,m)=0.0
        wbgc(2,m)=0.0
        wbgc(3,m)=0.0
        wbgc(4,m)=0.0
        wbgc(5,m)=0.0
     enddo

  else
! electromagnetic fields
!$omp parallel do private(m,igyro,e1,e2,e3,b1,b2,b3,b4,b5,wz1,wz0,wp0,wp1,wt00,wt10,wt01,wt11,ij)
     do m=1,mi
        e1=0.0
        e2=0.0
        e3=0.0
        b1=0.0
        b2=0.0
        b3=0.0
        b4=0.0
        b5=0.0
        wz1=wzion(m)                !weight for upper toroidal grid
        wz0=1.0-wz1                 !weight for lower toroidal grid

        do igyro=1,ngyroi
           wp1=wpion(igyro,m)       !outer flux surface
           wp0=1.0-wp1              !inner flux surface

           wt10=wp0*wtion0(igyro,m) !upper poloidal grid on inner flux surface
           wt00=wp0-wt10            !lower poloidal grid on inner flux surface

           wt11=wp1*wtion1(igyro,m) !upper poloidal grid on outer flux surface
           wt01=wp1-wt11            !lower poloidal grid on outer flux surface

           ij=jtion0(igyro,m)       !lower poloidal grid on inner flux surface
           e1=e1+wt00*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt00*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt00*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           b1=b1+wt00*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wt00*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wt00*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wt00*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
           b5=b5+wt00*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

           ij=ij+1                  !upper poloidal grid on inner flux surface
           e1=e1+wt10*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt10*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt10*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           b1=b1+wt10*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wt10*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wt10*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wt10*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
           b5=b5+wt10*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

           ij=jtion1(igyro,m)       !lower poloidal grid on outer flux surface
           e1=e1+wt01*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt01*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt01*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           b1=b1+wt01*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wt01*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wt01*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wt01*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
           b5=b5+wt01*(wz0*sapara(0,ij)+wz1*sapara(1,ij))

           ij=ij+1                  !upper poloidal grid on outer flux surface
           e1=e1+wt11*(wz0*gradphi(1,0,ij)+wz1*gradphi(1,1,ij))
           e2=e2+wt11*(wz0*gradphi(2,0,ij)+wz1*gradphi(2,1,ij))
           e3=e3+wt11*(wz0*gradphi(3,0,ij)+wz1*gradphi(3,1,ij))
           b1=b1+wt11*(wz0*gradapara(1,0,ij)+wz1*gradapara(1,1,ij))
           b2=b2+wt11*(wz0*gradapara(2,0,ij)+wz1*gradapara(2,1,ij))
           b3=b3+wt11*(wz0*gradapara(3,0,ij)+wz1*gradapara(3,1,ij))
           b4=b4+wt11*(wz0*gradphieff(3,0,ij)+wz1*gradphieff(3,1,ij))
           b5=b5+wt11*(wz0*sapara(0,ij)+wz1*sapara(1,ij))
        enddo

        wpgc(1,m)=gyrodum*e1
        wpgc(2,m)=gyrodum*e2
        wpgc(3,m)=gyrodum*e3
        wbgc(1,m)=gyrodum*b1
        wbgc(2,m)=gyrodum*b2
        wbgc(3,m)=gyrodum*b3
        wbgc(4,m)=gyrodum*b4
        wbgc(5,m)=gyrodum*b5
     enddo
  endif

! update GC position. Assuming psi0>spdpsi; will be relaxed later
!$omp parallel do private(m,psitmp,thetatmp,isp,dpx,dp2,jst,dtx,dt2,rg,q,g,gp,ri,rip,&
!$omp& b,dbdp,dbdt,ii,dedb,deni,gqi,upara,energy,dptdp,dptdt,dptdz,&
!$omp& dapdp,dapdt,dapdz,epara,vdr,vap,wdrive,wpara,wdrift,wdot,rdot,pdot,tdot,zdot,kappa,&
!$omp& lam,dlamdp,dlamdt,dlamdz,plampt,dp1,kappan,kappat,ti_inv,&
!$omp& upara0,dupara0dp,kappav,b_inv,majorr,dmajorrdp,dedb0,energy0,er,te,ne,kne,&
!$omp& uperp,angle,delvx,delvy,uperp0,vxtemp,vytemp,tandalpha,sindalpha,utemp,drgdp,&
!$omp& Ex,Ey,vx1,vy1,vx2,vy2)
  do m=1,mi
     psitmp=zion(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zion(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

! 1D spline in psi
     q=    qpsi(1,isp)   +qpsi(2,isp)*dpx   +qpsi(3,isp)*dp2
     g=    gpsi(1,isp)   +gpsi(2,isp)*dpx   +gpsi(3,isp)*dp2
     gp=                  gpsi(2,isp)       +gpsi(3,isp)*dpx*2.0
     ri=   cpsi(1,isp)   +cpsi(2,isp)*dpx   +cpsi(3,isp)*dp2
     rip=                 cpsi(2,isp)       +cpsi(3,isp)*dpx*2.0
     upara0=ropp(1,isp)  +ropp(2,isp)*dpx   +ropp(3,isp)*dp2
     dupara0dp=           ropp(2,isp)       +ropp(3,isp)*dpx*2.0
     er=   erpp(1,isp)   +erpp(2,isp)*dpx   +erpp(3,isp)*dp2

! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     drgdp=              rgpsi(2,isp)      +rgpsi(3,isp)*dpx*2.0

! 2D spline in (psi, theta)
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

     dbdt= bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dtx*2.0

     if(isp==1)then
        dbdp=(bsp(2,isp,jst)+bsp(5,isp,jst)*dtx+bsp(8,isp,jst)*dt2)*0.5/dpx+&
             +bsp(3,isp,jst)+bsp(6,isp,jst)*dtx+bsp(9,isp,jst)*dt2
     else
        dbdp= bsp(2,isp,jst)    +bsp(3,isp,jst)*dpx*2.0+ &
             (bsp(5,isp,jst)    +bsp(6,isp,jst)*dpx*2.0)*dtx+ &
             (bsp(8,isp,jst)    +bsp(9,isp,jst)*dpx*2.0)*dt2
     endif
     majorr= xsp(1,isp,jst)+xsp(2,isp,jst)*dpx+xsp(3,isp,jst)*dp2+ &
          (xsp(4,isp,jst)+xsp(5,isp,jst)*dpx+xsp(6,isp,jst)*dp2)*dtx+ &
          (xsp(7,isp,jst)+xsp(8,isp,jst)*dpx+xsp(9,isp,jst)*dp2)*dt2

     dmajorrdp= xsp(2,isp,jst)    +xsp(3,isp,jst)*dpx*2.0+ &
          (xsp(5,isp,jst)    +xsp(6,isp,jst)*dpx*2.0)*dtx+ &
          (xsp(8,isp,jst)    +xsp(9,isp,jst)*dpx*2.0)*dt2

     b_inv=1.0/b

! perturbed electric field
     dptdp=wpgc(1,m)
     dptdt=wpgc(2,m)
     dptdz=wpgc(3,m)-wpgc(2,m)/q

! perturbed magnetic field
     dapdp=wbgc(1,m)
     dapdt=wbgc(2,m)
     dapdz=wbgc(3,m)-wbgc(2,m)/q

! lambda=apara*b_inv
     lam=wbgc(5,m)*b_inv
     dlamdp=(dapdp-lam*dbdp)*b_inv
     dlamdt=(dapdt-lam*dbdt)*b_inv
     dlamdz=dapdz*b_inv

     upara0=upara0*majorr
     dupara0dp=dupara0dp*majorr+upara0*dmajorrdp

     upara=zion(4,m)*b*cmratio
     dedb=zion(4,m)*zion(4,m)*b*cmratio+cinv*zion(6,m)*zion(6,m)
     dedb0=cinv*((upara*upara-upara*upara0)*aion*b_inv+zion(6,m)*zion(6,m))
     deni=1.0/(g*q + ri + (zion(4,m)+lam)*(g*rip-ri*gp))
     gqi=1.0/(g*q+ri)
     energy=0.5*aion*upara*upara+zion(6,m)*zion(6,m)*b
     energy0=0.5*aion*(upara-upara0)*(upara-upara0)+zion(6,m)*zion(6,m)*b


     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     kappan=dp1*kapani(ii-1)+(1.0-dp1)*kapani(ii)
     kappat=dp1*kapati(ii-1)+(1.0-dp1)*kapati(ii)
     ti_inv=1.0/(dp1*meshti(ii-1)+(1.0-dp1)*meshti(ii))
     kappa=kappan+kappat*(energy0*ti_inv-1.5)
     kappav=-(upara-upara0)*dupara0dp*aion*ti_inv
     ne=dp1*meshne(ii-1)+(1.0-dp1)*meshne(ii)
     te=dp1*meshte(ii-1)+(1.0-dp1)*meshte(ii)
     kne=dp1*kapane(ii-1)+(1.0-dp1)*kapane(ii)

! parallel electric field epara, & dlambda/dt=-epara-gradphi(3,m)
     if (magnetic==0) then
        epara=-wpgc(3,m)*b*q*gqi
        plampt=0.0
     else
        epara=-wbgc(4,m)*b*q*gqi
        plampt=-(epara+wpgc(3,m)*q*b*gqi)*b_inv
     endif

! ExB drift in radial direction for w-dot and flux diagnostics
     vdr=(ri*dptdz-g*dptdt)*gqi
     vap=-(ri*dapdz-g*dapdt)*(upara-upara0)*gqi

     wdrive=(kappa+kappav)*(vdr+vap)
     wpara=epara*(upara-upara0)*qion*ti_inv
     wdrift=(g*dbdt*dptdp-g*dbdp*dptdt+ri*dbdp*dptdz)*gqi*dedb0*qion*ti_inv
     
     
     if(eqcurrent==1)then
       wdrift=wdrift+aion*(upara*upara-upara*upara0)*ti_inv*gqi*(gp*dptdt-rip*dptdz)
     endif

! self-consistent and external electric field for marker orbits
     dptdp=dptdp*perturb - er
     dptdt=dptdt*perturb
     dptdz=dptdz*perturb
     dlamdp=dlamdp*perturb
     dlamdt=dlamdt*perturb
     dlamdz=dlamdz*perturb
     lam=lam*perturb
     plampt=plampt*perturb

! particle velocity
      uperp=zion(6,m)*sqrt(b)*sqrt(2.0/real(aion))
      angle=zion(7,m)-zion(2,m)
      rdot=b_inv*dptdz
      zdot=upara
      pdot=uperp*cos(angle)/drgdp 
      tdot=uperp*sin(angle)/rg
      
! update particle position
     if(zion0(1,m) < paxis)then
! particles close to axis use (x,y) coordinates
       pdum=sqrt(zion0(1,m))
	   cost=cos(zion0(2,m))
       sint=sin(zion0(2,m))
       xdum   = pdum*cost
       ydum   = pdum*sint
       pdum=1.0/zion(1,m)
       xdot   = 0.5*pdot*xdum*pdum-ydum*tdot
       ydot   = 0.5*pdot*ydum*pdum+xdum*tdot
       pdum=sqrt(zion0(1,m))
       xdum   = pdum*cos(zion0(2,m)) + dtime*xdot
       ydum   = pdum*sin(zion0(2,m)) + dtime*ydot
       zion(1,m) = max(1.0e-8_lk*psi1,xdum*xdum+ydum*ydum)
       zion(2,m) = sign(1.0_lk,ydum)*acos(max(-1.0_lk,min(1.0_lk,xdum/sqrt(zion(1,m)))))
     else
       zion(1,m) = max(1.0e-8_lk*psi1,zion0(1,m)+dtime*pdot)
       zion(2,m) = zion0(2,m)+dtime*tdot
     endif

     zion(3,m) = zion0(3,m)+dtime*zdot
     zion(4,m) = zion0(4,m)+dtime*rdot

!Boris push, formulation from Buneman 1973
! First step
        Ex=-(wpgc(1,m)/drgdp*cos(zion(2,m))-wpgc(2,m)/rg*sin(zion(2,m)))
        Ey=-(wpgc(1,m)/drgdp*sin(zion(2,m))+wpgc(2,m)/rg*cos(zion(2,m)))
        delvx=cmratio*dtime*0.5*Ex*perturb
        delvy=cmratio*dtime*0.5*Ey*perturb
        uperp0=zion(7,m)*sqrt(b)*sqrt(2.0/real(aion))

!seperate the rotation from the acceleration from E field.
        vxtemp=uperp0*cos(zion0(7,m))+delvx
        vytemp=uperp0*sin(zion0(7,m))+delvy
        vx1=vxtemp
        vy1=vytemp

! Second step
!the rotated angle        
       tandalpha=cmratio*b*dtime*0.5
       sindalpha=2*tandalpha/(1+tandalpha*tandalpha)

       utemp=vxtemp+vytemp*tandalpha
       vytemp=vytemp-utemp*sindalpha
       vxtemp=utemp+vytemp*tandalpha

! Third step
!transform back to include E field
       vxtemp=vxtemp+delvx
       vytemp=vytemp+delvy
       vx2=vxtemp
       vy2=vytemp   
    
!get new zion(6,m) and zion(7,m)
        zion(6,m)=sqrt((vxtemp*vxtemp+vytemp*vytemp)/b)*sqrt(real(aion)/2.0)
       if (abs(vxtemp)/(abs(vytemp)+abs(vxtemp))<1e-6) then
           if(vytemp>0)then
                zion(7,m)=pi/2 
           else 
                zion(7,m)=-pi/2
           endif
       else
           zion(7,m)=atan(vytemp/vxtemp)
           if(vxtemp<0) zion(7,m)=zion(7,m)+pi
       endif

! Electric field and Velocity in X-Y coordinate
    
     wdot=0.5*(Ex*(vx1+vx2)+Ey*(vy1+vy2))*qion*ti_inv!+epara*upara*qion*ti_inv  ! For uniform plasma

     zion(5,m)=zion0(5,m)+wdot*dtime


     zion(2,m)=modulo(zion(2,m),pi2)
     zion(3,m)=modulo(zion(3,m),pi2)

! store GC information for flux measurements
     wpgc(1,m)=vdr
     wpgc(2,m)=energy
     wpgc(3,m)=upara*majorr
  enddo

  if(iload>99)zion(5,:)=1.0 ! for full-f simulation

  if(irk==2)then

! out of boundary particle
!$omp parallel do private(m)
     do m=1,mi
        if(zion(1,m) > psi1)then
              zion(7,m)=modulo(zion(2,m)*2.0-zion(7,m)+pi,pi2)
              zion(1,m)=2.0*psi1-zion(1,m)
              zion(2,m)=zion0(2,m)
              zion(3,m)=zion0(3,m)
              zion(4,m)=zion0(4,m)
              zion(5,m)=zion0(5,m)         

        elseif(zion(1,m) < psi0)then
              zion(7,m)=modulo(zion(2,m)*2.0-zion(7,m)+pi,pi2)
              zion(1,m)=2.0*psi0-zion(1,m)
              zion(2,m)=zion0(2,m)
              zion(3,m)=zion0(3,m)
              zion(4,m)=zion0(4,m)
              zion(5,m)=zion0(5,m)
        endif
     enddo

! Restore temperature profile for nonlinear delta-f simulation
     if(nonlinear==1 .and. iload==1)then

!$omp parallel do private(m,psitmp,isp,dpx,rg)
        do m=1,mi
           psitmp=zion(1,m)
           isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
           dpx=psitmp-spdpsi*real(isp-1)
           if(isp==1)dpx=sqrt(dpx)
           rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dpx*dpx
           iir(m)=max(0,min(mpsi,int((rg-rg0)*delr+0.5)))
        enddo

        dtem=0.0
        dden=0.0
        dmark=0.0
! Following loop is a gathering operation: should not use OpenMP
        do m=1,mi
           ii=iir(m)
           fullf=zion(7,m)
           deltaf=fullf*zion(5,m)
           dtem(ii)=dtem(ii)+wpgc(2,m)*deltaf
           dmark(ii)=dmark(ii)+wpgc(1,m)*fullf
           dden(ii)=dden(ii)+fullf
        enddo

        icount=mpsi+1
        ddum=0.0
        call MPI_ALLREDUCE(dtem,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dtem=ddum
        call MPI_ALLREDUCE(dmark,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dmark=ddum
        call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
        dden=ddum

        irsmooth=int(sqrt(real(mpsi)))
        dmark=dmark/max(1.0_lk,dden) !radial marker flux
        do i=1,irsmooth
           rdum=dmark(1)
           tdum=dmark(mpsi-1)
           dmark(1:mpsi-1)=0.5*dmark(1:mpsi-1)+0.25*(dmark(0:mpsi-2)+dmark(2:mpsi))
           dmark(0)=0.5*(dmark(0)+rdum)
           dmark(mpsi)=0.5*(dmark(mpsi)+tdum)
        enddo
        tdum=0.1
        pfluxi=(1.0-tdum)*pfluxi+tdum*dmark

! remove small scale temperature perturbation
        irsmooth=mpsi
        dtem=dtem/(meshti*max(1.0_lk,dden))
        do i=1,irsmooth
           rdum=dtem(1)
           tdum=dtem(mpsi-1)
           dtem(1:mpsi-1)=0.5*dtem(1:mpsi-1)+0.25*(dtem(0:mpsi-2)+dtem(2:mpsi))
           dtem(0)=0.5*(dtem(0)+rdum)
           dtem(mpsi)=0.5*(dtem(mpsi)+tdum)
        enddo
        tdum=0.01
        rdtemi=(1.0-tdum)*rdtemi+tdum*dtem

!c$omp parallel do private(m,ii)
!        do m=1,mi
!           ii=iir(m)
!           zion(5,m)=zion(5,m)-(wpgc(2,m)/meshti(ii)-1.5)*rdtemi(ii)
!        enddo
     endif
  endif

  if(idiag==0)then
! fluxes diagnose at irk=1
     diagion=0.0_lk
     dden=0.0_lk
     data1di=0.0_lk
     do m=1,mi
        psitmp=zion(1,m)
        isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
        dpx=psitmp-spdpsi*real(isp-1)
        if(isp==1)dpx=sqrt(dpx)
        dp2=dpx*dpx
        rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
       !! r= rpsi(1,isp)  +rpsi(2,isp)*dpx  +rpsi(3,isp)*dp2
       !! q=    qpsi(1,isp)   +qpsi(2,isp)*dpx   +qpsi(3,isp)*dp2

        ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
        dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
        ti_inv=1.0/(dp1*meshti(ii-1)+(1.0-dp1)*meshti(ii))

        fullf=zion(7,m)
        deltaf=fullf*zion0(5,m)
        energy=wpgc(2,m)*ti_inv-1.5
        vdr=wpgc(1,m)!!!*q/r
        angmom=wpgc(3,m)

! radial profile of particle and energy flux
        dden(ii-1)=dden(ii-1)+fullf*dp1
        dden(ii)  =dden(ii)+  fullf*(1.0-dp1)

        data1di(ii-1,1)=data1di(ii-1,1)+vdr*deltaf*dp1
        data1di(ii,  1)=data1di(ii,  1)+vdr*deltaf*(1.0-dp1)

        data1di(ii-1,2)=data1di(ii-1,2)+vdr*deltaf*energy*dp1
        data1di(ii,  2)=data1di(ii,  2)+vdr*deltaf*energy*(1.0-dp1)
! radial profiles of momentum flux
        data1di(ii-1,3)=data1di(ii-1,3)+vdr*deltaf*angmom*dp1
        data1di(ii,  3)=data1di(ii,  3)+vdr*deltaf*angmom*(1.0-dp1)

!!! ion diagnosis: density,entropy,flow,energy,fluxes of particle,momentum,heat
        diagion(1)=diagion(1)+deltaf
        diagion(2)=diagion(2)+deltaf*deltaf
        diagion(3)=diagion(3)+angmom
        diagion(4)=diagion(4)+angmom*deltaf
        diagion(5)=diagion(5)+energy
        diagion(6)=diagion(6)+energy*deltaf
        diagion(7)=diagion(7)+vdr*deltaf
        diagion(8)=diagion(8)+vdr*angmom*deltaf
        diagion(9)=diagion(9)+vdr*energy*deltaf
     enddo
     diagion(10)=real(mi)

! sum over all MPI processes
     icount=mpdiag
     diagtmp=0.0
     call MPI_REDUCE(diagion,diagtmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
     diagion=diagtmp

     icount=(mpsi+1)*mpdata1d
     data1dtmp=0.0
     call MPI_ALLREDUCE(data1di,data1dtmp,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

     icount=mpsi+1
     call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)

! radial profile data normalized by marker #
     do i=1,mpdata1d
        data1di(:,i)=data1dtmp(:,i)/max(1.0,ddum)
     enddo
  endif

end subroutine pushifk
