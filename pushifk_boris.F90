subroutine pushifk_boris
use global_parameters
use particle_array
use field_array
use equilibrium
use magnetic_island,only: gradalphaIs, alphaIs!, tornum,polnum ! tornum & polnum are obsolete
use spline_function, only: spx, dxdp, dxdt, dzdp, dzdt, dots, cross
implicit none


!local var
integer i,m,ii,ierror,ij,irsmooth,icount,iir(mimax),isp,jst,igyro
real(lk) rdum,tdum,zdum,dtime,q,q_inv,rg,b,dbdt,dbdp,dbdz,g,ri,gqi,deltaf,upara,energy,&
         dptdp,dptdt,dptdz,vdr,pdot,tdot,zdot,rdot,wdot,delr,xdum,ydum,xdot,ydot,&
         wz1,wz0,wp1,wp0,wt11,wt10,wt01,wt00,tp_inv,epara,wdrive,wdrift,wpara,wisland,&
         perturb,paxis,pdum,psitmp,thetatmp,zetatmp,e1,e2,e3,cmratio,&
         cinv,wpgc(3,mimax),dden(0:mpsi),ddum(0:mpsi),gp,rip,dedb0,model_fac,&
         diagtmp(mpdiag),data1dtmp(0:mpsi,mpdata1d),dpx,dp2,dtx,dt2,dzx,dx(27),&
         b1,b2,b3,b4,b5,wbgc(5,mimax),dapdp,dapdt,dapdz,dapidt,dapidz,&
         lam,dlamdp,dlamdt,dlamdz,plampt,delp(mpsi),dp1,fullf,energy0,term2,term3,&
         b_inv,majorr,dmajorrdp,angmom,er,upara0,kappa,kappan,kappat,kappav,ti_inv,vap,&
         i1,i2,i3,i4,wigc(4,mimax),te,ne,kne,cost,sint,gpp,gtt,gzz,dupara0dp,&
         gpt,gpz,gtz,gdelta,invgpp,invgtt,invgzz,invgpt,invgpz,invgtz,epsi0(3),etheta0(3),&
         ezeta0(3),delta0,delta1,delta2,invep(3),invet(3),invez(3),b_vec(3),&
         v_vec(3),temp_vec0(3),temp_vec1(3),temp_vec2(3),gammab,jacobian,upsi0,&
         utheta0,uzeta0,upsi1,utheta1,uzeta1,upsi2,utheta2,uzeta2,upsi3,utheta3,&
         uzeta3,upsi4,utheta4,uzeta4,dxdptmp,dxdttmp,dzdptmp,dzdttmp

paxis=0.0_lk
if(psi0<1.0e-8)paxis=psimesh(1)
delr=1.0_lk/deltar
cmratio=qion/aion
cinv=1.0_lk/qion
perturb=real(nonlinear)
delp=1.0_lk/deltap


!do not use Runge-Kutta method
dtime=0.5_lk*tstep
!$omp parallel do private(m)
do m=1,mi
  zion0(1:4,m)=zion(1:4,m)
  zion0(6,m)=zion(6,m)
  zion0(8,m)=zion(8,m)
enddo

if(irk==1)then
  !$omp parallel do private(m)
  do m=1,mi
    zion0(5,m)=zion(5,m)
  enddo
endif

  !magnetic island perturbation
  if(island==1)then
    !$omp parallel do private(m,igyro,i1,i2,i3,i4,wz1,wz0,wp1,wp0,wt00,wt10,wt01,wt11,ij)
    do m=1,mi
      i1=0.0  
      i2=0.0
      i3=0.0
      i4=0.0
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
        i1=i1+wt00*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt00*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt00*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt00*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))


        ij=ij+1                  !upper poloidal grid on inner flux surface
        i1=i1+wt10*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt10*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt10*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt10*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))

        ij=jtion1(igyro,m)       !lower poloidal grid on outer flux surface
        i1=i1+wt01*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt01*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt01*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt01*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))

        ij=ij+1                  !upper poloidal grid on outer flux surface
        i1=i1+wt11*(wz0*gradalphaIs(1,0,ij)+wz1*gradalphaIs(1,1,ij))
        i2=i2+wt11*(wz0*gradalphaIs(2,0,ij)+wz1*gradalphaIs(2,1,ij))
        i3=i3+wt11*(wz0*gradalphaIs(3,0,ij)+wz1*gradalphaIs(3,1,ij))
        i4=i4+wt11*(wz0*alphaIs(0,ij)+wz1*alphaIs(1,ij))
      enddo

      wigc(1,m)=i1
      wigc(2,m)=i2
      wigc(3,m)=i3
      wigc(4,m)=i4
    enddo
  endif
  

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

      wpgc(1,m)=e1
      wpgc(2,m)=e2
      wpgc(3,m)=e3
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

      wpgc(1,m)=e1
      wpgc(2,m)=e2
      wpgc(3,m)=e3
      wbgc(1,m)=b1
      wbgc(2,m)=b2
      wbgc(3,m)=b3
      wbgc(4,m)=b4
      wbgc(5,m)=b5
    enddo
  endif

!update particle position. Assuming psi0>spdpsi; will be relaxed later
!$omp parallel do private(m,psitmp,thetatmp,isp,dpx,dp2,jst,dtx,dt2,dzx,dx,rg,q,g,gp,rip,&
!$omp& ri,b,dbdp,dbdt,dbdz,ii,gqi,upara,energy,dptdp,dptdt,zetatmp,wdrive,wpara,wdrift,&
!$omp& dptdz,dapdp,dapdt,dapdz,epara,vdr,vap,wdot,dapidt,dapidz,upara0,wisland,term2,term3,&
!$omp& lam,dlamdp,dlamdt,dlamdz,plampt,dp1,kappan,kappat,ti_inv,b_inv,q_inv,majorr,er,gpp,gtt,gzz,&
!$omp& gpt,gpz,gtz,gdelta,invgpp,invgtt,invgzz,invgpt,invgpz,invgtz,epsi0,etheta0,kappa,kappav,&
!$omp& ezeta0,delta0,delta1,delta2,invep,invet,invez,b_vec,v_vec,dedb0,dupara0dp,model_fac,&
!$omp& temp_vec0,temp_vec1,temp_vec2,gammab,jacobian,upsi0,utheta0,uzeta0,dmajorrdp,&
!$omp& upsi1,utheta1,uzeta1,upsi2,utheta2,uzeta2,upsi3,utheta3,uzeta3,upsi4,energy0,&
!$omp& utheta4,uzeta4,dxdptmp,dxdttmp,dzdptmp,dzdttmp)
do m=1,mi
  psitmp=zion(1,m)
  isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
  dpx=psitmp-spdpsi*real(isp-1)
  dp2=dpx*dpx
  
  !when fieldmodel==1?
  thetatmp=zion(2,m)
  jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
  dtx=thetatmp-spdtheta*real(jst-1)
  dt2=dtx*dtx
  dzx=zion(3,m)-zeta0

  dx(1)=1.0
  dx(2)=dpx
  dx(3)=dp2
  dx(4:6)=dx(1:3)*dtx
  dx(7:9)=dx(1:3)*dt2
  dx(10:18)=dx(1:9)*dzx
  dx(19:27)=dx(1:9)*dzx*dzx
  
  zetatmp=zion(3,m)
  !1D spline in psi
  q=    qpsi(1,isp)   +qpsi(2,isp)*dpx    +qpsi(3,isp)*dp2
  g=    gpsi(1,isp)   +gpsi(2,isp)*dpx    +gpsi(3,isp)*dp2
  ri=   cpsi(1,isp)   +cpsi(2,isp)*dpx    +cpsi(3,isp)*dp2
  gp=                  gpsi(2,isp)       +gpsi(3,isp)*dpx*2.0
  rip=                 cpsi(2,isp)       +cpsi(3,isp)*dpx*2.0
  upara0=ropp(1,isp)  +ropp(2,isp)*dpx   +ropp(3,isp)*dp2
  dupara0dp=           ropp(2,isp)       +ropp(3,isp)*dpx*2.0
  er=   erpp(1,isp)   +erpp(2,isp)*dpx    +erpp(3,isp)*dp2
  majorr=spx(psitmp,thetatmp)
  dmajorrdp=dxdp(psitmp,thetatmp)
  
  model_fac=1.0
  if(fieldmodel==0 .and. numereq<100) model_fac=majorr
  ri=ri/model_fac
  
  ! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
  if(isp==1)dpx=sqrt(dpx)
  dp2=dpx*dpx
  rg= rgpsi(1,isp)    +rgpsi(2,isp)*dpx   +rgpsi(3,isp)*dp2
  
  !spline in(psi,theta,zeta)
  b=0.0
  do ii=1,spdim
    b=b+bsp(ii,isp,jst)*dx(ii)
  enddo
! zeta-derivative
  dx(10:18)=dx(1:9)*1.0
  dx(19:27)=dx(1:9)*dzx*2.0
  dbdz=0.0
  do ii = 10, spdim
    dbdz=dbdz +bsp(ii,isp,jst)*dx(ii)
  enddo
! theta-derivative
  dx(4:6)=dx(1:3)*1.0
  dx(7:9)=dx(1:3)*dtx*2.0
  dx(10:18)=dx(1:9)*dzx
  dx(19:27)=dx(1:9)*dzx*dzx
  dbdt=0.0
  do ii = 4, spdim
    dbdt=dbdt +bsp(ii,isp,jst)*dx(ii)
  enddo
  
! psi-derivative
  dx(1)=0.0
  dx(2)=1.0
  if(isp==1)dx(2)=0.5/dpx
  dx(3)=dpx*2.0
  if(isp==1)dx(3)=1.0
  dx(4:6)=dx(1:3)*dtx
  dx(7:9)=dx(1:3)*dt2
  dx(10:18)=dx(1:9)*dzx
  dx(19:27)=dx(1:9)*dzx*dzx
  dbdp=0.0
  do ii = 2, spdim
    dbdp=dbdp +bsp(ii,isp,jst)*dx(ii)
    dmajorrdp=dmajorrdp +xsp(ii,isp,jst)*dx(ii)
  enddo
  

  b_inv=1.0/b

#ifdef _FRC
  q_inv=0.0
#elif _CYLINDER
  q_inv=0.0
#else
  q_inv=1.0/q
#endif

  !perturbed electric field
  dptdp=wpgc(1,m)
  dptdt=wpgc(2,m)
  dptdz=wpgc(3,m)-wpgc(2,m)*q_inv
  
  !perturbed vector potential
  dapdp=wbgc(1,m)
  dapdt=wbgc(2,m)
  dapdt=wbgc(3,m)-wbgc(2,m)*q_inv
  
  !island perturbed magnetic field
  if(island==1)then
    dapidt=wigc(2,m)
    dapidz=wigc(3,m)
  endif
  !lambda=apara*b_inv
  lam=wbgc(5,m)*b_inv
  dlamdp=(dapdp-lam*dbdp)*b_inv
  dlamdt=(dapdt-lam*dbdt)*b_inv
  dlamdz=(dapdz-lam*dbdz)*b_inv
  
  upara=zion(4,m)*b*cmratio
  upara0=upara0*majorr
  dupara0dp=dupara0dp*majorr+upara0*dmajorrdp
  gqi=1.0/(g*q+ri)
  dedb0=cinv*((upara*upara-upara*upara0)*aion*b_inv+zion(6,m)*zion(6,m))
  energy=0.5*aion*upara+zion(6,m)**2*b
  energy0=0.5*aion*(upara-upara0)*(upara-upara0)+zion(6,m)*zion(6,m)*b
  
  ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
  dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
  kappan=dp1*kapani(ii-1)+(1.0-dp1)*kapani(ii)
  kappat=dp1*kapati(ii-1)+(1.0-dp1)*kapati(ii)
  ti_inv=1.0/(dp1*meshti(ii-1)+(1.0-dp1)*meshti(ii))
  kappa=kappan+kappat*(energy0*ti_inv-1.5)
  kappav=-(upara-upara0)*dupara0dp*aion*ti_inv
  
  !parallel electric field epara, & dlambda/dt=-epara-gradphi(3,m)
  if(magnetic==0)then
    epara=-(dptdt+q*dptdz)*b*gqi
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
  wdrift=(g*(dbdt*dptdp-dbdp*dptdt)+ri*(dbdp*dptdz-dbdz*dptdp))*gqi*dedb0*qion*ti_inv
  
  if(island==1)then
    wisland=-kappa*(ri*b*dapidz-g*b*dapidt)*(upara-upara0)*gqi
  endif
    
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
  
  ! island modification on lambda
    if(island==1)then
      lam=lam+wigc(4,m)
      dlamdp=dlamdp
      dlamdt=dlamdt+dapidt
      dlamdz=dlamdz+dapidz
    endif
    
  !particle velocity
  dxdptmp=dxdp(zion0(1,m),zion0(2,m))
  dxdttmp=dxdt(zion0(1,m),zion0(2,m))
  dzdptmp=dzdp(zion0(1,m),zion0(2,m))
  dzdttmp=dzdt(zion0(1,m),zion0(2,m))
  gpp=dxdptmp**2+dzdptmp**2
  gtt=dxdttmp**2+dzdttmp**2
  gzz=majorr**2
  gpt=dxdptmp*dxdttmp+dzdttmp*dzdptmp
  gpz=0
  gtz=0
  gdelta=gpp*gtt-gpt**2
  invgpp=gtt/gdelta
  invgtt=gpp/gdelta
  invgpt=-gpt/gdelta
  invgzz=1/gzz
  invgpz=0
  invgtz=0
!  covariant and contravariant variables in the first place
  epsi0 = (/sin(zetatmp)*dxdptmp,cos(zetatmp)*dxdptmp,dzdptmp/)
  etheta0=(/sin(zetatmp)*dxdttmp,cos(zetatmp)*dxdttmp,dzdttmp/)
  ezeta0= (/majorr*cos(zetatmp),-majorr*sin(zetatmp),0.0_lk/)
  delta0=dxdptmp*dzdttmp-dxdttmp*dzdptmp
  invep = (/dzdttmp*sin(zetatmp)/delta0,dzdttmp*cos(zetatmp)/delta0,&
            -dxdttmp/delta0/)
  invet = (/-dzdptmp*sin(zetatmp)/delta0,-dzdptmp*cos(zetatmp)/delta0,&
            dxdptmp/delta0/)
  invez = (/cos(zetatmp)/majorr,-sin(zetatmp)/majorr,0.0_lk/)
 
  b_vec=ri*invet+g*invez
  b=sqrt(dots(b_vec,b_vec,3))
  b_inv=1.0_lk/b
  gammab=1/(1+(b*cmratio*dtime)**2/4)
  jacobian=majorr*delta0

  upsi0=zion(4,m)
  utheta0=zion(6,m)
  uzeta0=zion(8,m)

  upsi1 = upsi0  + 0.5_lk*cmratio*dtime*(dptdp*invgpp+dptdt*invgpt)
  utheta1=utheta0+0.5_lk*cmratio*dtime*(dptdp*invgpt+dptdt*invgtt)
  uzeta1= uzeta0 + 0.5_lk*cmratio*dtime*dptdz*invgzz

  
  upsi2 = (1-0.5*gammab*(b*cmratio*dtime)**2)*upsi1+gammab/jacobian*&
             cmratio*dtime*(upsi1*g*gpt+utheta1*g*gtt-ri*uzeta1*gzz)
  utheta2=(1-0.5*gammab*(b*cmratio*dtime)**2)*utheta1+gammab/jacobian*&
             cmratio*dtime*(-upsi1*g*gpp-utheta1*g*gpt)+gammab*&
             (cmratio*dtime)**2/2*(ri*utheta1+g*uzeta1)*ri*invgtt
  uzeta2= (1-0.5*gammab*(b*cmratio*dtime)**2)*uzeta1+gammab/jacobian*&
             cmratio*dtime*(upsi1*ri*gpp+utheta1*ri*gpt)+gammab*&
            (cmratio*dtime)**2/2*(ri*utheta1+g*uzeta1)*g*invgzz
  
  upsi3 = upsi2 + 0.5*cmratio*dtime*(dptdp*invgpp+dptdt*invgpt)
  utheta3=utheta2+0.5*cmratio*dtime*(dptdp*invgpt+dptdt*invgtt)
  uzeta3= uzeta2+ 0.5*cmratio*dtime*dptdz*invgzz
  
  
  
  zion(1,m)=zion0(1,m)+upsi3*dtime/2
  zion(2,m)=zion0(2,m)+utheta3*dtime/2
  zion(3,m)=zion0(3,m)+uzeta3*dtime/2
  zion(2,m)=modulo(zion(2,m),pi2)
  zion(3,m)=modulo(zion(3,m),pi2)
  if(zion(1,m)<psi0 .or. zion(1,m)>psi1)then
    zion(1,m)=zion0(1,m)
    zion(2,m)=2*pi-zion0(2,m)
    zion(3,m)=zion0(3,m)
  endif
  
  dxdptmp=dxdp(zion(1,m),zion(2,m))
  dxdttmp=dxdt(zion(1,m),zion(2,m))
  dzdptmp=dzdp(zion(1,m),zion(2,m))
  dzdttmp=dzdt(zion(1,m),zion(2,m))
  zetatmp=zion(3,m)
  majorr=spx(zion(1,m),zion(2,m))
!  contravariant in the middle(approximate) point
  delta1=dxdptmp*dzdttmp-dxdttmp*dzdptmp
  invep=(/dzdttmp*sin(zetatmp)/delta1,dzdttmp*cos(zetatmp)/delta1,&
         -dxdttmp/delta1/)
  invet=(/-dzdptmp*sin(zetatmp)/delta1,-dzdptmp*cos(zetatmp)/delta1,&
         dxdptmp/delta1/)
  invez=(/cos(zetatmp)/majorr,-sin(zetatmp)/majorr,0.0_lk/)
  
  upsi4 = upsi3*dots(epsi0,invep,3)+utheta3*dots(etheta0,invep,3)+uzeta3*dots(ezeta0,invep,3)
  utheta4=upsi3*dots(epsi0,invet,3)+utheta3*dots(etheta0,invet,3)+uzeta3*dots(ezeta0,invet,3)
  uzeta4= upsi3*dots(epsi0,invez,3)+utheta3*dots(etheta0,invez,3)+uzeta3*dots(ezeta0,invez,3)

!for near-axis particles
  if(zion0(1,m)<paxis)then
    pdum=sqrt(zion0(1,m))
    cost=cos(zion0(2,m))
    sint=sin(zion0(2,m))
    xdum=pdum*cost
    ydum=pdum*sint
    pdum=1.0/zion0(1,m)
    xdot= 0.5*upsi4*xdum*pdum-ydum*utheta4
    ydot= 0.5*upsi4*ydum*pdum+xdum*utheta4
    pdum=sqrt(zion0(1,m))
    xdum= pdum*cos(zion0(2,m))+dtime*xdot
    ydum= pdum*sin(zion0(2,m))+dtime*ydot
    zion(1,m)=max(1.0e-8_lk*psi1,xdum*xdum+ydum*ydum)
    zion(2,m)=sign(1.0_lk,ydum)*acos(max(-1.0_lk,min(1.0_lk,xdum/sqrt(zion(1,m)))))
  else
    zion(1,m)=max(1.0e-8_lk*psi1,zion0(1,m)+upsi4*dtime)
    zion(2,m)=zion0(2,m)+utheta4*dtime
  endif
  zion(3,m)=zion0(3,m)+uzeta4*dtime
  zion(2,m)=modulo(zion(2,m),pi2)
  zion(3,m)=modulo(zion(3,m),pi2)
  if(zion(1,m)<psi0 .or. zion(1,m)>psi1)then
    zion(1,m)=zion0(1,m)
    zion(2,m)=2*pi-zion0(2,m)
    zion(3,m)=zion0(3,m)
  endif
  
  zetatmp=zion(3,m)
  majorr =spx(zion(1,m),zion(2,m))
  dxdptmp=dxdp(zion(1,m),zion(2,m))
  dxdttmp=dxdt(zion(1,m),zion(2,m))
  dzdptmp=dzdp(zion(1,m),zion(2,m))
  dzdttmp=dzdt(zion(1,m),zion(2,m))
  !  contravariant in the last point
  delta2=dxdptmp*dzdttmp-dxdttmp*dzdptmp
  invep=(/dzdttmp*sin(zetatmp)/delta2,dzdttmp*cos(zetatmp)/delta2,&
         -dxdttmp/delta2/)
  invet=(/-dzdptmp*sin(zetatmp)/delta2,-dzdptmp*cos(zetatmp)/delta2,&
         dxdptmp/delta2/)
  invez=(/cos(zetatmp)/majorr,-sin(zetatmp)/majorr,0.0_lk/)
  
  zion(4,m)=upsi3*dots(epsi0,invep,3)+utheta3*dots(etheta0,invep,3)+uzeta3*dots(ezeta0,invep,3)
  zion(6,m)=upsi3*dots(epsi0,invet,3)+utheta3*dots(etheta0,invet,3)+uzeta3*dots(ezeta0,invet,3)
  zion(8,m)=upsi3*dots(epsi0,invez,3)+utheta3*dots(etheta0,invez,3)+uzeta3*dots(ezeta0,invez,3) 

  CALL CROSS(v_vec,b_vec,temp_vec0)

  if(irk==2)then
    !from Kuley's paper, 2013
    wdot=-0.5*qion*ti_inv*((upsi0+upsi3)*dptdp+(utheta0+utheta3)*dptdt+(uzeta0+uzeta3)*dptdz)
    zion(5,m)=zion0(5,m)+wdot*dtime*2
  endif

  energy=0.5*aion*(upsi0**2*gpp+utheta0**2*gtt+uzeta0**2*gzz+2*upsi0*utheta0*gpt)

!  energyion(m)=energy
  wpgc(1,m)=vdr
  wpgc(2,m)=energy
  wpgc(3,m)=upara*majorr
enddo
if(iload>99)zion(5,:)=1.0 !for full-f simulation

if(idiag==0)then
  !fluxes diagnose at irk==1
  diagion=0.0_lk
  dden=0.0_lk
  data1di=0.0_lk
  do m=1,mi
    psitmp=zion(1,m)
    isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
    dpx=psitmp-spdpsi*real(isp-1)
    if (isp==1)dpx=sqrt(dpx)
    dp2=dpx*dpx
    rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
    ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
    dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
    tp_inv=1.0/(dp1*meshti(ii-1)+(1.0-dp1)*meshti(ii))
    fullf=zion(7,m)
    deltaf=fullf*zion0(5,m)
    energy=wpgc(2,m)*tp_inv-1.5
    vdr=wpgc(1,m)
    angmom=wpgc(3,m)
    
    !radial profile of particle and energy flux
    dden(ii-1)=dden(ii-1)+fullf*dp1
    dden(ii)  =dden(ii)+  fullf*(1.0-dp1)
    
    data1di(ii-1,1)=data1di(ii-1,1)+vdr*deltaf*dp1
    data1di(ii,  1)=data1di(ii,  1)+vdr*deltaf*(1.0-dp1)
    
    data1di(ii-1,2)=data1di(ii-1,2)+vdr*deltaf*energy*dp1
    data1di(ii,  2)=data1di(ii,  2)+vdr*deltaf*energy*(1.0-dp1)
    ! radial profiles of momentum flux
    data1di(ii-1,3)=data1di(ii-1,3)+vdr*deltaf*angmom*dp1
    data1di(ii,  3)=data1di(ii,  3)+vdr*deltaf*angmom*(1.0-dp1)
    
    !!! diagnosis: density,entropy,flow,energy,fluxes of particle,momentum,heat
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
  
  !radial profile data normalized by marker #
  do i=1,mpdata1d
    data1di(:,i)=data1dtmp(:,i)/max(1.0,ddum)
  enddo
endif
end subroutine pushifk_boris

