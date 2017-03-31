subroutine pushfield
  use global_parameters
  use field_array
  use particle_array
  use equilibrium
  use precision
  implicit none

  integer i,j,ij,isp,jst,ierror,iMHDxi_output
  real(lk) dtime,pdum,tdum,dpx,dp2,dtx,dt2,dbdt,dbdp,dbdz,&
           ri,g,b_inv,q_inv,&
           dptdp,dptdt,dptdz,wpara,wdrive,wstar,wdrift,wdamp,&
           ne,te,dnedp,dtedp,phipEq,wdrift0,web0,ddnedp,ddnedt,ddnedz,&
           ddpdp,ddpdt,ddpdz,dinddp,dinddt,dinddz,dpsidp,dpsidt,dpsidz,dne,&
           daparadp,daparadt,daparadz,dsdp,dphieffdp,dphieffdt,dphieffdz,&
           dkpparadp,dkpparadt,dkpparadz,dkpperpdp,dkpperpdt,dkpperpdz,&
           dqdp,dgdp,d2gdp2,dridp,d2ridp2,weqc,rbc(0:mpsi),&
           dudp,dudt,dudz,wstarnl,wperpnl,perturb,apara00tmp(0:mpsi),&
           dBdotdE00,JdBdotdE,adum(0:mpsi),gqi_inv

  perturb=real(nonlinear)
  !radial BC
  rbc(0:mpsi)=1.0
!  do i=0,nbound-1
!     rbc(i)=real(i)/real(nbound)
!  enddo

  if(irk==1)then

! 1st step of Runge-Kutta method
    dtime=0.5*tstep

!$omp parallel do private(ij)
    do ij=1,mgrid
        apara0(1,ij)=apara(1,ij)
        fluidne0(1,ij)=fluidne(1,ij)
        deltapsi0(1,ij)=deltapsi(1,ij)
    enddo

    apara00nl0=apara00nl

! 2nd step of Runge-Kutta method
  else
    dtime=tstep
  endif

! advance del_ne and a_parallel
!$omp parallel do private(i,j,ij,pdum,isp,dpx,dp2,g,ri,tdum,jst,dtx,dt2,&
!$omp& b_inv,q_inv,phipEq,dbdt,dbdp,dbdz,ne,te,dnedp,dtedp,&
!$omp& wpara,wstar,wdrift,wdrive,wdrift0,web0,dne,wdamp,&
!$omp& ddnedp,ddnedt,ddnedz,ddpdp,ddpdt,ddpdz,dptdp,dptdt,dptdz,dinddp,dinddt,dinddz,&
!$omp& daparadp,daparadt,daparadz,dpsidp,dpsidt,dpsidz,dphieffdp,dphieffdt,dphieffdz,&
!$omp& dkpparadp,dkpparadt,dkpparadz,dkpperpdp,dkpperpdt,dkpperpdz,&
!$omp& dqdp,dgdp,dridp,d2gdp2,d2ridp2,dsdp,weqc,dudp,dudt,dudz,wstarnl,wperpnl,&
!$omp& JdBdotdE,dBdotdE00,gqi_inv)
  do i=0,mpsi
     pdum=psimesh(i)
     q_inv=1/qmesh(i)
     dBdotdE00=0.0

     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     dp2=dpx*dpx

     g=   gpsi(1,isp)   +gpsi(2,isp)*dpx   +gpsi(3,isp)*dp2
     ri=  cpsi(1,isp)   +cpsi(2,isp)*dpx   +cpsi(3,isp)*dp2

     dqdp=               qpsi(2,isp)       +2.0*qpsi(3,isp)*dpx
     dgdp=               gpsi(2,isp)       +2.0*gpsi(3,isp)*dpx
     dridp=              cpsi(2,isp)       +2.0*cpsi(3,isp)*dpx

     gqi_inv=1.0/(g*qmesh(i)+ri)

! To make the 2nd derivatives of I and g continuous,
! currently we use linear interpolations as a temporary solution.
! In the future, higher order spline functions of I and g need to be introduced
     d2gdp2=2.0*(dpx*spdpsi_inv*(gpsi(3,isp+1)-gpsi(3,isp))+gpsi(3,isp))
     d2ridp2=2.0*(dpx*spdpsi_inv*(cpsi(3,isp+1)-cpsi(3,isp))+cpsi(3,isp))

     dsdp=(g*d2ridp2-ri*d2gdp2)/(g*qmesh(i)+ri) &
        -(g*dridp-ri*dgdp)*(dgdp*qmesh(i)+g*dqdp+dridp)*((gqi_inv)**2)

     ne=meshne(i)
     te=meshte(i)
     dnedp=-ne*kapane(i)
     dtedp=-te*kapate(i)
     phipEq=-mesher(i)! mesher = -dphi/dpsi

     do j=1,mtheta(i)
        ij=igrid(i)+j
        b_inv=1/bmesh(ij)

        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
        else
          tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
        endif
        jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
        dtx=tdum-spdtheta*real(jst-1)
        dt2=dtx*dtx

        dbdt = bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2+ &
             (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dtx*2.0
        dbdp = bsp(2,isp,jst)    +bsp(3,isp,jst)*dpx*2.0+ &
             (bsp(5,isp,jst)    +bsp(6,isp,jst)*dpx*2.0)*dtx+ &
             (bsp(8,isp,jst)    +bsp(9,isp,jst)*dpx*2.0)*dt2
        dbdz=0.0
        if (spdim==27) dbdz = bsp(10,isp,jst)+bsp(11,isp,jst)*dpx+bsp(12,isp,jst)*dp2+ &
             (bsp(13,isp,jst)+bsp(14,isp,jst)*dpx+bsp(15,isp,jst)*dp2)*dtx+ &
             (bsp(16,isp,jst)+bsp(17,isp,jst)*dpx+bsp(18,isp,jst)*dp2)*dt2

        dptdp=gradphi(1,1,ij)
        dptdt=gradphi(2,1,ij)
        dptdz=gradphi(3,1,ij)-gradphi(2,1,ij)*q_inv
        ddnedp=gradne(1,1,ij)
        ddnedt=gradne(2,1,ij)
        ddnedz=gradne(3,1,ij)-gradne(2,1,ij)*q_inv
        daparadp=gradapara(1,1,ij)
        daparadt=gradapara(2,1,ij)
        daparadz=gradapara(3,1,ij)-gradapara(2,1,ij)*q_inv
        dpsidp=gradpsi(1,1,ij)
        dpsidt=gradpsi(2,1,ij)
        dpsidz=gradpsi(3,1,ij)-gradpsi(2,1,ij)*q_inv

        dphieffdp=gradphieff(1,1,ij)
        dphieffdt=gradphieff(2,1,ij)
        dphieffdz=gradphieff(3,1,ij)-gradphieff(2,1,ij)*q_inv

        dne=sfluidne(1,ij)
        if(eta < machineEpsilon)then
           ! adiabatic pressure (either para or perp because they're equal)
           ddpdp=-qelectron*ne*dphieffdp+(te*dnedp+ne*dtedp+qelectron*ne*phipEq)*dpsidp
           ddpdt=-qelectron*ne*dphieffdt+(te*dnedp+ne*dtedp+qelectron*ne*phipEq)*dpsidt
           ddpdz=-qelectron*ne*dphieffdz+(te*dnedp+ne*dtedp+qelectron*ne*phipEq)*dpsidz
        else
           ! tearing mode
           ddpdp=te*ddnedp+dne*dtedp
           ddpdt=te*ddnedt
           ddpdz=te*ddnedz
        endif
   
        if(nhybrid>0)then
! parallel non-adiabatic pressure
          dkpparadp=gradpepara(1,1,ij)
          dkpparadt=gradpepara(2,1,ij)
          dkpparadz=gradpepara(3,1,ij)-gradpepara(2,1,ij)*q_inv
! perpendicular non-adiabatic pressure
          dkpperpdp=gradpeperp(1,1,ij)
          dkpperpdt=gradpeperp(2,1,ij)
          dkpperpdz=gradpeperp(3,1,ij)-gradpeperp(2,1,ij)*q_inv
        else
          dkpparadp=0.0
          dkpparadt=0.0
          dkpparadz=0.0
          dkpperpdp=0.0
          dkpperpdt=0.0
          dkpperpdz=0.0
        endif

        dinddp=dphieffdp-dptdp
        dinddt=dphieffdt-dptdt
        dinddz=dphieffdz-dptdz
        dudp=gradue(1,1,ij)
        dudt=gradue(2,1,ij)
        dudz=gradue(3,1,ij)-gradue(2,1,ij)*q_inv

! parallel compression
        wpara=jmesh(ij)*qmesh(i)*gradue(3,1,ij)
! wdrive: perturbed pressure gradient dot grad_B_0, interchange drive
        wdrive=jmesh(ij)*b_inv*b_inv*b_inv*(ri*dbdp*(ddpdz*2.0_lk+dkpparadz+dkpperpdz) &
              -ri*dbdz*(ddpdp*2.0_lk+dkpparadp+dkpperpdp) &
              +g*dbdt*(ddpdp*2.0_lk+dkpparadp+dkpperpdp)-g*dbdp*(ddpdt*2.0_lk+dkpparadt+dkpperpdt))
! wstar: perturbed ExB drift dot grad_n_0, omega_star term
        wstar=(ri*dptdz-g*dptdt)*dnedp*jmesh(ij)*b_inv*b_inv
! wdrift: perturbed ExB drift dot grad_B
        wdrift=2.0*(ne+perturb*dne)*jmesh(ij)*(g*dbdp*dptdt-g*dbdt*dptdp &
                    -ri*dbdp*dptdz+ri*dbdz*dptdp)*b_inv*b_inv*b_inv
! wdrift0: equilibrium radial ExB drift dot grad delta_n_e
        wdrift0=-2.0*dne*jmesh(ij)*(g*dbdt-ri*dbdz)*phipEq*b_inv*b_inv*b_inv
! web0: equilibrium radial ExB convection of perturbed density gradient
        web0=jmesh(ij)*(-ri*ddnedz+g*ddnedt)*phipEq*b_inv*b_inv
! weqc: equilibrium current terms (ion and fast ion flow not included yet)
        weqc=2.0/betae*rho0*rho0*jmesh(ij)*b_inv*dsdp*(ri*daparadz-g*daparadt) &
            -(dridp*(ddpdz+dkpparadz)-dgdp*(ddpdt+dkpparadt))/(g*qmesh(i)+ri) &
            +perturb*dne*(dridp*dptdz-dgdp*dptdt)/(g*qmesh(i)+ri)
! nonlinear compression along perturbed B-field (\delta B*\nabla(n_0u_\|/B_0))
        wperpnl=-jmesh(ij)*(ri*daparadz*dudp-ri*daparadp*dudz+g*daparadp*dudt-g*daparadt*dudp)*b_inv*perturb
! nonlinear exb convection
        wstarnl=jmesh(ij)*(ri*dptdz*ddnedp-ri*dptdp*ddnedz+g*dptdp*ddnedt-g*dptdt*ddnedp)*b_inv*b_inv*perturb

        if(iload==0)then ! ideal MHD: these two terms cancelled by ions in Poisson Eq.
           wstar=0.0
           wdrift=0.0
        else
           weqc=weqc+ne*(dridp*dptdz-dgdp*dptdt)*gqi_inv
        endif
        if(eqcurrent==0)weqc=0.0

!by zhs 41: without wdrive term
        !uncomment if necessary
        !wdrive=0.0
        wdamp=d4fluidne(ij)
        fluidne(1,ij) = fluidne0(1,ij) -dtime*rbc(i)*(wpara+wperpnl+wstarnl+wdrive+web0+wstar+wdrift+wdrift0+weqc-wdamp)
        apara(1,ij) = apara0(1,ij)+dtime*rbc(i)*(dinddt+qmesh(i)*dinddz)*b_inv*jmesh(ij)+&
          dtime*jmesh(ij)*qmesh(i)*gradext(3,1,ij)*b_inv
        deltapsi(1,ij) = deltapsi0(1,ij)-dtime*rbc(i)*dinddt*q_inv
!by Zhixuan: nonlinear genereation of zonal apara by dB\cdot \nabla \phi
        if(izonal>0) then
             JdBdotdE = daparadp*(ri*dptdz - g*dptdt) - dptdp*(ri*daparadz - g*daparadt) *b_inv*b_inv
             dBdotdE00 = dBdotdE00 +  JdBdotdE
        endif

!by Zhixuan: Gauge calculation to get MHD \xi_\perp, not needed for push particles, so no need to do Ronge-Kuta
!only ran when MHD displacement is needed for diagnosis (iMHDxi_output=1)
        iMHDxi_output=0
        if((irk==2) .and. (iMHDxi_output==1))then
            gradgaugef(1,1,ij)=gradgaugef(1,1,ij)+dinddp*dtime
            gradgaugef(2,1,ij)=gradgaugef(2,1,ij)+dinddt*dtime
            gradgaugef(3,1,ij)=gradgaugef(3,1,ij)+dinddz*dtime

!diagnosis only used when the step for snapshot
            if(mod(istep,mstep/msnap)==0 .or. istep*(1-irun)==ndiag)then

                MHDxiperp(1,1,ij)=(ri*gradgaugef(3,1,ij)-g*gradgaugef(2,1,ij))*gqi_inv
                MHDxiperp(2,1,ij)=g*gradgaugef(1,1,ij)*gqi_inv
                MHDxiperp(3,1,ij)=-ri*gradgaugef(1,1,ij)*gqi_inv
                MHDxi_mag(1,ij)=MHDxiperp(1,1,ij)*MHDxiperp(1,1,ij)*gdpp(ij)+&
                                MHDxiperp(2,1,ij)*MHDxiperp(2,1,ij)*gdtt(ij)+&
                                MHDxiperp(3,1,ij)*MHDxiperp(3,1,ij)*gdzz(ij)+&
                                MHDxiperp(1,1,ij)*MHDxiperp(2,1,ij)*gdpt(ij)*2.0
                MHDxi_mag(1,ij)=sqrt(MHDxi_mag(1,ij))
            endif
        endif

     enddo
     if(izonal>0)then
         apara00nl(i) = apara00nl0(i) + dtime*rbc(i)*dBdotdE00/jacobianpsi(i)/real(mtheta(i))
     endif
  enddo

  if (izonal>0)then
!calculate zonal <fluidne>
!$omp parallel do
     do i=0,mpsi
        do j=1,mtheta(i)
           ij=igrid(i)+j
           fluidne00(i)=fluidne00(i)+fluidne(0,ij)/jmesh(ij)
        enddo
        fluidne00(i)=fluidne00(i)/jacobianpsi(i)/real(mtheta(i))
     enddo

     call MPI_ALLREDUCE(apara00nl,adum,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierror)
     apara00nl=adum/real(mtoroidal)
     call MPI_ALLREDUCE(fluidne00,adum,mpsi+1,mpi_Rsize, MPI_SUM,toroidal_comm,ierror)
     fluidne00=adum/real(mtoroidal)

!calculate <A_para> and d<A_para>/dt
    do i=0,mpsi
       apara00(i)=aelectron*qion*meshni(i)*zonalci(i)/meshne(i)
       if(fload>0)then
          apara00(i)=apara00(i)+aelectron*qfast*meshnf(i)*zonalcf(i)/meshne(i)
       endif
       if(feload>0)then
          apara00(i)=apara00(i)+aelectron*qfaste*meshnfe(i)*zonalcfe(i)/meshnfe(i)
       endif
       if(nhybrid>0)then
          apara00(i)=apara00(i)+aelectron*qelectron*zonalce(i)
       endif
      apara00(i)=apara00(i)+apara00nl(i)*perturb
    enddo

!! (-0.0625 0.25 0.625 0.25 -0.0625) radial smoothing of (0,0) mode density apara00
    do i=1,2
       apara00tmp(0)=apara00(0)
       apara00tmp(mpsi)=apara00(mpsi)
       apara00tmp(1)=apara00(3)
       apara00tmp(mpsi-1)=apara00(mpsi-3)
       apara00tmp(2:mpsi-2)=apara00(0:mpsi-4)+apara00(4:mpsi)
       apara00tmp(1:mpsi-1)=0.625*apara00(1:mpsi-1)+0.25*(apara00(0:mpsi-2)+apara00(2:mpsi))-&
            0.0625*apara00tmp(1:mpsi-1)
       apara00=apara00tmp
    enddo
  else
    apara00=0.0
  endif

  call periodicity(fluidne)
  call periodicity(apara)
  call periodicity(deltapsi)
  call periodicity(MHDxi_mag)
  call periodicity(MHDxiperp(1,:,:))
  call periodicity(MHDxiperp(2,:,:))
  call periodicity(MHDxiperp(3,:,:))

  if(psi0==0.0_lk)then
    call continuousbc(fluidne)
    call continuousbc(apara)
    call continuousbc(deltapsi)
  endif

! smooth fields
  sfluidne=fluidne
  sdelapara=apara ! A|| without zonal current component
  sdeltapsi=deltapsi

  if (izonal>0) then
!$omp parallel do private(i,j,ij)
     do i=0,mpsi
        do j=1,mtheta(i)
           ij=igrid(i)+j
           sfluidne(1,ij)=sfluidne(1,ij)-fluidne00(i)
        enddo
     enddo
  endif

  call smooth(sfluidne)
  call smooth(sdelapara)
  call smooth(sdeltapsi)

  if(nfilter>0)then
    CALL FILTER(sfluidne)
    CALL FILTER(sdelapara)
    CALL FILTER(sdeltapsi)
  endif

  if(psi0==0.0_lk)then
    call continuousbc(sfluidne)
    call continuousbc(sdelapara)
    call continuousbc(sdeltapsi)
  endif

!add <A_para> with /delta A_para and expand d<A_para>/dt to poloidal plane
  if(nonlinear==1)then
    do i=0,mpsi
       sapara(0:1,igrid(i):igrid(i)+mtheta(i))=sdelapara(0:1,igrid(i):igrid(i)+mtheta(i))+apara00(i)
    enddo
  else
    sapara=sdelapara
  endif

end subroutine pushfield
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine continuousbc(farray)
  use global_parameters
  use field_array,only:igrid,mtheta
  implicit none

  integer ::i,j,ij
  real(lk)::farray(0:1,mgrid),ave

  i=1 
  ave=0.0_lk
  do j=1,mtheta(i)
    ij=igrid(i)+j
    ave = ave + farray(1,ij)/real(mtheta(i))
  enddo

  i=0
  do j=0,mtheta(i)
    ij=igrid(i)+j
    farray(1,ij)=ave
  enddo

end subroutine continuousbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine periodicity(fgrid)
  use global_parameters
  use field_array
  implicit none

  integer i,ii,jt,icount,idest,isource,isendtag,irecvtag,ierror,istatus(MPI_STATUS_SIZE)
  real(lk) sendr(mgrid),recvl(mgrid)
  real(lk),dimension(0:1,mgrid) :: fgrid

! toroidal end point, send E to right and receive from left
!$omp parallel do private(i)
  do i=1,mgrid
     sendr(i)=fgrid(1,i)
     recvl(i)=0.0
  enddo
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! unpack end point data for k=0
  if(myrank_toroidal==0)then
!$omp parallel do private(i,ii,jt)
     do i=0,mpsi
        ii=igrid(i)
        jt=mtheta(i)
        fgrid(0,ii+1:ii+jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
     enddo
  else
!$omp parallel do private(i)
     do i=1,mgrid
        fgrid(0,i)=recvl(i)
     enddo
  endif

! poloidal end point
!$omp parallel do private(i,ii,jt)
  do i=0,mpsi
     ii=igrid(i)
     jt=mtheta(i)
     fgrid(0:1,ii)=fgrid(0:1,ii+jt)
  enddo

end subroutine periodicity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
