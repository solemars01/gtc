#ifdef _hybrid
#define chargeParticle hybridChargeParticle
#else
#define chargeParticle gkChargeParticle
#endif

subroutine chargeParticle(zpart,wppart,wtpart0,wtpart1,density,flow,&
    jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,qpart,apart,&
    pload,ngyro,mp,ihybrid,nhybrid,pressurepara,pressureperp,&
    sfluidn,dnsave,switch)
  use precision
  use global_parameters,only: magnetic,npartdom,partd_comm,left_pe,&
    right_pe,myrank_toroidal,toroidal_comm,mtoroidal,irk,nbound,nboundR,&
    istep,mstep,nfilter,mpsi,mgrid,antenna,etemp0,eden0,r0,mype,izonal
  use field_array,only: igrid,mtheta,itran,dn_ext,rhom,guzz,nmodes,nmode,zeta0
  use equilibrium,only: lsp,spdpsi_inv,spdpsi,lst,spdtheta_inv,spdtheta,&
    bsp,xsp,spdim
  implicit none

  !TODO: check decleartions 
  !declaration of the dummy arguments
  integer pload,ngyro,mp
  integer,optional :: ihybrid,nhybrid
  integer,dimension(:,:) :: jtpart0,jtpart1
  real(lk) qpart,apart
  real(lk),dimension(:) :: wzpart,marker,markert
  real(lk),dimension(0:) :: zonal,zonalc,pmark
  real(lk),dimension(:,:) :: wppart,wtpart0,wtpart1
  real(lk),dimension(0:,:) :: density,flow
  real(lk),dimension(:,:) :: zpart
  real(lk),dimension(0:,:),optional :: pressurepara,pressureperp,sfluidn
  real(lk),dimension(0:,:,:),optional :: dnsave
  character(*),intent(in),optional :: switch

  !declaration of the local variables
  integer m,i,j,jt,ii,ierror,igyro,ij,isp,jst,&
    icount,idest,isource,isendtag,irecvtag,istatus(MPI_STATUS_SIZE),ia
  real(lk) gyrodum,weight,b,upara,cmratio,wz1,wz0,wp1,wp0,wt11,wt10,wt01,&
    wt00,adum(0:mpsi),dnitmp(0:1,mgrid),djitmp(0:1,mgrid),dpx,dp2,dzx,&
    dx(27),pdum,tdum,dtx,dt2,zf0,zc0,temp,majorr,dpitmppara(0:1,mgrid),&
    dpitmpperp(0:1,mgrid),energypara,energyperp
#ifdef _hybrid
  real(lk) sendl(mgrid,4),recvr(mgrid,4)
#else
  real(lk) sendl(mgrid,2),recvr(mgrid,2)
#endif
#ifdef _FRC
  real(8) dbesjn,zetagyravg,kperprhoi
#else
  real(lk) zetagyravg,kperprhoi
#endif
!  zetagyravg=1.0
!  kperprhoi=0.0  

  cmratio=qpart/apart
#ifdef _OPENACC
  !$acc parallel loop gang vector
#else
  !$omp parallel do private(ij)
#endif
  do ij=1,mgrid
    density(0:1,ij)=0.0_lk
    flow(0:1,ij)=0.0_lk
#ifdef _hybrid
    pressurepara(0:1,ij)=0.0_lk
    pressureperp(0:1,ij)=0.0_lk
#endif
  enddo
  !$acc end parallel
  
  if(magnetic==0)then
#ifdef _OPENACC
    !$acc parallel loop gang vector
#else
    ! scatter particle density for electrostatic simulation
    !$omp parallel do&
    !$omp& private(m,igyro,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij,dzx,dx,&
    !$omp& b,majorr,kperprhoi,zetagyravg,pdum,isp,dpx,tdum,jst,dtx,dt2,dp2)&
    !$omp& reduction(+: density)
#endif
     do m=1,mp

        pdum=zpart(1,m)
        isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
        dpx=pdum-spdpsi*real(isp-1)
        ! radial spline of b avoids sigularity near axis:
        ! y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
        if(isp==1)dpx=sqrt(dpx)
        dp2=dpx*dpx

        tdum=zpart(2,m)
        jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
        dtx=tdum-spdtheta*real(jst-1)
        dt2=dtx*dtx
        dzx=zpart(3,m)-zeta0

        dx(1)=1.0_lk
        dx(2)=dpx
        dx(3)=dp2
        dx(4:6)=dx(1:3)*dtx
        dx(7:9)=dx(1:3)*dt2
        dx(10:18)=dx(1:9)*dzx
        dx(19:27)=dx(1:9)*dzx*dzx

        b=0.0_lk
        majorr=0.0_lk
        do ii = 1, spdim
          b=b +bsp(ii,isp,jst)*dx(ii)
          majorr=majorr +xsp(ii,isp,jst)*dx(ii)
        enddo

        zetagyravg=1.0_lk
#ifdef _FRC
        !toroidal gyroavg by multiplying by besselj(0,k_perp*rho_i)
        !rho_i from sqrt(mu)*sqrt(2)/sqrt(b)/q->zpart(6,m)*sqrt(2)/sqrt(b)
        !k_perp from n/(Rmaj)--> n/majorr
        kperprhoi=(zpart(6,m)*sqrt(2.0_lk*apart/b))*(real(nmode)/majorr)
        zetagyravg=dbesjn(0,kperprhoi)
#endif
        weight=zpart(5,m)*zetagyravg
     
        wz1=weight*wzpart(m)      !weight for upper toroidal grid
        wz0=weight-wz1            !weight for lower toroidal grid
        do igyro=1,ngyro
          wp1=wppart(igyro,m)     !outer flux surface
          wp0=1.0_lk-wp1             !inner flux surface

          wt10=wp0*wtpart0(igyro,m) !upper poloidal grid on inner flux surface
          wt00=wp0-wt10           !lower poloidal grid on inner flux surface
           
          wt11=wp1*wtpart1(igyro,m) !upper poloidal grid on outer flux surface
          wt01=wp1-wt11           !lower poloidal grid on outer flux surface

          ! If no loop-level parallelism, write directly into array "density()"
          ij=jtpart0(igyro,m)     !lower poloidal grid on inner flux surface
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt00 !lower toroidal grid
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt00 !upper toroidal grid
           
          ij=ij+1                 !upper poloidal grid on inner flux surface
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt10
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt10
           
          ij=jtpart1(igyro,m)     !lower poloidal grid on outer flux surface
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt01
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt01
           
          ij=ij+1                 !upper poloidal grid on outer flux surface
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt11
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt11
        enddo
      enddo
    else
    ! scatter particle density and flow for electromagnetic simulation 
#ifdef _OPENACC
    !$acc parallel loop gang vector
#else
#ifdef _hybrid
    !$omp parallel do private(m,igyro,pdum,isp,dpx,dp2,tdum,jst,dtx,dt2,b,&
    !$omp& upara,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)&
    !$omp& reduction(+: density,flow)&
    !$omp& private(energypara,energyperp)&
    !$omp& reduction(+: pressurepara,pressureperp)
#else
    !$omp parallel do private(m,igyro,pdum,isp,dpx,dp2,tdum,jst,dtx,dt2,b,&
    !$omp& upara,weight,wz1,wz0,wp1,wp0,wt10,wt00,wt11,wt01,ij)&
    !$omp& reduction(+: density,flow)
#endif
#endif
     do m=1,mp

        ! 2D spline in (psi, theta)
        pdum=zpart(1,m)
        isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
        dpx=pdum-spdpsi*real(isp-1)

        tdum=zpart(2,m)
        jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
        dtx=tdum-spdtheta*real(jst-1)
        dt2=dtx*dtx

        ! radial spline of b avoids sigularity near axis:
        ! y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
        if(isp==1)dpx=sqrt(dpx)
        dp2=dpx*dpx

        b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
             (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
             (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
        
        upara=zpart(4,m)*b*cmratio !parallel velocity
#ifdef _hybrid
        energypara=apart*upara*upara
        energyperp=zpart(6,m)*zpart(6,m)*b
#endif
        weight=zpart(5,m)
        
        wz1=weight*wzpart(m)
        wz0=weight-wz1     
        !Todo: why reduction not work for the following loop?
        !!$acc loop reduction(+:density,flow)
        !$acc loop seq
        do igyro=1,ngyro
          wp1=wppart(igyro,m)       !outer flux surface
          wp0=1.0_lk-wp1               !inner flux surface
           
          wt10=wp0*wtpart0(igyro,m)
          wt00=wp0-wt10
           
          wt11=wp1*wtpart1(igyro,m)
          wt01=wp1-wt11

          ! If no loop-level parallelism, use original algorithm (write
          ! directly into array "density()".
          ij=jtpart0(igyro,m)
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt00
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt00
          !$acc atomic update
          flow(0,ij) = flow(0,ij) + wz0*wt00*upara
          !$acc atomic update
          flow(1,ij) = flow(1,ij) + wz1*wt00*upara
#ifdef _hybrid
          !$acc atomic update
          pressurepara(0,ij) = pressurepara(0,ij) + wz0*wt00*energypara
          !$acc atomic update
          pressurepara(1,ij) = pressurepara(1,ij) + wz1*wt00*energypara
          !$acc atomic update
          pressureperp(0,ij) = pressureperp(0,ij) + wz0*wt00*energyperp
          !$acc atomic update
          pressureperp(1,ij) = pressureperp(1,ij) + wz1*wt00*energyperp
#endif
         
          ij=ij+1
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt10
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt10
          !$acc atomic update
          flow(0,ij) = flow(0,ij) + wz0*wt10*upara
          !$acc atomic update
          flow(1,ij) = flow(1,ij) + wz1*wt10*upara
#ifdef _hybrid
          !$acc atomic update
          pressurepara(0,ij) = pressurepara(0,ij) + wz0*wt10*energypara
          !$acc atomic update
          pressurepara(1,ij) = pressurepara(1,ij) + wz1*wt10*energypara
          !$acc atomic update
          pressureperp(0,ij) = pressureperp(0,ij) + wz0*wt10*energyperp
          !$acc atomic update
          pressureperp(1,ij) = pressureperp(1,ij) + wz1*wt10*energyperp
#endif
    
          ij=jtpart1(igyro,m)
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt01
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt01
          !$acc atomic update
          flow(0,ij) = flow(0,ij) + wz0*wt01*upara
          !$acc atomic update
          flow(1,ij) = flow(1,ij) + wz1*wt01*upara
#ifdef _hybrid
          !$acc atomic update
          pressurepara(0,ij) = pressurepara(0,ij) + wz0*wt01*energypara
          !$acc atomic update
          pressurepara(1,ij) = pressurepara(1,ij) + wz1*wt01*energypara
          !$acc atomic update
          pressureperp(0,ij) = pressureperp(0,ij) + wz0*wt01*energyperp
          !$acc atomic update
          pressureperp(1,ij) = pressureperp(1,ij) + wz1*wt01*energyperp
#endif
           
          ij=ij+1
          !$acc atomic update
          density(0,ij) = density(0,ij) + wz0*wt11
          !$acc atomic update
          density(1,ij) = density(1,ij) + wz1*wt11
          !$acc atomic update
          flow(0,ij) = flow(0,ij) + wz0*wt11*upara
          !$acc atomic update
          flow(1,ij) = flow(1,ij) + wz1*wt11*upara
#ifdef _hybrid
          !$acc atomic update
          pressurepara(0,ij) = pressurepara(0,ij) + wz0*wt11*energypara
          !$acc atomic update
          pressurepara(1,ij) = pressurepara(1,ij) + wz1*wt11*energypara
          !$acc atomic update
          pressureperp(0,ij) = pressureperp(0,ij) + wz0*wt11*energyperp
          !$acc atomic update
          pressureperp(1,ij) = pressureperp(1,ij) + wz1*wt11*energyperp
#endif
        enddo
     enddo
!$acc end parallel
  endif

#ifndef GPU_UM
  !$acc update host(density,flow)
#endif
#ifdef _hybrid
#ifndef GPU_UM
  !$acc update host(pressurepara,pressureperp)
#endif
#endif
  ! If we have a particle decomposition on the toroidal domains, do a reduce
  ! operation to add up all the contributions to charge density on the grid
  if(npartdom>1)then
  !$omp parallel do private(ij)
    do ij=1,mgrid
      dnitmp(0:1,ij)=density(0:1,ij)
      density(0:1,ij)=0.0_lk
      djitmp(0:1,ij)=flow(0:1,ij)
      flow(0:1,ij)=0.0_lk
    enddo
    icount=2*mgrid
    call MPI_ALLREDUCE(dnitmp,density,icount,mpi_Rsize,MPI_SUM,partd_comm,ierror)
    call MPI_ALLREDUCE(djitmp,flow,icount,mpi_Rsize,MPI_SUM,partd_comm,ierror)

#ifdef _hybrid
    !$omp parallel do private(ij)
    do ij=1,mgrid
      dpitmppara(0:1,ij)=pressurepara(0:1,ij)
      pressurepara(0:1,ij)=0.0_lk
      dpitmpperp(0:1,ij)=pressureperp(0:1,ij)
      pressureperp(0:1,ij)=0.0_lk
    enddo
    call MPI_ALLREDUCE(dpitmppara,pressurepara,icount,mpi_Rsize,MPI_SUM,partd_comm,ierror)
    call MPI_ALLREDUCE(dpitmpperp,pressureperp,icount,mpi_Rsize,MPI_SUM,partd_comm,ierror)
#endif
  endif

  ! poloidal end cell, discard ghost cell j=0
  !$omp parallel do private(i)
  do i=0,mpsi
    density(:,igrid(i)+mtheta(i))=density(:,igrid(i)+mtheta(i))+density(:,igrid(i))
    flow(:,igrid(i)+mtheta(i))=flow(:,igrid(i)+mtheta(i))+flow(:,igrid(i))
#ifdef _hybrid
    pressurepara(:,igrid(i)+mtheta(i))=pressurepara(:,igrid(i)+mtheta(i))+pressurepara(:,igrid(i))
    pressureperp(:,igrid(i)+mtheta(i))=pressureperp(:,igrid(i)+mtheta(i))+pressureperp(:,igrid(i))
#endif
  enddo

  ! toroidal end cell
  !$omp parallel do private(i)
  do i=1,mgrid
    sendl(i,1)=density(0,i)
    sendl(i,2)=flow(0,i)
    recvr(i,1:2)=0.0_lk
#ifdef _hybrid
    sendl(i,3)=pressurepara(0,i)
    sendl(i,4)=pressureperp(0,i)
    recvr(i,3:4)=0.0_lk
#endif
  enddo

#ifdef _hybrid
  icount=4*mgrid
#else
  icount=2*mgrid
#endif
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource

  ! send density to left and receive from right
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
       recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)
     
  if(myrank_toroidal == mtoroidal-1)then
    ! B.C. at zeta=2*pi is shifted
    !$omp parallel do private(i,ii,jt)
    do i=0,mpsi
      ii=igrid(i)
      jt=mtheta(i)
      density(1,ii+1:ii+jt)=density(1,ii+1:ii+jt)+cshift(recvr(ii+1:ii+jt,1),itran(i))
      flow(1,ii+1:ii+jt)=flow(1,ii+1:ii+jt)+cshift(recvr(ii+1:ii+jt,2),itran(i))
#ifdef _hybrid
      pressurepara(1,ii+1:ii+jt)=pressurepara(1,ii+1:ii+jt)+&
           cshift(recvr(ii+1:ii+jt,3),itran(i))
      pressureperp(1,ii+1:ii+jt)=pressureperp(1,ii+1:ii+jt)+&
           cshift(recvr(ii+1:ii+jt,4),itran(i))
#endif
    enddo
  else
    ! B.C. at zeta<2*pi is continuous
    !$omp parallel do private(i,ii,jt)
    do i=1,mgrid
      density(1,i)=density(1,i)+recvr(i,1)
      flow(1,i)=flow(1,i)+recvr(i,2)
#ifdef _hybrid
      pressurepara(1,i)=pressurepara(1,i)+recvr(i,3)
      pressureperp(1,i)=pressureperp(1,i)+recvr(i,4)
#endif
    enddo
  endif
  
  !self-consistent density from external antenna potential
  if(present(switch) .and. switch=='density modification' .and. istep==1)then
    do ia=1,antenna
      do i=0,mpsi
        temp=(etemp0)/(eden0)*7.43e2*7.43e2/(r0*r0)
        do j=1,mtheta(i)
          ij=igrid(i)+j
          density(1,ij)=dn_ext(1,ij)*temp*real(pmark(i)/&
            (mtheta(i)*mtoroidal))+density(1,ij)
        enddo
      enddo
    enddo
  endif
  ! flux surface average and normalization  
  gyrodum=1.0_lk/real(ngyro)
  !$omp parallel do private(i,j,ij)
  do i=0,mpsi
    zonal(i)=0.0_lk
    zonalc(i)=0.0_lk
    do j=1,mtheta(i)
      ij=igrid(i)+j
      density(1,ij)=gyrodum*density(1,ij)
      zonal(i)=zonal(i)+density(1,ij)
      flow(1,ij)=gyrodum*flow(1,ij)
      zonalc(i)=zonalc(i)+flow(1,ij)
    enddo
  enddo

  ! time-averaged # of marker particles per cell ****This is currently not being used. Keep for future?
  !if(irk==2 .and. ( .not. present(nhybrid) .or. present(nhybrid) .and. ihybrid==nhybrid))markert=markert+density(1,:)

  !$omp parallel do private(i,j,ij)
  do i=0,mpsi
    do j=1,mtheta(i)
      ij=igrid(i)+j
      density(1,ij)=density(1,ij)*marker(ij)
      flow(1,ij)=flow(1,ij)*marker(ij)
#ifdef _hybrid
      pressurepara(1,ij)=pressurepara(1,ij)*marker(ij)
      pressureperp(1,ij)=pressureperp(1,ij)*marker(ij)
#endif
    enddo
  enddo
  
  ! global sum of zonal modes (flux surface averaged), broadcast to every toroidal PE
  call MPI_ALLREDUCE(zonal,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  zonal=adum/pmark
  zf0=sum(adum)

  call MPI_ALLREDUCE(zonalc,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  zonalc=adum/pmark
  zc0=sum(adum)

  !density subtracted by (0,0) mode
  if(magnetic==1)then
    !$omp parallel do private(i,j,ij)
    do i=0,mpsi
     if(izonal>0)then
        do j=1,mtheta(i)
          ij=igrid(i)+j
          density(1,ij)=density(1,ij)-zonal(i)
          flow(1,ij)=flow(1,ij)-zonalc(i)
        enddo
      endif
      ! poloidal BC condition
      density(1,igrid(i))=density(1,igrid(i)+mtheta(i))
      flow(1,igrid(i))=flow(1,igrid(i)+mtheta(i))
    enddo
  endif

  if(pload/=9)then         ! For GK calculation
    if(pload<100)then
      ! enforce charge/momentum conservation for zonal flow/field mode in delta-f simulation
      zc0=zc0/sum(pmark)
      zf0=zf0/sum(pmark)
      zonal(1:mpsi-1)=zonal(1:mpsi-1)-zf0
      zonalc(1:mpsi-1)=zonalc(1:mpsi-1)-zc0
    else
      ! full-f: zonal flow subtracted by equilibrium density
      zonal=zonal-1.0
    endif
    ! zero out zonal charge/current in radial boundary cell
    if(nbound>99)then
      call gaussbound(zonal)
      call gaussbound(zonalc)
      if(nboundR>99)then
        call gaussboundR(zonal)
        call gaussboundR(zonalc)
      endif
    else
      do i=0,nbound-1
        zonal(i)=zonal(i)*real(i)/real(nbound)
        zonal(mpsi-i)=zonal(mpsi-i)*real(i)/real(nbound)
        zonalc(i)=zonalc(i)*real(i)/real(nbound)
        zonalc(mpsi-i)=zonalc(mpsi-i)*real(i)/real(nbound)
      enddo
    endif
  else   ! For fully kinetic calculation
    ! full-f: zonal flow subtracted by equilibrium density
    zonal=zonal-1.0_lk
    ! This scheme is not correct for the delta-f simulation
  endif

  !if(irk==2 .and. istep==mstep .and. ( .not. present(nhybrid) .or. present(nhybrid) .and. ihybrid==nhybrid)) then
  !  markert=markert/real(mstep) ! # of marker particles per cell
  !endif

  ! smooth particle density and current  
  call smooth(density)
#ifdef _FRC
  if(nfilter>=2)call mfilter_frc(density)
#else
  if(nfilter>0)call filter(density)
#endif
  if(magnetic==1)then
    call smooth(flow)
#ifdef _hybrid
    call smooth(pressurepara)
    call smooth(pressureperp)
#endif
    if(nfilter>0)then
      call filter(flow)
#ifdef _hybrid
      call filter(pressurepara)
      call filter(pressureperp)
#endif
    endif
  endif
end subroutine chargeParticle
