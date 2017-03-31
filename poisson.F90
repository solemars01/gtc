subroutine poisson_solver(scenario)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use petsc_array,only: userb,userx,luserb,luserx
  use interfaces,only: axisExtrapolate
  use spline_function, only: spb, spx, spz
  implicit none

  integer :: i,ij,ia,j,mring,mindex,mtest,ierr,iteration,it
  integer,dimension(:),allocatable :: nindex
  integer,dimension(:,:),allocatable :: indexp
  real(lk),dimension(:,:),allocatable :: ring
  real(lk),dimension(mgrid) :: phi_ext
  real(lk) dentmp(mgrid),&
           wt,ave,cosave,sinave,tdum,rdum,ave2,sinave2,cosave2,ddum(0:mpsi)
  character(*),intent(in) :: scenario
  real(lk) :: temp,gam
  real(lk),dimension(:),allocatable :: tmp1,tmp2,phitmp,ptilde,perr
  save nindex,indexp,ring

  interface
     subroutine poisson_solver_initial(mring,mindex,nindex,indexp,ring)
       use precision
       use global_parameters,only: mgrid
       implicit none

       integer :: mring,mindex
       integer,dimension(mgrid) :: nindex
       integer,dimension(mindex,mgrid) :: indexp
       real(lk),dimension(mindex,mgrid) :: ring
     end subroutine poisson_solver_initial
  end interface

! number of gyro-ring
  mring=2
! number of summation: maximum is 32*mring+1
  mindex=32*mring+1

! initialize poisson solver
  if(istep==1 .and. irk==1 .and. scenario/='ES-hybrid')then
! if(istep==1 .and. irk==1 .and. scenario/='ES-hybrid-old')then  !ES-hybrid-old refers to resolving the poisson Eq per hybrid loop
    allocate(indexp(mindex,mgrid),ring(mindex,mgrid),nindex(mgrid),STAT=mtest)
    if (mtest /= 0) then
        write(0,*)mype,'*** Cannot allocate indexp: mtest=',mtest
        call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
    endif

! initialize
    call poisson_solver_initial(mring,mindex,nindex,indexp,ring)

  endif

!i.h. modification for Poisson's equation
  ddum=1.0_lk !ddum=1.0 for ilaplacian==0
  if (ilaplacian==1) then
!$omp parallel do private(i)
    do i=0,mpsi
      if (fload==0) then
        ddum(i)=qion*qion*meshni(i)/meshti(i)
      else
        ddum(i)=qion*qion*(meshni(i)+meshnf(i)*afast/aion)/meshti(i)
      endif
    enddo
  endif

  if(scenario=='ES-adiabatic')then
    !$omp parallel do private(i,j,ij)
    do i=0,mpsi
      do j=0,mtheta(i)
        ij=igrid(i)+j
        dentmp(ij)=qion*densityi(1,ij)*meshni(i)/(1.0_lk+real(ilaplacian, lk)*qelectron*qelectron*meshne(i)/(meshte(i)*ddum(i)))
      enddo
    enddo
  elseif(scenario=='ES-hybrid')then
    !$omp parallel do private(i,j,ij)
    do i=0,mpsi
      do j=0,mtheta(i)
        ij=igrid(i)+j
           phi(1,ij)=phi_zero(1,ij)+meshte(i)*densitye(1,ij)/qelectron 
      enddo
    enddo
  elseif(scenario=='ES-hybrid-old')then
    !$omp parallel do private(i,j,ij)
    do i=0,mpsi
      do j=0,mtheta(i)
        ij=igrid(i)+j
        dentmp(ij)=(qion*densityi(1,ij)*meshni(i)+qelectron*densitye(1,ij)*meshne(i))/&
                   (1.0_lk+real(ilaplacian,lk)*qelectron*qelectron*meshne(i)/(meshte(i)*ddum(i)))
      enddo
    enddo
  elseif(scenario=='ES-DKE')then
    !$omp parallel do private(i,j,ij)
    do i=0,mpsi
      do j=0,mtheta(i)
        ij=igrid(i)+j
        dentmp(ij)=qion*densityi(1,ij)*meshni(i)
      enddo
    enddo
  elseif(scenario=='EM-hybrid')then
    !$omp parallel do private(i,j,ij)
    do i=0,mpsi
      do j=0,mtheta(i)
        ij=igrid(i)+j
        dentmp(ij)=qion*densityi(1,ij)*meshni(i)+qelectron*sfluidne(1,ij)
      enddo
    enddo
  endif
if(iload==9)then
!$omp parallel do private(i,j,ij,temp)
     do i=0,mpsi
           temp=(etemp0/eden0)*7.43e2*7.43e2/(r0*r0)
           temp=1.0_lk
! normalized by (lambdaD_i)^2/(R_0)^2, lambdaDi, is the Debye-length, and R_0, is the major radius
       do j=1,mtheta(i)
           ij=igrid(i)+j
           dentmp(ij)=dentmp(ij)/temp
       enddo
     enddo
endif

! add fast ion contribution
  if(fload>0)then
!$omp parallel do private(i,j,ij,temp)
     do i=0,mpsi
        temp=1.0_lk+real(1-magnetic,lk)*real(ilaplacian,lk)*qelectron*qelectron*meshne(i)/(meshte(i)*ddum(i))
        do j=0,mtheta(i)
           ij=igrid(i)+j
           dentmp(ij)=dentmp(ij)+qfast*densityf(1,ij)*meshnf(i)/temp
        enddo
     enddo
  endif

!add fast elelctron contribution
   if(feload>0)then
!$omp parallel do private(i,j,ij,temp)
     do i=0,mpsi
        temp=1.0_lk+real(1-magnetic,lk)*real(ilaplacian,lk)*qelectron*qelectron*meshne(i)/(meshte(i)*ddum(i))
        do j=0,mtheta(i)
           ij=igrid(i)+j
           dentmp(ij)=dentmp(ij)+qfaste*densityfe(1,ij)*meshnfe(i)/temp
        enddo
     enddo
   endif

  if(scenario/='ES-hybrid')then  !Solve Poisson Eq for ion's and adiabatic electrons only
! Solve Poisson Equation using PETSc
#ifdef _PETSc
    if(ilaplacian==0 .and. fem==0)then
      !$omp parallel do private(i)
      do i=1,mgrid
        userb(i-1)=dentmp(i)  !!XY
      enddo
      call petsc_solver
      !$omp parallel do private(i)
      do i=1,mgrid
        phi(1,i)=userx(i-1)
      enddo
    elseif(ilaplacian==1 .and. fem==0)then
      if(psi0>0.0_lk)dentmp(igrid(0):igrid(0)+mtheta(0))=0.0_lk
      dentmp(igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0_lk
      !$omp parallel do private(i)
      do i=mgridlow,mgridhigh
        luserb(i-mgridlow)=dentmp(i)
      enddo
      call lapmat_pssolver
      !$omp parallel do private(i)
      do i=mgridlow,mgridhigh
        phi(1,i)=luserx(i-mgridlow)
      enddo
    elseif(fem>0)then
      call laplacian_fem(0,dentmp,phi(1,:)) 
    endif
! radial boundary
    if(psi0>0.0_lk)phi(1,igrid(0):igrid(0)+mtheta(0))=0.0_lk
    phi(1,igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0_lk

#else

! Solve Poisson Equation using iterative method (Lin & Lee, PRE95)
    allocate (tmp1(0:mpsi),tmp2(0:mpsi),perr(mgrid),ptilde(mgrid),phitmp(mgrid))
    gam=0.75 ! max. resolution for k=0.577
    iteration=5

! first iteration, first guess of phi. (1+T_i/T_e) phi - phi_title = n_i
!$omp parallel do private(i,j,ij)
    do i=0,mpsi
      if(fload>0)then
        tmp1(i)=qion*qion*(meshni(i)+meshnf(i)*afast/aion)*rho0*rho0/meshti(i)
      else
        tmp1(i)=qion*qion*meshni(i)*rho0*rho0/meshti(i)
      endif
      if(feload>0)then
        tmp2(i)=1.0_lk/(tmp1(i)*(1.0_lk-gam)+meshne(i)*rho0*rho0/meshte(i))+meshnfe(i)*rho0*rho0/meshte(i)
      else
        tmp2(i)=1.0_lk/(tmp1(i)*(1.0_lk-gam)+meshne(i)*rho0*rho0/meshte(i))
      endif

      do j=1,mtheta(i)
        ij=igrid(i)+j
        phitmp(ij)=dentmp(ij)*tmp2(i)
      enddo
    enddo

    do it=2,iteration
!$omp parallel do private(i,j)
      do i=1,mgrid
        ptilde(i)=0.0_lk
        do j=1,nindex(i)
          ptilde(i)=ptilde(i)+ring(j,i)*phitmp(indexp(j,i))
        enddo
      enddo

!$omp parallel do private(i,j,ij)
      do i=0,mpsi
        do j=1,mtheta(i)
          ij=igrid(i)+j
          perr(ij)=(ptilde(ij)-gam*phitmp(ij))*tmp1(i)
          phitmp(ij)=(dentmp(ij)+perr(ij))*tmp2(i)
        enddo
      enddo
! radial boundary
      phitmp(igrid(0):igrid(0)+mtheta(0))=0.0_lk
      phitmp(igrid(mpsi):igrid(mpsi)+mtheta(mpsi))=0.0_lk
    enddo

! store final results
!$omp parallel do private(i)
    do i=1,mgrid
      phi(1,i)=phitmp(i)
    enddo

#endif

    if (antenna>0) then
    ! antenna=1: launch external anntena, use standing wave in cos(m*theta-n*zeta)*cos(omega*t)
      phi_ext=0.0_lk
      do ia=1,antenna
        wt=cos(omega(ia)*tstep*(real(istep+mstep*irun)+0.5*irk)) !irun= # restart runs
!$omp parallel do private(i)
        do i=1,mgrid
          phi_ext(i)=phi_ext(i)+phiext(ia,i)*wt
        enddo

        ! write out driving source
        if(mype==0 .and. istep==1 .and. irk==1)then
          open(9991,file='time_ext.out',status='replace')
          open(9992,file='phi_ext.out',status='replace')
        endif
        if(mype==0 .and. irk==1)then
          write(9991,*)wt
          write(9992,*)dn_ext(1,igrid(mpsi/2)+1),phi_ext(igrid(mpsi/2)+1)
        endif
        if(mype==0 .and. istep==mstep .and. irk==1)then
          close(9991)
          close(9992)
        endif
      enddo

!$omp parallel do private(i)
      do i=1,mgrid
        phi(1,i)=phi(1,i)+phi_ext(i)
      enddo
    endif

    ! in GTC unit
    if (ilaplacian==1) then
!$omp parallel do private(i,j,ij)
      do i=0,mpsi
        do j=0,mtheta(i)
          ij=igrid(i)+j
          phi(1,ij)=phi(1,ij)*rho0*rho0+dentmp(ij)/ddum(i)
        enddo
      enddo
    else
!$omp parallel do private(i)
      do i=1,mgrid
        phi(1,i)=phi(1,i)*rho0*rho0
      enddo
    endif
  endif 
  if(bcond==1)call axisExtrapolate(phi(1,:))
! smooth potential
  CALL SMOOTH(phi)

! use one or a few toroidal modes
  if(nfilter>0)CALL FILTER(phi)

!!save phi adiabatic to be corrected in hybrid loop
  if(scenario=='ES-adiabatic')phi_zero=phi
end subroutine poisson_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine poisson_solver_initial(mring,mindex,nindex,indexp,ring)
  use precision
  use global_parameters
  use field_array
  use particle_array
  use equilibrium
  use petsc_array,only: userp,users,usera,userb,userx
  implicit none

  integer mring,mindex,ierr,nindex0max
  integer,dimension(mgrid) :: nindex
  integer,dimension(mindex,mgrid) :: indexp
  real(lk),dimension(mindex,mgrid) :: ring
  integer i,ii,ij,ij0,j,jt,j1,j0,kr,kp,nt,np,isp,ipjt,iflagboundary,icorrection
  real(lk) vring(3),fring(3),rg,ddelr,ddelt,wght,rr,rdum,wr,&
       wt1,wt0,tdum,zdum,b,pdum,pdum0,tdum0,t1,dpx,tempi,rhogyro

#ifdef _PETSc
  interface
     subroutine matrix(mindex,nindex,indexp,ring,nindex0max)
       use precision
       use global_parameters,only: mgrid
       implicit none

       integer :: mindex,nindex0max
       integer,dimension(mgrid) :: nindex
       integer,dimension(mindex,mgrid) :: indexp
       real(lk),dimension(mindex,mgrid) :: ring
     end subroutine matrix

     subroutine petsc_create(mindex)
       implicit none
       integer mindex
     end subroutine petsc_create
  end interface
#endif

  if(mring==1)then
! one ring, velocity in unit of rho0
     vring(1)=sqrt(2.0_lk)
     fring(1)=1.0_lk

  elseif(mring==2)then
! two rings good for up to k_perp=1.5
     vring(1)=0.9129713024553_lk
     vring(2)=2.233935334042_lk
     fring(1)=0.7193896325719_lk
     fring(2)=0.2806103674281_lk

  else
! three rings: exact(<0.8%) for up to k_perp=1.5
     vring(1)=0.388479356715_lk
     vring(2)=1.414213562373_lk
     vring(3)=2.647840808818_lk
     fring(1)=0.3043424333839_lk
     fring(2)=0.5833550690524_lk
     fring(3)=0.1123024975637_lk
  endif

  nindex=0
  ring=0.0_lk
  indexp=1

!OpenMP CANNOT be used here because of an gathering operation: ring(nt,ij0)=ring(nt,ij0)+rr
  do i=0,mpsi
     pdum0=psimesh(i)
     tempi=meshti(i)/(rho0*rho0)

     do j=1,mtheta(i)
        ij0=igrid(i)+j
        b=bmesh(ij0)
! grid position (pdum0,tdum0) in magnetic coordinates
        if(fielddir==1 .or. fielddir==3)then
          tdum0=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
        else
          tdum0=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
        endif
        jt=max(0,min(mtheta(i),floor(tdum0/deltat(i)+0.5)))
        ipjt=igrid(i)+jt

! 1st point having 1/4 weight is the original grid point
        nindex(ij0)=1
        indexp(1,ij0)=ij0
        ring(1,ij0)=0.25

        do kr=1,mring

! FLR from grid point and weight of 8-point for each ring. Counter-clockwise: 1,5,2,6,3,7,4,8
            rhogyro=vring(kr)*sqrt(aion*tempi*0.5/(b*qion*qion))


           do kp=1,8
              if(kp<5)then
                 ddelr=2.0_lk*pgyro(kp,ipjt)*rhogyro
                 ddelt=2.0_lk*tgyro(kp,ipjt)*rhogyro
                 wght=0.0625_lk*fring(kr)

              elseif(kp==5)then
                 ddelr=(pgyro(1,ipjt)+pgyro(2,ipjt))*rhogyro
                 ddelt=(tgyro(1,ipjt)+tgyro(2,ipjt))*rhogyro
                 wght=0.125_lk*fring(kr)
              elseif(kp==6)then
                 ddelr=(pgyro(2,ipjt)+pgyro(3,ipjt))*rhogyro
                 ddelt=(tgyro(2,ipjt)+tgyro(3,ipjt))*rhogyro
              elseif(kp==7)then
                 ddelr=(pgyro(3,ipjt)+pgyro(4,ipjt))*rhogyro
                 ddelt=(tgyro(3,ipjt)+tgyro(4,ipjt))*rhogyro
              elseif(kp==8)then
                 ddelr=(pgyro(4,ipjt)+pgyro(1,ipjt))*rhogyro
                 ddelt=(tgyro(4,ipjt)+tgyro(1,ipjt))*rhogyro
              endif


!!!!!!!!!!!Second order correction by Zhixuan!!!!!!!!!!
             icorrection=1
             if(icorrection==1)then
               if(kp<5)then
                  ddelr=ddelr+pgyro2(kp,ipjt)*4.0_lk*rhogyro*rhogyro
               elseif(kp>4)then
                  ddelr=ddelr+pgyro2(1,ipjt)*(1.0_lk+1.0_lk)*rhogyro*rhogyro
               endif

               if((kp==5).or.(kp==7))then
                  ddelt=ddelt-tgyro(2,ipjt)*tgyro(2,ipjt)*rhogyro*rhogyro
               elseif((kp==6).or.(kp==8))then
                  ddelt=ddelt+tgyro(2,ipjt)*tgyro(2,ipjt)*rhogyro*rhogyro
               endif
             endif


! particle position for each point with rho_i=2.0*vring/qion*sqrt(aion*T_i)
              pdum=pdum0+ddelr
              iflagboundary=0
              if(pdum < psi0)then
                     pdum=2.0_lk*psi0-pdum
                     iflagboundary=1
              endif
              if(pdum > psi1)then
                     pdum=2.0_lk*psi1-pdum
                     iflagboundary=1
              endif
              tdum=tdum0+ddelt


! linear interpolation in radial flux along constant theta (mangetic coordinates)
              isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
              dpx=pdum-spdpsi*real(isp-1)
! radaial spline of rg avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
              if(isp==1)dpx=sqrt(dpx)
              rg=rgpsi(1,isp)+rgpsi(2,isp)*dpx+rgpsi(3,isp)*dpx*dpx !radal coordinate rg
              ii=max(0,min(mpsi-1,floor((rg-rg0)/deltar))) ! label for inner flux surface
              wr=(pdum-psimesh(ii))/deltap(ii+1)           !weight for outer flux surface

! outer flux surface
              if(fielddir==1 .or. fielddir==3)then
                t1=modulo(tdum-(zeta1-pi2)*qtinv(ii+1),pi2)/deltat(ii+1)
              else
                t1=modulo(tdum-zeta1*qtinv(ii+1),pi2)/deltat(ii+1)
              endif
              j1=max(0,min(mtheta(ii+1)-1,floor(t1))) ! label for lower poloidal point
              wt1=t1-real(j1)                         !weight for upper poloidal point

! inner flux surface
              if(fielddir==1 .or. fielddir==3)then
                t1=modulo(tdum-(zeta1-pi2)*qtinv(ii),pi2)/deltat(ii)
              else
                t1=modulo(tdum-zeta1*qtinv(ii),pi2)/deltat(ii)
              endif
              j0=max(0,min(mtheta(ii)-1,floor(t1)))
              wt0=t1-real(j0)

! index and weight of each point
              do np=1,4
                 if(np==1)then
                    ij=igrid(ii+1)+j1+1 !upper poloidal grid on outer flux surface
                    rr=wght*wr*wt1
                 elseif(np==2)then
                    if(j1==0)j1=mtheta(ii+1) !use poloidal grid [1,mtheta]
                    ij=igrid(ii+1)+j1   !lower poloidal grid on outer flux surface
                    rr=wght*wr*(1.0_lk-wt1)
                 elseif(np==3)then
                    ij=igrid(ii)+j0+1   !upper poloidal grid on inner flux surface
                    rr=wght*(1.0_lk-wr)*wt0
                 else
                    if(j0==0)j0=mtheta(ii)
                    ij=igrid(ii)+j0     !lower poloidal grid on inner flux surface
                    rr=wght*(1.0_lk-wr)*(1.0_lk-wt0)
                 endif
                 if(iflagboundary.eq.1) rr=-1.0_lk*rr !if out of boundary then use anti-symmetric boundary condition.

                 do nt=1,nindex(ij0)
! redundant point
                    if(ij==indexp(nt,ij0))then
                       ring(nt,ij0)=ring(nt,ij0)+rr
                       goto 100
                    endif
                 enddo
! new point
                 nindex(ij0)=nindex(ij0)+1
                 nt=nindex(ij0)
                 indexp(nt,ij0)=ij
                 ring(  nt,ij0)=rr

100             continue
              enddo  !end of 4-point interpolation loop
           enddo     !end of 8-point-per-ring loop
        enddo        !end of ring loop
     enddo           !end of poloidal loop
  enddo              !end of radial loop

! check array size
  if(maxval(nindex)>mindex)then
     write(gtcout,*)'Poisson error',mype,maxval(nindex),' > ',mindex
     call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  endif

  rdum=0.0_lk
  tdum=0.0_lk
  zdum=1.0_lk
  do i=1,mgrid
     rdum=sum(ring(1:nindex(i),i))
     tdum=max(tdum,rdum)
     zdum=min(zdum,rdum)
  enddo
  if(mype==0)then
     write(gtcout,*)'poisson solver=',maxval(nindex),sum(nindex)/mgrid,tdum,zdum
     call FLUSH(gtcout)
  end if

#ifdef _PETSc
  call matrix(mindex,nindex,indexp,ring,nindex0max)

  call petsc_create(nindex0max)
#endif
end subroutine poisson_solver_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gyroinitial
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use spline_function, only: spx, spz, spb
  implicit none

  integer i,j,ij,ngyromax,igyro
  real(lk) bf,rhos,pdum0,pdum1,tdum0,tdum1,delx,delz,deltheta,delpsi

! for N-point gyro-averaging, N=1,4,8 for now.
! gyroaveraging for mesh points on regular magnetic coordinates, i.e., not on field-aligned mesh
  ngyromax=max(ngyroi,ngyrof)
  ngyromax=max(4,ngyromax)
  allocate(pgyro(ngyromax,mgrid),tgyro(ngyromax,mgrid),pgyro2(ngyromax,mgrid),tgyro2(ngyromax,mgrid))

! Starting from outermost point, counter-clockwise: 1,2,3,4.
! Calculate coeffecients here for rho_local=rho0*sqrt(2B_0/B), _0 represent quantities on axis.
! In subroutines poisson and locatei: rho_particle=rho_local*sqrt(m_i/m_proton)*sqrt(mu)/q_i.
!$omp parallel do private(i,j,ij,pdum0,pdum1,tdum0,tdum1,bf,rhos,delx,delz,deltheta,delpsi,igyro)
  do i=0,mpsi
     pdum0=psimesh(i)
     if(i==mpsi)then
        pdum1=psimesh(i-1)
     else
        pdum1=psimesh(i+1)
     endif

     do j=0,mtheta(i)
        ij=igrid(i)+j
        tdum0=deltat(i)*real(j)
        bf=spb(pdum0,tdum0)    !b-field on regular mesh of (pdum0,tdum0)
        rhos=rho0*sqrt(2.0_lk/bf) ! gyroradius of thermal proton with mu*B=T_e0

        tdum1=tdum0+deltat(i)
        tdum1=modulo(tdum1,pi2)
! distance between two poloidal grids
        delx=spx(pdum0,tdum1)-spx(pdum0,tdum0)
        delz=spz(pdum0,tdum1)-spz(pdum0,tdum0)
! change in poloidal angle
        deltheta=rhos*deltat(i)/sqrt(delx*delx+delz*delz)

! distance between two radial grids along constant theta
        delx=spx(pdum1,tdum0)-spx(pdum0,tdum0)
        delz=spz(pdum1,tdum0)-spz(pdum0,tdum0)
! change in radial flux
        delpsi=rhos*deltap(min(mpsi,i+1))/sqrt(delx*delx+delz*delz)

! first two points perpendicular to field line on poloidal surface
        pgyro(1,ij)=delpsi
        pgyro(3,ij)=0.0_lk-pgyro(1,ij)
! non-orthorgonality between psi and theta
        tgyro(1,ij)=0.0_lk
        tgyro(3,ij)=0.0_lk-tgyro(1,ij)

! the other two points tangential to field line
        tgyro(2,ij)=deltheta
        tgyro(4,ij)=0.0_lk-tgyro(2,ij)
        pgyro(2,ij)=0.0_lk
        pgyro(4,ij)=pgyro(2,ij)

! another 4-point for 8-point averaging; Starting from 1st quadrant, counter-clockwise: 5,6,7,8
        if(ngyroi==8)then
           pgyro(5,ij)=0.707_lk*pgyro(1,ij)
           tgyro(5,ij)=0.707_lk*tgyro(2,ij)
           pgyro(6,ij)=0.707_lk*pgyro(3,ij)
           tgyro(6,ij)=0.707_lk*tgyro(2,ij)
           pgyro(7,ij)=0.707_lk*pgyro(3,ij)
           tgyro(7,ij)=0.707_lk*tgyro(4,ij)
           pgyro(8,ij)=0.707_lk*pgyro(1,ij)
           tgyro(8,ij)=0.707_lk*tgyro(4,ij)
        endif

        do igyro=1,ngyromax
           pgyro2(igyro,ij)=0.5_lk*pgyro(igyro,ij)*pgyro(igyro,ij)*deltheta/delpsi+0.5_lk*tgyro(igyro,ij)*tgyro(igyro,ij)*delpsi/deltheta
           tgyro2(igyro,ij)=-pgyro(igyro,ij)*tgyro(igyro,ij)*deltheta/delpsi
        enddo

     enddo
  enddo
#ifndef GPU_UM
  !$acc update device(pgyro,tgyro,pgyro2,tgyro2)
#endif
end subroutine gyroinitial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix(mindex,nindex,indexp,ring,nindex0max)
!*Using nindex,indexp,ring generate matrices users, usera, and userp
  use precision
  use global_parameters
  use particle_array
  use field_array
  use petsc_array,only: usera,userx,userb,users,userp
  implicit none

  integer mindex,nnz,windexp,msort,nindex0max
  integer,dimension(mgrid) :: nindex,nindex0
  integer,dimension(mindex,mgrid) :: indexp,indexp0
  real(lk),dimension(mindex,mgrid) :: ring,ring0

  integer :: i,j,k,ij,ij0,ierror,mtest
  real(lk) :: wring,ddum,diagonal(0:mpsi)

!*local copy of variables
  nindex0=nindex
  indexp0=indexp
  ring0=ring

!*poloidal BC
  do i=0,mpsi
    nindex0(igrid(i))=nindex0(igrid(i)+mtheta(i))
    indexp0(:,igrid(i))=indexp0(:,igrid(i)+mtheta(i))
    indexp0(1,igrid(i))=indexp0(1,igrid(i))-mtheta(i)  !*"j=0" refer to self
    ring0(:,igrid(i))=ring0(:,igrid(i)+mtheta(i))
  enddo

  do i=0,mpsi
     if(fload>0)then
        ddum=qion*qion*(meshni(i)+meshnf(i)*afast/aion)*rho0*rho0/meshti(i)
     else
        ddum=qion*qion*meshni(i)*rho0*rho0/meshti(i)
     endif
     
     if(feload>0)then
        diagonal(i)=ddum+(1.0_lk-real(magnetic))*(meshne(i)+meshnfe(i))*rho0*rho0/meshte(i)
     else
        diagonal(i)=ddum+(1.0_lk-real(magnetic))*meshne(i)*rho0*rho0/meshte(i)
     endif

     do j=0,mtheta(i)
        ij=igrid(i)+j
        ring0(1,ij)=diagonal(i)-ring0(1,ij)*ddum
        do k=2,nindex0(ij)
           ring0(k,ij)=0.0_lk-ring0(k,ij)*ddum
        enddo
     enddo
  enddo

!*radial conditions.
  do i=0,mpsi,mpsi
     do j=0,mtheta(i)
        ij0=igrid(i)+j
        nindex0(ij0)=1
        indexp0(1,ij0)=ij0
        ring0(1,ij0)=diagonal(i)
     enddo
  enddo

!*sort indexp and ring
  do i=1,mgrid
    do msort=nindex0(i)-1,1,-1
      do j=1,msort
        if(indexp0(j,i)>indexp0(j+1,i)) then
          windexp=indexp0(j,i)
          wring=ring0(j,i)
          indexp0(j,i)=indexp0(j+1,i)
          ring0(j,i)=ring0(j+1,i)
          indexp0(j+1,i)=windexp
          ring0(j+1,i)=wring
        endif
      enddo
    enddo
  enddo

!*Count nonzero and allocate
  nindex0max=0
  nnz=0
  do i=1,mgrid
    if(nindex0(i).gt.nindex0max) nindex0max=nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
    enddo
  enddo

  allocate(usera(nnz),userp(nnz),userb(0:mgrid-1),users(-1:mgrid-1),&
    userx(0:mgrid-1),stat=mtest)
  if (mtest /= 0) then
    write(0,*)mype,'matrix: Cannot allocate userX'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

!*Make matrix usera,userp
  nnz=0
  users(-1)=0
  do i=1,mgrid
    users(i-1)=users(i-2)+nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
      usera(nnz)=ring0(j,i)
      userp(nnz)=indexp0(j,i)-1
    enddo
  enddo
end subroutine matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lapmat_initial
!*generate laplacian matrix
  use precision
  use global_parameters,only: mgrid,mpsi,mype,mpsilow,mpsihigh,mgridlow,mgridhigh,bcond
  use field_array,only: igrid,mtheta,wtp1,jtp1,deltat,psimesh,qmesh,bmesh,&
           gupp,gupt,gutt,guzz,qtinv,zeta1,phi,fluidne,mindexlap,&
           nindexlap,indexplap,lapmat
  use equilibrium,only:gpsi,cpsi,spdpsi,spdpsi_inv,lsp
  implicit none

  integer :: ierror,mtest,mgrid2
  real(lk) :: iguzz
  integer :: i,isp,ij,j,jt,jt00,jt10,jt20,jt01,jt11,jt21,jt31,jt30,jm,jc,jp,ibound
  real(lk) :: h0,h1,wt00,wt10,wt01,wt11,wp0,wp1,tdum,tdum0,tdum1
  real(lk) :: functiong,functioni,pdum,dpx,q,b
  real(lk) :: j1sdum1,j1sdum0,j1pdum1,j1pdum0,j3sdum1,j3sdum0,j3pdum1,j3pdum0
  real(lk) :: jdpgpp(mgrid),jdtgpt(mgrid),jdpgpt(mgrid),jdtgtt(mgrid),jacobtmp(mgrid)
 !! save jdpgpp,jdtgpt,jdpgpt,jdtgtt

#ifdef _PETSc
  interface
     subroutine lapmat_pscreate(mindexlap,matrixsize)
       implicit none
       integer mindexlap,matrixsize
     end subroutine lapmat_pscreate
  end interface
#endif

  iguzz=1.0_lk
  jacobtmp=1.0_lk
!compute the jacobian and related derivatives in Laplacian
! Jacobian=(gq+I)/B^2
!$omp parallel do private(i,pdum,isp,dpx,functiong,functioni,q,j,ij,b)
  do i=0,mpsi
     pdum=psimesh(i)
     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     functiong=gpsi(1,isp)+dpx*gpsi(2,isp)+dpx*dpx*gpsi(3,isp)
     functioni=cpsi(1,isp)+dpx*cpsi(2,isp)+dpx*dpx*cpsi(3,isp)
     q=qmesh(i)
     do j=1,mtheta(i)
        ij=igrid(i)+j
        b=bmesh(ij)
        jacobtmp(ij)=(functiong*q+functioni)/(b*b)
     enddo
     jacobtmp(igrid(i))=jacobtmp(igrid(i)+mtheta(i))
  enddo

  jdpgpp=0.0_lk
  jdpgpt=0.0_lk
  jdtgpt=0.0_lk
  jdtgtt=0.0_lk

!$omp parallel do private(i,j,ij,h0,h1,wp1,wp0,tdum,jm,jc,jp,wt10,wt00,wt11,&
!$omp& wt01,jt,jt00,jt01,jt11,jt21,jt31,jt10,jt20,jt30,j1sdum1,&
!$omp& j1sdum0,j1pdum1,j1pdum0,j3sdum1,j3sdum0,j3pdum1,j3pdum0)
  do i=mpsilow,mpsihigh
    h0=1.0_lk/(psimesh(i)-psimesh(i-1))
    h1=1.0_lk/(psimesh(i+1)-psimesh(i))
    wp1=h1/(h0+h1)
    wp0=1.0_lk-wp1

    tdum=1.0_lk/deltat(i)
    jm=igrid(i-1)
    jc=igrid(i)
    jp=igrid(i+1)

!!first computer dA/dpsi^2
    do j=1,mtheta(i)
      ij=igrid(i)+j

      wt10=wtp1(2,ij)   !upper poloidal grid on inner flux surface
      wt00=1.0_lk-wt10       !lower poloidal grid on inner flux surface
      wt11=wtp1(1,ij)   !upper poloidal grid on outer flux surface
      wt01=1.0_lk-wt11       !lower poloidal grid on outer flux surface

      jt11=jtp1(1,ij)   !!lower poloidal grid on outer flux surface
      jt21=jt11+1
      jt31=jt21+1-mtheta(i+1)*((jt21-jp)/mtheta(i+1))
      jt01=jp+mod(jt11-jp-1+mtheta(i+1),mtheta(i+1))

      jt=ij+1-mtheta(i)*(j/mtheta(i))

      jt10=jtp1(2,ij)   !!lower poloidal grid on the inner flux surface
      jt20=jt10+1
      jt30=jt20+1-mtheta(i-1)*((jt20-jm)/mtheta(i-1))
      jt00=jm+mod(jt10-jm-1+mtheta(i-1),mtheta(i-1))

      j1sdum1=wt01*jacobtmp(jt11)*gupp(jt11)+wt11*jacobtmp(jt21)*gupp(jt21)
      j1sdum0=wt00*jacobtmp(jt10)*gupp(jt10)+wt10*jacobtmp(jt20)*gupp(jt20)
      j1pdum1=j1sdum1+0.5_lk*wt01*wt11*(j1sdum1-((1.0_lk+wt11)*jacobtmp(jt31)*gupp(jt31)+(1.0_lk+wt01)*jacobtmp(jt01)*gupp(jt01))/3.0_lk)
      j1pdum0=j1sdum0+0.5_lk*wt00*wt10*(j1sdum0-((1.0_lk+wt10)*jacobtmp(jt30)*gupp(jt30)+(1.0_lk+wt00)*jacobtmp(jt00)*gupp(jt00))/3.0_lk)
      jdpgpp(ij)=(wp1*h1*j1pdum1-wp0*h0*j1pdum0+(wp0*h0-wp1*h1)*jacobtmp(ij)*gupp(ij))/jacobtmp(ij)

      j3sdum1=wt01*jacobtmp(jt11)*gupt(jt11)+wt11*jacobtmp(jt21)*gupt(jt21)
      j3sdum0=wt00*jacobtmp(jt10)*gupt(jt10)+wt10*jacobtmp(jt20)*gupt(jt20)
      j3pdum1=j3sdum1+0.5_lk*wt01*wt11*(j3sdum1-((1.0_lk+wt11)*jacobtmp(jt31)*gupt(jt31)+(1.0_lk+wt01)*jacobtmp(jt01)*gupt(jt01))/3.0_lk)
      j3pdum0=j3sdum0+0.5_lk*wt00*wt10*(j3sdum0-((1.0_lk+wt10)*jacobtmp(jt30)*gupt(jt30)+(1.0_lk+wt00)*jacobtmp(jt00)*gupt(jt00))/3.0_lk)
      jdpgpt(ij)=(wp1*h1*j3pdum1-wp0*h0*j3pdum0+(wp0*h0-wp1*h1)*jacobtmp(ij)*gupt(ij))/jacobtmp(ij)

      jdtgpt(ij)=0.5_lk*tdum*(jacobtmp(jt)*gupt(jt)-jacobtmp(ij-1)*gupt(ij-1))/jacobtmp(ij)
      jdtgtt(ij)=0.5_lk*tdum*(jacobtmp(jt)*gutt(jt)-jacobtmp(ij-1)*gutt(ij-1))/jacobtmp(ij)

    enddo
    jdpgpp(igrid(i))=jdpgpp(igrid(i)+mtheta(i))
    jdpgpt(igrid(i))=jdpgpt(igrid(i)+mtheta(i))
    jdtgpt(igrid(i))=jdtgpt(igrid(i)+mtheta(i))
    jdtgtt(igrid(i))=jdtgtt(igrid(i)+mtheta(i))

  enddo

!! 11 point for 2nd order finite difference
  mgrid2=mgridhigh-mgridlow+1
  mindexlap=11
!!allocate array
  allocate(nindexlap(mgridlow:mgridhigh),indexplap(mindexlap,mgridlow:mgridhigh),&
        lapmat(mindexlap,mgridlow:mgridhigh),stat=mtest)
  if (mtest /= 0) then
    write(0,*)mype,'lapmat_initial: Cannot allocate lapmat'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
  lapmat=0.0_lk
  nindexlap=0
  indexplap=0

!New laplacian operator for the inner boundary condition
!Only 7 points needed for the point after the boundary.
!lapmat(1-4,ij) represent 2nd flux surface, lapmat(5-7,ij) represent 3rd flux surface.
    i=mpsilow

    ibound=1
    h0=1.0_lk/(psimesh(i)-psimesh(i-1))
    h1=1.0_lk/(psimesh(i+1)-psimesh(i))
    wp1=h1/(h0+h1)
    wp0=1.0_lk-wp1

    tdum1=1.0_lk/deltat(i+1)
    tdum=1.0_lk/deltat(i)
    jm=igrid(i-1)
    jc=igrid(i)
    jp=igrid(i+1)

!!first computer dA/dpsi^2
    do j=1,mtheta(i)
      ij=igrid(i)+j
      wt10=wtp1(2,ij)   !upper poloidal grid on inner flux surface
      wt00=1.0_lk-wt10       !lower poloidal grid on inner flux surface
      wt11=wtp1(1,ij)   !upper poloidal grid on outer flux surface
      wt01=1.0_lk-wt11       !lower poloidal grid on outer flux surface

      jt11=jtp1(1,ij)   !!lower poloidal grid on outer flux surface
      jt21=jt11+1
      jt31=jt21+1-mtheta(i+1)*((jt21-jp)/mtheta(i+1))
      jt01=jp+mod(jt11-jp-1+mtheta(i+1),mtheta(i+1))

      jt=ij+1-mtheta(i)*(j/mtheta(i))

      nindexlap(ij)=mindexlap-4*ibound
!! current flux surface
      indexplap(1,ij)=ij-1
      lapmat(1,ij)=tdum*((-h0*wp0 + h1*wp1)*gupt(ij) + tdum*(gutt(ij)+iguzz*guzz(ij)/qmesh(i)/qmesh(i)))-&
           0.5_lk*tdum*(jdpgpt(ij) + jdtgtt(ij)) + 2.0_lk*gupt(ij)*(wp0*h0*tdum)*real(bcond)
      indexplap(2,ij)=ij
      lapmat(2,ij)=-2.0_lk*(h0*h1*gupp(ij) + tdum*tdum*(gutt(ij)+iguzz*guzz(ij)/qmesh(i)/qmesh(i)))+&
           (h0*wp0*(1.0_lk-2.0*real(bcond))-h1*wp1)*(jdpgpp(ij) + jdtgpt(ij)) +4.0_lk*gupp(ij)*wp0*h0*h1*real(bcond)
      indexplap(3,ij)=jt
      lapmat(3,ij)=tdum*((h0*wp0 - h1*wp1)*gupt(ij) + tdum*(gutt(ij)+iguzz*guzz(ij)/qmesh(i)/qmesh(i)))+&
           0.5_lk*tdum*(jdpgpt(ij) + jdtgtt(ij)) - 2.0*gupt(ij)*(wp0*h0*tdum)*real(bcond)

!!outer flux surface
      indexplap(4,ij)=jt01
      lapmat(4,ij)=-h1*wp1*wt01*(h0*(1.0_lk + wt01)*wt11*gupp(ij)/3.0_lk)*(1.0_lk-real(bcond)*wp0/wp1)-&
            h1*wp1*wt01*tdum1*gupt(ij)*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))-&
           (h1*wp1*wt01*(1.0_lk + wt01)*wt11*(jdpgpp(ij) + jdtgpt(ij)))/6.0_lk*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))
      indexplap(5,ij)=jt11
      lapmat(5,ij)=h1*wp1*h0*wt01*(2.0_lk + wt01*wt11)*gupp(ij)*(1.0_lk-real(bcond)*wp0/wp1)-&
           h1*wp1*tdum1*wt11*gupt(ij)*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))+&
           (h1*wp1*wt01*(2.0_lk + wt01*wt11)*(jdpgpp(ij) + jdtgpt(ij)))/2.0_lk*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))
      indexplap(6,ij)=jt21
      lapmat(6,ij)=h1*wp1*h0*wt11*(2.0_lk + wt01*wt11)*gupp(ij)*(1.0_lk-real(bcond)*wp0/wp1) +&
           h1*wp1*tdum1*wt01*gupt(ij)*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1)) +  &
           (h1*wp1*wt11*(2.0_lk + wt01*wt11)*(jdpgpp(ij) + jdtgpt(ij)))/2.0_lk*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))
      indexplap(7,ij)=jt31
      lapmat(7,ij)=-h1*wp1*wt11*(h0*wt01*(1.0_lk + wt11)*gupp(ij)/3.0_lk)*(1.0_lk-real(bcond)*wp0/wp1)+&
             h1*wp1*wt11*tdum1*gupt(ij)*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))-&
           (h1*wp1*wt01*wt11*(1.0_lk+wt11)*(jdpgpp(ij)+jdtgpt(ij)))/6.0_lk*(1.0_lk+real(bcond)*wp0*h0/(wp1*h1))
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!11
!*******
  !*2nd order derivative, in polar coordinates
  !* \partial/(\partial phi_t)(2 phi_t \partial #/\partial phi_t),
!$omp parallel do private(i,j,ij,h0,h1,wp1,wp0,tdum,tdum0,tdum1,jm,jc,jp,wt10,wt00,wt11,&
!$omp& wt01,jt,jt00,jt01,jt11,jt21,jt31,jt10,jt20,jt30,ibound)
  do i=mpsilow+1,mpsihigh

!ibound changes the nindexlap(ij) from 11 to 7 for boundary
    ibound=i/mpsihigh

    h0=1.0_lk/(psimesh(i)-psimesh(i-1))
    h1=1.0_lk/(psimesh(i+1)-psimesh(i))
    wp1=h1/(h0+h1)
    wp0=1.0_lk-wp1

   !! tdum1=deltat(i+1)/deltat(i)
   !! tdum=1.0_lk/(h1*deltat(i+1))
    tdum0=1.0_lk/deltat(i-1)
    tdum1=1.0_lk/deltat(i+1)
    tdum=1.0_lk/deltat(i)
    jm=igrid(i-1)
    jc=igrid(i)
    jp=igrid(i+1)

!!first computer dA/dpsi^2
    do j=1,mtheta(i)
      ij=igrid(i)+j

      wt10=wtp1(2,ij)   !upper poloidal grid on inner flux surface
      wt00=1.0_lk-wt10       !lower poloidal grid on inner flux surface
      wt11=wtp1(1,ij)   !upper poloidal grid on outer flux surface
      wt01=1.0_lk-wt11       !lower poloidal grid on outer flux surface

      jt11=jtp1(1,ij)   !!lower poloidal grid on outer flux surface
      jt21=jt11+1
      jt31=jt21+1-mtheta(i+1)*((jt21-jp)/mtheta(i+1))
      jt01=jp+mod(jt11-jp-1+mtheta(i+1),mtheta(i+1))

      jt=ij+1-mtheta(i)*(j/mtheta(i))

      jt10=jtp1(2,ij)   !!lower poloidal grid on the inner flux surface
      jt20=jt10+1
      jt30=jt20+1-mtheta(i-1)*((jt20-jm)/mtheta(i-1))
      jt00=jm+mod(jt10-jm-1+mtheta(i-1),mtheta(i-1))

      nindexlap(ij)=mindexlap-4*ibound
!! inner flux surface
      indexplap(1,ij)=jt00
      lapmat(1,ij)=-h0*wp0*wt00*(h1*(1.0_lk + wt00)*wt10*gupp(ij)/3.0_lk - tdum0*gupt(ij))+&
           (h0*wp0*wt00*(1.0_lk + wt00)*wt10*(jdpgpp(ij) + jdtgpt(ij)))/6.0_lk
      indexplap(2,ij)=jt10
      lapmat(2,ij)=h0*wp0*(h1*wt00*(2.0_lk + wt00*wt10)*gupp(ij) + tdum0*wt10*gupt(ij))-&
           (h0*wp0*wt00*(2.0_lk + wt00*wt10)*(jdpgpp(ij) + jdtgpt(ij)))/2.0_lk
      indexplap(3,ij)=jt20
      lapmat(3,ij)=h0*wp0*(h1*wt10*(2.0_lk + wt00*wt10)*gupp(ij) - tdum0*wt00*gupt(ij))-&
           (h0*wp0*wt10*(2.0_lk + wt00*wt10)*(jdpgpp(ij) + jdtgpt(ij)))/2.0_lk
      indexplap(4,ij)=jt30
      lapmat(4,ij)=-h0*wp0*wt10*(h1*wt00*(1.0_lk + wt10)*gupp(ij)/3.0_lk + tdum0*gupt(ij))+&
           (h0*wp0*wt00*wt10*(1.0_lk + wt10)*(jdpgpp(ij) + jdtgpt(ij)))/6.0_lk

!! current flux surface
      indexplap(5,ij)=ij-1
      lapmat(5,ij)=tdum*((-h0*wp0 + h1*wp1)*gupt(ij) + tdum*(gutt(ij)+iguzz*guzz(ij)/qmesh(i)/qmesh(i)))-&
           0.5_lk*tdum*(jdpgpt(ij) + jdtgtt(ij))
      indexplap(6,ij)=ij
      lapmat(6,ij)=-2.0_lk*(h0*h1*gupp(ij) + tdum*tdum*(gutt(ij)+iguzz*guzz(ij)/qmesh(i)/qmesh(i)))+&
           (h0*wp0 - h1*wp1)*(jdpgpp(ij) + jdtgpt(ij))
      indexplap(7,ij)=jt
      lapmat(7,ij)=tdum*((h0*wp0 - h1*wp1)*gupt(ij) + tdum*(gutt(ij)+iguzz*guzz(ij)/qmesh(i)/qmesh(i)))+&
           0.5_lk*tdum*(jdpgpt(ij) + jdtgtt(ij))

!!outer flux surface
      indexplap(8,ij)=jt01
      lapmat(8,ij)=-h1*wp1*wt01*(h0*(1.0_lk + wt01)*wt11*gupp(ij)/3.0_lk + tdum1*gupt(ij))-&
           (h1*wp1*wt01*(1.0_lk + wt01)*wt11*(jdpgpp(ij) + jdtgpt(ij)))/6.0_lk
      indexplap(9,ij)=jt11
      lapmat(9,ij)=h1*wp1*(h0*wt01*(2.0_lk + wt01*wt11)*gupp(ij) - tdum1*wt11*gupt(ij))+&
           (h1*wp1*wt01*(2.0_lk + wt01*wt11)*(jdpgpp(ij) + jdtgpt(ij)))/2.0_lk
      indexplap(10,ij)=jt21
      lapmat(10,ij)=h1*wp1*(h0*wt11*(2.0_lk + wt01*wt11)*gupp(ij) + tdum1*wt01*gupt(ij))+&
           (h1*wp1*wt11*(2.0_lk + wt01*wt11)*(jdpgpp(ij) + jdtgpt(ij)))/2.0_lk
      indexplap(11,ij)=jt31
      lapmat(11,ij)=-h1*wp1*wt11*(h0*wt01*(1.0_lk + wt11)*gupp(ij)/3.0_lk - tdum1*gupt(ij))-&
           (h1*wp1*wt01*wt11*(1.0_lk + wt11)*(jdpgpp(ij) + jdtgpt(ij)))/6.0_lk

    enddo
  enddo

!*poloidal BC
  do i=mpsilow,mpsihigh
!ibound treats shifted lapmat/indexplap indices for the radial boundary
    ibound=6
    if(i==mpsilow)ibound=2
    nindexlap(igrid(i))=nindexlap(igrid(i)+mtheta(i))
    indexplap(:,igrid(i))=indexplap(:,igrid(i)+mtheta(i))
    indexplap(ibound,igrid(i))=indexplap(ibound,igrid(i))-mtheta(i)  !*"j=0" refer to self
    lapmat(:,igrid(i))=lapmat(:,igrid(i)+mtheta(i))
  enddo

#ifdef _PETSc
  call lap2petsc

  call lapmat_pscreate(mindexlap,mgrid2)
#endif

end subroutine lapmat_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lap2petsc
!*Using nindex0,indexp0,ring0 generate matrices users, usera, and userp
  use precision
  use global_parameters,only: mgrid,mpsi,mype,rho0,magnetic,istep,gtcout,&
	    mpsilow,mpsihigh,mgridlow,mgridhigh
  use field_array,only: igrid,mtheta,mindexlap,nindexlap,indexplap,lapmat,rhom
  use particle_array,only: fload,feload,qion,aion,afast,afaste,qfaste,qelectron,meshti,meshte,meshni,meshne,meshnf,meshnfe,iload
  use petsc_array,only: lusera,lusers,luserp,luserb,luserx
  implicit none

  integer nnz,windexp,msort,nindex0max
  integer,dimension(mgridlow:mgridhigh) :: nindex0
  integer,dimension(mindexlap,mgridlow:mgridhigh) :: indexp0
  real(lk),dimension(mindexlap,mgridlow:mgridhigh) :: ring0

  integer :: ierror,mtest,mgrid2
  real(lk) :: wring,ddum(0:mpsi),diagonal(0:mpsi)
  integer :: i,ij,j,k
  integer,parameter :: idiagonal=1


!*local copy of variables
  nindex0=nindexlap
  indexp0=indexplap
  ring0=lapmat

!Normalization for gyrokinetic case
  if(iload/=9)then
! add diagonal term for poisson operator
     if(idiagonal==1)then
        do i=mpsilow,mpsihigh
            if(feload==0)then
              diagonal(i)=(1.0_lk-real(magnetic,lk))*meshne(i)*qelectron*qelectron*rho0*rho0/meshte(i)
           else
              diagonal(i)=(1.0_lk-real(magnetic,lk))*(meshne(i)*qelectron*qelectron+meshnfe(i)*qfaste*qfaste)*rho0*rho0/meshte(i)
           endif

           if (fload==0) then
             ddum(i)=qion*qion*meshni(i)*rho0*rho0/meshti(i)+diagonal(i)
           else
             ddum(i)=qion*qion*(meshni(i)+meshnf(i)*afast/aion)*rho0*rho0/meshti(i)+diagonal(i)
           endif
           do j=0,mtheta(i)
              ij=igrid(i)+j
              do k=1,nindex0(ij)
                 ring0(k,ij)=0.0_lk-ring0(k,ij)*ddum(i)*rhom(ij)*rhom(ij)
              enddo
              ring0(6,ij)=diagonal(i)+ring0(6,ij)
           enddo
        enddo
      endif
!Normalization for  fully kinetic case(None! already normalized in the source term)
   elseif(iload==9)then
      do i=mpsilow,mpsihigh
         do j=0,mtheta(i)
           ij=igrid(i)+j
           do k=1,nindex0(ij)
              ring0(k,ij)=0.0_lk-ring0(k,ij)
           enddo
         enddo
       enddo
   endif

!*sort indexp and ring
  do i=mgridlow,mgridhigh
    do msort=nindex0(i)-1,1,-1
      do j=1,msort
        if(indexp0(j,i)>indexp0(j+1,i)) then
          windexp=indexp0(j,i)
          wring=ring0(j,i)
          indexp0(j,i)=indexp0(j+1,i)
          ring0(j,i)=ring0(j+1,i)
          indexp0(j+1,i)=windexp
          ring0(j+1,i)=wring
        endif
      enddo
    enddo
  enddo

!*Count nonzero and allocate
  nindex0max=0
  nnz=0
  do i=mgridlow,mgridhigh
    if(nindex0(i).gt.nindex0max) nindex0max=nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
    enddo
  enddo
  mgrid2 = mgridhigh-mgridlow+1
  allocate(lusera(nnz),luserp(nnz),luserb(0:mgrid2-1),lusers(-1:mgrid2-1),&
    luserx(0:mgrid2-1),stat=mtest)


  if (mtest /= 0) then
    write(0,*)mype,'matrix: Cannot allocate userX'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

!*Make matrix lusera,luserp
  nnz=0
  lusers(-1)=0
  do i=mgridlow,mgridhigh
    lusers(i-mgridlow)=lusers(i-mgridlow-1)+nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
      lusera(nnz)=ring0(j,i)
      luserp(nnz)=indexp0(j,i)-mgridlow
    enddo
  enddo


end subroutine lap2petsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poisson_frc(scenario)
! Poisson solver on regular mesh in psi-zeta plane
  use global_parameters
  use field_array
  use particle_array
  use petsc_array
  use precision, only: lk
  implicit none

  integer i,j,k,ii,ij,ip,jt,kz,ipe,kpe,jtp,&
       icount,ierror,idest,isource,l,isp,maxmmode
  real(lk) wt,r,wz,zdum,tdum,ptemp,farray(0:1,mgrid),phiflux(mthetamax,mpsi),&
       ffilter(mtoroidal/2+1),xz(mtoroidal),&
       field3d(mtoroidal,mthetamax,mpsi),&
       ddum3(mpsi),ddum4(mpsi),ddum3i(mpsi),ddum4i(mpsi),temp1,temp2,&
       phi_zonal(mpsi),adum(mpsi)
  complex(lk) yz(mtoroidal/2+1),ddum5((mtoroidal/2+1),mpsi)
  character(*),intent(in) :: scenario

!!! for single n, no need to keep all other n's
  field3d=0.0_lk

!!! source term for poisson's eqn'
  if(scenario=='ES-adiabatic'.and.feload==0)then
!$omp parallel do private(i,j,ij)
    do i=1,mpsi-1
      do j=1,mtheta(i)
        ij=igrid(i)+j
        farray(1,ij)=qion*densityi(1,ij)*meshni(i)
      enddo
    enddo
  elseif(scenario=='ES-DKE')then
!$omp parallel do private(i,j,ij)
    do i=1,mpsi-1
      do j=1,mtheta(i)
        ij=igrid(i)+j
        farray(1,ij)=qion*densityi(1,ij)*meshni(i)+qfaste*densityfe(1,ij)*meshnfe(i)
      enddo
    enddo
  elseif(scenario=='ES-kinetic-hybrid')then
   !$omp parallel do private(i,j,ij)
    do i=1,mpsi-1
      do j=1,mtheta(i)
        ij=igrid(i)+j
        farray(1,ij)=phi_zero(1,ij)+meshte(i)*densitye(1,ij)/qelectron!
      enddo
    enddo
  endif

!2d theta-psi array
!$omp parallel do private(i,j,ij)
  do i=1,mpsi-1
     do j=1,mthetamax
        ij=igrid(i)+j
        phiflux(j,i)=farray(1,ij)
     enddo
  enddo

  field3d(myrank_toroidal+1,:,:)=phiflux(:,:)
  icount=mtoroidal*mthetamax*mpsi

  CALL MPI_ALLREDUCE(MPI_IN_PLACE,field3d,icount,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)

  ddum3=0.0_lk
  ddum4=0.0_lk
  ddum3i=0.0_lk
  ddum4i=0.0_lk
!!! for single n, no need to keep all other n's
  icount=mpsi
  yz=0.0_lk

!!transform from zeta-psi to n-psi
!!separate real and imaginary parts for petsc
!$omp parallel do private(i,yz,xz)
  do i=1,mpsi-1
     xz=field3d(:,myth0,i)
     call fftr1d(1,mtoroidal,1.0_lk,xz,yz,2)
     ddum3(i)=real(real(yz(nmodes(1)+1)))
     ddum3i(i)=real(aimag(yz(nmodes(1)+1)))
  enddo

! 0:fluxtube, kr-->0
! 1:semispectral, kr-->finite, partial torus 
if(solvermethod==0.and.(scenario=='ES-adiabatic'.or.scenario=='ES-DKE'))then 
  !$omp parallel do private(i)
  do i=1,mpsi-1
     ddum4(i)=spectraln(i)*ddum3(i)
     ddum4i(i)=spectraln(i)*ddum3i(i)
  enddo
 
elseif(solvermethod==1.and.scenario=='ES-adiabatic')then

!!!! constructing right-hand side !!!!!
  call laplacian2(ddum3,ddum4,icount)
  call laplacian2(ddum3i,ddum4i,icount)
!$omp parallel do private(i)
  do i=1,icount
     ddum3(i)=ddum3(i)-ddum4(i)
     ddum3i(i)=ddum3i(i)-ddum4i(i)
  enddo
!!! end constructing right-hand side !!!

!!!!!!!!!! solving poisson eqn using petsc !!!!!!!!!!!!!!!!!!!
!! petsc can only do real numbers so do !!!!!!!!!!!!!!!!!!!!!!
!! real and imaginary separately !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!
!!!solve real part
!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(i)
  do i=1,icount
     luserb2(i-1)=ddum3(i)
  enddo
#ifdef _PETSc
  call lapmat_pssolver2
#else
  if (mype==0) write(gtcout,*)'Poisson-FRC requires PETSc'
  stop
#endif

!$omp parallel do private(i)
  do i=1,icount
     ddum4(i)=luserx2(i-1)
  enddo

!!!!!!!!!!!!!!!!!!!!!!!!
!!!solve imag part
!!!!!!!!!!!!!!!!!!!!!!!!
!$omp parallel do private(i)
  do i=1,icount
     luserb2(i-1)=ddum3i(i)
  enddo
#ifdef _PETSc
  call lapmat_pssolver2
#else
  if (mype==0) write(gtcout,*)'Poisson-FRC requires PETSc'
  stop
#endif

!$omp parallel do private(i)
  do i=1,icount
     ddum4i(i)=luserx2(i-1)
  enddo
!!!!!!!!!! end solving poisson eqn using petsc !!!!!!!!!!!!!!!!
elseif(scenario=='ES-kinetic-hybrid')then
  ddum4=ddum3
  ddum4i=ddum3i
endif

!!! reconstructing back from n-psi to zeta-psi !!!
  field3d=0.0_lk
!$omp parallel do private(i,yz,xz)
  do i=1,mpsi-1
     xz=0.0_lk
     yz=0.0_lk
     yz(nmodes(1)+1)=cmplx(ddum4(i),ddum4i(i))
     call fftr1d(-1,mtoroidal,1.0_lk,xz,yz,2)
     field3d(:,myth0,i)=xz 
  enddo
  
  icount=mtoroidal*mthetamax*mpsi
  CALL MPI_ALLREDUCE(MPI_IN_PLACE,field3d,icount,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)

!$omp parallel do private(i,j)
  do i=1,mpsi-1
     do j=1,mthetamax
        phiflux(j,i)=field3d(myrank_toroidal+1,j,i)/real(mtoroidal/mthetamax)
     enddo
  enddo

  if(nfilter>=2)then
    ffilter=0.0_lk
    if(nfilter==2)then
      ffilter(mmodes+1)=1.0_lk
    endif
    yz=0.0_lk
    xz=0.0_lk
    !$omp parallel do private(i,xz,yz)
    do i=1,mpsi-1
       xz=phiflux(:,i)
       call fftr1d(1,mthetamax,1.0_lk,xz,yz,2)
       yz=yz*ffilter
       call fftr1d(-1,mthetamax,1.0_lk,xz,yz,2)
       phiflux(:,i)=xz
    enddo
  endif


!2d theta-psi array
!$omp parallel do private(i,j,ij)
  do i=1,mpsi-1
     do j=1,mthetamax
        ij=igrid(i)+j
        farray(1,ij)=phiflux(j,i)
     enddo
  enddo

  phi=farray
 !!for radially non-local solver, operator needs phi to be *rho0*rho0
 !!phi=farray*rho0*rho0
  
 !!in FRC version, is set to only poloidal smooth
  if(scenario=='ES-adiabatic'.or.scenario=='ES-DKE')call SMOOTH(phi)

 !! radially local, phi should be same on all flux surfaces
  !$omp parallel do private(j,i,ij)
     do j=1,mtheta(mpsi)
           i=igrid(1)+j
           ij=igrid(mpsi)+j
           phi(1,ij)=phi(1,i)
     enddo
  !$omp parallel do private(j,i,ij)
     do j=1,mtheta(0)
           i=igrid(1)+j
           ij=igrid(0)+j
           phi(1,ij)=phi(1,i)
     enddo
  call periodicity(phi)

!!save phi adiabatic to be corrected in hybrid loop
  if(scenario=='ES-adiabatic')phi_zero=phi
  
end subroutine poisson_frc
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine poisson_frc_initial
  use global_parameters
  use field_array
  use particle_array
  use petsc_array
  implicit none

  integer i,j,ii,jt,ij,mtest,icount,isp,msort,windex,nindex0max,nnz,ierr,&
          boundcond
  real(lk) wt,r,zdum,tdum,h0,h1,wring,rm1,rm2,temp1,temp2
  integer,dimension(:),allocatable :: nindexlap2loc
  integer,dimension(:,:),allocatable :: indexlap2loc
  real(lk),dimension(:),allocatable :: cdti,cdte,gpp2,gzz2
  real(lk),dimension(:,:),allocatable :: lapmat2loc

  !distribute mtheta across toroidal pe's
  myth0=ceiling((myrank_toroidal+1)/real(mtoroidal/mthetamax))

  icount=mpsi
  mindexlap2=3

  allocate(gpp2(mpsi),gzz2(mpsi),cdti(mpsi),cdte(mpsi),spectraln(mpsi),STAT=mtest)
  allocate(lapmat2(mindexlap2,icount),indexlap2(mindexlap2,icount),nindexlap2(icount),&
           lapmat2loc(mindexlap2,icount),indexlap2loc(mindexlap2,icount),nindexlap2loc(icount),STAT=mtest)

!!!!!! define laplacian operator matrix !!!!!!!
  gpp2=0.0_lk
  gzz2=0.0_lk
  cdti=0.0_lk
  cdte=0.0_lk
  spectraln=0.0_lk
  !$omp parallel do private(i,ij)
  do i=1,mpsi
     cdti(i)=qion*qion*meshni(i)/meshti(i)
     cdte(i)=qelectron*qelectron*meshne(i)/meshte(i)
     ij=igrid(i)+myth0
     gpp2(i)=gupp(ij)*rhom(ij)*rhom(ij) !!! includes rho_i^2 in the laplacian
     gzz2(i)=guzz(ij)*rhom(ij)*rhom(ij)
  enddo

  if(feload>0)cdte=0.0_lk

!! need to set solvermethod here and in poisson_frc(scenario)
!! 0:fluxtube: kr-->0
!! 1:semispectral: kr-->finite
solvermethod=0

!! need to set here and in charge/chargee
!! 0:free, can be set to zero boundary by nbound=0
!! 1:periodic, mpsi=1
boundcond=0

if(solvermethod==0)then
  if (mype==0) write(gtcout,*)'Solved in fluxtube geometry'
  !$omp parallel do private(i,temp1,temp2)
  do i=1,mpsi-1
     temp1=-gzz2(i)*real(nmode)*real(nmode)
     temp2=((cdti(i)+cdte(i))*(-1.0_lk)*temp1)
     spectraln(i)=(1.0_lk-temp1)/(temp2+cdte(i))
     if(mype==0)print *, spectraln(i)
  enddo
elseif(solvermethod==1)then
!!! psi-zeta Laplacian matrix initialization for theta belonging to individual PE (thpe)
  lapmat2=0.0_lk
  indexlap2=0.0_lk
  nindexlap2=0
  !$omp parallel do private(i,j,ii,h0,h1)
  do i=2,mpsi-1
     h0=psimesh(i)-psimesh(i-1)
     h1=psimesh(i+1)-psimesh(i)
     
     nindexlap2(i)=3
     indexlap2(1,i)=i-1
     indexlap2(2,i)=i
     indexlap2(3,i)=i+1

     lapmat2(1,i)=gpp2(i)*2.0_lk/(h0*(h1+h0))
     lapmat2(2,i)=-gpp2(i)*2.0_lk/(h0*h1)-gzz2(i)*real(nmode*nmode*nmodes(1)*nmodes(1))
     lapmat2(3,i)=gpp2(i)*2.0_lk/(h1*(h1+h0))
  enddo

  if(boundcond==0)then
    !!! inner boundary
     i=1
     h0=psimesh(i+1)-psimesh(i)
     h1=psimesh(i+2)-psimesh(i+1)
    
     nindexlap2(i)=3
     indexlap2(1,i)=i
     indexlap2(2,i)=i+1
     indexlap2(3,i)=i+2

     lapmat2(1,i)=gpp2(i)*h1/(h0*h0*(h0+h1))-gzz2(i)*real(nmode*nmode*nmodes(1)*nmodes(1))
     lapmat2(2,i)=-gpp2(i)/(h0*h0)
     lapmat2(3,i)=gpp2(i)/(h0*(h0+h1))
    
  !!! outer boundary
     i=mpsi
     h0=psimesh(i)-psimesh(i-1)
     h1=psimesh(i-1)-psimesh(i-2)
      
     nindexlap2(i)=3
     indexlap2(1,i)=i
     indexlap2(2,i)=i-1
     indexlap2(3,i)=i-2

     lapmat2(1,i)=gpp2(i)*h1/(h0*h0*(h0+h1))-gzz2(i)*real(nmode*nmode*nmodes(1)*nmodes(1))
     lapmat2(2,i)=-gpp2(i)/(h0*h0)
     lapmat2(3,i)=gpp2(i)/(h0*(h0+h1))
  elseif(boundcond==1)then
     h0=psimesh(mpsi)-psimesh(mpsi-1)
     h1=psimesh(2)-psimesh(1)

     i=1
     nindexlap2(i)=3
     indexlap2(1,i)=mpsi-1
     indexlap2(2,i)=1
     indexlap2(3,i)=2

     lapmat2(1,i)=gpp2(i)*2.0_lk/(h0*(h1+h0))
     lapmat2(2,i)=-gpp2(i)*2.0_lk/(h0*h1)-gzz2(i)*real(nmode*nmode*nmodes(1)*nmodes(1))
     lapmat2(3,i)=gpp2(i)*2.0_lk/(h1*(h1+h0))

     i=mpsi
     nindexlap2(i)=3
     indexlap2(1,i)=mpsi-1
     indexlap2(2,i)=mpsi
     indexlap2(3,i)=2

     lapmat2(1,i)=gpp2(i)*2.0_lk/(h0*(h1+h0))
     lapmat2(2,i)=-gpp2(i)*2.0_lk/(h0*h1)-gzz2(i)*real(nmode*nmode*nmodes(1)*nmodes(1))
     lapmat2(3,i)=gpp2(i)*2.0_lk/(h1*(h1+h0))
  endif


!!!!!! end defining laplacian operator matrix !!!!!!!
!!!!!! defining petsc matrix !!!!!!!

!local copy of variables
  nindexlap2loc=nindexlap2
  indexlap2loc=indexlap2
  lapmat2loc=lapmat2

!!!! left hand side operator
!!! multiply operator by rho^2, than the resulting phi muts be multiplied too
  do i=1,mpsi
     ii=i
     do jt=1,nindexlap2(i)
        !lapmat2loc(jt,ii)=rho0*rho0*(cdti(i)-cdte(i))*lapmat2loc(jt,ii)
        lapmat2loc(jt,ii)=rho0*rho0*(cdti(i)+cdte(i))*lapmat2loc(jt,ii)
     enddo
     !lapmat2loc(2,ii)=lapmat2loc(2,ii)-rho0*rho0*(-cdte(i))
     lapmat2loc(2,ii)=lapmat2loc(2,ii)-rho0*rho0*cdte(i)
  enddo

!*sort indexlap and lapmat
  do i=1,icount
    do msort=nindexlap2loc(i)-1,1,-1
      do j=1,msort
        if(indexlap2loc(j,i)>indexlap2loc(j+1,i)) then
          windex=indexlap2loc(j,i)
          wring=lapmat2loc(j,i)
          indexlap2loc(j,i)=indexlap2loc(j+1,i)
          lapmat2loc(j,i)=lapmat2loc(j+1,i)
          indexlap2loc(j+1,i)=windex
          lapmat2loc(j+1,i)=wring
        endif
      enddo
    enddo
  enddo
!*Count nonzero and allocate
  nindex0max=0
  nnz=0
  do i=1,icount
    if(nindexlap2loc(i).gt.nindex0max) nindex0max=nindexlap2loc(i)
    do j=1,nindexlap2loc(i)
      nnz=nnz+1
    enddo
  enddo
  allocate(lusera2(nnz),luserp2(nnz),luserb2(0:icount-1),lusers2(-1:icount-1),&
    luserx2(0:icount-1),stat=mtest)
  if (mtest /= 0) then
    write(0,*)mype,'Cannot allocate userX2'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierr)
  endif

!*Make matrix lusera2,luserp2
  nnz=0
  lusers2(-1)=0
  do i=1,icount
    lusers2(i-1)=lusers2(i-2)+nindexlap2loc(i)
    do j=1,nindexlap2loc(i)
      nnz=nnz+1
      lusera2(nnz)=lapmat2loc(j,i)
      luserp2(nnz)=indexlap2loc(j,i)-1
    enddo
  enddo

! Solve Poisson Equation using PETSc
#ifdef _PETSc
  call lapmat_pscreate2(mindexlap2,icount)
#else
  if (mype==0) write(gtcout,*)'Poisson-FRC requires PETSc'
  stop
#endif

  if (mype==0) write(gtcout,*)'End poisson_frc_initial'
!!!!!! end defining petsc matrix !!!!!!!
  if (mype==0) write(gtcout,*)'Solved using PETSc'
endif

  end subroutine poisson_frc_initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine laplacian2(scalar,scalar_out,icount)
  use precision
  use field_array,only: lapmat2,indexlap2,nindexlap2
  implicit none

  integer :: i,j,icount
  real(lk),dimension(icount),intent(in) :: scalar
  real(lk),dimension(icount),intent(out) :: scalar_out

!$omp parallel do private(i,j)
     do i=1,icount
        scalar_out(i)=0.0_lk
        do j=1,nindexlap2(i)
           scalar_out(i)=scalar_out(i)+lapmat2(j,i)*scalar(indexlap2(j,i))
        enddo
     enddo
end subroutine laplacian2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine laplacian_initial_fem
   use precision
   use global_parameters
   use field_array
   use equilibrium
   use spline_function, only: dxdt, dzdt
   implicit none

   integer :: ierror,mtest
   integer::trii,i,j,k,ki,kj,ij_fem
   real(lk)::deltae,b(3),c(3),lmatelem,dmatelem,dxdttmp(mgrid_fem),dzdttmp(mgrid_fem),qtmp(mgrid_fem),pdum,tdum,torcurv

#ifdef __PETSc
  interface
     subroutine lapmat_pscreate(mindexlap,matrixsize)
       implicit none
       integer mindexlap,matrixsize
     end subroutine lapmat_pscreate
  end interface
  interface
     subroutine lapmat_pscreate2(mindexlap,matrixsize)
       implicit none
       integer mindexlap,matrixsize
     end subroutine lapmat_pscreate2
  end interface
#endif

   mindex_fem=9
  
   allocate(nindex_fem(mgrid_fem),indexp_fem(mindex_fem,mgrid_fem),&
        lmat(mindex_fem,mgrid_fem),dmat(mindex_fem,mgrid_fem),stat=mtest)

   do i=0,mpsi
      pdum=psimesh(i)
      do j=1,mtheta(i)
         ij_fem=igrid_fem(i)+j
         if(fielddir==1 .or. fielddir==3)then
           tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
         else
           tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
         endif
         dxdttmp(ij_fem)=dxdt(pdum,tdum)
         dzdttmp(ij_fem)=dzdt(pdum,tdum)
         qtmp(ij_fem)=qmesh(i)
      enddo
   enddo
        
   dmat=0.0_lk
   lmat=0.0_lk
   nindex_fem=0
   indexp_fem=0


!loop over all triangles and add the three matrix elements of triangle to grid position	
   do trii=1,trilength
      i=trilist(trii,1)
      j=trilist(trii,2)
      k=trilist(trii,3)

!     deltae,b,c are same as defined in Huebner's Finite elements 
      deltae=0.5*(xygrid(i,1)-xygrid(k,1))*(xygrid(j,2)-xygrid(k,2)) -&
         0.5*(xygrid(k,1)-xygrid(j,1))*(xygrid(k,2)-xygrid(i,2))
      do ki=1,3
         i=trilist(trii,ki)
         j=trilist(trii,mod(ki,3)+1)
         k=trilist(trii,mod(ki+1,3)+1)
         b(ki) =-(xygrid(k,2)-xygrid(j,2))
         c(ki) =(xygrid(k,1)-xygrid(j,1))
      enddo

!     calculate each ki and kj in 3x3 matrix goes to i and j in full matrix 
!     if position (i,j) is not in full matrix -> increase size of array and set lapmat to matelem.
!     if position (i,j) is already in full matrix -> add matelem to lapmat
      do ki=1,3
         do kj=1,3
            i=trilist(trii,ki)
            j=trilist(trii,kj)   
            if(fem==1)then 
              torcurv = xygrid(i,1)!nabla^2=d/dx(x d/dx) + x d/dy^2 
            else
              torcurv = 1.0_lk!nabla^2 = d/dx^2 +d/dz^2  
            endif
            lmatelem= -0.25*(b(ki)*b(kj)+c(ki)*c(kj))/deltae*torcurv
!           diagonal elements are twice size of off diagonal elements
            if(ki==kj)then 
               dmatelem=deltae/6.0*torcurv
            else
               dmatelem=deltae/12.0*torcurv
            endif
!           if the matrix dmat/lmat already has data from previous triangle, add to matrix
            where(indexp_fem(:,i)==j)
               lmat(:,i)=lmat(:,i)+lmatelem
               dmat(:,i)=dmat(:,i)+dmatelem
            end where
!           if matrix location does not exist yet, create new point
            if(Count(indexp_fem(:,i)==j)==0)then
               nindex_fem(i)=nindex_fem(i)+1
               indexp_fem(nindex_fem(i),i)=j
               lmat(nindex_fem(i),i)=lmatelem
               dmat(nindex_fem(i),i)=dmatelem
            endif
         enddo
      enddo
   enddo

#ifdef _PETSc
  call lap2petsc_fem
   
  call lapmat_pscreate(mindex_fem,mgrid_fem)
  call lapmat_pscreate2(mindex_fem,mgrid_fem)
#endif


end subroutine laplacian_initial_fem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine lap2petsc_fem
!*Using nindex0,indexp0,ring0 generate matrices users, usera, and userp
  use precision
  use global_parameters,only: mgrid_fem,mpsi,mype,rho0,magnetic,istep,gtcout,psi0
  use field_array,only: igrid,igrid_fem,mtheta,mindex_fem,nindex_fem,indexp_fem,lmat,dmat,rhom
  use particle_array,only: fload,qion,aion,afast,qelectron,meshti,meshte,meshni,meshne,meshnf,iload
  use petsc_array,only: lusera,lusers,luserp,luserb,luserx,lusera2,lusers2,luserp2,luserb2,luserx2
  implicit none

  integer nnz,windexp,msort,nindex0max
  integer,dimension(mgrid_fem) :: nindex0
  integer,dimension(mindex_fem,mgrid_fem) :: indexp0
  real(lk),dimension(mindex_fem,mgrid_fem) :: ring0,ring0d

  integer :: ierror,mtest
  real(lk) :: wring,ddum(0:mpsi),diagonal(0:mpsi),zdum
  integer :: i,ij,j,k,ij_fem
  integer,parameter :: idiagonal=1



! igrid(i)    : represent grid points including poloidal BC repeated points
! igrid_fem(i): excludes poloidal BC repeated points


!*local copy of variables
  nindex0=nindex_fem
  indexp0=indexp_fem
  ring0=lmat
  ring0d=dmat

!Normalization for gyrokinetic case
  if(iload/=9)then
! add diagonal term for poisson operator
     if(idiagonal==1)then
        do i=0,mpsi
           diagonal(i)=(1.0_lk-real(magnetic,lk))*meshne(i)*qelectron*qelectron*rho0*rho0/meshte(i)
           if (fload==0) then
             ddum(i)=qion*qion*meshni(i)*rho0*rho0/meshti(i)+diagonal(i)
           else
             ddum(i)=qion*qion*(meshni(i)+meshnf(i)*afast/aion)*rho0*rho0/meshti(i)+diagonal(i)
           endif
           do j=1,mtheta(i)
              ij_fem=igrid_fem(i)+j
              ij=igrid(i)+j
              do k=1,nindex0(ij_fem)
                 ring0(k,ij_fem)=0.0_lk-ring0(k,ij_fem)*ddum(i)*rhom(ij)*rhom(ij)
                 if(indexp0(k,ij_fem)==ij_fem)ring0(k,ij_fem)=diagonal(i)+ring0(k,ij_fem)
              enddo
           enddo
        enddo
      endif
!Normalization for fully kinetic case(None! already normalized in the source term)
   elseif(iload==9)then
      do i=0,mpsi
        do j=1,mtheta(i)
           ij=igrid_fem(i)+j
           do k=1,nindex0(ij)
              ring0(k,ij)=0.0_lk-ring0(k,ij)
           enddo
         enddo
       enddo
   endif

!*sort indexp0 and ring0 for lmat matrix
  do i=1,mgrid_fem
    do msort=nindex0(i)-1,1,-1
      do j=1,msort
        if(indexp0(j,i)>indexp0(j+1,i)) then
          windexp=indexp0(j,i)
          wring=ring0(j,i)
          indexp0(j,i)=indexp0(j+1,i)
          ring0(j,i)=ring0(j+1,i)
          indexp0(j+1,i)=windexp
          ring0(j+1,i)=wring
        endif
      enddo
    enddo
  enddo
  
  indexp0=indexp_fem
!*sort again for indexp and ring0d
  do i=1,mgrid_fem
    do msort=nindex0(i)-1,1,-1
      do j=1,msort
        if(indexp0(j,i)>indexp0(j+1,i)) then
          windexp=indexp0(j,i)
          wring=ring0d(j,i)
          indexp0(j,i)=indexp0(j+1,i)
          ring0d(j,i)=ring0d(j+1,i)
          indexp0(j+1,i)=windexp
          ring0d(j+1,i)=wring
        endif
      enddo
    enddo
  enddo
  
!Applying BC to lapmat and dmat matrices
    do j=1,mtheta(mpsi-1)
      ij=igrid_fem(mpsi-1)+j
      do k=1,nindex0(ij)
        if(indexp0(k,ij)>igrid_fem(mpsi))then
          ring0d(k,ij)=0.0_lk
          ring0(k,ij)=0.0_lk
        endif
      enddo
    enddo
    do j=1,mtheta(mpsi)
      ij=igrid_fem(mpsi)+j
      do k=1,nindex0(ij)
          ring0d(k,ij)=0.0_lk
          ring0(k,ij)=0.0_lk
        if(indexp0(k,ij)==ij)then
          ring0d(k,ij)=1.0_lk
          ring0(k,ij)=1.0_lk
        endif
      enddo
    enddo
!Set inner BC if axis is not included in simulation  
  if(psi0>0.0_lk)then
    do j=1,mtheta(0)
      ij=igrid_fem(0)+j
      do k=1,nindex0(ij)
        ring0d(k,ij)=0.0_lk
        ring0(k,ij)=0.0_lk
        if(indexp0(k,ij)==ij)then
          ring0d(k,ij)=1.0_lk
          ring0(k,ij)=1.0_lk
        endif
      enddo
    enddo
    do j=1,mtheta(1)
      ij=igrid_fem(1)+j
      do k=1,nindex0(ij)
        if(indexp0(k,ij)<=igrid_fem(1))then
          ring0d(k,ij)=0.0_lk
          ring0(k,ij)=0.0_lk
        endif
      enddo
    enddo
  endif  

!*Count nonzero and allocate
  nindex0max=0
  nnz=0
  do i=1,mgrid_fem
    if(nindex0(i).gt.nindex0max) nindex0max=nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
    enddo
  enddo
  
  allocate(lusera(nnz),luserp(nnz),luserb(0:mgrid_fem-1),lusers(-1:mgrid_fem-1),luserx(0:mgrid_fem-1),&
  lusera2(nnz),luserp2(nnz),luserb2(0:mgrid_fem-1),lusers2(-1:mgrid_fem-1),luserx2(0:mgrid_fem-1),stat=mtest)

  if (mtest /= 0) then
    write(0,*)mype,'matrix: Cannot allocate userX'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

!*Make matrix lusera,luserp
  nnz=0
  lusers(-1)=0
  do i=1,mgrid_fem
    lusers(i-1)=lusers(i-2)+nindex0(i)
    do j=1,nindex0(i)
      nnz=nnz+1
      lusera(nnz)=ring0(j,i)
      lusera2(nnz)=ring0d(j,i)
      luserp(nnz)=indexp0(j,i)-1
    enddo
  enddo
! lusers and luserp are same for dmat and lmat
  lusers2=lusers
  luserp2=luserp 
   
  
end subroutine lap2petsc_fem

subroutine laplacian_fem(poissamp,scalar,scalar_out)
   use precision
   use global_parameters,only: mgrid,mgrid_fem,mpsi,mype,psi0
   use field_array,only: igrid,igrid_fem,mtheta,mindex_fem,nindex_fem,indexp_fem,lmat,dmat
   use petsc_array,only: luserb,luserx,luserp,luserb2,luserx2,luserp2
   implicit none

   integer::i,j,ij,k,ij_fem,nnz
   real(lk),dimension(mgrid_fem)::scalar_fem,rowsum
   integer,intent(in)::poissamp
   real(lk),dimension(mgrid),intent(in) :: scalar
   real(lk),dimension(mgrid),intent(out) :: scalar_out

   scalar_out=0.0_lk

#ifdef _PETSc
! igrid(i)    : represent grid points including poloidal BC repeated points
! igrid_fem(i): excludes poloidal BC repeated points

!convert scalar of size to scalar_fem of mgrid_fem size(remove redundant poloidal BC)
!$omp parallel do private(i,j,ij_fem,ij)
   do i=0,mpsi
       do j=1,mtheta(i)
         ij_fem=igrid_fem(i)+j
         ij=igrid(i)+j
         luserb(ij_fem-1)=0.0_lk
         luserb2(ij_fem-1)=0.0_lk
         scalar_fem(ij_fem)=scalar(ij)
       enddo
   enddo
!  solve ampere's law
   if(poissamp==1)then    
!     matrix multiply lmat and input scalar
!$omp parallel do private(i,j,ij,k)
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid_fem(i)+j
            do k=1,nindex_fem(ij) 
              luserb2(ij-1)=luserb2(ij-1)+lmat(k,ij)*scalar_fem(indexp_fem(k,ij))
            enddo
         enddo
      enddo 

      call lapmat_pssolver2 ! Call petsc to solve matrix

!     apply BCs
      if(psi0 > 0.0_lk)then
         do i=0,mpsi,mpsi
            do j=1,mtheta(i)
               ij=igrid_fem(i)+j
               do k=1,nindex_fem(ij) 
                 luserx2(ij-1)=0.0_lk
               enddo
            enddo
         enddo
!     no inner BC if psi0=0.0 
      else
         i=mpsi
         do j=1,mtheta(i)
            ij=igrid_fem(i)+j
            do k=1,nindex_fem(ij) 
              luserx2(ij-1)=0.0_lk
            enddo
         enddo
      endif
!     Convert back to array of mgrid size
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            ij_fem=igrid_fem(i)+j
            scalar_out(ij)=luserx2(ij_fem-1)
         enddo
      enddo
!  solve poisson's equation
   elseif(poissamp==0)then
!     matrix multiply dmat and input scalar
!$omp parallel do private(i,j,ij,k)
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid_fem(i)+j
            do k=1,nindex_fem(ij) 
               luserb(ij-1)=luserb(ij-1)+dmat(k,ij)*scalar_fem(indexp_fem(k,ij))
            enddo
         enddo
      enddo

      call lapmat_pssolver!call petsc to solve matrix equation

!     Convert back to array of mgrid size
!$omp parallel do private(i,j,ij,ij_fem)
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            ij_fem=igrid_fem(i)+j
            scalar_out(ij)=luserx(ij_fem-1)
         enddo
      enddo 
   endif
#endif
 
end subroutine laplacian_fem
