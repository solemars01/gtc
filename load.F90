! Load or restart all particle species
subroutine load
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use spline_function, only: dxdp, dxdt, dzdp, dzdt, spx, construct_spline1d, invert_spline0
  implicit none

  integer ierror,i,j,ij
  real(lk) dmarki(0:1,mgrid),pmarkit(0:mpsi),w_initial_ion,w_initial_fast,w_initial_electron,w_initial_faste

  interface
    subroutine gkLoad(species_name,w_initial)
      use precision
      implicit none
      character(*),intent(in) :: species_name
      real(lk) w_initial
    end subroutine
  end interface

!iload=0: ideal MHD with nonuniform pressure profile (eload=3) or uniform profile (eload=1)
!iload=1: Marker uniform temperature & density: n0,te0,ti0 (two-scale expansion)   
!iload=2: Marker non-uniform temperature & uniform density: n0,te(r),ti(r), marker weight corrected
!iload=3: Physical PDF for marker with non-uniform temperature & density: n(r),te(r),ti(r)
!iload>99: non-perturbative (full-f) simulation

! initialize particle position, velocity and weight: thermal ions, fast ions, electrons and fast electrons
  w_initial_ion=1.0e-8
  w_initial_fast=1.0e-8
  w_initial_electron=1.0e-8
  w_initial_faste=1.0e-8
  call gkLoad("thermal-ion",w_initial_ion)
  if(fload>0) call gkLoad("fast-ion",w_initial_fast)
  if(nhybrid>0)call gkLoad("thermal-electron",w_initial_electron)
  if(feload>0) call gkLoad("fast-electron",w_initial_faste)
  if(irun>0)call restart_load

! Tag particle for diagnostics
  if(track_particles==1 .and. irun==0)CALL TAG
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! For non-perturbative simulation:
! step 1: irun=0, nonlinear=0 for Monte-Carlo run to obtain PDF as constants of motion
! step 2: irun=1 Monte-Carlo run to calculate time-averaged # of marker particles per cell
! step 3: nonlinear=1: turbulence run by using time-averaged # of marker particles per cell 
! irun>1: continue turbulence run without changing markeri and pmarki

  if(irun==1 .and. nonlinear==1 .and. iload>99)then
! smoothing # of marker particles per cell
     dmarki(1,:)=markerit
     call smooth(dmarki)
     markeri=1.0/dmarki(1,:)

! # of marker particles per annulus
     pmarkit=0.0
!$omp parallel do private(i,j,ij)
     do i=0,mpsi
        do j=1,mtheta(i)
           ij=igrid(i)+j
           pmarkit(i)=pmarkit(i)+dmarki(1,ij)
        enddo
     enddo

! global sum for # of marker particles per cell
     call MPI_ALLREDUCE(pmarkit,pmarki,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)

! reset field and particle weight for perturbative (delta-f) simulation iload >1
!     zion(5,:)=0.0
!     phi=0.0
!     zonali=0.0
!     call random_number(zion(5,1:mi))
!     zion(5,1:mi)=0.001*zion(5,1:mi)

  endif

  markerit=0.0 ! re-set time averaging

end subroutine load

subroutine restart_load
  use global_parameters
  use particle_tracking,only:ptrackedi,ptrackede
  use particle_array,only:mimax,feload,mfemax
  implicit none

  integer iorestart,ierror
  character(len=9)cdum1
  if(track_particles==1) allocate(ptrackedi(nparami,max(mimax,1)))

#ifdef _FRC
  if(track_particles==1 .and. feload>0) allocate(ptrackede(nparamfe,max(mfemax,1)))
#endif
  if(mype==0)then
!!get restart directory info
     iorestart=345
     open(iorestart,file="FileExit.out",status="old")
     read(iorestart,"(A9,i1)")cdum1,FileExit
     read(iorestart,"(A9,i5)")cdum1,irest
     close(iorestart)
     irest=irest-1
  endif
  call mpi_bcast(irest,1,mpi_integer,0,mpi_comm_world,ierror)
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
  call restart_io("read")
  return

end subroutine restart_load

subroutine gkLoad(species_name,w_initial)
  use global_parameters,only: mype,gtcout,nparami,nparame,nparamf,nhybrid,&
                              ncyclee,nparamfe,ncyclefe,irun,mgrid,mpsi
  use particle_array
  use field_array,only: nmodes,nmode
  use equilibrium,only: nipp,tipp,nepp,tepp,nfpp,tfpp,nfepp,tfepp
  implicit none

  integer ierror,mtest
  real(lk) w_initial,mpartfactor
  character(*),intent(in) :: species_name

  interface
    subroutine gkLoadParticle(meshn,tppp,nppp,zpart,zpart0,wppart,&
        wtpart0,wtpart1,density,flow,data1d,jtpart0,jtpart1,wzpart,&
        marker,markert,pmark,pflux,rdtem,zonal,zonalc,mesht,qpart,apart,&
        pload,nparam,ngyro,mp,mpmax,w_initial,ncyclep,zpart1,&
        trap,trapfracn,trapfrac,hybrid,pressurepara,pressureperp,phit,dnt,&
        phisave,dnsave)

      use precision
      use spline_function, only: construct_spline1d, invert_spline0
      implicit none
      integer pload,nparam,ngyro,mp,mpmax
      integer,optional :: ncyclep,trap,hybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(lk) qpart,apart,w_initial
      real(lk),optional :: trapfracn
      real(lk),dimension(0:) :: meshn,mesht,pmark,pflux,rdtem,zonal,zonalc
      real(lk),dimension(:) :: wzpart,marker,markert
      real(lk),dimension(:),allocatable,optional :: trapfrac
      real(lk),dimension(:,:) :: tppp,nppp
      real(lk),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
      real(lk),dimension(0:,:) :: density,flow,data1d
      real(lk),dimension(:,:),optional :: zpart1
      real(lk),dimension(0:,:),optional :: pressurepara,&
        pressureperp,phit,dnt
      real(lk),dimension(0:,:,:),optional :: phisave,dnsave
    end subroutine gkLoadParticle
  end interface
 
#ifdef _TOROIDAL3D
  mpartfactor=3
#else
  mpartfactor=1
#endif

  if(species_name=='thermal-ion')then
    ! allocate particle memory
    mimax=mi+100*ceiling(sqrt(real(mi))) !particle array upper bound
    mimax=mimax*mpartfactor
    allocate(zion(nparami,mimax),zion0(nparami,mimax),wzion(mimax),&
      wpion(ngyroi,mimax),jtion0(ngyroi,mimax),jtion1(ngyroi,mimax),&
      wtion0(ngyroi,mimax),wtion1(ngyroi,mimax),STAT=mtest)
    allocate(pmarki(0:mpsi),markeri(mgrid),densityi(0:1,mgrid),&
      flowi(0:1,mgrid),pfluxi(0:mpsi),rdtemi(0:mpsi),zonali(0:mpsi),&
      zonalci(0:mpsi),data1di(0:mpsi,mpdata1d),markerit(mgrid),STAT=mtest)

    if(mtest /= 0) then
      write(0,*)mype,'*** Cannot allocate particle: mtest=',mtest,'@',&
        __FILE__, __LINE__
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
    endif

    call gkLoadParticle(meshni,tipp,nipp,zion,zion0,wpion,wtion0,wtion1,&
      densityi,flowi,data1di,jtion0,jtion1,wzion,markeri,markerit,pmarki,&
      pfluxi,rdtemi,zonali,zonalci,meshti,qion,aion,iload,&
      nparami,ngyroi,mi,mimax,w_initial)

    if(iload==9)call gk2fk

    if(mype==0 .and. irun==0)then
      write(gtcout,*)'mi=',mi
      call FLUSH(gtcout)
    endif

#ifndef GPU_UM
    !$acc update device(zion,zion0)
#endif
  elseif(species_name=='thermal-electron')then
    ! allocate particle memory
    memax=me+100*ceiling(sqrt(real(me))) !particle array upper bound
    memax=memax*mpartfactor
    allocate(zelectron(nparame,memax),zelectron0(nparame,memax),&
      wzelectron(memax),wpelectron(1,memax),jtelectron0(1,memax),&
      jtelectron1(1,memax),wtelectron0(1,memax),wtelectron1(1,memax),&
      STAT=mtest)
    allocate(pmarke(0:mpsi),markere(mgrid),densitye(0:1,mgrid),&
      flowe(0:1,mgrid),pfluxe(0:mpsi),rdteme(0:mpsi),zonale(0:mpsi),&
      zonalce(0:mpsi),data1de(0:mpsi,mpdata1d),markeret(mgrid),STAT=mtest)
    allocate(pressureepara(0:1,mgrid),pressureeperp(0:1,mgrid),&
      phisave(0:1,mgrid,2*nhybrid),phit(0:1,mgrid),&
      dnesave(0:1,mgrid,2*nhybrid),dnet(0:1,mgrid),STAT=mtest)

    if(ncyclee>0) allocate(zelectron1(nparame,memax),STAT=mtest)

    if(mtest /= 0) then
      write(0,*)mype,'*** Cannot allocate particle: mtest=',mtest,'@',&
        __FILE__, __LINE__
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
    endif

    call gkLoadParticle(meshne,tepp,nepp,zelectron,zelectron0,wpelectron,&
      wtelectron0,wtelectron1,densitye,flowe,data1de,jtelectron0,jtelectron1,&
      wzelectron,markere,markeret,pmarke,pfluxe,rdteme,zonale,zonalce,meshte,&
      qelectron,aelectron,eload,nparame,1,me,memax,&
      w_initial,ncyclep=ncyclee,zpart1=zelectron1,trap=etrap,&
      trapfracn=tfracn,trapfrac=tfrac,hybrid=nhybrid,&
      pressurepara=pressureepara,pressureperp=pressureeperp,phit=phit,&
      dnt=dnet,phisave=phisave,dnsave=dnesave)
    if(mype==0 .and. irun==0)then
       write(gtcout,*)'me=',me,',   trapped fraction',tfracn
       call FLUSH(gtcout)
    endif
#ifndef GPU_UM
    !$acc update device(zelectron,zelectron0)
#endif
  elseif(species_name=='fast-ion')then
    ! allocate particle memory
    mfmax=mf+100*ceiling(sqrt(real(mf))) !particle array upper bound
    mfmax=mfmax*mpartfactor
    allocate(zfast(nparamf,mfmax),zfast0(nparamf,mfmax),wzfast(mfmax),&
      wpfast(ngyrof,mfmax),jtfast0(ngyrof,mfmax),jtfast1(ngyrof,mfmax),&
      wtfast0(ngyrof,mfmax),wtfast1(ngyrof,mfmax),STAT=mtest)
    allocate(pmarkf(0:mpsi),markerf(mgrid),densityf(0:1,mgrid),&
      flowf(0:1,mgrid),pfluxf(0:mpsi),rdtemf(0:mpsi),zonalf(0:mpsi),&
      zonalcf(0:mpsi),data1df(0:mpsi,mpdata1d),markerft(mgrid),STAT=mtest)

    if(mtest /= 0) then
      write(0,*)mype,'*** Cannot allocate particle: mtest=',mtest,'@',&
        __FILE__, __LINE__
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
    endif

    call gkLoadParticle(meshnf,tfpp,nfpp,zfast,zfast0,wpfast,wtfast0,&
      wtfast1,densityf,flowf,data1df,jtfast0,jtfast1,wzfast,markerf,markerft,&
      pmarkf,pfluxf,rdtemf,zonalf,zonalcf,meshtf,qfast,afast,fload,&
      nparamf,ngyrof,mf,mfmax,w_initial)

    if(mype==0 .and. irun==0)then
      write(gtcout,*)'mf=',mf
      call FLUSH(gtcout)
    endif

#ifndef GPU_UM
    !$acc update device(zfast,zfast0)
#endif
  elseif(species_name=='fast-electron')then
    ! allocate particle memory
    mfemax=mfe+100*ceiling(sqrt(real(mfe))) !particle array upper bound
    mfemax=mfemax*mpartfactor
    allocate(zfaste(nparamfe,mfemax),zfaste0(nparamfe,mfemax),wzfaste(mfemax),&
      wpfaste(ngyrofe,mfemax),jtfaste0(ngyrofe,mfemax),&
      jtfaste1(ngyrofe,mfemax),wtfaste0(ngyrofe,mfemax),&
      wtfaste1(ngyrofe,mfemax),STAT=mtest)
    allocate(pmarkfe(0:mpsi),markerfe(mgrid),densityfe(0:1,mgrid),&
      flowfe(0:1,mgrid),pfluxfe(0:mpsi),rdtemfe(0:mpsi),zonalfe(0:mpsi),&
      zonalcfe(0:mpsi),data1dfe(0:mpsi,mpdata1d),markerfet(mgrid),STAT=mtest)

    if(ncyclefe>0) allocate(zfaste1(nparamfe,mfemax),STAT=mtest)

    if(mtest /= 0) then
      write(0,*)mype,'*** Cannot allocate particle: mtest=',mtest,'@',&
        __FILE__, __LINE__
      call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
    endif

    call gkLoadParticle(meshnfe,tfepp,nfepp,zfaste,zfaste0,wpfaste,wtfaste0,&
      wtfaste1,densityfe,flowfe,data1dfe,jtfaste0,jtfaste1,wzfaste,markerfe,&
      markerfet,pmarkfe,pfluxfe,rdtemfe,zonalfe,zonalcfe,meshtfe,qfaste,&
      afaste,feload,nparamfe,ngyrofe,mfe,mfemax,w_initial,&
      ncyclep=ncyclefe,zpart1=zfaste1,trap=fetrap,trapfracn=fetfracn,&
      trapfrac=fetfrac)

    if(mype==0 .and. irun==0)then
      write(gtcout,*)'mfe=',mfe,',   trapped fraction',fetfracn
      call FLUSH(gtcout)
    endif

#ifndef GPU_UM
    !$acc update device(zfaste,zfaste0)
#endif
    if(ncyclefe>0)then
#ifndef GPU_UM
      !$acc update device(zfaste1)
#endif
      continue
    endif
  else
    write(*,*)'load.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine gkLoad

subroutine gkLoadParticle(meshn,tppp,nppp,zpart,zpart0,wppart,wtpart0,&
    wtpart1,density,flow,data1d,jtpart0,jtpart1,wzpart,marker,markert,pmark,&
    pflux,rdtem,zonal,zonalc,mesht,qpart,apart,pload,nparam,ngyro,mp,mpmax,&
    w_initial,ncyclep,zpart1,trap,trapfracn,trapfrac,hybrid,&
    pressurepara,pressureperp,phit,dnt,phisave,dnsave)

  use precision
  use global_parameters,only: mgrid,mype,irun,npartdom,numberpe,pi,rg0,gtcout,&
    mpsi,psi0,psi1
  use particle_array,only: mpdata1d,sd_v0,sd_vc,sd_l0,sd_widthInv
  use field_array,only: deltap,psimesh,qmesh,mtheta,igrid,bmesh,deltat,deltar,&
    deltaz,zeta0,iflux,bmesh,nmodes,nmode,solvermethod
  use equilibrium,only: lsp,spdpsi_inv,spdpsi,gpsi,cpsi,spdtheta,spdtheta_inv,&
    bsp,rgpsi,lst,ropp,xsp
  use spline_function, only: invert_spline0, construct_spline1d
  implicit none

  !declaration of the dummy arguments
  !character(*),intent(in) :: w_init_drive
  integer pload,nparam,ngyro,mp,mpmax
  integer,optional :: ncyclep,trap,hybrid
  integer,dimension(:,:) :: jtpart0,jtpart1
  real(lk) qpart,apart
  real(lk),optional :: trapfracn
  real(lk),dimension(0:) :: meshn,mesht,pmark,pflux,rdtem,zonal,zonalc
  real(lk),dimension(:) :: wzpart,marker,markert
  real(lk),dimension(:),allocatable,optional :: trapfrac
  real(lk),dimension(:,:) :: tppp,nppp
  real(lk),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
  real(lk),dimension(0:,:) :: density,flow,data1d
  real(lk),dimension(:,:),optional :: zpart1
  real(lk),dimension(0:,:),optional :: pressurepara,&
    pressureperp,phit,dnt
  real(lk),dimension(0:,:,:),optional :: phisave,dnsave
  
  !declaration of the local variables
  integer i,isp,j,jst,ii,ij,m,mtest,ierror,ir,marknum(0:mpsi),&
    icount,jt,mmin,mmax,mptmp,perturbmarknum(0:mpsi),init
  real(lk) :: c0=2.515517,c1=0.802853,c2=0.010328,d1=1.432788,d2=0.189269,d3=0.001308
  real(lk) w_initial,energy,rg,delr,pdum,functiong,functioni,dp1,vth,&
    b,bm,q,upara,eperp,cmratio,zdum,tdum,umax,delpsi,delp(mpsi),&
    dtx,dt2,pmin,pmax,volt(3,lst),tvol(3,lst),spdt_inv,dpx,dp2,fullf,&
    xdum,dden(0:mpsi),ddum(0:mpsi),majorr,upara0,&
    markload(0:mpsi),tp_inv,np_inv,perturbpmark(0:mpsi),randvarv,&
    randvarf,randvart,sd_ftmp,fmax,zpartsign
  logical subcycle,trapped,fkhybrid
  mtest=0
  subcycle=.false.
  if(present(ncyclep))then
    if(ncyclep>0) subcycle=.true.
  endif

  trapped=.false.
  if(present(trap))then
    if(trap>0) trapped=.true.
  endif

  fkhybrid=.false.
  if(present(hybrid))then
    if(hybrid>0) fkhybrid=.true.
  endif

  if(trapped)then
    allocate(trapfrac(0:mpsi),STAT=mtest)
    dden=0.0_lk
    ddum=0.0_lk
    icount=mpsi+1
    trapfracn=1.0_lk
    trapfrac=1.0_lk
  endif

  if(mtest /= 0) then
    write(0,*)mype,'*** Cannot allocate particle: mtest=',mtest,'@',&
      __FILE__, __LINE__
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  pflux=0.0
  rdtem=0.0
  zonal=0.0
  zonalc=0.0
  pmark=0.0
  marker=0.0
  markert=0.0
  density=0.0
  flow=0.0
  data1d=0.0

  zpart=0.0
  zpart0=0.0
  if(subcycle) zpart1=0.0
  wzpart=0.0
  wppart=0.0
  jtpart0=0
  jtpart1=0
  wtpart0=0.0
  wtpart1=0.0

  if(fkhybrid)then
    pressurepara=0.0
    pressureperp=0.0
    phisave=0.0
    dnsave=0.0
  endif

  if(irun>0 .or. pload==0)return

  delp=1.0/deltap

  ! number of marker per grid, Jacobian=(gq+I)/B^2
  !$omp parallel do private(i,delpsi,pdum,isp,dpx,functiong,functioni,q,j,ij,b)
  do i=0,mpsi
    delpsi=0.5*(psimesh(min(mpsi,i+1))-psimesh(max(0,i-1)))
    pdum=psimesh(i)
    isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
    dpx=pdum-spdpsi*real(isp-1)
    functiong=gpsi(1,isp)+dpx*gpsi(2,isp)+dpx*dpx*gpsi(3,isp)
    functioni=cpsi(1,isp)+dpx*cpsi(2,isp)+dpx*dpx*cpsi(3,isp)
    q=qmesh(i)
    do j=1,mtheta(i)
      ij=igrid(i)+j
      b=bmesh(ij)
      marker(ij)=meshn(i)*delpsi*deltat(i)*(functiong*q+functioni)/(b*b)
    enddo
  enddo

  ! # of marker of per annulus
  ! cannot use OpenMP in this loop
  pmark=0.0
  do i=0,mpsi
    do j=1,mtheta(i)
      ij=igrid(i)+j
      pmark(i)=pmark(i)+marker(ij)
    enddo
  enddo

  ! normalization to marker #
  pdum=sum(pmark)
  tdum=real(mp)*real(npartdom)
  !$omp parallel do private(i,j,ij)
  do i=0,mpsi
    do j=1,mtheta(i)
      ij=igrid(i)+j
      marker(ij)=marker(ij)*tdum/pdum ! # of marker per cell
    enddo
    marker(igrid(i))=marker(igrid(i)+mtheta(i))
  enddo

#ifdef _FRC
  if(solvermethod==0)then
    do j=1,mtheta(1)
       ij=igrid(1)+j
       marker(ij)=marker(ij)+marker(igrid(0)+j)+marker(igrid(mpsi)+j)
    enddo
  endif
#endif

  !$omp parallel do private(i,j,ij)
  do i=0,mpsi
    do j=1,mtheta(i)
      ij=igrid(i)+j
      marker(ij)=1.0/marker(ij) !to avoid divide operation
    enddo
    marker(igrid(i))=marker(igrid(i)+mtheta(i))
  enddo

  marknum=int(real(mp)*pmark/sum(pmark)) !# of marker per annulus per MPI
  pmark=real(marknum*numberpe)  !# of marker per annulus

  if(pload==2)then !for loading uniform density
    markload=pmark/meshn
    !# of marker per annulus per MPI
    marknum=int(real(mp)*markload/sum(markload)) 
    pdum=sum(marknum*meshn)/sum(marknum)
    pmark=pmark*pdum
    marker=marker/pdum
  endif

  !initial density perturbation for fullk simulation
  init=0
  if(pload==9 .and. init==1)then
    do i=0,mpsi
      xdum=real(i)/real(mpsi)*PI
      xdum=1.0+1e-2*sin(xdum)
      perturbpmark(i)=pmark(i)*xdum
    enddo
    perturbmarknum=int(real(mp)*perturbpmark/sum(pmark))
    perturbpmark=real(perturbmarknum*numberpe)
    mp=sum(perturbmarknum)
  else
    perturbmarknum=marknum
  endif

  mp=sum(perturbmarknum)
  delr=1.0/deltar
  spdt_inv=1.0/spdtheta

  umax=4.0
  ! random # uniformly distributed between 0 and 1
  call random_number(zpart(1,1:mp))  !psi
  call random_number(zpart(2,1:mp))  !theta
  call random_number(zpart(3,1:mp))  !zeta
  call random_number(zpart(4,1:mp))  !rho_para
  call random_number(zpart(5,1:mp))  !marker weight
  call random_number(zpart(6,1:mp))  !sqrt mu

  mmin=0
  mmax=0
  do ir=0,mpsi
    ! range of psi
    pmin=psimesh(max(0,ir-1))
    pmax=psimesh(min(mpsi,ir+1))

    ! range of particle #
    mmin=mmax+1
    mmax=mmin-1+perturbmarknum(ir)

    ! load particles in radial direction: linear in psi within a grid cell
    !$omp parallel do private(m,pdum)
    do m=mmin,mmax
      pdum=zpart(1,m)
      zpart(1,m)=(1.0-pdum)*pmin+pdum*pmax
    enddo

    ! Jacobian as function of theta on each flux surface & its inversion
    volt=0.0
    tvol=0.0
    ! cannot use OpenMP in this loop
    do j=2,lst
      pdum=psimesh(ir)
      tdum=spdtheta*(real(j-1)-0.5)

      isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
      dpx=pdum-spdpsi*real(isp-1)
      if(isp==1)dpx=sqrt(dpx)
      dp2=dpx*dpx

      jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
      dtx=tdum-spdtheta*real(jst-1)
      dt2=dtx*dtx
      b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
           (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
           (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

      volt(1,j)=volt(1,j-1)+1.0/(b*b)
    enddo
    volt=2.0*pi*volt/volt(1,lst)
    call construct_spline1d(lst,spdtheta,volt,0)
    call invert_spline0(0,lst,spdtheta,spdtheta,volt,tvol)

    ! load particles in theta: uniform density in theta
    !$omp parallel do private(m,tdum,jt,dtx)
    do m=mmin,mmax
      tdum=2.0*pi*zpart(2,m)
      jt=max(1,min(lst-1,ceiling(tdum*spdt_inv)))
      dtx=tdum-spdtheta*real(jt-1)
      zpart(2,m)=tvol(1,jt)+dtx*tvol(2,jt)+dtx*dtx*tvol(3,jt)
    enddo
  enddo

  ! load particles uniformly in zeta
  !$omp parallel do private(m,tdum)
  do m=1,mp
    tdum=zpart(3,m)*deltaz
    zpart(3,m)=zeta0+tdum
  ! initialize random weight
    zpart(5,m)=w_initial*2.0*(zpart(5,m)-0.5)
#ifdef _FRC
  ! uncomment to set particular n mode 
    zpart(5,m)=zpart(5,m)*cos(zpart(3,m)*real(nmode))
#else
  ! set weight more heavily towards outer midplane
    zpart(5,m)=zpart(5,m)*(1.0+cos(zpart(2,m)))
#endif
  enddo

  if(pload==11)then ! Slowing down distribution
    fmax = sd_v0**2/(sd_v0**2+sd_vc**2)
  !$omp parallel do private(m,zdum,pdum,tdum,isp,dpx,dp2,jst,dtx,dt2,b)
    do m=1,mp
      randvarf = 1000.0_lk
      sd_ftmp  = 0.0001_lk
      do while(randvarf>sd_ftmp)
        call random_number(randvarv)
        call random_number(randvarf)
        call random_number(randvart)
        randvarf=fmax*randvarf
        randvarv=sd_v0*randvarv
        randvart=pi*randvart
 
        pdum=zpart(1,m)
        isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
        dpx=pdum-spdpsi*real(isp-1)
        dp2=dpx*dpx
      
        tdum=zpart(2,m)
        jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
        dtx=tdum-spdtheta*real(jst-1)
        dt2=dtx*dtx

        b =   bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
           (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
           (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
 
        sd_ftmp = exp(-((sin(randvart)**2/b-sd_l0)*sd_widthInv)**2)*randvarv**2/(sd_vc**3 +randvarv**3)*sin(randvart)
      enddo
        zpart(4,m)= randvarv*cos(randvart)
        zpart(6,m)= (randvarv*sin(randvart))**2
        zpart(4,m)=zpart(4,m)*apart/(qpart*b)
        zpart(6,m)=sqrt(0.5_lk*apart*zpart(6,m)/b)
        zpart0(4,m)=1.0_lk
        zpart0(6,m)=1.0_lk
    enddo
  else
  ! Maxwellian distribution in v_para, <v_para^2>=1.0, use zpartsign as temporary storage  
   tp_inv=1.0/tppp(1,1)
   np_inv=1.0/nppp(1,1)
  !$omp parallel do private(m,zdum,pdum,tdum,isp,dpx,dp2,jst,dtx,dt2,zpartsign,b,vth,rg,ii,dp1) 
   do m=1,mp
     vth=sqrt(tppp(1,1)/apart)
     if(pload==1)vth=sqrt(mesht(iflux)/apart)

     zdum=zpart(4,m)
     zpart(4,m)=zpart(4,m)-0.5
     zpartsign=sign(1.0_lk,zpart(4,m))
     zpart(4,m)=sqrt(max(1.0e-20_lk,log(1.0_lk/max(1.0e-20_lk,zpart(4,m)**2))))
     zpart(4,m)=zpart(4,m)-(c0+c1*zpart(4,m)+c2*zpart(4,m)**2)/&
       (1.0+d1*zpart(4,m)+d2*zpart(4,m)**2+d3*zpart(4,m)**3)
     if(zpart(4,m)>umax)zpart(4,m)=zdum
     zpart(4,m)=zpartsign*zpart(4,m)

    ! Maxwellian distribution in v_perp, <v_perp^2>=1.0
     zpart(6,m)=max(1.0e-20_lk,min(umax*umax,-log(max(1.0e-20_lk,zpart(6,m)))))

    ! zpart(4,:) to rho_para, zpart(6,:) to sqrt(mu)
     pdum=zpart(1,m)
     tdum=zpart(2,m)

     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

     jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
     dtx=tdum-spdtheta*real(jst-1)
     dt2=dtx*dtx
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
     majorr= xsp(1,isp,jst)+xsp(2,isp,jst)*dpx+xsp(3,isp,jst)*dp2+ &
          (xsp(4,isp,jst)+xsp(5,isp,jst)*dpx+xsp(6,isp,jst)*dp2)*dtx+ &
          (xsp(7,isp,jst)+xsp(8,isp,jst)*dpx+xsp(9,isp,jst)*dp2)*dt2
     upara0=ropp(1,isp)  +ropp(2,isp)*dpx   +ropp(3,isp)*dp2
     upara0=apart*majorr*upara0/(abs(qpart)*b) !convert upara0 to rho_para0
     
     zpart(4,m)=apart*zpart(4,m)*vth/(abs(qpart)*b)   !zpart4 is rho_para
     zpart(6,m)=vth*sqrt(apart*zpart(6,m)/b)          !zpart6 is sqrt(mu)
     zpart(7,m)=1.0

     !pload=2: non-uniform temperature & uniform density; marker weight zpart(7,:)=n_physical/n_marker
     !pload=3: non-uniform temperature & density; marker weight zpart(7,:)=1.0
     if(pload==2 .or. pload==3)then
       rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
       !radial grid on inner flux surface
       ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))
       dp1=(psimesh(ii)-pdum)*delp(ii) !weight for inner flux surface
       
       vth=sqrt((dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii))*tp_inv)
       zpart(4,m)=vth*zpart(4,m)  !use local temperature
       zpart(6,m)=vth*zpart(6,m)  !use local temperature
       if(pload==2)zpart(7,m)=(dp1*meshn(ii-1)+(1.0-dp1)*meshn(ii))*np_inv
     endif
     zpart(4,m)=zpart(4,m)+upara0 !transform into lab frame  
   enddo  
  endif


! trap=1: load trapped particles only 
  if(trapped .and. trap==1)then
     cmratio=qpart/apart
     mptmp=0
     trapfracn=real(mp,lk)
     trapfrac=0.0

! do not use omp here  
     do m=1,mp
        pdum=zpart(1,m)
        tdum=zpart(2,m)
        
        isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
        dpx=pdum-spdpsi*real(isp-1)
        if(isp==1)dpx=sqrt(dpx)
        dp2=dpx*dpx
        
        jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
        dtx=tdum-spdtheta*real(jst-1)
        dt2=dtx*dtx
        b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
             (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
             (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
        
        upara=zpart(4,m)*b*cmratio
        energy=0.5*apart*upara*upara+zpart(6,m)*zpart(6,m)*b
        
        rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
        ii=max(1,min(mpsi,int((rg-rg0)*delr+0.5)))
        bm=maxval(bmesh(igrid(ii):igrid(ii)+mtheta(ii))) !maximal B-field on flux-surface
        eperp=zpart(6,m)*zpart(6,m)*bm

! radial profile of particle and energy flux
        fullf=zpart(7,m)
        dp1=(psimesh(ii)-pdum)*delp(ii) !weight for inner flux surfac
        dden(ii-1)=dden(ii-1)+fullf*dp1
        dden(ii)  =dden(ii)+fullf*(1.0-dp1)
        
        if(eperp>energy)then
           mptmp=mptmp+1
           zpart(1,mptmp)=zpart(1,m)
           zpart(2,mptmp)=zpart(2,m)
           zpart(3,mptmp)=zpart(3,m)
           zpart(4,mptmp)=zpart(4,m)
           zpart(5,mptmp)=zpart(5,m)
           zpart(6,mptmp)=zpart(6,m)
           zpart(7,mptmp)=fullf !!zpart(7,m)
  
           trapfrac(ii-1)=trapfrac(ii-1)+fullf*dp1
           trapfrac(ii)  =trapfrac(ii)+fullf*(1.0-dp1)

        endif
     enddo
     mp=mptmp

     call MPI_ALLREDUCE(dden,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
     dden=ddum
     call MPI_ALLREDUCE(trapfrac,ddum,icount,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
     trapfrac=ddum/dden
     call MPI_ALLREDUCE(trapfracn,tdum,1,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)     
     zdum=real(mp,lk)
     call MPI_ALLREDUCE(zdum,trapfracn,1,mpi_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
     trapfracn=trapfracn/tdum
  endif
end subroutine gkLoadParticle

subroutine gk2fk
  use global_parameters
  use particle_array
  use field_array,only: nmodes,nmode
  use equilibrium
  use spline_function, only: dxdp, dxdt, dzdp, dzdt, spx
  implicit none


  integer m,isp,jst
  real(lk) pdum,tdum,dpx,dp2,dtx,dt2,gradt,gradz,gradp,b,ri,g,dxdptmp,dxdttmp,dzdptmp,dzdttmp,&
       pdot,tdot
  real(lk),dimension(:),allocatable:: alpha
  integer ierror
  
  interface
    subroutine gkLoadParticle(meshn,tppp,nppp,zpart,zpart0,wppart,&
        wtpart0,wtpart1,density,flow,data1d,jtpart0,jtpart1,wzpart,&
        marker,markert,pmark,pflux,rdtem,zonal,zonalc,mesht,qpart,apart,&
        pload,nparam,ngyro,mp,mpmax,w_initial,ncyclep,zpart1,&
        trap,trapfracn,trapfrac,hybrid,pressurepara,pressureperp,phit,dnt,&
        phisave,dnsave)

      use precision, only: lk
      use spline_function, only: construct_spline1d, invert_spline0
      implicit none
      integer pload,nparam,ngyro,mp,mpmax
      integer,optional :: ncyclep,trap,hybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(lk) qpart,apart,w_initial
      real(lk),optional :: trapfracn
      real(lk),dimension(0:) :: meshn,mesht,pmark,pflux,rdtem,zonal,zonalc
      real(lk),dimension(:) :: wzpart,marker,markert
      real(lk),dimension(:),allocatable,optional :: trapfrac
      real(lk),dimension(:,:) :: tppp,nppp
      real(lk),dimension(:,:) :: zpart,zpart0,wppart,wtpart0,wtpart1
      real(lk),dimension(0:,:) :: density,flow,data1d
      real(lk),dimension(:,:),optional :: zpart1
      real(lk),dimension(0:,:),optional :: pressurepara,&
        pressureperp,phit,dnt
      real(lk),dimension(0:,:,:),optional :: phisave,dnsave
    end subroutine gkLoadParticle
  end interface

  allocate(alpha(mimax))
  call random_number(alpha(mimax))
  alpha(:)=alpha(:)*2*pi

  !zion(1:3,m)->(psi,theta,zeta)
  !zion(4,m)->pdot
  !zion(5,m)->particle weight
  !zion(6,m)->(tdot)
  !zion(7,m)->
  !zion(8,m)->(zdot)
  !$omp parallel do private(m,pdum,isp,dpx,dp2,tdum,jst,dtx,dt2,&
  !$omp& gradt,gradz,gradp,b,ri,g,dxdptmp,dxdttmp,dzdptmp,dzdttmp)
  do m=1,mi
      pdum=zion(1,m)
      tdum=zion(2,m)
      dxdptmp=dxdp(pdum,tdum)
      dxdttmp=dxdt(pdum,tdum)
      dzdptmp=dzdp(pdum,tdum)
      dzdttmp=dzdt(pdum,tdum)
      isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
      dpx=pdum-spdpsi*real(isp-1)
      if(isp==1)dpx=sqrt(dpx)
      dp2=dpx*dpx
      jst=max(1,min(lst-1,ceiling(tdum*spdtheta_inv)))
      dtx=tdum-spdtheta*real(jst-1)
      dt2=dtx*dtx
      b=bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
        (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
        (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
      g=gpsi(1,isp)   +gpsi(2,isp)*dpx   +gpsi(3,isp)*dp2
      ri=cpsi(1,isp)   +cpsi(2,isp)*dpx   +cpsi(3,isp)*dp2
      gradp=sqrt(dxdttmp**2+dzdttmp**2)/abs(dxdptmp*dzdttmp-dxdttmp*dzdptmp)
      gradt=sqrt(dxdptmp**2+dzdptmp**2)/abs(dxdptmp*dzdttmp-dxdttmp*dzdptmp)
      gradz=1/spx(pdum,tdum)
      pdot=sqrt(zion(6,m)**2*b*2/aion)*sin(alpha(m))*gradp
      tdot=(zion(4,m)*qion/aion*ri*gradt*gradz+&
                sqrt(zion(6,m)**2*b*2/aion)*g*gradz/b*cos(alpha(m)))*gradt
      zion(8,m)=(zion(4,m)*qion/aion*g*gradz-&
                sqrt(zion(6,m)**2*b*2/aion)*ri*gradt*gradz/b*cos(alpha(m)))*gradz
      zion(4,m)=pdot
      zion(6,m)=tdot
  enddo

end subroutine gk2fk
