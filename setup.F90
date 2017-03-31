subroutine setup
 
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use particle_tracking
  use LAMYRIDGEeq
  implicit none

  integer mtheta0,mtheta1,ierror
  real(lk) delr,delt,randNumCheck
  namelist /key_parameters/ numberpe,nspecies,mgrid,mtheta0,mtheta1,mtdiag,&
    delr,delt,ulength,utime,rho0,betae,nparami,lsp,lst

! read MHD equilibrium data
  CALL EQDATA

if(numereq==14)then 
  CALL LR_initial    !initialize the LAMYRIDGE field aligned mesh and the equilibrium quantities on the field aligned mesh
endif
! allocate array for field quantites
  CALL fieldinitial

! gyro-averaging for sqrt(mu)=rho0 on grid of magnetic coordinates
  call gyroinitial

! initiate radial interpolation on poloidal grids for smooth.F90 and field.F90
  call rinterpolation

#ifdef _FRC
    call poisson_frc_initial
#else
    if(fem==0)then
    ! initiate laplacian matrix on poloidal grids for field calculate
      call lapmat_initial
    elseif(fem>0)then
      call createtriangles
      call laplacian_initial_fem
    endif
#endif

! initiate external field with non-zero boundary
  if (antenna>0) call phiext_initial

! first run plots equilibrium profile
  if(mype==0)then
     if(numereq<4) CALL EQPLOT
     if(numereq==11)CALL Eqplot_Cylinder_RZ
  endif     
! electromagnetic first run w/o fast ion, output Alfven continua data
! w/ fast ion, the mass density and the pressure contributing to the continua are not well-defined
! so in order to produce Alfven continua data, set fload==0
!  if(mype==0 .and. irun==0 .and. magnetic/=0 .and. fload==0)call alcondata
if(numereq==14)then  ! for FRC LR equilibrium in cylindrical coordinates
! number of particle  
  mi=int(real(mi)*real((numr-1)*(numz-1))/real(npartdom))   !# of ions per MPI process
  me=int(real(me)*real((numr-1)*(numz-1))/real(npartdom))   !# of electrons per MPI process
! number of particle  
  mf=int(real(mf)*real((numr-1)*(numz-1))/real(npartdom))   !# of fast ions per MPI process
  mfe=int(real(mfe)*real((numr-1)*(numz-1))/real(npartdom)) !# of fast electrons per MPI process
else    ! for tokamak equilibrium in Boozer coordinates
! number of particle
  mi=int(real(mi)*real(mgrid-mpsi)/real(npartdom))          !# of ions per MPI process
  me=int(real(me)*real(mgrid-mpsi)/real(npartdom))          !# of electrons per MPI process
  mf=int(real(mf)*real(mgrid-mpsi)/real(npartdom))          !# of fast ions per MPI process
  mfe=int(real(mfe)*real(mgrid-mpsi)/real(npartdom))        !# of fast electrons per MPI process 
endif
! # of species
  nspecies=1
  if(nhybrid>0)nspecies=nspecies+1
  if(fload>0 .and. feload >0)then
    write(0,*)'It is not availale when  feload and fload both >0'
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
  if(fload>0)nspecies=nspecies+1
  if(feload>0)nspecies=nspecies+1 

! Particles is tracked by tagging them with a number with an extra element to the particle array
  nparami=7
  nparame=7
  nparamf=7
  nparamfe=7
  if(track_particles == 1)then
    nparami=nparami+2
    nparame=nparame+2
  endif
  if(iload==9) nparami=nparami+1

#ifndef GPU_UM
  !$acc update device(nparami,nparame,nparamf,nparamfe)
#endif
!write test random number to check RNG initialization
  call random_number(randNumCheck)
  write(gtcout,*)'Rand Num Check = ',randNumCheck
! write out key parameter
  if(mype == 0) then
    delr=deltar*qion/sqrt(aion*meshti(iflux))
    delt=deltat(iflux)*(rg0+deltar*real(iflux))*qion/sqrt(aion*meshti(iflux))
    mtheta0=mtheta(0)
    mtheta1=mtheta(mpsi)
    write(gtcout,key_parameters)
    call FLUSH(gtcout)
  endif

  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

end subroutine setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fieldinitial
  use global_parameters
  use equilibrium
  use field_array
  use particle_array
  use magnetic_island
  use spline_function, only: spq,sppsirg,spb,sptorpsi,dtorpsi,sprgpsi,spgpsi,spcpsi,spzeff
  implicit none

  integer mtest,ierror,i,j,ii,ij,n,m,win0,win1,isp
  real(lk) envelope,rg,q,tdum,psi,sperr,delr,pleft,pright,pdum,dpx,dp2
  real(8) temp,dbesjn

! allocate memory
  allocate(qtinv(0:mpsi),qmesh(0:mpsi),meshze(0:mpsi),psimesh(0:mpsi),&
    tormesh(0:mpsi),itran(0:mpsi),mtheta(0:mpsi),deltap(mpsi),deltat(0:mpsi),&
    igrid(0:mpsi),phi00(0:mpsi),phip00(0:mpsi),apara00(0:mpsi),&
    apara00nl(0:mpsi),apara00nl0(0:mpsi),fluidne00(0:mpsi),igrid_fem(0:mpsi),STAT=mtest)
  if(mtest /= 0)then
    write(0,*)mype,'*** Cannot allocate qtinv: mtest=',mtest
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
  phi00=0.0
  phip00=0.0
  apara00=0.0
  apara00nl=0.0
  apara00nl0=0.0
  fluidne00=0.0

!analyticmesh provides temporary solution for inaccuracy of sqrt(x) spline functions near origin 
  !uncomment if necessary
  !call analyticmesh
! Define simulation domain and radial grid in psi, assume sprgpsi is an increasing function
  psi0=psi0*ped
  psi1=psi1*ped
  psimesh(0)=psi0
  psimesh(mpsi)=psi1
  tormesh(0)=sptorpsi(psi0)
  tormesh(mpsi)=sptorpsi(psi1)
  qmesh(0)=spq(psi0)
  qmesh(mpsi)=spq(psi1)
  rg0=sprgpsi(psi0)
  deltar=(sprgpsi(psi1)-rg0)/real(mpsi)
  sperr=0.0000001

  do i=1,mpsi-1
     rg=rg0+deltar*real(i) !target rg
     pleft=psimesh(i-1)
     pright=psimesh(mpsi)
     psi=0.5*(pleft+pright) !initial guess of psi
     delr=1.0-sprgpsi(psi)/rg
     ii=1

     do while (abs(delr)>sperr .and. ii<1000)
        ii=ii+1
        if(delr>0)then
           pleft=psi
        else
           pright=psi
        endif
        psi=0.5*(pleft+pright) !new value of psi
        delr=1.0-sprgpsi(psi)/rg
     enddo

!poloidal, toroidal, q flux on radial mesh
     psimesh(i)=psi
     tormesh(i)=sptorpsi(psi)
     qmesh(i)=spq(psi)

     if(mype==0)write(gtcout,*)i,sprgpsi(psimesh(i))/(rg0+deltar*real(i))-1.0
  enddo
  !endif

! radial grid size in psi
  do i=1,mpsi
     deltap(i)=psimesh(i)-psimesh(i-1)
  enddo

  if(mype==0)then
     write(gtcout,*)"rg0,rg1=",rg0/sprgpsi(ped),(rg0+deltar*mpsi)/sprgpsi(ped),&
        "psi0,spdpsi=",psimesh(0),spdpsi
        if(psimesh(0)<spdpsi)write(gtcout,*)"WARNING: psimesh(0)<spdpsi"

     write(gtcout,*)"spline consistency check: 4th=0; 5th=1"
     write(gtcout,*)"at last closed flux surface a_minor=",sprgpsi(ped),"at wall a_minor=",sprgpsi(psiw)

     do i=1,mpsi-1
        rg=rg0+deltar*real(i)
#ifdef _FRC
          write(gtcout,*)i,rg,qmesh(i),sprgpsi(psimesh(i))/rg-1.0
#else
          write(gtcout,*)i,rg/sprgpsi(ped),psimesh(i)/ped,qmesh(i),sprgpsi(psimesh(i))/rg-1.0,dtorpsi(psimesh(i))/qmesh(i)
#endif
     enddo
  endif

! particle temperature and density on radial mesh
  allocate(meshni(0:mpsi),kapani(0:mpsi),meshti(0:mpsi),kapati(0:mpsi),&
       meshte(0:mpsi),meshne(0:mpsi),kapane(0:mpsi),kapate(0:mpsi),&
       meshnf(0:mpsi),meshtf(0:mpsi),kapanf(0:mpsi),kapatf(0:mpsi),&
       meshnfe(0:mpsi),meshtfe(0:mpsi),kapanfe(0:mpsi),kapatfe(0:mpsi),&
       dtndpsi(0:mpsi),jacobianpsi(0:mpsi),mesher(0:mpsi))

! equilibrium quantities on radial grid points
!$omp parallel do private(i,pdum,isp,dpx,dp2)
  do i=0,mpsi
     pdum=psimesh(i)
     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     dp2=dpx*dpx

! 1D spline in psi
     meshze(i)=zepp(1,isp)+zepp(2,isp)*dpx+zepp(3,isp)*dp2

     meshni(i)=nipp(1,isp)+nipp(2,isp)*dpx+nipp(3,isp)*dp2
     kapani(i)=       0.0-(nipp(2,isp)+2.0*nipp(3,isp)*dpx)/meshni(i)
     meshti(i)=tipp(1,isp)+tipp(2,isp)*dpx+tipp(3,isp)*dp2
     kapati(i)=       0.0-(tipp(2,isp)+2.0*tipp(3,isp)*dpx)/meshti(i)
 
     meshnf(i)=nfpp(1,isp)+nfpp(2,isp)*dpx+nfpp(3,isp)*dp2
     kapanf(i)=       0.0-(nfpp(2,isp)+2.0*nfpp(3,isp)*dpx)/max(1.0e-10,meshnf(i))
     meshtf(i)=tfpp(1,isp)+tfpp(2,isp)*dpx+tfpp(3,isp)*dp2
     kapatf(i)=       0.0-(tfpp(2,isp)+2.0*tfpp(3,isp)*dpx)/meshtf(i)
     
     meshnfe(i)=nfepp(1,isp)+nfepp(2,isp)*dpx+nfepp(3,isp)*dp2
     kapanfe(i)=      0.0-(nfepp(2,isp)+2.0*nfepp(3,isp)*dpx)/max(1.0e-10,meshnfe(i))
     meshtfe(i)=tfepp(1,isp)+tfepp(2,isp)*dpx+tfepp(3,isp)*dp2
     kapatfe(i)=      0.0-(tfepp(2,isp)+2.0*tfepp(3,isp)*dpx)/meshtfe(i) 

#ifdef _FRC
     meshne(i)=meshnfe(i)
     kapane(i)=kapanfe(i)
     meshte(i)=meshtfe(i)
     kapate(i)=kapatfe(i)
#else
     meshne(i)=nepp(1,isp)+nepp(2,isp)*dpx+nepp(3,isp)*dp2
     kapane(i)=       0.0-(nepp(2,isp)+2.0*nepp(3,isp)*dpx)/meshne(i)
     meshte(i)=tepp(1,isp)+tepp(2,isp)*dpx+tepp(3,isp)*dp2
     kapate(i)=       0.0-(tepp(2,isp)+2.0*tepp(3,isp)*dpx)/meshte(i)
#endif
     mesher(i)=erpp(1,isp)+erpp(2,isp)*dpx+erpp(3,isp)*dp2
  enddo

!kill drive at simulation boundaries
  if(nbound>99)then  
    call gaussbound(meshze)
    call gaussbound(kapani)
    call gaussbound(kapati)
    call gaussbound(kapane)
    call gaussbound(kapate)
    call gaussbound(kapanfe)
    call gaussbound(kapatfe)
    call gaussbound(kapanf)
    call gaussbound(kapatf)
    call gaussbound(mesher)
    if(nboundR>99)then  
      call gaussboundR(meshze)
      call gaussboundR(kapani)
      call gaussboundR(kapati)
      call gaussboundR(kapane)
      call gaussboundR(kapate)
      call gaussboundR(kapanfe)
      call gaussboundR(kapatfe)
      call gaussboundR(kapanf)
      call gaussboundR(kapatf)
      call gaussboundR(mesher)
    endif
    if(mype==0)then
      if(irotation==0)then
        write(gtcout,*)'kappati,kappate'
        do i=0,mpsi
          write(gtcout,*)i,kapati(i),kapate(i)
        enddo
      else
        write(gtcout,*)'kappati,kappate,mesher'
        do i=0,mpsi
          write(gtcout,*)i,kapati(i),kapate(i),mesher(i)
        enddo
      endif
    endif
  endif
!by zhs set initialize iflux
!For Alven Eigenmode simulations, normalize via iflux values. 
  if(mype==0)write(gtcout,*)"iflux at i=",iflux  
  if(iload==0)then ! ideal MHD, only electron pressure (effectively T_i=0)
    nhybrid=0
    fload=0
    feload=0 
    if(eload==1)then !uniform marker profile if eload==1
      meshni=nipp(1,1)
      meshti=tipp(1,1)
    endif
  elseif(iload==1)then
#ifndef _FRC
    if(inorm==0)then  
      kapani=kapani*meshni/meshni(iflux)
      kapati=kapati*meshti/meshti(iflux)
    endif  
    meshni=meshni(iflux)
    meshti=meshti(iflux)
    if(fload>0)fload=1 ! make sure that when iload=1, fload=1
    if(nhybrid>0)eload=1 ! make sure that when iload=1, eload=1
    if(feload>0)feload=1 ! make sure that when iload=1, feload=1 
  endif
  if(eload==1)then
    if(inorm==0)then  
      kapane=kapane*meshne/meshne(iflux)
      kapate=kapate*meshte/meshte(iflux)
    endif
    meshne=meshne(iflux)
    meshte=meshte(iflux)
  endif
  if(fload==1)then
    if(inorm==0)then  
      kapanf=kapanf*meshnf/meshnf(iflux)
      kapatf=kapatf*meshtf/meshtf(iflux)
    endif
    meshnf=meshnf(iflux)
    meshtf=meshtf(iflux)
  endif
  if(feload==1)then
    if(inorm==0)then
      kapanfe=kapanfe*meshnfe/meshnfe(iflux)
      kapatfe=kapatfe*meshtfe/meshtfe(iflux)
    endif
    meshnfe=meshnfe(iflux)
    meshtfe=meshtfe(iflux)
  endif
#else
  meshni=nipp(1,1)
  meshti=tipp(1,1)
  meshne=nepp(1,1)
  meshte=tepp(1,1)
  meshnf=nfpp(1,1)
  meshtf=tfpp(1,1)
  meshnfe=nepp(1,1)
  meshtfe=tepp(1,1)  
#endif
  if(iload>1 .and. (fload==1 .or. eload==1 .or. feload==1))write(gtcout,*)"warning: iload, fload, eload=", iload,fload,eload
#ifdef _FRC
     ! numberpe must be a multiple of mthetamax
     if(mod(numberpe,mthetamax) /= 0)then
         write(gtcout,*)'Wrong PE # for FRC simulation; numberpe=',numberpe,'mthetamax=',mthetamax
         stop
     endif
     qtinv=0.0
     itran=0
     mtheta(0:mpsi)=mthetamax
     deltat(0:mpsi)=2.0_lk*pi/real(mthetamax)
#else
! Define poloidal grid: grid shift associated with fieldline following coordinates
     tdum=2.0*(rg0+deltar*real(mpsi/2))/(sqrt(meshti(mpsi/2))*real(mthetamax))
     do i=0,mpsi
        rg=rg0+deltar*real(i)
        mtheta(i)=2*max(1,int(rg/(tdum*sqrt(meshti(i)))+0.5_lk))
        deltat(i)=2.0_lk*pi/real(mtheta(i))
        q=qmesh(i)
        itran(i)=nint(real(mtheta(i))/(q*toroidaln))
        if (abs(itran(i))>0) then 
          qtinv(i)=real(mtheta(i))/(real(itran(i))*toroidaln) !q value for coordinate transformation
          qtinv(i)=1.0/qtinv(i) !inverse q to avoid divide operation
        else
        ! to avoid NaN in the case when itran(i)=0
        ! which is essentially the linear machine with only axial field case.
          qtinv(i)=0.0
        endif
        itran(i)=itran(i)-mtheta(i)*(itran(i)/mtheta(i))
        if(psi0==0.0_lk)mtheta(0)=1
     enddo
#endif
! un-comment the next two lines to use magnetic coordinate
!  qtinv=0.0
!  itran=0
 write(gtcout,*)"toroidaln =",toroidaln
! total number of grids on a poloidal plane
  mgrid=sum(mtheta+1)

! When doing filtering and diagnostics of toroidal mode, we need to switch from the 
! field-line following coordinates alpha-zeta to the magnetic coordinate in theta-zeta. 
! This requires a greater number of grid points in the zeta direction, which
! is mtdiag. Precisely, mtdiag should be mtheta/q but since mtheta changes
! from one flux surface to another, we use mtdiag=mtheta(mpsi)
  mthetamax=mtheta(mpsi)
  mtdiag=((mthetamax/numberpe)+1)*numberpe

! starting # in the label of poloidal grid on each flux surface
  igrid(0)=1
  do i=1,mpsi
     igrid(i)=igrid(i-1)+mtheta(i-1)+1
  enddo
  if(fem>0)then 
     igrid_fem(0)=0
     do i=1,mpsi
        igrid_fem(i)=igrid_fem(i-1)+mtheta(i-1)!does not include poloidal BC
     enddo
  endif 
 
!bcond 0: solution is zero at boundaries, 1: solution is assumed linear at inner boundary
!mpsilow is transition region between FD and exterpolation to origin with linear boundary condition
  if(rg0/sprgpsi(ped)<0.15 .and. bcond==1) then 
     mpsilow=ceiling(real(mpsi)*(0.15-rg0/sprgpsi(psiw))/3.0)
  else
     mpsilow=1
  endif 
  mpsihigh=mpsi-1
  mgridlow=igrid(mpsilow)
  mgridhigh=igrid(mpsihigh+1)-1
  if(fem>0)mgrid_fem=mgrid-(mpsi+1)!total grid points minus poloidal BC

!! xy toroidal averaging
  allocate(wtshift(0:mpsi),nshift(0:mpsi),STAT=mtest)
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate nshift: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  do i=0,mpsi
     if(fielddir==1 .or. fielddir==3)then
       tdum=(zeta1-pi2)*qtinv(i)
     else
       tdum=zeta1*qtinv(i)
     endif
     wtshift(i)=modulo(tdum,deltat(i))
     nshift(i)=modulo(floor(tdum/deltat(i)), mtheta(i))
  enddo
  wtshift=wtshift/deltat


! flux surface for mode diagnosis and filtering; default: 8 modes and on mpsi/2 flux surface
  modes=num_modes
  allocate(nmodes(modes),mmodes(modes))
! nfilter=0: all toroidal modes dynamics kept; diagonosis modes # k*delta_x=0.1,...,0.4
!  do i=1,modes
!     mmodes(i)=int(real(i*(mtheta(iflux)/2-1))/(10.0*pi))
!     nmodes(i)=int(0.5+real(mmodes(i))*abs(qtinv(iflux)))
!  enddo

!  if(nfilter>0)then
! or nfilter=1: n-modes explicitly specified
     nmodes(1:modes)=n_modes(1:modes)
     mmodes(1:modes)=m_modes(1:modes)
!  endif
  !if(mype==0)write(gtcout,"(a8,8i4)")"nmodes=",nmodes,"mmodes=",mmodes
  if(mype==0)write(gtcout,*)"nmodes=",nmodes
  if(mype==0)write(gtcout,*)"mmodes=",mmodes

! allocate memory
  allocate(phi(0:1,mgrid),phi_zero(0:1,mgrid),apara(0:1,mgrid),apara0(0:1,mgrid),&
       fluidne(0:1,mgrid),fluidne0(0:1,mgrid),fluidue(0:1,mgrid),gradue(3,0:1,mgrid),jmesh(mgrid),&
       gradphi(3,0:1,mgrid),gradapara(3,0:1,mgrid),bmesh(mgrid),gradne(3,0:1,mgrid),sapara(0:1,mgrid),sdelapara(0:1,mgrid),&
       deltapsi(0:1,mgrid),deltapsi0(0:1,mgrid),sdeltapsi(0:1,mgrid),gradpsi(3,0:1,mgrid),&
       gradphieff(3,0:1,mgrid),sfluidne(0:1,mgrid),gradext(3,0:1,mgrid),gradkne(3,0:1,mgrid),&
       gradpepara(3,0:1,mgrid),gradpeperp(3,0:1,mgrid),d4fluidne(mgrid),d2apara(mgrid),&
       gradgaugef(3,0:1,mgrid),MHDxiperp(3,0:1,mgrid),MHDxi_mag(0:1,mgrid),STAT=mtest)



! B-field & inverse Jacobian on poloidal mesh points
!$omp parallel do private(i,j,ij,psi,tdum)
  do i=0,mpsi
     psi=psimesh(i)
     jacobianpsi(i)=0.0
     do j=0,mtheta(i)
        ij=igrid(i)+j
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
        else
          tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
        endif
        bmesh(ij)=spb(psi,tdum)
        jmesh(ij)=bmesh(ij)**2/(qmesh(i)*spgpsi(psi)+spcpsi(psi))
        jacobianpsi(i)=jacobianpsi(i)+1.0/jmesh(ij)
     enddo
!delete the double counting first/last point on the fluxsurface
     jacobianpsi(i)=jacobianpsi(i)-1.0/jmesh(igrid(i))
     jacobianpsi(i)=jacobianpsi(i)/real(mtheta(i))
  enddo
! for calculating poisson zonal
  dtndpsi=(kapani-kapati)*meshti/meshni!!*meshni(0)/meshti(0)

! set up geometric tensor for zonal flow solver
  CALL SETGEOM

! inital adiabative electron density fluidne=cos(m*mtheta-n*zeta) with a radal envelope
! initial vector potential apara=0. All other dynamical quantities=0
  phi=0.0
  apara=0.0
  fluidne=0.0
  fluidue=0.0
  gradue=0.0
  gradphi=0.0
  gradapara=0.0
  gradne=0.0
  gradkne=0.0
  gradpepara=0.0
  gradpeperp=0.0
  sfluidne=0.0
  sapara=0.0
  sdelapara=0.0
  gradext=0.0
  deltapsi=0.0
  deltapsi0=0.0
  sdeltapsi=0.0
  gradpsi=0.0
  gradphieff=0.0
  d4fluidne=0.0
  d2apara=0.0

  gradgaugef=0.0
  MHDxiperp=0.0
  MHDxi_mag=0.0

! initial perturbation of electron current
  if(magnetic==1)then
     n=n_modes(1)
     m=m_modes(1)
     if(fielddir==1 .or. fielddir==3)n=-n
!by zhs use old initial profile
     win0=int(mpsi/4)
     win1=int(mpsi*3/4)
!    do i=win0,win1
!        envelope=1.0e-10 ! set this quantity to 0.0 for no initial field perturbation
!        envelope=envelope*sin(pi*real(i-win0)/real(win1-win0))

!~~~~~~~~~~~new smooth profile~~~~~~~~~~~~~~~~~~~~~
     do i=0,mpsi
!         rg=sqrt(psimesh(i)/psimesh(mpsi))
!         rg=20.0*(rg-0.5)
!         rg=rg*rg*rg*rg*rg*rg
!         rg=20.0*(rg-0.5)
!         temp=9.93611*rg
!         envelope=1.0e-6*dbesjn(m,temp)
          rg=real(i)/real(mpsi)
          rg=20.0*(rg-0.5)
          rg=rg*rg
          envelope=1.0e-5*exp(-rg)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         q=qmesh(i)

!$omp parallel do private(j,ij,tdum)
        do j=1,mtheta(i)
           ij=igrid(i)+j
           if(fielddir==1 .or. fielddir==3)then
             tdum=deltat(i)*real(j)+(zeta1-pi2)*qtinv(i)
           else
             tdum=deltat(i)*real(j)+zeta1*qtinv(i)
           endif
           !!change if necessary!!
           fluidne(1,ij)=1.0*envelope*(cos(real(m)*tdum-real(n)*zeta1)+&
                                       0.5*cos(real(m+1)*tdum-real(n)*zeta1)+&
                                       0.5*cos(real(m-1)*tdum-real(n)*zeta1))
           apara(1,ij)  = 0.0*envelope*cos(real(m)*tdum-real(n)*zeta1)
        enddo
     enddo
     sfluidne=fluidne
  endif

! define island structure
  if(island==1)then
    call islandinitial
  endif
! determines symmetric b-filed positions for particle bc
  if(numereq/=0) call bminmax

  ! copy grid parameters to GPU
#ifndef GPU_UM
  !$acc update device(mtheta,qtinv,psimesh)
  !$acc update device(igrid,deltap,deltat)
  !$acc update device(mesher)
  !$acc update device(meshti,meshni,kapani,kapati)
  !$acc update device(meshte,meshne,kapane,kapate)
  !$acc update device(meshtf,meshnf,kapanf,kapatf)
  !$acc update device(meshtfe,meshnfe,kapanfe,kapatfe)
#endif
end subroutine fieldinitial

! all the equations used in subroutine islandinitial can be found in the note on phoenix
subroutine islandinitial
  use global_parameters
  use equilibrium
  use field_array
  use particle_array
  use magnetic_island
  use spline_function, only: sprgpsi
  implicit none
  
  integer mtest,i,j,ij,n,m
  !integer ires,wi,wo,qs ! i of resonant surface, inner island width, outter island width, q of resonant surface
  real(lk) tdum
  real(lk),dimension(:),allocatable :: deltaq,hmesh,hmesh0,dhmesh0 ! delta q between different surfaces, psimesh-tormesh/qs Eq(14), hmesh normalized by the value of resonant surface Eq(15), delta hmesh0 between different surfaces
  real(lk),dimension(:,:),allocatable :: iarea,alpnm ! indicate whether a grid point is inside the island or not. if it is inside the island, the value is 1 otherwise 0
  
! allocate memeory
  allocate(hmesh(0:mpsi),hmesh0(0:mpsi),dhmesh0(0:mpsi),deltaq(0:mpsi),STAT=mtest)
  allocate(alphaIs(0:1,mgrid),hmesh_total(0:1,mgrid),hmesh_perturb(0:1,mgrid),hangle(0:1,mgrid),gradalphaIs(3,0:1,mgrid),iarea(0:1,mgrid),STAT=mtest)

  hmesh_perturb=0

  allocate(isn(l),ism(l))
  allocate(alpnm(l,0:mpsi)) 
  isn=(1,1)
  ism=(2,3)
  alphaIs=0.0
! q of the resonant surface
!  qs=real(polnum)/real(tornum)

! find the resonant surface
!  do i=0,mpsi
!     deltaq(i)=abs(qmesh(i)-real(polnum)/real(tornum))
!  enddo

 ! do i=1,mpsi-1
 !    if(deltaq(i-1)>deltaq(i).and.deltaq(i)<deltaq(i+1))then
 !      ires=i
 !    endif
  !enddo

  !do i=0,mpsi
  !   hmesh(i)=psimesh(i)-tormesh(i)/qs
  !enddo

! define equilibrium helical flux
! Use the condition: hmesh0(w)-hmesh0(ires)=2*alpha to find the island boundaries, Eq(17)
 ! do i=0,mpsi
 !    hmesh0(i)=hmesh(i)-hmesh(ires) !normalize the helical flux by that value of the resonant surface
 !    dhmesh0(i)=abs(hmesh0(i)-2.0*alp)
 ! enddo

! find the island boundaries in psi
  !do i=1,ires-1
  !   if(dhmesh0(i-1)>dhmesh0(i).and.dhmesh0(i)<dhmesh0(i+1))then
  !     wi=i
  !   endif
  !enddo
  
 ! do i=ires+1,mpsi-1
  !   if(dhmesh0(i-1)>dhmesh0(i).and.dhmesh0(i)<dhmesh0(i+1))then
  !     wo=i
  !   endif
  !enddo
!do i=1,mpsi
!if(abs(qmesh(i)-2)<abs(qmesh(i-1)-2))ires=i
!enddo
!do i=1,mpsi
!hmesh(i)=psimesh(i)-tormesh(i)*0.5-psimesh(ires)+tormesh(ires)*0.5-alp*2.0
!enddo
!do i=2,ires
!if(abs(hmesh(i))<abs(hmesh(i-1)))wi=i
!enddo
!do i=ires+2,mpsi
!if(abs(hmesh(i))<abs(hmesh(i-1)))wo=i
!enddo

!if(mype==0)then
!write(gtcout,*)"~^_^~,in,ires,out",wi,ires,wo
!write(gtcout,*)sprgpsi(psimesh(wi)),sprgpsi(psimesh(ires)),sprgpsi(psimesh(wo))
!endif
! un-comment the next two lines to use magnetic coordinate

! write out the ratio of island widths over minor radius   
!  if(mype==0)then
!  write(gtcout,*) "wi/a=",(sprgpsi(psimesh(ires))-sprgpsi(psimesh(wi)))/sprgpsi(ped),"wo/a=",(sprgpsi(psimesh(wo))-sprgpsi(psimesh(ires)))/sprgpsi(ped)
!  endif

! define magnetic perturbation inside the island
 ! do i=wi,wo
 !    do j=0,mtheta(i)
 !       ij=igrid(i)+j
  !      hmesh_perturb(1,ij)=alp
  !   enddo
  !enddo
 
  !do i=0,wi-1
  !   do j=0,mtheta(i)
  !      ij=igrid(i)+j
  !      hmesh_total(1,ij)=hmesh0(wi)+(hmesh0(i)-hmesh0(wi))/(hmesh0(0)-hmesh0(wi))*alp
        !hmesh_total(1,ij)=0
        !hfunction(1,ij)=(psimesh(i)-psimesh(ires))*(1-1.2**(-(abs(hmesh0(i)-abs(2*alp)))/abs(2*alp)))
  !   enddo
  !enddo

!define the island modification on lambda alpha_island & helical flux
!use the constant-psi condition on separatrix, which means the value of helical flux on separatrix is alpha
!if the helical flux of ij grid is smaller than alpha, this grid is outside the island, which means the island modification on lambda of this grid is zero

!!$omp parallel do private(i,j,ij,tdum)  
  !do i=wi,wo
 !    do j=0,mtheta(i)
 !       ij=igrid(i)+j
 !       if(fielddir==1 .or. fielddir==3)then
 !         tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
  !      else
   !       tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
   !     endif
    !    hangle(1,ij)=real(polnum)*tdum-real(tornum)*zeta1
    !    hmesh_total(1,ij)=hmesh0(i)-hmesh_perturb(1,ij)*cos(hangle(1,ij))
     !   if (hmesh_total(1,ij)<alp)then ! judgment condition  Eq(16)
     !      iarea(1,ij)=0.0
     !   else
      !     iarea(1,ij)=1.0
      !  endif
       ! hmesh_total(1,ij)=(hmesh0(i)-hmesh_perturb(1,ij)*cos(hangle(1,ij)))*iarea(1,ij)
        !ialpha(1,ij)=hmesh_perturb(1,ij)*iarea(1,ij)*cos(hangle(1,ij)) ! modification of lambda in equations of motion 
    ! enddo
  !enddo

 ! do i=wo+1,mpsi
 !    do j=0,mtheta(i)
 !       ij=igrid(i)+j
 !       hmesh_total(1,ij)=hmesh0(wo)+(hmesh0(i)-hmesh0(wo))/(hmesh0(mpsi)-hmesh0(wo))*alp
        !hmesh_total(1,ij)=0
        !hfunction(1,ij)=(psimesh(i)-psimesh(ires))*(1-1.2**(-(abs(hmesh0(i)-abs(2*alp)))/abs(2*alp)))
  !   enddo
!  enddo

   do n=1,l
     do i=0,mpsi
     alpnm(n,i)=-0.000008 !example profile of the island
     if(n==2)alpnm(n,i)=-0.00030 !example profile of the island
!$omp parallel do private(j,ij,tdum)
        do j=1,mtheta(i)
           ij=igrid(i)+j
             tdum=deltat(i)*real(j)+zeta1*qtinv(i)
           alphaIs(1,ij)=alphaIs(1,ij)+alpnm(n,i)*cos(real(ism(n))*tdum-real(isn(n))*zeta1)
         enddo
      enddo
    enddo
 ! call periodicity(hangle)
 ! call periodicity(hmesh_perturb)
 ! call periodicity(hmesh_total) 
  call periodicity(alphaIs)

end subroutine islandinitial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!XY: setup geometric tensors
subroutine setgeom
  use global_parameters
  use field_array
  use equilibrium
  use particle_array,only:aion,meshti,qion
  use spline_function, only: dxdp,dxdt,dxdz,dzdp,dzdt,dzdz,dfdp,dfdt,dfdz,spb,spx,sprpsi
  implicit none

  integer i,j,ij,isp,mtest,ierror
  real(lk) pdum,tdum,bdum,dpx,q,b,functioni,functiong,delpsi
  real(lk) b2m(mgrid),adum(0:mpsi),markeritmp(mgrid),pmarkitmp(0:mpsi)
  real(lk) xtmp,delmat,rdum
  real(lk) dxdptmp,dxdttmp,dxdztmp,dzdptmp,dzdttmp,dzdztmp,dfdptmp,dfdttmp,dfdztmp

  markeritmp=0.0
! Jacobian=(gq+I)/B^2
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
        markeritmp(ij)=delpsi*deltat(i)*(functiong*q+functioni)/(b*b)
     enddo
        markeritmp(igrid(i))=markeritmp(igrid(i)+mtheta(i))
  enddo

!jacobian flux surface average
! cannot use OpenMP in this loop
  pmarkitmp=0.0
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        pmarkitmp(i)=pmarkitmp(i)+markeritmp(ij)
     enddo
  enddo

!geometric tensor for zonal flow solver
! allocate memory
  allocate(gpsi200(0:mpsi),gupp(mgrid),gutt(mgrid),guzz(mgrid),&
     gupt(mgrid),gupz(mgrid),gutz(mgrid),&
     b2m00(0:mpsi),gdpp(mgrid),gdtt(mgrid),gdzz(mgrid),&
     gdpt(mgrid),gdpz(mgrid),gdtz(mgrid),STAT=mtest)
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate qtinv: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  gpsi200=0.0
  b2m00=0.0
  gupp=0.0
  gutt=0.0
  guzz=0.0
  gupt=0.0
  gupz=0.0
  gutz=0.0

  gdpp=0.0
  gdtt=0.0
  gdzz=0.0
  gdpt=0.0
  gdpz=0.0
  gdtz=0.0

  adum=0.0

  dxdptmp=0.0
  dxdttmp=0.0
  dxdztmp=0.0
  dzdptmp=0.0
  dzdttmp=0.0
  dzdztmp=0.0
  dfdptmp=0.0
  dfdttmp=0.0
  dfdztmp=1.0

  do i=0,mpsi
     pdum=psimesh(i)

     do j=1,mtheta(i)
        ij=igrid(i)+j
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
        else
          tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
        endif
        dxdptmp=dxdp(pdum,tdum)
        dxdttmp=dxdt(pdum,tdum)
        dzdptmp=dzdp(pdum,tdum)
        dzdttmp=dzdt(pdum,tdum)
        xtmp=spx(pdum,tdum)
       !!! 3D-case
#ifdef _TOROIDAL3D
          dxdztmp=dxdz(pdum,tdum)
          dzdztmp=dzdz(pdum,tdum)
          dfdptmp=dfdp(pdum,tdum)
          dfdttmp=dfdt(pdum,tdum)
          dfdztmp=dfdz(pdum,tdum)
#endif
        if(fielddir==2.or.fielddir==3)dfdztmp=-dfdztmp !direction of zeta is negative w.r.t. phi (R,phi,Z)-standard cylindrical coordinate 
        gdpp(ij)=dxdptmp*dxdptmp+dzdptmp*dzdptmp+dfdptmp*dfdptmp*xtmp*xtmp
        gdtt(ij)=dxdttmp*dxdttmp+dzdttmp*dzdttmp+dfdttmp*dfdttmp*xtmp*xtmp
        gdzz(ij)=dxdztmp*dxdztmp+dzdztmp*dzdztmp+dfdztmp*dfdztmp*xtmp*xtmp
        gdpt(ij)=dxdptmp*dxdttmp+dzdptmp*dzdttmp+dfdptmp*dfdttmp*xtmp*xtmp
        gdpz(ij)=dxdptmp*dxdztmp+dzdptmp*dzdztmp+dfdptmp*dfdztmp*xtmp*xtmp
        gdtz(ij)=dxdttmp*dxdztmp+dzdttmp*dzdztmp+dfdttmp*dfdztmp*xtmp*xtmp
        delmat=gdpp(ij)*gdtt(ij)*gdzz(ij)-gdpt(ij)*gdpt(ij)*gdzz(ij)-&
               gdtz(ij)*gdtz(ij)*gdpp(ij)-gdpz(ij)*gdpz(ij)*gdtt(ij)+&
               2.0*gdpt(ij)*gdpz(ij)*gdtz(ij)


        gupp(ij)=(gdtt(ij)*gdzz(ij)-gdtz(ij)*gdtz(ij))/delmat
        gupt(ij)=(gdpz(ij)*gdtz(ij)-gdpt(ij)*gdzz(ij))/delmat
        gutt(ij)=(gdpp(ij)*gdzz(ij)-gdpz(ij)*gdpz(ij))/delmat
        guzz(ij)=(gdpp(ij)*gdtt(ij)-gdpt(ij)*gdpt(ij))/delmat

#ifdef _CYLINDER
          guzz(ij)=1.0_lk
          gdzz(ij)=1.0_lk
#endif

        bdum=spb(pdum,tdum)
        b2m(ij)=bdum*bdum
     enddo
  enddo
!Setting g-operators to be zero for when psi0=0, avoids NaNs
  if(psi0==0.0)then
    i=0
    do j=1,mtheta(i)
       ij = igrid(i)+j
       gupp(ij)=0.0
       gutt(ij)=0.0
       guzz(ij)=0.0
       gupt(ij)=0.0
       gupz(ij)=0.0
       gutz(ij)=0.0
       gdpp(ij)=0.0
       gdtt(ij)=0.0
       gdzz(ij)=0.0
       gdpt(ij)=0.0
       gdpz(ij)=0.0
       gdtz(ij)=0.0
    enddo
  endif
!! periodicity and poloidal smoothness
  call polsmooth(gupp,1)
  call polsmooth(gutt,1)
  call polsmooth(guzz,2)
  call polsmooth(gupt,1)
  call polsmooth(gupz,1)
  call polsmooth(gutz,1)

  call polsmooth(gdpp,1)
  call polsmooth(gdtt,1)
  call polsmooth(gdzz,2)
  call polsmooth(gdpt,1)
  call polsmooth(gdpz,1)
  call polsmooth(gdtz,1)


!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j
        gpsi200(i)=gpsi200(i)+gupp(ij)*markeritmp(ij)
        b2m00(i)=b2m00(i)+b2m(ij)*markeritmp(ij)
     enddo
  enddo

! global sum of zonal modes (flux surface averaged), broadcast to every toroidal PE
  call MPI_ALLREDUCE(gpsi200,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  gpsi200=adum/(pmarkitmp*mtoroidal)

  call MPI_ALLREDUCE(b2m00,adum,mpsi+1,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  b2m00=adum/(pmarkitmp*mtoroidal)

  allocate(rhom(mgrid),STAT=mtest)
  if (mtest /= 0) then
     write(0,*)mype,'*** Cannot allocate rhom: mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
  rhom=1.0_lk
!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=0,mtheta(i)
        ij=igrid(i)+j
        rhom(ij)=sqrt(aion*meshti(i))/(qion*bmesh(ij))
     enddo
  enddo

end subroutine setgeom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=============================================================================
  Subroutine read_input_params
!=============================================================================

  use global_parameters
  use particle_array
  use field_array
  use particle_tracking
  use precision
  implicit none

  logical file_exist
  integer ierror,micell,mecell,mfcell,mfecell

  namelist /input_parameters/ irun,mstep,msnap,ndiag,toroidaln,nhybrid,nonlinear,paranl,&
      tstep,micell,mecell,mpsi,mthetamax,mtoroidal,ncyclee,psi0,psi1,&
      aion,qion,aelectron,qelectron,ngyroi,nbound,nboundR,iload,r0,b0,etemp0,eden0,&
      eload,etrap,icoll,ecoll,track_particles,ndata3d,magnetic,numereq,nfilter,&
      mfcell,afast,qfast,ngyrof,fload,neop,neot,neoz,ilaplacian,eqcurrent,&
      mfecell,afaste,qfaste,ngyrofe,feload,fetrap,ncyclefe,& 
      fieldmodel,bcond,fielddir,hypr1,hypr2,antenna,izonal,&
      n_modes,m_modes,island,iflux,omega_antenna,eta,fem,ismooth,inorm,sd_v0,sd_vc,&
      sd_l0,sd_widthInv,irestore,irotation, sdp_output

! Since it is preferable to have only one MPI process reading the input file,
! we choose the master process to set the default run parameters and to read
! the input file. The parameters will then be broadcast to the other processes.

  
  if(mype==0) then

    !default values for recently added parameters. 
    ilaplacian=1  ! 0: taking matrix element from iterative solver; 1: matrix element from finite difference
    eqcurrent=0   ! eqcurrent=0: drop curl B terms; eqcurrent=1: keep curl B terms
    fieldmodel=1  ! Analytic eq field model: 0: s-alpha like (cyclone) model; 1: first order (in r/R_0) model with parallel current
    bcond=0       ! 0: fixed zero boundary;  1:linear inner boundary in r 
    fielddir=0    ! Equilibrium magnetic field and current direction
    hypr1=0.0	  ! Parallel hyperviscosity
    hypr2=0.0     ! Perpendicular hyperviscosity 
    island=0      ! Magnetic island perturbation strength. if there is no island, alp=0
    antenna=0	  ! 0: no antenna
    izonal=2
    n_modes=(/1,1,1,1,1,1,1,1/)
    m_modes=(/1,2,3,4,5,6,7,8/)
    iflux=50
    ismooth=1
    fem=0
    sd_v0=0.00299123
    sd_vc=1.07*0.00299123
    sd_l0=0.5
    sd_widthInv=0
    nboundR=0
    fetrap=2  !all fast electrons should be loaded for the correct possion equation
    eta=0.0        ! resistivity, unit: Ohm*m. negative value: without tearing mode; others: with tearing mode.
    inorm=1
    irestore=1
    irotation=0
    toroidaln=1

! table below describes direction in the right poloidal plane
    ! fielddir  B_t  J
    ! 0         out  out
    ! 1         in   out
    ! 2         in   in
    ! 3         out  in

! Test if the input file gtc.in exists
    inquire(file='gtc.in',exist=file_exist)
    if (file_exist) then
      open(55,file='gtc.in',status='old')
      read(55,nml=input_parameters,iostat=ierror)
      if(ierror > 0)write(*,*)'WARNING: check input_parameters'
      close(55)
    else
      write(*,*)'Cannot find file gtc.in !!!'
      stop
    endif
! treat tearing mode

    if(eta .gt. machineEpsilon)then
       write(*,*)"Warning: set nhybrid = 0 to treat tearing mode"
       nhybrid=0
    endif

    mstep=max(ndiag,mstep)
    msnap=min(msnap,mstep/ndiag)

    if(magnetic==0)eqcurrent=0
    if(nonlinear==0)then
       paranl=0.0_lk
       izonal=0
    endif
    if(iload==9)then
       ngyroi=1
       ilaplacian=1
    endif
    if(nhybrid>0)ngyroe=1
#ifdef _CYLINDER
    fieldmodel=0
#endif

    if(psi0==0.0_lk .and. fem==0)bcond=1

! write parameters to standard output
    write(gtcout,input_parameters)

! read data for table lookup for Maxwell function
    if(icoll+ecoll >0)then
       inquire(file='maxwell.dat',exist=file_exist)
       if (file_exist) then
          open(66,file='maxwell.dat',status='old')
          read(66,*)maxwell
          close(66)
       else
          write(*,*)'Cannot find file maxwell.dat !!!'
          stop
       endif
    endif

    ! approximate # of particle per cell
    mi=micell
    me=mecell
    mf=mfcell
    mfe=mfecell  
 endif

! Now send the parameter values to all the other MPI processes
  call broadcast_input_params

end Subroutine read_input_params

!=============================================================================
  Subroutine broadcast_input_params
!=============================================================================

  use global_parameters
  use particle_array
  use particle_tracking
  use field_array,only: iflux
  implicit none

  integer,parameter :: n_integers=53+num_modes+num_modes, n_reals=25
  integer  :: ierror,integer_params(n_integers)
  real(lk) :: real_params(n_reals)

! The master process, mype=0, holds all the input parameters. We need
! to broadcast their values to the other processes. Instead of issuing
! an expensive MPI_BCAST() for each parameter, it is better to pack
! everything in a single vector, broadcast it, and unpack it.

! Pack all parameters in integer_params & real_params
  if(mype==0)then
    integer_params(1)=irun
    integer_params(2)=mstep
    integer_params(3)=msnap
    integer_params(4)=ndiag
    integer_params(5)=nhybrid
    integer_params(6)=izonal
    integer_params(7)=mi
    integer_params(8)=me
    integer_params(9)=mpsi
    integer_params(10)=mthetamax
    integer_params(11)=mtoroidal
    integer_params(12)=npartdom
    integer_params(13)=ncyclee
    integer_params(14)=ngyroi
    integer_params(15)=nbound
    integer_params(16)=iload
    integer_params(17)=track_particles
    integer_params(18)=ndata3d
    integer_params(19)=magnetic
    integer_params(20)=nonlinear
    integer_params(21)=numereq
    integer_params(22)=nfilter
    integer_params(23)=eload
    integer_params(24)=icoll
    integer_params(25)=ecoll
    integer_params(26)=neop
    integer_params(27)=neot
    integer_params(28)=neoz
    integer_params(29)=mf
    integer_params(30)=ngyrof
    integer_params(31)=fload
    integer_params(32)=etrap
    integer_params(33)=ilaplacian
    integer_params(34)=eqcurrent
    integer_params(35)=fieldmodel   
    integer_params(36)=bcond
    integer_params(37)=fielddir
    integer_params(38)=antenna
    integer_params(39)=ngyroe
    integer_params(40)=mfe
    integer_params(41)=feload
    integer_params(42)=ngyrofe
    integer_params(43)=fetrap
    integer_params(44)=ncyclefe
    integer_params(45)=ismooth
    integer_params(46)=fem
    integer_params(47)=nboundR
    integer_params(48)=irestore
    integer_params(49)=island
    integer_params(50)=inorm
    integer_params(51)=irotation
    integer_params(52)=iflux
    integer_params(53)=toroidaln
    integer_params(54:53+num_modes)=n_modes(1:num_modes)
    integer_params(54+num_modes:53+num_modes+num_modes)=m_modes(1:num_modes)

    real_params(1)=paranl
    real_params(2)=tstep
    real_params(3)=rg0
    real_params(4)=psi0
    real_params(5)=psi1
    real_params(6)=aion
    real_params(7)=qion
    real_params(8)=aelectron
    real_params(9)=qelectron
    real_params(10)=r0
    real_params(11)=b0
    real_params(12)=etemp0
    real_params(13)=eden0
    real_params(14)=afast
    real_params(15)=qfast
    real_params(16)=hypr1
    real_params(17)=hypr2
    real_params(18)=afaste 
    real_params(19)=qfaste
    real_params(20)=omega_antenna
    real_params(21)=eta
    real_params(22)=sd_v0
    real_params(23)=sd_vc
    real_params(24)=sd_l0
    real_params(25)=sd_widthInv
  endif

! Send input parameters to all processes
  call MPI_BCAST(integer_params,n_integers,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
  call MPI_BCAST(real_params,n_reals,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

! Broadcast maxwell function data
  call MPI_BCAST(maxwell,100001,mpi_Rsize,0,MPI_COMM_WORLD,ierror)

!   Unpack all parameters
  if(mype/=0)then
    irun=integer_params(1)
    mstep=integer_params(2)
    msnap=integer_params(3)
    ndiag=integer_params(4)
    nhybrid=integer_params(5)
    izonal=integer_params(6)
    mi=integer_params(7)
    me=integer_params(8)
    mpsi=integer_params(9)
    mthetamax=integer_params(10)
    mtoroidal=integer_params(11)
    npartdom=integer_params(12)
    ncyclee=integer_params(13)
    ngyroi=integer_params(14)
    nbound=integer_params(15)
    iload=integer_params(16)
    track_particles=integer_params(17)
    ndata3d=integer_params(18)
    magnetic=integer_params(19)
    nonlinear=integer_params(20)
    numereq=integer_params(21)
    nfilter=integer_params(22)
    eload=integer_params(23)
    icoll=integer_params(24)
    ecoll=integer_params(25)
    neop=integer_params(26)
    neot=integer_params(27)
    neoz=integer_params(28)
    mf=integer_params(29)
    ngyrof=integer_params(30)
    fload=integer_params(31)
    etrap=integer_params(32)
    ilaplacian=integer_params(33)
    eqcurrent=integer_params(34)
    fieldmodel=integer_params(35)
    bcond=integer_params(36)
    fielddir=integer_params(37)
    antenna=integer_params(38)
    ngyroe=integer_params(39)
    mfe=integer_params(40)
    feload=integer_params(41)
    ngyrofe=integer_params(42)
    fetrap=integer_params(43)
    ncyclefe=integer_params(44)
    ismooth=integer_params(45)
    fem=integer_params(46)
    nboundR=integer_params(47)
    irestore=integer_params(48)
    island=integer_params(49)
    inorm=integer_params(50)
    irotation=integer_params(51)
    iflux=integer_params(52)
    toroidaln=integer_params(53)
    n_modes(1:num_modes)=integer_params(54:53+num_modes)
    m_modes(1:num_modes)=integer_params(54+num_modes:53+num_modes+num_modes)

    paranl=real_params(1)
    tstep=real_params(2)
    rg0=real_params(3)
    psi0=real_params(4)
    psi1=real_params(5)
    aion=real_params(6)
    qion=real_params(7)
    aelectron=real_params(8)
    qelectron=real_params(9)
    r0=real_params(10)
    b0=real_params(11)
    etemp0=real_params(12)
    eden0=real_params(13)
    afast=real_params(14)
    qfast=real_params(15)
    hypr1=real_params(16)
    hypr2=real_params(17)
    afaste=real_params(18)  
    qfaste=real_params(19)
    omega_antenna=real_params(20)
    eta=real_params(21)
    sd_v0=real_params(22)
    sd_vc=real_params(23)
    sd_l0=real_params(24)
    sd_widthInv=real_params(25)
  endif

end subroutine broadcast_input_params

!=============================================================================
    Subroutine set_particle_decomp
!=============================================================================

  use global_parameters
  use field_array
  implicit none

  integer  :: i,j,k,pe_number,mtest,ierror
  logical file_exist 


! ----- First we verify the consistency of mtoroidal and npartdom -------
! The number of toroidal domains (mtoroidal) times the number of particle
! "domains" (npartdom) needs to be equal to the number of processor "numberpe".
! numberpe cannot be changed since it is given on the command line.

! numberpe must be a multiple of mtoroidal
  if(mod(numberpe,mtoroidal) /= 0)then
    write(gtcout,*)'Wrong PE #; PE=',numberpe,'Mtoroidal=',mtoroidal
    stop
  endif

!number of particle decomposition.
  npartdom=numberpe/mtoroidal
  if(mype==0)then
    write(gtcout,*)'*******************************************************'
    write(gtcout,*)'  Using npartdom=',npartdom,' and mtoroidal=',mtoroidal
    write(gtcout,*)'*******************************************************'
    write(gtcout,*)
  endif

! We now give each PE (task) a unique domain identified by 2 numbers: the
! particle and toroidal domain numbers.
!    particle_domain_location = rank of the particle domain holding mype
!    toroidal_domain_location = rank of the toroidal domain holding mype
!
! On the IBM SP, the MPI tasks are distributed in an orderly fashion to each
! node unless the LoadLeveler instruction "#@ blocking = unlimited" is used.
! On Seaborg for example, the first 16 tasks (mype=0-15) will be assigned to
! the first node that has been allocated to the job, then the next 16
! (mype=16-31) will be assigned to the second node, etc. When not using the
! OpenMP, we want the particle domains to sit on the same node because
! communication is more intensive. To achieve this, successive PE numbers are
! assigned to the particle domains first.
! It is easy to achieve this ordering by simply using mype/npartdom for
! the toroidal domain and mod(mype,npartdom) for the particle domain.
!
!  pe_number=0
!  do j=0,mtoroidal-1
!     do i=0,npartdom-1
!        pe_grid(i,j)=pe_number
!        particle_domain_location(pe_number)=i
!        toroidal_domain_location(pe_number)=j
!        pe_number=pe_number+1
!     enddo
!  enddo

  particle_domain_location=mod(mype,npartdom)
  toroidal_domain_location=mype/npartdom

!  write(0,*)'mype=',mype,"  particle_domain_location =",&
!            particle_domain_location,' toroidal_domain_location =',&
!            toroidal_domain_location,' pi=',pi

! Domain decomposition in toroidal direction.
  if(toroidaln<1)then
     write(gtcout,*)'Please check gtc.in for toroidaln (must be an integer equal to or larger than 1)!!!',&
                    'You put in toroidaln=',toroidaln,' but now it is toroidaln=1.'
     toroidaln=1
  endif     
  torbound=real(pi2)/real(toroidaln)

  zeta0=torbound*real(toroidal_domain_location)/real(mtoroidal)
  zeta1=torbound*real(toroidal_domain_location+1)/real(mtoroidal)

! grid spacing in the toroidal direction, one toroidal grip per toroidal domain
  deltaz=zeta1-zeta0

  if(mype==0)then
     write(gtcout,*)'Toroidal domain : [0,2pi/',toroidaln,'], deltaz=',deltaz
  endif

! ---- Create particle domain communicator and toroidal communicator -----
! We now need to create a new communicator which will include only the
! processes located in the same toroidal domain. The particles inside
! each toroidal domain are divided equally between "npartdom" processes.
! Each one of these processes will do a charge deposition on a copy of
! the same grid, requiring a toroidal-domain-wide reduction after the
! deposition. The new communicator will allow the reduction to be done
! only between those processes located in the same toroidal domain.
!
! We also need to create a purely toroidal communicator so that the
! particles with the same particle domain id can exchange with their
! toroidal neighbors.
!
! Form 2 subcommunicators: one that includes all the processes located in
! the same toroidal domain (partd_comm), and one that includes all the
! processes part of the same particle domain (toroidal_comm).
! Here is how to create a new communicator from an old one by using
! the MPI call "MPI_COMM_SPLIT()".
! All the processes passing the same value of "color" will be placed in
! the same communicator. The "rank_in_new_comm" value will be used to
! set the rank of that process on the communicator.
!  call MPI_COMM_SPLIT(old_comm,color,rank_in_new_comm,new_comm,ierror)

! particle domain communicator (for communications between the particle
! domains WITHIN the same toroidal domain)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,toroidal_domain_location,&
                      particle_domain_location,partd_comm,ierror)

! toroidal communicator (for communications BETWEEN toroidal domains of same
! particle domain number)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD,particle_domain_location,&
                      toroidal_domain_location,toroidal_comm,ierror)

  call mpi_comm_size(partd_comm,nproc_partd,ierror)
  call mpi_comm_rank(partd_comm,myrank_partd,ierror)

  call mpi_comm_size(toroidal_comm,nproc_toroidal,ierror)
  call mpi_comm_rank(toroidal_comm,myrank_toroidal,ierror)

!  write(0,*)'mype=',mype,'  nproc_toroidal=',nproc_toroidal,&
!       ' myrank_toroidal=',myrank_toroidal,'  nproc_partd=',nproc_partd,&
!       ' myrank_partd=',myrank_partd

  if(nproc_partd/=npartdom)then
    write(0,*)'*** nproc_partd=',nproc_partd,' NOT EQUAL to npartdom=',npartdom
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  if(nproc_toroidal/=mtoroidal)then
    write(*,*)'*** nproc_toroidal=',nproc_toroidal,' NOT EQUAL to mtoroidal=',&
              mtoroidal
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! We now find the toroidal neighbors of the current toroidal domain and
! store that information in 2 easily accessible variables. This information
! is needed several times inside the code, such as when particles cross
! the domain boundaries. We will use the toroidal communicator to do these
! transfers so we don't need to worry about the value of myrank_partd.
! We have periodic boundary conditions in the toroidal direction so the
! neighbor to the left of myrank_toroidal=0 is (mtoroidal-1).

  left_pe=mod(myrank_toroidal-1+mtoroidal,mtoroidal)
  right_pe=mod(myrank_toroidal+1,mtoroidal)

end subroutine set_particle_decomp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initial(cdate,timecpu)
  use global_parameters
  use particle_array
  implicit none

  integer ierror,i,m,nsize
  real(doubleprec) timecpu(10)
  integer,dimension(:),allocatable :: mget,mput
  character(len=16) cdate(4)
#ifdef _OPENACC
  integer :: numDevice,myDevice,ierr
#endif


! MPI initialize, total # of PE, and rank of PE
  call mpi_init(ierror)
  call mpi_comm_size(mpi_comm_world,numberpe,ierror)
  call mpi_comm_rank(mpi_comm_world,mype,ierror)

! open standard output file, record program starting time
  if(mype==0)then
     open(gtcout,file='gtc.out',form='formatted',status='replace')
     call date_and_time(cdate(1),cdate(2))
     write(gtcout,*)"GTC Version 3"
     write(gtcout,*)"Compiled on ",cdate(1)," at ",cdate(2)

#ifdef _FRC
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Compiled for field-reversed configuration equilibrium.'
     write(gtcout,'("===================================",/)')
#elif _TOKAMAK
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Compiled for tokamak equilibrium.'
     write(gtcout,'("===================================",/)')
#elif _CYLINDER
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Compiled for cylinder equilibrium.'
     write(gtcout,'("===================================",/)')
#elif _TOROIDAL3D
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Compiled for 3-D equilibrium.'
     write(gtcout,'("===================================",/)')
#else
     write(gtcout,*)'Please compile with a specific configuration for ICONFIG in Makefile.'
     stop
#endif

#ifdef _PETSc
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Compiled with PETSc.'
     write(gtcout,'("===================================",/)')
#else 
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Using Iterative Solver.'
     write(gtcout,'("===================================",/)')
#endif

! check OPENMP
#ifdef _OPENMP
!$omp parallel
!$omp master
     nthreads=omp_get_num_threads()
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Number of OpenMP threads = ',nthreads
     write(gtcout,'("===================================",/)')
!$omp end master
!$omp end parallel
#else
     write(gtcout,'(/,"===================================")')
     write(gtcout,*)' Run without OpenMP threads'
     write(gtcout,'("===================================",/)')
#endif
  endif

! initialize GPU
#ifdef _OPENACC
  call acc_init(acc_device_nvidia)
  numDevice=acc_get_num_devices(acc_device_nvidia)
  if(numDevice/=0)then
    myDevice=mod(mype,numDevice)
    call acc_set_device_num(myDevice,acc_device_nvidia)
    myDevice=acc_get_device_num(acc_device_nvidia)
    ! print *, 'pe=',mype,'using ',idevice
    !call setupgpu(1)
    if(mype==0)then
      write(gtcout,'(/,"===================================")')
      write(gtcout,*)' Number of OpenACC Devices = ',numDevice 
      write(gtcout,'("===================================",/)')
    endif
  else
    call acc_set_device_type(acc_device_host)
    if(mype==0)then
      write(gtcout,'(/,"===================================")')
      write(gtcout,*)' Run without OpenACC Devices'
      write(gtcout,'("===================================",/)')
    endif
  endif
#endif

! initiate adios
#if ADIOS
  CALL adios_init ("config.xml"//char(0), MPI_COMM_WORLD, MPI_COMM_SELF, MPI_INFO_NULL)
#endif

!!XY rolling restart, initialize restart directory and control var
  irest=0
  FileExit=0

! initialize timer
  istep=0
  timecpu=0.0
  timecpu(9)=mpi_wtime()
  timecpu(10)=timecpu(9)

! numerical constant
  pi=4.0_lk*atan(1.0_lk)
  pi2=2.0_lk*pi
  pi2_inv=1.0_lk/pi2

! **** Use the intrinsic F90 random number generator *****
! initialize f90 random number generator
!    call random_seed !comemt out this line for identitical random number
  call random_seed(size=nsize)
  allocate(mget(nsize),mput(nsize))
  call random_seed(get=mget)
  do i=1,nsize
     !TODO: system_clock causes problem on Titan with pgi compiler

     !call system_clock(m)   !random initialization for collision

     if(irun==0)m=1         !same initialization
     mput(i)=111111*(mype+1)+m+mget(i)
  enddo
  call random_seed(put=mput)
  deallocate(mget,mput)

! Read the input file that contains the run parameters
  call read_input_params

! Set up the particle decomposition within each toroidal domain
  call set_particle_decomp

! PETSc initialize
#ifdef _PETSc
call petsc_init
#endif

end subroutine initial
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine final(cdate,timecpu)
  use global_parameters
  implicit none

  integer ierror,iorestart,i
  real(doubleprec) timecpu(10)
  character(len=16) ic(9)
  character(len=16) cdate(4)

!total cpu and wall clock time
  timecpu(10)=mpi_wtime()
  timecpu(9)=timecpu(10)-timecpu(9)
  ic(1)='pusher:'
  ic(2)='shift:'
  ic(3)='charge:'
  ic(4)='electron:'
  ic(5)='fast:'
  ic(6)='Poisson:'
  ic(7)='initial:'
  ic(8)='IO:'
  ic(9)='total:'
  if(mype==0)then
     write(gtcout,*)'CPU TIME USAGE (in SEC):'
     !write(gtcout,*)ic
     !write(gtcout,'(8(1pe8.1),/)')timecpu(1:8)
     write(gtcout,'(A16,1pe12.5,1pe16.5)')&
        (ic(i),timecpu(i),timecpu(i)/timecpu(9), i=1,9)

! Restart file info
     FileExit=1
     iorestart=345
     open(iorestart,file="FileExit.out",status="replace")
     write(iorestart,"(A9,i1)")"FileExit=",FileExit
     write(iorestart,"(A9,i5)")"irest   =",irest
     if(mod(irest+1,2)==0)then
        write(iorestart,"(A12)")"restart_dir1"
     else
        write(iorestart,"(A12)")"restart_dir2"
     endif
     close(iorestart)

! record program end time
     call date_and_time(cdate(3),cdate(4))
     write(gtcout,*)'Program starts at DATE=', cdate(1), 'TIME=', cdate(2)
     write(gtcout,*)'Program ends at   DATE=', cdate(3), 'TIME=', cdate(4)

! close standard output file
     if(gtcout /= 6 .and. gtcout /= 0)close(gtcout)
  endif

! PETSc finalize
#ifdef _PETSc
#ifdef _FRC
     call petsc2_final
#else
     call petsc_final
#endif
#endif

#if ADIOS
  CALL adios_finalize (mype)
#endif

! MPI finalize
  call mpi_finalize(ierror)

end subroutine final
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!=========================================
subroutine timer(t0,timecpu)
!=========================================
  use precision
  implicit none
  real(doubleprec) t0,t1,dt,timecpu

! Get cpu usage time since the begine of run and subtract value from the previous call
  t1=mpi_wtime()
  dt=t1-t0
  t0=t1
  timecpu=timecpu+dt

! Get wall clock time and subtract value from the previous call
!  t1wc=mpi_wtime()
!  dtwc=t1wc-t0wc
!  t0wc=t1wc
!  timewc=timewc+dtwc

end subroutine timer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine polsmooth(farray,ip)
  use precision
  use global_parameters,only:mthetamax,mgrid,mpsi
  use field_array,only:mtheta,igrid
  implicit none

  integer i,j,ii,jt,ip
  real(lk) pright(mthetamax),pright0(mthetamax),pright1(mthetamax),&
       pright2(mthetamax),pright3(mthetamax),farray(mgrid)
  if(ip==1)then
! poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
    do j=1,4
!$omp parallel do private(i)
       do i=1,mpsi-1
          farray(igrid(i))=farray(igrid(i)+mtheta(i)) ! poloidal BC
       enddo

!$omp parallel do private(i,ii,jt,pright,pright0,pright1,pright2,pright3)
       do i=0,mpsi
          ii=igrid(i)
          jt=mtheta(i)
          pright(1:jt)=farray(ii+1:ii+jt)
          pright0(1:jt)=cshift(pright(1:jt),-1)
          pright1(1:jt)=cshift(pright(1:jt),1)
          pright2(1:jt)=cshift(pright0(1:jt),-1)
          pright3(1:jt)=cshift(pright1(1:jt),1)
          farray(ii+1:ii+jt)=0.20_lk*pright(1:jt)+0.20_lk*(pright0(1:jt)+pright1(1:jt))+&
                  0.20_lk*(pright2(1:jt)+pright3(1:jt))
       enddo
    enddo
  endif
!$omp parallel do private(i,ii,jt)
  do i=0,mpsi
     ii=igrid(i)
     jt=mtheta(i)
     farray(ii)=farray(ii+jt)
  enddo

end subroutine polsmooth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! if fload>0, fast ion mass density is taken into account but not fast ion pressure
! in that case this calculation may not make sense
subroutine alcondata
  use global_parameters
  use equilibrium
  use particle_array,only:aion,afast,betae,fload!,afaste,feload
  use spline_function, only:sptorpsi,spq,spgpsi,spcpsi,spb,dbdt,dxdp,dxdt,dzdp,dzdt 
  implicit none

! gamma for ions and electrons
  real(lk),parameter :: gammai=7.0_lk/4.0_lk, gammae=1.0_lk

! number of radial grid points
  integer,parameter :: nrad=200

! ratio of psi at first radial grid point to ped
  real(lk),parameter :: psifrac=1e-6_lk

! number of poloidal (theta) grid points in the sampling for fft
  integer,parameter :: npol=2560

! m-harmonics from m=0 to m=nfft will be written out for the equilibrium quantities
! (should be much smaller than npol)
  integer,parameter :: nfft=50

! internal parameters, no need to change in most cases
! file number for writing out continua data
  integer,parameter :: iacd=240
  integer,parameter :: nprofiledata=5,nfftdata=5
  real(lk) acdprofile(nprofiledata,nrad)
  complex(lk) acdfft(0:nfft,nfftdata,nrad)
  real(lk),dimension(0:npol-1) :: hsamp,jsamp,ksamp,lsamp,nsamp
  complex(lk),dimension(0:npol/2) :: hsampfft,jsampfft,ksampfft,lsampfft,nsampfft
  integer i,j,isp
  real(lk) psi,dpx,dp2
  real(lk) theta
  real(lk) ne,ni,nf,te,ti,tf,nfe,tfe
  real(lk) dpdx,dpdz
  real(lk) b,g

  open(iacd,file='alcon.dat',status='replace')
  write(iacd,*)nrad,nfft,nprofiledata,nfftdata

  acdprofile=0.0_lk
  acdfft=0.0_lk
  do i=1,nrad
!   uniform in psi, very not uniform in rho
!    psi=ped*(psifrac+(1.0_lk-psifrac)*(i-1)/(nrad-1))

!   nonuniform in psi, more uniform in rho
    psi=sqrt(psifrac)+(1.0_lk-sqrt(psifrac))*(i-1)/(nrad-1)
    psi=psi*psi*ped

    isp=max(1,min(lsp-1,ceiling(psi*spdpsi_inv)))
    dpx=psi-spdpsi*real(isp-1)
    dp2=dpx*dpx

    acdprofile(1,i)=sqrt(sptorpsi(psi)/sptorpsi(ped))

    acdprofile(2,i)=abs(spq(psi))
    g=abs(spgpsi(psi))
    acdprofile(3,i)=g*acdprofile(2,i)+spcpsi(psi)
    ne=         nepp(1,isp)   +nepp(2,isp)*dpx   +nepp(3,isp)*dp2
    ni=         nipp(1,isp)   +nipp(2,isp)*dpx   +nipp(3,isp)*dp2
    nf=         nfpp(1,isp)   +nfpp(2,isp)*dpx   +nfpp(3,isp)*dp2
   !nfe=        nfepp(1,isp)  +nfepp(2,isp)*dpx  +nfepp(3,isp)*dp2
    te=         tepp(1,isp)   +tepp(2,isp)*dpx   +tepp(3,isp)*dp2
    ti=         tipp(1,isp)   +tipp(2,isp)*dpx   +tipp(3,isp)*dp2
    tf=         tfpp(1,isp)   +tfpp(2,isp)*dpx   +tfpp(3,isp)*dp2
   !tfe=        tfepp(1,isp)  +tfepp(2,isp)*dpx  +tfepp(3,isp)*dp2
    acdprofile(4,i)=betae/(2.0_lk*rho0*rho0)*(gammai*ni*ti+gammae*ne*te)
    acdprofile(5,i)=aion*ni
    if(fload>0)acdprofile(5,i)=acdprofile(5,i)+afast*nf
   !if(feload>0)acdprofile(5,i)=acdprofile(5,i)+afaste*nfe

    do j=0,npol-1
      theta=pi2/real(npol)*real(j)
      b=spb(psi,theta)
      dpdx=1.0_lk/(dxdp(psi,theta)-dxdt(psi,theta)*dzdp(psi,theta)/dzdt(psi,theta))
      dpdz=1.0_lk/(dzdp(psi,theta)-dzdt(psi,theta)*dxdp(psi,theta)/dxdt(psi,theta))
      hsamp(j)=(dpdx*dpdx+dpdz*dpdz)/acdprofile(3,i)
      jsamp(j)=(dpdx*dpdx+dpdz*dpdz)*acdprofile(3,i)/(b**4)
      ksamp(j)=2.0_lk*g*dbdt(psi,theta)/(b**3)
      lsamp(j)=(acdprofile(4,i)+b*b)/(b**4)*acdprofile(3,i)
      nsamp(j)=acdprofile(4,i)*ksamp(j)*ksamp(j)/lsamp(j)
    enddo
    call fftr1d(1,npol,pi2_inv,hsamp,hsampfft,2)
    call fftr1d(1,npol,pi2_inv,jsamp,jsampfft,2)
    call fftr1d(1,npol,pi2_inv,ksamp,ksampfft,2)
    call fftr1d(1,npol,pi2_inv,lsamp,lsampfft,2)
    call fftr1d(1,npol,pi2_inv,nsamp,nsampfft,2)
    do j=0,nfft
      acdfft(j,1,i)=(hsampfft(j))/real(npol)
      acdfft(j,2,i)=(jsampfft(j))/real(npol)
      acdfft(j,3,i)=(ksampfft(j))/real(npol)
      acdfft(j,4,i)=(lsampfft(j))/real(npol)
      acdfft(j,5,i)=(nsampfft(j))/real(npol)
    enddo
  enddo

  write(iacd,*)acdprofile
  write(iacd,*)acdfft
  close(iacd)
end subroutine alcondata
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!initialize the external phi with non-zero boundary
subroutine phiext_initial
  use global_parameters
  use equilibrium
  use field_array
  use particle_array,only:betae,aion,qion
  implicit none

  integer i,j,ij,ia,win0,win1,m,n,ierror,mtest
  real(lk) adum,envelope,tdum,qmin
  real(8) temp,dbesjn

  interface
    subroutine laplacian(scalar,scalar_out,none0bound)
      use precision
      use global_parameters,only:mgrid
      implicit none

      real(lk),dimension(mgrid),intent(in) :: scalar
      real(lk),dimension(mgrid),intent(out) :: scalar_out
      integer, optional, intent (in) :: none0bound
    end subroutine laplacian
  end interface
  
  allocate(phiext(0:antenna,mgrid),omega(antenna),dn_ext(0:antenna,mgrid),STAT=mtest)
  if(mtest /=0) then
     write(0,*)mype,'*** Cannot allocate phiext:mtest=',mtest
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  phiext=0.0
  omega=0.0

  do ia=1,antenna
    !if(ia==1)then
         n=n_modes(ia)
         m=m_modes(ia)
    !endif

! radial domain applying external source
!by zhs use old antenna profile
     win0=int(mpsi/4)
     win1=int(3*mpsi/4)
     adum=pi/real(win1-win0)
     omega(ia)=omega_antenna !set your own frequency

     if(mype==0)then
        write(gtcout,*)'Antenna=',antenna
        write(gtcout,*)'Antenna',ia,' freq=',omega(ia),'=',omega(ia)/(pi2*utime),'Hz'
        call FLUSH(gtcout)
     endif

!~~~~~~~~~~~antenna profile~~~~~~~~~~~~~~~~~~~~~
     do i=0,mpsi
         envelope=1.0e-5*(psimesh(i)-psimesh(0))/(psimesh(mpsi)-psimesh(0))
         !rg=sqrt((psimesh(i)-psimesh(0))/(psimesh(mpsi)-psimesh(0)))
         !temp=rg*38.1599
         !envelope=1.e6*dbesjn(m,temp)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!$omp parallel do private(j,ij,tdum)
        do j=1,mtheta(i)
           ij=igrid(i)+j
#ifdef _FRC
             tdum=deltat(i)*real(j)
#else
             tdum=deltat(i)*real(j)+zeta1*qtinv(i)
#endif
           phiext(ia,ij)=envelope*cos(real(m)*tdum-real(n)*zeta1)
        enddo
        phiext(ia,igrid(i))=phiext(ia,igrid(i)+mtheta(i))
     enddo


!The third argurment, 1, here indicates that phiext generally has non-zero boundary
     call laplacian(phiext(ia,:),dn_ext(ia,:), 1)

  enddo
end subroutine phiext_initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!analyticmesh provides temporary solution for inaccuracy of sqrt(x) spline functions near origin 
subroutine analyticmesh
       
  use global_parameters
  use equilibrium,only:psiw,ped
  use field_array
  use spline_function, only: sptorpsi
  implicit none

  logical file_exist
  integer i,ii
  real(lk) pleft,pright,psi,rdum,tordum,delr,sperr,rg,qion,qelectron,&
       q(3),ti(3),tf(3),tfe(3),te(3),ne(3),ni(3),nf(3),nfe(3),ze(3),itemp0,ftemp0,fden0,fetemp0,feden0,er(3),&
       psiw_analytic,ped_analytic,q_analytic(3),&
       ze_analytic(3),er_analytic(3),itemp0_analytic,ftemp0_analytic,&
       fden0_analytic,ne_analytic(3),te_analytic(3),ti_analytic(3),&
       tf_analytic(3),nf_analytic(3),fetemp0_analytic,feden0_analytic,&
       tfe_analytic(3),nfe_analytic(3)

  namelist /equilibrium_parameters/ psiw_analytic,ped_analytic,q_analytic,&
       ze_analytic,er_analytic,itemp0_analytic,ftemp0_analytic,&
       fden0_analytic,ne_analytic,te_analytic,ti_analytic,&
       tf_analytic,nf_analytic,fetemp0_analytic,feden0_analytic,&
       tfe_analytic,nfe_analytic

! Test if the input file gtc.in exists
  inquire(file='gtc.in',exist=file_exist)
  if (file_exist) then
    open(55,file='gtc.in',status='old')
    read(55,nml=equilibrium_parameters)
    close(55)
  else
    write(*,*)'Cannot find file gtc.in !!!'
    stop
  endif
  psiw=psiw_analytic
  ped=ped_analytic
  q(1:3)=q_analytic(1:3)
  ze(1:3)=ze_analytic(1:3)
  er(1:3)=er_analytic(1:3)
  itemp0=itemp0_analytic
  ftemp0=ftemp0_analytic
  fden0=fden0_analytic
  fetemp0=fetemp0_analytic
  feden0=feden0_analytic
  ne(1:3)=ne_analytic(1:3)
  te(1:3)=te_analytic(1:3)
  ti(1:3)=ti_analytic(1:3)
  tf(1:3)=tf_analytic(1:3)
  nf(1:3)=nf_analytic(1:3)
  tfe(1:3)=tfe_analytic(1:3)
  nfe(1:3)=nfe_analytic(1:3)

  qmesh(0)=q(1)+q(2)*psi0+q(3)*psi0*psi0
  qmesh(mpsi)=q(1)+q(2)*psi1+q(3)*psi1*psi1

  psi0=psi0*ped
  psi1=psi1*ped

  psimesh(0)=psi0
  psimesh(mpsi)=psi1
  tormesh(0)=sptorpsi(psi0)
  tormesh(mpsi)=sptorpsi(psi1)

  tordum=q(1)*psi0+0.5*q(2)*psi0*psi0/ped+q(3)*psi0*psi0*psi0/(ped**2)/3.0
  rdum=sqrt(2.0*tordum)
  rg0=rdum
  tordum=q(1)*psi1+0.5*q(2)*psi1*psi1/ped+q(3)*psi1*psi1*psi1/(ped**2)/3.0
  rdum=sqrt(2.0*tordum)
  deltar=(rdum-rg0)/real(mpsi)
  sperr=0.0000001

  do i=1,mpsi-1
     rg=rg0+deltar*real(i) !target rg
     pleft=psimesh(i-1)
     pright=psimesh(mpsi)
     psi=0.5*(pleft+pright) !initial guess of psi
     tordum=q(1)*psi+0.5*q(2)*psi*psi/ped+q(3)*psi*psi*psi/(ped**2)/3.0
     rdum=sqrt(2.0*tordum)
     delr=1.0-rdum/rg
     ii=1

     do while (abs(delr)>sperr .and. ii<1000)
        ii=ii+1
        if(delr>0)then
           pleft=psi
        elseif(delr<0)then
           pright=psi
        endif
        psi=0.5*(pleft+pright) !new value of psi
	tordum=q(1)*psi+0.5*q(2)*psi*psi/ped+q(3)*psi*psi*psi/(ped*ped)/3.0
     	rdum=sqrt(2.0*tordum)
        delr=1.0-rdum/rg
     enddo
!poloidal, toroidal, q flux on radial mesh
     psimesh(i)=psi
     tormesh(i)=tordum
     qmesh(i)=q(1)+q(2)*psi/ped+q(3)*(psi/ped)*(psi/ped)
  enddo

end subroutine analyticmesh


subroutine createtriangles

   use global_parameters
   use equilibrium
   use field_array
   use spline_function, only: spx, spz
   implicit none 
  
   integer:: i,j,ij,ij0,ij1,ij2,trii1,trii2
   real(lk)::pdum,tdum

!  calculate the number of triangles for system  
!  the are mtheta(i)+mtheta(i+1) for each flux surface    
   trilength=mtheta(1)+mtheta(0)
   if(psi0==0.0_lk)trilength=mtheta(1)
   do i=1,mpsi-1
      trilength=trilength+mtheta(i)+mtheta(i+1)
   enddo
   allocate(trilist(trilength,3),xygrid(mgrid_fem,2))
  
!Define xygrid for FEM solver
   do i=0,mpsi
      pdum=psimesh(i)
      do j=1,mtheta(i)
         ij=igrid_fem(i)+j
         if(fielddir==1 .or. fielddir==3)then
            tdum=modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)
         else
            tdum=modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)
         endif
         xygrid(ij,1)=spx(pdum,tdum)
         xygrid(ij,2)=spz(pdum,tdum)
      enddo
   enddo

!   if(mype==0)then
!      open(129,file='xygrid.txt')
!         do i=0,mpsi
!            do j=1,mtheta(i)
!               ij=igrid_fem(i)+j
!               write(129,*) xygrid(ij,1),xygrid(ij,2),i
!            enddo
!         enddo
!      close(129)
!   endif
  
!FEM triangle grid setup  
   trii2=0
   trilist = -1
   do i=0,mpsi-1
!     triangle algorithm doesn't work when psi0=0.0, add triangles manually
      if(i==0 .and. psi0==0.0_lk)then
         do j=1,mtheta(1)
               trii2=trii2+1
               trilist(trii2,1) = 1
               trilist(trii2,2) = j+1
               trilist(trii2,3) = modulo(j,mtheta(1))+2
            enddo
      else
!	 ij0 to ij1 -> lower flux surface points(i),ij1 to ij2 -> upper flux surface points(i+1)			
         ij0 =igrid_fem(i)+1
         ij1 = igrid_fem(i)+mtheta(i)
         ij2 = igrid_fem(i+1)+mtheta(i+1)
 
!	 calculate triangles position in array			
         trii1=trii2+1
         trii2=trii1+ij2-ij0
!        Solve triangular mesh for flux surfaces i and i+1    		
         call addtriangles(xygrid(ij0:ij1,:),mtheta(i),xygrid(ij1+1:ij2,:),mtheta(i+1),&
         igrid_fem(i),trilist(trii1:trii2,:))
      endif
   enddo

!   if(mype==0)then
!      open(130,file='triangles.txt')
!      do i=1,trilength
!         write(130,*) trilist(i,:)
!      enddo
!      close(130)
!   endif
    
end subroutine createtriangles

subroutine addtriangles(isurface,ilength,osurface,olength,igrid_fem,triangles)

   use precision
   implicit none
   
   integer,intent(in)::ilength,olength,igrid_fem
   real(lk),intent(in)::isurface(ilength,2),osurface(olength,2)
   integer,intent(inout)::triangles(ilength+olength,3)
   integer:: i,i1,j1,i2,j,j2,k,k2,ij,ij2,trii
   real(lk)::neardist,dist,dist1(3),dist2(3)

!  locale triangle number between two flux surfaces   
   trii=1

!  Calculate initial i1 and j1 points to contruct first triangle
   i1=1
   neardist=1.0_lk
   do j=1,olength
      dist = (isurface(i1,1)-osurface(j,1))*(isurface(i1,1)-osurface(j,1)) &
          +(isurface(i1,2)-osurface(j,2))*(isurface(i1,2)-osurface(j,2)) 

      if(dist <= neardist)then
         neardist = dist
         j1=j
      endif
   enddo

!     calculate the next i-point and j-point, for triangle, (i1,j1,i2) or (i1,j1,j2)
!     whichever triangle is better formed, add to triangle list
   do i=1,ilength+olength
      i2=modulo(i1,ilength)+1
      j2=modulo(j1,olength)+1

!     calculate distances between each point of triangle with i2 gridpoint
      dist1(1)=(isurface(i1,1)-osurface(j1,1))*(isurface(i1,1)-osurface(j1,1)) &
            +(isurface(i1,2)-osurface(j1,2))*(isurface(i1,2)-osurface(j1,2))    
      dist1(2)=(isurface(i1,1)-isurface(i2,1))*(isurface(i1,1)-isurface(i2,1)) &
            +(isurface(i1,2)-isurface(i2,2))*(isurface(i1,2)-isurface(i2,2))
      dist1(3)= (isurface(i2,1)-osurface(j1,1))*(isurface(i2,1)-osurface(j1,1)) &
            +(isurface(i2,2)-osurface(j1,2))*(isurface(i2,2)-osurface(j1,2))

!     calculate distances between each point of triangle with i2 gridpoint
      dist2(1)=(isurface(i1,1)-osurface(j1,1))*(isurface(i1,1)-osurface(j1,1)) &
            +(isurface(i1,2)-osurface(j1,2))*(isurface(i1,2)-osurface(j1,2))
      dist2(2)=(isurface(i1,1)-osurface(j2,1))*(isurface(i1,1)-osurface(j2,1)) &
            +(isurface(i1,2)-osurface(j2,2))*(isurface(i1,2)-osurface(j2,2))
      dist2(3)=(osurface(j2,1)-osurface(j1,1))*(osurface(j2,1)-osurface(j1,1)) &
            +(osurface(j2,2)-osurface(j1,2))*(osurface(j2,2)-osurface(j1,2))

!     Choose best triangle. 1,2,3 grid points are written in counter clockwise order
      if(sum(dist1)>=sum(dist2))then
         triangles(trii,1)=i1
         triangles(trii,2)=ilength+j1
         triangles(trii,3)=ilength+j2
         trii = trii+1
         j1=modulo(j1,olength)+1
      else
         triangles(trii,1)=i1
         triangles(trii,2)=ilength+j1
         triangles(trii,3)=i2
         trii = trii+1
         i1=modulo(i1,ilength)+1
      endif
   enddo

!  add igrid to triangles for correct position in total set of grid points
   triangles(:,:)= triangles(:,:)+igrid_fem

end subroutine addtriangles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bminmax
  use global_parameters
  use equilibrium
  use field_array
  use spline_function, only: spb
  implicit none

  integer i,j,nb,ii,ifs
  real(lk) bdum,testb,tdum,pdum,deltab,tdum0
  real(lk),dimension(:),allocatable :: test01,test02

  minbfield(0:1)=100.0
  maxbfield(0:1)=0.0
  nb=mthetamax/2 ! number of b-field bins between bmin and bmax
  allocate(thetaupp(0:1,nb),thetadwn(0:1,nb)) ! poloidal theta indexes for each b-field value in upp and down halves

  do ifs=0,1
! ifs labels outer (1) and inner (0) flux surfaces  
    if(ifs==0) ii=0
    if(ifs==1) ii=mpsi

    pdum=psimesh(ii)
! goes over theta to find minimuma nd maximum b-field (assuming it's monotonic)
    do j=0,mtheta(ii)
       tdum=real(j)*deltat(ii)
       bdum=spb(pdum,tdum)
       if (bdum<minbfield(ifs)) then
          minbfield(ifs)=bdum ! min b-filed value
          thetabmin(ifs)=tdum ! poloidal angle of minimum b-field value
       endif
       if (bdum>maxbfield(ifs)) then 
          maxbfield(ifs)=bdum ! max b-filed value
          thetabmax(ifs)=tdum ! poloidal angle of maximum b-field value
       endif
    enddo
! finds location of same b-field values at the opposite theta half-plane
    deltab=(maxbfield(ifs)-minbfield(ifs))/real(nb-1)
    do i=1,nb
       testb=minbfield(ifs)+deltab*real(i-1)

       bdum=minbfield(ifs)
       j=0
       tdum=thetabmin(ifs)
       do while(bdum<testb)
          j=j+1
          tdum=mod(thetabmin(ifs)+real(j)*deltat(ii),pi2)
          bdum=spb(pdum,tdum)
       enddo
       thetadwn(ifs,i)=floor(tdum/deltat(ii))
  
       bdum=minbfield(ifs)
       j=0
       tdum=thetabmin(ifs)
       do while(bdum<testb)
          j=j+1
          tdum=mod(pi2+thetabmin(ifs)-real(j)*deltat(ii),pi2)
          bdum=spb(pdum,tdum)
       enddo
       thetaupp(ifs,i)=floor(tdum/deltat(ii))
    enddo
    thetabmax(ifs)=mod(pi2+thetabmax(ifs)-thetabmin(ifs),pi2) ! w.r.t. thetabmin
  enddo
end subroutine bminmax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussbound(y)
  use precision
  use global_parameters,only: nbound,mpsi,nboundR
  implicit none
  real(lk) sigma
  real(lk) y(0:mpsi)
  integer i
  
  sigma=((nbound-101)/2)/(2*sqrt(2*alog(2.0))) !magnitude is halved at fourth of nbound domain.FWHM=nbound/2
  do i=0,nbound-101
    y(i)=y(i)*exp(-(((i-(nbound-101))/sigma)**2)/2)
    if(nboundR==0)y(mpsi-i)=y(mpsi-i)*exp(-(((i-(nbound-101))/sigma)**2)/2)
  enddo
end subroutine gaussbound
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gaussboundR(y)
  use precision
  use global_parameters,only: nboundR,mpsi
  implicit none
  real(lk) sigma
  real(lk) y(0:mpsi)
  integer i
  
  sigma=((nboundR-101)/2)/(2*sqrt(2*alog(2.0))) !magnitude is halved at fourth of nboundR domain.FWHM=nboundR/2
  do i=0,nboundR-101
    y(mpsi-i)=y(mpsi-i)*exp(-(((i-(nboundR-101))/sigma)**2)/2)
  enddo
end subroutine gaussboundR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Field-aligned grids of LAMIRIDGE_FRC initialization
subroutine LR_initial
      use precision
      use global_parameters,only:gtcout,mype
      use cylindricalRZ,only:lsr,lsz,psirz_Eq,spdR,spdZ
      use LAMYRIDGEeq
      use spline_function, only: spline2d,sppsiR1d
      implicit none

      integer i,j,ij,ij1

! simulation domain in Z direction
      Zsimulation0=-9.3058_lk
      Zsimulation1=9.3058_lk

! simulation domain in psi direction
      psi1_frc=2.0_lk
      psi0_frc=0.5_lk
      min_r=0.001_lk ! minimal r allowed in the simulation, normalized by xnorm

! simulation grid number along psi and Z direction
      numr=48
      numz=96
      nrzgrid=numr*numz

      write(gtcout,*)"numr =",numr
      write(gtcout,*)"numz =",numz

! grid size along Z direction
      lz=Zsimulation1-Zsimulation0
      delta_Z=lz/real(numz-1)
!create grid index arrays
      allocate(igrid_frc(numr),nz_frc(numr))
      nz_frc(:)=numz  ! the number of the grid on each psi flux
      igrid_frc(1)=0
      do i=2,numr
        igrid_frc(i)=igrid_frc(i-1)+nz_frc(i-1)
      enddo

      write(gtcout,*)"igrid_frc=",igrid_frc(numr)
      write(gtcout,*)"nz_frc =",nz_frc(numr)

! create the array to store the R,Z position of field aligned mesh
      allocate(RZ_reg(nrzgrid,2))
      do i=1,numr
        do j=1,nz_frc(i)
         ij=igrid_frc(i)+j
         RZ_reg(ij,2) = (j-1)*delta_Z + Zsimulation0
        enddo
      enddo

      allocate(psiEq0(nrzgrid),meshni_cy(nrzgrid),meshti_cy(nrzgrid),meshte_cy(nrzgrid),meshne_cy(nrzgrid),psi_frc(numr),psifrc_tmp(nrzgrid),deltar_frc(nrzgrid),dpsifrcdr(nrzgrid),dpsifrcdz(nrzgrid))

! difference of psi between each flux surfaces
      deltapsi_frc=(psi1_frc-psi0_frc)/real(numr-1)

      write(gtcout,*)"deltapsi_frc =",deltapsi_frc
! each flux surface psi value 
      do i=1,numr
        psi_frc(i) = psi0_frc+real(i-1)*deltapsi_frc
      enddo
! find the R position of each grids
!$omp parallel do private(i,j,ij) 
      do i=1,numr
        do j=1,nz_frc(i)
          ij=igrid_frc(i)+j
          call sppsiR1d(0,lsr,lsz,spdR,spdZ,r_Eq,z_Eq,psirz_Eq,RZ_reg(ij,2)+zmax,psi_frc(i),RZ_reg(ij,1),psifrc_tmp(ij))
        enddo
      enddo

!write(gtcout,*)"RZ_reg1",RZ_reg(:,1)
!write(gtcout,*)"RZ_reg2",RZ_reg(:,2)

! calculate the distance between neighbor grid in R direction
!$omp parallel do private(i,j,ij,ij1) 
      do i=1,numr-1
        do j=1,nz_frc(i)
          ij=igrid_frc(i)+j
          ij1=igrid_frc(i+1)+j
          deltar_frc(ij)=RZ_reg(ij1,1)-RZ_reg(ij,1)
        enddo
      enddo

!$omp parallel do private(i,j,ij) 
      do j=1,nz_frc(numr)
         ij=igrid_frc(numr)+j
         deltar_frc(ij)=0.0
      enddo

! caulculate the equilibrium quantites on the field aligned mesh
!$omp parallel do private(i) 
      do i=1,nrzgrid
        psiEq0(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,0,lsr,lsz,spdR,spdZ,psirz_Eq,0,0)
        dpsifrcdr(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,1,lsr,lsz,spdR,spdZ,psirz_Eq,0,0)
        dpsifrcdz(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,2,lsr,lsz,spdR,spdZ,psirz_Eq,0,0)
      enddo

!$omp parallel do private(i) 
      do i=1,nrzgrid
        meshne_cy(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,0,lsr,lsz,spdR,spdZ,lrne_Eq,0,0)
        meshni_cy(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,0,lsr,lsz,spdR,spdZ,lrni_Eq,0,0)
        meshte_cy(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,0,lsr,lsz,spdR,spdZ,lrte_Eq,0,0)
        meshti_cy(i) = spline2d(RZ_reg(i,1),RZ_reg(i,2)+zmax,0,lsr,lsz,spdR,spdZ,lrti_Eq,0,0)
      enddo

end subroutine LR_initial

