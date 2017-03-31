module system_env
  use mpi
  use omp_lib
#ifdef _OPENACC
  use openacc
  use cudafor
#endif
end module system_env

module precision
  use system_env
  integer, parameter :: doubleprec=kind(1.0d0),&
    singleprec=kind(1.0e0),defaultprec=kind(0.0)
	
	! lk stands for Length of real Kind
#ifdef DOUBLE_PRECISION
  integer, parameter :: lk=doubleprec,mpi_Rsize=MPI_DOUBLE_PRECISION,&
                        mpi_Csize=MPI_DOUBLE_COMPLEX
#else
  integer, parameter :: lk=singleprec,mpi_Rsize=MPI_REAL,&
                        mpi_Csize=MPI_COMPLEX
#endif
  real,parameter:: machineEpsilon=10.0_lk*epsilon(1.0_lk)
  save
end module precision

module global_parameters
  use precision
  integer,parameter :: gtcout=11, num_modes=8

! control parameters
  integer track_particles,ndata3d,mgrid,mpsi,mthetamax,mtoroidal,iodiag,iodata1d,&
       nspecies,istep,ndiag,msnap,mstep,mstepall,izonal,nbound,nboundR,irun,irk,idiag,&
       ncyclee,ncyclefe,mtdiag,nhybrid,ihybrid,nparami,magnetic,nonlinear,nfilter,numereq,&
       neop,neot,neoz,eqcurrent,fielddir,nparamf,nparamfe,nparame,mpsilow,mgridlow,&
       mpsihigh,mgridhigh,bcond,fieldmodel,ilaplacian,antenna,tearingmode,myth0,&
       n_modes(num_modes),m_modes(num_modes),island,fem,mgrid_fem,trilength,ismooth,inorm,irestore,&
       irotation,idiagonal,sdp_output
  real(lk) paranl,psi0,psi1,rg0,pi,pi2,pi2_inv,tstep,ulength,utime,rho0,maxwell(100001),&
       r0,b0,etemp0,eden0,hypr1,hypr2,omega_antenna,eta,toroidaln
#ifdef _FRC
  real(16) torbound
#else
  real(lk) torbound
#endif
#ifndef GPU_UM
  !$acc declare create(nparami,nparame,nparamf,nparamfe)
#endif
! MPI toroidal and particle decomposion
  integer mype,numberpe,npartdom,toroidal_comm,partd_comm,nproc_partd,myrank_partd,&
       nproc_toroidal,myrank_toroidal,left_pe,right_pe,&
       toroidal_domain_location,particle_domain_location

#ifdef _OPENMP
  integer nthreads
#endif

!!XY rolling restart
  integer irest,FileExit
  character(len=10) restart_dir1,restart_dir2
  save
end module global_parameters

module equilibrium
  use precision
  integer lsp,lst
#ifdef _TOROIDAL3D
  integer,parameter :: spdim=27 ! non-axisymmetric equilibrium
#else
  integer,parameter :: spdim=9
#endif
  real(lk) psiw,ped,spdpsi,spdtheta,spdrg,spdtor,spdpsi_inv,spdtheta_inv,spdrg_inv,spdtor_inv
  real(lk),dimension(:),allocatable :: stpp,mipp,mapp,mesher
  real(lk),dimension(:,:),allocatable :: qpsi,gpsi,ppsi,rpsi,torpsi,tpsi,npsi,nepp,tepp,&
       tipp,nipp,tfpp,nfpp,tfepp,nfepp,zepp,ropp,erpp,spcos,spsin,rgpsi,psirg,psitor,cpsi,xygrid
  real(lk),dimension(:,:,:),allocatable :: bsp,xsp,zsp,gsp,jsp,fsp,rd,nu,dl,ha,hb
#ifndef GPU_UM
  !$acc declare create(qpsi,gpsi,ropp,erpp,rgpsi,cpsi)
  !$acc declare create(bsp,xsp,mesher)
#endif
  save
end module equilibrium

module magnetic_island
  use precision   
  integer,parameter :: l=1 !island number  
  integer wi,wo,ires,qs
  integer,dimension(:),allocatable :: isn,ism
  real(lk),dimension(:,:),allocatable :: hmesh_total,hmesh_perturb,hangle,alphaIs !helical flux , alpha_zero, ksi in Eq(15), alpha in Eq(11)
  real(lk),dimension(:,:,:),allocatable :: gradalphaIs !gradient of alphais
#ifndef GPU_UM
  !$acc declare create(gradalphaIs,alphaIs)
#endif
  save
end module magnetic_island

module cylindricalRZ
  use precision
  integer  lsr,lsz,nbbbs,limitr 
  real(lk) rdim,zdim,rleft,psi_sep,spdR,spdZ
  real(lk),dimension(:),allocatable :: fpol,rbbbs,zbbbs,rlim,zlim
  real(lk),dimension(:,:),allocatable :: psirz
  real(lk),dimension(:,:,:),allocatable :: psirz_eq,psieq_p
  save
end module cylindricalRZ

module particle_array
  use precision
! particle diagnostics: # of quantities per species in history.out and data1d.out
  integer,parameter :: mpdiag=10,mpdata1d=3

! electron
  integer me,me1,memax,ngyroe,eload,etrap,ecoll
  real(lk) qelectron,aelectron,betae,tauee,tauei
  real(lk),dimension(mpdiag) :: diagelectron
  integer,dimension(:,:),allocatable :: jtelectron0,jtelectron1
  real(lk),dimension(:),allocatable :: wzelectron,meshte,meshne,kapane,kapate,&
    dtndpsi,pmarke,zonale,zonalce,markere,rdteme,pfluxe,markeret
  real(lk),dimension(:,:),allocatable :: zelectron,zelectron0,zelectron1,&
    wpelectron,wtelectron0,wtelectron1,densitye,flowe,phit,data1de,dnet,&
    pressureepara,pressureeperp
#ifndef GPU_UM
  !$acc declare create(diagelectron,jtelectron0,jtelectron1)
  !$acc declare create(wzelectron,meshte,meshne,kapane,kapate,rdteme)
  !$acc declare create(meshte,meshne,kapane,kapate)
  !$acc declare create(zelectron,zelectron0,zelectron1,wpelectron,wtelectron0)
  !$acc declare create(wtelectron1,densitye,flowe,data1de,phit,dnet,pressureepara,pressureeperp)
#endif
  real(lk),dimension(:,:,:),allocatable :: phisave,dnesave
  real(lk) tfracn
  real(lk),dimension(:),allocatable :: tfrac

! themal ion (main ion)
  integer mi,mimax,ngyroi,iload,icoll
  real(lk) qion,aion,tauii
  real(lk),dimension(mpdiag) :: diagion
  integer,dimension(:,:),allocatable :: jtion0,jtion1
  real(lk),dimension(:),allocatable :: wzion,meshti,meshni,kapani,kapati,&
       jacobianpsi,pmarki,zonali,zonalci,markeri,rdtemi,pfluxi,markerit
  real(lk),dimension(:,:),allocatable :: zion,zion0,wpion,wtion0,wtion1,&
    densityi,flowi,data1di
#ifndef GPU_UM
  !$acc declare create(diagion,jtion0,jtion1)
  !$acc declare create(wzion,meshti,meshni,kapani,kapati,rdtemi)
  !$acc declare create(zion,zion0,wpion,wtion0,wtion1,densityi,flowi,data1di)
#endif

! fast ion
  integer mf,mfmax,ngyrof,fload
  real(lk) qfast,afast,sd_v0,sd_vc,sd_l0,sd_widthInv
  real(lk),dimension(mpdiag) :: diagfast
  integer,dimension(:,:),allocatable :: jtfast0,jtfast1
  real(lk),dimension(:),allocatable :: wzfast,meshtf,meshnf,kapanf,kapatf,&
       pmarkf,zonalf,zonalcf,markerf,rdtemf,pfluxf,markerft
  real(lk),dimension(:,:),allocatable :: zfast,zfast0,wpfast,wtfast0,wtfast1,densityf,flowf,data1df
#ifndef GPU_UM
  !$acc declare create(diagfast,jtfast0,jtfast1)
  !$acc declare create(wzfast,meshtf,meshnf,kapanf,kapatf,rdtemf)
  !$acc declare create(zfast,zfast0,wpfast,wtfast0,wtfast1,densityf,flowf,data1df)
#endif

!fast electron
  integer mfe,mfe1,mfemax,ngyrofe,feload,fetrap
  real(lk) qfaste,afaste
  real(lk),dimension(mpdiag) :: diagfaste
  integer,dimension(:,:),allocatable :: jtfaste0,jtfaste1
  real(lk),dimension(:),allocatable :: wzfaste,meshtfe,meshnfe,kapanfe,kapatfe,&
       pmarkfe,zonalfe,zonalcfe,markerfe,rdtemfe,pfluxfe,markerfet
  real(lk),dimension(:,:),allocatable :: zfaste,zfaste0,zfaste1,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,data1dfe
  real(lk) fetfracn
  real(lk),dimension(:),allocatable :: fetfrac
#ifndef GPU_UM
  !$acc declare create(diagfaste,jtfaste0,jtfaste1)
  !$acc declare create(wzfaste,meshtfe,meshnfe,kapanfe,kapatfe,rdtemfe)
  !$acc declare create(zfaste,zfaste0,zfaste1,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,data1dfe)
#endif
  save
end module particle_array

module particle_tracking
  use precision
  real(lk),dimension(:,:),allocatable :: ptrackedi,ptrackede,ptrackedf
  integer,dimension(3) :: ntrackp
  save
end module particle_tracking

module field_array
  use precision
! PIC global fieldline following mesh
  real(lk) deltar,deltaz,zeta1,zeta0
  integer,dimension(:),allocatable :: itran,igrid,mtheta,nshift,igrid_fem
  real(lk),dimension(:),allocatable :: deltap,deltat,qtinv,qmesh,bmesh,jmesh,psimesh,tormesh,meshze
#ifndef GPU_UM
  !$acc declare create(igrid,mtheta)
  !$acc declare create(deltap,deltat,qtinv,psimesh)
#endif
  real(lk),dimension(:),allocatable :: gpsi200,b2m00,wtshift
  real(lk),dimension(:),allocatable :: gupp,gupt,gutt,guzz,gupz,gutz,rhom
  real(lk),dimension(:),allocatable :: gdpp,gdpt,gdtt,gdzz,gdpz,gdtz
  real(lk),dimension(:),allocatable :: spectraln

! fields on mesh: phi, apara, fluidne, fluidue, zonal and gradient
  real(lk),dimension(:),allocatable :: phi00,phip00,apara00,apara00nl,apara00nl0,fluidne00,d4fluidne,d2apara
  real(lk),dimension(:,:),allocatable :: phi,apara,fluidne,fluidue,apara0,fluidne0,&
                                         deltapsi,deltapsi0,sdeltapsi,&
                                         sapara,sfluidne,sdelapara,MHDxi_mag,phi_zero
  real(lk),dimension(:,:,:),allocatable :: gradphi,gradapara,gradue,gradne,gradext,gradpsi,gradphieff,&
                                           gradkne,gradpepara,gradpeperp,gradgaugef,MHDxiperp
#ifndef GPU_UM
  !$acc declare create(phip00,sapara)
  !$acc declare create(gradphi,gradapara,gradpsi,gradphieff)
#endif

! external field
  real(lk),dimension(:),allocatable :: omega
  real(lk),dimension(:,:),allocatable :: phiext,dn_ext

! diagnostics and filtering
  integer iflux,modes,solvermethod
  integer,dimension(:),allocatable :: nmodes,mmodes
  integer nmode

! radial interpolation
  integer,dimension(:,:),allocatable :: jtp1,jtp2
  real(lk),dimension(:,:),allocatable :: wtp1,wtp2

! laplacian coefficients
  integer mindexlap,mindexlap2,mindex_fem
  integer,dimension(:),allocatable :: nindexlap,nindexlap2,nindex_fem
  integer,dimension(:,:),allocatable :: indexplap,indexlap2,indexp_fem
  integer,dimension(:,:),allocatable :: trilist
  real(lk),dimension(:,:),allocatable :: lapmat,lapmat2,lmat,dmat

! theta of max and min b-field for particle baoundary cond.
  real(lk) maxbfield(0:1),minbfield(0:1),thetabmin(0:1),thetabmax(0:1)
  integer,dimension(:,:),allocatable :: thetaupp,thetadwn

! gyro averaging
  real(lk),dimension(:,:),allocatable :: pgyro,tgyro,pgyro2,tgyro2
#ifndef GPU_UM
  !$acc declare create(pgyro,tgyro,pgyro2,tgyro2)
#endif
  save
end module field_array

module petsc_array
  use precision
  integer newcomm,nproc_newcomm,myrank_newcomm
  integer,dimension(:),allocatable :: userp,users,luserp,lusers,luserp2,lusers2
  real(lk),dimension(:),allocatable :: usera,userb,userx,lusera,luserb,luserx,lusera2,luserb2,luserx2
  save
end module petsc_array

module data_type
  integer,parameter :: bp_char=0,bp_short=1,&
                bp_int=2,bp_long=3,bp_longlong=4,bp_float=5,&
                bp_double=6,bp_longdouble=7,bp_pointer=8,&
                bp_string=9,bp_uchar=50,bp_ushort=51,bp_uint=52,&
                bp_ulong=53,bp_ulonglong=54
  save
end module data_type

module LAMYRIDGEeq
  use precision
  integer numr,numz,nrzgrid
  integer,dimension(:),allocatable :: igrid_frc,nz_frc

  real(lk) zmax,psi0_frc,psi1_frc,Zsimulation0,Zsimulation1,min_r,&
           lz,delta_Z,deltapsi_frc           
  real(lk),dimension(:),allocatable :: r_Eq,z_Eq,psiEq0,meshni_cy,meshne_cy,&
                                       meshti_cy,meshte_cy,psi_frc,psifrc_tmp,&
                                       deltar_frc,dpsifrcdr,dpsifrcdz
  real(lk),dimension(:,:),allocatable ::RZ_reg
  real(lk),dimension(:,:,:),allocatable :: B0mag_Eq,B0R_Eq,B0Z_Eq,lrne_Eq,lrni_Eq,lrti_Eq,lrte_Eq
  save
end module LAMYRIDGEeq

module interfaces
  interface
    subroutine push(species_name,icycle,irksub,ihybrid)
      implicit none
      character(*),intent(in) :: species_name
      integer,intent(in),optional :: icycle,irksub,ihybrid
    end subroutine push

    subroutine axisExtrapolate(farray)
      use precision
      use global_parameters,only: mgrid

      real(lk),dimension(mgrid) :: farray
    end subroutine axisExtrapolate

  end interface
end module interfaces
