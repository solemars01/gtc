subroutine shift(species_name)
  use global_parameters,only: nparami,nparame,nparamf,nparamfe,&
    mpi_comm_world,mtoroidal,mype,left_pe,right_pe,myrank_toroidal,&
    toroidal_comm,numberpe,pi
  use particle_array,only: mimax,mi,zion,zion0,memax,me,zelectron,zelectron0,&
    mfmax,mf,zfast,zfast0,mfemax,mfe,zfaste,zfaste0
  use field_array,only: zeta0,zeta1

  implicit none
  integer ierror
  character(*),intent(in) :: species_name

#ifdef _OPENACC
  interface
    subroutine shift_cuda(zpart,zpart0,mpmax,mp,nparam,mtoroidal,mype,&
        left_pe,right_pe,myrank_toroidal,toroidal_comm0,numberpe,pi,zeta0,&
        zeta1) bind(C, name='shift_cuda')
      use iso_c_binding
      real(4),device :: zpart(nparam,*),zpart0(nparam,*)
      integer :: mp
      integer,value :: mpmax,nparam,mtoroidal,mype,left_pe,right_pe,&
        myrank_toroidal,toroidal_comm0,numberpe
      real(4),value :: pi,zeta0,zeta1
    end subroutine shift_cuda
  end interface
#define uniShift shift_cuda
#else
  interface
    subroutine shiftParticle(zpart,zpart0,mpmax,mp,nparam,mtoroidal,&
        mype,left_pe,right_pe,myrank_toroidal,toroidal_comm0,numberpe,pi,&
        zeta0,zeta1)
      use precision
      implicit none
      
      integer mpmax,mp,nparam,mtoroidal,mype,left_pe,right_pe,myrank_toroidal,&
        toroidal_comm0,numberpe
      real(lk),dimension(:,:) :: zpart,zpart0
      real(lk) pi,zeta0,zeta1
    end subroutine shiftParticle
  end interface
#define uniShift shiftParticle
#endif

  if(species_name=='thermal-ion')then
    !$acc host_data use_device(zion,zion0)
    call uniShift(zion,zion0,mimax,mi,nparami,mtoroidal,mype,&
      left_pe,right_pe,myrank_toroidal,toroidal_comm,numberpe,pi,zeta0,zeta1)
    !$acc end host_data
  elseif(species_name=='thermal-electron')then
    !$acc host_data use_device(zelectron,zelectron0)
    call uniShift(zelectron,zelectron0,memax,me,nparame,mtoroidal,mype,&
      left_pe,right_pe,myrank_toroidal,toroidal_comm,numberpe,pi,zeta0,zeta1)
    !$acc end host_data
  elseif(species_name=='fast-ion')then
    !$acc host_data use_device(zfast,zfast0)
    call uniShift(zfast,zfast0,mfmax,mf,nparamf,mtoroidal,mype,&
      left_pe,right_pe,myrank_toroidal,toroidal_comm,numberpe,pi,zeta0,zeta1)
    !$acc end host_data
  elseif(species_name=='fast-electron')then
    !$acc host_data use_device(zfaste,zfaste0)
    call uniShift(zfaste,zfaste0,mfemax,mfe,nparamfe,mtoroidal,mype,&
      left_pe,right_pe,myrank_toroidal,toroidal_comm,numberpe,pi,zeta0,zeta1)
    !$acc end host_data
  else
    write(*,*)'shift.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
#ifdef _OPENACC
  call acc_async_wait_all()
#endif
end subroutine shift

subroutine shiftParticle(zpart,zpart0,mpmax,mp,nparam,mtoroidal,mype,&
    left_pe,right_pe,myrank_toroidal,toroidal_comm0,numberpe,pi,zeta0,zeta1)
  use precision
  use global_parameters,only: torbound
  implicit none
  
  integer i,m,msendleft(2),msendright(2),mrecvleft(2),mrecvright(2),&
    mtop,m0,msend,msendtmp,mrecvtmp,idest,isource,isendtag,irecvtag,nzphase,&
    kzpart(mpmax),iright(mpmax),ileft(mpmax),isendcount,irecvcount,&
    istatus(MPI_STATUS_SIZE),ierror,iteration,lasth,ompswitch
  real(lk) zetaright,zetaleft,pi2_inv
#ifdef _FRC
  real(16) torbound_inv
#else
  real(lk) torbound_inv
#endif
  real(lk),dimension(:,:),allocatable :: recvleft,recvright,sendleft,sendright

  integer mpmax,mp,nparam,mtoroidal,mype,left_pe,right_pe,myrank_toroidal,&
    toroidal_comm0,numberpe
  real(lk) pi,zeta0,zeta1
  real(lk),dimension(:,:) :: zpart,zpart0
#ifdef _OPENMP
  integer msleft(32,0:15),msright(32,0:15)
  integer nthreads,gnthreads,iam,delm,mbeg,mend
#endif

  !For developing/debug mode using single mpi, and you want to skip shift
  if(mtoroidal==1)return

  torbound_inv=1.0_lk/torbound
  pi2_inv=0.5_lk/pi
  nzphase=2*nparam ! nzphase=14 if track_particles=1, =12 otherwise
  m0=1
  iteration=0
  
100 iteration=iteration+1
  !For debuging shift using single mpi
  if(iteration>max(mtoroidal,2))then
     write(*,*)'endless particle sorting loop at PE=',mype
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

  msend=0
  msendright=0
  msendleft=0

  if(m0 <= mp)then
!$omp parallel do private(m)
     do m=1,mp
        kzpart(m)=0
     enddo

#ifdef _OPENMP
! This section of code (down to #else) is included by the preprocessor when
! the compilation is performed with OpenMP support. We must then use a few
! temporary arrays and add some work distribution code.
! First, we initialize the shared (and temporary) arrays msleft and msright
! to zero.
     msleft=0
     msright=0

! Then we start the parallel region with
!$omp parallel private(nthreads,iam,delm,i,mbeg,mend,m,zetaright,zetaleft) &
!$omp& shared(gnthreads,msleft,msright)
     nthreads=omp_get_num_threads()    !Get the number of threads ready to work
     iam=omp_get_thread_num()          !Get my thread number (position)
     delm=(mp-m0+1)/nthreads       !Calculate the number of steps per thread
     i=mod((mp-m0+1),nthreads)
!$omp single                 !Put nthread in global variable for later use.
     gnthreads=nthreads      !nthread is the same for all threads so only one
!$omp end single nowait      !of them needs to copy the value in gnthreads

! We now distribute the work between the threads. The loop over the particles
! is distributed equally (as much as possible) between them.
     mbeg=m0+min(iam,i)*(delm+1)+max(0,(iam-i))*delm
     mend=mbeg+delm+(min((iam+1),i)/(iam+1))-1

! label particle to be moved
     do m=mbeg,mend
        zetaright=min(torbound,zpart(3,m))-zeta1
        zetaleft=zpart(3,m)-zeta0
        
        if(zetaright*zetaleft > 0.0_lk)then
           zetaright=zetaright*torbound_inv
           zetaright=zetaright-real(floor(zetaright))
           msright(3,iam)=msright(3,iam)+1
           kzpart(mbeg+msright(3,iam)-1)=m
           
           if(zetaright < 0.5_lk)then
! particle to move right
              msright(1,iam)=msright(1,iam)+1
              iright(mbeg+msright(1,iam)-1)=m

! particle to move left
           else
              msleft(1,iam)=msleft(1,iam)+1
              ileft(mbeg+msleft(1,iam)-1)=m
           endif
        endif
     enddo
! End of the OpenMP parallel region
!$omp end parallel

! Now that we are out of the parallel region we need to gather and rearrange
! the results of the multi-thread calculation. We need to end up with the
! same arrays as for the sequential (single-threaded) calculation.
     do m=0,gnthreads-1
        delm=(mp-m0+1)/gnthreads
        i=mod((mp-m0+1),gnthreads)
        mbeg=m0+min(m,i)*(delm+1)+max(0,(m-i))*delm
        if( msleft(2,m) /= 0 )msendleft(2)=msendleft(1)+msleft(2,m)
        do i=1,msleft(1,m)
           ileft(msendleft(1)+i)=ileft(mbeg+i-1)
        enddo
        msendleft(1)=msendleft(1)+msleft(1,m)
        if( msright(2,m) /= 0 )msendright(2)=msendright(1)+msright(2,m)
        do i=1,msright(1,m)
           iright(msendright(1)+i)=iright(mbeg+i-1)
        enddo
        msendright(1)=msendright(1)+msright(1,m)
        do i=1,msright(3,m)
           kzpart(msend+i)=kzpart(mbeg+i-1)
        enddo
        msend=msend+msright(3,m)
     enddo

#else
!  This section of code replaces the section above when the compilation does
!  NOT include the OpenMP support option. Temporary arrays msleft and msright
!  are not needed as well as the extra code for thread work distribution.

     do m=m0,mp
        zetaright=min(torbound,zpart(3,m))-zeta1
        zetaleft=zpart(3,m)-zeta0
        
        if( zetaright*zetaleft > 0.0_lk )then
           zetaright=zetaright*torbound_inv
           zetaright=zetaright-real(floor(zetaright))
           msend=msend+1
           kzpart(msend)=m
           
           if( zetaright < 0.5_lk )then
! # of particle to move right
              msendright(1)=msendright(1)+1
              iright(msendright(1))=m

! # of particle to move left
           else
              msendleft(1)=msendleft(1)+1
              ileft(msendleft(1))=m
           endif
        endif
     enddo
     
#endif

  endif

  if (msend /= (msendleft(1)+msendright(1))) then
     write(*,*)'mype=',mype,'  msend NOT equal to msendleft+msendright'
     msend=msendleft(1)+msendright(1)
  endif

  if(iteration>1)then

! test: at least 1 particle needs to be shifted.
     mrecvtmp=0
     msendtmp=0
     if(msend>0)msendtmp=1

     call MPI_ALLREDUCE(msendtmp,mrecvtmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierror)

! no particle to be shifted, return
     if ( mrecvtmp == 0 ) then
!        write(0,*)istep,irk,mype,mp,m0,iteration
        return
     endif
  endif

! an extra space to prevent zero size when msendright(1)=msendleft(1)=0
  allocate(sendright(1:nzphase,max(1,msendright(1))),sendleft(1:nzphase,max(1,msendleft(1))))

! pack particle to move right
!$omp parallel do private(m)
  do m=1,msendright(1)
     sendright(1:nparam,m)=zpart(1:nparam,iright(m))
     sendright(nparam+1:2*nparam,m)=zpart0(1:nparam,iright(m))
  enddo

! pack particle to move left
!$omp parallel do private(m)
  do m=1,msendleft(1)    
     sendleft(1:nparam,m)=zpart(1:nparam,ileft(m))
     sendleft(nparam+1:2*nparam,m)=zpart0(1:nparam,ileft(m))
  enddo

!"fill the holes" and "send/recv" are two independent jobs. So if we have more than 1 OpenMP threads we can use 2 threads to do these 2 jobs at the same time.
!In the cases when we use only OpenMP thread or we do not use OpenMP,
!ompswitch=1, so these 2 part of code will be executed in series.
!In the cases when we use more than 2 OpenMP threads. ompswitch=1 for the second
!thread, ompswitch=0 for all other threads. So only the second thread will
!repack the particles left and the first thread (master thread) will send and
!receive particles from other MPI tasks at the same time.

!$OMP PARALLEL PRIVATE(iam,ompswitch,nthreads)
  ompswitch=1
#ifdef _OPENMP
  iam=omp_get_thread_num()
  nthreads=omp_get_num_threads()
  if (iam.ne.1) ompswitch=0
  if (nthreads.eq.1) ompswitch=1
#endif

! independent job 1: fill the holes. 
! This part will be executed only when ompswitch.eq.1 
  if (ompswitch.eq.1) then
  mtop=mp
! # of particles remain on local PE
  mp=mp-msendleft(1)-msendright(1)
  lasth=msend
  do i=1,msend
     m=kzpart(i)
     if (m > mp) exit  !Break out of the DO loop if m > mp
     do while(mtop == kzpart(lasth))
        mtop=mtop-1
        lasth=lasth-1
     enddo
     zpart(1:nparam,m)=zpart(1:nparam,mtop)
     zpart0(1:nparam,m)=zpart0(1:nparam,mtop)
     mtop=mtop-1
     if (mtop == mp) exit  !Break out of the DO loop if mtop=mp
  enddo
  endif


! independent job 2: send/recv,  only the master thread will do this part.
!$OMP MASTER

! send # of particle to move right
  mrecvleft=0
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendright,2,MPI_INTEGER,idest,isendtag,&
       mrecvleft,2,MPI_INTEGER,isource,irecvtag,toroidal_comm0,istatus,ierror)
  allocate(recvleft(1:nzphase,max(1,mrecvleft(1))))
 
! send particle to right and receive from left
  recvleft=0.0
  isendcount=max(1,msendright(1))*nzphase
  irecvcount=max(1,mrecvleft(1))*nzphase
  call MPI_SENDRECV(sendright,isendcount,mpi_Rsize,idest,isendtag,recvleft,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm0,istatus,ierror)
  
! send # of particle to move left
  mrecvright=0
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(msendleft,2,MPI_INTEGER,idest,isendtag,&
       mrecvright,2,MPI_INTEGER,isource,irecvtag,toroidal_comm0,istatus,ierror)
  
  allocate(recvright(1:nzphase,max(1,mrecvright(1))))

! send particle to left and receive from right
  recvright=0.0
  isendcount=max(1,msendleft(1))*nzphase
  irecvcount=max(1,mrecvright(1))*nzphase
  call MPI_SENDRECV(sendleft,isendcount,mpi_Rsize,idest,isendtag,recvright,&
       irecvcount,mpi_Rsize,isource,irecvtag,toroidal_comm0,istatus,ierror)
!$OMP END MASTER
! end of 2 independent jobs
!$OMP END PARALLEL
  
! need extra particle array
  if(mp+mrecvleft(1)+mrecvright(1) > mpmax)then
     write(*,*)"need bigger particle array",mype,mpmax,mp+mrecvleft(1)+mrecvright(1)
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif

! unpack particle, particle moved from left
!$omp parallel do private(m)
  do m=1,mrecvleft(1)
     zpart(1:nparam,m+mp)=recvleft(1:nparam,m)
     zpart0(1:nparam,m+mp)=recvleft(nparam+1:2*nparam,m)
  enddo

! particle moved from right
!$omp parallel do private(m)
  do m=1,mrecvright(1)
     zpart(1:nparam,m+mp+mrecvleft(1))=recvright(1:nparam,m)
     zpart0(1:nparam,m+mp+mrecvleft(1))=recvright(nparam+1:2*nparam,m)
  enddo
  
  mp=mp+mrecvleft(1)+mrecvright(1)
  
  deallocate(sendleft,sendright,recvleft,recvright)
  m0=mp-mrecvright(1)-mrecvleft(1)+1
  goto 100
  
end subroutine shiftParticle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine locate(species_name)
  use particle_array
  implicit none

  integer ierror
  character(*),intent(in) :: species_name

  interface
     subroutine locateParticle(zpart,wzpart,wppart,wtpart0,wtpart1,&
        jtpart0,jtpart1,qpart,apart,ngyro,mp)
      use precision
      implicit none

      real(lk) qpart,apart
      real(lk),dimension(:) :: wzpart
      real(lk),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(lk),dimension(:,:) :: zpart
      integer,dimension(:,:) :: jtpart0,jtpart1
      integer ngyro,mp
    end subroutine locateParticle
  end interface
  
  if(species_name=='thermal-ion')then
    call locateParticle(zion,wzion,wpion,wtion0,wtion1,jtion0,jtion1,qion,&
      aion,ngyroi,mi)
  elseif(species_name=='thermal-electron')then
    call locateParticle(zelectron,wzelectron,wpelectron,wtelectron0,&
      wtelectron1,jtelectron0,jtelectron1,qelectron,aelectron,ngyroe,me)
  elseif(species_name=='fast-ion')then
    call locateParticle(zfast,wzfast,wpfast,wtfast0,wtfast1,jtfast0,jtfast1,&
      qfast,afast,ngyrof,mf)
  elseif(species_name=='fast-electron')then
    call locateParticle(zfaste,wzfaste,wpfaste,wtfaste0,wtfaste1,jtfaste0,&
      jtfaste1,qfaste,afaste,ngyrofe,mfe)
  else
    write(*,*)'push.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine locate

! find ion location for interpolation in both gathering (chargei.F90) and scattering (pushi.F90)
subroutine locateParticle(zpart,wzpart,wppart,wtpart0,wtpart1,jtpart0,&
    jtpart1,qpart,apart,ngyro,mp)
  use precision
  use global_parameters,only: rho0,mpsilow,rg0,psi0,psi1,pi2,fielddir,mpsi,&
    mype,istep
  use field_array,only: deltar,deltat,deltap,deltaz,psimesh,zeta0,mtheta,&
    igrid,pgyro,pgyro2,tgyro,tgyro2,qmesh,qtinv,solvermethod
  use equilibrium,only: lsp,spdpsi_inv,spdpsi,rgpsi
  implicit none

  !declaration of the dummy arguments
  real(lk) qpart,apart
  real(lk),dimension(:) :: wzpart
  real(lk),dimension(:,:) :: wppart,wtpart0,wtpart1
  real(lk),dimension(:,:) :: zpart
  integer,dimension(:,:) :: jtpart0,jtpart1
  integer ngyro,mp

  integer m,igyro,ip,jt,ij0,ii,im,j00,j01,isp
  real(lk) pdum0,tdum0,zdum,rho,rg,tflr,pdum,tdum,dpx,delr,delt(0:mpsi),delp(mpsi),delz,rhoc,xdum,ydum,sqrtq0_inv,paxis
#ifndef GPU_UM
  !$acc declare create(delt(0:mpsi),delp(mpsi))
#endif

  delr=1.0_lk/deltar
  delt=1.0_lk/deltat
  delz=1.0_lk/deltaz
  delp=1.0_lk/deltap
#ifndef GPU_UM
  !$acc update device(delt(0:mpsi),delp(mpsi))
#endif
  rhoc=sqrt(apart)/(rho0*qpart) !apart & qpart are needed to calculate ion gyroradius
  if(ngyro==1)rhoc=0.0_lk

  sqrtq0_inv=1.0_lk/sqrt(qmesh(0))
  paxis=0.0_lk
  if(psi0==0.0_lk)paxis=12.5_lk*rho0*rho0*sqrtq0_inv!Use x-y locate for r < 5 rho0


#ifdef _OPENACC
  !$acc parallel loop gang vector
#else
  !$omp parallel do private(m,igyro,pdum0,tdum0,zdum,rho,isp,dpx,&
  !$omp& rg,ip,jt,ij0,pdum,ii,tflr,im,tdum,j00,j01)
#endif
  do m=1,mp
     pdum0=zpart(1,m)
     tdum0=zpart(2,m)
     zdum=zpart(3,m)
     rho=zpart(6,m)*rhoc
     wzpart(m)=(zdum-zeta0)*delz  !weight for upper toroidal grid

! Guiding Center location
     isp=max(1,min(lsp-1,ceiling(pdum0*spdpsi_inv)))
     dpx=pdum0-spdpsi*real(isp-1)
! radaial spline of rg avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dpx*dpx
     ip=max(0,min(mpsi,floor((rg-rg0)*delr+0.5_lk)))       !radial grid closest to GC
     jt=max(0,min(mtheta(ip),floor(tdum0*delt(ip)+0.5_lk))) !poloidal grid closest to GC
     ij0=igrid(ip)+jt                                   !GC grid index in magnetic coordinates

     do igyro=1,ngyro
!        pdum=max(psimesh(0),min(psimesh(mpsi),pdum0+rho*pgyro(igyro,ij0))) !particle position
        pdum=pdum0+rho*pgyro(igyro,ij0)+rho*rho*pgyro2(igyro,ij0) !particle position
        !if(pdum<0)print *,"pdum,bounded pdum:",pdum,modulo(pdum-psi0,psi1-psi0)+psi0
        !pdum=modulo(pdum-psi0,psi1-psi0)+psi0 !periodic BC
        if(pdum < psi0)pdum=2.0_lk*psi0-pdum !reflective BC
        if(pdum > psi1)pdum=2.0_lk*psi1-pdum
! particle position in theta
        tflr=tdum0+rho*tgyro(igyro,ij0)+rho*rho*tgyro2(igyro,ij0)
! Near axis gyroaveraging
        if(pdum0<paxis)then
          pdum = sqrt(pdum0)
          xdum = pdum*cos(tdum0)+rho*rho0*sqrtq0_inv*cos(tdum0+pi2*real(igyro-1)/4.0_lk) 
          ydum = pdum*sin(tdum0)+rho*rho0*sqrtq0_inv*sin(tdum0+pi2*real(igyro-1)/4.0_lk) 
          pdum = max(1.0e-8_lk*psi1,xdum*xdum+ydum*ydum)
          tflr = sign(1.0_lk,ydum)*acos(max(-1.0_lk,min(1.0_lk,xdum/&
            sqrt(zpart(1,m)))))
        endif
        isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
        dpx=pdum-spdpsi*real(isp-1)
! radaial spline of rg avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
        if(isp==1)dpx=sqrt(dpx)
        rg=rgpsi(1,isp)+rgpsi(2,isp)*dpx+rgpsi(3,isp)*dpx*dpx
        ii=max(0,min(mpsi-1,floor((rg-rg0)*delr)))    !radial grid on inner flux surface
        !solvermethod for fluxtube simulations:
        !0:fluxtube
        !1:semispectral
#ifdef _FRC
        if(solvermethod==0)then
        ! assuming mpsi=2, for fluxtube, then flux surfaces = [0,1,2]
           if(ii<1)then
             wppart(igyro,m)=1.0_lk !all weight on outer flux surface between [0,1]
           else
             wppart(igyro,m)=0.0_lk !all weight on inner flux surface between [1,2]
           endif
        else
           wppart(igyro,m)=(pdum-psimesh(ii))*delp(ii+1) !weight for outer flux surface
        endif
#else
        wppart(igyro,m)=(pdum-psimesh(ii))*delp(ii+1) !weight for outer flux surface
#endif
! inner flux surface
        im=ii
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(tflr-(zdum-pi2)*qtinv(im),pi2)*delt(im)
        else
          tdum=modulo(tflr-zdum*qtinv(im),pi2)*delt(im)
        endif
        j00=max(0,min(mtheta(im)-1,floor(tdum))) 
        jtpart0(igyro,m)=igrid(im)+j00  !lower poloidal grid on inner flux surface
        wtpart0(igyro,m)=tdum-real(j00) !weight for upper poloidal grid on inner flux surface
        if(ii == 0 .and. psi0== 0.0_lk)wtpart0(igyro,m)=1.0_lk
! outer flux surface
        im=ii+1
        if(fielddir==1 .or. fielddir==3)then
          tdum=modulo(tflr-(zdum-pi2)*qtinv(im),pi2)*delt(im)
        else
          tdum=modulo(tflr-zdum*qtinv(im),pi2)*delt(im)
        endif
        j01=max(0,min(mtheta(im)-1,floor(tdum)))
        jtpart1(igyro,m)=igrid(im)+j01  !lower poloidal grid on outer flux surface
        wtpart1(igyro,m)=tdum-real(j01) !weight for upper poloidal grid on outer flux surface
     enddo
  enddo
  !$acc end parallel
end subroutine locateParticle
