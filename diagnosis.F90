subroutine diagnosis
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer,parameter :: nfield=3,mfdiag=4,mfdata1d=2
  integer i,j,k,ii,ip,icount,ierror,mtgrid,nf,mtest
  real(lk) fieldmode(2,modes,nfield),fieldtime(mfdiag,nfield),partdata(mpdiag,nspecies),&
       field00(0:mpsi,nfield),fieldrms(0:mpsi,nfield),rmstmp(0:mpsi,nfield)
  real(lk),dimension(:,:,:),allocatable :: fieldgrid
  real(lk),dimension(:),allocatable :: kapatmti,kapatmte,kapatmtf,kapatmtfe,kapatmni,kapatmne,kapatmnf,kapatmnfe
  save kapatmti,kapatmte,kapatmtf,kapatmtfe,kapatmni,kapatmne,kapatmnf,kapatmnfe

! open output files history.out and data1d.out
  if(mype==0 .and. istep==ndiag)call opendiag(mpdiag,nfield,mfdiag,modes,mpdata1d,mfdata1d)
  call MPI_BCAST(mstepall,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
!!! particle diagnosis: main ion, electron and/or EP, imputities, ...
! record the initial temperature/density gradient for normalization
  if(istep==ndiag)then
     allocate(kapatmti(0:mpsi),kapatmte(0:mpsi),kapatmtf(0:mpsi),kapatmtfe(0:mpsi),kapatmni(0:mpsi),&
     kapatmne(0:mpsi),kapatmnf(0:mpsi),kapatmnfe(0:mpsi),STAT=mtest)
     if (mtest /= 0) then
       write(0,*)mype,'*** Cannot allocate kapatmti: mtest=',mtest
       call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
     endif

     do i=0,mpsi
        kapatmti(i)=max(abs(kapati(i)),1.0e-2_lk)
        kapatmte(i)=max(abs(kapate(i)),1.0e-2_lk)
        kapatmtf(i)=max(abs(kapatf(i)),1.0e-2_lk)
        kapatmtfe(i)=max(abs(kapatfe(i)),1.0e-2_lk)
        kapatmni(i)=max(abs(kapani(i)),1.0e-2_lk)
        kapatmne(i)=max(abs(kapane(i)),1.0e-2_lk)
        kapatmnf(i)=max(abs(kapanf(i)),1.0e-2_lk)
        kapatmnfe(i)=max(abs(kapanfe(i)),1.0e-2_lk)
     enddo

  endif

!!convert heat and particle flux to Bohm unit with normalization on middle flux surface
  diagion(7)=diagion(7)/(rho0*rho0*kapatmni(iflux)*gpsi200(iflux))
  diagion(9)=diagion(9)/(rho0*rho0*kapatmti(iflux)*gpsi200(iflux))
  if(nhybrid>0)then
     diagelectron(7)=diagelectron(7)*tfracn/(rho0*rho0*kapatmne(iflux)*gpsi200(iflux))
     diagelectron(9)=diagelectron(9)*tfracn/(rho0*rho0*kapatmte(iflux)*gpsi200(iflux))
  endif
  if(fload>0)then
     diagfast(7)=diagfast(7)/(rho0*rho0*kapatmnf(iflux)*gpsi200(iflux))
     diagfast(9)=diagfast(9)/(rho0*rho0*kapatmtf(iflux)*gpsi200(iflux))
  endif
 
  if(feload>0)then
     diagfaste(7)=diagfaste(7)*fetfracn/(rho0*rho0*kapatmnfe(iflux)*gpsi200(iflux))
     diagfaste(9)=diagfaste(9)*fetfracn/(rho0*rho0*kapatmtfe(iflux)*gpsi200(iflux))
  endif

  partdata=0.0_lk
  do ip=1,nspecies
     if(ip==1)then
        partdata(:,ip)=diagion/diagion(10) !first species is thermal ion

     elseif(ip==2)then
        if(nhybrid>0)then
           partdata(:,ip)=diagelectron/diagelectron(10) !second species is electron
        elseif(feload>0)then
           partdata(:,ip)=diagfaste/diagfaste(10) !second species is fast electron
        else
           partdata(:,ip)=diagfast/diagfast(10) !second species is fast ion if electron and fast electron not loaded
        endif

     elseif(ip==3)then
        if(feload>0)then
          partdata(:,ip)=diagfaste/diagfaste(10) !third species is fast electron
        else
          partdata(:,ip)=diagfast/diagfast(10) !third species is fast ion
        endif
     endif
  enddo

!!! field diagnosis: phi, a_para, fluid_ne, ...
  ii=igrid(iflux)
  mtgrid=mtheta(iflux)
  fieldtime=0.0_lk 
  fieldmode=0.0_lk
  allocate(fieldgrid(0:1,0:mtgrid,nfield))
  
  do nf=1,nfield
     if(nf==1)then
        fieldgrid(:,:,nf)=phi(:,ii:ii+mtgrid)/(rho0*rho0)
!time history of field quantity at theta=zeta=0 & i=iflux
        fieldtime(1,nf)=phi(0,ii)/(rho0*rho0) 
        fieldtime(2,nf)=phip00(iflux)/rho0 
        fieldtime(3,nf)=sqrt(sum(phip00(0:mpsi)**2)/real(mpsi+1))/rho0
!$omp parallel do private(i)
        do i=0,mpsi
           field00(i,nf)=phi00(i)/(rho0*rho0)
           fieldrms(i,nf)=sum(phi(0,igrid(i):igrid(i)+mtheta(i)-1)**2)/(rho0**4)
        enddo

     elseif(nf==2)then
! apara is renormalized in such a way that it has the same
! amplitude as phi in ideal shear Alfven waves
        fieldgrid(:,:,nf)=sapara(:,ii:ii+mtgrid)/(rho0*sqrt(betae*aion))
        fieldtime(1,nf)=sapara(0,ii)/(rho0*sqrt(betae*aion))
        fieldtime(2,nf)=apara00(iflux)/(rho0*sqrt(betae*aion))
        fieldtime(3,nf)=sqrt(sum(apara00(0:mpsi)**2)/real(mpsi+1))/(rho0*sqrt(betae*aion)) 
!$omp parallel do private(i)
        do i=0,mpsi
           field00(i,nf)=apara00(i)/(rho0*sqrt(betae*aion))
           fieldrms(i,nf)=sum(sapara(0,igrid(i):igrid(i)+mtheta(i)-1)**2)/(rho0*rho0*betae*aion)
        enddo

     elseif(nf==3)then
        fieldgrid(:,:,nf)=sfluidne(:,ii:ii+mtgrid)
        fieldtime(1,nf)=sfluidne(0,ii)
        fieldtime(2,nf)=fluidne00(iflux)
        fieldtime(3,nf)=sqrt(sum(fluidne00(0:mpsi)**2)/real(mpsi+1))
!$omp parallel do private(i)
        do i=0,mpsi
           field00(i,nf)=fluidne00(i)
           fieldrms(i,nf)=sum(sfluidne(0,igrid(i):igrid(i)+mtheta(i)-1)**2)
        enddo
     endif
  enddo

!Volume-averaged RMS
  icount=nfield*(mpsi+1)
  call MPI_REDUCE(fieldrms,rmstmp,icount,mpi_Rsize,MPI_SUM,0,toroidal_COMM,ierror)
  fieldrms=rmstmp
  do nf=1,nfield
     fieldtime(4,nf)=sqrt(sum(fieldrms(:,nf))/real(mtoroidal*sum(mtheta)))
!$omp parallel do private(i)
     do i=0,mpsi
        fieldrms(i,nf)=sqrt(fieldrms(i,nf)/real(mtoroidal*mtheta(i)))
     enddo
  enddo

  call spectrum(nfield,modes,nmodes,mmodes,iflux,mtgrid,fieldgrid,fieldmode)

  deallocate(fieldgrid)

  if(mype==0)then
! write particle and field data to history file
     do i=1,nspecies
        do j=1,mpdiag
           write(iodiag,102)partdata(j,i)
        enddo
     enddo
     do i=1,nfield
        do j=1,mfdiag
           write(iodiag,102)fieldtime(j,i)
        enddo
     enddo
     do i=1,nfield
        do j=1,modes
           write(iodiag,102)fieldmode(1,j,i),fieldmode(2,j,i)
        enddo
     enddo
     call FLUSH(iodiag)

! convert 1d heat and particle diffusivity to Bohm unit
     do i=1,mpdata1d
         if(i==1)data1di(:,i)=data1di(:,i)/(rho0*rho0*kapatmni*gpsi200)
         if(i==2)data1di(:,i)=data1di(:,i)/(rho0*rho0*kapatmti*gpsi200)
         if(nspecies>1)then
            if(nhybrid>0)then
               if(i==1)data1de(:,i)=data1de(:,i)*tfrac/(rho0*rho0*kapatmne*gpsi200)
               if(i==2)data1de(:,i)=data1de(:,i)*tfrac/(rho0*rho0*kapatmte*gpsi200)
               if(i==3)data1de(:,i)=data1de(:,i)/sqrt(gpsi200)
               if(i==4)data1de(:,i)=data1de(:,i)/sqrt(gpsi200)
            endif

            if(fload>0)then
               if(i==1)data1df(:,i)=data1df(:,i)/(rho0*rho0*kapatmnf*gpsi200)
               if(i==2)data1df(:,i)=data1df(:,i)/(rho0*rho0*kapatmtf*gpsi200)
               if(i==3)data1df(:,i)=data1df(:,i)/sqrt(gpsi200)
               if(i==4)data1df(:,i)=data1df(:,i)/sqrt(gpsi200)
            endif
             
              if(feload>0)then
               if(i==1)data1dfe(:,i)=data1dfe(:,i)*fetfrac/(rho0*rho0*kapatmnfe*gpsi200)
               if(i==2)data1dfe(:,i)=data1dfe(:,i)*fetfrac/(rho0*rho0*kapatmtfe*gpsi200)
               if(i==3)data1dfe(:,i)=data1dfe(:,i)/sqrt(gpsi200)
               if(i==4)data1dfe(:,i)=data1dfe(:,i)/sqrt(gpsi200)
            endif

         endif
     enddo

! write radial-time data to data1d.out
     write(iodata1d,102)data1di
     if(nspecies>1)then
        if(nhybrid>0)write(iodata1d,102)data1de
        if(fload>0)write(iodata1d,102)data1df
        if(feload>0)write(iodata1d,102)data1dfe
     endif

     write(iodata1d,102)field00
     write(iodata1d,102)fieldrms
     call FLUSH(iodata1d)

! monitor: rms of phi, apara, fluidne, heatflux of ion & electron
     write(gtcout,"(I7,6e15.6)")istep+mstepall,(fieldtime(4,i),i=1,nfield),&
          (partdata(9,i),i=1,nspecies)
     call FLUSH(gtcout)
  endif

101 format(i6)
102 format(e13.6)

end subroutine diagnosis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine frc_diagnosis
  use global_parameters
  use field_array
  use particle_array
  use equilibrium
  use cylindricalRZ
  use LAMYRIDGEeq
  use spline_function,only: spline2d
  implicit none

  integer,parameter :: nfield=1,mfdiag=7
  integer :: i,j,m,ij,frcdiag,frcdata1d,frctracking,ndstep,icount,ierror,nf,nrsp,nzsp,ii
  real(lk) :: fieldtime(mfdiag,nfield),rdum,zdum,sdr,sdz,sdr_inv,sdz_inv,etotdum,bdum,psidum
  101 format(' ',4g16.4)
  102 format(' ',9g16.8)
  103 format(' ',8g16.4)

  frcdiag=654
  frcdata1d=987
  frctracking=624
  ndstep=mstep/ndiag
  mstepall=0
  if(mype==0 .and. istep==ndiag .and. irun==0)then
! open time history file in new run
     open(frcdiag,file='frc_history.out',status='replace')
     write(frcdiag,102)ndstep,ndiag
! open output file for radial-time data in new run
 !    open(iodata1d,file='data1d.out',status='replace')
 !    write(iodata1d,102)ndstep,mpsi+1,nspecies,nhybrid,mpdata1d,nfield,mfdata1d
  endif

  if(mype==0 .and. istep==ndiag .and. irun==0)then
! open time history file in new run
     open(frctracking,file='frc_tracking.out',status='replace')
  endif

  if(mype==0)then
    rdum=zelectron(1,1)
    zdum=zelectron(3,1)
    bdum = spline2d(rdum,zdum+zmax,0,lsr,lsz,spdR,spdZ,B0Mag_Eq,0,0)
    write(frctracking,102)istep,zelectron(1,1),zelectron(2,1),zelectron(3,1),zelectron(6,1),zelectron(8,1),zelectron(7,1),bdum
    call FLUSH(frctracking)
  endif

end subroutine frc_diagnosis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! calculate (n,m) mode amplitude for diagnostics
subroutine spectrum(nfield,modes,nmodes,mmodes,iflux,mtgrid,fieldgrid,fieldmode)
  use global_parameters
  use field_array,only: zeta0,deltaz,qtinv
  implicit none

  integer nfield,modes,nmodes(modes),mmodes(modes),iflux,mtgrid
  real(lk) fieldgrid(0:1,0:mtgrid,nfield),fieldmode(2,modes,nfield)
  integer i,j,k,jt,kz,mzeach,jpe,indp,indt,indp1,indt1,mteach,icount,ierror,nf
  real(lk) wt,dt,wz,zdum,tdum,xz(mtdiag),fieldflux(mtdiag/mtoroidal,mtdiag,nfield),&
       eachzeta((mtdiag/mtoroidal)*(mtdiag/mtoroidal)*nfield),&
       allzeta((mtdiag/mtoroidal)*mtdiag*nfield)
  complex(lk) yzeta(nfield,(mtdiag/mtoroidal)*modes),ytheta(nfield,mtdiag*modes),&
       yz(mtdiag/2+1),ye(mtdiag) 

! mesh in magnetic coordinates
  mzeach=mtdiag/mtoroidal
  mteach=mzeach
  dt=pi2/real(mtdiag)

! Interpolate on a flux surface from fieldline coordinates to magnetic coordinates. 
! Use mtdiag for both poloidal and toroidal grid points.
!$omp parallel do private(kz,nf,j,tdum,zdum,jt,wt,wz)           
  do kz=1,mzeach
     wz=real(kz)/real(mzeach)
     zdum=zeta0+wz*deltaz
     do nf=1,nfield
        do j=1,mtdiag
           if(fielddir==1 .or. fielddir==3)then
             tdum=pi2_inv*modulo(dt*real(j)-(zdum-pi2)*qtinv(iflux),pi2)*real(mtgrid)
           else
             tdum=pi2_inv*modulo(dt*real(j)-zdum*qtinv(iflux),pi2)*real(mtgrid)
           endif
           jt=max(0,min(mtgrid-1,int(tdum)))
           wt=tdum-real(jt)
           
           fieldflux(kz,j,nf)=wz*((1.0_lk-wt)*fieldgrid(1,jt,nf)+wt*fieldgrid(1,jt+1,nf))+&
                (1.0_lk-wz)*((1.0_lk-wt)*fieldgrid(0,jt,nf)+wt*fieldgrid(0,jt+1,nf))
        enddo
     enddo
  enddo
 
! transpose 2-d matrix from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  allzeta=0.0_lk
  icount=mteach*mzeach*nfield  
  do jpe=0,mtoroidal-1

!$omp parallel do private(j,nf,k,jt,indt,indp1,indp)
     do j=1,mteach
        jt=jpe*mteach+j
        indt=(j-1)*mzeach
        do nf=1,nfield
           indp1=indt+(nf-1)*mteach*mzeach
           do k=1,mzeach
              indp=indp1+k
              eachzeta(indp)=fieldflux(k,jt,nf)
           enddo
        enddo
     enddo
     
     call MPI_GATHER(eachzeta,icount,mpi_Rsize,allzeta,icount,mpi_Rsize,jpe,toroidal_comm,ierror)
  enddo
    
! transform to k space
  yz=0.0_lk
  do j=1,mteach
     indt1=(j-1)*mzeach
     do nf=1,nfield
        indt=indt1+(nf-1)*mteach*mzeach

!$omp parallel do private(kz,k,indp)      
        do kz=0,mtoroidal-1
           do k=1,mzeach
              indp=kz*icount+indt+k
              xz(kz*mzeach+k)=allzeta(indp)
           enddo
        enddo
      
        call fftr1d(1,mtdiag,1.0_lk,xz,yz,2)

! record toroidal mode amplitude for diagnostic
        do kz=1,modes
           yzeta(nf,j+mteach*(kz-1))=yz(nmodes(kz)+1)
        enddo

     enddo
  enddo

! gather toroidal mode amplitude for calculation of poloidal mode amplitudes
  ytheta=0.0_lk
  icount=mteach*modes*nfield
  call MPI_GATHER(yzeta,icount,mpi_Csize,ytheta,icount,mpi_Csize,0,toroidal_comm,ierror)

  icount=mteach*modes
  if(myrank_toroidal == 0)then
     do nf=1,nfield
        do kz=1,modes

!$omp parallel do private(i,j)     
           do i=0,mtoroidal-1
              do j=1,mteach
                 ye(j+i*mteach)=ytheta(nf,j+(kz-1)*mteach+i*icount)
              enddo
           enddo

           call fftc1d(1,mtdiag,1.0_lk,ye)
           
           if(fielddir==1 .or. fielddir==3)then
             fieldmode(1,kz,nf)=real(ye(mmodes(kz)+1))/(real(mtdiag)**2)
             fieldmode(2,kz,nf)=aimag(ye(mmodes(kz)+1))/(real(mtdiag)**2)
           else
             fieldmode(1,kz,nf)=real(ye(mtdiag-mmodes(kz)+1))/(real(mtdiag)**2)
             fieldmode(2,kz,nf)=aimag(ye(mtdiag-mmodes(kz)+1))/(real(mtdiag)**2)
           endif
        enddo
     enddo
  endif

end subroutine spectrum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! filter out some poloidal modes in linear simulation
subroutine mfilter_frc(farray)

  use global_parameters
  use field_array
  implicit none

  integer i,j,k,ii,ij,ip,jt,kz,ipe,kpe,jtp,&
       icount,ierror,idest,isource,l,isp,maxmmode
  real(lk) farray(0:1,mgrid),farray2d(mthetamax,mpsi),&
       ffilter(mtoroidal/2+1),xz(mtoroidal)
  complex(lk) yz(mtoroidal/2+1)

 !2d theta-psi array
 !$omp parallel do private(i,j,ij)
  do i=1,mpsi-1
     do j=1,mthetamax
        ij=igrid(i)+j
        farray2d(j,i)=farray(1,ij)
     enddo
  enddo

  ffilter=0.0_lk
  if(nfilter==2)then
      ffilter(mmodes+1)=1.0_lk
  endif

  yz=0.0_lk
  xz=0.0_lk
 
  !$omp parallel do private(i,xz,yz)
  do i=1,mpsi-1
       xz=farray2d(:,i)
       call fftr1d(1,mthetamax,1.0_lk,xz,yz,2)
       yz=yz*ffilter
       call fftr1d(-1,mthetamax,1.0_lk,xz,yz,2)
       farray2d(:,i)=xz
  enddo

  !2d theta-psi array
  !$omp parallel do private(i,j,ij)
  do i=1,mpsi-1
     do j=1,mthetamax
        ij=igrid(i)+j
        farray(1,ij)=farray2d(j,i)
     enddo
  enddo
 
  !$omp parallel do private(i)
  do i=1,mpsi-1
        farray(1,igrid(i))=farray(1,igrid(i)+mtheta(i)) ! poloidal BC
  enddo


end subroutine mfilter_frc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! filter out some toroidal (and optionally poloidal) modes in linear simulation
subroutine filter(farray)

  use global_parameters
  use field_array
  implicit none
  
  integer i,j,k,ii,ij,ip,jt,kz,jpe,indp,indt,indp1,indt1,meachtheta,mteach,mzeach,jtp,&
       icount,ierror,idest,isource,isendtag,l,irecvtag,joffset,mteachpd,&
       istatus(MPI_STATUS_SIZE)
  real(lk) wt,r,dt,wz,zdum,tdum,ptemp,sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),&
       phism(mgrid),pleft(mthetamax),pright(mthetamax),phitmp(mgrid),ffilter(mtdiag/2+1),&
       xz(mtdiag),farray(0:1,mgrid),&
       phifluxpartd(mtdiag/mtoroidal,mpsi,mtdiag/npartdom),&
       phifluxlz(mpsi,mtdiag/numberpe,mtdiag),&
       phiflux(mpsi,mtdiag),&
       allzeta((mtdiag/numberpe)*mtdiag*mpsi)
  complex(lk) yz(mtdiag/2+1)
     
! periodic BC in poloidal and toroidal
  CALL PERIODICITY(farray)
!by zhs nfilter > 2, use zhs's filter
  if(nfilter > 2)then  
     call mnfilter(farray)
     return
  endif

! mesh in magnetic coordinates
  mzeach=mtdiag/mtoroidal
  mteach=mzeach/npartdom
  mteachpd=(mtdiag/npartdom)
  dt=2.0_lk*pi/real(mtdiag)
  
  ffilter=0.0_lk
  ffilter(nmodes+1)=1.0_lk
  allzeta=0.0_lk

  joffset=mteachpd*myrank_partd

! Interpolate on a flux surface from fieldline coordinates to magnetic
! coordinates. Use mtdiag for both poloidal and toroidal grid points.
!$omp parallel do private(kz,i,j,wz,ii,zdum,tdum,jt,wt)           
  do kz=1,mzeach
     wz=real(kz)/real(mzeach)
     zdum=zeta0+wz*deltaz
     do i=1,mpsi
        ii=igrid(i)
        do j=1,mteachpd
           if(fielddir==1 .or. fielddir==3)then
             tdum=pi2_inv*modulo(dt*real(j+joffset)-(zdum-pi2)*qtinv(i),pi2)*real(mtheta(i))
           else
             tdum=pi2_inv*modulo(dt*real(j+joffset)-zdum*qtinv(i),pi2)*real(mtheta(i))
           endif
           jt=max(0,min(mtheta(i)-1,int(tdum)))
           wt=tdum-real(jt)
           phifluxpartd(kz,i,j)=((1.0_lk-wt)*farray(1,ii+jt)+wt*farray(1,ii+jt+1))*&
                wz+(1.0_lk-wz)*((1.0_lk-wt)*farray(0,ii+jt)+wt*farray(0,ii+jt+1))
        enddo
     enddo
  enddo

! transpose 2-d matrix from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  allzeta=0.0_lk
  icount=mteach*mzeach*mpsi
 !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp) 
   do jpe=0,mtoroidal-1
      do j=1,mteach
         jt=jpe*mteach+j
         indt=(j-1)*mpsi
         do k=1,mzeach
            indp1=indt+(k-1)*mpsi*mteach
            do i=1,mpsi
               indp=indp1+i
               allzeta(indp+jpe*icount)=phifluxpartd(k,i,jt)
            enddo
         enddo
      enddo
   enddo

  call MPI_ALLTOALL(allzeta,icount,mpi_Rsize,phifluxlz,icount,&
          mpi_Rsize,toroidal_comm,ierror)
  
! transform to k space
  yz=0.0_lk
!$omp parallel do private(j,i,kz,k,indt1,indt,indp,xz,yz)
  do j=1,mteach
     !indt1=(j-1)*mzeach
     do i=1,mpsi
        xz=phifluxlz(i,j,:)
        call fftr1d(1,mtdiag,1.0_lk,xz,yz,2)
! linear run only keep a few modes
        yz=ffilter*yz
! transform back to real space
        call fftr1d(-1,mtdiag,1.0_lk,xz,yz,2)
        phifluxlz(i,j,:)=xz
     enddo
  enddo

  call MPI_ALLTOALL(phifluxlz,icount,mpi_Rsize,allzeta,&
          icount,mpi_Rsize,toroidal_comm,ierror)

!$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
  do jpe=0,mtoroidal-1
    do j=1,mteach
      jt=jpe*mteach+j
      indt=(j-1)*mpsi
      do k=1,mzeach
        indp1=indt+(k-1)*mpsi*mteach
        do i=1,mpsi
          indp=indp1+i
          phifluxpartd(k,i,jt)=allzeta(indp+jpe*icount)
        enddo
      enddo
    enddo
  enddo
  
  call MPI_ALLGATHER(phifluxpartd(mzeach,:,:),mpsi*mteachpd,mpi_Rsize,phiflux,mpsi*mteachpd,mpi_Rsize,partd_comm,ierror)

  if(nfilter==2)then ! filter poloidal modes
     ffilter=0.0_lk
     ffilter(mmodes+1)=1.0_lk
  
     do i=1,mpsi
        do j=1,mtdiag
           xz(j)=phiflux(i,j)
        enddo
        call fftr1d(1,mtdiag,1.0_lk,xz,yz,2)
        yz=ffilter*yz
        call fftr1d(-1,mtdiag,1.0_lk,xz,yz,2)
        do j=1,mtdiag
           phiflux(i,j)=xz(j)
        enddo
     enddo
  endif

! interpolate field from magnetic coordinates to fieldline coordinates
!$omp parallel do private(i,j,ii,tdum,jt,wt,jtp)           
  do i=1,mpsi
     ii=igrid(i)              
     do j=1,mtheta(i)
        if(fielddir==1 .or. fielddir==3)then
          tdum=pi2_inv*modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)*real(mtdiag)
        else
          tdum=pi2_inv*modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)*real(mtdiag)
        endif
        jt=max(0,min(mtdiag-1,int(tdum)))
        wt=tdum-real(jt)
        jtp=jt+1
        if(jt==0)jt=mtdiag
           
        farray(1,ii+j)=wt*phiflux(i,jtp)+(1.0-wt)*phiflux(i,jt)
     enddo
  enddo

  call periodicity(farray)

end subroutine filter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mnfilter(farray)

  use global_parameters
  use field_array
  implicit none

  integer i,j,k,ii,ij,ip,jt,kz,jpe,indp,indt,indp1,indt1,meachtheta,mteach,mzeach,jtp,&
       icount,ierror,idest,isource,isendtag,l,irecvtag,joffset,mteachpd,&
       istatus(MPI_STATUS_SIZE)
  real(lk) wt,r,dt,wz,zdum,tdum,ptemp,sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),&
       phism(mgrid),pleft(mthetamax),pright(mthetamax),phitmp(mgrid),ffilter(mtdiag/2+1),&
       xz(mtdiag),farray(0:1,mgrid),&
       phifluxpartd(mtdiag/mtoroidal,mpsi,mtdiag/npartdom),&
       phifluxlz(mpsi,mtdiag/numberpe,mtdiag),&
       phiflux(mpsi,mtdiag),&
       allzeta((mtdiag/numberpe)*mtdiag*mpsi)

!by zhs
  integer localpe,localmt
  real(lk) field2d(mtdiag,mtdiag,mpsi/mtoroidal),ffiltern(mtdiag/2+1),ffilterm(mtdiag),ffiltern0(mtdiag)

  complex(lk) yz(mtdiag/2+1),zz(mtdiag),cfield2d1(mpsi,mtdiag/numberpe,mtdiag),&
  cfield2d2(mpsi,mtdiag/numberpe,mtdiag),cphifluxpartd(mtdiag/mtoroidal,mpsi,mtdiag/npartdom),&
  allcfield(mtdiag*mtdiag*mpsi/numberpe)

! mesh in magnetic coordinates
  mzeach=mtdiag/mtoroidal
  mteach=mzeach/npartdom
  mteachpd=(mtdiag/npartdom)
  dt=2.0_lk*pi/real(mtdiag)

  ffilter=0.0_lk
  ffiltern=0.0_lk
  ffilterm=0.0_lk
  ffiltern0=0.0_lk
  ffilter(nmodes+1)=1.0_lk
  ffiltern(nmodes(1)+1)=1.0_lk
!  ffiltern(1)=1.0    !keep the n=0 mode
!  ffiltern0(mtdiag-1)=1.0 !keep a few m mode for n=0
!  ffiltern0(3)=1.0
  ffilterm(mtdiag+1-mmodes(1))=1.0_lk
  if(nfilter==4)ffilterm((mtdiag-mmodes(1)):(mtdiag+2-mmodes(1)))=1.0_lk
  if(nfilter==5)ffilterm((mtdiag-1-mmodes(1)):mtdiag+3-mmodes(1))=1.0_lk

!by zhs
  localmt=mod(nmodes(1)+1,mteach)
  localpe=(nmodes(1)+1-localmt)/mteach
  if(localmt == 0)then
     localmt=mteach
     localpe=localpe-1
  endif

  allzeta=0.0_lk

  joffset=mteachpd*myrank_partd

! Interpolate on a flux surface from fieldline coordinates to magnetic
! coordinates. Use mtdiag for both poloidal and toroidal grid points.
!$omp parallel do private(kz,i,j,wz,ii,zdum,tdum,jt,wt)
  do kz=1,mzeach
     wz=real(kz)/real(mzeach)
     zdum=zeta0+wz*deltaz
     do i=1,mpsi
        ii=igrid(i)
        do j=1,mteachpd
           if(fielddir==1 .or. fielddir==3)then
             tdum=pi2_inv*modulo(dt*real(j+joffset)-(zdum-pi2)*qtinv(i),pi2)*real(mtheta(i))
           else
             tdum=pi2_inv*modulo(dt*real(j+joffset)-zdum*qtinv(i),pi2)*real(mtheta(i))
           endif
           jt=max(0,min(mtheta(i)-1,int(tdum)))
           wt=tdum-real(jt)
           phifluxpartd(kz,i,j)=((1.0_lk-wt)*farray(1,ii+jt)+wt*farray(1,ii+jt+1))*&
                wz+(1.0_lk-wz)*((1.0_lk-wt)*farray(0,ii+jt)+wt*farray(0,ii+jt+1))
        enddo
     enddo
  enddo

! transpose 2-d matrix from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  allzeta=0.0_lk
  icount=mteach*mzeach*mpsi
 !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
   do jpe=0,mtoroidal-1
      do j=1,mteach
         jt=jpe*mteach+j
         indt=(j-1)*mpsi
         do k=1,mzeach
            indp1=indt+(k-1)*mpsi*mteach
            do i=1,mpsi
               indp=indp1+i
               allzeta(indp+jpe*icount)=phifluxpartd(k,i,jt)
            enddo
         enddo
      enddo
   enddo

  call MPI_ALLTOALL(allzeta,icount,mpi_Rsize,phifluxlz,icount,&
          mpi_Rsize,toroidal_comm,ierror)

! transform to k space
  yz=0.0_lk
  cfield2d1=0.0_lk
!$omp parallel do private(j,i,kz,k,indt1,indt,indp,xz,yz)
  do j=1,mteach
     do i=1,mpsi
        xz=phifluxlz(i,j,:)
        call fftr1d(1,mtdiag,1.0_lk,xz,yz,2)
! linear run only keep a few modes
        yz=ffiltern*yz
        do k=1,mtdiag
           if(k<mtdiag/2+2)then
              cfield2d1(i,j,k)=yz(k)
           else
              cfield2d1(i,j,k)=0.0_lk
           endif
        enddo
     enddo
  enddo

  allcfield=0.0_lk  !similar to allzeta structure
  call MPI_ALLTOALL(cfield2d1,icount,mpi_complex,allcfield,&
          icount,mpi_complex,toroidal_comm,ierror)

!transform back to (mtoroidal,mzetach)
  cphifluxpartd=0.0_lk  !similar to phifluxpartd
  !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
  do jpe=0,mtoroidal-1
    do j=1,mteach
      jt=jpe*mteach+j
      indt=(j-1)*mpsi
      do k=1,mzeach
        indp1=indt+(k-1)*mpsi*mteach
        do i=1,mpsi
          indp=indp1+i
          cphifluxpartd(k,i,jt)=allcfield(indp+jpe*icount)
        enddo
      enddo
    enddo
  enddo

! transpose 2-d matrix from (mtoroidal,mzeach) to (mtoroidal*mzeach,1)
  allcfield=0.0_lk
  icount=mteach*mteachpd*mpsi !here replace mzeach by mteachpd
 !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
   do jpe=0,npartdom-1      ! zeta in one pe is divided by jpe(kpe)
      do j=1,mteach         ! zeta in one jpe is mteach(j as k)
         jt=jpe*mteach+j  
         indt=(j-1)*mpsi
         do k=1,mteachpd    ! theta in one partdom (k as j)
            indp1=indt+(k-1)*mpsi*mteach
            do i=1,mpsi
               indp=indp1+i
               allcfield(indp+jpe*icount)=cphifluxpartd(jt,i,k)
            enddo
         enddo
      enddo
   enddo

  cfield2d1=0.0_lk
  call MPI_ALLTOALL(allcfield,icount,mpi_complex,cfield2d1,&
          icount,mpi_complex,partd_comm,ierror)

! transform from k to k space
  zz=0.0_lk
  cfield2d2=0.0_lk
!~~~~~~~~~~~~~~~~~~new scheme for fft in 2d,valid for multi mode~~~~~~~~~~~~~
  if(mype==localpe)then
     j=localmt
!$omp parallel do private(i,k,zz)
     do i=1,mpsi
        do k=1,mtdiag
           zz(k)=cfield2d1(i,j,k)
        enddo
!fft in m direction
        call fftc1d(1,mtdiag,1.0_lk,zz)
        do k=1,mtdiag
           zz(k)=ffilterm(k)*zz(k)
        enddo
! transform back to k_n space
        call fftc1d(-1,mtdiag,1.0_lk,zz)
        do k=1,mtdiag
           cfield2d2(i,j,k)=zz(k)
        enddo
     enddo
  endif
!to control the n=0,m/=0 harmonics in NL simulation
  if(mype==0)then
    j=1          
!$omp parallel do private(i,k,zz)
     do i=1,mpsi
        do k=1,mtdiag
           zz(k)=cfield2d1(i,j,k)
        enddo
!fft in m direction
        call fftc1d(1,mtdiag,1.0_lk,zz)
        do k=1,mtdiag
           zz(k)=ffiltern0(k)*zz(k)
        enddo
! transform back to k_n space
        call fftc1d(-1,mtdiag,1.0_lk,zz)
        do k=1,mtdiag
           cfield2d2(i,j,k)=zz(k)
        enddo
     enddo
  endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!~~~~~~~~~~~~~~~original scheme for fft in 2d,valid for single n mode~~~~~~
!!$omp parallel do private(j,i,k,zz)
!  do j=1,mteach
!     do i=1,mpsi
!        do k=1,mtdiag
!           zz(k)=cfield2d1(i,j,k)
!        enddo
!!fft in m direction
!        call fftc1d(1,mtdiag,1.0,zz)
!        do k=1,mtdiag
!           zz(k)=ffilterm(k)*zz(k)
!        enddo
!! transform back to k_n space
!        call fftc1d(-1,mtdiag,1.0,zz)
!        do k=1,mtdiag
!           cfield2d2(i,j,k)=zz(k)
!        enddo
!     enddo
!  enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!transform back from (mtoroidal*mzeach,1) to (mtoroidal,mzeach)
  allcfield=0.0_lk
  call MPI_ALLTOALL(cfield2d2,icount,mpi_complex,allcfield,&
          icount,mpi_complex,partd_comm,ierror)

 !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
   do jpe=0,npartdom-1      ! zeta in one pe is divided by jpe(kpe)
      do j=1,mteach         ! zeta in one jpe is mteach(j as k)
         jt=jpe*mteach+j
         indt=(j-1)*mpsi
         do k=1,mteachpd    ! theta in one partdom (k as j)
            indp1=indt+(k-1)*mpsi*mteach
            do i=1,mpsi
               indp=indp1+i
               cphifluxpartd(jt,i,k)=allcfield(indp+jpe*icount)
            enddo
         enddo
      enddo
   enddo

!transform back from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  icount=mteach*mzeach*mpsi
  !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
  do jpe=0,mtoroidal-1
    do j=1,mteach
      jt=jpe*mteach+j
      indt=(j-1)*mpsi
      do k=1,mzeach
        indp1=indt+(k-1)*mpsi*mteach
        do i=1,mpsi
          indp=indp1+i
          allcfield(indp+jpe*icount)=cphifluxpartd(k,i,jt)
        enddo
      enddo
    enddo
  enddo

  cfield2d1=0.0_lk
  call MPI_ALLTOALL(allcfield,icount,mpi_complex,cfield2d1,&
          icount,mpi_complex,toroidal_comm,ierror)

! transform back to real space
  xz=0.0_lk
!$omp parallel do private(j,i,kz,k,xz,yz)
  do j=1,mteach
     do i=1,mpsi
        do k=1,mtdiag/2+1
           yz(k)=cfield2d1(i,j,k)
        enddo
        call fftr1d(-1,mtdiag,1.0_lk,xz,yz,2)
        phifluxlz(i,j,:)=xz
     enddo
  enddo

  call MPI_ALLTOALL(phifluxlz,icount,mpi_Rsize,allzeta,icount,&
          mpi_Rsize,toroidal_comm,ierror)

! transpose 2-d matrix from (1,mtoroidal*mzeach) to (mtoroidal,mzeach)
 !$omp parallel do private(jpe,j,i,k,jt,indt,indp1,indp)
   do jpe=0,mtoroidal-1
      do j=1,mteach
         jt=jpe*mteach+j
         indt=(j-1)*mpsi
         do k=1,mzeach
            indp1=indt+(k-1)*mpsi*mteach
            do i=1,mpsi
               indp=indp1+i
               phifluxpartd(k,i,jt)=allzeta(indp+jpe*icount)
            enddo
         enddo
      enddo
   enddo

  call MPI_ALLGATHER(phifluxpartd(mzeach,:,:),mpsi*mteachpd,mpi_Rsize,phiflux,mpsi*mteachpd,mpi_Rsize,partd_comm,ierror)

! interpolate field from magnetic coordinates to fieldline coordinates
!$omp parallel do private(i,j,ii,tdum,jt,wt,jtp)
  do i=1,mpsi
     ii=igrid(i)
     do j=1,mtheta(i)
        if(fielddir==1 .or. fielddir==3)then
          tdum=pi2_inv*modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)*real(mtdiag)
        else
          tdum=pi2_inv*modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)*real(mtdiag)
        endif
        jt=max(0,min(mtdiag-1,int(tdum)))
        wt=tdum-real(jt)
        jtp=jt+1
        if(jt==0)jt=mtdiag

        farray(1,ii+j)=wt*phiflux(i,jtp)+(1.0_lk-wt)*phiflux(i,jt)
     enddo
  enddo

  call periodicity(farray)

end subroutine mnfilter

subroutine mfilter(farray)
  use global_parameters
  use field_array
  implicit none

  integer i,j,k,ii,ij,ip,jt,kz,ipe,kpe,indp,indt,indp1,indt1,meachtheta,mteach,mzeach,jtp,&
       icount,ierror,idest,isource,isendtag,l,irecvtag,istatus(MPI_STATUS_SIZE)
  real(lk) wt,r,dt,wz,zdum,tdum,ptemp,sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),&
       phism(mgrid),pleft(mthetamax),pright(mthetamax),phitmp(mgrid),ffilter(mtdiag/2+1),&
       xz(mtdiag),farray(0:1,mgrid),phiflux(mtdiag/mtoroidal,mtdiag,mpsi),&
       allzeta2(mtdiag*mtdiag*mpsi/mtoroidal),allzeta(mtdiag*mtdiag*mpsi/mtoroidal)
       
!by zhs
  integer mpeach
  real(lk) field2d(mtdiag,mtdiag,mpsi/mtoroidal),ffilters(mtdiag/2+1,mtdiag)

  complex(lk) yz(mtdiag/2+1),zz(mtdiag),cfield2d(mtdiag/2+1,mtdiag)

! mesh in magnetic coordinates
  mzeach=mtdiag/mtoroidal
  mteach=mzeach
  mpeach=mpsi/mtoroidal

  dt=2.0_lk*pi/real(mtdiag)

  ffilter=0.0_lk
  ffilters=0.0_lk
!by zhs use the first mode. when nfilter>3, mult m
  ffilter(nmodes(1)+1)=1.0_lk
  ffilters(nmodes(1)+1,mtdiag+1-mmodes(1))=1.0_lk
  if(nfilter>3)then  !nfilter>3, include m-1,m+1 in the simulation
    ffilters(nmodes(1)+1,mtdiag-mmodes(1))=1.0_lk
    ffilters(nmodes(1)+1,mtdiag+2-mmodes(1))=1.0_lk
  endif
!  ffilters(nmodes(1)+1,mtdiag-1-mmodes(1))=1.0
!  ffilters(nmodes(1)+1,mtdiag+3-mmodes(1))=1.0

  field2d=0.0_lk
  cfield2d=0.0_lk
! Interpolate on a flux surface from fieldline coordinates to magnetic
! coordinates. Use mtdiag for both poloidal and toroidal grid points.
!$omp parallel do private(kz,i,j,wz,ii,zdum,tdum,jt,wt)
  do kz=1,mzeach
     wz=real(kz)/real(mzeach)
     zdum=zeta0+wz*deltaz
     do i=1,mpsi
        ii=igrid(i)
        do j=1,mtdiag
           if(fielddir==1 .or. fielddir==3)then
             tdum=pi2_inv*modulo(dt*real(j)-(zdum-pi2)*qtinv(i),pi2)*real(mtheta(i))
           else
             tdum=pi2_inv*modulo(dt*real(j)-zdum*qtinv(i),pi2)*real(mtheta(i))
           endif
           jt=max(0,min(mtheta(i)-1,int(tdum)))
           wt=tdum-real(jt)

           phiflux(kz,j,i)=((1.0_lk-wt)*farray(1,ii+jt)+wt*farray(1,ii+jt+1))*&
                wz+(1.0_lk-wz)*((1.0_lk-wt)*farray(0,ii+jt)+wt*farray(0,ii+jt+1))
        enddo
     enddo
  enddo

! transpose 2-d matrix from (mtoroidal,mzeach) to (1,mtoroidal*mzeach)
  allzeta=0.0_lk
  icount=mtdiag*mzeach*mpeach

!$omp parallel do private(ipe,j,i,k,ip,indp,indp1,indt)
  do ipe=0,mtoroidal-1
     do i=1,mpeach
        ip=mpeach*ipe+i
        indp1=(i-1)*mtdiag*mzeach
        do j=1,mtdiag
           indt=indp1+(j-1)*mzeach
           do k=1,mzeach
              indp=indt+k
              allzeta(indp+ipe*icount)=phiflux(k,j,ip)
           enddo
       enddo
     enddo
  enddo
     call MPI_ALLTOALL(allzeta,icount,mpi_Rsize,allzeta2,icount,&
          mpi_Rsize,toroidal_comm,ierror)


!transpose 1-d array (mtdiag*mtdiag,mpeach) allzeta to 2-d array field2d(mtdiag,mtdiag,mpeach)
!$omp parallel do private(j,i,kz,k,kpe,indp,indp1,indt)
  do kpe=0,mtoroidal-1
     do i=1,mpeach
        indp1=kpe*icount+(i-1)*mtdiag*mzeach
        do j=1,mtdiag
           indt=indp1+(j-1)*mzeach
           do k=1,mzeach
              indp=indt+k
              kz=kpe*mzeach+k
              field2d(kz,j,i)=allzeta2(indp)
           enddo
        enddo
     enddo
  enddo

! transform to k space
  yz=0.0_lk
  zz=0.0_lk

!filter n and m mode
!$omp parallel do private(j,i,k,xz,yz,zz,cfield2d)
     do i=1,mpeach
        do j=1,mtdiag
           do k=1,mtdiag
                 xz(k)=field2d(k,j,i)
           enddo
! fft in n direction
           call fftr1d(1,mtdiag,1.0_lk,xz,yz,2)
           do k=1,mtdiag/2+1
              cfield2d(k,j)=yz(k)
           enddo
        enddo
        do k=1,mtdiag/2+1
           do j=1,mtdiag
              zz(j)=cfield2d(k,j)
           enddo
! fft in m direction
           call fftc1d(1,mtdiag,1.0_lk,zz)
           do j=1,mtdiag
              cfield2d(k,j)=zz(j)
           enddo
        enddo

        cfield2d(:,:)=ffilters(:,:)*cfield2d(:,:)
! tranform back to real space
        do k=1,mtdiag/2+1
           do j=1,mtdiag
              zz(j)=cfield2d(k,j)
           enddo
           call fftc1d(-1,mtdiag,1.0_lk,zz)
           do j=1,mtdiag
              cfield2d(k,j)=zz(j)
           enddo
        enddo

        do j=1,mtdiag
           do k=1,mtdiag/2+1
              yz(k)=cfield2d(k,j)
           enddo
           call fftr1d(-1,mtdiag,1.0_lk,xz,yz,2)
           do k=1,mtdiag
              field2d(k,j,i)=xz(k)
           enddo
        enddo
     enddo

! tranpose from field2d back to allzeta
!$omp parallel do private(j,i,kz,k,kpe,indp,indt,indp1)
  do kpe=0,mtoroidal-1
     do i=1,mpeach
        indp1=kpe*icount+(i-1)*mtdiag*mzeach
        do j=1,mtdiag
           indt=indp1+(j-1)*mzeach
           do k=1,mzeach
              indp=indt+k
              kz=kpe*mzeach+k
              allzeta2(indp)=field2d(kz,j,i)
           enddo
        enddo
     enddo
  enddo

  call MPI_ALLTOALL(allzeta2,icount,mpi_Rsize,allzeta,&
          icount,mpi_Rsize,toroidal_comm,ierror)

!$omp parallel do private(ipe,j,i,k,ip,indp,indt,indp1)
  do ipe=0,mtoroidal-1
     do i=1,mpeach
        ip=mpeach*ipe+i
        indp1=(i-1)*mtdiag*mzeach
        do j=1,mtdiag
           indt=indp1+(j-1)*mzeach
           do k=1,mzeach
              indp=indt+k
              phiflux(k,j,ip)=allzeta(indp+ipe*icount)
           enddo
        enddo
     enddo
  enddo

! interpolate field from magnetic coordinates to fieldline coordinates
!$omp parallel do private(i,j,ii,tdum,jt,wt,jtp)
  do i=1,mpsi
     ii=igrid(i)
     do j=1,mtheta(i)
        if(fielddir==1 .or. fielddir==3)then
          tdum=pi2_inv*modulo(deltat(i)*real(j)+(zeta1-pi2)*qtinv(i),pi2)*real(mtdiag)
        else
          tdum=pi2_inv*modulo(deltat(i)*real(j)+zeta1*qtinv(i),pi2)*real(mtdiag)
        endif
        jt=max(0,min(mtdiag-1,int(tdum)))
        wt=tdum-real(jt)
        jtp=jt+1
        if(jt==0)jt=mtdiag

        farray(1,ii+j)=wt*phiflux(mzeach,jtp,i)+(1.0_lk-wt)*phiflux(mzeach,jt,i)
     enddo
  enddo

  call periodicity(farray)

end subroutine mfilter
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine opendiag(mpdiag,nfield,mfdiag,modes,mpdata1d,mfdata1d)
  use global_parameters
  implicit none

  integer mpdiag,nfield,mfdiag,modes,mpdata1d,mfdata1d
  integer i,j,iotmp,ndstep,ndata
  real(lk) tdum
  character(len=100) cdum

  iodiag=123
  iotmp=456
  iodata1d=789

! # of time steps and # of data added to history file
  ndstep=mstep/ndiag
  ndata=(nspecies*mpdiag+nfield*(2*modes+mfdiag))
  if(irun==0)then
     mstepall=0
! open time history file in new run
     open(iodiag,file='history.out',status='replace')
     write(iodiag,101)ndstep,nspecies,mpdiag,nfield,modes,mfdiag
     write(iodiag,102)tstep*ndiag

! open output file for radial-time data in new run
     open(iodata1d,file='data1d.out',status='replace')            
     write(iodata1d,101)ndstep,mpsi+1,nspecies,nhybrid,mpdata1d,nfield,mfdata1d
  else

! find history file from previous run
     irest=irest-1
     if(mod(irest,2)==0)then
        cdum="restart_dir1/history_restart.out"
     else
        cdum="restart_dir2/history_restart.out"
     endif

     open(iotmp,file=trim(cdum),status='old')
     rewind(iotmp)

! open time history file for restart run
     open(iodiag,file='history.out',status='replace')

! # of time steps
     read(iotmp,101)mstepall
     write(iodiag,101)mstepall+ndstep
     do i=1,5
        read(iotmp,101)j
        write(iodiag,101)j
     enddo
     
! copy restart history data file
     do i=0,ndata*mstepall
        read(iotmp,102)tdum
        write(iodiag,102)tdum
     enddo
     close(iotmp)

!!copy restart data1d.out file
     open(iodata1d,file='data1d.out',status='replace')
     if(mod(irest,2)==0)then
        cdum="restart_dir1/data1d_restart.out"
     else
        cdum="restart_dir2/data1d_restart.out"
     endif
     
     irest=irest+1
     
     open(iotmp,file=trim(cdum),status='old')
     rewind(iotmp)
     read(iotmp,101)j
     write(iodata1d,101)mstepall+ndstep
     do i=1,6
        read(iotmp,101)j
        write(iodata1d,101)j
     enddo

     ndata=(mpsi+1)*(nspecies*mpdata1d+nfield*mfdata1d)
     do i=1,ndata*mstepall
        read(iotmp,102)tdum
        write(iodata1d,102)tdum
     enddo
     close(iotmp)

!!calculate current total time step        
     mstepall=mstepall*ndiag
     
  endif

101 format(i6)
102 format(e13.6)

end subroutine opendiag
