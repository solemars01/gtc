!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!the unformatted bin file
subroutine restart_native(iop)
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer ierror
  character(*),intent(in)::iop
  character(len=80) :: restart_fname
  character(len=18) cdum

  if(mype < 10)then
     write(cdum,'("DATA_RESTART.0000",i1)')mype
  elseif(mype < 100)then
     write(cdum,'("DATA_RESTART.000",i2)')mype
  elseif(mype < 1000)then
     write(cdum,'("DATA_RESTART.00",i3)')mype
  elseif(mype < 10000)then
     write(cdum,'("DATA_RESTART.0",i4)')mype
  else
     write(cdum,'("DATA_RESTART.",i5)')mype
  endif

!!XY rolling restart: save two copies of restart data
  if(mod(irest,2)==0)then
     restart_fname="restart_dir1/"//trim(cdum)
  else
     restart_fname="restart_dir2/"//trim(cdum)
  endif
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

! record particle information for future restart run
  if(iop=="write")then
     open(222,file=trim(restart_fname),status='replace',form='unformatted')
     write(222)mi,me,rdtemi,pfluxi,phi,phip00,phi00,zonali,zonalci,meshti,meshte,meshni,meshne,pmarki,markeri,markerit,&
          zion(:,1:mi)
     if(nhybrid>0)write(222)rdteme,pfluxe,zonale,zonalce,phisave,dnesave,zelectron(:,1:me),&
          pmarke,markere,tfrac,tfracn,&
          densitye,pressureepara,pressureeperp
     if(fload>0)write(222)mf,rdtemf,zonalf,zonalcf,meshtf,meshnf,pfluxf,zfast(:,1:mf),&
          pmarkf,markerf
     if(feload>0)write(222)mfe,rdtemfe,zonalfe,zonalcfe,meshtfe,meshnfe,pfluxfe,zfaste(:,1:mfe),&
          pmarkfe,markerfe,fetfrac,fetfracn
     if(magnetic>0)write(222)apara,sapara,sdelapara,fluidne,sfluidne,fluidue,deltapsi,sdeltapsi,apara00,apara00nl,fluidne00,gradgaugef
     close(222)
     if(mype==0)write(gtcout,*)'write to ',trim(restart_fname)
     call restart_hist

! read particle information to restart previous run
  else
     if(mype==0)write(gtcout,*)'read in ',trim(restart_fname)
     open(333,file=trim(restart_fname),status='old',form='unformatted')
     read(333)mi,me,rdtemi,pfluxi,phi,phip00,phi00,zonali,zonalci,meshti,meshte,meshni,meshne,pmarki,markeri,markerit,&
          zion(:,1:mi)
     if(nhybrid>0)read(333)rdteme,pfluxe,zonale,zonalce,phisave,dnesave,zelectron(:,1:me),&
          pmarke,markere,tfrac,tfracn,&
          densitye,pressureepara,pressureeperp
     if(fload>0)read(333)mf,rdtemf,zonalf,zonalcf,meshtf,meshnf,pfluxf,zfast(:,1:mf),&
          pmarkf,markerf
     if(feload>0)read(333)mfe,rdtemfe,zonalfe,zonalcfe,meshtfe,meshnfe,pfluxfe,zfaste(:,1:mfe),&
          pmarkfe,markerfe,fetfrac,fetfracn
     if(magnetic>0)read(333)apara,sapara,sdelapara,fluidne,sfluidne,fluidue,deltapsi,sdeltapsi,apara00,apara00nl,fluidne00,gradgaugef
     close(333)
     irest=irest+1 !!XY rolling rstart
     if(mype==0)write(gtcout,*)restart_fname,'read over'

  endif

end subroutine restart_native
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine restart_io(iop)
  use global_parameters
  use particle_array
  use field_array
  use data_type
  use particle_tracking
  implicit none
  integer ierror, mrequest, mfmode
  integer i,j,k,m,subsize,startidx,endidx,ndata
  integer mquantity,mflx,n_mode,mstepfinal,noutputs
  real(lk) dum
  character(*),intent(in)::iop
  character(len=80) :: restart_fname
  character(len=18) cdum
  integer :: save_restart_files
  INTEGER(KIND=MPI_OFFSET_KIND) mype_filesize, sum_filesize

#if ADIOS
#define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
#define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
#define ADIOS_READ(a,b) call adios_read(a,'b'//char(0),b)
  integer*8 buf_id, group_id
  real(lk),dimension(:),allocatable::zion0_read,zelectron0_read
  character(len=80) :: dirstr
#endif

  if(iop/="read" .and. iop/="write")then
     write(gtcout,*)'*** subroutine restart_io (iop <> "read" or "write")',iop
     call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
     return
  endif

#if ADIOS
!!xy rolling restart
  if(mod(irest,2)==0)then
     write(restart_fname,'(a,i5.5,".bp")')"restart_dir1/restart_",myrank_toroidal
  else
     write(restart_fname,'(a,i5.5,".bp")')"restart_dir2/restart_",myrank_toroidal
  endif

! setup the element path for this node
  write(dirstr,'("/node",i5.5,"/param")')mype
  dirstr=trim(dirstr)//char(0)
!     call MPI_BARRIER(MPI_COMM_WORLD,ierror)
!     start_time = mpi_wtime()
  call adios_get_group (group_id, "restart"//char(0))
! set the path for all vars in the type for proper sizing
  call adios_set_path (group_id,dirstr//char(0));
  restart_fname=trim(restart_fname)//char(0)
  if (iop=="read") then
     call adios_open_read (buf_id, group_id, restart_fname)
  else
     call adios_open (buf_id, group_id, restart_fname)
  endif
! write the sizing paramters for both reading and writing
  ADIOS_WRITE(buf_id,partd_comm)
  ADIOS_WRITE(buf_id,mpsi+1)
     !!ADIOS_WRITE(buf_id,mzeta+1)
  ADIOS_WRITE(buf_id,nparami)
  ADIOS_WRITE(buf_id,mimax)
  ADIOS_WRITE(buf_id,mgrid)
  ADIOS_WRITE(buf_id,memax)
  ADIOS_WRITE(buf_id,nhybrid)
  ADIOS_WRITE(buf_id,2*nhybrid)

  if(iop=="write")then
     !!ADIOS_WRITE(buf_id,mzeta)
     ADIOS_WRITE(buf_id,mi)
     ADIOS_WRITE(buf_id,me)
     ADIOS_WRITE(buf_id,rdtemi)
     ADIOS_WRITE(buf_id,pfluxi)
     ADIOS_WRITE(buf_id,phi00)
     ADIOS_WRITE(buf_id,phip00)
     ADIOS_WRITE(buf_id,zonali)
     ADIOS_WRITE(buf_id,meshti)
     ADIOS_WRITE(buf_id,meshte)
     ADIOS_WRITE(buf_id,meshni)
     ADIOS_WRITE(buf_id,meshne)
     ADIOS_WRITE(buf_id,phi)
     ADIOS_WRITE(buf_id,zion)

     ADIOS_WRITE(buf_id,pmarki)
     ADIOS_WRITE(buf_id,markeri)
     ADIOS_WRITE(buf_id,markerit)


     if(nhybrid>0)then
        ADIOS_WRITE(buf_id,zelectron)
        ADIOS_WRITE(buf_id,phisave)
        ADIOS_WRITE(buf_id,rdteme)
        ADIOS_WRITE(buf_id,pfluxe)
        ADIOS_WRITE(buf_id,zonale)

        ADIOS_WRITE(buf_id,pmarke)
        ADIOS_WRITE(buf_id,markere)

        ADIOS_WRITE(buf_id,tfrac)
        ADIOS_WRITE(buf_id,tfracn)

     endif

     if(magnetic>0)then
        ADIOS_WRITE(buf_id,apara)
        ADIOS_WRITE(buf_id,fluidne)
     endif
     if(fload>0)then
        ADIOS_WRITE(buf_id,mfmax)
        ADIOS_WRITE(buf_id,mf)
        ADIOS_WRITE(buf_id,rdtemf)
        ADIOS_WRITE(buf_id,zonalf)
        ADIOS_WRITE(buf_id,meshtf)
        ADIOS_WRITE(buf_id,meshnf)
        ADIOS_WRITE(buf_id,pfluxf)
        ADIOS_WRITE(buf_id,zfast)

        ADIOS_WRITE(buf_id,pmarkf)
        ADIOS_WRITE(buf_id,markerf)

     endif
      if(feload>0)then
        ADIOS_WRITE(buf_id,mfemax)
        ADIOS_WRITE(buf_id,mfe)
        ADIOS_WRITE(buf_id,rdtemfe)
        ADIOS_WRITE(buf_id,zonalfe)
        ADIOS_WRITE(buf_id,meshtfe)
        ADIOS_WRITE(buf_id,meshnfe)
        ADIOS_WRITE(buf_id,pfluxfe)
        ADIOS_WRITE(buf_id,zfaste)

        ADIOS_WRITE(buf_id,pmarkfe)
        ADIOS_WRITE(buf_id,markerfe)
  
        ADIOS_WRITE(buf_id,fetfrac)
        ADIOS_WRITE(buf_id,fetfracn)      
  
     endif
 
!     call adios_get_data_size (buf_id, mype_filesize)
     call adios_close (buf_id)
     if(mype==0)write(gtcout,*)'write to ',trim(restart_fname)

     call restart_hist
!!***********************************************

  else
!!     if(mype==0)write(gtcout,*)"param::",mpsi,mzeta
!     if(mype==0)write(gtcout,*)'read in',restart_fname
     allocate(zion0_read(mimax))
!!     ADIOS_READ(buf_id,mzeta)
     ADIOS_READ(buf_id,mi)
     ADIOS_READ(buf_id,me)
     ADIOS_READ(buf_id,rdtemi)
     ADIOS_READ(buf_id,pfluxi)
     ADIOS_READ(buf_id,phi00)
     ADIOS_READ(buf_id,phip00)
     ADIOS_READ(buf_id,zonali)
     ADIOS_READ(buf_id,meshti)
     ADIOS_READ(buf_id,meshte)
     ADIOS_READ(buf_id,meshni)
     ADIOS_READ(buf_id,meshne)
     call adios_read(buf_id,"zion0"//char(0),zion0_read)
     ADIOS_READ(buf_id,phi)
     ADIOS_READ(buf_id,zion)

     ADIOS_READ(buf_id,pmarki)
     ADIOS_READ(buf_id,markeri)
     ADIOS_READ(buf_id,markerit)

     if(nhybrid>0)then
        allocate(zelectron0_read(memax))
        call adios_read(buf_id,"zelectron0"//char(0),zelectron0_read)
        ADIOS_READ(buf_id,zelectron)
        ADIOS_READ(buf_id,phisave)
        ADIOS_READ(buf_id,rdteme)
        ADIOS_READ(buf_id,pfluxe)
        ADIOS_READ(buf_id,zonale)

        ADIOS_READ(buf_id,pmarke)
        ADIOS_READ(buf_id,markere)
        ADIOS_READ(buf_id,tfrac)
        ADIOS_READ(buf_id,tfracn)


     endif
     if(magnetic>0)then
        ADIOS_READ(buf_id,apara)
        ADIOS_READ(buf_id,fluidne)
     endif

     if(fload>0)then
        ADIOS_READ(buf_id,mfmax)
        ADIOS_READ(buf_id,mf)
        ADIOS_READ(buf_id,rdtemf)
        ADIOS_READ(buf_id,zonalf)
        ADIOS_READ(buf_id,meshtf)
        ADIOS_READ(buf_id,meshnf)
        ADIOS_READ(buf_id,pfluxf)
        ADIOS_READ(buf_id,zfast)

        ADIOS_READ(buf_id,pmarkf)
        ADIOS_READ(buf_id,markerf)

     endif
     
     if(feload>0)then
        ADIOS_READ(buf_id,mfemax)
        ADIOS_READ(buf_id,mfe)
        ADIOS_READ(buf_id,rdtemfe)
        ADIOS_READ(buf_id,zonalfe)
        ADIOS_READ(buf_id,meshtfe)
        ADIOS_READ(buf_id,meshnfe)
        ADIOS_READ(buf_id,pfluxfe)
        ADIOS_READ(buf_id,zfaste)
 
        ADIOS_READ(buf_id,pmarkfe)
        ADIOS_READ(buf_id,markerfe)
 
        ADIOS_READ(buf_id,fetfrac)
        ADIOS_READ(buf_id,fetfracn)
 
      endif

!     call adios_get_data_size (buf_id, mype_filesize)
     call adios_close (buf_id)
     if(mype==0)write(gtcout,*)'read in ',trim(restart_fname)
     irest=irest+1
!     if(myrank_partd<nproc_partd-1)call MPI_Wait(mrequest,mstatus,ierror)
!     if(mype==0)write(gtcout,*)restart_fname,'read over'
  endif  ! end of read

#else
  call restart_native(iop)
#endif

end subroutine restart_io
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine restart_hist
  use global_parameters
  implicit none
  integer i,ii(6),mstepfinal,ndstep,ndata
  real dum
  character(len=80) :: restart_fname

!! write out two copies of history.out and data1d.out for restart from a crash
     if(mype==0 .and. istep<=mstep)then
        if(mod(irest,2)==0)then
           restart_fname="restart_dir1/"
        else
           restart_fname="restart_dir2/"
        endif

! history.out
        open(333,file=trim(restart_fname)//'history_restart.out',status='replace')
        rewind(iodiag)
        read(iodiag,101)mstepfinal
        ndstep=mstepfinal-mstep/ndiag+istep/ndiag
        write(333,101)ndstep
        do i=1,5
           read(iodiag,101)ii(i)
           write(333,101)ii(i)
        enddo

        ndata=ii(1)*ii(2)+ii(3)*(2*ii(4)+ii(5))
        do i=0,ndata*ndstep
           read(iodiag,102)dum
           write(333,102)dum
        enddo
        close(333)
        write(gtcout,*)'write to ',trim(restart_fname)//'history_restart.out'
! data1d.out
        open(444,file=trim(restart_fname)//'data1d_restart.out',status='replace')
        rewind(iodata1d)
        read(iodata1d,101)mstepfinal
        ndstep=mstepfinal-mstep/ndiag+istep/ndiag
        write(444,101)ndstep
        do i=1,6
           read(iodata1d,101)ii(i)
           write(444,101)ii(i)
        enddo


        ndata=ii(1)*(ii(2)*ii(4)+ii(5)*ii(6))
        do i=1,ndata*ndstep
           read(iodata1d,102)dum
           write(444,102)dum
        enddo
        close(444)
        write(gtcout,*)'write to ',trim(restart_fname)//'data1d_restart.out'

! save current location of restart files
        open(555,file="FileExit.out",status="replace")
        write(555,"(A9,i1)")"FileExit=",FileExit
        write(555,"(A9,i5)")"irest   =",irest+1
        if(mod(irest,2)==0)then
           write(555,"(A12)")"restart_dir1"
        else
           write(555,"(A12)")"restart_dir2"
        endif
        close(555)

        if(istep==mstep)then
           close(iodiag)
           close(iodata1d)
        endif

     endif

     irest=irest+1 !!XY rolling rstart
101  format(i6)
102  format(e13.6)

end subroutine restart_hist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
