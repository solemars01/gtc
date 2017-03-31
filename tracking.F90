subroutine tag
  use global_parameters
  use particle_array
  use particle_tracking
  implicit none

! Each MPI process has its own "ptracked" array to accumulate the tracked
! particles that may or may not reside in its subdomain.
! The vector "ntrackp" contains the number of tracked particles currently residing on a MPI process

  ntrackp=0
  allocate(ptrackedi(nparami,max(mimax,1))) !!array to store ion tracked particles
  ptrackedi=0.0
  if(nhybrid>0)then
     allocate(ptrackede(nparame,max(memax,1)))  !!array to store tracked electrons
     ptrackede=0.0
  endif

! If irun==0 and if particle tracking is "on", tag each particle with a unique number
  call tag_particles

end subroutine tag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine partout

  use global_parameters
  implicit none
  call locate_tracked_particles
  call write_tracked_particles

end subroutine partout

!========================================
!! originally written by s. either, modified by Y. Xiao April 2008
subroutine tag_particles

!========================================

  use global_parameters
  use particle_array
  use field_array
  use particle_tracking
  implicit none
  integer :: m,np
  real(lk) :: ainside,aoutside,thickratio,a0,a1
! We keep track of the particles by tagging them with a unique number.
! We add an extra element to the particle array, which holds the
! particle tag, i.e. just a number.
! The input parameter "nptrack" is the total number of particles that
! we track.

! We tag each particle with a unique number that will be carried along
! with it as it moves from one processor/domain to another. To facilitate
! the search for the tracked particles, we give the non-tracked particles
! a negative value.
! CAVEAT: We are storing the tag in a floating point array, which means that
!         for a large number of particles, the truncation error associated
!         with the precision may make some particles have the same tag.
!         This is particularly true in single precision where the maximum
!         number of particles that can be distinguished with tags differing
!         by unity is 2**24-1 = 16,777,215 (16.7 million). The most important
!         is to distinguish the tracked particles between them. The other
!         particles just need to have a negative tag. In order to minimize
!         the effect of the precision truncation, we retag the tracked
!         particles with positive numbers starting at 1.

  do m=1,mi
     zion(nparami-1,m)=-real(m,lk) !!!-real(m+mype*mi)
     zion(nparami,m)=real(mype+1,lk)
     zion0(nparami-1,m)=-real(m,lk)
     zion0(nparami,m)=real(mype+1,lk)
  enddo    

  if (nhybrid>0) then
     do m=1,me
        zelectron(nparame-1,m)=-real(m,lk) !!!-real(m+mype*mi)
        zelectron(nparame,m)=real(mype+1,lk)
        zelectron0(nparame-1,m)=-real(m,lk)
        zelectron0(nparame,m)=real(mype+1,lk)
     enddo

  endif

  a0=rg0
  a1=rg0+deltar*real(mpsi)
  thickratio=0.000002
  ainside=0.5*(1.0-thickratio)*a1+0.5*(1.0+thickratio)*a0
  ainside=ainside*ainside*0.5

  aoutside=0.5*(1.0+thickratio)*a1+0.5*(1.0-thickratio)*a0
  aoutside=aoutside*aoutside*0.5

  if (nhybrid>0) then
      np=0
      do m=1,me
        
         if(zelectron(1,m)>ainside .and. zelectron(1,m)<aoutside) then
            np=np+1
            zelectron(nparame-1,m)=real(np,lk)
            zelectron0(nparame-1,m)=real(np,lk)
         endif
      enddo
  endif

  if(mype==0)write(*,*)'np=',np,'ainside=',ainside,'aoutside=',aoutside
!!  if(mype==0)write(*,*)'zelectron(1,1:20)=',zelectron(1,1:20),'zelectron(nparame,1)=',zelectron(nparame,1)
!! we only condiser those particles originally satify a certain condition
! namely within an annulous ain<r<aout
! We divide the number of tracked particles equally between all processors
! as much as possible. If nptrack is not a multiple of numberpe, we add
! 1 particle to some of the the processors until we get to nptrack. We also
! start at mi/2 so that the lower numbers fall around r=0.5 (see subroutine 
! load for the initial particle distribution in r).
!  acentre=0.5*(a0+a1)
!  awidth=0.1*(a1-a0)*0.5
!  ainside= acentre-awidth
!  aoutside=acentre+awidth

 np=0
! this part is for tgaaing a single particle 
 ! if(mype==0)then
	! np=np+1
	! zion(nparami-1,1)=real(np,wp)
        ! zion0(nparami-1,1)=real(np,wp)
 ! endif

! use this part for selecting a group of particles

  do m=1,mi
!      r=sqrt(2.0*zion(1,m))
     if(zion(1,m)>ainside .and. zion(1,m)<aoutside) then
       np=np+1
       zion(nparami-1,m)=real(np,lk)
       zion0(nparami-1,m)=real(np,lk)
     endif
  enddo
!! if (nhybrid>0) then
!!      np=0
!!      do m=1,me
!!         r=sqrt(2.0*zelectron(1,m))
!!         if(r>ainside .and. r<aoutside) then
!!            np=np+1
!!            zelectron(nparami-1,m)=real(np,wp)
!!            zelectron0(nparami-1,m)=real(np,wp)
!!         endif
!!      enddo
!! endif

!  if(mype==0)write(78,*)'np=',np,'  npp=',npp,'  zion(nparami-1,mi/2)=',zion(nparami-1,mi/2),&
!                        '  zion(nparami-1,1)=',zion(nparami-1,1)

!  if(mype==0)then
!  ! On the master process (mype=0), we pick "nptrack" particles that will
!  ! be followed at every time step. To facilitate the search of those
!  ! particles among all the others in the zion arrays of each processor,
!  ! we give them a positive number from 1 to nptrack. The particles are
!  ! picked around r=0.5, which is mi/2.
!    do m=(mi-nptrack)/2,(mi+nptrack)/2-1
!       zion(nparami-1,m)=-zion(nparami-1,m)
!       zion0(nparami-1,m)=-zion(nparami-1,m)
!    enddo
!  endif

end subroutine tag_particles

!========================================

subroutine locate_tracked_particles

!========================================

  use global_parameters
  use particle_array
  use particle_tracking
  implicit none
  integer :: m,npp

! Check if tracked particles are located on this processor. Particles that
! are tracked at every time step have a positive zion(nparami-1,m) value. All the
! others have a negative value.
  ntrackp=0
  ptrackedi=0.0

  npp=0

  do m=1,mi
     if(zion(nparami-1,m)>0.0)then
       npp=npp+1
       ptrackedi(1:nparami,npp)=zion(1:nparami,m)
     endif
  enddo
  ntrackp(1)=npp
  npp=0

  if (nhybrid>0) then
     ptrackede=0.0
     do m=1,me
        if(zelectron(nparame-1,m)>0.0)then
          npp=npp+1
          ptrackede(1:nparame,npp)=zelectron(1:nparame,m)
        endif
     enddo
  endif
  ntrackp(2)=npp
!!  if(mype==0)write(*,*)'me=',me,'np=',npp,'zelectron(1,1:6)=',zelectron(nparame-1,1:6),&
!!                        'zelectron(nparame,1)= ',zelectron(nparame,1)

end subroutine locate_tracked_particles

!========================================

subroutine write_tracked_particles

!========================================

  use global_parameters
  use particle_tracking
  implicit none
 
  integer :: j
  character(len=30) :: cdum
#if ADIOS
    character(len=50),SAVE:: fname
    character(len=50)::dirstr
    integer*8 :: group_id, buf_id,comm
    #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
#endif


#if ADIOS
    comm=MPI_COMM_WORLD
    write(fname,'("trackp_dir/TRACKP_",i5.5,".bp")')mstepall+istep
    !!!if(mype==0)write(*,*)"tracking filename ",fname
    fname=trim(fname)//char(0)
    write(dirstr,'("node_",i5.5)')mype!,mstepall+istep !!!(mstepall+istep)/ndiag
    dirstr=trim(dirstr)//char(0)
    call adios_get_group (group_id, "particles"//char(0))
    call adios_set_path (group_id,dirstr//char(0))
    call adios_open (buf_id, group_id, fname)!!call adios_open_append(buf_id, group_id, fname)
    ADIOS_WRITE(buf_id,comm)
    ADIOS_WRITE(buf_id,mype)
    ADIOS_WRITE(buf_id,nparami)
    ADIOS_WRITE(buf_id,nspecies)
!!    call adios_write(buf_id,'ntracki'//char(0),ntrackp(1))
!!    call adios_write(buf_id,'ptrackedi'//char(0),ptrackedi(:,1:ntrackp(1)))
    if(nhybrid>0)then
        call adios_write(buf_id,'ntracke'//char(0),ntrackp(2))
        call adios_write(buf_id,'ptrackede'//char(0),ptrackede(:,1:ntrackp(2)))
    endif

    !ADIOS_WRITE(buf_id,ptracked)
    !call adios_get_data_size (buf_id, mype_filesize)
    !write(*,*)"npp,nparam,data_size",npp,nparam,filesize
    call adios_close (buf_id)
    if(mype==0)write(*,*)"tracking filename ",fname
#else
  
     write(cdum,'("trackp_dir/TRACKP.",i5.5)')mype
     cdum=trim(cdum)//char(0)
     if(irun==0 .and. istep==ndiag)then
        open(57,file=cdum,status='replace',position='append')
     else
        open(57,file=cdum,status='old',position='append')
     endif
     write(57,*)istep+irun*mstep
     write(57,*)ntrackp(1:nspecies)

     do j=1,ntrackp(1)
        write(57,*)ptrackedi(1:nparami,j)
     enddo
     if(nhybrid>0)then
        do j=1,ntrackp(2)
           write(57,*)ptrackede(1:nparame,j)
        enddo
     endif

     close(57)
#endif
end subroutine write_tracked_particles
!========================================
