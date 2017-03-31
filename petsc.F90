! fix some PETSc versioning problems
#  define PETSC_DEFAULT_REAL PETSC_DEFAULT_DOUBLE_PRECISION

#ifdef _PETSc30ANDBEFORE
#  define _PETSc31P8ANDBEFORE
#endif

#ifdef _PETSc31P8ANDBEFORE
#  define _PETSc35BEFORE
#  define PetscBool PetscTruth
#  define MatCreateAIJ MatCreateMPIAIJ
#endif

#ifdef _PETSc35BEFORE
#  define KSPSetOperators(ksp,Apetsc,Bpetsc,ierr) KSPSetOperators(ksp,Apetsc,Bpetsc,DIFFERENT_NONZERO_PATTERN,ierr)
#endif

subroutine petsc_init
  use global_parameters,only:mype,partd_comm,gtcout,psi0
  use petsc_array,only: newcomm,nproc_newcomm,myrank_newcomm

  implicit none 

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  PetscErrorCode ierr
  integer ierror

#ifndef _PETSc35BEFORE
  if(psi0==0.0)then
    write(gtcout,*)'ERROR: Known bug for using psi0=0.0 with PETSc version 3.5 and later.'
    write(gtcout,*)'       Try loading a module for an older version of PETSc and setting'
    write(gtcout,*)'       the -D_PETSc35BEFORE flag as a PETSC_OPT in the Makefile.'
    write(gtcout,*)'       PETSc version 3.4.4.0 is known to work.'
    write(gtcout,*)'       Stopping the code...'
    stop
  endif
#endif

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_Comm_dup(partd_comm,newcomm,ierror)

  call mpi_comm_size(newcomm,nproc_newcomm,ierror)

  call mpi_comm_rank(newcomm,myrank_newcomm,ierror)

  if(ierror /= 0) then
    write(*,*)mype,'main_petsc: newcomm error',ierror
    call MPI_ABORT(MPI_COMM_WORLD,1,ierror)
  endif
        
  return
end subroutine petsc_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine petsc_final
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  Vec   x,b,x_recv,x2,b2,x_recv2
  Mat   Apetsc,Apetsc2
  KSP   ksp,ksp2
  PC    pc,pc2
  IS    is_from,is_to,is_from2,is_to2
  VecScatter     gather,gather2

  common/petsc_common/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather
  common/petsc_common2/Apetsc2,x2,b2,ksp2,pc2,x_recv2,is_from2,is_to2,gather2

  PetscErrorCode ierr

  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)
  call VecDestroy(x_recv,ierr)
  call VecScatterDestroy(gather,ierr)
  call ISDestroy(is_from,ierr)
  call ISDestroy(is_to,ierr)
  call MatDestroy(Apetsc,ierr)
  call KSPDestroy(ksp,ierr)

  call VecDestroy(x2,ierr)
  call VecDestroy(b2,ierr)
  call VecDestroy(x_recv2,ierr)
  call VecScatterDestroy(gather2,ierr)
  call ISDestroy(is_from2,ierr)
  call ISDestroy(is_to2,ierr)
  call MatDestroy(Apetsc2,ierr)
  call KSPDestroy(ksp2,ierr)

  call PetscFinalize(ierr)

  return
end subroutine petsc_final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine petsc_create(mindex)
  use global_parameters,only: mgrid
  use petsc_array, only:newcomm,usera,userp,users
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  integer mindex
  double precision tol,value(mindex)
  PetscErrorCode ierr
  PetscInt n,i,j,ij,ilow,ihigh,nz
  PetscInt idx(1)
  PetscBool flg
  integer col(mindex)

  n=mgrid

  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
  call VecCreate(newcomm,x,ierr)
  call VecSetSizes(x,PETSC_DECIDE,n,ierr)
  call VecSetFromOptions(x,ierr)
  call VecDuplicate(x,b,ierr)

  call VecCreate(PETSC_COMM_SELF,x_recv,ierr)
  call VecSetSizes(x_recv,n,n,ierr)
  call VecSetType(x_recv,VECSEQ,ierr)
  call VecSetFromOptions(x_recv,ierr)

  call ISCreateStride(newcomm,n,0,1,is_from,ierr)
  call ISCreateStride(PETSC_COMM_SELF,n,0,1,is_to,ierr)

  call VecScatterCreate(x,is_from,x_recv,is_to,gather,ierr)


!     call MatCreate(newcomm,Apetsc,ierr)
!     CHKERRQ(ierr)
!     call MatSetSizes(Apetsc,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!     CHKERRQ(ierr)
!     call MatSetType(Apetsc,MATMPIAIJ,ierr)
!     CHKERRQ(ierr)
!*optimized matrix creation
  call MatCreateAIJ(newcomm,PETSC_DECIDE,PETSC_DECIDE,n,n, &
    mindex,PETSC_NULL_INTEGER,mindex,PETSC_NULL_INTEGER,Apetsc,ierr)
  call MatSetFromOptions(Apetsc,ierr)
  call MatGetOwnershipRange(Apetsc,ilow,ihigh,ierr)

  ij=users(ilow-1)
  do i=ilow,ihigh-1
    nz=users(i)-users(i-1)
    do j=1,nz
      ij=ij+1
      col(j) = userp(ij)
      value(j) = usera(ij)
    enddo
    idx(1)=i
    call MatSetValues(Apetsc,1,idx,nz,col,value,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(Apetsc,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Apetsc,MAT_FINAL_ASSEMBLY,ierr)
  !     call MatTranspose(Apetsc,Bpetsc1,ierr)
  !     call MatAXPY(Bpetsc1,1.d0,Apetsc,DIFFERENT_NONZERO_PATTERN,ierr)
  !     call MatPartitioningSetAdjacency(mpart, Bpetsc1, ierr )
  !     call MatCopy(Apetsc,Bpetsc1,DIFFERENT_NONZERO_PATTERN,ierr)

  call KSPCreate(newcomm,ksp,ierr)
  call KSPSetOperators(ksp,Apetsc,Apetsc,ierr)
  !     call KSPGetPC(ksp,pc,ierr)
  !     call PCSetType(pc,PCJACOBI,ierr)
  ! 1.d-4 seems to be adequate
  !      tol=1.d-4
        tol=1.d-6
  !     tol=1.d-8
  call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL, &
    PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
  call KSPSetFromOptions(ksp,ierr)

  return
end subroutine petsc_create

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine petsc_solver
  use global_parameters,only:istep,irk,mype
  use petsc_array,only:userx,userb
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif
 
  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  PetscScalar temp
  PetscErrorCode ierr
  PetscInt i,ilow,ihigh,isize,ilo,ihi
!  PetscInt isi
!  PetscOffset idx
!  PetscScalar vec_x(1)
  PetscScalar, pointer :: vec_x(:)

  call VecGetLocalSize(b,isize,ierr)
  call VecGetOwnershipRange(b,ilow,ihigh,ierr)
  do i=ilow,ihigh-1
    temp=userb(i)
    call VecSetValue(b,i,temp,INSERT_VALUES,ierr)
  enddo
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

!     if(istep==1 .and. irk==1 .and. mype==0)then
!       write(*,*)'checking vector b:'
!       call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
!       CHKERRQ(ierr)
!     endif
!     if(istep==1 .and. irk==1 .and. mype==0)then
!       write(*,*)'checking matrix a:'
!       call MatView(Apetsc,PETSC_VIEWER_STDOUT_SELF,ierr)
!       CHKERRQ(ierr)
!     endif

  call KSPSetInitialGuessNonzero(ksp,Petsc_True,ierr)
  call KSPSolve(ksp,b,x,ierr)

!* 
  call VecScatterBegin(gather,x,x_recv,INSERT_VALUES, &
    SCATTER_FORWARD,ierr)
  call VecScatterEnd(gather,x,x_recv,INSERT_VALUES, &
    SCATTER_FORWARD,ierr)
!*
  call VecGetOwnershipRange(x_recv,ilo,ihi,ierr)
!  call VecGetLocalSize(x_recv,isi,ierr)
!  call VecGetArray(x_recv,vec_x,idx,ierr)
  call VecGetArrayF90(x_recv,vec_x,ierr)
  do i=ilo,ihi-1
    userx(i)=vec_x(i-ilo+1)
  enddo
!  call VecRestoreArray(x,vec_x,idx,ierr)
  call VecRestoreArrayF90(x_recv,vec_x,ierr)

  return
end subroutine petsc_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine laplacian_core
  use global_parameters,only:istep,irk,mype
  use petsc_array,only:userx,userb
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif
 
  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  PetscScalar temp
  PetscErrorCode ierr
  PetscInt i,ilow,ihigh,isize,ilo,ihi
!  PetscInt isi
!  PetscOffset idx
!  PetscScalar vec_x(1)
  PetscScalar, pointer :: vec_x(:)

  call VecGetLocalSize(b,isize,ierr)
  call VecGetOwnershipRange(b,ilow,ihigh,ierr)
  do i=ilow,ihigh-1
    temp=userb(i)
    call VecSetValue(b,i,temp,INSERT_VALUES,ierr)
  enddo
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

!     if(istep==1 .and. irk==1 .and. mype==0)then
!       write(*,*)'checking vector b:'
!       call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
!       CHKERRQ(ierr)
!     endif
!     if(istep==1 .and. irk==1 .and. mype==0)then
!       write(*,*)'checking matrix a:'
!       call MatView(Apetsc,PETSC_VIEWER_STDOUT_SELF,ierr)
!       CHKERRQ(ierr)
!     endif

  call MatMult(Apetsc,b,x,ierr)

!* 
  call VecScatterBegin(gather,x,x_recv,INSERT_VALUES, &
    SCATTER_FORWARD,ierr)
  call VecScatterEnd(gather,x,x_recv,INSERT_VALUES, &
    SCATTER_FORWARD,ierr)
!*
  call VecGetOwnershipRange(x_recv,ilo,ihi,ierr)
!  call VecGetLocalSize(x_recv,isi,ierr)
!  call VecGetArray(x_recv,vec_x,idx,ierr)
  call VecGetArrayF90(x_recv,vec_x,ierr)
  do i=ilo,ihi-1
    userx(i)=vec_x(i-ilo+1)
  enddo
!  call VecRestoreArray(x,vec_x,idx,ierr)
  call VecRestoreArrayF90(x_recv,vec_x,ierr)

  return
end subroutine laplacian_core

!********************************YX following

subroutine lapmat_pscreate(mindexlap,matrixsize)
  use global_parameters,only: mgridlow,mgridhigh,mype,istep,gtcout
  use petsc_array, only:newcomm,lusera,luserp,lusers
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common2/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  integer mindexlap,matrixsize
  double precision tol,value(mindexlap)
  PetscErrorCode ierr
  PetscInt n,i,j,ij,ilow,ihigh,nz
  PetscInt idx(1)
  PetscBool flg
  integer col(mindexlap)

  n=matrixsize

  call PetscOptionsGetInt(PETSC_NULL_CHARACTER,'-n',n,flg,ierr)
  call VecCreate(newcomm,x,ierr)
  call VecSetSizes(x,PETSC_DECIDE,n,ierr)
  call VecSetFromOptions(x,ierr)
  call VecDuplicate(x,b,ierr)

!*
  call VecCreate(PETSC_COMM_SELF,x_recv,ierr)
  call VecSetSizes(x_recv,n,n,ierr)
  call VecSetType(x_recv,VECSEQ,ierr)
  call VecSetFromOptions(x_recv,ierr)

  call ISCreateStride(newcomm,n,0,1,is_from,ierr)
  call ISCreateStride(PETSC_COMM_SELF,n,0,1,is_to,ierr)

  call VecScatterCreate(x,is_from,x_recv,is_to,gather,ierr)
!*

!     call MatCreate(newcomm,Apetsc,ierr)
!     CHKERRQ(ierr)
!     call MatSetSizes(Apetsc,PETSC_DECIDE,PETSC_DECIDE,n,n,ierr)
!     CHKERRQ(ierr)
!     call MatSetType(Apetsc,MATMPIAIJ,ierr)
!     CHKERRQ(ierr)
!*optimized matrix creation
  call MatCreateAIJ(newcomm,PETSC_DECIDE,PETSC_DECIDE,n,n, &
     mindexlap,PETSC_NULL_INTEGER,mindexlap,PETSC_NULL_INTEGER, &
     Apetsc,ierr)

  call MatSetFromOptions(Apetsc,ierr)
  call MatGetOwnershipRange(Apetsc,ilow,ihigh,ierr)

  ij=lusers(ilow-1)
  do i=ilow,ihigh-1
    nz=lusers(i)-lusers(i-1)
    do j=1,nz
      ij=ij+1
      col(j) = luserp(ij)
      value(j) = lusera(ij)
    enddo
    idx(1)=i
    call MatSetValues(Apetsc,1,idx,nz,col,value,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(Apetsc,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Apetsc,MAT_FINAL_ASSEMBLY,ierr)
!     call MatTranspose(Apetsc,Bpetsc1,ierr)
!     call MatAXPY(Bpetsc1,1.d0,Apetsc,DIFFERENT_NONZERO_PATTERN,ierr)
!     call MatPartitioningSetAdjacency(mpart, Bpetsc1, ierr )
!     call MatCopy(Apetsc,Bpetsc1,DIFFERENT_NONZERO_PATTERN,ierr)

  call KSPCreate(newcomm,ksp,ierr)
  call KSPSetOperators(ksp,Apetsc,Apetsc,ierr)
!     call KSPGetPC(ksp,pc,ierr)
!     call PCSetType(pc,PCJACOBI,ierr)
! 1.d-4 seems to be adequate
!  tol=1.d-4
  tol=1.d-6
!     tol=1.d-8
  call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL, &
       PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
  call KSPSetFromOptions(ksp,ierr)

  if(mype==0 .and. istep==1)then
    write(gtcout,*)'Finish create lap-petsc matrix'
    call flush(gtcout)
  endif


  return
end subroutine lapmat_pscreate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lapmat_pssolver
  use global_parameters,only:istep,irk,mype
  use petsc_array,only:luserx,luserb
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif
 
  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common2/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  PetscScalar temp
  PetscErrorCode ierr
  PetscInt i,ilow,ihigh,isize,ilo,ihi
!  PetscInt isi
!  PetscOffset idx
!  PetscScalar vec_x(1)
  PetscScalar, pointer :: vec_x(:)

  call VecGetLocalSize(b,isize,ierr)
  call VecGetOwnershipRange(b,ilow,ihigh,ierr)
  do i=ilow,ihigh-1
    temp=luserb(i)
    call VecSetValue(b,i,temp,INSERT_VALUES,ierr)
  enddo
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

!!!!!!!!!!!!!!!!!!!diag petsc  !!!!!!!
!!!!     if(istep==1 .and. mype==0)then 
!!!   write(*,*)'checking vector b:'
!   call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
!   CHKERRQ(ierr)
!!!     endif
!!!     if(istep==1 .and. mype==0)then
!!!   write(*,*)'checking matrix a:'
!!!   call MatView(Apetsc,PETSC_VIEWER_STDOUT_SELF,ierr)
!!!   CHKERRQ(ierr)
!!     endif

  call KSPSetInitialGuessNonzero(ksp,Petsc_True,ierr)
  call KSPSolve(ksp,b,x,ierr)

!* 
  call VecScatterBegin(gather,x,x_recv,INSERT_VALUES, &
       SCATTER_FORWARD,ierr)
  call VecScatterEnd(gather,x,x_recv,INSERT_VALUES, &
       SCATTER_FORWARD,ierr)
!*
  call VecGetOwnershipRange(x_recv,ilo,ihi,ierr)
!  call VecGetLocalSize(x_recv,isi,ierr)
!  call VecGetArray(x_recv,vec_x,idx,ierr)
  call VecGetArrayF90(x_recv,vec_x,ierr)
  do i=ilo,ihi-1
    luserx(i)=vec_x(i-ilo+1)
  enddo
!  call VecRestoreArray(x,vec_x,idx,ierr)
  call VecRestoreArrayF90(x_recv,vec_x,ierr)

  return
end subroutine lapmat_pssolver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lapmat_psmult
  use global_parameters,only:istep,irk,mype,gtcout
  use petsc_array,only:luserx,luserb
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif
 
  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common2/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  PetscScalar temp
  PetscErrorCode ierr
  PetscInt i,ilow,ihigh,isize,ilo,ihi
!  PetscInt isi
!  PetscOffset idx
!  PetscScalar vec_x(1)
  PetscScalar, pointer :: vec_x(:)

  call VecGetLocalSize(b,isize,ierr)
  call VecGetOwnershipRange(b,ilow,ihigh,ierr)
  do i=ilow,ihigh-1
    temp=luserb(i)
    call VecSetValue(b,i,temp,INSERT_VALUES,ierr)
  enddo
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)

!     if(istep==1 .and. irk==1 .and. mype==0)then
!   write(*,*)'checking vector b:'
!   call VecView(b,PETSC_VIEWER_STDOUT_SELF,ierr)
!   CHKERRQ(ierr)
!     endif
!     if(istep==1 .and. irk==1 .and. mype==0)then
!   write(*,*)'checking matrix a:'
!   call MatView(Apetsc,PETSC_VIEWER_STDOUT_SELF,ierr)
!   CHKERRQ(ierr)
!     endif

  call MatMult(Apetsc,b,x,ierr)

!* 
  call VecScatterBegin(gather,x,x_recv,INSERT_VALUES, &
       SCATTER_FORWARD,ierr)
  call VecScatterEnd(gather,x,x_recv,INSERT_VALUES, &
       SCATTER_FORWARD,ierr)
!*
  call VecGetOwnershipRange(x_recv,ilo,ihi,ierr)
!  call VecGetLocalSize(x_recv,isi,ierr)
!  call VecGetArray(x_recv,vec_x,idx,ierr)
  call VecGetArrayF90(x_recv,vec_x,ierr)
  do i=ilo,ihi-1
    luserx(i)=vec_x(i-ilo+1)
  enddo
!  call VecRestoreArray(x,vec_x,idx,ierr)
  call VecRestoreArrayF90(x_recv,vec_x,ierr)
  if(mype==0 .and. istep==0)then
    write(gtcout,*)'Finish laplacian multiplication'
    call flush(gtcout)
  endif

  return

end subroutine lapmat_psmult !!lap_psoperator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lapmat_pscreate2(mindexlap,matrixsize)
  use global_parameters,only: mype,istep,gtcout
  use petsc_array, only:newcomm,lusera2,luserp2,lusers2
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common3/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  integer mindexlap,matrixsize
  double precision tol,value(mindexlap)
  PetscErrorCode ierr
  PetscInt n,i,j,ij,ilow,ihigh,nz
  PetscInt idx(1)
  PetscBool flg
  integer col(mindexlap)

  n=matrixsize

  call VecCreate(newcomm,x,ierr)
  call VecSetSizes(x,PETSC_DECIDE,n,ierr)
  call VecSetFromOptions(x,ierr)
  call VecDuplicate(x,b,ierr)

!*
  call VecCreate(PETSC_COMM_SELF,x_recv,ierr)
  call VecSetSizes(x_recv,n,n,ierr)
  call VecSetType(x_recv,VECSEQ,ierr)
  call VecSetFromOptions(x_recv,ierr)

  call ISCreateStride(newcomm,n,0,1,is_from,ierr)
  call ISCreateStride(PETSC_COMM_SELF,n,0,1,is_to,ierr)

  call VecScatterCreate(x,is_from,x_recv,is_to,gather,ierr)

!*optimized matrix creation
  call MatCreateAIJ(newcomm,PETSC_DECIDE,PETSC_DECIDE,n,n, &
     mindexlap,PETSC_NULL_INTEGER,mindexlap,PETSC_NULL_INTEGER, &
     Apetsc,ierr)

  call MatSetFromOptions(Apetsc,ierr)
  call MatGetOwnershipRange(Apetsc,ilow,ihigh,ierr)

  ij=lusers2(ilow-1)
  do i=ilow,ihigh-1
    nz=lusers2(i)-lusers2(i-1)
    do j=1,nz
      ij=ij+1
      col(j) = luserp2(ij)
      value(j) = lusera2(ij)
    enddo
    idx(1)=i
    call MatSetValues(Apetsc,1,idx,nz,col,value,INSERT_VALUES,ierr)
  enddo

  call MatAssemblyBegin(Apetsc,MAT_FINAL_ASSEMBLY,ierr)
  call MatAssemblyEnd(Apetsc,MAT_FINAL_ASSEMBLY,ierr)

  call KSPCreate(newcomm,ksp,ierr)
  call KSPSetOperators(ksp,Apetsc,Apetsc,ierr)

! 1.d-4 seems to be adequate
!  tol=1.d-4
  tol=1.d-6
!     tol=1.d-8
  call KSPSetTolerances(ksp,tol,PETSC_DEFAULT_REAL, &
       PETSC_DEFAULT_REAL,PETSC_DEFAULT_INTEGER,ierr)
  call KSPSetFromOptions(ksp,ierr)
  return
end subroutine lapmat_pscreate2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lapmat_pssolver2
  use global_parameters,only:istep,irk,mype
  use petsc_array,only:luserx2,luserb2
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common3/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather

  PetscScalar temp
  PetscErrorCode ierr
  PetscInt i,ilow,ihigh,isize,ilo,ihi
  PetscScalar, pointer :: vec_x(:)

  call VecGetLocalSize(b,isize,ierr)
  call VecGetOwnershipRange(b,ilow,ihigh,ierr)
  do i=ilow,ihigh-1
    temp=luserb2(i)
    call VecSetValue(b,i,temp,INSERT_VALUES,ierr)
  enddo
  call VecAssemblyBegin(b,ierr)
  call VecAssemblyEnd(b,ierr)
!  call KSPSetInitialGuessNonzero(ksp,Petsc_True,ierr)
  call KSPSolve(ksp,b,x,ierr)

!*
  call VecScatterBegin(gather,x,x_recv,INSERT_VALUES, &
       SCATTER_FORWARD,ierr)
  call VecScatterEnd(gather,x,x_recv,INSERT_VALUES, &
       SCATTER_FORWARD,ierr)
!*
  call VecGetOwnershipRange(x_recv,ilo,ihi,ierr)
  call VecGetArrayF90(x_recv,vec_x,ierr)
  do i=ilo,ihi-1
    luserx2(i)=vec_x(i-ilo+1)
  enddo
  call VecRestoreArrayF90(x_recv,vec_x,ierr)

  return
end subroutine lapmat_pssolver2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine petsc2_final
  implicit none

#ifdef _PETSc30ANDBEFORE
#include "finclude/petscall.h90"
#else
#include "petsc/finclude/petsc.h90"
#endif

  Vec   x,b,x_recv
  Mat   Apetsc
  KSP   ksp
  PC    pc
  IS    is_from,is_to
  VecScatter     gather

  common/petsc_common3/Apetsc,x,b,ksp,pc,x_recv,is_from,is_to,gather
  

  PetscErrorCode ierr

  call VecDestroy(x,ierr)
  call VecDestroy(b,ierr)
  call VecDestroy(x_recv,ierr)
  call VecScatterDestroy(gather,ierr)
  call ISDestroy(is_from,ierr)
  call ISDestroy(is_to,ierr)
  call MatDestroy(Apetsc,ierr)
  call KSPDestroy(ksp,ierr)

  call PetscFinalize(ierr)

  return
end subroutine petsc2_final
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
