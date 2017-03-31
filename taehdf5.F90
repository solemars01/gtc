!----------------------------------------------------------------------
! taehdf5.f                 -  routines to save data in HDF5
!
! Test examples: see test_hdf5_dumpread.f90 and test_hdf5_wr.f90 for usage.
!                See also dump_h5 routine in LamyRidge.
!
! This module consists of simplifying wrappers to HDF5 I/O routines:
!
! 1) h5dump_attr - write a scalar attribute    of type INTEGER, DOUBLE, or CHARACTER(len=*)
! 2) h5read_attr - read  a scalar attribute    of type INTEGER, DOUBLE, or CHARACTER(len=*)
! 3) h5dump_data - write a 0D, 1D, or 2D array of type INTEGER, DOUBLE, or CHARACTER(len=*)
! 4) h5read_data - read  a 0D, 1D, or 2D array of type INTEGER, DOUBLE, or CHARACTER(len=*)
!
! Also there are the WR wrappers to the above routines, which allow one to select between
! WRITE/READ using WR=0 or 1:
!
! 5) h5wr_attr  - write (WR=0) or read (WR=1) an attribute, i.e. call h5dump_attr or h5read_attr
! 6) h5wr_data  - write (WR=0) or read (WR=1) an array,     i.e. call h5dump_data or h5read_data
! 7) h5wr_file  - Create a new file  (WR=0) or Open an existing file  (WR=1)
! 8) h5wr_group - Create a new group (WR=0) or Open an existing group (WR=1)
!
! Note on scalars: I usually store scalars as attributes  (*_attr wrappers above) but you can 
!                  just as well store them in a dataspace (*_data wrappers above).
!
! Note on files: h5wr_file(0,...) will create a file as new, DELETING any existing contents.
!                h5wr_file(1,...) will open a new file in WRITE/READ mode - you can add new
!                                 data to an old file.
!
! Note on IDL: the above data can be quickly accessed with the idl/gui/get_h5_*pro routines.
!
! Also if you have some 2D data associated with a 2D mesh which you would like to display
! in the IDL GUI, then there is a convenient wrapper for that too (OUTPUT only at present):
!
! 9) printh5 - write structured HDF5 output - see LamyRidge/src/plot_hdf5.f for example.
!
! In the PRINTH5 interface, the first three arguments are REQUIRED:
!
!    FILE - HDF5 data file name
!    NAME - name of data group created in the file
!    F    - 1D or 2D data array created in the group
!
! The other arguments are OPTIONAL, and may be specified by KEYWORD in arbitrary order:
!
!    X      - 1D array of X-axis points (should have SIZE(F,DIM=1))
!    Y      - 1D array of Y-axis points (should have SIZE(F,DIM=2))
!    TITLE  - Title to appear on plot in IDL GUI
!    XTITLE - X-axis Title to appear on plot in IDL GUI
!    YTITLE - Y-axis Title to appear on plot in IDL GUI
!
! E.g. with ordered arguments:
!
!    call printh5('data_123.h5','Er', er)
!    call printh5('data_123.h5','Ez', ez)
!    call printh5('data_123.h5','T', tm ,r,z,'Temperature (J)','R (m)','Z (m)')
!
! or in arbitrary order by using keyword arguments:
!
!    call printh5('data_123.h5','Psi',psi,x=radius,xtitle='R (m)',y=z,ytitle='Z (m)')
!    call printh5('data_123.h5','Psi slice',psi(:,nz/2),title='Slice through midplane',x=radius,xtitle='R' )
!
! The matching IDL GUI code in taepic/idl will substitute values 
! for any OPTIONAL arguments that were omitted.
!
!
! See also            http://www.hdfgroup.org/HDF5/doc1.6/
! and in particular,  http://www.hdfgroup.org/HDF5/doc1.6/RM_H5Front.html
!----------------------------------------------------------------------
! Sean Dettrick, Oct '08.
!
! Questions: sean@trialphaenergy.com
!----------------------------------------------------------------------

module taehdf5
  use hdf5
  !use precision,only:wp
!  use numeric_kinds
  implicit none

  integer*1,parameter,private:: rn=8,cn=8 !!,longlong = selected_int_kind( 16 )

  ! store between printh5 statements:

  CHARACTER(len=200) :: last_file ! Previous filename
  
  logical,public:: taehdf5_quiet=.false.

  ! Access statements:

  public::printh5
  private::data_to_group

  ! Generic interfaces:

  interface printh5
     module procedure printh5_3D_group
     module procedure printh5_3D_file
     module procedure printh5_2D_group
     module procedure printh5_2D_file
     module procedure printh5_1D_group
     module procedure printh5_1D_file
  end interface

  interface get_h5_type
     module procedure get_h5_type_path
     module procedure get_h5_type_group
  end interface
  interface get_h5_dims
     module procedure get_h5_dims_group
     module procedure get_h5_dims_path
  end interface
  interface data_to_group
     module procedure data_to_group_2D
     module procedure data_to_group_1D
     module procedure data_to_group_1Dint
     module procedure data_to_group_string
  end interface
  interface h5wr_attr
     module procedure h5wr_attr_int
     module procedure h5wr_attr_double
     module procedure h5wr_attr_string
  end interface
  interface h5dump_attr
     module procedure h5dump_attr_int
     module procedure h5dump_attr_double
     module procedure h5dump_attr_string
  end interface
  interface h5read_attr
     module procedure h5read_attr_int
     module procedure h5read_attr_double
     module procedure h5read_attr_string
  end interface
  interface h5wr_data 
     module procedure h5wr_data_logical_0d
     module procedure h5wr_data_logical_1d
     module procedure h5wr_data_logical_2d
     module procedure h5wr_data_int_0d
     module procedure h5wr_data_int_1d
     module procedure h5wr_data_int_2d
     module procedure h5wr_data_double_0d
     module procedure h5wr_data_double_1d
     module procedure h5wr_data_double_2d
     module procedure h5wr_data_double_3d
     module procedure h5wr_data_string_0d
     module procedure h5wr_data_string_1d
     module procedure h5wr_data_string_2d
  end interface
  interface h5dump_data 
     module procedure h5dump_data_logical_0d
     module procedure h5dump_data_int_0d
     module procedure h5dump_data_int_1d
     module procedure h5dump_data_int_2d
!!$     module procedure h5dump_data_longlong_1d
     module procedure h5dump_data_double_0d
     module procedure h5dump_data_double_1d
     module procedure h5dump_data_double_2d
     module procedure h5dump_data_double_3d
     module procedure h5dump_data_string_0d_auto
     module procedure h5dump_data_string_0d
     module procedure h5dump_data_string_1d
     module procedure h5dump_data_string_2d
     module procedure h5dump_data_double_0d_auto
     module procedure h5dump_data_double_1d_auto
     module procedure h5dump_data_double_2d_auto
     module procedure h5dump_data_double_3d_auto
     module procedure h5dump_data_int_0d_auto
     module procedure h5dump_data_int_1d_auto
     module procedure h5dump_data_int_2d_auto
!!$     module procedure h5dump_data_longlong_1d_auto
  end interface
  interface h5read_data 
     module procedure h5read_data_logical_0d
     module procedure h5read_data_int_0d
     module procedure h5read_data_int_1d
     module procedure h5read_data_int_2d
     module procedure h5read_data_double_0d
     module procedure h5read_data_double_1d
     module procedure h5read_data_double_2d
     module procedure h5read_data_double_3d
     module procedure h5read_data_double_4d
     module procedure h5read_data_string_0d
     module procedure h5read_data_string_1d
     module procedure h5read_data_string_2d
     module procedure h5read_data_double_0d_auto
     module procedure h5read_data_double_1d_auto
     module procedure h5read_data_double_2d_auto
     module procedure h5read_data_double_3d_auto
     module procedure h5read_data_double_4d_auto
     module procedure h5read_data_int_0d_auto
     module procedure h5read_data_int_1d_auto
     module procedure h5read_data_int_2d_auto
  end interface
  interface h5append_data
     module procedure h5append_data_int_0d_fixed
     module procedure h5append_data_int_0d_auto_fixed
     module procedure h5append_data_double_0d_fixed
     module procedure h5append_data_double_0d_auto_fixed
     module procedure h5append_data_double_1d_fixed
     module procedure h5append_data_double_1d_auto_fixed
     module procedure h5append_data_double_0d_auto
     module procedure h5append_data_int_0d_auto
     module procedure h5append_data_double_0d
     module procedure h5append_data_int_0d
!!$     module procedure h5append_data_double_1d_auto
!!$     module procedure h5append_data_int_1d_auto
     module procedure h5append_data_double_1d
     module procedure h5append_data_int_1d
  end interface

  interface h5dataset_exists
     module procedure h5dataset_exists_id
     module procedure h5dataset_exists_path
  end interface
  interface h5group_exists
     module procedure h5group_exists_id
     module procedure h5group_exists_path
  end interface

contains

  !-----------------------------------------------------------------------
  ! Attempt to get HDF5 errors to abort the program, so it can be debugged
  !
  ! Sadly, an error function can't be set in Fortran, so we try to do
  ! it instead using the C interface.
  !
  ! See
  !    http://www.hdfgroup.org/HDF5/doc1.6/RM_H5E.html#Error-SetAuto
  !
  subroutine set_h5errors_abort
  end subroutine set_h5errors_abort

  !-----------------------------------------------------------------------
  ! Get the dimensions of a dataset by reading from its dataspace
  !
  ! See
  !    http://www.hdfgroup.org/HDF5/doc1.6/RM_H5D.html   for dataset API
  !    http://www.hdfgroup.org/HDF5/doc1.6/RM_H5S.html   for dataspace API
  !
  subroutine get_h5_dims_group(group_id,name,dims,num_dims)
    ! Arguments
    INTEGER(HID_T),               intent(IN) :: group_id ! Group identifier
    CHARACTER(len=*),             intent(IN) :: name     ! Name of dataspace
    INTEGER(HSIZE_T),DIMENSION(*),intent(OUT):: dims     ! Array of dims
    INTEGER                      ,intent(OUT):: num_dims ! Number of dims
    ! Local
    INTEGER                        :: ierr            ! Error code 
    INTEGER(HID_T)                 :: dataspace_id    ! Dataspace identifier 
    INTEGER(HID_T)                 :: dataset_id      ! Data set identifier 
    INTEGER(HSIZE_T), DIMENSION(4) :: maxdims         ! Array to store max dimension sizes 

    call h5dopen_f(group_id,name,dataset_id,ierr)        ! open dataset of NAME
    call h5dget_space_f(dataset_id, dataspace_id, ierr)  ! get dataspace
    call h5dclose_f(dataset_id,ierr)                   ! close dataset

    call h5sget_simple_extent_ndims_f(dataspace_id,num_dims,ierr)       ! get rank of dataset
    call h5sget_simple_extent_dims_f(dataspace_id, dims, maxdims, ierr) ! get dims of dataset
    return
  end subroutine get_h5_dims_group

  !---------------
  ! path interface

  subroutine get_h5_dims_path(file,path,dims,num_dims)
    ! Arguments
    CHARACTER(len=*),             intent(IN) :: file,path ! filename and internal HDF5 path input
    INTEGER(HSIZE_T),DIMENSION(*),intent(OUT):: dims      ! Array of dims
    INTEGER                      ,intent(OUT):: num_dims  ! Number of dims
    ! Local
    !INTEGER                        :: ierr            ! Error code 
    !INTEGER(HID_T)                 :: dataspace_id    ! Dataspace identifier 
    !INTEGER(HID_T)                 :: dataset_id      ! Data set identifier 
    !INTEGER(HSIZE_T), DIMENSION(4) :: maxdims         ! Array to store max dimension sizes 
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0 .AND. h5dataset_exists(ids(last),name)) then
       call get_h5_dims(ids(last),name,dims,num_dims)
    else
       if (.not.taehdf5_quiet) print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
       num_dims=0
    endif
    call h5close_path(ids)

    return
  end subroutine get_h5_dims_path

  !-----------------------------------------------------------------------
  ! Get the type of a dataset by reading from its dataspace
  !
  !    call get_h5_type('dump000.bin','mesh/v_kind',type_id)

  subroutine get_h5_type_path(file,path,type_id)
    ! Arguments
    CHARACTER(len=*),intent(IN) :: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),  intent(OUT):: type_id   ! matching Native data type
    ! Local
    INTEGER        :: ierr            ! Error code 
    INTEGER(HID_T) :: etype_id        ! Existing data type
    !INTEGER(HID_T) :: dataspace_id    ! Dataspace identifier 
    !INTEGER(HID_T) :: dataset_id      ! Data set identifier 
    !logical::equal
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    call h5read_path(file,path,ids,last,name)
    if (ids(last)<=0 .OR. .NOT.h5dataset_exists(ids(last),name)) then
       if (.not.taehdf5_quiet) print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
       type_id=0
       call h5close_path(ids)
       return
    endif

    call get_h5_type_group(ids(last),TRIM(name),type_id)

    call h5close_path(ids)

    return
  end subroutine get_h5_type_path

  !-----------------------------------------------------------------------
  ! Get the type of a dataset by reading from its dataspace
  !
  !    call get_h5_type(mesh_id, 'v_kind',type_id)

  subroutine get_h5_type_group(loc_id,name,type_id)
    ! Arguments
    INTEGER(HID_T),  intent(IN):: loc_id     ! Location identifier
    CHARACTER(len=*),intent(IN):: name       ! Name of object
    INTEGER(HID_T),  intent(OUT):: type_id   ! matching Native data type
    ! Local
    INTEGER        :: ierr            ! Error code 
    INTEGER(HID_T) :: etype_id        ! Existing data type
    !INTEGER(HID_T) :: dataspace_id    ! Dataspace identifier 
    INTEGER(HID_T) :: dataset_id      ! Data set identifier 
    logical::equal
    !INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    !integer :: last

    call h5dopen_f(loc_id,name,dataset_id,ierr)
    call h5dget_type_f(dataset_id, etype_id, ierr)  ! get existing datatype
    call h5dclose_f(dataset_id,ierr)

    call h5tequal_f(etype_id,H5T_NATIVE_DOUBLE,equal,ierr)
    if (equal) then
       call h5tclose_f(etype_id,ierr)
       type_id=H5T_NATIVE_DOUBLE
       return
    endif

    call h5tequal_f(etype_id,H5T_NATIVE_INTEGER,equal,ierr)
    if (equal) then
       call h5tclose_f(etype_id,ierr)
       type_id=H5T_NATIVE_INTEGER
       return
    endif

    ! unknown type:

    type_id=0

    return
  end subroutine get_h5_type_group

!-----------------------------------------------------------------------
! Does an HDF5 object exist?
! At the moment this can be used to test for the existence of a 
!    FILE      - see h5file_exists    in taehdf5.cpp,
!    GROUP     - see h5group_exists   in taehdf5.cpp, or 
!    DATASET   - see h5dataset_exists in taehdf5.cpp
!    ATTRIBUTE - see h5attribute_exists

  FUNCTION h5object_exists(object_type,loc_id,name) RESULT(exists)
    INTEGER,         intent(IN):: object_type ! H5I_<type>_F
    INTEGER(HID_T),  intent(IN):: loc_id      ! Location identifier
    CHARACTER(len=*),intent(IN):: name        ! Name of object
    LOGICAL                    :: exists
    INTEGER(HID_T):: obj_id
    INTEGER       :: ierr

    exists=.FALSE.

    ! Turn off error printing

    call H5Eset_auto_f(0,ierr)

    ! Try to open existing file, group, or dataset:

    if     (object_type==H5I_FILE_F)   then ; CALL h5fopen_f(name,H5F_ACC_RDWR_F,obj_id,ierr)
    elseif (object_type==H5I_GROUP_F)  then ; CALL h5gopen_f(loc_id,name,obj_id,ierr)
    elseif (object_type==H5I_DATASET_F)then ; CALL h5dopen_f(loc_id,name,obj_id,ierr)
    elseif (object_type==H5I_ATTR_F)   then ; CALL h5aopen_name_f(loc_id,name,obj_id,ierr)
    else
       print*,'H5OBJECT_EXISTS - UNKNOWN OBJECT TYPE'
       exists=.FALSE.
       obj_id=-1
    endif

    if (ierr<0) then ! object does not exist
       exists=.FALSE.
    else               ! object exists - open was successful, so close it again
       exists=.TRUE.

       if     (object_type==H5I_FILE_F)   then ; CALL h5fclose_f(obj_id,ierr)
       elseif (object_type==H5I_GROUP_F)  then ; CALL h5gclose_f(obj_id,ierr)
       elseif (object_type==H5I_DATASET_F)then ; CALL h5dclose_f(obj_id,ierr)
       elseif (object_type==H5I_ATTR_F)   then ; CALL h5aclose_f(obj_id,ierr)
       endif
    endif

    ! Turn error printing back on

    call H5Eset_auto_f(1,ierr)

    RETURN
  END FUNCTION h5object_exists

!-----------------------------------------------------------------------
! Use h5object_exists to test for a FILE, GROUP, DATASET or ATTRIBUTE

  FUNCTION h5file_exists(name) RESULT(exists)
    CHARACTER(len=*),intent(IN):: name        ! Name of object
    LOGICAL                    :: exists

    exists=h5object_exists(H5I_FILE_F,0,name)

    RETURN
  END FUNCTION h5file_exists


  FUNCTION h5group_exists_id(loc_id,name) RESULT(exists)
    INTEGER(HID_T),  intent(IN):: loc_id      ! Location identifier
    CHARACTER(len=*),intent(IN):: name        ! Name of object
    LOGICAL                    :: exists

    exists=h5object_exists(H5I_GROUP_F,loc_id,name)
!!$    if (exists) then
!!$       print*,'___'//name//'___ group exists'
!!$    else
!!$       print*,'___'//name//'___ group non-existent'
!!$    endif

    RETURN
  END FUNCTION h5group_exists_id


  FUNCTION h5dataset_exists_id(loc_id,name) RESULT(exists)
    INTEGER(HID_T),  intent(IN):: loc_id      ! Location identifier
    CHARACTER(len=*),intent(IN):: name        ! Name of object
    LOGICAL                    :: exists

    exists=h5object_exists(H5I_DATASET_F,loc_id,name)

    RETURN
  END FUNCTION h5dataset_exists_id


  FUNCTION h5attribute_exists(loc_id,name) RESULT(exists)
    INTEGER(HID_T),  intent(IN):: loc_id      ! Location identifier
    CHARACTER(len=*),intent(IN):: name        ! Name of object
    LOGICAL                    :: exists

    exists=h5object_exists(H5I_ATTR_F,loc_id,name)

    RETURN
  END FUNCTION h5attribute_exists

  !-----------------------------------------------------------------------
  ! Does dataset exist, by filename and path
  !

  FUNCTION h5dataset_exists_path(file,path) RESULT(exists)
    CHARACTER(len=*),intent(IN):: file,path ! filename and internal HDF5 path input
    LOGICAL                    :: exists
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0.AND.h5dataset_exists(ids(last),name)) then
       exists=.TRUE.
    else
       exists=.FALSE.
    endif

    RETURN
  END FUNCTION h5dataset_exists_path

  FUNCTION h5group_exists_path(file,path) RESULT(exists)
    CHARACTER(len=*),intent(IN):: file,path ! filename and internal HDF5 path input
    LOGICAL                    :: exists
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0.AND.h5group_exists(ids(last),name)) then
       exists=.TRUE.
    else
       exists=.FALSE.
    endif

    RETURN
  END FUNCTION h5group_exists_path

  !-----------------------------------------------------------------------
  ! Open a "path" in a file, i.e. open/create a file and then 
  ! recursively open/create an array of groups inside the file,
  ! by parsing the PATH argument.
  !
  ! Example:
  ! INPUT:  ('test.h5','grand/parent/child/X',...
  ! OUTPUT:  ids=loc_ids of [test.h5, grand, parent, child]
  !          leaf=4,
  !          name='X'
  !
  subroutine h5open_path(file,path,ids,leaf,name)
    CHARACTER(len=*)            ,intent(IN) ::file,path  ! filename and internal HDF5 path
    INTEGER(HID_T),dimension(10),intent(OUT)::ids        ! loc_ids of each branch in path
    integer                     ,intent(OUT)::leaf       ! leaf position on path
    CHARACTER(len=*)            ,intent(OUT)::name       ! name of leaf at end of path
    INTEGER::ierr,in,i,islash,ng
    CHARACTER(len=80)::gpath

    ! use -1 as an end marker:

    ids=-1

    ! open interface

    CALL h5open_f(ierr)

    ! open or create file

    if (h5file_exists(file)) then
       CALL h5fopen_f(file,H5F_ACC_RDWR_F,ids(1),ierr)
    else
       CALL h5fcreate_f(file,H5F_ACC_TRUNC_F,ids(1),ierr)
    endif

    ! parse path - for every '/' character, open or create a subgroup

    in=1

    if (INDEX(path,'/').eq.1) in=2  ! ignore leading '/' in group path

    do i=2,size(ids)
       
       islash=INDEX(path(in:),'/')

       if (islash > 0) then

          ! open or create a group path(in:in+islash-1)

          gpath=path(in:in+islash-2)

          if (h5group_exists(ids(i-1),TRIM(gpath))) then
             call h5gopen_f(ids(i-1),TRIM(gpath),ids(i),ierr) ! Open existing group
          else
             call h5gcreate_f(ids(i-1),TRIM(gpath),ids(i),ierr) ! Create the group
          endif

          in=in+islash
       else
          exit ! exit loop
       endif
    enddo

    ! name and position of leaf at end of path

    name=path(in:)
    leaf=i-1

    return
  end subroutine h5open_path

  !-----------------------------------------------------------------------
  ! Read a "path" in a file, i.e. open a file and then 
  ! recursively open an array of groups inside the file,
  ! by parsing the PATH argument.
  !
  ! This is the read-only version of h5open_path, above.
  !
  ! Example:
  ! INPUT:  ('test.h5','grand/parent/child/X',...
  ! OUTPUT:  ids=loc_ids of [test.h5, grand, parent, child]
  !          leaf=4,
  !          name='X'
  !
  subroutine h5read_path(file,path,ids,leaf,name)
    CHARACTER(len=*)            ,intent(IN) ::file,path  ! filename and internal HDF5 path
    INTEGER(HID_T),dimension(10),intent(OUT)::ids        ! loc_ids of each branch in path
    integer                     ,intent(OUT)::leaf       ! leaf position on path
    CHARACTER(len=*)            ,intent(OUT)::name       ! name of leaf at end of path
    INTEGER::ierr,in,i,islash,ng
    CHARACTER(len=80)::gpath

    ! use -1 as an end marker:

    ids=-1

    ! open interface

    CALL h5open_f(ierr)

    ! open or create file

    if (.NOT.h5file_exists(file)) then
       ids=-1
       name=''
       leaf=0
       print*,'H5READ_PATH - FILE ',file,' NOT FOUND'
       return
    else
       CALL h5fopen_f(file,H5F_ACC_RDWR_F,ids(1),ierr)
    endif

    ! parse path - for every '/' character, open or create a subgroup

    in=1

    if (INDEX(path,'/').eq.1) in=2  ! ignore leading '/' in group path

    do i=2,size(ids)
       
       islash=INDEX(path(in:),'/')

       if (islash > 0) then

          ! open a group path(in:in+islash-1)

          gpath=path(in:in+islash-2)

          if (.NOT.h5group_exists(ids(i-1),TRIM(gpath))) then
             ids(i-1)=-1
             name=path(in:)
             leaf=i-1
             if (.not.taehdf5_quiet) print*,'H5READ_PATH - PATH ',path,' NOT FOUND IN FILE ',file
             return
          else
             call h5gopen_f(ids(i-1),TRIM(gpath),ids(i),ierr) ! Open existing group
          endif

          in=in+islash
       else
          exit ! exit loop
       endif
    enddo

    ! name and position of leaf at end of path

    name=path(in:)
    leaf=i-1

    return
  end subroutine h5read_path

  !-----------------------------------------------------------------------
  ! Close an array of ids which were opened/created by h5open_path.
  !
  subroutine h5close_path(ids)
    INTEGER(HID_T),dimension(10),intent(INOUT)::ids
    INTEGER::i,ierr

    ! close groups, file, and interface

    do i=size(ids),2,-1
       if (ids(i) > 0) call h5gclose_f(ids(i),ierr)
    enddo
    if (ids(1)>0) call h5fclose_f(ids(1),ierr)
    call h5close_f(ierr)

    return
  end subroutine h5close_path

  !=============================================================
  ! H5COPY_DATA* wrappers
  !
  ! Copy data from one file to another
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Copy data at end of path from src file to dest file
  !
  ! example arguments:  ('dump000.bin','dump002.h5','vacuum/psi')
  ! or                  ('dump000.bin','dump002.h5','Deepak/dens_n','neutral/en_n')
  subroutine h5copy_data(src,dest,srcpath,destpath_in)
    CHARACTER(len=*),intent(IN) :: src,dest,srcpath,destpath_in  ! filenames and internal HDF5 path input
    OPTIONAL::destpath_in
    CHARACTER(len=80)::destpath
    INTEGER(HSIZE_T),DIMENSION(4):: dims          ! Array of dims
    INTEGER                      :: num_dims      ! Number of dims
    INTEGER(HID_T)               :: type
    INTEGER::I0
    INTEGER,dimension(:),    allocatable::I1
    INTEGER,dimension(:,:),  allocatable::I2
!!$    INTEGER,dimension(:,:,:),allocatable::I3
    REAL(rn)::D0
    REAL(rn),dimension(:),    allocatable::D1
    REAL(rn),dimension(:,:),  allocatable::D2
    REAL(rn),dimension(:,:,:),allocatable::D3

    if (present(destpath_in)) then
       destpath=destpath_in
    else
       destpath=srcpath
    endif

    call get_h5_dims(src,srcpath,dims,num_dims)
    call get_h5_type(src,srcpath,type)

    if (num_dims<0.OR.num_dims>3.OR.(type/=H5T_NATIVE_DOUBLE.AND.type/=H5T_NATIVE_INTEGER)) then
       print*,'FAILED h5copy_data(',src,',',dest,',',srcpath,')'
       return
    endif

    select case(num_dims)
    case(0)
       if (type==H5T_NATIVE_DOUBLE) then
          call h5read_data(src, D0,srcpath)
          call h5dump_data(dest,D0,destpath)
       else if (type==H5T_NATIVE_INTEGER) then
          call h5read_data(src, I0,srcpath)
          call h5dump_data(dest,I0,destpath)
       endif
    case(1)
       if (type==H5T_NATIVE_DOUBLE) then
          allocate(D1(dims(1)))
          call h5read_data(src, D1,srcpath)
          call h5dump_data(dest,D1,destpath)
          deallocate(D1)
       else if (type==H5T_NATIVE_INTEGER) then
          allocate(I1(dims(1)))
          call h5read_data(src, I1,srcpath)
          call h5dump_data(dest,I1,destpath)
          deallocate(I1)
       endif
    case(2)
       if (type==H5T_NATIVE_DOUBLE) then
          allocate(D2(dims(1),dims(2)))
          call h5read_data(src, D2,srcpath)
          call h5dump_data(dest,D2,destpath)
          deallocate(D2)
       else if (type==H5T_NATIVE_INTEGER) then
          allocate(I2(dims(1),dims(2)))
          call h5read_data(src, I2,srcpath)
          call h5dump_data(dest,I2,destpath)
          deallocate(I2)
       endif
    case(3)
       if (type==H5T_NATIVE_DOUBLE) then
          allocate(D3(dims(1),dims(2),dims(3)))
          call h5read_data(src, D3,srcpath)
          call h5dump_data(dest,D3,destpath)
          deallocate(D3)
       else if (type==H5T_NATIVE_INTEGER) then
!!$          allocate(I3(dims(1),dims(2),dims(3)))
!!$          call h5read_data(src, I3,srcpath)
!!$          call h5dump_data(dest,I3,destpath)
!!$          deallocate(I3)
          print*,'taehdf5 - do not have 3D integer interface yet.'
          stop
       endif
    end select

    return
  end subroutine h5copy_data

  !=============================================================
  ! H5DUMP_*_AUTO wrappers
  !
  ! Open/Create File, Open/Create group path, dump data at end of path
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_double_0d_auto(file,F,path)
    REAL(kind=rn),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_double_0d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_double_0d_auto
  !-----------------------------------------------------------------------
  ! Dump F into file at end of path, with description argument
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue','Put a desription here')
  !
  subroutine h5dump_data_double_0d_auto_descrip(file,F,path,descrip)
    REAL(kind=rn),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*),intent(IN):: file,path,descrip  ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids              ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name             ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_double_0d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_double_0d_auto_descrip

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_double_1d_auto(file,F,path)
    REAL(kind=rn),dimension(:),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_double_1d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_double_1d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_double_2d_auto(file,F,path)
    REAL(kind=rn),dimension(:,:),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_double_2d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_double_2d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_double_3d_auto(file,F,path)
    REAL(kind=rn),dimension(:,:,:),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_double_3d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_double_3d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_int_0d_auto(file,F,path)
    integer,intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_int_0d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_int_0d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_int_1d_auto(file,F,path)
    integer,dimension(:),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_int_1d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_int_1d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
!!$  subroutine h5dump_data_longlong_1d_auto(file,F,path)
!!$    integer(longlong),dimension(:),intent(IN):: F         ! value to put on leaf
!!$    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
!!$    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
!!$    CHARACTER(len=80)            :: name                ! name of leaf at end of path
!!$    integer :: last
!!$
!!$    ! open path
!!$    ! dump the data to the leaf at the end of the path
!!$    ! close path
!!$
!!$    call h5open_path(file,path,ids,last,name)
!!$    call h5dump_data_longlong_1d(ids(last),F,name)
!!$    call h5close_path(ids)
!!$
!!$    return
!!$  end subroutine h5dump_data_longlong_1d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_int_2d_auto(file,F,path)
    integer,dimension(:,:),intent(IN):: F         ! value to put on leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_int_2d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_int_2d_auto

  !=============================================================
  ! H5WR_ file/group wrappers
  ! 
  ! Allow Create/Open of a File or Group from the same interface
  !    WR=0 -> WRITE
  !    WR=1 -> READ
  !=============================================================

  !-----------------------------------------------------------------------
  ! WR=0 -> Create an HDF5 file
  ! WR=1 -> Open an existing HDF5 file
  !
  subroutine h5wr_file(WR,fid,name)
    INTEGER,         intent(IN)   :: WR    ! =1 or 0
    INTEGER(HID_T),  intent(INOUT):: fid   ! File identifier
    CHARACTER(len=*),intent(IN)   :: name  ! Name of file

    INTEGER       :: ierr  ! Error code 

    select case(WR)
    case(0) ! Create the HDF5 file
       CALL h5fcreate_f(name,H5F_ACC_TRUNC_F,fid,ierr)
    case(1) ! Open existing HDF5 file
       CALL h5fopen_f(name,H5F_ACC_RDWR_F,fid,ierr)
    case(-2:-1) ! ASCII Debug output
       open(FILE=name,UNIT=fid,STATUS='unknown') 
    end select
    return
  end subroutine h5wr_file

  !-----------------------------------------------------------------------
  ! WR=0 -> Create an HDF5 group in file fid
  ! WR=1 -> Open an existing HDF5 group in file fid
  !
  subroutine h5wr_group(WR,fid,gid,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: fid   ! File identifier
    INTEGER(HID_T),  intent(INOUT):: gid   ! Group identifier
    CHARACTER(len=*),intent(IN)   :: name  ! Name of group

    INTEGER       :: ierr  ! Error code 

    select case(WR)
    case(0:1)
       if (h5group_exists(fid,name)) then
          call h5gopen_f(fid,name,gid,ierr) ! Open existing group
       else
          call h5gcreate_f(fid,name,gid,ierr) ! Create the group
       endif
!!$    case(0)
!!$       call h5gcreate_f(fid,name,gid,ierr) ! Create the group
!!$    case(1)
!!$       call h5gopen_f(fid,name,gid,ierr) ! Open existing group
    case(-2:-1) ! for ASCII debug output
       gid=fid
    end select
    return
  end subroutine h5wr_group

  !=============================================================
  ! H5WR_ATTR wrappers
  ! 
  ! Allow Write/Read of a Scalar int/double/string ATTRIBUTE from the same interface.
  !    WR=0 -> WRITE
  !    WR=1 -> READ
  !
  ! These have a common interface:
  !
  !    call h5wr_attr(WR,group_id,attribute_value,attribute_name)
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! WR=0 -> Write an HDF5 integer attribute in   group gid
  ! WR=1 -> Read  an HDF5 integer attribute from group gid
  !
  subroutine h5wr_attr_int(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    INTEGER,         intent(INOUT):: F     ! Value of attribute
    CHARACTER(len=*),intent(IN)   :: name  ! Name of attribute

    select case(WR)
    case(0)
       call h5dump_attr_int(gid,F,name) ! Write the attribute
    case(1)
       call h5read_attr_int(gid,F,name) ! Read the attribute
    case(-1) ! Ascii debug
       write(gid,*) F,name
    case(-2)
       read(gid,*) F
    end select
    return
  end subroutine h5wr_attr_int
  !-----------------------------------------------------------------------
  ! WR=0 -> Write an HDF5 double attribute in   group gid
  ! WR=1 -> Read  an HDF5 double attribute from group gid
  !
  subroutine h5wr_attr_double(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    REAL(rn)        ,intent(INOUT):: F     ! Value of attribute
    CHARACTER(len=*),intent(IN)   :: name  ! Name of attribute

    select case(WR)
    case(0)
       call h5dump_attr_double(gid,F,name) ! Write the attribute
    case(1)
       call h5read_attr_double(gid,F,name) ! Read the attribute
    case(-1) ! Ascii debug
       write(gid,*) F,name
    case(-2)
       read(gid,*) F
    end select
    return
  end subroutine h5wr_attr_double
  !-----------------------------------------------------------------------
  ! WR=0 -> Write an HDF5 string attribute in   group gid
  ! WR=1 -> Read  an HDF5 string attribute from group gid
  !
  subroutine h5wr_attr_string(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    CHARACTER(len=*),intent(INOUT):: F     ! Value of attribute
    CHARACTER(len=*),intent(IN)   :: name  ! Name of attribute

    select case(WR)
    case(0)
       call h5dump_attr_string(gid,F,name) ! Write the attribute
    case(1)
       call h5read_attr_string(gid,F,name) ! Read the attribute
    case(-1) ! Ascii debug
       write(gid,*) F,name
    case(-2)
       read(gid,*) F
    end select
    return
  end subroutine h5wr_attr_string


  !=============================================================
  ! H5READ_*_AUTO wrappers
  !
  ! Open File, Open group path, read data at end of path
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_double_0d_auto(file,F,path)
    REAL(kind=rn),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    !INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

!!$    call get_h5_type(file,path,type)
!!$    if (type.ne.H5T_NATIVE_DOUBLE) then
!!$       print*,'WRONG TYPE PASSED TO h5read_data_double_0d_auto(',TRIM(file),",F,'",TRIM(path),"')"
!!$       print*,'maybe ',path,' is an integer or string?'
!!$       stop
!!$    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0.AND.h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else if (ids(last)>0.AND.h5attribute_exists(ids(last),name)) then
       call h5read_attr(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    return
  end subroutine h5read_data_double_0d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_double_1d_auto(file,F,path)
    REAL(kind=rn),dimension(:),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    !INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

!!$    call get_h5_type(file,path,type)
!!$    if (type.ne.H5T_NATIVE_DOUBLE) then
!!$       print*,'WRONG TYPE PASSED TO h5read_data_double_1d_auto(',TRIM(file),",F,'",TRIM(path),"')"
!!$       print*,'maybe ',path,' is an integer or string?'
!!$       stop
!!$    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0 .AND. h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    return
  end subroutine h5read_data_double_1d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_double_2d_auto(file,F,path,found)
    REAL(kind=rn),dimension(:,:),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)            ,intent(IN) :: file,path ! filename and internal HDF5 path input
    logical,            optional,intent(OUT):: found
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

    call get_h5_type(file,path,type)
    if (type.ne.H5T_NATIVE_DOUBLE) then
       if (present(found)) then
          FOUND=.FALSE.
          return
       else
          print*,'WRONG TYPE PASSED TO h5read_data_double_2d_auto(',TRIM(file),",F,'",TRIM(path),"')"
          print*,'maybe ',path,' is an integer or string?'
          stop
       endif
    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0 .AND. h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    if (present(found)) then
       FOUND=.TRUE.
    endif
    return
  end subroutine h5read_data_double_2d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_double_3d_auto(file,F,path,found)
    REAL(kind=rn),dimension(:,:,:),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)              ,intent(IN) :: file,path ! filename and internal HDF5 path input
    logical,              optional,intent(OUT):: found
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

    call get_h5_type(file,path,type)
    if (type.ne.H5T_NATIVE_DOUBLE) then
       if (present(found)) then
          FOUND=.FALSE.
          return
       else
          print*,'WRONG TYPE PASSED TO h5read_data_double_3d_auto(',TRIM(file),",F,'",TRIM(path),"')"
          print*,'maybe ',path,' is an integer or string?'
          stop
       endif
    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0 .AND. h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    if (present(found)) then
       FOUND=.TRUE.
    endif
    return
  end subroutine h5read_data_double_3d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_double_4d_auto(file,F,path,found)
    REAL(kind=rn),dimension(:,:,:,:),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)                ,intent(IN) :: file,path ! filename and internal HDF5 path input
    logical,                optional,intent(OUT):: found
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

    call get_h5_type(file,path,type)
    if (type.ne.H5T_NATIVE_DOUBLE) then
       if (present(found)) then
          FOUND=.FALSE.
          return
       else
          print*,'WRONG TYPE PASSED TO h5read_data_double_4d_auto(',TRIM(file),",F,'",TRIM(path),"')"
          print*,'maybe ',path,' is an integer or string?'
          stop
       endif
    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0 .AND. h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    if (present(found)) then
       FOUND=.TRUE.
    endif
    return
  end subroutine h5read_data_double_4d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_int_0d_auto(file,F,path)
    INTEGER,intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    !INTEGER(HID_T)               :: type

    ! open path
    ! read the data OR ATTRIBUTE from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0.AND.h5dataset_exists(ids(last),name)) then

!!$       ! we can check type on datasets, so check that user has asked for the right type
!!$       call get_h5_type(file,path,type)
!!$       if (type.ne.H5T_NATIVE_INTEGER) then
!!$          print*,'WRONG TYPE PASSED TO h5read_data_int_0d_auto(',TRIM(file),",F,'",TRIM(path),"')"
!!$          print*,'maybe ',path,' is an REAL or string?'
!!$          stop
!!$       endif

       call h5read_data(ids(last),F,name)
    else if (ids(last)>0.AND.h5attribute_exists(ids(last),name)) then
       ! don't think we can check type on attributes, so just read it:
       call h5read_attr(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    return
  end subroutine h5read_data_int_0d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_int_1d_auto(file,F,path)
    INTEGER,dimension(:),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

    call get_h5_type(file,path,type)
    if (type.ne.H5T_NATIVE_INTEGER) then
       print*,'WRONG TYPE PASSED TO h5read_data_int_1d_auto(',TRIM(file),",F,'",TRIM(path),"')"
       print*,'maybe ',path,' is an REAL or string?'
       stop
    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0.AND.h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    return
  end subroutine h5read_data_int_1d_auto

  !-----------------------------------------------------------------------
  ! Read F from path in file
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5read_data_int_2d_auto(file,F,path)
    INTEGER,dimension(:,:),intent(OUT):: F         ! value read from leaf
    CHARACTER(len=*)            ,intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last
    INTEGER(HID_T)               :: type

    ! check that user has asked for the right type

    call get_h5_type(file,path,type)
    if (type.ne.H5T_NATIVE_INTEGER) then
       print*,'WRONG TYPE PASSED TO h5read_data_int_2d_auto(',TRIM(file),",F,'",TRIM(path),"')"
       print*,'maybe ',path,' is an REAL or string?'
       stop
    endif

    ! open path
    ! read the data from the leaf at the end of the path
    ! close path

    call h5read_path(file,path,ids,last,name)
    if (ids(last)>0 .AND. h5dataset_exists(ids(last),name)) then
       call h5read_data(ids(last),F,name)
    else
       print*,'NOT FOUND: H5READ_DATA('//TRIM(file)//',F,'//TRIM(path)//')'
    endif
    call h5close_path(ids)

    return
  end subroutine h5read_data_int_2d_auto

  !=============================================================
  ! H5WR_DATA wrappers
  ! 
  ! Allow Write/Read of 0D/1D/2D int/double/string DATA from the same interface.
  !    WR=0 -> WRITE
  !    WR=1 -> READ
  !
  ! These have a common interface:
  !
  !    call h5wr_data(WR,group_id,data_value(s),data_name)
  !
  !=============================================================


  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 0D logical to   HDF5 group gid
  ! WR=1 -> Read  a 0D logical from HDF5 group gid
  !
  subroutine h5wr_data_logical_0d(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    LOGICAL,         intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_logical_0d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_logical_0d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) F,name
    case(-2)
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_logical_0d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 1D logical to   HDF5 group gid
  ! WR=1 -> Read  a 1D logical from HDF5 group gid
  !
  subroutine h5wr_data_logical_1d(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    LOGICAL,DIMENSION(:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_logical_1d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_logical_1d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*)
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_logical_1d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 2D logical to   HDF5 group gid
  ! WR=1 -> Read  a 2D logical from HDF5 group gid
  !
  subroutine h5wr_data_logical_2d(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    LOGICAL,DIMENSION(:,:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_logical_2d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_logical_2d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_logical_2d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 0D integer to   HDF5 group gid
  ! WR=1 -> Read  a 0D integer from HDF5 group gid
  !
  subroutine h5wr_data_int_0d(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    INTEGER,         intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_int_0d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_int_0d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) F,name
    case(-2)
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_int_0d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 1D integer to   HDF5 group gid
  ! WR=1 -> Read  a 1D integer from HDF5 group gid
  !
  subroutine h5wr_data_int_1d(WR,gid,F,name)
    INTEGER,             intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),      intent(IN)   :: gid   ! Group identifier
    INTEGER,DIMENSION(:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),    intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_int_1d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_int_1d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_int_1d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 2D integer to   HDF5 group gid
  ! WR=1 -> Read  a 2D integer from HDF5 group gid
  !
  subroutine h5wr_data_int_2d(WR,gid,F,name)
    INTEGER,               intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),        intent(IN)   :: gid   ! Group identifier
    INTEGER,DIMENSION(:,:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),      intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_int_2d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_int_2d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_int_2d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 0D double to   HDF5 group gid
  ! WR=1 -> Read  a 0D double from HDF5 group gid
  !
  subroutine h5wr_data_double_0d(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    REAL(rn),        intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_double_0d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_double_0d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_double_0d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 1D double to   HDF5 group gid
  ! WR=1 -> Read  a 1D double from HDF5 group gid
  !
  subroutine h5wr_data_double_1d(WR,gid,F,name)
    INTEGER,              intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),       intent(IN)   :: gid   ! Group identifier
    REAL(rn),DIMENSION(:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),     intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_double_1d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_double_1d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_double_1d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 2D double to   HDF5 group gid
  ! WR=1 -> Read  a 2D double from HDF5 group gid
  !
  subroutine h5wr_data_double_2d(WR,gid,F,name)
    INTEGER,                intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),         intent(IN)   :: gid   ! Group identifier
    REAL(rn),DIMENSION(:,:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),       intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_double_2d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_double_2d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_double_2d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 3D double to   HDF5 group gid
  ! WR=1 -> Read  a 3D double from HDF5 group gid
  !
  subroutine h5wr_data_double_3d(WR,gid,F,name)
    INTEGER,                  intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),           intent(IN)   :: gid   ! Group identifier
    REAL(rn),DIMENSION(:,:,:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),         intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_double_3d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_double_3d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_double_3d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 0D string to   HDF5 group gid
  ! WR=1 -> Read  a 0D string from HDF5 group gid
  !
  subroutine h5wr_data_string_0d(WR,gid,F,name)
    INTEGER,         intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),  intent(IN)   :: gid   ! Group identifier
    CHARACTER(len=*),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_string_0d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_string_0d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_string_0d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 1D string to   HDF5 group gid
  ! WR=1 -> Read  a 1D string from HDF5 group gid
  !
  subroutine h5wr_data_string_1d(WR,gid,F,name)
    INTEGER,                      intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),               intent(IN)   :: gid   ! Group identifier
    CHARACTER(len=*),DIMENSION(:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),             intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_string_1d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_string_1d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_string_1d
  !-----------------------------------------------------------------------
  ! WR=0 -> Write a 2D string to   HDF5 group gid
  ! WR=1 -> Read  a 2D string from HDF5 group gid
  !
  subroutine h5wr_data_string_2d(WR,gid,F,name)
    INTEGER,                        intent(IN)   :: WR    ! W or R
    INTEGER(HID_T),                 intent(IN)   :: gid   ! Group identifier
    CHARACTER(len=*),DIMENSION(:,:),intent(INOUT):: F     ! Value of input
    CHARACTER(len=*),               intent(IN)   :: name  ! Name of input

    select case(WR)
    case(0)
       call h5dump_data_string_2d(gid,F,name) ! Write the data
    case(1)
       call h5read_data_string_2d(gid,F,name) ! Read the data
    case(-1) ! Ascii debug
       write(gid,*) name
       write(gid,*) F
    case(-2)
       read(gid,*) 
       read(gid,*) F
    end select
    return
  end subroutine h5wr_data_string_2d

  !=============================================================
  ! H5DUMP_ATTR int/double/string wrappers
  ! 
  ! These are Shortcuts to WRITE attributes to a group or object.
  ! These have a uniform interface:
  !
  !    call h5dump_attr(group_id,attribute_value,attribute_name)
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Store an integer attribute

  subroutine h5dump_attr_int(loc_id,f,name)
    ! arguments
    INTEGER         ,intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: loc_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of attribute
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier

    dims=1

    ! If the named attribute already exists at this location then delete it.

    if (h5attribute_exists(loc_id,name)) call H5Adelete_F(loc_id, name,hdferr)

    ! Store F in a space in loc_id

    call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
    call h5acreate_f(loc_id,name,H5T_NATIVE_INTEGER,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_attr_int


  !-----------------------------------------------------------------------
  ! Store an integer attribute

  subroutine h5dump_attr_double(loc_id,f,name)
    ! arguments
    REAL(kind=rn)   ,intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: loc_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of attribute
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier

    dims=1

    ! If the named attribute already exists at this location then delete it.

    if (h5attribute_exists(loc_id,name)) call H5Adelete_F(loc_id, name,hdferr)

    ! Store F in a space in loc_id

    call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
    call h5acreate_f(loc_id,name,H5T_NATIVE_DOUBLE,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_attr_double

  !-----------------------------------------------------------------------
  ! Store a CHARACTER string attribute

  subroutine h5dump_attr_string(loc_id,f,name)
    ! arguments
    CHARACTER(len=*) ,intent(IN):: f,name
    INTEGER(HID_T),intent(IN)   :: loc_id
    ! local variables
    INTEGER :: hdferr         ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims
    INTEGER(HID_T):: type_id        ! ID of variable length string type
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id     ! Dataset identifier
    INTEGER(HID_T):: attr_id        ! Attribute identifier

    dims=1

    if (LEN(f) < 1) then
          print*,'IGNORING h5dump_attr_string(',loc_id,',f,',name,') - LEN(F) must be > 0'
    endif

    ! If the named attribute already exists at this location then delete it.

    if (h5attribute_exists(loc_id,name)) call H5Adelete_F(loc_id, name,hdferr)

    ! Create new type_id variable length string data type for F

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN_TRIM(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Store F in a space in loc_id

    call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
    call h5acreate_f(loc_id,name,type_id,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,type_id, TRIM(F), dims,hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5dump_attr_string

  !=============================================================
  ! H5READ_ATTR int/double/string wrappers
  ! 
  ! These are Shortcuts to READ attributes from a group or object.
  ! These have a uniform interface:
  !
  !    call h5read_attr(loc_id,attribute_value,attribute_name)
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Store an integer attribute

  subroutine h5read_attr_int(loc_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: loc_id
    INTEGER         ,intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of attribute
    INTEGER(HID_T):: int_id             ! ID of integer data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier

    dims=1

    ! Open the named attribute and read F

    call h5aopen_name_f(loc_id,name,attr_id,hdferr)
    call h5aread_f(attr_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
    call h5aclose_f(attr_id,hdferr)
    return
  end subroutine h5read_attr_int


  !-----------------------------------------------------------------------
  ! Store an integer attribute

  subroutine h5read_attr_double(loc_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: loc_id
    REAL(kind=rn)   ,intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of attribute
    INTEGER(HID_T):: double_id          ! ID of double data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier

    dims=1

    ! Open the named attribute and read F

    call h5aopen_name_f(loc_id,name,attr_id,hdferr)
    call h5aread_f(attr_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5aclose_f(attr_id,hdferr)
    return
  end subroutine h5read_attr_double

  !-----------------------------------------------------------------------
  ! Read a CHARACTER string attribute

  subroutine h5read_attr_string(loc_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN)   :: loc_id
    CHARACTER(len=*),intent(INOUT):: f        ! (needs IN to get LEN(F))
    CHARACTER(len=*),intent(IN)   :: name
    ! local variables
    INTEGER :: hdferr         ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims
    INTEGER(HID_T):: attr_id        ! Attribute identifier
    INTEGER(HID_T):: type_id        ! Data Type identifier

    dims=1

    ! Create a type of LEN(F)*H5T_NATIVE_CHARACTER

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Open the named attribute and read F

    call h5aopen_name_f(loc_id,name,attr_id,hdferr)
    call h5aread_f(attr_id,type_id, F, dims,hdferr)
    call h5aclose_f(attr_id,hdferr)

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5read_attr_string

  !=============================================================
  ! H5DUMP_DATA 0D/1D/2D int/double/string wrappers
  ! 
  ! These are Shortcuts to WRITE 0D/1D/2D data of type int/double/string
  ! to a group or object.
  ! These have a uniform interface:
  !
  !    call h5dump_data(group_id,data_value,data_name)
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Store an logical as int as there is no bool native type (?)

  subroutine h5dump_data_logical_0d(group_id,f,name)
    ! arguments
    LOGICAL         ,intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: Fint ! write an int as there is no bool native type

    ! Map the bool to an int

    if (F) then 
       Fint=1
    else
       Fint=0
    end if

    call h5dump_data(group_id,Fint,name)

    return
  end subroutine h5dump_data_logical_0d
  !-----------------------------------------------------------------------
  ! Store an logical as int as there is no bool native type (?)

  subroutine h5dump_data_logical_1d(group_id,f,name)
    ! arguments
    LOGICAL,DIMENSION(:),intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER,DIMENSION(size(F)) :: Fint ! write an int as there is no bool native type

    ! Map the bool to an int

    where (F)
       Fint=1
    elsewhere
       Fint=0
    end where

    call h5dump_data(group_id,Fint,name)

    return
  end subroutine h5dump_data_logical_1d
  !-----------------------------------------------------------------------
  ! Store an logical as int as there is no bool native type (?)

  subroutine h5dump_data_logical_2d(group_id,f,name)
    ! arguments
    LOGICAL,DIMENSION(:,:),intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER,DIMENSION(size(F,DIM=1),size(F,DIM=2)) :: Fint ! write an int as there is no bool native type

    ! Map the bool to an int

    where (F)
       Fint=1
    elsewhere
       Fint=0
    end where

    call h5dump_data(group_id,Fint,name)

    return
  end subroutine h5dump_data_logical_2d
  !-----------------------------------------------------------------------
  ! Store an integer scalar

  subroutine h5dump_data_int_0d(group_id,f,name)
    ! arguments
    INTEGER         ,intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of integer data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=1

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if (num_dims > 1 .OR. sdims(1) > dims(1)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new dataspace and dataset
          call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new dataspace and dataset
       call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    ! Write to the dataset

    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_data_int_0d

  !-----------------------------------------------------------------------
  ! Store a 1D integer array

  subroutine h5dump_data_int_1d(group_id,f,name)
    ! arguments
    INTEGER,DIMENSION(:),intent(IN):: f
    CHARACTER(len=*)    ,intent(IN):: name
    INTEGER(HID_T)      ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of integer data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=SHAPE(F)
    if(.NOT.ALL(dims>0)) return

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if (num_dims > 1 .OR. sdims(1) .NE. dims(1)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new dataspace and dataset
          call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new dataspace and dataset
       call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    ! Write to the dataset

    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_data_int_1d

  !-----------------------------------------------------------------------
  ! Store a 1D integer array

!!$  subroutine h5dump_data_longlong_1d(group_id,f,name)
!!$    ! arguments
!!$    INTEGER(longlong),DIMENSION(:),intent(IN):: f
!!$    CHARACTER(len=*)    ,intent(IN):: name
!!$    INTEGER(HID_T)      ,intent(IN):: group_id
!!$    ! local variables
!!$    INTEGER :: hdferr                   ! Error code 
!!$    INTEGER(HID_T):: int_id             ! ID of integer data type
!!$    INTEGER(HID_T):: space_id           ! Dataspace identifier
!!$    INTEGER(HID_T):: dataset_id         ! Dataset identifier
!!$    INTEGER(HID_T):: attr_id            ! Attribute identifier
!!$    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
!!$    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
!!$    INTEGER                      :: num_dims ! Number of dims
!!$
!!$    dims=SHAPE(F)
!!$    if(.NOT.ALL(dims>0)) return
!!$
!!$    ! Open or Create a dataset of the right size in group_id
!!$
!!$    if (h5dataset_exists(group_id,name)) then
!!$       ! Check that the existing dataspace and dataset is the right size
!!$       ! - if not, then unlink it and create new one.
!!$
!!$       call get_h5_dims(group_id,name,sdims,num_dims)
!!$       if (num_dims > 1 .OR. sdims(1) .NE. dims(1)) then
!!$          ! Nuke this dataset
!!$          call H5Gunlink_f(group_id,name,hdferr)
!!$
!!$          ! Create a new dataspace and dataset
!!$          call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
!!$          call h5dcreate_f(group_id,name,H5T_NATIVE_LLONG,space_id,dataset_id,hdferr) ! create dataset_id
!!$       else
!!$          ! Everything cool - just open the dataset
!!$          call h5dopen_f(group_id,name,dataset_id,hdferr)
!!$          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
!!$       endif
!!$    else
!!$       ! Create a new dataspace and dataset
!!$       call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
!!$       call h5dcreate_f(group_id,name,H5T_NATIVE_LLONG,space_id,dataset_id,hdferr) ! create dataset_id
!!$    endif
!!$
!!$    ! Write to the dataset
!!$
!!$    call h5dwrite_f(dataset_id,H5T_NATIVE_LLONG, F, dims,hdferr)
!!$    call h5dclose_f(dataset_id,hdferr)
!!$    call h5sclose_f(space_id,hdferr)
!!$    return
!!$  end subroutine h5dump_data_longlong_1d

  !-----------------------------------------------------------------------
  ! Store a 2D integer array

  subroutine h5dump_data_int_2d(group_id,f,name)
    ! arguments
    INTEGER,DIMENSION(:,:),intent(IN):: f
    CHARACTER(len=*)      ,intent(IN):: name
    INTEGER(HID_T)        ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of integer data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(2):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=SHAPE(F)
    if(.not.ALL(dims>0)) return

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if ( num_dims > 2        .OR. &
            sdims(1).NE.dims(1) .OR. &
            sdims(2).NE.dims(2)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new 2D dataspace and dataset
          call h5screate_simple_f(2,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new 2D dataspace and dataset
       call h5screate_simple_f(2,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    ! Write to the dataset

    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)

    return
  end subroutine h5dump_data_int_2d

  !-----------------------------------------------------------------------
  ! Store a double scalar

  subroutine h5dump_data_double_0d(group_id,f,name)
    ! arguments
    REAL(rn)        ,intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=1

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if (num_dims > 1 .OR. sdims(1).NE.dims(1)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new 1D dataspace and dataset
          call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new 1D dataspace and dataset
       call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    ! Store F in a space in group_id

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_data_double_0d

  !-----------------------------------------------------------------------
  ! Store a 1D double array

  subroutine h5dump_data_double_1d(group_id,f,name)
    ! arguments
    REAL(rn),DIMENSION(:),intent(IN):: f
    CHARACTER(len=*)     ,intent(IN):: name
    INTEGER(HID_T)       ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=SHAPE(F)
    if(dims(1)<1) return

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if (num_dims > 1 .OR. sdims(1).NE.dims(1)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new 1D dataspace and dataset
          call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new 1D dataspace and dataset
       call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_data_double_1d

  !-----------------------------------------------------------------------
  ! Store a 2D double array

  subroutine h5dump_data_double_2d(group_id,f,name)
    ! arguments
    REAL(rn),DIMENSION(:,:),intent(IN):: f
    CHARACTER(len=*)       ,intent(IN):: name
    INTEGER(HID_T)         ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(2):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=SHAPE(F)
    if(.not.ALL(dims>0)) return

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if ( num_dims > 2        .OR. &
            sdims(1).NE.dims(1) .OR. &
            sdims(2).NE.dims(2)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new 2D dataspace and dataset
          call h5screate_simple_f(2,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new 2D dataspace and dataset
       call h5screate_simple_f(2,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_data_double_2d


  !-----------------------------------------------------------------------
  ! Store a 3D double array

  subroutine h5dump_data_double_3d(group_id,f,name)
    ! arguments
    REAL(rn),DIMENSION(:,:,:),intent(IN):: f
    CHARACTER(len=*)         ,intent(IN):: name
    INTEGER(HID_T)           ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: attr_id            ! Attribute identifier
    INTEGER(HSIZE_T),dimension(3):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims

    dims=SHAPE(F)
    if(.not.ALL(dims>0)) return

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if ( num_dims > 3        .OR. &
            sdims(1).NE.dims(1) .OR. &
            sdims(2).NE.dims(2) .OR. &
            sdims(3).NE.dims(3)) then

          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)

          ! Create a new 3D dataspace and dataset
          call h5screate_simple_f(3,dims,space_id,hdferr)                               ! create space_id
          call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
       else
          ! Everything cool - just open the dataset
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
       endif
    else
       ! Create a new 3D dataspace and dataset
       call h5screate_simple_f(3,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
    endif

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine h5dump_data_double_3d


  !-----------------------------------------------------------------------
  ! Store a scalar CHARACTER string F into file at end of path
  !
  ! example arguments:  ('test.h5', 'test string','grand/parent/child/Fvalue')
  !
  subroutine h5dump_data_string_0d_auto(file,F,path)
    CHARACTER(len=*),intent(IN):: file,F,path ! filename, value, and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5dump_data_string_0d(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5dump_data_string_0d_auto

  !-----------------------------------------------------------------------
  ! Store a scalar CHARACTER string

  subroutine h5dump_data_string_0d(group_id,f,name)
    ! arguments
    CHARACTER(len=*),intent(IN):: f,name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: type_id,etype_id   ! ID of data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims
    INTEGER(SIZE_T) :: esize
    LOGICAL :: create_new

    create_new=.FALSE.

    dims=1

    if (LEN(f) < 1) then
       print*,'IGNORING h5dump_data_string_0d(',group_id,',f,',name,') - LEN(F) must be > 0'
       print*,'xxx'//f//'xxx'
    endif

    ! Open or Create a dataset of the right size in group_id

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if (num_dims > 1 .OR. sdims(1).NE.dims(1)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)
          create_new=.TRUE.
       else
          ! Existing data space is right size
          create_new=.FALSE.
          ! But, check that the string length is correct
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_type_f(dataset_id, etype_id, hdferr)  ! get existing datatype
          call h5tget_size_f(etype_id,esize,hdferr)         ! get existing type size in bytes
          call h5tclose_f(etype_id,hdferr)
          call h5dclose_f(dataset_id,hdferr)
          if (esize.NE.int(LEN(F),SIZE_T)) then
             ! existing string type is wrong size - nuke this dataset
             call H5Gunlink_f(group_id,name,hdferr)
             create_new=.TRUE.
          endif
       endif
    else
       create_new=.TRUE.
    endif

    ! Create new type_id variable length string data type for F

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    if (create_new) then
       ! Create a new dataspace and dataset
       call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
       call h5dcreate_f(group_id,name,type_id,space_id,dataset_id,hdferr)
    else
       ! Everything cool - just open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    endif

    ! Store F in a space in group_id

    call h5dwrite_f(dataset_id,type_id, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5dump_data_string_0d

  !-----------------------------------------------------------------------
  ! Store a 1D CHARACTER string array

  subroutine h5dump_data_string_1d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN)             :: group_id
    CHARACTER(len=*),intent(IN),dimension(:):: f
    CHARACTER(len=*),intent(IN)             :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: type_id,etype_id   ! ID of data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HSIZE_T),dimension(1):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims
    LOGICAL :: create_new
    INTEGER,dimension(size(F))::lens    ! Length of each F(i)
    INTEGER :: i,length
    INTEGER(SIZE_T)::esize

    dims=SHAPE(F)
    if(ANY(dims==0)) return

    length=LEN(F(1))

    ! Check that all F(i) have the same LEN:

    if (size(F) > 0) then
       do i=1,size(F)
          lens(i)=LEN(F(i))
       enddo
       if (length < 1 .OR. ANY(lens/=lens(1))) then
          print*,'h5dump_data_string_1d - all LEN(F) must be the same and > 0'
          print*,lens
          return
       endif
    endif

    ! Decide whether to Open or Create a dataset

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if (num_dims > 1 .OR. sdims(1).NE.dims(1)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)
          create_new=.TRUE.
       else
          ! Existing data space is right size
          create_new=.FALSE.
          ! But, check that the string length is correct
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_type_f(dataset_id, etype_id, hdferr)  ! get existing datatype
          call h5tget_size_f(etype_id,esize,hdferr)         ! get existing type size in bytes
          call h5tclose_f(etype_id,hdferr)
          call h5dclose_f(dataset_id,hdferr)
          if (esize.NE.int(LEN(F(1)),SIZE_T)) then
             ! existing string type is wrong size - nuke this dataset
             call H5Gunlink_f(group_id,name,hdferr)
             create_new=.TRUE.
          endif
       endif
    else
       create_new=.TRUE.
    endif

    ! Create data type for CHARACTER(len=length) string

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,int(LEN(F(1)),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Open or Create the dataset

    if (create_new) then
       ! Create a new dataspace and dataset
       call h5screate_simple_f(1,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,type_id,space_id,dataset_id,hdferr) ! create dataset_id
    else
       ! Everything cool - just open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    endif

    ! Write to the dataset

    call h5dwrite_f(dataset_id,type_id, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5dump_data_string_1d

  !-----------------------------------------------------------------------
  ! Store a 2D CHARACTER string array

  subroutine h5dump_data_string_2d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN)               :: group_id
    CHARACTER(len=*),intent(IN),dimension(:,:):: f
    CHARACTER(len=*),intent(IN)               :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: type_id,etype_id   ! ID of double data type
    INTEGER(HID_T):: space_id           ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HSIZE_T),dimension(2):: dims ! Dimensions of input
    INTEGER(HSIZE_T),DIMENSION(4):: sdims    ! Dims of existing space
    INTEGER                      :: num_dims ! Number of dims
    LOGICAL :: create_new
    INTEGER,dimension(size(F,DIM=1),size(F,DIM=2))::lens    ! Length of each F(i)
    INTEGER :: i,j,length
    INTEGER(SIZE_T)::esize

    dims=SHAPE(F)
    if(ANY(dims==0)) return

    length=LEN(F(1,1))

    ! Check that all F(i) have the same LEN:

    if (size(F) > 0) then
       do i=1,size(F,DIM=1)
          do j=1,size(F,DIM=2)
             lens(i,j)=LEN(F(i,j))
          enddo
       enddo
       if (length < 1 .OR. ANY(lens/=lens(1,1))) then
          print*,'h5dump_data_string_2d - all LEN(F) must be the same and > 0'
          print*,lens
          return
       endif
    endif

    ! Decide whether to Open or Create a dataset

    if (h5dataset_exists(group_id,name)) then
       ! Check that the existing dataspace and dataset is the right size
       ! - if not, then unlink it and create new one.

       call get_h5_dims(group_id,name,sdims,num_dims)
       if ( num_dims > 2        .OR. &
            sdims(1).NE.dims(1) .OR. &
            sdims(2).NE.dims(2)) then
          ! Nuke this dataset
          call H5Gunlink_f(group_id,name,hdferr)
          create_new=.TRUE.
       else
          ! Existing data space is right size
          create_new=.FALSE.
          ! But, check that the string length is correct
          call h5dopen_f(group_id,name,dataset_id,hdferr)
          call h5dget_type_f(dataset_id, etype_id, hdferr)  ! get existing datatype
          call h5tget_size_f(etype_id,esize,hdferr)         ! get existing type size in bytes
          call h5tclose_f(etype_id,hdferr)
          call h5dclose_f(dataset_id,hdferr)
          if (esize.NE.int(LEN(F(1,1)),SIZE_T)) then
             ! existing string type is wrong size - nuke this dataset
             call H5Gunlink_f(group_id,name,hdferr)
             create_new=.TRUE.
          endif
       endif
    else
       create_new=.TRUE.
    endif

    ! Create data type for CHARACTER(len=length) string

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN(F(1,1)),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Open or Create the dataset

    if (create_new) then
       ! Create a new 2D dataspace and dataset
       call h5screate_simple_f(2,dims,space_id,hdferr)                               ! create space_id
       call h5dcreate_f(group_id,name,type_id,space_id,dataset_id,hdferr) ! create dataset_id
    else
       ! Everything cool - just open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    endif

    ! Write to the dataset

    call h5dwrite_f(dataset_id,type_id, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5dump_data_string_2d

  !=============================================================
  ! H5READ_DATA 0D/1D/2D int/double/string wrappers
  ! 
  ! These are Shortcuts to READ 0D/1D/2D data of type int/double/string
  ! from a group or object.
  ! These have a uniform interface:
  !
  !    call h5read_data(group_id,data_value,data_name)
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Read a logical scalar

  subroutine h5read_data_logical_0d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: group_id
    LOGICAL         ,intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER :: Fint ! read an int as there is no bool native type

    call h5read_data(group_id,Fint,name)

    if (Fint==0) then 
       F=.FALSE.
    else
       F=.TRUE.
    end if

    return
  end subroutine h5read_data_logical_0d
  !-----------------------------------------------------------------------
  ! Read a logical scalar

  subroutine h5read_data_logical_1d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: group_id
    LOGICAL,DIMENSION(:),intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER,DIMENSION(size(F)) :: Fint ! read an int as there is no bool native type

    call h5read_data(group_id,Fint,name)

    where (Fint==0) 
       F=.FALSE.
    elsewhere
       F=.TRUE.
    end where

    return
  end subroutine h5read_data_logical_1d
  !-----------------------------------------------------------------------
  ! Read a logical scalar

  subroutine h5read_data_logical_2d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: group_id
    LOGICAL,DIMENSION(:,:),intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER,DIMENSION(size(F,DIM=1),size(F,DIM=2)) :: Fint ! read an int as there is no bool native type

    call h5read_data(group_id,Fint,name)

    where (Fint==0) 
       F=.FALSE.
    elsewhere
       F=.TRUE.
    end where

    return
  end subroutine h5read_data_logical_2d

  !-----------------------------------------------------------------------
  ! Read an integer scalar

  subroutine h5read_data_int_0d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: group_id
    INTEGER         ,intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=1

    ! Open the named dataset and read F

    if (h5dataset_exists(group_id,name)) then
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       if (hdferr == 0) then
          call h5dread_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
          call h5dclose_f(dataset_id,hdferr)
          return
       else
          print*,'FAILED h5d_open_f in h5read_data_int_0d(',group_id,'f,',name,')'
          return
       endif
    else
       ! Dataset not found - Look for an attribute instead

       if (h5attribute_exists(group_id,name)) then
          call h5read_attr(group_id,F,name)
          return
       endif

       print*,'NOT FOUND: H5READ_DATA_INT_0D(',group_id,name,')'
    endif

    return
  end subroutine h5read_data_int_0d

  !-----------------------------------------------------------------------
  ! Read a 1D integer array

  subroutine h5read_data_int_1d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)      ,intent(IN) :: group_id
    INTEGER,DIMENSION(:),intent(OUT):: f
    CHARACTER(len=*)    ,intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_int_1d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_int_1d(',group_id,'f,',name,')'
    endif
    return
  end subroutine h5read_data_int_1d

  !-----------------------------------------------------------------------
  ! Read a 2D integer array

  subroutine h5read_data_int_2d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)        ,intent(IN) :: group_id
    INTEGER,DIMENSION(:,:),intent(OUT):: f
    CHARACTER(len=*)      ,intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(2)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_int_2d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_int_2d(',group_id,'f,',name,')'
    endif
    return
  end subroutine h5read_data_int_2d

  !-----------------------------------------------------------------------
  ! Read a double scalar

  subroutine h5read_data_double_0d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN) :: group_id
    REAL(rn)        ,intent(OUT):: f
    CHARACTER(len=*),intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=1

    ! Open the named dataset in the group and read F

    if (h5dataset_exists(group_id,name)) then
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       if (hdferr == 0) then
          call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
          call h5dclose_f(dataset_id,hdferr)
          return
       else
          print*,'FAILED h5d_open_f in h5read_data_double_0d(',group_id,'f,',name,')'
          return
       endif
    else
       ! Dataset not found - Look for an attribute instead

       if (h5attribute_exists(group_id,name)) then
          call h5read_attr(group_id,F,name)
          return
       endif

       print*,'NOT FOUND: H5READ_DATA_DOUBLE_0D(',group_id,name,')'
    endif

    return
  end subroutine h5read_data_double_0d

  !-----------------------------------------------------------------------
  ! Read a 1D double array

  subroutine h5read_data_double_1d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)       ,intent(IN) :: group_id
    REAL(rn),DIMENSION(:),intent(OUT):: f
    CHARACTER(len=*)     ,intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_double_1d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_double_1d(',group_id,'f,',name,')'
    endif
    return
  end subroutine h5read_data_double_1d

  !-----------------------------------------------------------------------
  ! Read a 2D double array

  subroutine h5read_data_double_2d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)         ,intent(IN) :: group_id
    REAL(rn),DIMENSION(:,:),intent(OUT):: f
    CHARACTER(len=*)       ,intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(2)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_double_2d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_double_2d(',group_id,'f,',name,')'
    endif
    return
  end subroutine h5read_data_double_2d


  !-----------------------------------------------------------------------
  ! Read a 3D double array

  subroutine h5read_data_double_3d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)           ,intent(IN) :: group_id
    REAL(rn),DIMENSION(:,:,:),intent(OUT):: f
    CHARACTER(len=*)         ,intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(3)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_double_2d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_double_3d(',group_id,'f,',name,')'
    endif
    return
  end subroutine h5read_data_double_3d


  !-----------------------------------------------------------------------
  ! Read a 4D double array

  subroutine h5read_data_double_4d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)             ,intent(IN) :: group_id
    REAL(rn),DIMENSION(:,:,:,:),intent(OUT):: f
    CHARACTER(len=*)           ,intent(IN) :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(4)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_double_2d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_double_3d(',group_id,'f,',name,')'
    endif
    return
  end subroutine h5read_data_double_4d


  !-----------------------------------------------------------------------
  ! Read a scalar CHARACTER string

  subroutine h5read_data_string_0d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)  ,intent(IN)   :: group_id
    CHARACTER(len=*),intent(INOUT):: f        ! (needs IN to get LEN(F))
    CHARACTER(len=*),intent(IN)   :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: type_id            ! Data Type identifier

    dims=1

    ! Create a type of LEN(F)*H5T_NATIVE_CHARACTER

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,type_id, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_string_0d(',group_id,'f,',name,')'
    endif

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5read_data_string_0d

  !-----------------------------------------------------------------------
  ! Read a 1D CHARACTER string array

  subroutine h5read_data_string_1d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)   ,intent(IN)               :: group_id
    CHARACTER(len=*),intent(INOUT),dimension(:):: f        ! (needs IN to get LEN(F))
    CHARACTER(len=*) ,intent(IN)               :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code
    INTEGER(HSIZE_T),dimension(1)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: type_id            ! Data Type identifier
    INTEGER,dimension(size(F))::lens    ! Length of each F(i)
    INTEGER :: i

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_string_1d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Check that all F(i) have the same LEN:

    do i=1,size(F)
       lens(i)=LEN(F(i))
    enddo
    if (ANY(lens/=lens(1))) then
       print*,'h5read_data_string_1d - all LEN(F) must be the same'
       print*,lens
       return
    endif

    ! Create a type of LEN(F)*H5T_NATIVE_CHARACTER

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,type_id, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_string_1d(',group_id,'f,',name,')'
    endif

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5read_data_string_1d

  !-----------------------------------------------------------------------
  ! Read a 2D CHARACTER string array

  subroutine h5read_data_string_2d(group_id,f,name)
    ! arguments
    INTEGER(HID_T)   ,intent(IN)                 :: group_id
    CHARACTER(len=*),intent(INOUT),dimension(:,:):: f        ! (needs IN to get LEN(F))
    CHARACTER(len=*) ,intent(IN)                 :: name
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HSIZE_T),dimension(2)::dims ! Dimensions of input
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: type_id            ! Data Type identifier
    INTEGER,dimension(size(F,DIM=1),size(F,DIM=2))::lens    ! Length of each F(i)
    INTEGER :: i,j

    dims=SHAPE(F)
    if(ANY(dims==0)) then
       print*,'h5read_data_string_1d on ',name,' ZERO-SIZED ARRAY ARGUMENT, IGNORING !!'
       print*,dims
       return
    endif

    ! Check that all F(i) have the same LEN:

    do i=1,size(F,DIM=1)
       do j=0,size(F,DIM=2)
          lens(i,j)=LEN(F(i,j))
       enddo
    enddo
    if (ANY(lens/=lens(1,1))) then
       print*,'h5read_data_string_1d - all LEN(F) must be the same'
       print*,lens
       return
    endif

    ! Create a type of LEN(F)*H5T_NATIVE_CHARACTER

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Open the named dataset in the group and read F

    call h5dopen_f(group_id,name,dataset_id,hdferr)
    if (hdferr == 0) then
       call h5dread_f(dataset_id,type_id, F, dims,hdferr)
       call h5dclose_f(dataset_id,hdferr)
    else
       print*,'FAILED h5d_open_f in h5read_data_string_1d(',group_id,'f,',name,')'
    endif

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine h5read_data_string_2d


  !=======================================================================
  ! Old, highly structured printh5 routines.
  ! These used to be necessary to use the IDL GUI.
  ! Now they aren't necessary, but they're still recommended.
  !=======================================================================

  !-----------------------------------------------------------------------
  ! Output a 3D array F to HDF5 file
  !
  ! Optional x,y,title,xtitle,ytitle values are also stored
  ! in the HDF5 file if they are present.
  !
  subroutine printh5_3D_file(file,path,f,x,y,z,title,xtitle,ytitle,ztitle)
    ! arguments
    real(rn),intent(IN),dimension(:,:,:)     :: f
    real(rn),intent(IN),dimension(:),optional:: x,y,z
    CHARACTER(len=*),intent(IN)         :: file
    CHARACTER(len=*),intent(IN)         :: path
    CHARACTER(len=*),intent(IN),optional:: title,xtitle,ytitle,ztitle
    ! for using path:
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path

    call h5open_path(file,TRIM(path),ids,last,name)

    ! printh5 to the group at the leaf at the end of the path

    call printh5_3d_group(ids(last),name,f,x,y,z,title,xtitle,ytitle,ztitle)

    ! close path

    call h5close_path(ids)

    return
  end subroutine printh5_3D_file

  !-----------------------------------------------------------------------
  ! Output a 3D array F to HDF5 group
  !
  ! Optional x,y,title,xtitle,ytitle values are also stored
  ! in the HDF5 file if they are present.
  !
  subroutine printh5_3D_group(parent_id,name,f,x,y,z,title,xtitle,ytitle,ztitle)
    ! arguments
    INTEGER(HID_T),intent(IN)           :: parent_id
    CHARACTER(len=*),intent(IN)         :: name
    real(rn),intent(IN),dimension(:,:,:):: f
    real(rn),intent(IN),dimension(:),optional:: x,y,z
    CHARACTER(len=*),intent(IN),optional:: title,xtitle,ytitle,ztitle
    ! local variables
    integer :: nx,ny,nz
    INTEGER :: hdferr         ! Error code 
    INTEGER(HID_T):: group_id       ! Group identifier
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id     ! Dataset identifier

    nx=size(f,DIM=1)
    ny=size(f,DIM=2)
    nz=size(f,DIM=3)

    ! check that x,y sizes match shape of f:

    if (present(x)) then
       if (size(x)/=nx) then
          print*,'printh5_3D: size(x)=',size(x),'/=size(f,DIM=1)=',nx
          stop
       endif
    endif
    if (present(y)) then
       if (size(y)/=ny) then
          print*,'printh5_3D: size(y)=',size(y),'/=size(f,DIM=2)=',ny
          stop
       endif
    endif
    if (present(z)) then
       if (size(z)/=nz) then
          print*,'printh5_3D: size(z)=',size(z),'/=size(f,DIM=3)=',nz
          stop
       endif
    endif

    ! open/create a group inside the parent_id

    if (h5group_exists(parent_id,name)) then
       call h5gopen_f(parent_id,name,group_id,hdferr)
    else
       call h5gcreate_f(parent_id,name,group_id,hdferr)
    endif

    ! Store array data and metadata

    call h5dump_data( group_id, F,        'data')
    call h5dump_data( group_id, shape(F), 'dims')

    if (present(x)) call h5dump_data( group_id, x, 'X')
    if (present(y)) call h5dump_data( group_id, y, 'Y')
    if (present(z)) call h5dump_data( group_id, z, 'Z')

    ! Store name and title

    call h5dump_data( group_id, TRIM(name), 'name')
    if (present( title)) call h5dump_data( group_id,  title,  'title')
    if (present(xtitle)) call h5dump_data( group_id, xtitle, 'xtitle')
    if (present(ytitle)) call h5dump_data( group_id, ytitle, 'ytitle')
    if (present(ztitle)) call h5dump_data( group_id, ytitle, 'ztitle')

    ! Close group

    CALL h5gclose_f(group_id,hdferr)

    return
  end subroutine printh5_3D_group

  !-----------------------------------------------------------------------
  ! Output a 2D array F to HDF5 file
  !
  ! Optional x,y,title,xtitle,ytitle values are also stored
  ! in the HDF5 file if they are present.
  !
  subroutine printh5_2D_file(file,path,f,x,y,title,xtitle,ytitle)
    ! arguments
    real(rn),intent(IN),dimension(:,:)       :: f
    real(rn),intent(IN),dimension(:),optional:: x
    real(rn),intent(IN),dimension(:),optional:: y
    CHARACTER(len=*),intent(IN)         :: file
    CHARACTER(len=*),intent(IN)         :: path
    CHARACTER(len=*),intent(IN),optional:: title,xtitle,ytitle
    ! for using path:
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path

    call h5open_path(file,TRIM(path),ids,last,name)

    ! printh5 to the group at the leaf at the end of the path

    call printh5_2d_group(ids(last),name,f,x,y,title,xtitle,ytitle)

    ! close path

    call h5close_path(ids)

    return
  end subroutine printh5_2D_file

  !-----------------------------------------------------------------------
  ! Output a 2D array F to HDF5 group
  !
  ! Optional x,y,title,xtitle,ytitle values are also stored
  ! in the HDF5 file if they are present.
  !
  subroutine printh5_2D_group(parent_id,name,f,x,y,title,xtitle,ytitle)
    ! arguments
    INTEGER(HID_T),intent(IN)           :: parent_id
    CHARACTER(len=*),intent(IN)         :: name
    real(rn),intent(IN),dimension(:,:)       :: f
    real(rn),intent(IN),dimension(:),optional:: x
    real(rn),intent(IN),dimension(:),optional:: y
    CHARACTER(len=*),intent(IN),optional:: title,xtitle,ytitle
    ! local variables
    integer :: nx,ny
    INTEGER :: hdferr         ! Error code 
    INTEGER(HID_T):: group_id       ! Group identifier
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id     ! Dataset identifier

    nx=size(f,DIM=1)
    ny=size(f,DIM=2)

    ! check that x,y sizes match shape of f:

    if (present(x)) then
       if (size(x)/=nx) then
          print*,'printh5_2D: size(x)=',size(x),'/=size(f,DIM=1)=',nx
          stop
       endif
    endif
    if (present(y)) then
       if (size(y)/=ny) then
          print*,'printh5_2D: size(y)=',size(y),'/=size(f,DIM=1)=',ny
          stop
       endif
    endif

    ! open/create a group inside the parent_id

    if (h5group_exists(parent_id,name)) then
       call h5gopen_f(parent_id,name,group_id,hdferr)
    else
       call h5gcreate_f(parent_id,name,group_id,hdferr)
    endif

    ! Store array data and metadata

    call h5dump_data( group_id, F,        'data')
    call h5dump_data( group_id, shape(F), 'dims')

    if (present(x)) call h5dump_data( group_id, x, 'X')
    if (present(y)) call h5dump_data( group_id, y, 'Y')

    ! Store name and title

    call h5dump_data( group_id, TRIM(name), 'name')
    if (present( title)) call h5dump_data( group_id,  title,  'title')
    if (present(xtitle)) call h5dump_data( group_id, xtitle, 'xtitle')
    if (present(ytitle)) call h5dump_data( group_id, ytitle, 'ytitle')

    ! Close group

    CALL h5gclose_f(group_id,hdferr)

    return
  end subroutine printh5_2D_group

  !-----------------------------------------------------------------------
  ! Output a 1D array F to HDF5 file
  !
  ! Optional x,y,title,xtitle,ytitle values are also stored
  ! in the HDF5 file if they are present.
  !
  subroutine printh5_1D_file(file,path,f,x,title,xtitle)
    ! arguments
    real(rn),intent(IN),dimension(:)         :: f
    real(rn),intent(IN),dimension(:),optional:: x
    CHARACTER(len=*),intent(IN)         :: file
    CHARACTER(len=*),intent(IN)         :: path
    CHARACTER(len=*),intent(IN),optional:: title,xtitle
    ! for using path:
    INTEGER(HID_T),dimension(10) :: ids                 ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name                ! name of leaf at end of path
    integer :: last

    ! open path

    call h5open_path(file,TRIM(path),ids,last,name)

    ! printh5 to the group at the leaf at the end of the path

    call printh5_1D_group(ids(last),name,f,x,title,xtitle)

    ! close path

    call h5close_path(ids)

    return
  end subroutine printh5_1D_file

  !-----------------------------------------------------------------------
  ! Output a 1D array F to HDF5 group
  !
  ! Optional x,y,title,xtitle,ytitle values are also stored
  ! in the HDF5 file if they are present.
  !
  subroutine printh5_1D_group(parent_id,name,f,x,title,xtitle)
    ! arguments
    INTEGER(HID_T),intent(IN)           :: parent_id
    CHARACTER(len=*),intent(IN)         :: name
    real(rn),intent(IN),dimension(:)         :: f
    real(rn),intent(IN),dimension(:),optional:: x
    CHARACTER(len=*),intent(IN),optional:: title,xtitle
    ! local variables
    integer :: nx
    INTEGER :: hdferr         ! Error code 
    INTEGER(HID_T):: group_id       ! Group identifier
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id     ! Dataset identifier

    nx=size(f)

    ! check that x,y sizes match shape of f:

    if (present(x)) then
       if (size(x)/=nx) then
          print*,'printh5_1D: size(x)=',size(x),'/=size(f,DIM=1)=',nx
          stop
       endif
    endif

    ! open/create a group inside the parent_id

    if (h5group_exists(parent_id,name)) then
       call h5gopen_f(parent_id,name,group_id,hdferr)
    else
       call h5gcreate_f(parent_id,name,group_id,hdferr)
    endif

    ! Store array data and metadata

    call h5dump_data( group_id, F,        'data')
    call h5dump_data( group_id, shape(F), 'dims')

    if (present(x)) call h5dump_data( group_id, x, 'X')

    ! Store name and title

    call h5dump_data( group_id, TRIM(name), 'name')
    if (present( title)) call h5dump_data( group_id,  title,  'title')
    if (present(xtitle)) call h5dump_data( group_id, xtitle, 'xtitle')

    ! Close group

    CALL h5gclose_f(group_id,hdferr)

    return
  end subroutine printh5_1D_group

  !-----------------------------------------------------------------------
  ! Shortcut to store a 2D array in an HDF5 group
  !
  subroutine data_to_group_2D(F,name,group_id)
    ! arguments
    real(rn),intent(IN),dimension(:,:):: f
    CHARACTER(len=*),intent(IN)       :: name
    INTEGER(HID_T),intent(IN)         :: group_id
    ! local variables
    INTEGER :: hdferr         ! Error code 
    INTEGER(HSIZE_T),dimension(2)::dims
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    INTEGER(HID_T):: dataset_id     ! Dataset identifier

    ! Store array data in a dataset in a dataspace

    dims(1)=SIZE(F,DIM=1)
    dims(2)=SIZE(F,DIM=2)
    call h5screate_simple_f(2,dims,space_id,hdferr)                          ! create space_id
    call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine data_to_group_2D

  !-----------------------------------------------------------------------
  ! Shortcut to store a 1D array in an HDF5 group
  !
  subroutine data_to_group_1D(F,name,group_id)
    ! arguments
    real(rn),intent(IN),dimension(:):: f
    CHARACTER(len=*),intent(IN)     :: name
    INTEGER(HID_T),intent(IN)       :: group_id
    ! local variables
    INTEGER :: hdferr         ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    INTEGER(HID_T):: dataset_id     ! Dataset identifier

    ! Store array data in a dataset in a dataspace

    dims(1)=SIZE(F,DIM=1)
    call h5screate_simple_f(1,dims,space_id,hdferr)                          ! create space_id
    call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr) ! create dataset_id
    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine data_to_group_1D

  !-----------------------------------------------------------------------
  ! Shortcut to store a 1D array in an HDF5 group
  !
  subroutine data_to_group_1Dint(F,name,group_id)
    ! arguments
    INTEGER,intent(IN),dimension(:) :: f
    CHARACTER(len=*),intent(IN)     :: name
    INTEGER(HID_T),intent(IN)       :: group_id
    ! local variables
    INTEGER :: hdferr         ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    INTEGER(HID_T):: dataset_id     ! Dataset identifier

    ! Store array data in a dataset in a dataspace

    dims=SHAPE(F)
    call h5screate_simple_f(1,dims,space_id,hdferr)                          ! create space_id
    call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr) ! create dataset_id
    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, dims,hdferr)
    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    return
  end subroutine data_to_group_1Dint

  !-----------------------------------------------------------------------
  ! Add a CHARACTER string attribute to a group
  !
  subroutine data_to_group_string(F,name,group_id)
    ! arguments
    CHARACTER(len=*) ,intent(IN):: f,name
    INTEGER(HID_T),intent(IN)   :: group_id
    ! local variables
    INTEGER :: hdferr         ! Error code 
    INTEGER(HSIZE_T),dimension(1)::dims
    INTEGER(HID_T):: type_id        ! ID of variable length string type
    INTEGER(HID_T):: space_id       ! Dataspace identifier
    !INTEGER(HID_T):: dataset_id     ! Dataset identifier
    INTEGER(HID_T):: attr_id        ! Attribute identifier

    dims=1

    ! Create new type_id variable length string data type for F

    call h5tcopy_f(H5T_NATIVE_CHARACTER,type_id,hdferr)
    call h5tset_size_f(type_id,INT(LEN_TRIM(F),SIZE_T),hdferr)
    call h5tset_strpad_f(type_id,H5T_STR_SPACEPAD_F,hdferr)

    ! Store F in a space in group_id

    call h5screate_f(H5S_SCALAR_F,space_id,hdferr)
    call h5acreate_f(group_id,name,type_id,space_id,attr_id,hdferr)
    call h5awrite_f(attr_id,type_id, TRIM(F), dims,hdferr)
    call h5aclose_f(attr_id,hdferr)
    call h5sclose_f(space_id,hdferr)

    ! Release the datatype

    call h5tclose_f(type_id,hdferr)
    return
  end subroutine data_to_group_string

  !=============================================================
  ! H5APPEND* wrappers
  !
  ! Append a scalar to an array - useful for orbits, time traces, etc
  !
  !=============================================================

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',numout,iout,F,'grand/parent/child/Fvalue')
  !
  subroutine h5append_data_int_0d_auto_fixed(file,fnum,fpos,F,path)
    INTEGER,         intent(IN):: F           ! value to append to array on leaf
    INTEGER,         intent(IN):: fnum,fpos   ! max number of F values to write, position of this value in array
    CHARACTER(len=*),intent(IN):: file,path   ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids       ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name      ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5append_data(ids(last),fnum,fpos,F,name)
    call h5close_path(ids)

    return
  end subroutine h5append_data_int_0d_auto_fixed

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',numout,iout,F,'grand/parent/child/Fvalue')
  !
  subroutine h5append_data_double_0d_auto_fixed(file,fnum,fpos,F,path)
    REAL(kind=rn),   intent(IN):: F           ! value to append to array on leaf
    INTEGER,         intent(IN):: fnum,fpos   ! max number of F values to write, position of this value in array
    CHARACTER(len=*),intent(IN):: file,path   ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids       ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name      ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5append_data(ids(last),fnum,fpos,F,name)
    call h5close_path(ids)

    return
  end subroutine h5append_data_double_0d_auto_fixed

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',numout,iout,F,'grand/parent/child/Fvalue')
  !
  subroutine h5append_data_double_1d_auto_fixed(file,fnum,fpos,F,path)
    REAL(rn),DIMENSION(:),intent(IN):: F           ! value to append to array on leaf
    INTEGER,         intent(IN):: fnum,fpos   ! max number of F values to write, position of this value in array
    CHARACTER(len=*),intent(IN):: file,path   ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids       ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name      ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5append_data(ids(last),fnum,fpos,F,name)
    call h5close_path(ids)

    return
  end subroutine h5append_data_double_1d_auto_fixed

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5append_data_double_0d_auto(file,F,path)
    REAL(kind=rn),   intent(IN):: F         ! value to append to array on leaf
    CHARACTER(len=*),intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids     ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name    ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5append_data(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5append_data_double_0d_auto

  !-----------------------------------------------------------------------
  ! Dump F into file at end of path
  !
  ! example arguments:  ('test.h5',F,'grand/parent/child/Fvalue')
  !
  subroutine h5append_data_int_0d_auto(file,F,path)
    INTEGER,         intent(IN):: F         ! value to append to array on leaf
    CHARACTER(len=*),intent(IN):: file,path ! filename and internal HDF5 path input
    INTEGER(HID_T),dimension(10) :: ids     ! loc_ids of each branch in path
    CHARACTER(len=80)            :: name    ! name of leaf at end of path
    integer :: last

    ! open path
    ! dump the data to the leaf at the end of the path
    ! close path

    call h5open_path(file,path,ids,last,name)
    call h5append_data(ids(last),F,name)
    call h5close_path(ids)

    return
  end subroutine h5append_data_int_0d_auto

  !-----------------------------------------------------------------------
  ! Append scalar data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_double_0d(group_id,f,name)
    ! arguments
    REAL(rn),        intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: cparms             ! Property identifier
    INTEGER(HSIZE_T),DIMENSION(1)::dims,chunksize,maxdims,addsize,start,newdims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size,written_size
    INTEGER:: intsize

    chunksize=100  ! minimum extension chunk
    rank=1        ! will write to a 1D array

    ! Open or Create a CHUNKED dataset in group_id

    maxdims=H5S_UNLIMITED_F
    if (h5dataset_exists(group_id,name)) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial size 0, and "unlimited" maximum extent
       dims=0
       call h5screate_simple_f(rank,dims,space_id,hdferr,maxdims)   ! create space_id
       ! Create a CHUNKED dataset, rank 1, and chunk size chunksize
       call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,hdferr)
       call h5pset_chunk_f(cparms,rank,chunksize,hdferr)
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr,cparms)
       call h5pclose_f(cparms,hdferr)
    endif

    ! get size of existing dataspace, and if necessary,
    ! extend the dataspace with extra chunks

    if (h5attribute_exists(dataset_id,'size')) then
       call h5read_attr(dataset_id,intsize,'size')
       written_size=intsize
    else
       written_size=0
    endif
    call h5dget_space_f(dataset_id,oldspace_id,hdferr)
    call h5sget_simple_extent_npoints_f(oldspace_id,old_size,hdferr)
    call h5sclose_f(oldspace_id,hdferr)

    if (written_size+1 > old_size) then
       newdims=old_size+chunksize
       call h5dextend_f(dataset_id,newdims,hdferr)
    endif

    ! finally, append the data via a new dataspace which starts at the
    ! end of the old data

    addsize=1
    start=written_size
    
    call h5dget_space_f(dataset_id,newspace_id,hdferr) ! get the file space
    call h5sget_simple_extent_npoints_f(newspace_id,new_size,hdferr)

    call h5sselect_hyperslab_f(newspace_id,H5S_SELECT_SET_F,start, &!!! OR npoints+1 (start)
                               addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, addsize,hdferr, memspace_id,newspace_id,H5P_DEFAULT_F)

    ! write the CURRENT SIZE attribute to the dataset as an attribute

    intsize=start(1)+addsize(1)
    call h5dump_attr(dataset_id,intsize,'size')

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(newspace_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_double_0d

  !-----------------------------------------------------------------------
  ! Append scalar data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_double_0d_fixed(group_id,fnum,fpos,f,name)
    ! arguments
    REAL(rn),        intent(IN):: f
    INTEGER,         intent(IN):: fnum,fpos
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HSIZE_T),DIMENSION(1)::dims,maxdims,addsize,start,newdims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size,written_size
    INTEGER:: intsize

    if (fpos < 1.or.fpos > fnum) then
       print*,'h5append_data_double_0d - illegal values of fnum,fpos:',fnum,fpos
       print*,'ignoring'
       return
    endif

    rank=1        ! will write to a 1D array
    maxdims=fnum
    dims=1

    ! Open or Create a dataset in group_id

!!$    if (h5dataset_exists(group_id,name)) then
    if (fpos > 1) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial fnum
       call h5screate_simple_f(rank,maxdims,space_id,hdferr)   ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr)
    endif

    ! insert the data into the dataspace using fpos as the offset

    start  =fpos-1
    addsize=1
    
    call h5sselect_hyperslab_f(space_id,H5S_SELECT_SET_F,start,addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr) !,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, addsize,hdferr, memspace_id,space_id) !,H5P_DEFAULT_F)

    ! write the CURRENT SIZE attribute to the dataset as an attribute

    call h5dump_attr(dataset_id,fpos,'size')

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_double_0d_fixed

  !-----------------------------------------------------------------------
  ! Append 1D data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_double_1d_fixed(group_id,fnum,fpos,f,name)
    ! arguments
    REAL(rn),DIMENSION(:),intent(IN):: f
    INTEGER,         intent(IN):: fnum,fpos
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HSIZE_T),DIMENSION(1)::maxdims,addsize,start,newdims!,dims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size,written_size
    INTEGER:: intsize

    if (fpos < 1.or.fpos > fnum) then
       print*,'h5append_data_double_0d - illegal values of fnum,fpos:',fnum,fpos
       print*,'ignoring'
       return
    endif

    rank=1        ! will write to a 1D array
    maxdims=fnum

    ! Open or Create a dataset in group_id

!!$    if (h5dataset_exists(group_id,name)) then
    if (fpos > 1) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial fnum
       call h5screate_simple_f(rank,maxdims,space_id,hdferr)   ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr)
    endif

    ! insert the data into the dataspace using fpos as the offset

    start  =fpos-1
    addsize=size(f)
    
    call h5sselect_hyperslab_f(space_id,H5S_SELECT_SET_F,start,addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr) !,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, addsize,hdferr, memspace_id,space_id) !,H5P_DEFAULT_F)

    ! write the CURRENT SIZE attribute to the dataset as an attribute

    intsize=fpos+size(f)-1
    call h5dump_attr(dataset_id,intsize,'size')

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_double_1d_fixed

  !-----------------------------------------------------------------------
  ! Append scalar data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_int_0d_fixed(group_id,fnum,fpos,f,name)
    ! arguments
    INTEGER,         intent(IN):: f
    INTEGER,         intent(IN):: fnum,fpos
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HSIZE_T),DIMENSION(1)::dims,maxdims,addsize,start,newdims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size,written_size
    INTEGER:: intsize

    if (fpos < 1.or.fpos > fnum) then
       print*,'h5append_data_double_0d - illegal values of fnum,fpos:',fnum,fpos
       print*,'ignoring'
       return
    endif

    rank=1        ! will write to a 1D array
    maxdims=fnum
    dims=1

    ! Open or Create a dataset in group_id

!!$    if (h5dataset_exists(group_id,name)) then
    if (fpos > 1) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial fnum
       call h5screate_simple_f(rank,maxdims,space_id,hdferr)   ! create space_id
       call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr)
    endif

    ! insert the data into the dataspace using fpos as the offset

    start  =fpos-1
    addsize=1
    
    call h5sselect_hyperslab_f(space_id,H5S_SELECT_SET_F,start,addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr) !,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, addsize,hdferr, memspace_id,space_id) !,H5P_DEFAULT_F)

    ! write the CURRENT SIZE attribute to the dataset as an attribute

    call h5dump_attr(dataset_id,fpos,'size')

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_int_0d_fixed

  !-----------------------------------------------------------------------
  ! Append 1D array data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_double_1d(group_id,f,name)
    ! arguments
    REAL(rn),dimension(:),intent(IN):: f
    CHARACTER(len=*),     intent(IN):: name
    INTEGER(HID_T)  ,     intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: cparms             ! Property identifier
    INTEGER(HSIZE_T),DIMENSION(1)::dims,chunksize,maxdims,addsize,start,newdims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size

    chunksize=20  ! minimum extension chunk
    rank=1        ! will write to a 1D array

    ! Open or Create a CHUNKED dataset in group_id

    maxdims=H5S_UNLIMITED_F
    if (h5dataset_exists(group_id,name)) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial size 0, and "unlimited" maximum extent
       dims=0
       call h5screate_simple_f(rank,dims,space_id,hdferr,maxdims)   ! create space_id
       ! Create a CHUNKED dataset, rank 1, and chunk size chunksize
       call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,hdferr)
       call h5pset_chunk_f(cparms,rank,chunksize,hdferr)
       call h5dcreate_f(group_id,name,H5T_NATIVE_DOUBLE,space_id,dataset_id,hdferr,cparms)
       call h5pclose_f(cparms,hdferr)
    endif

    ! get size of existing dataspace, and if necessary,
    ! extend the dataspace with extra chunks

    call h5dget_space_f(dataset_id,oldspace_id,hdferr)
    call h5sget_simple_extent_npoints_f(oldspace_id,old_size,hdferr)
    newdims=old_size+size(F)
    call h5dextend_f(dataset_id,newdims,hdferr)

    ! finally, append the data via a new dataspace which starts at the
    ! end of the old data

    addsize=size(F)
    start=old_size
    
    call h5dget_space_f(dataset_id,newspace_id,hdferr) ! get the file space
    call h5sget_simple_extent_npoints_f(newspace_id,new_size,hdferr)

    call h5sselect_hyperslab_f(newspace_id,H5S_SELECT_SET_F,start, &!!! OR npoints+1 (start)
                               addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE, F, addsize,hdferr, memspace_id,newspace_id,H5P_DEFAULT_F)

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(oldspace_id,hdferr)
    call h5sclose_f(newspace_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_double_1d

  !-----------------------------------------------------------------------
  ! Append scalar data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_int_0d(group_id,f,name)
    ! arguments
    INTEGER,         intent(IN):: f
    CHARACTER(len=*),intent(IN):: name
    INTEGER(HID_T)  ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: cparms             ! Property identifier
    INTEGER(HSIZE_T),DIMENSION(1)::dims,chunksize,maxdims,addsize,start,newdims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size

    chunksize=20  ! minimum extension chunk
    rank=1        ! will write to a 1D array

    ! Open or Create a CHUNKED dataset in group_id

    maxdims=H5S_UNLIMITED_F
    if (h5dataset_exists(group_id,name)) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial size 0, and "unlimited" maximum extent
       dims=0
       call h5screate_simple_f(rank,dims,space_id,hdferr,maxdims)   ! create space_id
       ! Create a CHUNKED dataset, rank 1, and chunk size chunksize
       call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,hdferr)
       call h5pset_chunk_f(cparms,rank,chunksize,hdferr)
       call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr,cparms)
       call h5pclose_f(cparms,hdferr)
    endif

    ! get size of existing dataspace, and if necessary,
    ! extend the dataspace with extra chunks

    call h5dget_space_f(dataset_id,oldspace_id,hdferr)
    call h5sget_simple_extent_npoints_f(oldspace_id,old_size,hdferr)
    newdims=old_size+1
    call h5dextend_f(dataset_id,newdims,hdferr)

    ! finally, append the data via a new dataspace which starts at the
    ! end of the old data

    addsize=1
    start=old_size
    
    call h5dget_space_f(dataset_id,newspace_id,hdferr) ! get the file space
    call h5sget_simple_extent_npoints_f(newspace_id,new_size,hdferr)

    call h5sselect_hyperslab_f(newspace_id,H5S_SELECT_SET_F,start, &!!! OR npoints+1 (start)
                               addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, addsize,hdferr, memspace_id,newspace_id,H5P_DEFAULT_F)

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(oldspace_id,hdferr)
    call h5sclose_f(newspace_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_int_0d

  !-----------------------------------------------------------------------
  ! Append 1D array data to a 1D sequence.
  ! Create the sequence if it doesn't exist.
  !

  subroutine h5append_data_int_1d(group_id,f,name)
    ! arguments
    INTEGER,dimension(:),intent(IN):: f
    CHARACTER(len=*)    ,intent(IN):: name
    INTEGER(HID_T)      ,intent(IN):: group_id
    ! local variables
    INTEGER :: hdferr                   ! Error code 
    INTEGER(HID_T):: int_id             ! ID of double data type
    INTEGER(HID_T):: space_id,oldspace_id,memspace_id,newspace_id   ! Dataspace identifier
    INTEGER(HID_T):: dataset_id         ! Dataset identifier
    INTEGER(HID_T):: cparms             ! Property identifier
    INTEGER(HSIZE_T),DIMENSION(1)::dims,chunksize,maxdims,addsize,start,newdims
    INTEGER:: rank !,npoints,size,addsize
    INTEGER(HSIZE_T):: size,old_size,new_size

    chunksize=20  ! minimum extension chunk
    rank=1        ! will write to a 1D array

    ! Open or Create a CHUNKED dataset in group_id

    maxdims=H5S_UNLIMITED_F
    if (h5dataset_exists(group_id,name)) then
       ! open the dataset
       call h5dopen_f(group_id,name,dataset_id,hdferr)
       call h5dget_space_f(dataset_id, space_id, hdferr)  ! get dataspace
    else
       ! Create a new 1D dataspace,
       ! with rank 1, initial size 0, and "unlimited" maximum extent
       dims=0
       call h5screate_simple_f(rank,dims,space_id,hdferr,maxdims)   ! create space_id
       ! Create a CHUNKED dataset, rank 1, and chunk size chunksize
       call h5pcreate_f(H5P_DATASET_CREATE_F,cparms,hdferr)
       call h5pset_chunk_f(cparms,rank,chunksize,hdferr)
       call h5dcreate_f(group_id,name,H5T_NATIVE_INTEGER,space_id,dataset_id,hdferr,cparms)
       call h5pclose_f(cparms,hdferr)
    endif

    ! get size of existing dataspace, and if necessary,
    ! extend the dataspace with extra chunks

    call h5dget_space_f(dataset_id,oldspace_id,hdferr)
    call h5sget_simple_extent_npoints_f(oldspace_id,old_size,hdferr)
    newdims=old_size+size(F)
    call h5dextend_f(dataset_id,newdims,hdferr)

    ! finally, append the data via a new dataspace which starts at the
    ! end of the old data

    addsize=size(F)
    start=old_size
    
    call h5dget_space_f(dataset_id,newspace_id,hdferr) ! get the file space
    call h5sget_simple_extent_npoints_f(newspace_id,new_size,hdferr)

    call h5sselect_hyperslab_f(newspace_id,H5S_SELECT_SET_F,start, &!!! OR npoints+1 (start)
                               addsize,hdferr)

    call h5screate_simple_f(rank,addsize,memspace_id,hdferr,maxdims)   ! define memory space to write

    call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER, F, addsize,hdferr, memspace_id,newspace_id,H5P_DEFAULT_F)

    call h5dclose_f(dataset_id,hdferr)
    call h5sclose_f(space_id,hdferr)
    call h5sclose_f(oldspace_id,hdferr)
    call h5sclose_f(newspace_id,hdferr)
    call h5sclose_f(memspace_id,hdferr)

    return
  end subroutine h5append_data_int_1d


end module taehdf5
