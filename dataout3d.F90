Subroutine dataout3d
  use global_parameters
  use field_array
  implicit none
#if ADIOS
  #define ADIOS_WRITE_PATH(a,b,c) call adios_write_path(a,'b'//char(0),b,c//char(0))
  #define ADIOS_WRITE(a,b) call adios_write(a,'b'//char(0),b)
  integer*8 group_handle,type_id
  integer next,previous,current
#else
  include 'netcdf.inc'
#endif

  integer i,j,n,ij,istatus,ncid,dataid2(3),dimid(2),dataid1(2)
  real(lk) data3d(mgrid),radial(mpsi+1),&
       r,theta,zeta,theta0
  character(len=1) cdum(3)
  character(len=9) vdum
  character(len=100) fdum
! HDF5 declarations

!!first of all, write down the dimension information

#ifdef ADIOS
! Since we have several MPI processes on each plane, only one of these will
! participate in writing out the data for that plane. Let's pick the processes
! with myrank_partd=0 on each plane.
  if(myrank_partd==0)then
     if(myrank_toroidal==0)then
        previous=-1
     else
        previous=myrank_toroidal-1
     endif
     current = myrank_toroidal
     if(myrank_toroidal==(nproc_toroidal-1))then
        next=-1
     else
        next=myrank_toroidal+1
     endif
     if(istep==ndiag)then
        fdum='phi_dir/RUNcoord.bp'//char(0)
        write(dirstr,'("/Coordinate_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
        dirstr=trim(dirstr)//char(0)
        
        call adios_get_group(type_id,"output3d.0"//char(0))
        call adios_set_path (type_id,dirstr)
        call adios_open(group_handle,type_id,fdum)
! call adios_group_by(type_id,"mype"//char(0),toroidal_comm,previous,current,next)
        do i=0,mpsi
           radial(i+1)=(rg0+deltar*real(i))/rho0
        enddo
! Mapping of magnetic coordinates onto cartesian grid and write into file
! We explicitely add the zeta=0 poloidal plane since it will be used by
! the vizualization software to connect the torus. myrank_toroidal=0 takes care of
! this, having one extra plane to write into the output file (kp=1).
        cdum(1)='X'
        cdum(2)='Y'
        cdum(3)='Z'
!! if(myrank_toroidal.eq.0)then
!!   kp=1   !add the plane zeta=0 to myrank_toroidal=0 (zeta0=0 for myrank_toroidal=0)
!! else
!!   kp=0   !no changes for myrank_toroidal > 0
!! endif
        kp=0
        do k=1,1+kp
           zeta=zeta0+deltaz*real(k-kp)
           do i=0,mpsi
              r=rg0+deltar*real(i)
              theta_start(i+1,k)=zeta*qtinv(i)
              do j=0,mtheta(i)
                 ij=igrid(i)+j
                 theta=deltat(i)*real(j)+zeta*qtinv(i)
                 theta0=theta
                 X(ij,k)=cos(zeta)*(1.0+r*cos(theta0))
                 Y(ij,k)=-sin(zeta)*(1.0+r*cos(theta0))
                 Z(ij,k)=r*sin(theta0)
              enddo
           enddo
        enddo

!!write dimension
        ADIOS_WRITE(group_handle,toroidal_comm)
        ADIOS_WRITE(group_handle,mpsi+1)
        ADIOS_WRITE(group_handle,kp+1)
!!write data
        ADIOS_WRITE(group_handle,myrank_toroidal)
        ADIOS_WRITE(group_handle,mpsi)
        ADIOS_WRITE(group_handle,kp)
        ADIOS_WRITE(group_handle,nproc_toroidal)
        ADIOS_WRITE(group_handle,radial)
        call adios_write(group_handle,"mtheta"//char(0),(mtheta+1))
        ADIOS_WRITE(group_handle,itran)
        ADIOS_WRITE(group_handle,mgrid)
        ADIOS_WRITE(group_handle,X)
        ADIOS_WRITE(group_handle,Y)
        ADIOS_WRITE(group_handle,Z)
        call adios_close(group_handle)
     endif !! end if if(istep==ndiag)

! potential data file
     kp=0
     do k=1,1+kp
        do i=0,mpsi
           do j=0,mtheta(i)
              ij=igrid(i)+j
              dataout(ij,k)=phi(k-kp,ij)
           enddo
        enddo
     enddo
     
     write(output_fname,'("phi_dir/PHI_",i5.5,".bp")')mstepall+istep
     output_fname=trim(output_fname)//char(0)
     write(dirstr,'("/Potential_tor",i5.5)')myrank_toroidal !!!(mstepall+istep)/ndiag
     dirstr=trim(dirstr)//char(0)
     call adios_get_group (type_id, "output3d.1"//char(0))
     call adios_set_path (type_id,dirstr)
     call adios_open (group_handle, type_id, output_fname)
     ADIOS_WRITE(group_handle,toroidal_comm)
     ADIOS_WRITE(group_handle,kp+1)  ! TODO check this for correctness

     ADIOS_WRITE(group_handle,myrank_toroidal)
     ADIOS_WRITE(group_handle,nproc_toroidal)
     ADIOS_WRITE(group_handle,mgrid)
!!ADIOS_WRITE(group_handle,mzeta)
     ADIOS_WRITE(group_handle,kp)
     ADIOS_WRITE(group_handle,mpsi)
     ADIOS_WRITE(group_handle,mpsi+1)
!!ADIOS_WRITE(group_handle,mype)
     call adios_write(group_handle,"mtheta"//char(0),(mtheta+1))
     call adios_write (group_handle, "phi"//char(0), dataout)
     call adios_close (group_handle)

  endif !!end myrank_partd

! if(myrank_partd==0)then
! The processes not participating in writing to file wait at the following
! barrier. We use the partd_comm communicator since it does not require a
! fully global synchronization.
  call MPI_BARRIER(partd_comm,ierr)

#else

!!!netcdf 3d data
  if(istep==ndiag)then
     if(mype==0)then
        do i=0,mpsi
           radial(i+1)=(rg0+deltar*real(i))/rho0
        enddo
! write run-dimension
        fdum='phi_dir/RUNdimen.ncd'
! open netcdf data file
        istatus=nf_create(trim(fdum),nf_clobber,ncid)

! define data array dimension
        istatus=nf_def_dim(ncid,'scalar',1,dimid(1))
        istatus=nf_def_dim(ncid,'flux-surfaces',mpsi+1,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,'PE-number',nf_int,1,dimid(1),dataid1(1))
        istatus=nf_def_var(ncid,'flux-surface-number',nf_int,1,dimid(1),dataid1(2))
        istatus=nf_def_var(ncid,'radial-grids',nf_real,1,dimid(2),dataid2(1))
        istatus=nf_def_var(ncid,'poloidal-grids',nf_int,1,dimid(2),dataid2(2))
        istatus=nf_def_var(ncid,'index-shift',nf_int,1,dimid(2),dataid2(3))
! end of netcdf definition
        istatus=nf_enddef(ncid)

! write data
        istatus=nf_put_var_int(ncid,dataid1(1),numberpe)
        istatus=nf_put_var_int(ncid,dataid1(2),mpsi+1)
        istatus=nf_put_var_real(ncid,dataid2(1),radial)
        istatus=nf_put_var_int(ncid,dataid2(2),mtheta+1)
        istatus=nf_put_var_int(ncid,dataid2(3),itran)
! check error
        if (istatus .ne. nf_noerr) then
           print *, nf_strerror(istatus)
        endif
! close netcdf data file
        istatus=nf_close(ncid)
     endif

! grid coordinates,each PE writes to a separate file 
     cdum(1)='X'
     cdum(2)='Y'
     cdum(3)='Z'
     do n=1,3
        zeta=zeta1
        do i=0,mpsi
           r=rg0+deltar*real(i)
           do j=0,mtheta(i)
              ij=igrid(i)+j
              theta=deltat(i)*real(j)+zeta*qtinv(i)
              theta0=theta
! grid coordinates (x,y,z), use geometry center as the origin 
              if(n==1)then
                 data3d(ij)=cos(zeta)*(1.0+r*cos(theta0))
              elseif(n==2)then
                 data3d(ij)=-sin(zeta)*(1.0+r*cos(theta0))
              else
                 data3d(ij)=r*sin(theta0)
              endif
           enddo
        enddo
        
! coordinate data file name
        if(myrank_toroidal < 10)then
           write(fdum,'("phi_dir/NCD",a1,"coor.00",i1)')cdum(n),myrank_toroidal
        elseif(myrank_toroidal < 100)then
           write(fdum,'("phi_dir/NCD",a1,"coor.0",i2)')cdum(n),myrank_toroidal
        else
           write(fdum,'("phi_dir/NCD",a1,"coor.",i3)')cdum(n),myrank_toroidal
        endif
        
! variable name
        write(vdum,'(a1,"coordina")')cdum(n)

! open netcdf data file
        istatus=nf_create(trim(fdum),nf_clobber,ncid)
! define data array dimension
        istatus=nf_def_dim(ncid,'poloidal_grid',mgrid,dimid(1))
        istatus=nf_def_dim(ncid,'toroidal_grid',1,dimid(2))
! define data array id
        istatus=nf_def_var(ncid,vdum,nf_real,2,dimid,dataid1(1))
! end of netcdf definition
        istatus=nf_enddef(ncid)
! write data
        istatus=nf_put_var_real(ncid,dataid1(1),data3d)
        istatus=nf_close(ncid)        
     enddo
  endif

! potential data file
  if(myrank_partd==0)then
     do i=0,mpsi
        do j=0,mtheta(i)
           ij=igrid(i)+j
           data3d(ij)=phi(1,ij)
        enddo
     enddo
     
     write(fdum,'("phi_dir/PHI_",i0,"_",i0,".ncd")')(mstepall+istep),myrank_toroidal

   ! variable name
     write(vdum,'("Potential")')

   ! open netcdf data file
     istatus=nf_create(trim(fdum),nf_clobber,ncid)
   ! define data array dimension
     istatus=nf_def_dim(ncid,'poloidal_grid',mgrid,dimid(1))
     istatus=nf_def_dim(ncid,'toroidal_grid',1,dimid(2))
   ! define data array id
     istatus=nf_def_var(ncid,vdum,nf_real,2,dimid,dataid1(1))
   ! end of netcdf definition
     istatus=nf_enddef(ncid)
   ! write data
     istatus=nf_put_var_real(ncid,dataid1(1),data3d)
     istatus=nf_close(ncid)
  endif
  
#endif
      
end subroutine dataout3d
