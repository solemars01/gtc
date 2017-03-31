subroutine smooth(farray)

  use global_parameters
  use field_array
  use interfaces,only: axisExtrapolate
  use particle_array, only: iload
  use precision
  implicit none

  integer i,j,ii,ij,ip,jt,&
       icount,ierror,idest,isource,isendtag,&
       irecvtag,istatus(MPI_STATUS_SIZE)
  real(lk) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),phism(mgrid),&
       ptemp,pleft(mthetamax),pright(mthetamax),pright0(mthetamax),pright1(mthetamax),&
       pright2(mthetamax),pright3(mthetamax),farraytmp(mgrid),&
       farray(0:1,mgrid),sigma


if(iload/=9)then
! zero out sources in radial boundary cells
    if(nbound>99)then
      sigma=((nbound-101)/2)/(2*sqrt(2*alog(2.0)))
      do i=0,nbound-101
        j=mpsi-i
        if(psi0 .gt. 0.0_lk)then
          farray(1,igrid(i):igrid(i)+mtheta(i))=farray(1,igrid(i):igrid(i)+mtheta(i))*exp(-(((i-(nbound-101))/sigma)**2)/2)
          if(nboundR==0)farray(1,igrid(j):igrid(j)+mtheta(j))=farray(1,igrid(j):igrid(j)+mtheta(j))*exp(-(((i-(nbound-101))/sigma)**2)/2)
        endif
      enddo
      if(nboundR>99)then !right boundary uses different number of points in boundary condition
        sigma=((nboundR-101)/2)/(2*sqrt(2*alog(2.0)))
        do i=0,nboundR-101
          j=mpsi-i
          if(psi0 .gt. 0.0_lk)then
            farray(1,igrid(j):igrid(j)+mtheta(j))=farray(1,igrid(j):igrid(j)+mtheta(j))*exp(-(((i-(nboundR-101))/sigma)**2)/2)
          endif
        enddo
      endif 
    else 
      do i=0,nbound-1
        j=mpsi-i
        if(psi0 .gt. 0.0_lk)then
          farray(1,igrid(i):igrid(i)+mtheta(i))=farray(1,igrid(i):igrid(i)+mtheta(i))*real(i)/real(nbound)
          farray(1,igrid(j):igrid(j)+mtheta(j))=farray(1,igrid(j):igrid(j)+mtheta(j))*real(i)/real(nbound)
       endif
      enddo
    endif
     

! Remove n=0, m/=0 mode for nonlinear sims with all n's
  if(nfilter==0 .and. nonlinear>0) call removn0(farray)
endif


!$omp parallel do private(i)
  do i=1,mgrid
     farraytmp(i)=farray(1,i)
  enddo

  do ip=1,ismooth
#ifdef _FRC
     !$omp parallel do private(i)
     do i=1,mpsi-1
        farraytmp(igrid(i))=farraytmp(igrid(i)+mtheta(i)) ! poloidal BC
     enddo
     ! poloidal smoothing (-0.0625 0.25 0.625 0.25 -0.0625)
     !$omp parallel do private(i,ii,jt,pright,pright0,pright1,pright2,pright3)
     do i=1,mpsi-1
        ii=igrid(i)
        jt=mtheta(i)
        pright(1:jt)=farraytmp(ii+1:ii+jt)
        pright0(1:jt)=cshift(pright(1:jt),-1)
        pright1(1:jt)=cshift(pright(1:jt),1)
        pright2(1:jt)=cshift(pright0(1:jt),-1)
        pright3(1:jt)=cshift(pright1(1:jt),1)
        farraytmp(ii+1:ii+jt)=0.625*pright(1:jt)+0.25*(pright0(1:jt)+pright1(1:jt))-&
             0.0625*(pright2(1:jt)+pright3(1:jt))
     enddo
#else
     if(fem>0.and.eta<machineEpsilon)then
       call smooth_perp_fem(farraytmp)
     else
       call smooth_perp(farraytmp)
     endif
     call smooth_para(farraytmp)
#endif
  enddo

!$omp parallel do private(i)
  do i=1,mgrid
     farray(1,i)=farraytmp(i)
  enddo

  if(bcond==1)call axisExtrapolate(farray(1,:))

  call periodicity(farray)

end subroutine smooth

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine smooth_perp(phitmp)

  use global_parameters
  use field_array
  use particle_array,only:iload
  implicit none

  integer i,j,ii,ij,ip,jt
  real(lk) phism(mgrid),ptemp,phitmp(mgrid),&
       farray(0:1,mgrid)

! radial smoothing
!$omp parallel do private(i)
     do i=1,mpsi-1
        phitmp(igrid(i))=phitmp(igrid(i)+mtheta(i)) ! poloidal BC
     enddo

     phism=0.0
!$omp parallel do private(i,j,ij)
     do i=1,mpsi-1
        do j=1,mtheta(i)
           ij=igrid(i)+j
           phism(ij)=0.25*((1.0-wtp1(1,ij))*phitmp(jtp1(1,ij))+wtp1(1,ij)*phitmp(jtp1(1,ij)+1)+&
                (1.0-wtp1(2,ij))*phitmp(jtp1(2,ij))+wtp1(2,ij)*phitmp(jtp1(2,ij)+1))-&
                0.0625*((1.0-wtp2(1,ij))*phitmp(jtp2(1,ij))+wtp2(1,ij)*phitmp(jtp2(1,ij)+1)+&
                (1.0-wtp2(2,ij))*phitmp(jtp2(2,ij))+wtp2(2,ij)*phitmp(jtp2(2,ij)+1))
        enddo
     enddo

   !todo
   if(.false.)then
!reflective outer boundary
!$omp parallel do private(j,ij)
     do j=1,mtheta(mpsi-1)
         ij=igrid(mpsi-1)+j
       if(iload/=9)then
         phism(ij)=phism(ij)+0.125*((1.0-wtp2(1,ij))*phitmp(jtp2(1,ij))+wtp2(1,ij)*phitmp(jtp2(1,ij)+1))
       else !Non-zero outer boundary
         phism(ij)=phism(ij)+0.125*((1.0-wtp2(1,ij))*phitmp(jtp2(1,ij))+wtp2(1,ij)*phitmp(jtp2(1,ij)+1))-&
                   0.125*((1.0-wtp1(1,ij))*phitmp(jtp1(1,ij))+wtp1(1,ij)*phitmp(jtp1(1,ij)+1))
       endif
     enddo
!reflective inner boundary
!$omp parallel do private(j,ij)
     do j=1,mtheta(1)
         ij=igrid(1)+j
       if(iload/=9)then
         phism(ij)=phism(ij)+0.125*((1.0-wtp2(2,ij))*phitmp(jtp2(2,ij))+wtp2(2,ij)*phitmp(jtp2(2,ij)+1))
       else !Non-zero inner boundary
         phism(ij)=phism(ij)+0.125*((1.0-wtp2(2,ij))*phitmp(jtp2(2,ij))+wtp2(2,ij)*phitmp(jtp2(2,ij)+1))-&
                   0.125*((1.0-wtp1(1,ij))*phitmp(jtp1(2,ij))+wtp1(2,ij)*phitmp(jtp1(2,ij)+1))
       endif
     enddo
   endif

!$omp parallel do private(i,j,ij)
     do i=1,mpsi-1
        do j=1,mtheta(i)
           ij=igrid(i)+j
           phitmp(ij)=0.625*phitmp(ij)+phism(ij)
        enddo
     enddo

end subroutine smooth_perp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine smooth_perp_fem(phitmp)

  use global_parameters
  use field_array
  use particle_array,only:iload
  implicit none

  integer i,j,k,ii,ij,ij_fem,ip,jt
  real(lk) phism(mgrid),ptemp,phitmp(mgrid),&
       smat(mindex_fem,mgrid_fem),rowsum(mgrid_fem),phitmp_fem(mgrid_fem),phism_fem(mgrid_fem)


   phism=0.0

   do i=0,mpsi
      do j=1,mtheta(i)
         ij=igrid(i)+j
         ij_fem=igrid_fem(i)+j
         phitmp_fem(ij_fem)=phitmp(ij)
      enddo
   enddo

      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid_fem(i)+j
            rowsum(ij)=real(nindex_fem(ij)-1) 
         enddo
      enddo 

      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid_fem(i)+j
            do k=1,nindex_fem(ij)
              if(ij==indexp_fem(k,ij))then 
                smat(k,ij)=0.6
              else
                smat(k,ij)=0.4/rowsum(ij)
              endif
            enddo
         enddo
      enddo
     
      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid_fem(i)+j
            phism_fem(ij)=0.0_lk
            do k=1,nindex_fem(ij)
               phism_fem(ij)=phism_fem(ij)+smat(k,ij)*phitmp_fem(indexp_fem(k,ij))
            enddo
         enddo
      enddo 

      do i=0,mpsi
         do j=1,mtheta(i)
            ij=igrid(i)+j
            ij_fem=igrid_fem(i)+j
            phism(ij)=phism_fem(ij_fem)
         enddo
      enddo

!reflective outer boundary
!$omp parallel do private(j,ij)
     do j=1,mtheta(mpsi-1)
         ij=igrid(mpsi-1)+j
       if(iload/=9)then
         phism(ij)=phism(ij)+0.125*((1.0-wtp2(1,ij))*phitmp(jtp2(1,ij))+wtp2(1,ij)*phitmp(jtp2(1,ij)+1))
       else !Non-zero outer boundary
         phism(ij)=phism(ij)+0.125*((1.0-wtp2(1,ij))*phitmp(jtp2(1,ij))+wtp2(1,ij)*phitmp(jtp2(1,ij)+1))-&
                   0.125*((1.0-wtp1(1,ij))*phitmp(jtp1(1,ij))+wtp1(1,ij)*phitmp(jtp1(1,ij)+1))
       endif
     enddo

!$omp parallel do private(i,j,ij)
     do i=1,mpsi-1
        do j=1,mtheta(i)
           ij=igrid(i)+j
           phitmp(ij)=0.625*phitmp(ij)+0.375*phism(ij)
        enddo
     enddo

end subroutine smooth_perp_fem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smooth_para(phitmp)

  use global_parameters
  use field_array
  implicit none

  integer i,j,ii,ij,ip,jt,&
       icount,ierror,idest,isource,isendtag,&
       irecvtag,istatus(MPI_STATUS_SIZE)
  real(lk) sendl(mgrid),sendr(mgrid),recvl(mgrid),recvr(mgrid),phism(mgrid),&
       ptemp,pleft(mthetamax),pright(mthetamax),pright0(mthetamax),pright1(mthetamax),&
       pright2(mthetamax),pright3(mthetamax),phitmp(mgrid)



! parallel smoothing
! send phi to right and receive from left
     sendr=phitmp
     recvl=0.0
     icount=mgrid
     idest=right_pe
     isource=left_pe
     isendtag=myrank_toroidal
     irecvtag=isource
     call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,recvl,icount,&
          mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! send phi to left and receive from right
     sendl=phitmp
     recvr=0.0
     idest=left_pe
     isource=right_pe
     isendtag=myrank_toroidal
     irecvtag=isource
     call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,recvr,icount,&
          mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

!$omp parallel do private(i,ii,jt,j,ij,ptemp,pleft,pright)
     do i=1,mpsi-1
        ii=igrid(i)
        jt=mtheta(i)
        if(myrank_toroidal==0)then !down-shift for zeta=0
           pleft(1:jt)=cshift(recvl(ii+1:ii+jt),-itran(i))
           pright(1:jt)=recvr(ii+1:ii+jt)
        elseif(myrank_toroidal==mtoroidal-1)then !up-shift for zeta=2*pi
           pright(1:jt)=cshift(recvr(ii+1:ii+jt),itran(i))
           pleft(1:jt)=recvl(ii+1:ii+jt)
        else
           pleft(1:jt)=recvl(ii+1:ii+jt)
           pright(1:jt)=recvr(ii+1:ii+jt)
        endif

        do j=1,mtheta(i)
           ij=igrid(i)+j
           ptemp=phitmp(ij)
           phitmp(ij)=0.5*ptemp+0.25*(pleft(j)+pright(j))
        enddo
     enddo

end subroutine smooth_para

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! radial interpolation: find poloidal grid index across flux surfaces & along constant theta
subroutine rinterpolation
  use global_parameters
  use field_array
  use particle_array,only:iload
  implicit none

  integer i,ip,j,indp,indt,ij,jt,ibin,ibout
  real(lk) tdum,wt

! allocate memory
  allocate(jtp1(2,mgrid),jtp2(2,mgrid),wtp1(2,mgrid),wtp2(2,mgrid))

     if(iload==9)then
         ibin=0
         ibout=mpsi
     else
         ibin=1
         ibout=mpsi-1
     endif
!$omp parallel do private(i,ip,j,indp,indt,ij,tdum,jt,wt)
  do i=ibin,ibout
     do ip=1,2
        indp=i+ip
        indt=i-ip

!reflective boundary condition
        if (indp>mpsi) then
            indp=2*mpsi-indp
        endif
        if (indt<0) then
            indt=-indt
        endif

        do j=1,mtheta(i)
           ij=igrid(i)+j
! upward
           if(fielddir==1 .or. fielddir==3)then
             tdum=(real(j)*deltat(i)+(zeta1-pi2)*(qtinv(i)-qtinv(indp)))/deltat(indp)
           else
             tdum=(real(j)*deltat(i)+zeta1*(qtinv(i)-qtinv(indp)))/deltat(indp)
           endif
           jt=floor(tdum)
           wt=tdum-real(jt)
           jt=mod(jt+mtheta(indp),mtheta(indp))
           if(ip==1)then
              wtp1(1,ij)=wt
              jtp1(1,ij)=igrid(indp)+jt
           else
              wtp2(1,ij)=wt
              jtp2(1,ij)=igrid(indp)+jt
           endif
! downward
           if(fielddir==1 .or. fielddir==3)then
             tdum=(real(j)*deltat(i)+(zeta1-pi2)*(qtinv(i)-qtinv(indt)))/deltat(indt)
           else
             tdum=(real(j)*deltat(i)+zeta1*(qtinv(i)-qtinv(indt)))/deltat(indt)
           endif
           jt=floor(tdum)
           wt=tdum-real(jt)
           jt=mod(jt+mtheta(indt),mtheta(indt))
           if(ip==1)then
              wtp1(2,ij)=wt
              jtp1(2,ij)=igrid(indt)+jt
           else
              wtp2(2,ij)=wt
              jtp2(2,ij)=igrid(indt)+jt
           endif
        enddo
     enddo
  enddo
if(iload/=9)then
!!****end points
!! upward
  i=0
!$omp parallel do private(ip,j,indp,ij,tdum,jt,wt)
  do ip=1,2
     indp=min(mpsi,i+ip)
     do j=1,mtheta(i)
           ij=igrid(i)+j
           if(fielddir==1 .or. fielddir==3)then
             tdum=(real(j)*deltat(i)+(zeta1-pi2)*(qtinv(i)-qtinv(indp)))/deltat(indp)
           else
             tdum=(real(j)*deltat(i)+zeta1*(qtinv(i)-qtinv(indp)))/deltat(indp)
           endif
           jt=floor(tdum)
           wt=tdum-real(jt)
           jt=mod(jt+mtheta(indp),mtheta(indp))
           if(ip==1)then
              wtp1(1,ij)=wt
              jtp1(1,ij)=igrid(indp)+jt
           else
              wtp2(1,ij)=wt
              jtp2(1,ij)=igrid(indp)+jt
           endif
     enddo
 enddo
!! downward
  i=mpsi
!$omp parallel do private(ip,j,indt,ij,tdum,jt,wt)
  do ip=1,2
     indt=max(0,i-ip)
     do j=1,mtheta(i)
           ij=igrid(i)+j
           if(fielddir==1 .or. fielddir==3)then
             tdum=(real(j)*deltat(i)+(zeta1-pi2)*(qtinv(i)-qtinv(indt)))/deltat(indt)
           else
             tdum=(real(j)*deltat(i)+zeta1*(qtinv(i)-qtinv(indt)))/deltat(indt)
           endif
           jt=floor(tdum)
           wt=tdum-real(jt)
           jt=mod(jt+mtheta(indt),mtheta(indt))
           if(ip==1)then
              wtp1(2,ij)=wt
              jtp1(2,ij)=igrid(indt)+jt
           else
              wtp2(2,ij)=wt
              jtp2(2,ij)=igrid(indt)+jt
           endif
     enddo
  enddo
endif
end subroutine rinterpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! d_phieff/d_t for fluid-kinetic hybrid electron model
subroutine dphieff
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer ij,j,i,jt
  real(lk) deltime

  deltime=max(1.0e-10,tstep)
  jt=irk+2*(ihybrid-1)

  if (magnetic==0) then
  !$omp parallel do private(ij)
    do ij=1,mgrid
      !corresponding ihybrid method for dphi/dt
       phit(:,ij)=(phi(:,ij)-phisave(:,ij,jt))/deltime   
       phisave(:,ij,jt)=phi(:,ij) !phisave is stored per ihybrid iteration
      !converged nhybrid method for dphi/dt
!       phit(:,ij)=(phi(:,ij)-phisave(:,ij,irk))/deltime
!       if(ihybrid==nhybrid) phisave(:,ij,irk)=phi(:,ij) !phisave is stored at ihybrid = nhybrid, needed for etrap=2
    enddo
    if(istep==1 .and. mype==0 .and. irk==1 .and. ihybrid==nhybrid)write(gtcout,*)'Using corresponding dphi/dt'
!    if(istep==1 .and. mype==0 .and. irk==1 .and. ihybrid==nhybrid)write(gtcout,*)'Using converged dphi/dt'
  else
!$omp parallel do private(i,j,ij)
    do i=0,mpsi
       do j=0,mtheta(i)
          ij=igrid(i)+j
         !corresponding ihybrid method for dne/dt
          if(ihybrid==1)then
            dnet(:,ij)=(sfluidne(:,ij)-dnesave(:,ij,jt))/deltime
            dnesave(:,ij,jt)=sfluidne(:,ij)
            !dnet(:,ij)=(sfluidne(:,ij)-dnesave(:,ij,irk))/deltime
          else
            dnet(:,ij)=(sfluidne(:,ij)-densitye(:,ij)*meshne(i)-dnesave(:,ij,jt))/deltime
            dnesave(:,ij,jt)=sfluidne(:,ij)-densitye(:,ij)*meshne(i) !dnesave is stored per ihybrid iteration
       
           !converged nhybrid method for dne/dt
           !dnesave is stored at ihybrid = nhybrid, needed for etrap=2
            !dnet(:,ij)=(sfluidne(:,ij)-dnesave(:,ij,irk))/deltime-&
            !           (densitye(:,ij)-dnesave(:,ij,irk+2))*meshne(i)/deltime
            !if(ihybrid==nhybrid)then
            !    dnesave(:,ij,irk)=sfluidne(:,ij)
            !    dnesave(:,ij,irk+2)=densitye(:,ij)
            !endif
          endif
       enddo
    enddo
 endif

end subroutine dphieff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine removn0(farray)
  use global_parameters
  use field_array,only:deltat,qtinv,zeta1,mtheta,igrid,nshift,wtshift
  implicit none
  integer i,j,ij,ierror
  real(lk) farray(0:1,mgrid),phitmp(mgrid),wtmp

!! map alpha grid to theta grid
!$omp parallel do private(i,wtmp,j,ij)
  do i=0,mpsi
     phitmp(igrid(i)+1:igrid(i)+mtheta(i))=cshift(farray(1,igrid(i)+1:igrid(i)+mtheta(i)),-nshift(i))
     phitmp(igrid(i))=phitmp(igrid(i)+mtheta(i))
     wtmp=wtshift(i)
! linear interpolation to theta grid
     do j=1,mtheta(i)
        ij=igrid(i)+j
        farray(1,ij)=phitmp(ij-1)*wtmp+phitmp(ij)*(1.0_lk-wtmp)
     enddo
     farray(1,igrid(i))=farray(1,igrid(i)+mtheta(i))
  enddo

!! remove n=0 m\=0 mode
  call MPI_ALLREDUCE(farray(1,:),phitmp,mgrid,mpi_Rsize,MPI_SUM,toroidal_comm,ierror)
  farray(1,:)=farray(1,:)-phitmp/real(mtoroidal)

!! map theta grid back to alpha grid
!$omp parallel do private(i,wtmp,j,ij)
  do i=0,mpsi
     wtmp=wtshift(i)
     do j=0,mtheta(i)-1
        ij=igrid(i)+j
        phitmp(ij)=farray(1,ij+1)*wtmp+farray(1,ij)*(1.0_lk-wtmp)
     enddo
     phitmp(igrid(i)+mtheta(i))=phitmp(igrid(i))

     farray(1,igrid(i)+1:igrid(i)+mtheta(i))=cshift(phitmp(igrid(i)+1:igrid(i)+mtheta(i)),nshift(i))
     farray(1,igrid(i))=farray(1,igrid(i)+mtheta(i))
  enddo

end subroutine removn0


subroutine axisExtrapolate(farray)
  use precision
  use global_parameters,only: mgrid,rg0,mpsilow
  use field_array,only: mtheta,igrid,deltar,deltat,zeta1,qtinv
  implicit none

  integer :: i,j,ij
  real(lk),dimension(mgrid) :: farray
  real(lk) ave,cosave,sinave,tdum,rdum

  !Extrapolate to the magnetic axis.
  ave=0.0
  cosave=0.0
  sinave=0.0

  do i=mpsilow+1, mpsilow+2
    rdum=(rg0+deltar*real(mpsilow+1))/(real(mtheta(i))*(rg0+deltar*real(i)))
    do j=1, mtheta(i)
      ij=igrid(i)+j
      tdum=deltat(i)*real(j)+zeta1*qtinv(i)
      ave=ave+farray(ij)*rdum
      cosave=cosave+2.0*farray(ij)*cos(tdum)*rdum
      sinave=sinave+2.0*farray(ij)*sin(tdum)*rdum
    enddo
  enddo

  sinave=0.5_lk*sinave
  cosave=0.5_lk*cosave
  ave=0.5_lk*ave

!Exterpolation mpsilow flux surface to origin
  do i=0, mpsilow
    rdum=(rg0+deltar*real(i))/(rg0+deltar*real(mpsilow+1))
    do j=1, mtheta(i)
      ij=igrid(i)+j
      tdum=deltat(i)*real(j)+(zeta1)*qtinv(i)
      farray(ij)=ave*rdum !!!
      farray(ij)=farray(ij)+cosave*cos(tdum)*rdum
      farray(ij)=farray(ij)+sinave*sin(tdum)*rdum
    enddo
  enddo

!theta boundary condition
  do i=0, mpsilow+1
      farray(igrid(i))=farray(igrid(i)+mtheta(i))
  enddo
end subroutine axisExtrapolate
