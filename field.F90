! Solve fields: phi_eff, phi_ind, fluidue
subroutine field_solver
  use global_parameters
  use particle_array
  use field_array
  use equilibrium, only: psiw
  implicit none

  integer i,j,ij,n,m,win0,win1,intantenna,ierror
  real(lk) gperp2,g2_inv,adum,wa,envelope,tdum,wt,qmin,d2fluidnetmp(0:1,1:mgrid),&
           d4fluidne1(mgrid),d4fluidne2(mgrid),c1,c2
  real(lk),dimension(mgrid) :: phi_ext

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
  
  d2fluidnetmp=0.0_lk
  d4fluidne=0.0_lk
  d4fluidne1=0.0_lk
  d4fluidne2=0.0_lk
! faked k_perp for benchmarking petsc
  gperp2=0.1_lk !k_perp rho0 value
  gperp2=(gperp2/rho0)**2 !k_perp^2 in basic unit
  g2_inv=1.0_lk/gperp2

#ifdef _PETSc
  if((istep==1).and.(mype==0).and.(irk==1))write(gtcout,*)'ILAPLACIAN=',ilaplacian
  if(ilaplacian==0) then
     call laplacian_petsc(sdelapara(1,:),d2apara) !use Poisson Eq's metric, valid for k_perp rho <<1
  elseif(ilaplacian==1 .and. fem==0)then
     call laplacian(sdelapara(1,:),d2apara)! finite-difference for Laplacian
  elseif(ilaplacian==1 .and. fem>0)then 
     call laplacian_fem(1,sdelapara(1,:),d2apara)!
  endif
#else
  if(mype==0)then
    write(gtcout,*) "Please re-compile GTC with PETSc to do electromagnetic simulations."
    call flush(gtcout)	
  endif
  call MPI_FINALIZE(ierror)
  call exit(1)
#endif

! calculate 4-th order derivative of fluid density in parallel direction for hyperviscosity
! hyperviscosity coefficients (c1 for parallel and c2 for perpendicular)
  c1=-hypr1/(mtoroidal*mtoroidal*mtoroidal*mtoroidal*tstep)
  c2=-hypr2/(mthetamax*mthetamax*mthetamax*mthetamax*tstep)
  if (abs(c1)+abs(c2)>1.0e-10) then
   d4fluidne1(:)=sfluidne(1,:)
   d4fluidne2(:)=sfluidne(1,:)
    do n = 1, 2
      call grad2para(d4fluidne1,d2fluidnetmp(:,:))
      d4fluidne1(:)=d2fluidnetmp(1,:)
      call grad2psi(d4fluidne2(:),d2fluidnetmp(:,:))
      d4fluidne2(:)=d2fluidnetmp(1,:)
    enddo
!$omp parallel do private(i,j,ij)
    do i=0,mpsi
     do j=0,mtheta(i)   ! j begins from 0 so that we don't need poloidal BC
        ij=igrid(i)+j
        d4fluidne(ij)=c1*d4fluidne1(ij)+c2*d4fluidne2(ij)*1.0e-12
     enddo
    enddo
 endif

! intantenna=1: launch external anntena, use standing wave in cos(m*theta-n*zeta)*cos(omega*t)
! if using intantenna, must set your own n & m & wa
  intantenna=0
  phi_ext=0.0_lk
  if(intantenna==1)then
     n=0
     m=0

! radial domain applying external source
!by zhs use old intantenna profile
     win0=int(mpsi/4)
     win1=int(3*mpsi/4)
     adum=pi/real(win1-win0)
     wa=0.0_lk

     if(mype==0 .and. istep==1 .and. irk==1)then
        write(gtcout,*)'Antenna freq=',wa,'=',wa/(pi2*utime),'Hz'
        call FLUSH(gtcout)
     endif

! phase of driving source
     wt=cos(wa*tstep*(real(istep+mstep*irun)+0.5*irk)) !irun = # of restart runs
!     do i=win0,win1
!        envelope=1.0e-10*sin(adum*real(i-win0))
!~~~~~~~~~~~new smooth profile~~~~~~~~~~~~~~~~~~~~~
     do i=0,mpsi
         adum=real(i)/real(mpsi)
!         adum=20.0*(adum-0.5)
!         adum=adum*adum*adum*adum
         adum=20.0_lk*(adum-0.5_lk)
         adum=adum*adum
         envelope=1.0e-10*exp(-adum)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!$omp parallel do private(j,ij,tdum)
        do j=1,mtheta(i)
           ij=igrid(i)+j
           if(fielddir==1 .or. fielddir==3)then
             tdum=deltat(i)*real(j)+(zeta1-pi2)*qtinv(i)
           else
             tdum=deltat(i)*real(j)+zeta1*qtinv(i)
           endif
           phi_ext(ij)=envelope*cos(real(m)*tdum-real(n)*zeta1)*wt
        enddo
     enddo
     call gradient(phi_ext(:),gradext)

! write out driving source
     if(mype==0 .and. istep==1 .and. irk==1)then
        open(9991,file='time_ext.out',status='replace')
        open(9992,file='phi_ext.out',status='replace')
     endif
     if(mype==0 .and. irk==1)then
        write(9991,*)wt
        write(9992,*)phi_ext(igrid(mpsi/2)+1)
     endif
     if(mype==0 .and. istep==mstep .and. irk==1)then
        close(9991)
        close(9992)
     endif
  endif

! solver fluid_ue in EM simulation
!$omp parallel do private(i,j,ij)
  do i=0,mpsi
     do j=1,mtheta(i)
        ij=igrid(i)+j

! Electrostatic potential; , phi=[ B_0(R)/B_0(R_0)]^2 (delta_ne)
! Reduced MHD: k_perp rho_s <<1; GK Poisson Eq. becomes Laplace's Eq.
!        phi(1,ij)=-sfluidne(1,ij)*g2_inv*b*b !use constant k_perp; overide petsc solution

! Electron parallel flow from Ampere's law (Holod09: Eq. 11)
! d2apara is the perpendicular laplacian in GTC unit
        fluidue(1,ij)=2.0_lk*d2apara(ij)*rho0*rho0/(betae*bmesh(ij))+&
             meshni(i)*qion*flowi(1,ij)/bmesh(ij)

        if(fload>0)then
           fluidue(1,ij)=fluidue(1,ij)+meshnf(i)*qfast*flowf(1,ij)/bmesh(ij)
        endif

         if(feload>0)then
           fluidue(1,ij)=fluidue(1,ij)+meshnfe(i)*qfaste*flowfe(1,ij)/bmesh(ij)
        endif

! use constant k_perp in Ampere's law for benchmarking Laplacian operator
!        fluidue(1,ij)=(-gperp2*sapara(1,ij)*rho0*rho0/betae)/bmesh(ij)
     enddo
  enddo

! use one or a few toroidal modes
  if(nfilter>0)then
     CALL FILTER(fluidue)
  endif

! smooth potential
  CALL SMOOTH(fluidue)

end subroutine field_solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine field_gradient
  use global_parameters
  use field_array
  use particle_array
  use magnetic_island
  use equilibrium,only: mesher 
  implicit none

  integer i,j,ij,ibin,ibout
  real(lk) rg,q,b_inv,q_inv,ne,te,dnedp,dtedp,deffdp,deffdt,deffdz,ddnedp,ddnedt,ddnedz,dlamdt,dlamdz,&
           dpsidp,dpsidt,dpsidz,ddpdp,ddpdt,ddpdz,phipEq
! field gradients
  if(magnetic==1)then
     call gradient(fluidue(1,:),gradue)
     call gradient(sapara(1,:),gradapara)
     call gradient(sfluidne(1,:),gradne)
     call gradient(sdeltapsi(1,:),gradpsi)
     if (nhybrid>0) then
       call gradient(densitye(1,:),gradkne)
       call gradient(pressureepara(1,:),gradpepara)
       call gradient(pressureeperp(1,:),gradpeperp)
     endif


     if(iload==0 .and. eta .le. 0.0_lk)then
!       ideal MHD 
        gradphieff(:,1,ij)=0.0_lk
     else
!       gyrokinetic/tearing
! gradient of phi_eff
!$omp parallel do private(i,j,ij,b_inv,q_inv,ne,te,dnedp,dtedp,deffdp,deffdt,deffdz,&
!$omp& ddnedp,ddnedt,ddnedz,dlamdt,dlamdz,dpsidp,dpsidt,dpsidz,ddpdp,ddpdt,ddpdz)
        do i=0,mpsi
           q_inv=1.0_lk/qmesh(i)
           if(qmesh(i)==0.0_lk)q_inv=0.0_lk

           do j=1,mtheta(i)
              ij=igrid(i)+j
              b_inv=1.0_lk/bmesh(ij)

              ddnedp=gradne(1,1,ij)
              ddnedt=gradne(2,1,ij)
              ddnedz=gradne(3,1,ij)-gradne(2,1,ij)*q_inv

              ne=meshne(i)
              te=meshte(i)
              phipEq=-mesher(i)
              dnedp = -ne*kapane(i)
              dtedp = -te*kapate(i)

              dpsidp=gradpsi(1,1,ij)
              dpsidt=gradpsi(2,1,ij)
              dpsidz=gradpsi(3,1,ij)-gradpsi(2,1,ij)*q_inv

              deffdp=-(ddnedp-dpsidp*dnedp)*qelectron*te/ne+dpsidp*phipEq
              deffdt=-(ddnedt-dpsidt*dnedp)*qelectron*te/ne+dpsidt*phipEq
              deffdz=-(ddnedz-dpsidz*dnedp)*qelectron*te/ne+dpsidz*phipEq

              gradphieff(1,1,ij)=deffdp+gradkne(1,1,ij)*qelectron*te
              gradphieff(2,1,ij)=deffdt+gradkne(2,1,ij)*qelectron*te
              gradphieff(3,1,ij)=deffdz+deffdt*q_inv+gradkne(3,1,ij)*qelectron*te+&
                                 eta*2.0_lk*d2apara(ij)*rho0*rho0/betae
           enddo
        enddo
        if(nfilter>0)then
           CALL FILTER(gradphieff(1,:,:))
           CALL FILTER(gradphieff(2,:,:))
           CALL FILTER(gradphieff(3,:,:))
        else
           call periodicity(gradphieff(1,:,:))
           call periodicity(gradphieff(2,:,:))
           call periodicity(gradphieff(3,:,:))
        endif
     endif
  endif

! grad phi_electrostatic
  call gradient(phi(1,:),gradphi)

! grad alpha_island
  if(island==1)then
    call gradient(alphaIs(1,:),gradalphaIs)
  endif

! add (0,0) mode to d phi/d psi
  if(izonal==0)then
     phip00=0.0_lk
     phi00=0.0_lk
  elseif(izonal==1)then
! solve zonal flow use new method 
     call zonal  !zonal flow solver for general geometry

     if(iload==9)then
         ibin=0
         ibout=mpsi
     else
         ibin=1
         ibout=mpsi-1
     endif
!$omp parallel do private(i)
     do i=ibin,ibout
        gradphi(1,0:1,igrid(i):igrid(i)+mtheta(i))=gradphi(1,0:1,igrid(i):igrid(i)+mtheta(i))+&
             phip00(i)
     enddo

  else
  ! solve zonal flow
     call zonal_old    !old verson solver
     !$omp parallel do private(i,rg,q)
     do i=1,mpsi-1
        rg=rg0+deltar*real(i)
        q=qmesh(i)
        gradphi(1,0:1,igrid(i):igrid(i)+mtheta(i))=gradphi(1,0:1,igrid(i):igrid(i)+mtheta(i))+&
           phip00(i)*q/rg
     enddo
  endif
end subroutine field_gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! old version of zonal flow solver
subroutine zonal_old
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i
  real(lk) den00(0:mpsi),rg

! zonal density
  phip00=qion*zonali*meshni+qelectron*fluidne00*real(magnetic)
  if(nhybrid>0)phip00=phip00+qelectron*zonale*meshne
  if(fload>0)phip00=phip00+qfast*zonalf*meshnf
  if(feload>0)phip00=phip00+qfaste*zonalfe*meshnfe

! (-0.0625 0.25 0.625 0.25 -0.0625) radial smoothing of (0,0) mode density phip00
  do i=1,2
     den00(0)=phip00(0)
     den00(mpsi)=phip00(mpsi)
     den00(1)=phip00(3)
     den00(mpsi-1)=phip00(mpsi-3)
     den00(2:mpsi-2)=phip00(0:mpsi-4)+phip00(4:mpsi)
     den00(1:mpsi-1)=0.625_lk*phip00(1:mpsi-1)+0.25_lk*(phip00(0:mpsi-2)+phip00(2:mpsi))-&
          0.0625_lk*den00(1:mpsi-1)
     phip00=den00
  enddo

! phi00=r*E_r, E_r(a0)=0. Trapezoid rule
  den00=phip00
  phip00=0.0_lk
  do i=1,mpsi
     rg=rg0+deltar*real(i)
     phip00(i)=phip00(i-1)+0.5_lk*deltar*((rg-deltar)*den00(i-1)+rg*den00(i))
  enddo

! subtract net momentum
!     phip00=phip00-sum(phip00)/real(mpsi+1)

! d phi/dr, in equilibrium unit
  do i=0,mpsi
     rg=rg0+deltar*real(i)
     phip00(i)=-phip00(i)/rg
  enddo

! add FLR contribution using Pade approximation: b*<phi>=(1+b)*<n>
  phi00=den00*meshti/(qion*qion) !!rho0*rho0 ! potential in equilibrium unit
  do i=1,mpsi-1
     phip00(i)=(phip00(i)/aion+0.5_lk*(phi00(i+1)-phi00(i-1))/deltar)/meshni(i)
  enddo

! (0,0) mode potential store in phi00
  phi00=0.0_lk
  do i=1,mpsi
     phi00(i)=phi00(i-1)+0.5_lk*deltar*(phip00(i-1)+phip00(i))
  enddo
  if(izonal==0)phip00=0.0_lk

end subroutine zonal_old

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!XY new zonal solver

subroutine zonal
  use global_parameters
  use particle_array
  use field_array
  implicit none

  integer i
  real(lk) den00(0:mpsi),h1,h2,delpsi,temp,psidum(0:mpsi)

! zonal density
  phip00=qion*zonali*meshni+qelectron*fluidne00*real(magnetic)
  if(nhybrid>0)phip00=phip00+qelectron*zonale*meshne
  if(fload>0)phip00=phip00+qfast*zonalf*meshnf
  if(feload>0)phip00=phip00+qfaste*zonalfe*meshnfe

! (-0.0625 0.25 0.625 0.25 -0.0625) radial smoothing of (0,0) mode density phip00
  do i=1,2
     den00(0)=phip00(0)
     den00(mpsi)=phip00(mpsi)
     den00(1)=phip00(3)
     den00(mpsi-1)=phip00(mpsi-3)
     den00(2:mpsi-2)=phip00(0:mpsi-4)+phip00(4:mpsi)
     den00(1:mpsi-1)=0.625_lk*phip00(1:mpsi-1)+0.25_lk*(phip00(0:mpsi-2)+phip00(2:mpsi))-&
          0.0625_lk*den00(1:mpsi-1)
     phip00=den00
  enddo
  if(iload/=9)then
!   den00=<e dn>, phip00= d_<e dn>/d_psip, nonuniform psip grids
    den00=phip00
! initial point
    h1=psimesh(1)-psimesh(0)
    h2=psimesh(2)-psimesh(1)
    phip00(0)=(den00(1)-den00(0))*(h1+h2)/(h1*h2)+(den00(0)-den00(2))*h1/(h2*(h1+h2))
    do i=1,mpsi-1
       h1=psimesh(i)-psimesh(i-1)
       h2=psimesh(i+1)-psimesh(i)
       phip00(i)=(den00(i+1)-den00(i))*h1/(h2*(h1+h2))+(den00(i)-den00(i-1))*h2/(h1*(h1+h2))
    enddo
! final point
    h1=psimesh(mpsi-1)-psimesh(mpsi-2)
    h2=psimesh(mpsi)-psimesh(mpsi-1)
    phip00(mpsi)=(den00(mpsi)-den00(mpsi-1))*(h1+h2)/(h1*h2)+(den00(mpsi-2)-den00(mpsi))*h1/(h2*(h1+h2))

! for integration
    phi00=den00*b2m00/(aion*meshni)+gpsi200*dtndpsi*phip00/(qion*qion)
    phi00=phi00*pmarki

! first term in d_<phi>/d_psip
    den00=phip00*meshti/(qion*qion*meshni)

! phi00=r*E_r, E_r(a0)=0. Trapezoid rule
    phip00=0.0_lk
    do i=1,mpsi
       phip00(i)=phip00(i-1)+(phi00(i)+phi00(i-1))*0.5_lk
    enddo

! d phi/dr, in equilibrium unit
    do i=0,mpsi
       delpsi=0.5_lk*(psimesh(min(mpsi,i+1))-psimesh(max(0,i-1)))
       phip00(i)=den00(i)-phip00(i)*delpsi/(gpsi200(i)*pmarki(i))
    enddo

  else ! Fully kinetic ion

    den00=phip00*jacobianpsi
    phip00=0.0_lk
    do i=0,mpsi
! normalized by (lambdaD_i)^2/(R_0)^2, lambdaDi, is the Debye-length, and R_0, is the major radius
       temp=(etemp0)/(eden0)*7.43e2*7.43e2/(r0*r0)
       temp=1.0_lk
       den00(i)=-den00(i)/temp

       if(i>0)then
         phip00(i)=phip00(i-1)+(den00(i)+den00(i-1))*(psimesh(i)-psimesh(i-1))*0.5_lk
       endif
       psidum(i)=0.5_lk*(psimesh(min(mpsi,i+1))-psimesh(max(0,i-1)))
    enddo
    do i=0,mpsi
       phip00(i)=phip00(i)/(gpsi200(i)*jacobianpsi(i))
    enddo

! zonal potential are zero at both boundaries
! should be checked (animesh?)
    phip00=phip00-sum(phip00*psidum)/sum(psidum/(gpsi200*jacobianpsi))/(gpsi200*jacobianpsi)
  endif

! (0,0) mode potential store in phi00
  phi00=0.0_lk
  do i=1,mpsi
     phi00(i)=phi00(i-1)+(psimesh(i)-psimesh(i-1))*(phip00(i)+phip00(i-1))*0.5_lk
  enddo
  if(iload==9)then
    do i=0,mpsi
       phip00(i)=phip00(i)*rho0*rho0
       phi00(i)=phi00(i)*rho0*rho0
    enddo
  endif

end subroutine zonal


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gradient(scalar,vector)
  use global_parameters
  use field_array
  use particle_array,only:iload
  implicit none

  integer i,ii,ij,j,k,icount,idest,isource,isendtag,irecvtag,ierror,jt,&
      istatus(MPI_STATUS_SIZE)
  real(lk) difft(0:mpsi),diffz,pleft(mthetamax),pright(mthetamax),&
      sendl(mgrid),recvr(mgrid),sendr(mgrid),recvl(mgrid),&
      q,delq
  real(lk),dimension(mgrid),intent(in) :: scalar
  real(lk),dimension(3,0:1,mgrid),intent(out) :: vector
  real(lk) h0,h1,wp0,wp1,wt10,wt00,wt11,wt01,sdum0,sdum1,pdum0,pdum1
  integer jm,jp,jt11,jt21,jt31,jt01,jt10,jt20,jt30,jt00,iout,iin,ibin,ibout

! finite difference for e-field in equilibrium unit
  difft=0.5_lk/deltat
  diffz=0.5_lk/deltaz
!$omp parallel do private(i,j,k)
  do i=1,mgrid
     do j=0,1
        do k=1,3
          vector(k,j,i)=0.0_lk
        enddo
     enddo
  enddo

!! d_scalar/d_psito higher order
       if(iload==9)then
         ibin=0
         ibout=mpsi
       else
         ibin=1
         ibout=mpsi-1
       endif
!$omp parallel do private(i,j,ij,h0,h1,wp1,wp0,jm,jp,wt10,wt00,wt11,wt01,&
!$omp& jt11,jt21,jt31,jt01,jt10,jt20,jt30,jt00,sdum1,sdum0,pdum1,pdum0,iin,iout)
  do i=ibin,ibout

     iout=i+1
     iin=i-1

     if(i==0)iin=1
     if(i==mpsi)iout=mpsi-1

     h0=1.0_lk/(psimesh(i)-psimesh(iin))
     h1=1.0_lk/(psimesh(iout)-psimesh(i))

     if(i==0)h0=-h0
     if(i==mpsi)h1=-h1

     wp1=h1/(h0+h1)
     wp0=1.0_lk-wp1

     jm=igrid(iin)
     jp=igrid(iout)

     do j=1,mtheta(i)
        ij=igrid(i)+j
        wt10=wtp1(2,ij)   !upper poloidal grid on inner flux surface
        wt00=1.0_lk-wt10       !lower poloidal grid on inner flux surface
        wt11=wtp1(1,ij)   !upper poloidal grid on outer flux surface
        wt01=1.0_lk-wt11       !lower poloidal grid on outer flux surface

        jt11=jtp1(1,ij)   !!lower poloidal grid on outer flux surface
        jt21=jt11+1
        jt31=jt21+1-mtheta(iout)*((jt21-jp)/mtheta(iout))
        jt01=jp+mod(jt11-jp-1+mtheta(iout),mtheta(iout))

        jt10=jtp1(2,ij)   !!lower poloidal grid on the innoer flux surface
        jt20=jt10+1
        jt30=jt20+1-mtheta(iin)*((jt20-jm)/mtheta(iin))
        jt00=jm+mod(jt10-jm-1+mtheta(iin),mtheta(iin))

        sdum1=wt01*scalar(jt11)+wt11*scalar(jt21)
        sdum0=wt00*scalar(jt10)+wt10*scalar(jt20)
        pdum1=sdum1+0.5_lk*wt01*wt11*(sdum1-((1.0_lk+wt11)*scalar(jt31)+(1.0_lk+wt01)*scalar(jt01))/3.0_lk)
        pdum0=sdum0+0.5_lk*wt00*wt10*(sdum0-((1.0_lk+wt10)*scalar(jt30)+(1.0_lk+wt00)*scalar(jt00))/3.0_lk)
        if(i==0)then
          vector(1,1,ij)=h1*(pdum1-scalar(ij))
        elseif(i==mpsi)then
          vector(1,1,ij)=h0*(scalar(ij)-pdum0)
        else
        vector(1,1,ij)=wp1*h1*pdum1-wp0*h0*pdum0+scalar(ij)*(wp0*h0-wp1*h1)
        endif
     enddo
  enddo

! d_scalar/d_theta

!$omp parallel do private(i,j,ij,jt)
  do i=ibin,ibout
     do j=1,mtheta(i)
        ij=igrid(i)+j
        jt=j+1-mtheta(i)*(j/mtheta(i))
        vector(2,1,ij)=difft(i)*(scalar(igrid(i)+jt)-scalar(igrid(i)+j-1))
     enddo
  enddo

! send scalar to right and receive from left
  sendr=scalar(:)
  recvl=0.0_lk
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! send scalar to left and receive from right
  sendl=scalar(:)
  recvr=0.0_lk
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
      recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! unpack scalar_boundary and calculate E_zeta at boundaries toroidal
!$omp parallel do private(i,j,ii,jt,ij,pleft,pright)
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

! d_scalar/d_zeta
     do j=1,mtheta(i)
        ij=igrid(i)+j
        vector(3,1,ij)=(pright(j)-pleft(j))*diffz
     enddo
  enddo

! adjust the difference between safety factor q and qtinv for fieldline coordinate
!$omp parallel do private(i,j,q,delq,ij)
  do i=1,mpsi-1
     q=qmesh(i)
#ifdef _FRC
        delq=0.0_lk
#else
        delq=(1.0_lk/q-qtinv(i))
#endif
     
     do j=1,mtheta(i)
        ij=igrid(i)+j
        vector(3,:,ij)=vector(3,:,ij)+delq*vector(2,:,ij)
     enddo
  enddo

  call periodicity(vector(1,:,:))
  call periodicity(vector(2,:,:))
  call periodicity(vector(3,:,:))

end subroutine gradient
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grad2para(scalarin,scalarout)
  use global_parameters
  use field_array
  implicit none

  integer i,ii,ij,j,k,icount,idest,isource,isendtag,irecvtag,ierror,jt,&
      istatus(MPI_STATUS_SIZE)
  real(lk) difft(0:mpsi),diffz,pleft(mthetamax),pright(mthetamax),&
      sendl(mgrid),recvr(mgrid),sendr(mgrid),recvl(mgrid),q,delq,d2dt2(mgrid)
  real(lk),dimension(mgrid),intent(in) :: scalarin
  real(lk),dimension(0:1,mgrid),intent(out) :: scalarout

! finite difference for e-field in equilibrium unit
  difft=1.0_lk/deltat
  diffz=1.0_lk/deltaz
  d2dt2=0.0_lk
!$omp parallel do private(i,j,k)
  do i=1,mgrid
     do j=0,1
          scalarout(j,i)=0.0_lk
     enddo
  enddo

! d2_scalar/d_theta2
!$omp parallel do private(i,j,ij,jt)
  do i=1,mpsi-1
     do j=1,mtheta(i)
        ij=igrid(i)+j
        jt=j+1-mtheta(i)*(j/mtheta(i))
        d2dt2(ij)=difft(i)*difft(i)*(scalarin(igrid(i)+jt)+scalarin(igrid(i)+j-1)-2.0_lk*scalarin(ij))
     enddo
  enddo


! send scalar to right and receive from left
  sendr=scalarin(:)
  recvl=0.0_lk
  icount=mgrid
  idest=right_pe
  isource=left_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendr,icount,mpi_Rsize,idest,isendtag,&
      recvl,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! send scalar to left and receive from right
  sendl=scalarin(:)
  recvr=0.0_lk
  idest=left_pe
  isource=right_pe
  isendtag=myrank_toroidal
  irecvtag=isource
  call MPI_SENDRECV(sendl,icount,mpi_Rsize,idest,isendtag,&
      recvr,icount,mpi_Rsize,isource,irecvtag,toroidal_comm,istatus,ierror)

! unpack scalar_boundary and calculate E_zeta at boundaries toroidal
!$omp parallel do private(i,j,ii,jt,ij,pleft,pright,q,delq)
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

! d2_scalar/d_zeta2
     do j=1,mtheta(i)
        ij=igrid(i)+j
        scalarout(1,ij)=(pright(j)+pleft(j)-2.0_lk*scalarin(ij))*diffz*diffz
     enddo
  enddo

! adjust the difference between safety factor q and qtinv for fieldline coordinate
!$omp parallel do private(i,j,q,delq,ij)
  do i=1,mpsi-1
     q=qmesh(i)
     delq=(1.0_lk/(q*q)-qtinv(i)*qtinv(i))
     do j=1,mtheta(i)
        ij=igrid(i)+j
        scalarout(1,ij)=scalarout(1,ij)+delq*d2dt2(ij)
     enddo
  enddo
  call periodicity(scalarout(:,:))

end subroutine grad2para
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grad2psi(scalarin,scalarout)
  use global_parameters
  use field_array
  implicit none

  integer i,ij,j
  real(lk) h0,h1,wp0,wp1,wt10,wt00,wt11,wt01,sdum0,sdum1,pdum0,pdum1
  integer jm,jp,jt11,jt21,jt31,jt01,jt10,jt20,jt30,jt00
  real(lk),dimension(mgrid),intent(in) :: scalarin
  real(lk),dimension(0:1,mgrid),intent(out) :: scalarout
!! d2_scalar/d_psi2 higher order
!$omp parallel do private(i,j,ij,h0,h1,wp1,wp0,jm,jp,wt10,wt00,wt11,wt01,&
!$omp& jt11,jt21,jt31,jt01,jt10,jt20,jt30,jt00,sdum1,sdum0,pdum1,pdum0)
  do i=1,mpsi-1

     h0=1.0_lk/(psimesh(i)-psimesh(i-1))
     h1=1.0_lk/(psimesh(i+1)-psimesh(i))
     wp1=h1/(h0+h1)
     wp0=1.0_lk-wp1

     jm=igrid(i-1)
     jp=igrid(i+1)

     do j=1,mtheta(i)
        ij=igrid(i)+j
        wt10=wtp1(2,ij)   !upper poloidal grid on inner flux surface
        wt00=1.0_lk-wt10       !lower poloidal grid on inner flux surface
        wt11=wtp1(1,ij)   !upper poloidal grid on outer flux surface
        wt01=1.0_lk-wt11       !lower poloidal grid on outer flux surface

        jt11=jtp1(1,ij)   !!lower poloidal grid on outer flux surface
        jt21=jt11+1
        jt31=jt21+1-mtheta(i+1)*((jt21-jp)/mtheta(i+1))
        jt01=jp+mod(jt11-jp-1+mtheta(i+1),mtheta(i+1))

        jt10=jtp1(2,ij)   !!lower poloidal grid on the innoer flux surface
        jt20=jt10+1
        jt30=jt20+1-mtheta(i-1)*((jt20-jm)/mtheta(i-1))
        jt00=jm+mod(jt10-jm-1+mtheta(i-1),mtheta(i-1))

        sdum1=wt01*scalarin(jt11)+wt11*scalarin(jt21)
        sdum0=wt00*scalarin(jt10)+wt10*scalarin(jt20)
        pdum1=sdum1+0.5_lk*wt01*wt11*(sdum1-((1.0_lk+wt11)*scalarin(jt31)+(1.0_lk+wt01)*scalarin(jt01))/3.0_lk)
        pdum0=sdum0+0.5_lk*wt00*wt10*(sdum0-((1.0_lk+wt10)*scalarin(jt30)+(1.0_lk+wt00)*scalarin(jt00))/3.0_lk)
        scalarout(1,ij)=2.0_lk*h0*h1*(wp1*pdum1+wp0*pdum0-scalarin(ij))
     enddo
  enddo
  call periodicity(scalarout(:,:))
end subroutine grad2psi
!!!!!!!!!

subroutine laplacian_petsc(scalar,scalar_out)
  use precision
  use global_parameters,only: mgrid,mpsi,rho0
  use particle_array,only: aion,qion,meshti,meshni
  use field_array,only: igrid,mtheta,bmesh
  use petsc_array,only: userb,userx
  implicit none

  integer :: i,ij,j
  !real(lk),dimension(mpsi+1) :: ddum
  real(lk),dimension(mgrid),intent(in) :: scalar
  real(lk),dimension(mgrid),intent(out) :: scalar_out

  scalar_out=0.0_lk

!$omp parallel do private(i)
  do i=1,mgrid
     userb(i-1)=scalar(i)
  enddo

#ifdef _PETSc
  call laplacian_core
#endif

!$omp parallel do private(i,j,ij)
  do i=1,mpsi-1
     do j=0,mtheta(i)
        ij=igrid(i)+j
        scalar_out(ij)=-userx(ij-1)*bmesh(ij)*bmesh(ij)/(aion*rho0*rho0*meshni(i))
     enddo
  enddo

!*Boundary conditions at the inner wall and at the outer wall
!$omp parallel do private(i,j,ij)
  do i=0,mpsi,mpsi
    do j=0,mtheta(i)
      ij=igrid(i)+j
      scalar_out(ij)=0.0_lk
    enddo
  enddo

end subroutine laplacian_petsc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine laplacian(scalar,scalar_out,none0bound)
  use precision
  use global_parameters,only: mgrid,mpsilow,mpsihigh,mgridlow,mgridhigh,mpsi,mype,istep,&
   pi2,bcond,rg0!,rho0
 !! use particle_array,only: aion,qion,meshti
  use field_array,only: igrid,mtheta,lapmat,indexplap,nindexlap,&
   zeta1,qtinv,deltat,psimesh,deltar
  use interfaces,only: axisExtrapolate
  implicit none

  integer, optional, intent (in) :: none0bound


  integer :: i,ij,j,ii,jj,inone0bound
  real(lk),dimension(mgrid),intent(in) :: scalar
  real(lk),dimension(mgrid),intent(out) :: scalar_out
  real(lk) tdum,rdum,ave,sinave,cosave,ave2,sinave2,cosave2


  if(PRESENT(none0bound))then
      inone0bound=none0bound
  else
      inone0bound=0
  endif

!$omp parallel do private(i,j)
     do i=mgridlow,mgridhigh
        scalar_out(i)=0.0_lk
        do j=1,nindexlap(i)
           scalar_out(i)=scalar_out(i)+lapmat(j,i)*scalar(indexplap(j,i))
        enddo
     enddo

!For linear boundary conditions, Fourier decomposes mpsilow+1 flux surface and linear exterpolates
!solution to the magnetic axis.
  if(bcond==1)then
    call axisExtrapolate(scalar_out(:))
!Zero inner radial boundary condition for bcond=0
  else
    do j=0,mtheta(0)
      ij= igrid(0)+j
      scalar_out(ij)=0.0_lk
    enddo
  endif

!Outer radial boundary Condition
   do j=0,mtheta(mpsi)
      ij=igrid(mpsi)+j
      scalar_out(ij)=0.0_lk
   enddo


!correction for the mpihigh flux surface for nonzero outer boundary
   if(inone0bound > 0)then
      do j=0, mtheta(mpsihigh)
         ii = igrid(mpsihigh) + j
         do jj=nindexlap(ii)+1,11
            scalar_out(ii)=scalar_out(ii)+lapmat(jj,ii)*scalar(indexplap(jj,ii))
         enddo
      enddo
   endif

end subroutine laplacian
