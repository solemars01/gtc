subroutine snapshot
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use magnetic_island
  use sdp, only: sdp_perturbation_output
  use spline_function, only: spx, spz
  implicit none
  
  integer,parameter :: nfield=3,nvgrid=65,iosnap=222,iophi=223
  integer mtgrid,nf,ns,m,ip,jenergy,kpitch,nsnap,i,ierror,j,ij,jt,icount,isp,jst
  real(lk) upara,delf,rmu,fullf,ai_inv,vi_inv,ae_inv,ve_inv,afe_inv,vfe_inv,emax_inv,delr,rg,b,&
       energy,pitch,wt,tdum,pdum,dpx,dp2,dtx,dt2,psitmp,thetatmp,&
       profile(0:mpsi,6,nspecies),pdf(nvgrid,4,nspecies),&
       dfield(mgrid),proftmp(0:mpsi,6,nspecies),pdftmp(nvgrid,4,nspecies),niflux(0:mpsi),eqmeshni(0:mpsi)
  real(lk),dimension(:),allocatable :: eachflux
  real(lk),dimension(:,:),allocatable :: allflux
  real(lk),dimension(:,:,:),allocatable :: poloidata,fluxdata
  character(len=64) cdum0,cdum1,cdum2

! particle species: 1=ion, 2=electron, 3=EP
! radial profiles: density, flow, energy
  profile=0.0_lk

! distribution function: energy, pitch angle
  pdf=0.0_lk

! species ns=1: thermal ion
  ns=1
  ai_inv=1.0_lk/aion
  vi_inv=sqrt(aion)/rho0 !parallel velocity normalized by c_s
  delr=1.0_lk/deltar
  emax_inv=0.2_lk


!!$omp parallel do private(i,pdum,isp,dpx,dp2)
  do i=0,mpsi
     pdum=psimesh(i)
     isp=max(1,min(lsp-1,ceiling(pdum*spdpsi_inv)))
     dpx=pdum-spdpsi*real(isp-1)
     dp2=dpx*dpx

     eqmeshni(i)=nipp(1,isp)+nipp(2,isp)*dpx+nipp(3,isp)*dp2
 enddo
  if(iload/=9)then  
  do m=1,mi
     psitmp=zion(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zion(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
    

     fullf=zion(7,m)
     delf=fullf*zion(5,m)
     rmu=zion(6,m)
     upara=zion(4,m)*b*qion*ai_inv*vi_inv 
     energy=0.5_lk*upara*upara+rmu*rmu*b*ai_inv*vi_inv*vi_inv ! energy is normalized by etemp0
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tipp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5_lk*(pitch+1.0_lk))))! pitch bin

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo
  else
  do m=1,mi
     psitmp=zion(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zion(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2


     fullf=zion(7,m)
     delf=fullf*zion(5,m)
     rmu=zion(6,m)
     upara=zion(4,m)*b*qion*ai_inv*vi_inv
     energy=0.5_lk*upara*upara+rmu*rmu*b*ai_inv*vi_inv*vi_inv ! energy is normalized by etemp0
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tipp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5*(pitch+1.0))))! pitch bin

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo

  endif
! species ns=2: electron
  if(nhybrid>0)then
  ns=2
  ae_inv=1.0/aelectron
  ve_inv=sqrt(aelectron)/rho0 ! electron velocity normalized by v_the
  do m=1,me

     psitmp=zelectron(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zelectron(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

     fullf=zelectron(7,m)
     delf=fullf*zelectron(5,m)
     rmu=zelectron(6,m)
     upara=zelectron(4,m)*b*qelectron*ae_inv*ve_inv 
     energy=0.5_lk*upara*upara+rmu*rmu*b*ae_inv*ve_inv*ve_inv ! energy is normalized by etemp0
     pitch=upara/sqrt(2.0_lk*energy)

     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv))) !energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5_lk*(pitch+1.0_lk))))  !pitch bin

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo
  endif

! species ns=3: fast ion
  if(fload>0)then
  ns=ns+1
  ai_inv=1.0/afast
  vi_inv=sqrt(afast)/rho0
  
  do m=1,mf
     psitmp=zfast(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zfast(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

     fullf=zfast(7,m)
     delf=fullf*zfast(5,m)
     rmu=zfast(6,m)
     upara=zfast(4,m)*b*qfast*ai_inv*vi_inv 
     energy=0.5_lk*upara*upara+rmu*rmu*b*ai_inv*vi_inv*vi_inv
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tfpp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5_lk*(pitch+1.0))))! pitch bin

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0_lk-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0_lk-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0_lk-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0_lk-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0_lk-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0_lk-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo
  endif
  
! species ns=4: fast electron
  if(feload>0)then
  ns=ns+1
  afe_inv=1.0_lk/afaste
  vfe_inv=sqrt(afaste)/rho0

  do m=1,mfe
     psitmp=zfaste(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of rg & b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline is regular
     thetatmp=zfaste(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

     rg= rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     b=    bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
          (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
          (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

     fullf=zfaste(7,m)
     delf=fullf*zfaste(5,m)
     rmu=zfaste(6,m)
     upara=zfaste(4,m)*b*qfaste*afe_inv*vfe_inv
     energy=0.5_lk*upara*upara+rmu*rmu*b*afe_inv*vfe_inv*vfe_inv
     pitch=upara/sqrt(2.0_lk*energy)

! particles sorted into bins in radius, energy, and pitch
     ip=max(1,min(mpsi,ceiling((rg-rg0)*delr)))! radial bin
     jenergy=max(1,min(nvgrid,ceiling(real(nvgrid)*energy*emax_inv*rho0*rho0/tfepp(1,1))))! energy bin
     kpitch=max(1,min(nvgrid,ceiling(real(nvgrid)*0.5*(pitch+1.0))))! pitch bin

     pdum=(psimesh(ip)-psitmp)/deltap(ip)
     profile(ip-1,1,ns)=profile(ip-1,1,ns)+pdum*fullf
     profile(ip,1,ns)=profile(ip,1,ns)+(1.0-pdum)*fullf
     profile(ip-1,2,ns)=profile(ip-1,2,ns)+pdum*delf
     profile(ip,2,ns)=profile(ip,2,ns)+(1.0-pdum)*delf
     profile(ip-1,3,ns)=profile(ip-1,3,ns)+pdum*fullf*upara
     profile(ip,3,ns)=profile(ip,3,ns)+(1.0-pdum)*fullf*upara
     profile(ip-1,4,ns)=profile(ip-1,4,ns)+pdum*delf*upara
     profile(ip,4,ns)=profile(ip,4,ns)+(1.0-pdum)*delf*upara
     profile(ip-1,5,ns)=profile(ip-1,5,ns)+pdum*fullf*energy
     profile(ip,5,ns)=profile(ip,5,ns)+(1.0-pdum)*fullf*energy
     profile(ip-1,6,ns)=profile(ip-1,6,ns)+pdum*delf*energy
     profile(ip,6,ns)=profile(ip,6,ns)+(1.0-pdum)*delf*energy

     pdf(jenergy,1,ns)=pdf(jenergy,1,ns)+fullf
     pdf(jenergy,2,ns)=pdf(jenergy,2,ns)+delf
     pdf(kpitch,3,ns)=pdf(kpitch,3,ns)+fullf
     pdf(kpitch,4,ns)=pdf(kpitch,4,ns)+delf
  enddo
  endif

! sum over MPI tasks
  icount=nspecies*6*(mpsi+1)
  call MPI_REDUCE(profile,proftmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  profile=proftmp

  icount=nspecies*4*(nvgrid)
  call MPI_REDUCE(pdf,pdftmp,icount,mpi_Rsize,MPI_SUM,0,MPI_COMM_WORLD,ierror)
  pdf=pdftmp
  
! normalization by marker #
  do ns=1,nspecies
     do ip=0,mpsi
        profile(ip,2:6,ns)=profile(ip,2:6,ns)/profile(ip,1,ns)
     enddo
     if(ns==1)then
        profile(:,1,ns)=profile(:,1,ns)/pmarki

     elseif(ns==2)then
        if(nhybrid>0)then
           profile(:,1,ns)=profile(:,1,ns)/pmarke
        else
           if(fload>0)then
              profile(:,1,ns)=profile(:,1,ns)/pmarkf
           else
               if(feload>0)then
                 profile(:,1,ns)=profile(:,1,ns)/pmarkfe
              endif
           endif 
        endif
     elseif(ns==3)then
          if(nhybrid>0)then
           if(fload>0)then
              profile(:,1,ns)=profile(:,1,ns)/pmarkf
           else
              if(feload>0)then
                 profile(:,1,ns)=profile(:,1,ns)/pmarkfe
              endif
           endif
        else
           profile(:,1,ns)=profile(:,1,ns)/pmarkfe
        endif
     elseif(nf==4)then
        profile(:,1,ns)=profile(:,1,ns)/pmarkfe
     endif

     do j=1,nvgrid
        pdf(j,2,ns)=pdf(j,2,ns)/max(1.0,pdf(j,1,ns))
        pdf(j,4,ns)=pdf(j,4,ns)/max(1.0,pdf(j,3,ns))
     enddo
  enddo

! poloidal resolution=poloidal grid on iflux surface
  mtgrid=mtheta(iflux)

  allocate(poloidata(0:mtgrid,0:mpsi,nfield+2),fluxdata(0:mtgrid,mtoroidal,nfield),&
       eachflux(mtgrid),allflux(mtgrid,mtoroidal))
  poloidata=0.0_lk
  fluxdata=0.0_lk

! field quantities: phi, a_para, fluidne. Last two coloumn of poloidal for coordinates
  do nf=1,nfield
     if(nf==1)then
        dfield=phi(0,:)
     endif
     if(nf==2)then
#ifdef _FRC
             dfield=densityi(1,:)
#else
             dfield=sapara(0,:)
#endif
     endif
     if(nf==3)then
#ifdef _FRC
        if(feload>0)then
            dfield=densityfe(1,:)
        elseif(nhybrid>0)then
            dfield=densitye(1,:)
        endif
#else       
            dfield=sfluidne(0,:)
#endif
    endif 

! gather potential on flux surface
     allflux=0.0
     do j=1,mtgrid
        eachflux(j)=dfield(igrid(iflux)+j)
     enddo

     icount=mtgrid
     call MPI_GATHER(eachflux,icount,mpi_Rsize,allflux,icount,mpi_Rsize,0,toroidal_comm,ierror) 
     fluxdata(1:mtgrid,:,nf)=allflux

! poloidal BC
     fluxdata(0,:,nf)=fluxdata(mtgrid,:,nf)

! potential data on poloidal plain uses polar coordinates
     do j=0,mtgrid
        tdum=pi2*real(j)/real(mtgrid)
        do i=0,mpsi
           if(fielddir==1 .or. fielddir==3)then
             jt=max(0,min(mtheta(i),1+int(modulo(tdum+pi2*qtinv(i),pi2)/deltat(i))))
             wt=modulo(tdum+pi2*qtinv(i),pi2)/deltat(i)-real(jt-1)
           else
             jt=max(0,min(mtheta(i),1+int(tdum/deltat(i))))
             wt=tdum/deltat(i)-real(jt-1)
           endif
           poloidata(j,i,nf)=(wt*dfield(igrid(i)+jt)+(1.0-wt)*dfield(igrid(i)+jt-1)) !/rho0**2
           if(nf==1)then
              poloidata(j,i,nf)=poloidata(j,i,nf)/(rho0*rho0)
           elseif(nf==2)then
! apara is renormalized in such a way that it has the same
! amplitude as phi in ideal shear Alfven waves
#ifndef _FRC
              poloidata(j,i,nf)=poloidata(j,i,nf)/(rho0*sqrt(betae*aion))
#endif
! no need to renormalize fluidne
           endif
        enddo
     enddo
  enddo

! poloidal grid position in polar coordinates
  do j=0,mtgrid
     tdum=2.0_lk*pi*real(j)/real(mtgrid)
     do i=0,mpsi
        pdum=psimesh(i)
        poloidata(j,i,nfield+1)=spx(pdum,tdum)
        poloidata(j,i,nfield+2)=spz(pdum,tdum)
     enddo
  enddo


  niflux=0.0_lk
  do i=0,mpsi
     do j=0,mtheta(i)
        ij=igrid(i)+j
        niflux(i)=niflux(i)+densityi(0,ij)
     enddo
     niflux(i)=niflux(i)/(mtheta(i)+1)
  enddo

! open snapshot output file
  if(mype==0)then
     nsnap=mstepall+istep
     write(cdum0,'(i7.7,".out")')nsnap
     cdum1='snap'//trim(cdum0)
     open(iosnap,file=cdum1,status='replace')

! parameters: # of species, fields, and grids in velocity, radius, poloidal, toroidal; T_up
     write(iosnap,101)nspecies,nfield,nvgrid,mpsi+1,mtgrid+1,mtoroidal
     write(iosnap,102)1.0/emax_inv

! write out particle radial profile and pdf, and 2D field
     write(iosnap,102)profile,pdf,poloidata,fluxdata

! close snapshot file
     close(iosnap)
     
     if(island==1)then
     cdum2='density'//trim(cdum0)
     open(6570,file=cdum2,status='replace') 
     do i=0, mpsi
        ij=igrid(i)
        !write(6570,*) i,eqmeshni(i),densityi(0,ij+mtheta(i)/2-3),densityi(0,ij+mtheta(i)/2-2),densityi(0,ij+mtheta(i)/2-1),densityi(0,ij+mtheta      (i)/2),densityi(0,ij+mtheta(i)/2+1),densityi(0,ij+mtheta(i)/2+2),densityi(0,ij+mtheta(i)/2+3)
!        write(6570,*) i,eqmeshni(i),eqmeshni(i)+densityi(0,ij+mtheta(i)/2),eqmeshni(i)+densityi(0,ij),niflux(i)
        write(6570,*) i,(densityi(1,ij+mtheta(i)/2)+1.0)*eqmeshni(i),eqmeshni(i)*(densityi(1,ij)+1.0),(densitye(1,ij)+1)*eqmeshni(i),(densitye(1,ij+mtheta(i)/2)+1)*eqmeshni(i),eqmeshni(i)*(1.0+pressureepara(1,ij))
        enddo
     close(6570)
     endif
  endif

! Write out files for Synthetic Diagnostics Platform
  if(mype==0) then
     if(sdp_output==1) then
        call sdp_perturbation_output
     endif
  endif

101 format(i6)
102 format(e13.6)
  
end subroutine snapshot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine phishot

use global_parameters
use field_array
use particle_array, only:zonali

  integer,parameter :: iophi=223
  integer i,nsnap
  character(len=64) cdum

  if (mype==0) then
      nsnap=mstepall+istep
      write(cdum,'(i7.7,".out")')nsnap
      cdum='phi'//trim(cdum)
      open(iophi,file=cdum,status='replace') 
      do i=0,mpsi
        write(iophi,*)i,phi(1,igrid(i)),phi00(i),zonali(i),phi(1,igrid(i))+phi00(i)
      enddo
      close(iophi)
   endif

end subroutine phishot
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine frc_snapshot
  use precision
  use global_parameters
  use field_array
  use particle_array!,only:densityi,densityfe
  use LAMYRIDGEeq
  use cylindricalRZ
  implicit none

  integer :: i,j,ij,frcsnap,icount,ierror,nsnap
  character(len=99) :: fname,format_string,snapshot_fname
  101 format(' ',4g16.4)
  102 format(' ',9g16.8)
  103 format(' ',8g16.4)
if(mype==0)then
    frcsnap=321 !open frc equilibrium plot output file

    format_string = "(A16,I0.5)" ! 7 characters for "restart"
                                 ! 5 characters for mpi id with padded zeros
    nsnap=mstepall+istep
    write(fname,'(i7.7,".out")')nsnap
    snapshot_fname='frc_CYL_snapshot'//trim(fname)
    open(frcsnap,file=snapshot_fname,status='replace')
 !    write(frcsnap,102)r_Eq,z_Eq,psirz_Eq(1,:,:),B0Mag_Eq(1,:,:),B0R_Eq(1,:,:),B0Z_Eq(1,:,:)
    write(frcsnap,102)RZ_reg(:,1),RZ_reg(:,2),psiEq0,psifrc_tmp,psi_frc,meshni_cy,meshne_cy,meshti_cy,meshte_cy,deltar_frc
    close(frcsnap)
endif
end subroutine frc_snapshot
