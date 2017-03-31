!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lorentz(species_name)
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use precision
  implicit none
  
  character(*),intent(in) :: species_name
  
  interface
    subroutine collision_pitch(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                               ncyclep)
      use precision
      !declaration of the dummy arguments
      integer mp,pcoll
      real(lk) qpart,apart,taupp
      real(lk),dimension(0:) :: mesht,meshn
      real(lk),dimension(:,:) :: tppp
      real(lk),dimension(:,:) :: zpart
      integer,optional :: ncyclep
    end subroutine collision_pitch
  end interface

  !call collision_pitch(zion,meshti,meshni,qion,aion,tipp,mi,icoll,tauii)
  if(species_name=="fast-electron")then
    call collision_pitch(zfaste,meshtfe,meshnfe,qfaste,afaste,tfepp,mfe,ecoll,tauei)!,ncyclefe)
  elseif(species_name=="thermal-electron")then  
    call collision_pitch(zelectron,meshte,meshne,qelectron,aelectron,tepp,me,ecoll,tauei)!,ncyclee)  
  endif

end subroutine lorentz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine collision_pitch(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                        ncyclep)
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  implicit none

  !declaration of the dummy arguments
  integer mp,pcoll
  real(lk) qpart,apart,taupp
  real(lk),dimension(0:) :: mesht,meshn
  real(lk),dimension(:,:) :: tppp
  real(lk),dimension(:,:) :: zpart
  integer,optional :: ncyclep
  
  !declaration of the local arguments
  integer m,isp,jst,ii
  real(lk) v,dele,zv,ap_inv,ve,nuei,freq,upara,vperp2,pitch,b,randome(mp),rg,dp1,&
       psitmp,dpx,dp2,thetatmp,dtx,dt2,zeff,eden,etemp,ve0,ptmp,tp_inv,delr,delp(mpsi)

  ap_inv=1.0_lk/apart
  tp_inv=1.0_lk/tppp(1,1)
  delr=1.0_lk/deltar
  delp=1.0_lk/deltap
  ve0=rho0*sqrt(2.0*ap_inv)
  
  if(present(ncyclep))then
    nuei=1.88_lk*real(pcoll)*tstep/(taupp*2.0_lk*real(ncyclep))
  else
    nuei=1.88_lk*real(pcoll)*tstep/(taupp*2.0_lk)
  endif
  
  call random_number(randome(1:mp))

!$omp parallel do private(m,psitmp,isp,dpx,dp2,thetatmp,jst,dtx,dt2,rg,ii,dp1,zeff,eden,etemp,b,&
!$omp& upara,vperp2,v,pitch,ve,zv,freq,ptmp)
  do m=1,mp

! radial spline grid
     psitmp=zpart(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
!     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline grid
     thetatmp=zpart(2,m)
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx

! 1D spline in radial
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     zeff=dp1*meshze(ii-1)+(1.0_lk-dp1)*meshze(ii)
     eden=dp1*meshn(ii-1)+(1.0_lk-dp1)*meshn(ii)
     etemp=(dp1*mesht(ii-1)+(1.0_lk-dp1)*mesht(ii))*tp_inv

! 2D spline in (psi, theta) for B-field
     b= bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
       (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
       (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2
     
     upara=zpart(4,m)*b*abs(qpart)*ap_inv
     vperp2=zpart(6,m)*zpart(6,m)*2.0_lk*b*ap_inv
     v=sqrt(upara*upara+vperp2)
     pitch=upara/v
     
     ve=ve0*sqrt(etemp)
     zv=max(0.1_lk,min(10.0_lk,v/ve))
     
! nui_ei for electron-ion collision
     freq=nuei*zeff*eden/(zv*zv*zv*etemp**1.5_lk)

! uniform square Montle-Carlo method
     pitch=pitch*(1.0_lk-freq)+(randome(m)-0.5_lk)*sqrt(12.0_lk*freq*max(1.0e-10,1.0_lk-pitch*pitch))
     ptmp=aint(pitch)
     if(abs(ptmp) > 0.5_lk)pitch=sign(1.0_lk,pitch)-(pitch-ptmp)

     upara=v*pitch
     zpart(4,m)=upara*apart/(abs(qpart)*b)
     vperp2=v*v-upara*upara
     zpart(6,m)=sqrt(0.5_lk*apart*(vperp2)/b)
  enddo
  
  return
end subroutine collision_pitch
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fokkerplanck(species_name)
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  use precision
  implicit none
  character(*),intent(in) :: species_name
  
  interface
    subroutine collision_fokkerplanck(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                                        ncyclep,z_eff)
        use precision
        !declaration of the dummy arguments
        integer mp,pcoll
        real(lk) qpart,apart,taupp
        real(lk),dimension(0:) :: mesht,meshn
        real(lk),dimension(:,:) :: tppp
        real(lk),dimension(:,:) :: zpart
        integer,optional :: ncyclep,z_eff
    end subroutine collision_fokkerplanck
  end interface
  
  if(species_name=='fast-electron')then
#ifdef _FRC
    call collision_fokkerplanck(zfaste,meshtfe,meshnfe,qfaste,afaste,tfepp,mfe,ecoll,tauee)
#else 
    call collision_fokkerplanck(zfaste,meshtfe,meshnfe,qfaste,afaste,tfepp,mfe,ecoll,tauee)!,ncyclep=ncyclefe)
#endif
  elseif(species_name=='thermal_electron')then
    call collision_fokkerplanck(zelectron,meshte,meshne,qelectron,aelectron,tepp,me,ecoll,tauee)!,ncyclep=ncyclee)  
  elseif(species_name=='thermal_ion')then
    call collision_fokkerplanck(zion,meshti,meshni,qion,aion,tipp,mi,icoll,tauii,z_eff=1)
  endif

end subroutine fokkerplanck
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine collision_fokkerplanck(zpart,mesht,meshn,qpart,apart,tppp,mp,pcoll,taupp,&
                                  ncyclep,z_eff)
  use precision
  use global_parameters
  use particle_array
  use field_array
  use equilibrium
  implicit none

  !declaration of the dummy arguments
  integer mp,pcoll
  real(lk) qpart,apart,taupp
  real(lk),dimension(0:) :: mesht,meshn
  real(lk),dimension(:,:) :: tppp
  real(lk),dimension(:,:) :: zpart
  integer,optional :: ncyclep,z_eff

  integer m,k,icount,ierror,ip,jt,kz,mcell(mp),isp,jst,ii
  real(lk) upara,vperp2,delu,delv2,dneop,dneot,dneoz,psitmp,thetatmp,zetatmp,dpx,dp2,dtx,dt2,b,&
       v,xv,zv,vmin,freq,phix,dphi,f,g,h,sp,sn,dp,dn,dpn,ap_inv,vp,vp0,nupp,pden,ptemp,zeff,&
       random1(mp),random2(mp),bgc(mp),delm(neop*neot*neoz),dele(neop*neot*neoz),meshz(0:mpsi),&
       marker(neop*neot*neoz),adum(neop*neot*neoz),uflow(neop*neot*neoz),tp_inv,rg,dp1,delr,delp(mpsi)

  ap_inv=1.0_lk/apart
  tp_inv=1.0_lk/tppp(1,1)
  delr=1.0_lk/deltar
  delp=1.0_lk/deltap
  vp0=sqrt(2.0_lk*tppp(1,1)*ap_inv)
  if(present(ncyclep))then
    nupp=1.88_lk*real(pcoll)*tstep/(taupp*2.0_lk*real(ncyclep))
  else
    nupp=1.88_lk*real(pcoll)*tstep/taupp
  endif
  
  meshz=1.0_lk
  if(present(z_eff))meshz=meshze
          
  vmin=1.0e-10*vp0*vp0
  dneop=real(neop)/(psi1-psi0)
  dneot=real(neot)/pi2
  dneoz=real(neoz)/pi2
 
! GC cell location, B-field, & velocity
!$omp parallel do private(m,psitmp,thetatmp,zetatmp,ip,jt,kz,isp,dpx,dp2,jst,dtx,dt2,&
!$omp& b,upara,vperp2)
  do m=1,mp
     psitmp=zpart(1,m)
     thetatmp=zpart(2,m)
     zetatmp=zpart(3,m)
     
! neoclassical cell in psi,theta, zeta
     ip=max(1,min(neop,ceiling((psitmp-psi0)*dneop)))
     jt=max(1,min(neot,ceiling(thetatmp*dneot)))           
     kz=max(1,min(neoz,ceiling(zetatmp*dneoz)))         
! GC cell
     mcell(m)=(kz-1)*neop*neot+(jt-1)*neop+ip

! radial spline grid
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
! radaial spline of b avoids sigularity near axis: y=y1+y2*sqrt(x)+y3*x for psi<spdpsi
!     if(isp==1)dpx=sqrt(dpx)
     dp2=dpx*dpx

! poloidal spline grid
     jst=max(1,min(lst-1,ceiling(thetatmp*spdtheta_inv)))
     dtx=thetatmp-spdtheta*real(jst-1)
     dt2=dtx*dtx
     
! B-field 2D spline in (psi, theta)
     b= bsp(1,isp,jst)+bsp(2,isp,jst)*dpx+bsp(3,isp,jst)*dp2+ &
       (bsp(4,isp,jst)+bsp(5,isp,jst)*dpx+bsp(6,isp,jst)*dp2)*dtx+ &
       (bsp(7,isp,jst)+bsp(8,isp,jst)*dpx+bsp(9,isp,jst)*dp2)*dt2

     upara=zpart(4,m)*b*qpart*ap_inv
     vperp2=zpart(6,m)*zpart(6,m)*b*2.0_lk*ap_inv
     
! temporal storage of GC B-field, parallel & perpendicular velocity
     bgc(m)=b
     zpart(4,m)=upara
     zpart(6,m)=vperp2
  enddo

! center-of-mass coordinats; required OpenMP vector parallelization
  uflow=0.0_lk
  marker=0.0_lk
  do m=1,mp
     upara=zpart(4,m)
     ip=mcell(m)
     uflow(ip)=uflow(ip)+upara
     marker(ip)=marker(ip)+1.0_lk
  enddo

! global sum
  icount=neop*neot*neoz
  call MPI_ALLREDUCE(uflow,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  uflow=adum
  call MPI_ALLREDUCE(marker,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  marker=max(1.0_lk,adum)
  uflow=uflow/marker

! random number for ion-ion collision
  call random_number(random1(1:mp))
  call random_number(random2(1:mp))

! velocity changes due to ion-ion collision
!$omp parallel do private(m,psitmp,isp,dpx,dp2,rg,ii,dp1,pden,ptemp,ip,upara,vperp2,v,vp,zv,xv,&
!$omp& freq,k,phix,dphi,f,g,h,sp,sn,dp,dn,dpn,delu,delv2,zeff)
  do m=1,mp

! local ion density & temporature
     psitmp=zpart(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
     dp2=dpx*dpx  
     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2

     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     pden=dp1*meshn(ii-1)+(1.0_lk-dp1)*meshn(ii)
     ptemp=(dp1*mesht(ii-1)+(1.0_lk-dp1)*mesht(ii))*tp_inv
     zeff=dp1*meshz(ii-1)+(1.0_lk-dp1)*meshz(ii)

! transform to center-of-mass for like-species collision
     ip=mcell(m)
     zpart(4,m)=zpart(4,m)-uflow(ip)
     upara=zpart(4,m)
     vperp2=zpart(6,m)
     v=sqrt(upara*upara+vperp2)

! velocity normalized by local temporature
     vp=vp0*sqrt(ptemp)
     zv=max(0.1_lk,min(10.0_lk,v/vp))
     xv=zv*zv
     
! ion-ion collision frequency
     freq=zeff*nupp*pden/(zv*xv*ptemp**1.5)

! Maxwellian integral by table look-up for range of (0.0001,10.0)
     k = min(100000,int(xv*10000.0_lk) + 1)
     phix = (real(k)-10000.0_lk*xv)*maxwell(k) + (10000.0_lk*xv-real(k-1))*maxwell(k+1)
     if(xv < 0.025_lk)phix=4.0_lk/3.0_lk*sqrt(xv/pi)*xv*(1.0_lk-0.6_lk*xv+3.0_lk/14.0_lk*xv*xv)
     if(xv > 10.0_lk)phix=1.0_lk-2.0_lk/exp(xv)*sqrt(xv/pi)*(1.0_lk+1.0_lk/(2.0_lk*xv)-1.0_lk/(4.0_lk*xv*xv))
     dphi = 2.0_lk*sqrt(xv/pi)/exp(xv)

! coefficients for like-species collisions
     f=freq*2.0_lk*phix
     g=freq*(phix-0.5_lk*phix/xv+dphi)
     h=freq*phix/xv
     
     sp=upara*f
     sn=vperp2*(2.0_lk*f-h-g)-2.0_lk*upara*upara*g
     dp=upara*upara*h+vperp2*g
     dn=4.0_lk*vperp2*v*v*v*v*g*h/dp
     dpn=2.0_lk*vperp2*upara*(h-g)
     
! parallel and perpendicular drag and diffusion
     delu= (random1(m)-0.5_lk)*sqrt(12.0_lk*dp)-sp
     delv2=(random1(m)-0.5_lk)*dpn*sqrt(12.0_lk/dp)+(random2(m)-0.5_lk)*sqrt(12.0_lk*dn)-sn

! temporal storage of velocity changes
     random1(m)=delu
     random2(m)=delv2
  enddo

! momentum and energy changes due to collisions, require OpenMP vector parallelization
  dele=0.0_lk
  delm=0.0_lk
  do m=1,mp
     ip=mcell(m)
     delm(ip)=delm(ip)+zpart(5,m)*random1(m)
     dele(ip)=dele(ip)+zpart(5,m)*(random2(m)+random1(m)*(2.0_lk*zpart(4,m)+random1(m)))
  enddo

! global sum
  call MPI_ALLREDUCE(delm,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  delm=adum
  call MPI_ALLREDUCE(dele,adum,icount,MPI_Rsize,MPI_SUM,MPI_COMM_WORLD,ierror)
  dele=adum
  delm=sqrt(4.5_lk*pi)*2.0_lk*delm/(vp0*vp0*marker)
  dele=sqrt(4.5_lk*pi)*dele/(1.5_lk*vp0*vp0*marker)

! local conservation of momentum and energy
!$omp parallel do private(m,upara,vperp2,psitmp,isp,dpx,dp2,rg,ii,dp1,ptemp,vp,v,zv,xv,k,phix,dphi,ip)
  do m=1,mp
! new parallel and perpendicular velocity
     upara=zpart(4,m)+random1(m)
     vperp2=max(vmin,zpart(6,m)+random2(m))

! local ion temporature
     psitmp=zpart(1,m)
     isp=max(1,min(lsp-1,ceiling(psitmp*spdpsi_inv)))
     dpx=psitmp-spdpsi*real(isp-1)
     dp2=dpx*dpx     

     rg=  rgpsi(1,isp)  +rgpsi(2,isp)*dpx  +rgpsi(3,isp)*dp2
     ii=max(1,min(mpsi,ceiling((rg-rg0)*delr)))    !radial grid on inner flux surface
     dp1=(psimesh(ii)-psitmp)*delp(ii) !weight for inner flux surface
     ptemp=(dp1*mesht(ii-1)+(1.0-dp1)*mesht(ii))*tp_inv

     vp=vp0*sqrt(ptemp)
     v=sqrt(upara*upara+vperp2)
     zv=max(0.1_lk,min(10.0_lk,v/vp))
     xv=zv*zv
    
     k = min(100000,int(xv*10000.0_lk) + 1)
     phix = (k-10000.0_lk*xv)*maxwell(k) + (10000.0_lk*xv-k+1)*maxwell(k+1)
     if(xv .lt. 0.025_lk)phix=4.0_lk/3.0_lk*sqrt(xv/pi)*xv*(1.0_lk-0.6_lk*xv+3.0_lk/14.0_lk*xv*xv)
     if(xv .gt. 10.0_lk)phix=1.0_lk-2.0_lk/exp(xv)*sqrt(xv/pi)*(1.0_lk+1.0_lk/(2.0_lk*xv)-1.0_lk/(4.0_lk*xv*xv))
     dphi = 2.0_lk*sqrt(xv/pi)/exp(xv)
     
     ip=mcell(m)     
     zpart(5,m)=zpart(5,m)-phix/(xv*zv)*upara*delm(ip)-(phix-dphi)/zv*dele(ip)
     zpart(4,m)=(upara+uflow(ip))*apart/(qpart*bgc(m))
     zpart(6,m)=sqrt(0.5_lk*vperp2*apart/bgc(m))
  enddo

end subroutine collision_fokkerplanck
