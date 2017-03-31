subroutine charge(species_name)
  use global_parameters,only: ihybrid,nhybrid
  use particle_array
  use field_array,only: sfluidne
  implicit none

  interface
    subroutine gkChargeParticle(zpart,wppart,wtpart0,wtpart1,density,&
        flow,jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,&
        qpart,apart,pload,ngyro,mp,ihybrid,nhybrid,&
        pressurepara,pressureperp,sfluidn,dnsave,&
        switch)
      use precision
      implicit none

      integer pload,ngyro,mp
      integer,optional :: ihybrid,nhybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(lk) qpart,apart
      real(lk),dimension(:) :: wzpart,marker,markert
      real(lk),dimension(0:) :: zonal,zonalc,pmark
      real(lk),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(lk),dimension(0:,:) :: density,flow
      real(lk),dimension(:,:) :: zpart
      real(lk),dimension(:,:),optional :: pressurepara,pressureperp,sfluidn
      real(lk),dimension(:,:,:),optional :: dnsave
      character(*),intent(in),optional :: switch
    end subroutine gkChargeParticle

    subroutine hybridChargeParticle(zpart,wppart,wtpart0,wtpart1,density,&
        flow,jtpart0,jtpart1,wzpart,zonal,zonalc,marker,markert,pmark,&
        qpart,apart,pload,ngyro,mp,ihybrid,nhybrid,&
        pressurepara,pressureperp,sfluidn,dnsave,&
        switch)
      use precision
      implicit none

      integer pload,ngyro,mp
      integer,optional :: ihybrid,nhybrid
      integer,dimension(:,:) :: jtpart0,jtpart1
      real(lk) qpart,apart
      real(lk),dimension(:) :: wzpart,marker,markert
      real(lk),dimension(0:) :: zonal,zonalc,pmark
      real(lk),dimension(:,:) :: wppart,wtpart0,wtpart1
      real(lk),dimension(0:,:) :: density,flow
      real(lk),dimension(:,:) :: zpart
      real(lk),dimension(:,:),optional :: pressurepara,pressureperp,sfluidn
      real(lk),dimension(:,:,:),optional :: dnsave
      character(*),intent(in),optional :: switch
    end subroutine hybridChargeParticle
  end interface

  integer ierror
  character(*),intent(in) :: species_name

  if(species_name=='thermal-ion')then
    !ideal MHD
    if(iload==0)then
      densityi=0.0
      flowi=0.0
    else
      call gkChargeParticle(zion,wpion,wtion0,wtion1,densityi,flowi,jtion0,&
        jtion1,wzion,zonali,zonalci,markeri,markerit,pmarki,qion,aion,&
        iload,ngyroi,mi,switch="w/o density modification")
    endif
  elseif(species_name=='fast-ion')then
    call gkChargeParticle(zfast,wpfast,wtfast0,wtfast1,densityf,flowf,&
      jtfast0,jtfast1,wzfast,zonalf,zonalcf,markerf,markerft,pmarkf,&
      qfast,afast,fload,ngyrof,mf)
  elseif(species_name=='fast-electron')then
    call gkChargeParticle(zfaste,wpfaste,wtfaste0,wtfaste1,densityfe,flowfe,&
      jtfaste0,jtfaste1,wzfaste,zonalfe,zonalcfe,markerfe,markerfet,pmarkfe,&
      qfaste,afaste,feload,ngyrofe,mfe)
  elseif(species_name=='thermal-electron')then
    call hybridChargeParticle(zelectron,wpelectron,wtelectron0,wtelectron1,&
      densitye,flowe,jtelectron0,jtelectron1,wzelectron,zonale,zonalce,&
      markere,markeret,pmarke,qelectron,aelectron,eload,ngyroe,me,&
      ihybrid=ihybrid,nhybrid=nhybrid,pressurepara=pressureepara,&
      pressureperp=pressureeperp,sfluidn=sfluidne,dnsave=dnesave)
  else
    write(*,*)'push.F90: wrong choice'
    call mpi_barrier(mpi_comm_world,ierror)
  endif
end subroutine charge
