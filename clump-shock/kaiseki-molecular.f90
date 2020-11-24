program main
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm,time1,timeloop,jump,timestep,cnt10,ijk,ijkl,loopcem,icem,icn
  integer nx,ny,nz,rgn,zero,sloop,smesh,initime,iii,im,ixpi,jypj,kzpk,ixmi,jymj,kzmk,ixdi,iydi,izdi
  integer,allocatable,dimension(:,:,:,:) ::tag
  !integer,allocatable,dimension(:,:,:,:,:) ::tagr
  integer,allocatable,dimension(:,:,:,:) ::tagr
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num,timech
  character(1) num1
  double precision Lbox,Msun,dl,Mshd,ganmim1,Mtot,Mtotcl,mu,Mh2tot,Mtotushd,Mh2totushd,Mtotclushd
  integer,allocatable,dimension(:) :: cas,count1,count1S
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1
  double precision ,allocatable,dimension(:) :: MM,S,lgM,phi,tff1,vdsp,vdspx,vdspy,vdspz,csmean,csmean1
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg,rshn,thhcem
  double precision ,allocatable,dimension(:) ::emag,ekin,ethm,egrv,div,Tmean,mrho,Mh2,Mh
  double precision ,allocatable,dimension(:,:) :: Massr,rrr
  double precision ,allocatable,dimension(:,:,:) :: jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv,Tn
  double precision ,allocatable,dimension(:,:,:) :: vdotf,vdotfg,vdotfp,vdotfb,Ncem
  integer,allocatable,dimension(:):: lgMcnt
  double precision cs,nsisu,ff,vv,ffb,ffg,ffp,lx,ly,lz,Mstarcum,Mh2totall,Mstarcumushd,eloc,Meloc,Mhfra,Melochfra
  double precision :: kb=8.63359d0,h2fra,Melocstar,Mhfrastar,Melochfrastar,nth100,Mth100,Mtot10,Mh2totall10,Mtot100,Mh2totall100
  DOUBLE PRECISION, parameter :: G=1.11142d-4
  !character(46) :: dir='/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm/'
character(54) :: dir='/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-hiden10/'
!character(55) :: dir='/glv0/maedarn/clst-form-HIcol/cnv100-wb-wg-sm-optthick/'
!  character(53) :: dir='/glv0/maedarn/clst-form-HIcol/cnv100-wsb-wg-sm-200pc/'

  !----parameter----
  nx=512
  ny=512
  nz=512
  val=18
  loop=1
  Lbox=100.d0
 ! Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
  Msun=2.4d-2 !1pc * 1pc * 1pc * 1m_p/cc
  rgn=300
  Mshd=1.d2
  ganmim1=3.d0/2.d0 !1.d0/((5.d0/3.d0-1.d0)
  sloop=30
  smesh=5
  initime=8
  timeloop=8
  !initime=98/4
  !timeloop=98/4
  timestep=1
  jump=1
  mu=1.27d0
  loopcem=5
  icn=1
  nth100=100.0*mu
  !----parameter----

  dl=Lbox/dble(nx)

  allocate(U(-1:nx+2,-1:ny+2,-1:nz+2,val))
  !allocate(tag(1:nx,1:ny,1:nz,loop))
  allocate(tag(0:nx+1,0:ny+1,0:nz+1,loop))
  allocate(lgM(0:sloop),lgMcnt(sloop))
  allocate(x(nx),y(ny),z(nz),xg(0:rgn),yg(0:rgn),zg(0:rgn),MM(0:rgn),S(0:rgn),tff1(0:rgn),vdsp(0:rgn),&
       vdspx(0:rgn),vdspy(0:rgn),vdspz(0:rgn))
  allocate(vx(0:rgn),vy(0:rgn),vz(0:rgn),vxg(0:rgn),vyg(0:rgn),vzg(0:rgn),csmean(0:rgn),csmean1(0:rgn))
  allocate(emag(0:rgn),ekin(0:rgn),ethm(0:rgn),egrv(0:rgn))
  allocate(cas(0:rgn),count1(0:rgn),div(0:rgn),ig(0:rgn),jg(0:rgn),kg(0:rgn),count1S(0:rgn),phi(0:rgn),Tmean(0:rgn))
  allocate(jx(1:nx,1:ny,1:nz),jy(1:nx,1:ny,1:nz),jz(1:nx,1:ny,1:nz),fx(1:nx,1:ny,1:nz),fy(1:nx,1:ny,1:nz),fz(1:nx,1:ny,1:nz)&
       ,fgx(1:nx,1:ny,1:nz),fgy(1:nx,1:ny,1:nz),fgz(1:nx,1:ny,1:nz),fpx(1:nx,1:ny,1:nz),fpy(1:nx,1:ny,1:nz),fpz(1:nx,1:ny,1:nz)&
       ,fbx(1:nx,1:ny,1:nz),fby(1:nx,1:ny,1:nz),fbz(1:nx,1:ny,1:nz),divv(1:nx,1:ny,1:nz),Tn(1:nx,1:ny,1:nz))
  allocate(vdotf(1:nx,1:ny,1:nz),vdotfg(1:nx,1:ny,1:nz),vdotfp(1:nx,1:ny,1:nz),vdotfb(1:nx,1:ny,1:nz))
  allocate(mrho(0:rgn),Mh2(0:rgn),Mh(0:rgn))
  allocate(Ncem(0:8,0:rgn,loopcem),thhcem(loopcem))

!-------------------
  allocate(rshn(5))
  allocate(Massr(0:rgn,5))
  allocate(rrr(0:rgn,3))
  !allocate(tagr(1:nx,1:ny,1:nz,rgn,5))
  !tagr(:,:,:,:,:)=0
  allocate(tagr(1:nx,1:ny,1:nz,5))
  tagr(:,:,:,:)=0
  Massr(:,:)=0.d0
  rrr(:,:)=0.d0
  rshn(1)=1.d0
  rshn(2)=2.d0
  rshn(3)=3.d0
  rshn(4)=4.d0
  rshn(5)=5.d0
  thhcem(1)=1.d1
  thhcem(2)=1.d2
  thhcem(3)=1.d3
  thhcem(4)=1.d4
  thhcem(5)=1.d5
!-------------------
  !U(:,:,:,:)=0.e0
  tag(:,:,:,:)=0



  do time1=initime,timeloop,timestep
     write(timech,'(i3.3)') time1
     

     U(:,:,:,:)=0.e0
     !open(110,file='/home/maedarn/cnv100wbg-cluster/Allhirshighden098.DAT',FORM='FORMATTED')
     open(110,file=dir//'All/All'//timech//'.DAT',access='stream',FORM='UNFORMATTED')
     do k=1,nz
        do j=1,ny
           do i=1,nx
              read(110) (U(i,j,k,n),n=1,val)
           end do
        end do
        write(*,*) k, U(5,5,k,2),'read-phyval'
     end do
     close(110)


        write(num,'(i3.3)') m
        Mtot10=0.d0
        Mh2totall10=0.d0
        Mtot100=0.d0
        Mh2totall100=0.d0
        do k=1,nz
           do j=1,ny
              do i=1,nz
                 if(U(i,j,k,1)>10.d0*mu) then
                 Mtot10=Mtot10+dble(U(i,j,k,1))*(dl**3)
                 Mh2totall10=Mh2totall10+dble(U(i,j,k,11))*dl*dl*dl*2.d0
                 end if
                 if(U(i,j,k,1)>100.d0*mu) then
                 Mtot100=Mtot100+dble(U(i,j,k,1))*(dl**3)
                 Mh2totall100=Mh2totall100+dble(U(i,j,k,11))*dl*dl*dl*2.d0
                 end if
              end do
           end do
           write(*,*) k,n,'mass'
        end do
        !tff1(n)=tff1(n)*Msun
        !Massr(:,:)=Massr(:,:)*Msun


        !end do
        Mtot10=Mtot10*Msun
        Mh2totall10=Mh2totall10*Msun
        Mtot100=Mtot100*Msun
        Mh2totall100=Mh2totall100*Msun
  
    
        open(850,file=dir//'Clmpanly/molecular.DAT',FORM='FORMATTED',position='append')
        write(850,*) dble(time1)*1.d0/dble(timestep)+dble(initime-1)*0.25 ,Mtot10,Mh2totall10,Mtot100,Mh2totall100
        close(850)

     write(*,*)'kokoha?'
  end do

  write(*,*)
!  deallocate(U,lgM,lgMcnt,phi,tff1,vdsp,vdspx,vdspy,vdspz,csmean,csmean1)
!  deallocate(tag,MM,S)
!  deallocate(x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg)
!  deallocate(emag,ekin,ethm,egrv)
!  deallocate(cas,count1,div,ig,jg,kg,count1S)
!  deallocate(jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv,Tn,Tmean)
!  deallocate(vdotf,vdotfg,vdotfp,vdotfb)

  deallocate(U)
  !deallocate(tag(1:nx,1:ny,1:nz,loop))
  deallocate(tag)
  deallocate(lgM,lgMcnt)
  deallocate(x,y,z,xg,yg,zg,MM,S,tff1,vdsp,vdspx,vdspy,vdspz)
  deallocate(vx,vy,vz,vxg,vyg,vzg,csmean,csmean1)
  deallocate(emag,ekin,ethm,egrv)
  deallocate(cas,count1,div,ig,jg,kg,count1S,phi,Tmean)
  deallocate(jx,jy,jz,fx,fy,fz&
       ,fgx,fgy,fgz,fpx,fpy,fpz&
       ,fbx,fby,fbz,divv,Tn)
  deallocate(vdotf,vdotfg,vdotfp,vdotfb)
  deallocate(mrho)


  deallocate(rshn)
  deallocate(Massr,Mh2,Mh)
  deallocate(rrr)
  deallocate(tagr,Ncem,thhcem)
  !U(:,:,:,:)=0.e0
  !tag(:,:,:,:)=0


end program main
