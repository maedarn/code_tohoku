program main
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm,time1,timeloop,jump,timestep,cnt10,cc,icn
  integer nx,ny,nz,rgn,zero,sloop,smesh
  integer,allocatable,dimension(:,:,:,:) ::tag
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num,timech
  character(1) num1
  double precision Lbox,Msun,dl,Mshd,ganmim1,Mtot,Mass1
  integer,allocatable,dimension(:) :: cas,count1,count1S
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1,shd
  double precision ,allocatable,dimension(:) :: MM,S,lgM,phi,tff,Mass,ttt
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg
  double precision ,allocatable,dimension(:) :: emag,ekin,ethm,egrv,div
  double precision ,allocatable,dimension(:,:,:) :: jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv
  double precision ,allocatable,dimension(:,:,:) :: vdotf,vdotfg,vdotfp,vdotfb
  integer,allocatable,dimension(:):: lgMcnt
  double precision cs,nsisu,ff,vv,ffb,ffg,ffp,tff1
  DOUBLE PRECISION, parameter :: G=1.11142d-4

  !----parameter----
  nx=512
  ny=512
  nz=512
  val=9
  loop=5
  Lbox=100.d0
  Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
  rgn=1000
  Mshd=1.d2
  ganmim1=3.d0/2.d0 !1.d0/((5.d0/3.d0-1.d0)
  sloop=30
  smesh=5
  timeloop=7
  timestep=4
  jump=4
  !----parameter----

  dl=Lbox/dble(nx)

  allocate(U(-1:nx+2,-1:ny+2,-1:nz+2,val))
  !allocate(tag(1:nx,1:ny,1:nz,loop))
  allocate(tag(0:nx+1,0:ny+1,0:nz+1,loop))
  allocate(lgM(0:sloop),lgMcnt(sloop))
  allocate(x(nx),y(ny),z(nz),xg(rgn),yg(rgn),zg(rgn),MM(rgn),S(rgn))
  allocate(vx(rgn),vy(rgn),vz(rgn),vxg(rgn),vyg(rgn),vzg(rgn))
  allocate(emag(rgn),ekin(rgn),ethm(rgn),egrv(rgn))
  allocate(cas(0:rgn),count1(rgn),div(rgn),ig(rgn),jg(rgn),kg(rgn),count1S(rgn),phi(rgn))
  allocate(jx(1:nx,1:ny,1:nz),jy(1:nx,1:ny,1:nz),jz(1:nx,1:ny,1:nz),fx(1:nx,1:ny,1:nz),fy(1:nx,1:ny,1:nz),fz(1:nx,1:ny,1:nz)&
       ,fgx(1:nx,1:ny,1:nz),fgy(1:nx,1:ny,1:nz),fgz(1:nx,1:ny,1:nz),fpx(1:nx,1:ny,1:nz),fpy(1:nx,1:ny,1:nz),fpz(1:nx,1:ny,1:nz)&
       ,fbx(1:nx,1:ny,1:nz),fby(1:nx,1:ny,1:nz),fbz(1:nx,1:ny,1:nz),divv(1:nx,1:ny,1:nz))
  allocate(vdotf(1:nx,1:ny,1:nz),vdotfg(1:nx,1:ny,1:nz),vdotfp(1:nx,1:ny,1:nz),vdotfb(1:nx,1:ny,1:nz))
  allocate(tff(175),Mass(175),ttt(0:175))

  !U(:,:,:,:)=0.e0
  tag(:,:,:,:)=0



 !do time1=1,timeloop,timestep
  !time1=175
  open(unit=350,file='cnt3.dat')
  read(350,*) cnt10
  close(350)
  if(cnt10==1)then
     tff(:)=0.d0
     Mass(:)=0.d0
     ttt(:)=0.d0
     ttt(0)=-1.d0
  end if
  if(cnt10.ne.1)then
     open(920,file='tff.DAT',FORM='FORMATTED')
     do i=1,175
     read(920,*) tff(i),Mass(i),ttt(i)
     end do
     close(920)
  end if
  icn=(cnt10-1)/timestep+1
  !write(cntc,'(i3.3)') cnt10
  write(timech,'(i3.3)') cnt10
  cnt10=cnt10+timestep
  open(unit=350,file='cnt3.dat')
  write(350,*) cnt10
  close(350)
  !write(timech,'(i3.3)') time1

  U(:,:,:,:)=0.e0
  open(110,file='Allnewbigtime'//timech//'.DAT',FORM='UNFORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nz
           read(110) U(i,j,k,1), U(i,j,k,2), U(i,j,k,3), U(i,j,k,4), U(i,j,k,5),&
                U(i,j,k,6), U(i,j,k,7), U(i,j,k,8), U(i,j,k,9)
        end do
     end do
     write(*,*) k, U(5,5,k,2),'read-phyval'
  end do
  close(110)

  shd=1.e3
  cc=0
  tff1=0.d0
  Mass1=0.d0

  do k=1,nz
     do j=1,ny
        do i=1,nz
           if(U(i,j,k,1)>shd)then
              !tff1=1/dsqrt(G*dble(U(i,j,k,1)))+tff1
              tff1=tff1+dble(U(i,j,k,1))/dsqrt(G*dble(U(i,j,k,1)))*dl*dl*dl
              cc=cc+1
              Mass1=dble(U(i,j,k,1))*dl*dl*dl+Mass1
           end if
        end do
     end do
     write(*,*) k,'caltff'
  end do
  if(cc>0) then
   tff1=tff1/dble(cc)
   Mass1=Mass1*Msun
   tff(icn)=tff1
   Mass(icn)=Mass1
 end if
!end do

  ttt(icn)=ttt(icn-1)+1.d0
  open(910,file='tff.DAT',FORM='FORMATTED')
  do i=1,175
     write(910,*) tff(i),Mass(i),ttt(i)
  end do
  close(910)

  deallocate(U,lgM,lgMcnt,phi)
  deallocate(tag,MM,S)
  deallocate(x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg)
  deallocate(emag,ekin,ethm,egrv)
  deallocate(cas,count1,div,ig,jg,kg,count1S)
  deallocate(jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv)
  deallocate(vdotf,vdotfg,vdotfp,vdotfb,tff,Mass,ttt)
end program mai
