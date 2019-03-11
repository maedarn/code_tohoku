program main
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm,time1,timeloop,jump,timestep,cnt10,cc,icn
  integer nx,ny,nz,rgn,zero,sloop,smesh
  integer,allocatable,dimension(:,:,:,:) ::tag
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num,timech
  character(1) num1
  double precision Lbox,Msun,dl,Mshd,ganmim1,Mtot,Mass1,tffr
  integer,allocatable,dimension(:) :: cas,count1,count1S
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1,shd
  double precision ,allocatable,dimension(:) :: MM,S,lgM,phi,tff,Mass,ttt,tffint,Massint
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg
  double precision ,allocatable,dimension(:) :: emag,ekin,ethm,egrv,div
  double precision ,allocatable,dimension(:,:,:) :: jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv
  double precision ,allocatable,dimension(:,:,:) :: vdotf,vdotfg,vdotfp,vdotfb
  integer,allocatable,dimension(:):: lgMcnt
  double precision cs,nsisu,ff,vv,ffb,ffg,ffp,tff1,Massr
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


  allocate(tff(175),Mass(175),ttt(0:175),tffint(175),Massint(175))

  open(910,file='tff.DAT',FORM='FORMATTED')
  do i=1,175
     read(910,*) tff(i),Mass(i),ttt(i)
  end do
  close(910)

  tffint(:)=0.d0
  Massint(:)=0.d0
  tffr=0.d0
  Massr=0.d0
  do i=1,175
     if(tff(i)>0.d0) then
        tffr=tffr+1.d0/tff(i)
        Massr=Massr+1.d0/tff(i)*Mass(i)
        tffint(i)=tffr
        Massint(i)=Massr
     end if
  end do
  open(920,file='tff.DAT',FORM='FORMATTED')
  do i=1,175
     write(920,*) tffint(i),Massint(i),ttt(i),1.d5/Massint(i)
  end do
  close(920)

  deallocate(tff,Mass,ttt,tffint,Massint)
end program main
