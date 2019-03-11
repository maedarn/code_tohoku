module comvar
  integer parameter :: nx=512,ny=512,nz=512,val=9
  real(4),dimension(nx,ny,nz,val) ::U
  real(4),dimension(nx/2,ny/2,nz/2,val) ::Un
end module comvar
program main
  use comvar
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm,time1,timeloop,jump,timestep,cnt10
  integer nx,ny,nz,rgn,zero,sloop,smesh,smg
  integer,allocatable,dimension(:,:,:,:) ::tag
  !real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num,timech
  character(1) num1
  character(9) ch1
  character(9) ch2
  double precision Lbox,Msun,dl,Mshd,ganmim1,phi,Mtot
  integer,allocatable,dimension(:) :: cas,count1,count1S
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1
  double precision ,allocatable,dimension(:) :: MM,S,lgM
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg
  double precision ,allocatable,dimension(:) :: emag,ekin,ethm,egrv,div
  double precision ,allocatable,dimension(:,:,:) :: jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv
  double precision ,allocatable,dimension(:,:,:) :: vdotf,vdotfg,vdotfp,vdotfb
  integer,allocatable,dimension(:):: lgMcnt
  double precision cs,nsisu,ff,vv,ffb,ffg,ffp

  !----parameter----
  !nx=512
  !ny=512
  !nz=512
  !val=9
  loop=1
  Lbox=100.d0
  Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
  rgn=1000
  Mshd=1.d2
  ganmim1=3.d0/2.d0 !1.d0/((5.d0/3.d0-1.d0)
  sloop=20
  smesh=5
  timeloop=7
  timestep=15
  jump=4
  smg=2
  ch1='Allnewbig'
  ch2='Allbigsm1'
  !----parameter----

  dl=Lbox/dble(nx)

  !allocate(U(1:nx,1:ny,1:nz,val))
  !allocate(Un(1:nx/2,1:ny/2,1:nz/2,val))
  allocate(tag(1:nx,1:ny,1:nz,loop))
  allocate(lgM(0:sloop),lgMcnt(sloop))
  allocate(x(nx),y(ny),z(nz),xg(rgn),yg(rgn),zg(rgn),MM(rgn),S(rgn))
  allocate(vx(rgn),vy(rgn),vz(rgn),vxg(rgn),vyg(rgn),vzg(rgn))
  allocate(emag(rgn),ekin(rgn),ethm(rgn),egrv(rgn))
  allocate(cas(0:rgn),count1(rgn),div(rgn),ig(rgn),jg(rgn),kg(rgn),count1S(rgn))
  allocate(jx(1:nx,1:ny,1:nz),jy(1:nx,1:ny,1:nz),jz(1:nx,1:ny,1:nz),fx(1:nx,1:ny,1:nz),fy(1:nx,1:ny,1:nz),fz(1:nx,1:ny,1:nz)&
       ,fgx(1:nx,1:ny,1:nz),fgy(1:nx,1:ny,1:nz),fgz(1:nx,1:ny,1:nz),fpx(1:nx,1:ny,1:nz),fpy(1:nx,1:ny,1:nz),fpz(1:nx,1:ny,1:nz)&
       ,fbx(1:nx,1:ny,1:nz),fby(1:nx,1:ny,1:nz),fbz(1:nx,1:ny,1:nz),divv(1:nx,1:ny,1:nz))
  allocate(vdotf(1:nx,1:ny,1:nz),vdotfg(1:nx,1:ny,1:nz),vdotfp(1:nx,1:ny,1:nz),vdotfb(1:nx,1:ny,1:nz))

  !U(:,:,:,:)=0.e0
!  tag(:,:,:,:)=0


!  nsisu=2
!  lgM(:)=0.d0
!  lgM(0)=10.0d0**nsisu
!  do j=1,sloop
!     nsisu=2
!     nsisu =nsisu + 1.d0/dble(smesh)*dble(j)
!     lgM(j) = 10.0d0**nsisu
!     write(*,*) lgM(j)
!  end do


  do time1=1,timeloop
     open(unit=350,file='cnt5.dat')
     read(350,*) cnt10
     close(350)
     !write(cntc,'(i3.3)') cnt10
     write(timech,'(i3.3)') cnt10
     cnt10=cnt10+timestep
     open(unit=350,file='cnt5.dat')
     write(350,*) cnt10
     close(350)
     !write(timech,'(i3.3)') time1
     U(:,:,:,:)=0.e0
     Un(:,:,:,:)=0.e0
  open(110,file=ch1//timech//'.DAT',FORM='UNFORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nz
           read(110) U(i,j,k,1), U(i,j,k,2), U(i,j,k,3), U(i,j,k,4), U(i,j,k,5),&
                U(i,j,k,6), U(i,j,k,7), U(i,j,k,8), U(i,j,k,9)
        end do
     end do
     write(*,*) k, U(5,5,5,2),'read-phyval'
  end do
  close(110)

  call sm(U(1,1,1,1),Un(1,1,1,1),nx,ny,nz)
  call sm(U(1,1,1,9),Un(1,1,1,9),nx,ny,nz)

open(100,file=ch2//timech//'.DAT',FORM='UNFORMATTED')
do k=1,nz/2
   do j=1,ny/2
      do i=1,nx/2
         write(100,*) Un(i,j,k,1),Un(i,j,k,9)
      end do
   end do
   write(*,*) k
end do
close(100)

end do
  deallocate(lgM,lgMcnt)
  deallocate(tag,MM,S)
  deallocate(x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg)
  deallocate(emag,ekin,ethm,egrv)
  deallocate(cas,count1,div,ig,jg,kg,count1S)
  deallocate(jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv)
  deallocate(vdotf,vdotfg,vdotfp,vdotfb)

end program main

subroutine sm(U1,Un1,nx1,ny1,nz1)
  !real(4),allocatable,dimension(:,:,:) :: U1
  !real(4),allocatable,dimension(:,:,:) :: Un1
  real(4,dimension(nx1,ny1,nz1) :: U1
  real(4),dimension(nx1/2,ny1/2,nz1/2) :: Un1
  integer nx1,ny1,nz1,i,j,k

  !allocate(U(1:nx,1:ny,1:nz,val))
  !allocate(Un(1:nx/2,1:ny/2,1:nz/2,val))
  do k=1,nz/2
     do j=1,ny/2
        do i=1,nz/2
           Un(i,j,k)=(U(2*i-1,2*j-1,2*k-1)+U(2*i,2*j-1,2*k-1)+U(2*i,2*j,2*k-1)+U(2*i,2*j-1,2*k)&
                +U(2*i-1,2*j,2*k-1)+U(2*i-1,2*j-1,2*k)+U(2*i-1,2*j,2*k)+U(2*i,2*j,2*k))/8.d0
        end do
     end do
     !write(*,*) k, U(5,5,5,2),'read-phyval'
  end do
end subroutine sm
