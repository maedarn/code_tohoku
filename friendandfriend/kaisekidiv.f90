program main
  implicit none
  integer i,j,k,val,m,loop
  integer nx,ny,nz
  integer,allocatable,dimension(:,:,:,:) ::tag
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num
  double precision Lbox,M

  !----parameter----
  nx=512
  ny=512
  nz=512
  val=8
  loop=5
  Lbox=100.d0
  !----parameter----

  allocate(U(-1:nx+2,,-1:ny+2,-1:nz+2,val))
  allocate(tag(1:nx,,1:ny,1:nz,loop))

  U(:,:,:,:)=0.e0
  tag(:,:,:)=0

  do m=1,loop
     write(num,'(i3.3)') m
     open(150,file='/Users/maeda/Desktop/kaiseki/image/tag'//num//'.DAT',FORM='FORMATTED')
     do k=1,nz
        do j=1,ny
           do i=1,nz
              read(150,*) tag(i,j,k,m)
           end do
        end do
        write(*,*) k,m,'read-tag'
     end do
     close(150)
  end do

  open(110,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/all.DAT',FORM='FORMATTED')
  do k=1,nz
     do j=1,ny
        do i=1,nz
           read(110,*) U(i,j,k,1), U(i,j,k,2), U(i,j,k,3), U(i,j,k,4), U(i,j,k,5),&
                U(i,j,k,6), U(i,j,k,7), U(i,j,k,8)
        end do
     end do
     write(*,*) k,'read-phyval'
  end do
  close(110)

  
end program main
