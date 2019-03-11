program main
  implicit none
  real(4),allocatable,dimension(:,:,:,:) :: U
  real(4),allocatable,dimension(:) :: shd
  integer i,j,k,nx,ny,nz,val,time,lstime,is,maxshd,timejump
  integer intime
  character(3) timech
  real(4) ,allocatable,dimension(:) :: ttt
  real(4) ,allocatable,dimension(:,:) :: tff,Mass,inttff!,intMass
  real(4)  Lbox,Mass1,tff1,dl,Msun,dt
  real(4), parameter :: G=1.11142e-4

  nx=512
  ny=512
  nz=512
  val=9
  lstime=102
  Lbox=100.e0
  Msun=1.473e-2
  maxshd=1
  timejump=4
  intime=94
  dt=1.e0

  dl=Lbox/real(nx)

  allocate(U(nx,ny,nz,val))
  allocate(tff(lstime,maxshd),Mass(lstime,maxshd),ttt(lstime),inttff(0:lstime,maxshd))
  allocate(shd(maxshd))
  tff(:,:)=0.e0
  Mass(:,:)=0.e0
  inttff(:,:)=0.e0
  !intMass(:,:)=0.e0
  ttt(:)=0.e0

  do i=1,maxshd
     shd(i)=1.e4*real(i)
  end do

  do time=intime,lstime,timejump
     write(*,*) time,'timestep'
     ttt(time)=0.25e0*real(time)
     write(timech,'(i3.3)') time
     open(110,file='Allnewbigtime'//timech//'.DAT',FORM='UNFORMATTED')
     U(:,:,:,:)=0.e0
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

     open(120,file='Allnewbigtimeden'//timech//'.DAT',FORM='UNFORMATTED')
     do k=1,nz
        do j=1,ny
           do i=1,nz
              write(120) U(i,j,k,1)
           end do
        end do
        write(*,*) k, U(5,5,k,2),'read-phyval'
     end do
     close(120)
  end do

  deallocate(U)
  deallocate(tff,Mass,inttff,ttt,shd)
end program main
