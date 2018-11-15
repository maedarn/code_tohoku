program computeraw
 implicit none
 !integer :: nx, ny, nz
 integer(kind=8) :: rec_len
 !integer :: i,j,k
 !parameter(nx = 100, ny = 20, nz = 50)
 !real, dimension(nx,ny,nz) :: array
 real(4), allocatable :: u(:,:,:,:) !, v(:,:,:,:)
 real(4), allocatable :: v(:,:,:,:)
 real(8), allocatable ::  pbr(:,:),pbl(:,:)
 real(8), allocatable ::  pbrall(:,:),pblall(:,:)
 double precision , allocatable :: x(:),y(:),z(:),inputlength(:)
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times,mesh,IST,KST,JST,lengthx,lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname,isat,iend,iii
 integer :: count0=0,count1=1,count2=2,count3=3
 character(50) filename,data
 character(4) NPE,count0c,count1c,count2c,count3c!,time
 character(6) time
 character(2) xxx

  !******parameter********
  core=64
  mesh=32+2
  val=4
  lasttime=300
  initime=0
  Ncellx=32
  Ncelly=32
  Ncellz=32
  !***********************

if(core.eq.4)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 1; end if
if(core.eq.8)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 2; end if
if(core.eq.16)   then; NSPLTx = 4; NSPLTy = 2; NSPLTz = 2; end if
if(core.eq.32)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 2; end if
if(core.eq.64)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 4; end if
if(core.eq.128)  then; NSPLTx = 8; NSPLTy = 4; NSPLTz = 4; end if
if(core.eq.256)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 4; end if
if(core.eq.512)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 8; end if
if(core.eq.1024) then; NSPLTx = 8; NSPLTy = 8; NSPLTz =16; end if

   lengthx=(mesh-2)*NSPLTx
   lengthy=(mesh-2)*NSPLTy
   lengthz=(mesh-2)*NSPLTz

   !ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(u(-1:mesh,1:lengthy,1:lengthz,1:val))
   !ALLOCATE(v(-1:mesh,-1:mesh,-1:mesh,1:val))
   ALLOCATE(v(-1:lengthx+2,1:lengthy,1:lengthz,1:val))

   count0 = 0 +4
   count1 = 1 +4
   count2 = 2 +4
   count3 = 3 +4
do iii = 2,mesh-1

   count0 = count0 +4
   count1 = count1 +4
   count2 = count2 +4
   count3 = count3 +4
   if(count0 > core-1) then
      count0 = count0 - core
   end if
   if(count1 > core-1) then
      count1 = count1 - core
   end if
   if(count2 > core-1) then
      count2 = count2 - core
   end if
   if(count3 > core-1) then
      count3 = count3 - core
   end if
   write(xxx,'(I2.2)') iii
   write(*,*) iii,count0,count1,count2,count3
   write(count0c,'(I4.4)') count0
   write(count1c,'(I4.4)') count1
   write(count2c,'(I4.4)') count2
   write(count3c,'(I4.4)') count3
   open(unit=150,file='/Users/maeda/Desktop/kaiseki/testcode9/final'//xxx//count0c//'.dat')
   open(unit=151,file='/Users/maeda/Desktop/kaiseki/testcode9/final'//xxx//count1c//'.dat')
   open(unit=152,file='/Users/maeda/Desktop/kaiseki/testcode9/final'//xxx//count2c//'.dat')
   open(unit=153,file='/Users/maeda/Desktop/kaiseki/testcode9/final'//xxx//count3c//'.dat')
   do k = 1,lengthz
      do j = 1,lengthy
         ! do i =  -1, mesh
         read(150,*) u(kz,j,k,1) !,i=-1,mesh)
         read(151,*) u(kz,j,k,2)
         read(152,*) u(kz,j,k,3)
         read(153,*) u(kz,j,k,4)
      end do
   end do
   close(150)
   close(151)
   close(152)
   close(153)
   do k = 1,lengthz
      do j = 1,lengthy
         !do i = 0, mesh+1
         v(iii+32*0-1,j,k,1) = u(kz,j,k,1)
         v(iii+32*1-1,j,k,1) = u(kz,j,k,2)
         v(iii+32*2-1,j,k,1) = u(kz,j,k,3)
         v(iii+32*3-1,j,k,1) = u(kz,j,k,4)
         !end do
      end do
   end do
end do


count0 = -4
do iii = 0,1
   count0 = count0 +4
   write(count0c,'(I4.4)') count0
   write(xxx,'(I2.2)') iii
   open(unit=150,file='/Users/maeda/Desktop/kaiseki/testcode9/final'//xxx//count0c//'.dat')
   write(*,*) iii,count0c
   do k = 1,lengthz
      do j = 1,lengthy
         ! do i =  -1, mesh
         read(150,*) u(kz,j,k,1) !,i=-1,mesh)
      end do
   end do
   close(150)
   do k = 1,lengthz
      do j = 1,lengthy
         !do i = 0, mesh+1
         v(iii+32*0-1,j,k,1) = u(kz,j,k,1)
         !end do
      end do
   end do
end do

count3 = 3+4
do iii = mesh,mesh+1
   count3 = count3 +4
   write(count3c,'(I4.4)') count3
   write(xxx,'(I2.2)') iii
   open(unit=153,file='/Users/maeda/Desktop/kaiseki/testcode9/final'//xxx//count3c//'.dat')
   do k = 1,lengthz
      do j = 1,lengthy
         ! do i =  -1, mesh
         read(153,*) u(kz,j,k,4) !,i=-1,mesh)
      end do
   end do
   close(153)
   do k = 1,lengthz
      do j = 1,lengthy
         !do i = 0, mesh+1
         v(iii+32*3-1,j,k,1) = u(kz,j,k,4)
         !end do
      end do
   end do
end do



open(100,file='/Users/maeda/Desktop/kaiseki/testcode9/INIPHI2.DAT')!,FORM='UNFORMATTED')
!open(100,file='INIPHI'//time//'.DAT',FORM='UNFORMATTED')
do k=1,lengthz
   do j=1,lengthy
      do i=-1,lengthx+1
         write(100,*) v(i,j,k,1)
      end do
   end do
   !write(*,*) k,times
end do
close(100)


DEALLOCATE(v)
DEALLOCATE(u)

end program
