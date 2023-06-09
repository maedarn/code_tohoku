program computeraw
 implicit none
 !integer :: nx, ny, nz
 integer(kind=8) :: rec_len
 !integer :: i,j,k
 !parameter(nx = 100, ny = 20, nz = 50)
 !real, dimension(nx,ny,nz) :: array
 real(4), allocatable :: u(:,:,:,:) !, v(:,:,:,:)
 real(4), allocatable ::  v(:,:,:,:)
 real(8), allocatable ::  pbr(:,:),pbl(:,:)
 real(8), allocatable ::  pbrall(:,:),pblall(:,:)
 double precision , allocatable :: x(:),y(:),z(:),inputlength(:)
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times,mesh,IST,KST,JST,lengthx,lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname
 character(50) filename,data
 character(3) NPE!,time
 character(6) time

  !******parameter********
  core=64
  mesh=16+2
  val=5
  lasttime=300
  initime=0
  Ncellx=16
  Ncelly=16
  Ncellz=16
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

   !lengthx=mesh*NSPLTx
   !lengthy=mesh*NSPLTy
   !lengthz=mesh*NSPLTz
   length=lengthx+lengthy+lengthz
   ALLOCATE(inputlength(1:length))
   ALLOCATE(x(1:lengthx))
   ALLOCATE(y(1:lengthy))
   ALLOCATE(z(1:lengthz))


!400 format(D25.17)
!   open(140,file='cdnt.DAT')
!   read(140,400) ( x(i) , i=1, Ncellx*NSPLTx )
!   Read(140,400) ( y(j) , j=1, Ncelly*NSPLTy )
!   read(140,400) ( z(k) , k=1, Ncellz*NSPLTz )

   goto 2090
   do i=1,length
      !read(140,'(D10.3)') inputlength(i)
      read(140,*) inputlength(i)
   end do
   do i=1,lengthx
      x(i)=inputlength(i)
   end do
   do i=1,lengthy
      y(i)=inputlength(i+lengthx)
   end do
   do i=1,lengthz
      z(i)=inputlength(i+lengthx+lengthy)
   end do
   2090 continue

   close(140)

   lengthx=(mesh-2)*NSPLTx
   lengthy=(mesh-2)*NSPLTy
   lengthz=(mesh-2)*NSPLTz

   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(u(-1:mesh,-1:mesh,-1:mesh,1:val))
   ALLOCATE(pbr(0:mesh-1+2,0:mesh-1+2))
   ALLOCATE(pbl(0:mesh-1+2,0:mesh-1+2))
   ALLOCATE(pbrall(1:lengthy,1:lengthz))
   ALLOCATE(pblall(1:lengthy,1:lengthz))


   !   goto 800




   do times=initime,lasttime

      write(time,'(I6.6)') times

      open(100,file='/Users/maeda/Desktop/kaiseki/testcode17/INIPHI'//time//'.dat',FORM='UNFORMATTED')
      do k=1,lengthz
         do j=1,lengthy
            do i=1,lengthx
               read(100) v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5)
            end do
         end do
         write(*,*) k,times
      end do
      close(100)
      open(110,file='/Users/maeda/Desktop/kaiseki/testcode17/INIPHI2'//time//'.DAT',FORM='FORMATTED')
      do k=1,lengthz
         do j=1,lengthy
            do i=1,lengthx
               write(110,*) v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5)
            end do
         end do
         write(*,*) k,times
      end do
      close(110)
   end do

DEALLOCATE(inputlength)
DEALLOCATE(x)
DEALLOCATE(y)
DEALLOCATE(z)

DEALLOCATE(v)
DEALLOCATE(u)

end program computeraw
