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

   open(unit=150,file='/Users/maeda/Desktop/kaiseki/testcode9/INIPHI.DAT')
   open(unit=151,file='/Users/maeda/Desktop/kaiseki/testcode9/INIPHI2.DAT')
   do k = 1,lengthz
      do j = 1,lengthy
         do i =  -1, lengthx+1
            read(150,*) v(i,j,k,1)
            !read(151,*) v(i,j,k,2)
         end do
      end do
   end do
   do k = 1,lengthz
      do j = 1,lengthy
         do i =  -1, lengthx+1
            !read(150,*) v(i,j,k,1)
            read(151,*) v(i,j,k,2)
         end do
      end do
   end do
   close(150)
   close(151)
   do k = 1,lengthz
      do j = 1,lengthy
         do i =  -1, lengthx+1
            v(i,j,k,3)=v(i,j,k,1)-v(i,j,k,2)
         end do
      end do
   end do
   open(100,file='/Users/maeda/Desktop/kaiseki/testcode9/INIPHIdif.DAT')
   do k=1,lengthz
      do j=1,lengthy
         do i=-1,lengthx+1
            write(100,*) v(i,j,k,3)
            if( v(i,j,k,3) .ne. 0.0d0) then
               write(*,*) i,j,k,'.ne.0'
            end if
         end do
      end do
   end do
   close(100)


DEALLOCATE(v)
DEALLOCATE(u)

end program
