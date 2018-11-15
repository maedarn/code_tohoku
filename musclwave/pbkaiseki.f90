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
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname,isat,iend
 character(50) filename,data
 character(4) NPE!,time
 character(6) time

  !******parameter********
  core=64
  mesh=32+2
  val=1
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
   ALLOCATE(u(-1:mesh,-1:mesh,-1:mesh,1:val))
   !ALLOCATE(v(-1:mesh,-1:mesh,-1:mesh,1:val))
   ALLOCATE(v(-1:lengthx+2,1:lengthy,1:lengthz,1:val))

do usecore=0,core-1,1
   IST = mod(usecore,NSPLTx)
   KST = usecore/(NSPLTx*NSPLTy)
   JST = usecore/NSPLTx-NSPLTy*KST
   !KST = mod((k-corestart)/4,NSPLTy)
   write(*,*) IST,JST,KST,usecore
   !write(time,'(I6.6)') times
   write(NPE,'(I4.4)') usecore
   open(unit=150,file='/Users/maeda/Desktop/kaiseki/testcode9/init'//NPE//'.DAT')!,FORM='UNFORMATTED') !バイナリ
   !open(unit=150,file='INIPHI'//NPE//time//'.DAT',FORM='UNFORMATTED') !バイナリ
   !     open(unit=150,file=trim(data)//NPE//'.dat')
   !open(unit=150,file='2D000031.dat',FORM='UNFORMATTED') !バイナリ
   !open(unit=250,file='NEW'//time//NPE//'.dat') !バイナリ
   do i =  -1, mesh
   do k = 1, mesh-2 !*************************************************** +1 いる？ **********************************
      do j = 1, mesh-2
         ! do i =  -1, mesh
         read(150,*) u(i,j,k,1) !,i=-1,mesh)
      end do
   end do
   !write(*,*) 'ok11'
   read (150,'(a)') time
   read (150,'(a)') time
   !read(150,*) time
   !read(150,*) time
   end do
   close(150)
   isat=1
   iend=Ncellx
   if(IST==0)then
      isat = -1
   end if
   if(IST==NSPLTx-1)then
      iend = Ncellx+2
   end if
   do k = 1, Ncellz
      do j = 1, Ncelly
         do i = isat, iend
            !write(*,*)i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),isat, iend,IST,JST,KST,usecore
            v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),1)=u(i,j,k,1)
            !write(*,*)i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz)
         end do
      end do
   end do
end do
open(100,file='/Users/maeda/Desktop/kaiseki/testcode9/INIPHI.DAT')!,FORM='UNFORMATTED')
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
