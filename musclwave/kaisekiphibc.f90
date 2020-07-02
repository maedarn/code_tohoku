program computeraw
 implicit none
 !integer :: nx, ny, nz
 integer(kind=8) :: rec_len
 integer cnt10,jump,ts,timejump
 !integer :: i,j,k
 !parameter(nx = 100, ny = 20, nz = 50)
 !real, dimension(nx,ny,nz) :: array
 real(8), allocatable :: u(:,:,:) !, v(:,:,:,:)
 real(8), allocatable ::  v(:,:,:)
 double precision , allocatable :: x(:),y(:),z(:),inputlength(:)
 real(4) mass,nth
 real(8) mass1
 integer(4), allocatable :: x1shock(:,:),x2shock(:,:)
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times=0,mesh,IST,KST,JST,lengthx,&
      lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname
 character(50) filename,data
 character(3) NPE,time,cntctj
 character(6) cntc
 character(34) :: dir='/glv0/maedarn/test-grvwave/PHIINI/'
 
  !******parameter********
  core=64
  mesh=32+2
  val=2
  jump=1
  ts=1
  initime=0
  lasttime=2
  Ncellx=32
  Ncelly=32
  Ncellz=32
  timejump=1
  nth=10.e0
  !***********************

  mass=0.e0
  mass1=0.d0
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

   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(Tn(1:lengthx,1:lengthy,1:lengthz))
   ALLOCATE(u(-1:mesh,-1:mesh,-1:mesh,1:val))
   ALLOCATE(x1shock(1:lengthy,1:lengthz),x2shock(1:lengthy,1:lengthx))

  
write(*,*)'1',initime,lasttime,timejump
do times=initime,lasttime,timejump
     write(*,*)times
       write(cntc,'(i1.1)') times
       do usecore=0,core-1,1
         IST = mod(usecore,NSPLTx)
         KST = usecore/(NSPLTx*NSPLTy)
         JST = usecore/NSPLTx-NSPLTy*KST
         write(*,*) IST,JST,KST,time,usecore
         write(NPE,'(I3.3)') usecore/NSPLTx
         open(unit=150,file=dir//'bcsave'//cntc//NPE//'.DAT',FORM='UNFORMATTED')
         do k = -1, Ncellz+2
            do j = -1, Ncelly+2
                  read(150) u(j,k,1) ,u(j,k,2)
            end do
         end do
         close(150)
         do k = 1, Ncellz
            do j = 1, Ncelly
                  v(j+JST*(mesh-2),k+KST*(mesh-2),1)=u(i,j,k,1)
                  v(j+JST*(mesh-2),k+KST*(mesh-2),2)=u(i,j,k,2)
         end do
!write(*,*)'Cmetal',u(k,30,30,14),u(k,30,30,15)
      end do
   end do

   write(cntctj,'(i3.3)') times/timejump
   !open(100,file=dir//'All/All'//cntctj//'.DAT',access='stream',FORM='UNFORMATTED')
   open(100,file=dir//'Allbc'//cntc//'.DAT',FORM='FORMATTED')
   do k=1,lengthz,jump
      do j=1,lengthy,jump
         write(100,*) v(j,k,1) ,v(j,k,2)
      end do
      write(*,*) k,u(64,64,64,1)!,data
      write(100,*)
   end do
   close(100)
end do


DEALLOCATE(v)
DEALLOCATE(u)
DEALLOCATE(x1shock,x2shock)

end program computeraw
