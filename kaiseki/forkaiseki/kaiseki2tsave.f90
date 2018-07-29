program main
  implicit none
  real(4), allocatable :: u(:,:,:,:)
  double precision , allocatable :: v(:,:,:,:),x(:),y(:),z(:),inputlength(:),t(:)
  integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times,mesh,IST,KST,JST,lengthx,lengthy,lengthz,corexy,midle,corestart,lasttime,length
  integer :: i,j,k,kz=1,tsave,rate,it
  character(50) filename
  character(3) NPE,time
  character(5) timename

  !******parameter********
  core=64
  mesh=64
  val=17
  lasttime=19
  tsave
  rate=4
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

   lengthx=mesh*NSPLTx
   lengthy=mesh*NSPLTy
   lengthz=mesh*NSPLTz
   length=lengthx+lengthy!+lengthz
   ALLOCATE(inputlength(1:length))
   ALLOCATE(x(1:lengthx))
   ALLOCATE(y(1:lengthy))
   ALLOCATE(t(1:tsave*4))

   open(140,file='cdnt.DAT')
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
   close(140)

   !ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   !ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val))
   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val))

   !midle=core-core/NSPLTx/NSPLETy/2
   !corexy=core/NSPLTz
   corexy=NSPLTx*NSPLTy
   midle=NSPLTz/2
   corestart=core-corexy*midle
   !IST = mod(core,NSPLTx); KST = core/(NSPLTx*NSPLTy); JST = core/NSPLTx-NSPLTy*KST

!100 format(D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3)
   open(230,file='time.dat')
   do i=1,tsave
      it=1+i/rate
      read(250,*) t(it)
   end do


   do times=0,lasttime
      do k=corestart,corestart+corexy-1,1
         IST = mod(k,NSPLTx) !;KST = k/(NSPLTx*NSPLTy)
         KST = mod((k-corestart)/4,NSPLTy)
         write(*,*) IST,KST,time,k
         write(time,'(I3.3)') times
         write(NPE,'(I3.3)') k
         open(unit=150,file='2D'//time//NPE//'.dat',FORM='UNFORMATTED') !バイナリ
         !open(unit=150,file='2D000031.dat',FORM='UNFORMATTED') !バイナリ
         open(unit=250,file='NEW'//time//NPE//'.dat') !バイナリ
         do j=1,mesh
            !do i=1,mesh
               !write(*,*) 'ok1'
            read(150) (u(i,j,kz,1),u(i,j,kz,2),u(i,j,kz,3),u(i,j,kz,4),u(i,j,kz,5),u(i,j,kz,6), &
                 u(i,j,kz,7),u(i,j,kz,8), u(i,j,kz,9),u(i,j,kz,10),u(i,j,kz,11),u(i,j,kz,12),u(i,j,kz,13),&
                 u(i,j,kz,14),u(i,j,kz,15),u(i,j,kz,16),u(i,j,kz,17),i=1,mesh)
             ! write(*,*) 'ok2'
            !end do
            write(250,*)  (u(i,j,kz,1),u(i,j,kz,2),u(i,j,kz,3),u(i,j,kz,4),u(i,j,kz,5),u(i,j,kz,6), &
                 u(i,j,kz,7),u(i,j,kz,8), u(i,j,kz,9),u(i,j,kz,10),u(i,j,kz,11),u(i,j,kz,12),u(i,j,kz,13),&
                 u(i,j,kz,14),u(i,j,kz,15),u(i,j,kz,16),u(i,j,kz,17),i=1,mesh)
         end do
         close(250)
         close(150)
         do j=1,mesh
            do i=1,mesh

               v(i+KST*mesh,j+IST*mesh,kz,1)=dble(u(i,j,kz,1))
               v(i+KST*mesh,j+IST*mesh,kz,2)=dble(u(i,j,kz,2))
               v(i+KST*mesh,j+IST*mesh,kz,3)=dble(u(i,j,kz,3))
               v(i+KST*mesh,j+IST*mesh,kz,4)=dble(u(i,j,kz,4))
               v(i+KST*mesh,j+IST*mesh,kz,5)=dble(u(i,j,kz,5))
               v(i+KST*mesh,j+IST*mesh,kz,6)=dble(u(i,j,kz,6))
               v(i+KST*mesh,j+IST*mesh,kz,7)=dble(u(i,j,kz,7))
               v(i+KST*mesh,j+IST*mesh,kz,8)=dble(u(i,j,kz,8))
               v(i+KST*mesh,j+IST*mesh,kz,9)=dble(u(i,j,kz,9))
               v(i+KST*mesh,j+IST*mesh,kz,10)=dble(u(i,j,kz,10))
               v(i+KST*mesh,j+IST*mesh,kz,11)=dble(u(i,j,kz,11))
               v(i+KST*mesh,j+IST*mesh,kz,12)=dble(u(i,j,kz,12))
               v(i+KST*mesh,j+IST*mesh,kz,13)=dble(u(i,j,kz,13))
               v(i+KST*mesh,j+IST*mesh,kz,14)=dble(u(i,j,kz,14))
               v(i+KST*mesh,j+IST*mesh,kz,15)=dble(u(i,j,kz,15))
               v(i+KST*mesh,j+IST*mesh,kz,16)=dble(u(i,j,kz,16))
               v(i+KST*mesh,j+IST*mesh,kz,17)=dble(u(i,j,kz,17))

              ! u(:,:,:,:) =0.0e0
            end do
         end do
      end do

      open(160,file='2Dall'//time//'.DAT',form='formatted')
      !write(timename,'(f5.4)') t(times+1)
      !open(160,file=//timename//'.DAT',form='formatted')
!101   format(D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3)
      do j=1,lengthy
            do i=1,lengthx

               write(160,*) x(i),y(j),v(i,j,kz,1),v(i,j,kz,2),v(i,j,kz,3),v(i,j,kz,4),v(i,j,kz,5),v(i,j,kz,6),v(i,j,kz,7),&
                    v(i,j,kz,8),v(i,j,kz,9),v(i,j,kz,10),v(i,j,kz,11),v(i,j,kz,12),v(i,j,kz,13),v(i,j,kz,14),v(i,j,kz,15), &
                    v(i,j,kz,16),v(i,j,kz,17)

            end do
            write(160,*)
         end do
         close(160)
      end do
    end program main

