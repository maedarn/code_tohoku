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
 character(3) NPE,time

  !******parameter********
  core=512
  mesh=64+1
  val=17
  !lasttime=800
  !initime=800
  Ncellx=64+1
  Ncelly=64+1
  Ncellz=64+1
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

   lengthx=(mesh-1)*NSPLTx
   lengthy=(mesh-1)*NSPLTy
   lengthz=(mesh-1)*NSPLTz

   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val))
   ALLOCATE(pbr(0:mesh-1+2,0:mesh-1+2))
   ALLOCATE(pbl(0:mesh-1+2,0:mesh-1+2))
   ALLOCATE(pbrall(1:lengthy,1:lengthz))
   ALLOCATE(pblall(1:lengthy,1:lengthz))

   !lengthx=mesh*NSPLTx
   !lengthy=mesh*NSPLTy
   !lengthz=mesh*NSPLTz

!   midle=core-core/NSPLTx/NSPLETy/2
!   corexy=core/NSPLTz
!   corexy=NSPLTx*NSPLTy
!   midle=NSPLTz/2
!   corestart=core-corexy*midle

   !   goto 800

   !do bptime=1,2
   !   if(bptime==1) then
   !      ALLOCATE(pbr(0:mesh-1+2,0:mesh-1+2))
   !      ALLOCATE(pbl(0:mesh-1+2,0:mesh-1+2))
   !      ALLOCATE(pbrall(1:lengthy,1:lengthz))
   !      ALLOCATE(pblall(1:lengthy,1:lengthz))
   !   end if
   !   if(bptime==2) then
   !      mesh=16
   !      ALLOCATE(pbr(1:mesh-1+2,1:mesh-1+2))
   !      ALLOCATE(pbl(1:mesh-1+2,1:mesh-1+2))
   !      ALLOCATE(pbrall(1:lengthy,1:lengthz))
   !      ALLOCATE(pblall(1:lengthy,1:lengthz))
   !   end if
   do usecore=0,core-1,1
      IST = mod(usecore,NSPLTx)
      KST = usecore/(NSPLTx*NSPLTy)
      JST = usecore/NSPLTx-NSPLTy*KST
      !KST = mod((k-corestart)/4,NSPLTy)
      write(*,*) IST,JST,KST,time,usecore
      !      write(time,'(I3.3)') times
      write(NPE,'(I3.3)') usecore
      if(IST==7) then
         open(unit=150,file='bpl'//NPE//'.DAT')

         do k = 0, Ncellz
            do j = 0, Ncelly
               read(150,*) pbl(j,k)
            end do
         end do
         close(150)

         do k = 1, Ncellz-1
            do j = 1, Ncelly-1
               pblall(j+JST*(mesh-1),k+KST*(mesh-1))=pbl(j,k)
            end do
         end do
      end if
      if(IST==NSPLTx-1) then
         open(unit=160,file='bpr'//NPE//'.DAT')
         do k = 0, Ncellz
            do j = 0, Ncelly
               read(160,*) pbr(j,k)
            end do
         end do
         close(160)

         do k = 1, Ncellz-1
            do j = 1, Ncelly-1
               pbrall(j+JST*(mesh-1),k+KST*(mesh-1))=pbr(j,k)
            end do
         end do
      end if
   end do
    !  800 continue

   open(unit=250,file='bpl.dat')
   open(unit=350,file='bpr.dat')

   do k = 1,lengthz
      do j = 1,lengthy
         write(250,*) dble(j),dble(k),pblall(j,k)
      end do
      write(250,*)
   end do
   do k = 1,lengthz
      do j = 1,lengthy
         write(350,*) dble(j),dble(k),pbrall(j,k)
      end do
      write(350,*)
   end do
   close(250)
   close(350)

   DEALLOCATE(pbr)
   DEALLOCATE(pbl)
   DEALLOCATE(pbrall)
   DEALLOCATE(pblall)

   !800 continue
!end do
!goto 870
   do dataname=1,5
   !do dataname=1,1
      if(dataname==1) then
         data='final'
      end if
      if(dataname==2) then
         data='exact'
      end if
      if(dataname==3) then
         data='d'
      end if
      if(dataname==4) then
         data='dpr'
      end if
      if(dataname==5) then
         data='INIT'
      end if
   !do times=initime,lasttime
      do usecore=0,core-1,1
         IST = mod(usecore,NSPLTx)
         KST = usecore/(NSPLTx*NSPLTy)
         JST = usecore/NSPLTx-NSPLTy*KST
         !KST = mod((k-corestart)/4,NSPLTy)
         write(*,*) IST,JST,KST,time,usecore
   !      write(time,'(I3.3)') times
         write(NPE,'(I3.3)') usecore
   !      open(unit=150,file=time//NPE//'.dat',FORM='UNFORMATTED') !バイナリ
         open(unit=150,file=trim(data)//NPE//'.dat')
         !open(unit=150,file='2D000031.dat',FORM='UNFORMATTED') !バイナリ
         !open(unit=250,file='NEW'//time//NPE//'.dat') !バイナリ
         do k = 1, Ncellz !*************************************************** +1 いる？ **********************************
            do j = 1, Ncelly
               do i= 1, Ncellx
                  read(150,*) u(i,j,k,1)! &!,u(i,j,k,2),u(i,j,k,3),u(i,j,k,4),u(i,j,k,5), &
                    !u(i,j,k,6),u(i,j,k,7),u(i,j,k,8), &
                    !u(i,j,k,9),u(i,j,k,10),u(i,j,k,11),u(i,j,k,12), &
                    !u(i,j,k,13),u(i,j,k,14),u(i,j,k,15),u(i,j,k,16), u(i,j,k,17) &
                    !, i=1,Ncellx )
                    end do
            end do
         end do
         close(150)
         do k = 1, Ncellz-1
            do j = 1, Ncelly-1
               do i = 1, Ncellx-1


                  v(i+IST*(mesh-1),j+JST*(mesh-1),k+KST*(mesh-1),1)=u(i,j,k,1)
                  goto 200
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,2)=u(i,j,k,2)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,3)=u(i,j,k,3)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,4)=u(i,j,k,4)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,5)=u(i,j,k,5)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,6)=u(i,j,k,6)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,7)=u(i,j,k,7)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,8)=u(i,j,k,8)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,9)=u(i,j,k,9)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,10)=u(i,j,k,10)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,11)=u(i,j,k,11)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,12)=u(i,j,k,12)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,13)=u(i,j,k,13)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,14)=u(i,j,k,14)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,15)=u(i,j,k,15)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,16)=u(i,j,k,16)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,17)=u(i,j,k,17)
200 continue
              ! u(:,:,:,:) =0.0e0
            end do
         end do
      end do
   end do
   !inquire (IOLENGTH=rec_len) v
   !open(1,file='Dall'//time//'.raw',status = 'unknown', &
   !     form='unformatted', access='direct',recl=rec_len)
   open(100,file=trim(data)//'.DAT')
   do k=1,lengthz
      do j=1,lengthy
         do i=1,lengthx
            !write(1,rec=1)  x(i),y(j),z(k),v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7),&
   write(100,*) v(i,j,k,1)!,v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7),&
                   ! v(i,j,k,8)!,v(i,j,k,9),v(i,j,k,10),v(i,j,k,11),v(i,j,k,12),v(i,j,k,13),v(i,j,k,14),v(i,j,k,15), &
                    !v(i,j,k,16),v(i,j,k,17)

         end do
      end do
      write(*,*) k,data
   end do
   close(100)

end do
!870 continue
!end do
!   open(140,file='cdnt.DAT')
!   do i= 1,lengthx
!      write(140,*) sngl(x(i))
!   end do
!   close(140)


!open(1,file='myarray.raw',status = 'unknown', &
!form='unformatted', access='direct',recl=rec_len)
!write(1,rec=1) array
!close(1)
!end do
!goto 900
DEALLOCATE(inputlength)
DEALLOCATE(x)
DEALLOCATE(y)
DEALLOCATE(z)

DEALLOCATE(v)
DEALLOCATE(u)

!DEALLOCATE(pbr)
!DEALLOCATE(pbl)
!DEALLOCATE(pbrall)
!DEALLOCATE(pblall)
!900 continue
end program computeraw
