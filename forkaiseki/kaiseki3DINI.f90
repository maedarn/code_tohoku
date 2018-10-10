program computeraw
 implicit none
 !integer :: nx, ny, nz
 integer(kind=8) :: rec_len
 !integer :: i,j,k
 !parameter(nx = 100, ny = 20, nz = 50)
 !real, dimension(nx,ny,nz) :: array
 real(8), allocatable :: u(:,:,:,:) !, v(:,:,:,:)
 real(4), allocatable ::  v(:,:,:,:)
 double precision , allocatable :: x(:),y(:),z(:),inputlength(:)
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times,mesh,IST,KST,JST,lengthx,lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore
 character(50) filename
 character(3) NPE,time

  !******parameter********
  core=512
  mesh=32
  val=17
  !lasttime=800
  !initime=800
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
!   read(140,400) ( y(j) , j=1, Ncelly*NSPLTy )
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

   lengthx=(mesh)*NSPLTx
   lengthy=(mesh)*NSPLTy
   lengthz=(mesh)*NSPLTz

   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val))

   !lengthx=mesh*NSPLTx
   !lengthy=mesh*NSPLTy
   !lengthz=mesh*NSPLTz

!   midle=core-core/NSPLTx/NSPLETy/2
!   corexy=core/NSPLTz
!   corexy=NSPLTx*NSPLTy
!   midle=NSPLTz/2
!   corestart=core-corexy*midle


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
         open(unit=150,file='INIT'//NPE//'.dat')
         !open(unit=150,file='2D000031.dat',FORM='UNFORMATTED') !バイナリ
         !open(unit=250,file='NEW'//time//NPE//'.dat') !バイナリ

         !li=0; IF(IST.eq.0) li=1
         !lj=0; IF(JST.eq.0) lj=1
         !lk=0; IF(KST.eq.0) lk=1

         do k = 1, Ncellz !*************************************************** +1 いる？ **********************************
            do j = 1, Ncelly
               do i=1, Ncellx
                read(150,*) u(i,j,k,1)! &!,u(i,j,k,2),u(i,j,k,3),u(i,j,k,4),u(i,j,k,5), &
                    !u(i,j,k,6),u(i,j,k,7),u(i,j,k,8), &
                    !u(i,j,k,9),u(i,j,k,10),u(i,j,k,11),u(i,j,k,12), &
                    !u(i,j,k,13),u(i,j,k,14),u(i,j,k,15),u(i,j,k,16), u(i,j,k,17) &
                    !, i=1,Ncellx )
                    end do
            end do
         end do
         close(150)
         do k = 1, Ncellz
            do j = 1, Ncelly
               do i = 1, Ncellx


                  v(i+IST*(mesh),j+JST*(mesh),k+KST*(mesh),1) = sngl(u(i,j,k,1) )
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
   open(100,file='INIT.DAT')
   do k=1,lengthz
      do j=1,lengthy
         do i=1,lengthx
            !write(1,rec=1)  x(i),y(j),z(k),v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7),&
   write(100,*) v(i,j,k,1)!,v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7),&
                   ! v(i,j,k,8)!,v(i,j,k,9),v(i,j,k,10),v(i,j,k,11),v(i,j,k,12),v(i,j,k,13),v(i,j,k,14),v(i,j,k,15), &
                    !v(i,j,k,16),v(i,j,k,17)

         end do
      end do
      write(*,*) k
   end do
   close(100)
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
end program computeraw
