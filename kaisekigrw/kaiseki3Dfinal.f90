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
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname,iiii
 character(50) filename,data
 character(3) NPE!,time
 character(6) time
 double precision dua,dub,duc
 character(250) pwd
 real(4) rms,rms2,rms3,rms4,dx
 

  !******parameter********
  core=512
  mesh=64+2
  val=7
  lasttime=0
  initime=0
  Ncellx=64
  Ncelly=64
  Ncellz=64
  dx=100.e0/real(64*8)
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

   write(*,*) 'pwd?'
   write(*,*)
   read(*,*) pwd
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
      rms=0.0e0
      rms2=0.0e0
      rms3=0.0e0
      rms4=0.0e0
      do usecore=0,core-1,1
         IST = mod(usecore,NSPLTx)
         KST = usecore/(NSPLTx*NSPLTy)
         JST = usecore/NSPLTx-NSPLTy*KST
         !KST = mod((k-corestart)/4,NSPLTy)
         write(*,*) IST,JST,KST,time,usecore
         write(time,'(I6.6)') times
         write(NPE,'(I3.3)') usecore
         !open(unit=150,file='/Users/maeda/Desktop/kaiseki/'//trim(pwd)//NPE//'.DAT')
         open(unit=150,file='/Users/maeda/Desktop/kaiseki/multigrid-plane-512(512core)-gosa/final000'//NPE//'.DAT')
         !open(unit=160,file='/Users/maeda/Desktop/kaiseki/multigrid-plane-512(512core)-gosa/phiexact'//NPE//'.DAT')
         !open(unit=150,file=trim(pwd)//NPE//'.DAT')!,FORM='UNFORMATTED') !バイナリ
    !     open(unit=150,file=trim(data)//NPE//'.dat')
         !open(unit=150,file='2D000031.dat',FORM='UNFORMATTED') !バイナリ
         !open(unit=250,file='NEW'//time//NPE//'.dat') !バイナリ
         do k = 1, mesh-2 !*************************************************** +1 いる？ **********************************
            do j = 1, mesh-2
               do i= 1, mesh-2
                  ! do k = 1, mesh-1 !*************************************************** +1 いる？ **********************************
                  ! do j = 1, mesh-1
                  !  do i= 1, mesh-1
               read(150,*) u(i,j,k,1),u(i,j,k,2),u(i,j,k,3),u(i,j,k,4),u(i,j,k,5),u(i,j,k,6),u(i,j,k,7)!(u(i,j,k,1) ,i=-1,mesh)
            end do
            end do
         end do
        ! write(*,*)'aa'
         !do k = -1, mesh !*************************************************** +1 いる？ **********************************
         !   do j = -1, mesh
         !      do i= -1, mesh
         !         read(160,*) u(i,j,k,2)
         !         u(i,j,k,3)=u(i,j,k,1)-u(i,j,k,2)
         !   end do
         !   end do
         !end do
         close(150)
         !close(160)
         do k = 1, Ncellz
            do j = 1, Ncelly
               do i = 1, Ncellx
                  !v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),1)=(u(i,j,k,1)+u(i+1,j,k,1)+u(i,j+1,k,1)+u(i,j,k+1,1)&
                  !     +u(i,j+1,k+1,1)+u(i+1,j+1,k,1)+u(i+1,j,k+1,1)+u(i+1,j+1,k+1,1))/8.0d0
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),1)=u(i,j,k,1)
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),2)=u(i,j,k,2)
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),3)=u(i,j,k,3)
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),4)=u(i,j,k,4)
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),5)=u(i,j,k,5)
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),6)=u(i,j,k,6)
                  v(i+IST*(Ncellx),j+JST*(Ncelly),k+KST*(Ncellz),7)=u(i,j,k,7)
               end do
            end do
         end do
          !write(*,*)'aa'
      end do
      do k=1,lengthz
         do j=1,lengthy-2
            do i=1,lengthx
               rms=v(i,j,k,3)**2+rms
            end do
         end do
      end do
      !rms=rms/real(lengthz*lengthy*(lengthx))
      write(*,*) 'rms1'

      do k=1,lengthz
         do j=1,lengthy
            do i=2,lengthx-1
               rms2=((-v(i-1,j,k,1)+v(i+1,j,k,1)-v(i-1,j,k,2)+v(i+1,j,k,2))*0.5e0)**2+rms2
            end do
         end do
      end do
      !rms2=rms2/real(lengthz*lengthy*(lengthx-2))
      write(*,*) 'rms2'

      do k=1,lengthz
         do j=2,lengthy-1
            do i=1,lengthx
               rms3=((-v(i,j-1,k,1)+v(i,j+1,k,1)-v(i,j-1,k,2)+v(i,j+1,k,2))*0.5e0)**2+rms3
            end do
         end do
      end do
      write(*,*) 'rms3'
      !rms3=rms3/real(lengthz*(lengthy-2)*(lengthx))

      do k=2,lengthz-1
         do j=1,lengthy
            do i=1,lengthx
               rms4=((-v(i,j,k-1,1)+v(i,j,k+1,1)-v(i,j,k-1,2)+v(i,j,k+1,2))*0.5e0)**2+rms4
            end do
         end do
      end do
      write(*,*) 'rms4'
      !rms3=rms3/real(lengthz*(lengthy-2)*(lengthx))
      !rms4=rms4/real((lengthz-2)*(lengthy)*(lengthx))
      !v(0,:,:,1)= v(0,:,:,2)
      !v(:,0,:,1)= v(:,0,:,2)
      !v(:,:,0,1)= v(:,:,0,2)
      !v(lengthx-1,:,:,1)= v(lengthx-1,:,:,2)
      !v(:,lengthy-1,:,1)= v(:,lengthy-1,:,2)
      !v(:,:,lengthz-1,1)= v(:,:,lengthz-1,2)
!      v(0,:,:,1)= v(1,:,:,1)
!      v(:,0,:,1)= v(:,1,:,1)
!      v(:,:,0,1)= v(:,:,1,1)
!      v(lengthx-1,:,:,1)= v(lengthx-2,:,:,1)
!      v(:,lengthy-1,:,1)= v(:,lengthy-2,:,1)
!      v(:,:,lengthz-1,1)= v(:,:,lengthz-2,1)

!      dx=1.e0/dx

!      do k=1,lengthz
!         do j=1,lengthy
!            do i=1,lengthx
!               v(i,j,k,4)=((-v(i-1,j,k,2)+v(i+1,j,k,2))*dx*0.5e0)**2/sqrt(((-v(i,j,k-1,2)&
!                    +v(i,j,k+1,2))*dx*0.5e0)**2+((-v(i-1,j,k,2)+v(i+1,j,k,2))*dx*0.5e0)**2&
!                    +((-v(i,j-1,k,2)+v(i,j+1,k,2))*dx*0.5e0)**2)
!            end do
!         end do
!      end do
      open(100,file='/Users/maeda/Desktop/kaiseki/multigrid-plane-512(512core)-gosa/final.DAT')
      do k=1,lengthz
         do j=1,lengthy
            do i=1,lengthx
               write(100,*) v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7)&
                    ,v(i,j,k,6)/v(i,j,k,5)
               !,(-v(i-1,j,k,1)+v(i+1,j,k,1))*dx*0.5e0 &
                 !   ,(-v(i-1,j,k,2)+v(i+1,j,k,2))*dx*0.5e0,(-v(i-1,j,k,1)+v(i+1,j,k,1)+v(i-1,j,k,2)-v(i+1,j,k,2))*dx*0.5e0&
                 !   ,v(i,j,k,4)
            end do
         end do
         write(*,*) k,times
      end do
      close(100)
      rms=rms/real(lengthz*lengthy*(lengthx))
      rms2=rms2/real(lengthz*lengthy*(lengthx-2))
      rms3=rms3/real(lengthz*(lengthy-2)*(lengthx))
      rms4=rms4/real((lengthz-2)*(lengthy)*(lengthx))
      write(*,*) times,rms,rms2,rms3,rms4
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

end program computeraw
