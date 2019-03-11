program computeraw
 implicit none
 !integer :: nx, ny, nz
 integer(kind=8) :: rec_len
 integer cnt10,jump,ts,timejump,icn,cc,timestep
 !integer :: i,j,k
 !parameter(nx = 100, ny = 20, nz = 50)
 !real, dimension(nx,ny,nz) :: array
 real(4), allocatable :: u(:,:,:,:) !, v(:,:,:,:)
 real(4), allocatable ::  v(:,:,:,:)
 real(8), allocatable ::  pbr(:,:),pbl(:,:)
 real(8), allocatable ::  pbrall(:,:),pblall(:,:)
 double precision , allocatable :: x(:),y(:),z(:),inputlength(:)
 !real(4) mass
 !real(8) mass1
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times,mesh,IST,KST,JST,lengthx,lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname
 character(50) filename,data
 character(3) NPE,time,cntc
 double precision ,allocatable,dimension(:) :: tff,Mass,ttt
 double precision Lbox,Mass1,tff1,dl,Msun
 real(4) shd,mxrho1
 real(4),allocatable,dimension(:) :: mxrho
 DOUBLE PRECISION, parameter :: G=1.11142d-4

  !******parameter********
  core=512
  mesh=64+2
  val=17
  jump=1
  ts=1
  !lasttime=800
  !initime=800
  Ncellx=64+1
  Ncelly=64+1
  Ncellz=64+1
  timejump=15
  Lbox=100.d0
  timestep=4
  Msun=1.473d-2
  !***********************

  !mass=0.e0
  !mass1=0.d0
  !dl=Lbox/dble(nx)
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
   dl=Lbox/dble(lengthx)

   ALLOCATE(v(1:lengthx,1:lengthy,1:lengthz,1:val))
   ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val))
   allocate(tff(175),Mass(175),ttt(0:175),mxrho(175))


   open(unit=350,file='cnt3.dat')
   read(350,*) cnt10
   close(350)
   write(cntc,'(i3.3)') cnt10
   if(cnt10==1)then
      tff(:)=0.d0
      Mass(:)=0.d0
      ttt(:)=0.d0
      ttt(0)=-1.d0
      mxrho(:)=0.e0
   end if
   if(cnt10.ne.1)then
      open(920,file='tff.DAT',FORM='FORMATTED')
      do i=1,175
         read(920,*) tff(i),Mass(i),ttt(i),mxrho(i)
      end do
      close(920)
   end if
   icn=(cnt10-1)/timestep+1
   !write(cntc,'(i3.3)') cnt10
   !write(cntc,'(i3.3)') cnt10
   cnt10=cnt10+timestep
   open(unit=350,file='cnt3.dat')
   write(350,*) cnt10
   close(350)



   !do times=initime,lasttime
      do usecore=0,core-1,1
         IST = mod(usecore,NSPLTx)
         KST = usecore/(NSPLTx*NSPLTy)
         JST = usecore/NSPLTx-NSPLTy*KST
         write(*,*) IST,JST,KST,time,usecore
   !      write(time,'(I3.3)') times
         write(NPE,'(I3.3)') usecore
   !      open(unit=150,file=time//NPE//'.dat',FORM='UNFORMATTED') !バイナリ
open(unit=150,file=cntc//NPE//'.dat',FORM='UNFORMATTED')
         !open(unit=150,file='2D000031.dat',FORM='UNFORMATTED') !バイナリ
!write(*,*)'kakuknin'
!open(unit=250,file='NEW'//time//NPE//'.dat') !バイナリ
         do k = 1, Ncellz !*************************************************** +1 いる？ **********************************
            do j = 1, Ncelly
               !do i= 1, Ncellx
                  read(150) (u(i,j,k,1) ,u(i,j,k,2),u(i,j,k,3),u(i,j,k,4),u(i,j,k,5), &
                    u(i,j,k,6),u(i,j,k,7),u(i,j,k,8), &
                    u(i,j,k,9),u(i,j,k,10),u(i,j,k,11),u(i,j,k,12), &
                    u(i,j,k,13),u(i,j,k,14),u(i,j,k,15),u(i,j,k,16), u(i,j,k,17) &
                    , i=1,Ncellx)
                    !end do
            end do
         end do
         close(150)
         do k = 1, Ncellz-1
            do j = 1, Ncelly-1
               do i = 1, Ncellx-1

    !              if(u(i,j,k,1)>1.d3 ) then

                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)=u(i,j,k,1)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),2)=u(i,j,k,2)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),3)=u(i,j,k,3)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),4)=u(i,j,k,4)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),5)=u(i,j,k,5)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),6)=u(i,j,k,6)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),7)=u(i,j,k,7)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),8)=u(i,j,k,8)
             !     v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),17)=u(i,j,k,17)
    !              else
    !                 v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)=0.d0
    !              end if
!                  mass=mass+v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)
!                  mass1=mass1+dble(v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1))

                  goto 200
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,2)=u(i,j,k,2)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,3)=u(i,j,k,3)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,4)=u(i,j,k,4)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,5)=u(i,j,k,5)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,6)=u(i,j,k,6)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,7)=u(i,j,k,7)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,8)=u(i,j,k,8)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,9)=u(i,j,k,9)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,10)=u(i,j,k,10)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,11)=u(i,j,k,11)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,12)=u(i,j,k,12)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,13)=u(i,j,k,13)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,14)=u(i,j,k,14)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,15)=u(i,j,k,15)
            !   v(i+IST*mesh,j+JST*mesh,k+KST*mesh,16)=u(i,j,k,16)
               v(i+IST*mesh,j+JST*mesh,k+KST*mesh,17)=u(i,j,k,17)
200 continue
              ! u(:,:,:,:) =0.0e0
            end do
         end do
      end do
   end do

   shd=1.e3
   cc=0
   tff1=0.d0
   Mass1=0.d0
   mxrho1=0.e0

   do k=1,lengthz
      do j=1,lengthy
         do i=1,lengthx
            if(v(i,j,k,1)>shd)then
               !tff1=1/dsqrt(G*dble(U(i,j,k,1)))+tff1
               tff1=tff1+dble(v(i,j,k,1))*dsqrt(G*dble(v(i,j,k,1)))*dl*dl*dl
               cc=cc+1
               Mass1=dble(v(i,j,k,1))*dl*dl*dl+Mass1
               mxrho1=amax1(mxrho1,v(i,j,k,1))
            end if
         end do
      end do
      write(*,*) k,'caltff'
   end do
   !if(cc>0) then
      !tff1=tff1/dble(cc)
      tff1=tff1*Msun
      Mass1=Mass1*Msun
      tff(icn)=tff1
      Mass(icn)=Mass1
      mxrho(icn)=mxrho1
   !end if
   !end do

   ttt(icn)=ttt(icn-1)+1.d0
   open(910,file='tff.DAT',FORM='FORMATTED')
   do i=1,175
      write(910,*) tff(i),Mass(i),ttt(i),mxrho(i)
   end do
   close(910)

DEALLOCATE(v)
DEALLOCATE(u)
deallocate(tff,Mass,ttt)

end program computeraw
