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
 real(4) mass
 real(8) mass1
 integer :: val,core=0,NSPLTx,NSPLTy,NSPLTz,times,mesh,IST,KST,JST,lengthx,lengthy,lengthz,corexy,midle,corestart,lasttime,length
 integer :: i,j,k,kz=1,initime, Ncellx, Ncelly, Ncellz,usecore,dataname
 character(50) filename,data
 character(3) NPE,time

  !******parameter********
  core=512
  mesh=64+2
  val=17
  !lasttime=800
  !initime=800
  Ncellx=64+1
  Ncelly=64+1
  Ncellz=64+1
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
   ALLOCATE(u(1:mesh,1:mesh,1:mesh,1:val))
!   ALLOCATE(pbr(0:mesh-1+2,0:mesh-1+2))
!   ALLOCATE(pbl(0:mesh-1+2,0:mesh-1+2))
!   ALLOCATE(pbrall(1:lengthy,1:lengthz))
!   ALLOCATE(pblall(1:lengthy,1:lengthz))

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

!goto 870

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
open(unit=150,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/135'//NPE//'.dat',FORM='UNFORMATTED')
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
!         if((JST<1)) then
           ! (JST>7)) then
!            u(:,:,:,1)=0.d0
!         end if
!         if((JST>5)) then
           ! (JST>7)) then
!            u(:,:,:,1)=0.d0
!         end if
!         if((KST<4)) then
           ! (JST>7)) then
!            u(:,:,:,1)=0.d0
!         end if
         close(150)
         write(*,*) u(5,5,5,1),u(5,5,5,2),u(5,5,5,17),u(64,64,64,17)
         do k = 1, Ncellz-1
            do j = 1, Ncelly-1
               do i = 1, Ncellx-1

    !              if(u(i,j,k,1)>1.d3 ) then

                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)=u(i,j,k,1)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),2)=u(i,j,k,2)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),3)=u(i,j,k,3)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),4)=u(i,j,k,4)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),5)=u(i,j,k,5)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),6)=u(i,j,k,6)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),7)=u(i,j,k,7)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),8)=u(i,j,k,8)
                  v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),17)=u(i,j,k,17)
    !              else
    !                 v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)=0.d0
    !              end if
                  mass=mass+v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1)
                  mass1=mass1+dble(v(i+IST*(mesh-2),j+JST*(mesh-2),k+KST*(mesh-2),1))

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

   mass=mass*100.e0*100.e0*100.e0/real(512)/real(512)/real(512)
   mass1=mass1*100.d0*100.d0*100.d0/dble(512)/dble(512)/dble(512)
   write(*,*) mass,mass1
   !inquire (IOLENGTH=rec_len) v
   !open(1,file='Dall'//time//'.raw',status = 'unknown', &
   !     form='unformatted', access='direct',recl=rec_len)
   open(100,file='/Users/maeda/Desktop/kaiseki/cnv100wbwg/all.DAT',FORM='FORMATTED')
   do k=1,lengthz
      do j=1,lengthy
         do i=1,lengthx
            !write(1,rec=1)  x(i),y(j),z(k),v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7),&
   write(100,*) v(i,j,k,1),v(i,j,k,2),v(i,j,k,3),v(i,j,k,4),v(i,j,k,5),v(i,j,k,6),v(i,j,k,7),&
                    v(i,j,k,8)&!,v(i,j,k,9),v(i,j,k,10),v(i,j,k,11),v(i,j,k,12),v(i,j,k,13),v(i,j,k,14),v(i,j,k,15), &
                    !v(i,j,k,16)
                    ,v(i,j,k,17)
   !write(*,*)v(i,j,k,1)
end do
      end do
      write(*,*) k,data
   end do
   close(100)
    write(*,*) mass,mass1

!end do

!DEALLOCATE(inputlength)
!DEALLOCATE(x)
!DEALLOCATE(y)
!DEALLOCATE(z)

DEALLOCATE(v)
DEALLOCATE(u)

end program computeraw
