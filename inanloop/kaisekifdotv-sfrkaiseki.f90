program main
  implicit none
  integer i,j,k,val,m,loop,n,n1,n2,n3,dm,time1,timeloop,jump,timestep,cnt10
  integer nx,ny,nz,rgn,zero,sloop,smesh,err
  integer,allocatable,dimension(:,:,:,:) ::tag
  real(4),allocatable,dimension(:,:,:,:) ::U
  character(3) num,timech,ch1
  character(1) num1
  double precision Lbox,Msun,dl,Mshd,ganmim1,Mtot,mu
  integer,allocatable,dimension(:) :: cas,count1,count1S
  integer,allocatable,dimension(:) :: ig,jg,kg
  real(4) dm1
  double precision ,allocatable,dimension(:) :: MM,S,lgM,phi,tff1,vdsp,vdspx,vdspy,vdspz,csmean,csmean1
  double precision ,allocatable,dimension(:) :: x,y,z,xg,yg,zg,vx,vy,vz,vxg,vyg,vzg
  double precision ,allocatable,dimension(:) :: emag,ekin,ethm,egrv,div,Tmean
  double precision ,allocatable,dimension(:,:,:) :: jx,jy,jz,fx,fy,fz,fgx,fgy,fgz,fpx,fpy,fpz,fbx,fby,fbz,divv,Tn
  double precision ,allocatable,dimension(:,:,:) :: vdotf,vdotfg,vdotfp,vdotfb
  integer,allocatable,dimension(:):: lgMcnt
  double precision cs,nsisu,ff,vv,ffb,ffg,ffp!,tt1
  double precision :: kb=8.63359d0
  DOUBLE PRECISION, parameter :: G=1.11142d-4
  double precision ,allocatable,dimension(:,:,:) :: UU
  double precision ,allocatable,dimension(:,:,:) :: UU2
  double precision ,allocatable,dimension(:,:,:,:) :: UU3
  double precision ,allocatable,dimension(:) :: tt1
  double precision a1,a2,a3,a4

  !----parameter----
  nx=512
  ny=512
  nz=512
  val=9
  loop=5
  Lbox=100.d0
  Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
  !rgn=3000
  rgn=30
  Mshd=1.d2
  ganmim1=3.d0/2.d0 !1.d0/((5.d0/3.d0-1.d0)
  sloop=30
  smesh=5
  timeloop=50
  timestep=4
  jump=4
  mu=1.27d0
  !----parameter----

  dl=Lbox/dble(nx)


  allocate(UU(loop,rgn,21),UU2(loop,rgn,8),UU3(loop,rgn,8,50),tt1(50))
  UU(:,:,:)=0.d0
  UU2(:,:,:)=0.d0
  UU3(:,:,:,:)=0.d0
  tt1(:)=0.d0

  do time1=1,timeloop
     open(unit=350,file='cnt3.dat')
     read(350,*) cnt10
     close(350)
     !write(cntc,'(i3.3)') cnt10
     write(timech,'(i3.3)') cnt10
     cnt10=cnt10+timestep
     open(unit=350,file='cnt3.dat')
     write(350,*) cnt10
     close(350)
     !write(timech,'(i3.3)') time1
     tt1(time1)=0.25d0*dble(cnt10)
     do i=1,loop
        write(ch1,'(i3.3)') i
        open(unit=200,file='mdiv'//timech//ch1//'.DAT')
        do j=1,rgn
           read(200,*,iostat=err) UU(i,j,1),UU(i,j,2),UU(i,j,3),UU(i,j,4),UU(i,j,5),UU(i,j,6),UU(i,j,7),UU(i,j,8),&
                UU(i,j,9),UU(i,j,10),UU(i,j,11),UU(i,j,12),UU(i,j,13),UU(i,j,14),UU(i,j,15),UU(i,j,16),UU(i,j,17),UU(i,j,18),&
                UU(i,j,19),UU(i,j,20),UU(i,j,21)

           !write(*,*)'ok',i,j
           !read(*, *, iostat=err) i
           !if(err == 0) then
              ! エラーがないときの処理
           !   print *, "  success. ", i
           !else if(err > 0) then
              ! エラーが起こったときの処理
           !   print *, "  error."
           !   exit
           !end if
           if(err>0)then
!520           UU(i,j,:)=0.d0
               UU(i,j,:)=0.d0
            end if
        end do
        close(200)
     end do

     !(0.05*$12*$13/$18)
     do i=1,loop
        do j=1,rgn
           if((UU(i,j,18).ne.0.d0).and.(UU(i,j,20).ne.0.d0).and.(UU(i,j,12).ne.0.d0)) then
              UU2(i,j,1)=UU(i,j,2)
              UU2(i,j,2)=UU(i,j,13)
              UU2(i,j,3)=0.05d0*UU(i,j,12)*UU(i,j,13)/UU(i,j,18)
              UU2(i,j,4)=0.01d0*UU(i,j,12)*UU(i,j,13)/UU(i,j,18)
              UU2(i,j,5)=0.05d0*UU(i,j,12)*UU(i,j,13)/UU(i,j,20)
              UU2(i,j,6)=0.01d0*UU(i,j,12)*UU(i,j,13)/UU(i,j,20)
              UU2(i,j,7)=UU(i,j,2)/(0.05d0*UU(i,j,12))
              UU2(i,j,8)=UU(i,j,1)
           end if
        end do
     end do
     UU3(:,:,:,time1)=UU2(:,:,:)
  end do


  open(unit=360,file='mass1.dat')
  open(unit=370,file='mass2.dat')
  open(unit=380,file='mass3.dat')
  open(unit=390,file='mass4.dat')
!  open(unit=760,file='mass1-1.dat')
!  open(unit=770,file='mass2-1.dat')
!  open(unit=780,file='mass3-1.dat')
!  open(unit=790,file='mass4-1.dat')
  open(unit=400,file='mass.dat')

  do time1=1,timeloop
     do i=1,loop
        do j=1,rgn
           if((UU3(i,j,3,time1).ne.0.d0).and.(UU3(i,j,4,time1).ne.0.d0).and.(UU3(i,j,5,time1).ne.0.d0)&
                .and.(UU3(i,j,6,time1).ne.0.d0)) then
              a1=UU3(i,j,1,time1)/UU3(i,j,3,time1)
              a2=UU3(i,j,1,time1)/UU3(i,j,4,time1)
              a3=UU3(i,j,1,time1)/UU3(i,j,5,time1)
              a4=UU3(i,j,1,time1)/UU3(i,j,6,time1)

              if(a1<1.d0)then
                 write(360,*) tt1(time1),UU3(i,j,1,time1),UU3(i,j,2,time1),UU3(i,j,3,time1),UU3(i,j,4,time1)&
                      ,UU3(i,j,5,time1),UU3(i,j,6,time1),UU3(i,j,7,time1),UU3(i,j,8,time1),i,j,a1,a2,a3,a4
!                 write(760,*) tt1(time1-1),UU3(i,j,1,time1-1),UU3(i,j,2,time1-1),UU3(i,j,3,time1-1),UU3(i,j,4,time1-1)&
!                      ,UU3(i,j,5,time1-1),UU3(i,j,6,time1-1),UU3(i,j,7,time1-1),UU3(i,j,8,time1-1),i,j!,a1,a2,a3,a4
              end if
              if(a2<1.d0)then
                 write(370,*)tt1(time1),UU3(i,j,1,time1),UU3(i,j,2,time1),UU3(i,j,3,time1),UU3(i,j,4,time1)&
                      ,UU3(i,j,5,time1),UU3(i,j,6,time1),UU3(i,j,7,time1),UU3(i,j,8,time1),i,j,a1,a2,a3,a4
!                 write(770,*) tt1(time1-1),UU3(i,j,1,time1-1),UU3(i,j,2,time1-1),UU3(i,j,3,time1-1),UU3(i,j,4,time1-1)&
!                      ,UU3(i,j,5,time1-1),UU3(i,j,6,time1-1),UU3(i,j,7,time1-1),UU3(i,j,8,time1-1),i,j!,a1,a2,a3,a4
              end if
              if(a3<1.d0)then
                 write(380,*)tt1(time1),UU3(i,j,1,time1),UU3(i,j,2,time1),UU3(i,j,3,time1),UU3(i,j,4,time1)&
                      ,UU3(i,j,5,time1),UU3(i,j,6,time1),UU3(i,j,7,time1),UU3(i,j,8,time1),i,j,a1,a2,a3,a4
!                 write(780,*) tt1(time1-1),UU3(i,j,1,time1-1),UU3(i,j,2,time1-1),UU3(i,j,3,time1-1),UU3(i,j,4,time1-1)&
!                      ,UU3(i,j,5,time1-1),UU3(i,j,6,time1-1),UU3(i,j,7,time1-1),UU3(i,j,8,time1-1),i,j!,a1,a2,a3,a4
              end if
              if(a4<1.d0)then
                 write(390,*)tt1(time1),UU3(i,j,1,time1),UU3(i,j,2,time1),UU3(i,j,3,time1),UU3(i,j,4,time1)&
                      ,UU3(i,j,5,time1),UU3(i,j,6,time1),UU3(i,j,7,time1),UU3(i,j,8,time1),i,j,a1,a2,a3,a4
!                 write(790,*) tt1(time1-1),UU3(i,j,1,time1-1),UU3(i,j,2,time1-1),UU3(i,j,3,time1-1),UU3(i,j,4,time1-1)&
!                      ,UU3(i,j,5,time1-1),UU3(i,j,6,time1-1),UU3(i,j,7,time1-1),UU3(i,j,8,time1-1),i,j!,a1,a2,a3,a4
              end if
              write(400,*)tt1(time1),UU3(i,j,1,time1),UU3(i,j,2,time1),UU3(i,j,3,time1),UU3(i,j,4,time1)&
                      ,UU3(i,j,5,time1),UU3(i,j,6,time1),UU3(i,j,7,time1),UU3(i,j,8,time1),i,j,a1,a2,a3,a4
           end if
        end do
        write(360,*)
        write(370,*)
        write(380,*)
        write(390,*)
        write(400,*)
!        write(760,*)
!        write(770,*)
!        write(780,*)
!        write(790,*)
     end do
  end do
  close(360)
  close(370)
  close(380)
  close(390)
!  close(760)
!  close(770)
!  close(780)
!  close(790)
  close(400)


  deallocate(UU,UU2,UU3,tt1)



end program main
