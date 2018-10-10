program main
  implicit none
  integer i,j,k,kwmax, NSPLTx, NSPLTy, NSPLTz,IS,JS,KS,NPE,MRANK
  double precision , dimension(:,:,:) , allocatable :: U
  double precision , dimension(:) , allocatable :: x,y,z
  !integer(kind=8), allocatable, dimension(:) :: seed
  double precision rho,n,kw,pi,Lbox,dxx,dyy,dzz,rho0,Linv,Klm,delta,deltax,dsp
  integer dx,dy,dz,Ncellx,Ncelly,Ncellz,a,b,c
  !integer(kind=8) :: clock
  !integer(kind=8) :: nrand
  double precision , allocatable :: rand(:,:,:)
  double precision , allocatable :: plane(:,:)
  character(3) NPENUM

  !******parameter******
  Lbox=20.0d0
  dx=256
  dy=256
  dz=256
  kwmax=50
  rho0= 1.2998844461310257d0
  n=1.2998844461310257d-1/2.0d0
  NSPLTx=4
  NSPLTy=4
  NSPLTz=4
  NPE=64
  delta=Lbox/40.d0
  !******parameter******

  Linv=1.0d0/Lbox
  Klm=-11.0d0/3.0d0
  Ncellx=dx/NSPLTx
  Ncelly=dy/NSPLTy
  Ncellz=dz/NSPLTz

  ALLOCATE(x(1:dx),y(1:dy),z(1:dz))
  ALLOCATE(U(1:dx,1:dy,1:dz))
  ALLOCATE(plane(-1:dy+2,-1:dz+2))
  ALLOCATE(rand(-dx/2:dx/2,-dy/2:dy/2,-dz/2:dz/2))
  call random_seed
  call random_number(rand(-dx/2:dx/2,-dy/2:dy/2,-dz/2:dz/2))

  pi=acos(-1.0d0)
  x(:)=0.0d0
  y(:)=0.0d0
  z(:)=0.0d0
  U(:,:,:)=0.0d0
  plane(:,:)=Lbox/2.0d0

  dxx=Lbox/dble(dx)
  dyy=Lbox/dble(dy)
  dzz=Lbox/dble(dz)

  do i=2,dx
     x(i)=dxx+x(i-1)
  end do
  do i=2,dx
     y(i)=dyy+y(i-1)
  end do
  do i=2,dx
     z(i)=dzz+z(i-1)
  end do

  !call random_seed(size=nrand)
  !allocate(seed(nrand))
  !call system_clock(count=clock)
  !seed = clock
  !call random_seed(put=seed)

  dsp=0.0d0

  do k= -kwmax,kwmax,1
     write(*,*) k,'point1'
     do j= -kwmax,kwmax,1
        kw=dsqrt(dble(j**2+k**2))
        if(k==0 .and. j==0) then
           goto 2018
        end if
        if(kw+1 > dble(kwmax)) then
           goto 2019
        end if
        !if(kw+1 < 20) then
        !   goto 2020
        !end if
        !do i= -kwmax,-1,1
           call random_number(rand(i,j,k))
           !write(*,*) rand(i,j,k),j,k
           do c=1,dz
              do b=1,dy
                 !do a=1,dx
                 deltax=delta*(kw**(Klm))*dsin(2*pi*(j*y(b)+k*z(c))*Linv+2*pi*rand(i,j,k))
                 plane(b,c)=deltax+plane(b,c)
                ! end do
              end do
           end do
2018       continue
2019       continue
!2020       continue
        !end do
     end do
  end do
  plane(-1,:)=plane(dy-1,:)
  plane(0,:)=plane(dy,:)
  plane(dy+1,:)= plane(1,:)
  plane(dy+2,:)= plane(2,:)

  plane(:,-1)=plane(:,dz-1)
  plane(:,0)=plane(:,dz)
  plane(:,dz+1)=plane(0,1)
  plane(:,dz+2)=plane(:,2)

  do c=1,dz
     do b=1,dy
        dsp=dsp+(plane(b,c)-Lbox/2.0d0)**2
     end do
  end do
  dsp=dsp/dble(dy*dz)
  dsp=dsqrt(dsp)
write(*,*) 'dsp=',dsp

!  do k= 1,kwmax,1
!     write(*,*) k,'point2'
!     do j= 1,kwmax,1
!       ! do i= 1,kwmax,1
!           call random_number(rand(i,j,k))
!           do c=1,dz
!              do b=1,dy
!                 !do a=1,dx
!                    kw=dsqrt(dble(j**2+k**2))
!                    deltax=delta*kw**(Klm)*dsin(2*pi*(j*y(b)+k*z(c))*Linv+2*pi*rand(i,j,k))
!                    plane(b,c)=deltax+plane(b,c)
!                 !end do
!              end do
!           end do
!        !end do
!     end do
!  end do

open(unit=18,file='dxplane3.dat')
do j=1,dz
   do i=1,dy
      write(18,*) y(i),z(j),plane(i,j)
   end do
   write(18,*)
end do
close(18)

open(unit=28,file='delta3.dat',FORM='UNFORMATTED')
do j=-1,dz+2
   do i=-1,dy+2
      write(28) plane(i,j)
   end do
   !write(18,*)
end do
close(28)

  !***************************************
  !do k=1,dz
  !  do j=1,dy
  !    do i=1,dx
  !      U(i,j,k)=dble(i+k+j)
  !   end do
  !end do
  !end do
  !***************************************

!  do MRANK = 0, NPE-1
!     IS = mod(MRANK,NSPLTx); KS = MRANK/(NSPLTx*NSPLTy); JS = MRANK/NSPLTx-NSPLTy*KS
!     write(*,*) IS,KS,JS
!    !if((JS.eq.JST).and.(KS.eq.KST)) then
!      WRITE(NPENUM,'(I3.3)') MRANK
!      open(unit=8,file='Dplane'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
!      do k = Ncellz*KS+1,  Ncellz*KS+Ncellz
!         do j = Ncelly*JS+1,  Ncelly*JS+Ncelly
!            !write(*,*) IS,KS,JS
!           write(8) (U(i,j,k),i=Ncellx*IS+1,Ncellx*IS+Ncellx)
!         end do
!      end do
!      close(8)



      !***************************************
      ! end if
!      open(11,FILE='2D'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
!      k=1;do j = Ncelly*JS+1,  Ncelly*JS+Ncelly
!      write(11) (sngl(U(i,j,k)), i=Ncellx*IS+1,Ncellx*IS+Ncellx)
!      end do
!      close(11)
      !***************************************



 !  end do


end program main

