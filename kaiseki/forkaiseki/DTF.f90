program main
  implicit none
  integer i,j,k,kwmax, NSPLTx, NSPLTy, NSPLTz,IS,JS,KS,NPE,MRANK
  double precision , dimension(:,:,:) , allocatable :: U
  double precision , dimension(:) , allocatable :: x,y,z
  !integer(kind=8), allocatable, dimension(:) :: seed
  double precision rho,n,kw,pi,Lbox,dxx,dyy,dzz,rho0,Linv,Klm,dsp
  integer dx,dy,dz,Ncellx,Ncelly,Ncellz,a,b,c,IST,JST,KST
  !integer(kind=8) :: clock
  !integer(kind=8) :: nrand
  double precision , allocatable :: rand(:,:,:)
  character(3) NPENUM

  !******parameter******
  Lbox=20.0d0
  dx=128
  dy=128
  dz=128
  kwmax=16
  rho0= 1.2998844461310d0
  n=(1.29988d-1)/2.0d0
  NSPLTx=4
  NSPLTy=4
  NSPLTz=4
  NPE=64
  !******parameter******
  write(*,*) n
  Linv=1.0d0/Lbox
  Klm=-11.0d0/3.0d0
  Ncellx=dx/NSPLTx
  Ncelly=dy/NSPLTy
  Ncellz=dz/NSPLTz

  ALLOCATE(x(1:dx),y(1:dy),z(1:dz))
  ALLOCATE(U(-1:dx+2,1:dy,1:dz))
  ALLOCATE(rand(-dx/2:dx/2,-dy/2:dy/2,-dz/2:dz/2))
  call random_seed
  call random_number(rand(-dx/2:dx/2,-dy/2:dy/2,-dz/2:dz/2))

  pi=3.1415926535d0!acos(-1.0d0)
  x(:)=0.0d0
  y(:)=0.0d0
  z(:)=0.0d0
  U(:,:,:)=rho0

  dxx=Lbox/dble(dx)
  dyy=Lbox/dble(dy)
  dzz=Lbox/dble(dz)

  x(1)=0.5d0*dxx
  y(1)=0.5d0*dyy
  z(1)=0.5d0*dzz
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


  do k= -kwmax,kwmax,1
     write(*,*) k,'point1'
     do j= -kwmax,kwmax,1
        do i= -kwmax,kwmax,1
           kw=dsqrt(dble(i)**2 + dble(j)**2 + dble(k)**2)
           if(k==0 .and. j==0 .and. i==0) then
              goto 2018
           end if
           if(kw > kwmax) then
              goto 2019
           end if
           call random_number(rand(i,j,k))
           write(*,*) i,j,k,kw,'point2',dble(k)
           do c=1,dz
              do b=1,dy
                 do a=1,dx
                    rho=n*(kw**(Klm))*dsin(2*pi*(dble(i)*x(a)+dble(j)*y(b)+dble(k)*z(c))*Linv+2*pi*rand(i,j,k))
                    U(a,b,c)=rho+U(a,b,c)
                    if(U(a,b,c)<0) then
                       write(*,*) 'ERR'
                       goto 2020
                    end if
                    if(U(a,b,c)>1.0d3) then
                       write(*,*) 'ERR'
                       goto 2021
                    end if
                 end do
              end do
           end do
2018       continue
2019       continue
        end do
     end do
  end do

  U(-1,:,:)=U(dx-1,:,:)
  U(0,:,:)=U(dx,:,:)
  U(dx+1,:,:)=U(1,:,:)
  U(dx+2,:,:)=U(2,:,:)

  dsp=0.0d0
 do c=1,dz
    do b=1,dy
       do a=1,dx
        dsp=dsp+(U(a,b,c)-rho0)**2
     end do
  end do
end do
  dsp=dsp/dble(dy*dz*dx)
  dsp=dsqrt(dsp)
write(*,*) 'dsp=',dsp


  !***************************************
  !do k=1,dz
  !  do j=1,dy
  !    do i=1,dx
  !      U(i,j,k)=dble(i+k+j)
  !   end do
  !end do
  !end do
  !***************************************

  do MRANK = 0, NPE-1
     IS = mod(MRANK,NSPLTx); KS = MRANK/(NSPLTx*NSPLTy); JS = MRANK/NSPLTx-NSPLTy*KS
     write(*,*) IS,KS,JS
     !if((JS.eq.JST).and.(KS.eq.KST)) then
     WRITE(NPENUM,'(I3.3)') MRANK
     open(unit=8,file='D'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     do k = Ncellz*KS+1,  Ncellz*KS+Ncellz
        do j = Ncelly*JS+1,  Ncelly*JS+Ncelly
           !write(*,*) IS,KS,JS
           write(8) (sngl(U(i,j,k)),i=Ncellx*IS+1,Ncellx*IS+Ncellx)
        end do
     end do
     close(8)



     !***************************************
     ! end if
     open(11,FILE='2D'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     k=1;do j = Ncelly*JS+1,  Ncelly*JS+Ncelly
     write(11) (sngl(U(i,j,k)), i=Ncellx*IS+1,Ncellx*IS+Ncellx)
     end do
     close(11)
     !***************************************



   end do

2020 continue
2021 continue
end program main
