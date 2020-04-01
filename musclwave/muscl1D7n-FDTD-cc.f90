module comvar
  implicit none
  integer , parameter :: ndx=66,ndy=66,itrstp=20000,svnum=1000,dim=2
  double precision :: cg=1.d0 , Lbox=1.d2 , dx , dy , h=10.0d0 , hcen=50.0d0
  double precision , allocatable , dimension(:,:) :: Phi  , rho , Phidt , Phidt2
  double precision , allocatable , dimension(:,:,:) :: vg , vgdt
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0
  DOUBLE PRECISION , allocatable , dimension(:) :: x , y
  character(56) :: dir='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testwave2D/'
end module comvar


program main
  use comvar
  implicit none
  integer :: i,j,k,l,m,n,mode,sv=0
  double precision :: dt

  allocate(Phi(-1:ndx,-1:ndy),rho(-1:ndx,-1:ndy))
  allocate(Phidt(-1:ndx,-1:ndy),Phidt2(-1:ndx,-1:ndy))
  allocate(vg(-1:ndx,-1:ndy,dim))
  allocate(vgdt(-1:ndx,-1:ndy,dim))
  allocate(x(-1:ndx),y(-1:ndy))
  dx=Lbox/dble(ndx)
  dy=Lbox/dble(ndy)

  call INIT()
  dt = cg/dx*coeff

  call source(dt*0.5d0)
  Phidt(:,:) = Phi(:,:)
  call stepphi(dt*0.5d0)
  call BC(1)

  do i = 1,itrstp
     Phidt2(:,:)=Phidt(:,:)
     Phidt(:,:) =Phi(:,:)
     vgdt(:,:,:)  =vg(:,:,:)

     call stepvg(dt)
     call BC(3)

     call source(dt)
     Phidt2(:,:)=Phidt(:,:)
     Phidt(:,:) = Phi(:,:)
     call stepphi(dt)
     call BC(1)

     if(mod(i,svnum)==1) then
        call saveu(sv)
     end if
  end do
  call saveu(sv)

  deallocate(Phi,rho)
  deallocate(Phidt,Phidt2)
  deallocate(vg)
  deallocate(vgdt)
  deallocate(x,y)

end program main

subroutine stepvg(dt)
  use comvar
  integer :: i,j,k,l,m,n
  double precision :: dt
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: slopx,slopy,ur,ul
  !vgdt(:,:,:) = vg(:,:,:)

  !call vanalbadaphi(slopx,slopy)
!  do i = 0,ndx-1
!     ul(i,j) = Phi(i+1,j) + 0.25d0 * ep * slop(i+1) &
!          * ((1.0d0-slop(i+1,j)*kappa)*(Phi(i+1,j)-Phi(i+1-1,j)) + (1.0d0+slop(i+1)*kappa)*(Phi(i+1+1,j) - Phi(i+1,j))) !i+1/2
!     ur(i,j) = Phi(i,j) - 0.25d0 * ep * slop(i) &
!          * ((1.0d0+slop(i)*kappa)*(Phi(i,j)-(i-1,j)) + (1.0d0-slop(i,j)*kappa)*(Phi(i+1,j) - Phi(i,j))) !i-1/2
!  end do


  
  do j = 0 , ndy-2
     do i = 0 , ndx-2
        !vg(i,j,1) = vgdt(i,j,1) -cg*dt/dx*(Phi(i+1,j) - Phi(i,j))*slopx(i+1,j)
        !vg(i,j,2) = vgdt(i,j,2) -cg*dt/dy*(Phi(i,j+1) - Phi(i,j))*slopy(i,j+1)
        !vg(i,j,1) = vgdt(i,j,1) -cg*dt/dx*(Phi(i-1,j)+(Phi(i+1,j) - Phi(i-1,j))*0.5d0)!*slopx(i+1,j)
        !vg(i,j,2) = vgdt(i,j,2) -cg*dt/dy*(Phi(i,j-1)+(Phi(i,j+1) - Phi(i,j-1))*0.5d0)!*slopy(i,j+1)
     end do
  end do
end subroutine stepvg

subroutine stepphi(dt)
  use comvar
  integer :: i,j,k,l,m,n
  double precision :: dt , r , ep=1.0d-20 , pr
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: slopx,slopy
  !Phidt(:,:)=Phi(:,:)
  !call vanalbadavg(slopx,1)
  !call vanalbadavg(slopy,2)

  do j = 1 , ndy-2
     do i = 1 , ndx-2
        !Phi(i,j) = Phidt(i,j) - cg*dt/dx*(vg(i,j,1) - vg(i-1,j,1))*slopx(i-1,j) &
        !     - cg*dt/dy*(vg(i,j,2) - vg(i,j-1,2))*slopy(i,j-1)
        !Phi(i,j) = Phidt(i,j) - cg*dt/dx*(vg(i-1,j,1)+(vg(i,j+1,1) - vg(i-1,j,1))*0.5d0) &
        !     - cg*dt/dy*(vg(i,j-1,2)+(vg(i,j+1,2) - vg(i,j-1,2))*0.5d0)!*slopy(i,j-1)
        !r=(vg(i,j,1) - vg(i-1,j,1)+ep)/(vg(i+1,j,1) - vg(i,j,1)+ep)
        !pr=(r**2+r)/(r**2+1.d0)
        Phi(i,j) = Phidt(i,j) - cg*dt/dx*(vg(i,j,1) - vg(i-1,j,1))*pr &
             - cg*dt/dy*(vg(i,j,2) - vg(i,j-1,2))*slopy(i,j-1)
     end do
  end do
end subroutine stepphi

subroutine source(dt)
  use comvar
  integer :: i,j,k,l,m,n
  double precision :: dt
  do j = 1 , ndy
     do i = 1 , ndx
        Phi(i,j) = 2.d0 * Phidt(i,j) - Phidt2(i,j) - G4pi * cg * cg * dt * dt * rho(i,j)
     end do
  end do
end subroutine source

subroutine INIT()
  use comvar
  integer :: i,j,k,l,m,n
  double precision :: dinit1,meanrho
  double precision :: dt

  !----------x--------------
  dx = Lbox/dble(ndx-2)
  x(1) = dx/2.0d0
  x(0) = x(1) - dx
  x(-1) = x(0) - dx
  do i=2,ndx
     x(i) = x(i-1) + dx
  end do
  !----------x--------------

  !----------y--------------
  dy = Lbox/dble(ndy-2)
  y(1) = dy/2.0d0
  y(0) = y(1) - dy
  y(-1) = y(0) - dy
  do i=2,ndy
     y(i) = y(i-1) + dy
  end do
  !----------y--------------

  !------rho------
  dinit1 = 2.0d0/G4pi/90.d0
  do j = -1,ndy
     do i = -1,ndx
        if( dabs(x(i) - hcen) .le. h) then
           rho(i,j) = dinit1
           !rho(i) = 0.0d0
        else
           rho(i,j) = 0.0d0
           !rho(i) = dinit1
           !rho(i) = dinit1*1.d-2
        end if
     end do
  end do

  meanrho=0.d0
  do j = 1,ndy-2
     do i = 1,ndx-2
        meanrho=meanrho+rho(i,j)
     end do
  end do
  meanrho=meanrho/dble(ndx-2)/dble(ndy-2)

  do j = -1,ndy
     do i = -1,ndx
        rho(i,j)=rho(i,j)-meanrho
     end do
  end do
  !------rho------


  Phidt2(:,:) = 0.d0
  Phidt(:,:) = 0.d0
  Phi(:,:) = 0.d0
  vgdt(:,:,:) = 0.d0
  vg(:,:,:) = 0.d0

end subroutine INIT

subroutine BC(mode)
  use comvar
  integer :: mode
  integer :: i,j,k,l,m,n

  if(mode==1) then

     !--------Phi*x--------
     Phi(-1,:)=Phi(ndx-3,:)
     Phi( 0,:)=Phi(ndx-2,:)

     Phi(ndx-1,:)=Phi(1,:)
     Phi(ndx  ,:)=Phi(2,:)
     !--------Phi*x--------

     !--------Phi*y--------
     Phi(:,-1)=Phi(:,ndx-3)
     Phi(:, 0)=Phi(:,ndx-2)

     Phi(:,ndx-1)=Phi(:,1)
     Phi(:,ndx  )=Phi(:,2)
     !--------Phi*y--------

  end if

  if(mode==3) then

     !--------vg*x--------
     vg(-1,:,1)=vg(ndx-3,:,1)
     vg( 0,:,1)=vg(ndx-2,:,1)

     vg(ndx-1,:,1)=vg(1,:,1)
     vg(ndx  ,:,1)=vg(2,:,1)
     !--------vg*x--------

     !--------vg*y--------
     vg(:,-1,2)=vg(:,ndx-3,2)
     vg(:, 0,2)=vg(:,ndx-2,2)

     vg(:,ndx-1,2)=vg(:,1,2)
     vg(:,ndx  ,2)=vg(:,2,2)
     !--------vg*y--------

  end if
end subroutine BC

subroutine saveu(in1)
  use comvar
  integer :: i,in1
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file=dir//'phi'//name//'.dat')
  do j=-1,ndy
     do i=-1,ndx
        write(21,*) x(i),y(j),Phi(i,j),vg(i,j,1),vg(i,j,2),rho(i,j)
     end do
     write(21,*)
  end do
  close(21)

  open(22,file=dir//'phix'//name//'.dat')
  !do j=-1,ndy
  j=ndy/2.d0
  do i=-1,ndx
     write(22,*) x(i),y(j),Phi(i,j),vg(i,j,1),vg(i,j,2),rho(i,j)
  end do
  write(22,*)
  !end do
  close(22)
  in1=in1+1
end subroutine saveu

subroutine vanalbadaphi(Phigrdx,Phigrdy)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  integer :: i , j , ip , im , xi , yj , jp , jm
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phigrdx,Phigrdy

  do j = 0 , ndy-1
     do i = 0 , ndx-1
        ip=i+1
        im=i-1

        delp = Phi(ip,j)-Phi(i,j)
        delm = Phi(i,j)-Phi(im,j)
        flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
        Phigrdx(i,j) = flmt
     end do
  end do

  do j = 0 , ndy-1
     do i = 0 , ndx-1
        jp=j+1
        jm=j-1

        delp = Phi(i,jp)-Phi(i,j)
        delm = Phi(i,j)-Phi(i,jm)
        flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
        Phigrdy(i,j) = flmt
     end do
  end do
end subroutine vanalbadaphi

subroutine vanalbadavg(Phigrd,mode)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  integer :: i , j , ip , im , xi , yj , jp , jm , mode
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phigrd

  if(mode==1)then
     do j = 0 , ndy-1
        do i = 0 , ndx-1
           ip=i+1
           im=i-1

           delp = vg(ip,j,1)-vg(i,j,1)
           delm = vg(i,j,1)-vg(im,j,1)
           flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
           Phigrd(i,j) = flmt
        end do
     end do
  end if

  if(mode==3)then
     do j = 0 , ndy-1
        do i = 0 , ndx-1
           jp=j+1
           jm=j-1

           delp = vg(i,jp,2)-vg(i,j,2)
           delm = vg(i,j,2)-vg(i,jm,2)
           flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
           Phigrd(i,j) = flmt
        end do
     end do
  end if
end subroutine vanalbadavg
