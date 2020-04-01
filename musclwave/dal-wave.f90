module comvar
  implicit none
  integer, parameter :: ndx=66,ndy=66,laststep=520000,istx=1,ienx=2,isty=1,ieny=2,svnum=1000 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  integer :: iwx,iwy,iwz
  DOUBLE PRECISION :: cg = 1.0d0 , dx,dy != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0 ,meanrho!,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=66 , ndy2=66 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x
  DOUBLE PRECISION , dimension(-1:ndy2) :: y
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2) ::   Phicgm ,rho, Phi1step , Phi2step ,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2) :: Phidt,Phigrd,Phiexa
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,2) :: source ,sourcedt,sourcedt2
end module grvvar

program muscl1D
  !implicit none :: まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2
  integer :: n,m


  call INITIAL()
  call BC(1)
  !call muslcslv1D(Phi,Phi1step,dt,13)

  do i=1,laststep
     call time(dt)
     write(*,*) i ,'step'



!     call timesource(Phicgp,Phi1step,dt,2)
!     call timesource(Phi1step,rho,dt,1)

!     call timesource(Phicgm,Phi2step,dt,2)
!     call timesource(Phi2step,rho,dt,1)

     !call BC(1)
     !call cnbn(Phicgp,Phicgm)
     !call cnbn(Phi1step,Phi2step)
     !call BC(1)

     !call BC()
     !call muslcslv1D(Phi1step,rho,dt,1)
     !call pbstep(dt)

     !do l=-1,ndz
     do m=-1,ndy
        do n=-1,ndx
           !rho(n,m,l) = U(n,m,l,1)-rhomean
           source(n,m,1)=-source(n,m,1)+G4pi*rho(n,m)
           source(n,m,2)=-source(n,m,2)+G4pi*rho(n,m)
     !      source(n,m,3)=-source(n,m,3)+G4pi*rho(n,m)
     !source(n,m,l,1)=source(n,m,l,1)-G4pi*rho(n,m,l)
     !source(n,m,l,2)=source(n,m,l,2)-G4pi*rho(n,m,l)
     !source(n,m,l,3)=source(n,m,l,3)-G4pi*rho(n,m,l)
        end do
     end do!;end do

     call BC(4)
     call BC(3)
     iwx=1;iwy=0
     call muslcslv1D(Phi1step,source(-1,-1,1),dt*0.5d0,3,2)
     iwx=0;iwy=1
     call muslcslv1D(Phi2step,source(-1,-1,2),dt*0.5d0,3,2)
     !call muslcslv1D(Phi1step,rho,0.5d0*dt,3)
     call BC(4)
     call BC(3)
     iwx=1;iwy=0
     call muslcslv1D(Phi1step,source(-1,-1,1),dt*0.5d0,2,2)
     iwx=0;iwy=1
     call muslcslv1D(Phi2step,source(-1,-1,2),dt*0.5d0,1,2)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,1)

     call BC(4)
     call BC(3)
     call BC(1)
     !Phidt(:)=Phi(:)
     iwx=1;iwy=0
     call muslcslv1D(Phicgp,Phi1step,dt,4,2)
     iwx=0;iwy=1
     call muslcslv1D(Phicgm,Phi2step,dt,4,2)
     call BC(1)
     iwx=1;iwy=0
     call muslcslv1D(Phicgp,Phi1step,dt,1,2)
     iwx=0;iwy=1
     call muslcslv1D(Phicgm,Phi2step,dt,2,2)
     !call BC(1)

     call BC(4)
     call BC(3)
     iwx=1;iwy=0
     !call muslcslv1D(Phi1step,rho,dt,3)
     call muslcslv1D(Phi1step,source(-1,-1,1),dt*0.5d0,3,2)
     iwx=0;iwy=1
     call muslcslv1D(Phi2step,source(-1,-1,2),0.5d0*dt,3,2)
     call BC(4)
     call BC(3)
     iwx=1;iwy=0
     !call muslcslv1D(Phi1step,rho,dt,1)
     call muslcslv1D(Phi1step,source(-1,-1,1),dt*0.5d0,2,2)
     iwx=0;iwy=1
     call muslcslv1D(Phi2step,source(-1,-1,2),dt*0.5d0,1,2)
     !call BC(3)
     !call BC(4)
     !call BC(1)

     if(mod(i,svnum)==1) then
        call saveu(sv)
     end if
  end do
  call BC(3)
  call BC(4)
  call BC(1)
  call saveu(sv)
end program muscl1D

subroutine cnbn(phi1,phi2)
  use comvar
  !double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx) :: phin,phi1,phi2
  integer :: i
  do i = 1 , ndx-2
     phin(i)=(Phi1(i)+phi2(i))*0.5d0
  end do
  phi1(:)=phin(:)
  phi2(:)=phin(:)
end subroutine cnbn

subroutine INITIAL()
  use comvar
  use grvvar
  integer :: i,j
  double precision :: amp,pi=3.1415926535d0,haba,meanphi

  dinit1 = 2.0d0/G4pi/90.d0

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


  !---------Phi-------------
  Phicgp(:,:)=0.0d0
  Phicgm(:,:)=0.0d0
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:,:)=0.0d0
  Phi2step(:,:)=0.0d0
  !Phi1step(:)=+G4pi*meanrho*cg*Lbox
  !Phi2step(:)=0.0d0
  !-------Phi1step-----------

  source(:,:,:)=0.d0
  sourcedt(:,:,:)=0.d0
  sourcedt2(:,:,:)=0.d0

  !-------Phidt-----------
  Phidt(:,:)=0.0d0
  !-------Phdt-----------




  !---------rho-------------
  do i = -1,ndx
     if( dabs(x(i) - hcen) .le. h) then
        rho(i,:) = dinit1
        !rho(i) = 0.0d0
     else
        rho(i,:) = 0.0d0
        !rho(i) = dinit1
        !rho(i) = dinit1*1.d-2
     end if
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


  Phi1step(:,:)=0.d0 !+G4pi*meanrho*cg*Lbox
  Phi2step(:,:)=0.d0 !+G4pi*meanrho*cg*Lbox
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  !dinit1=dinit1-meanrho
  meanphi=0.d0
  open(142,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscle1/phiexact.DAT')
  open(143,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscle1/INIden.DAT')
  open(144,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscle1/phigrd.DAT')
  do j= -1,ndy
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i,j) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
        meanphi=meanphi+G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
     else
        Phiexa(i,j) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
        meanphi=meanphi+G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
     end if
     write(143,*) sngl(rho(i,1))
  end do
  end do
  meanphi=meanphi/dble(ndx+2)/dble(ndy+2)

  !Phiexa(:,:)=Phiexa(:,:)-meanphi

  do j=-1,ndx
  do i=0,ndx-1
     Phigrd(i,j)=(-Phiexa(i-1,j)+Phiexa(i+1,j))*0.5d0/dx
     !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
  end do
  Phigrd(-1,j)=(-Phiexa(0,j)+Phiexa(1,j))/dx
  Phigrd(ndx,j)=(Phiexa(ndx-1,j)-Phiexa(ndx-2,j))/dx
  end do
  !Phigrd(-1,:)=(-Phiexa(0,:)+Phiexa(1,:))/dx
  !Phigrd(ndx,:)=(Phiexa(ndx-1,:)-Phiexa(ndx-2,:))/dx

  !do i=0,ndx-1
  !   Phigrd(i)=-(-Phiexa(i-1)+Phiexa(i+1))*0.5d0/dx
     !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
  !end do
  !Phigrd(-1)=-(-Phiexa(0)+Phiexa(1))/dx
  !Phigrd(ndx)=-(Phiexa(ndx-1)-Phiexa(ndx-2))/dx

!  do i=-1,ndx
!     write(144,*) sngl(x(i)) , Phigrd(i) !, Phiexa(i-1),Phiexa(i+1)
!  end do

!  bcphi1(1) = G4pi * dinit1 * h * dabs(x(1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
!  bcphi2(1) = G4pi * dinit1 * h * dabs(x(ndx-2) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

!  bcphi1(2) = G4pi * dinit1 * h * dabs(x(0) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
!  bcphi2(2) = G4pi * dinit1 * h * dabs(x(ndx-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

!  bcphi1(3) = G4pi * dinit1 * h * dabs(x(-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
!  bcphi2(3) = G4pi * dinit1 * h * dabs(x(ndx) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
  close(142)
  close(143)
  close(144)
  !200 continue
  !--------Phiexa-----------


  !---------wave--------
  goto 201
  !do i = -1, ndx
  !   amp = 1.d-3
  !   Phi(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  !   Phi1step(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  !end do


!  do i = -1, ndx
!     amp = 1.d-3
!     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!     Phi1step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!  end do
  201 continue
  !---------wave--------


  !-------Phidt------------
  goto 302
  !Phi(1)= bcphi1(1)
  !Phi(0)= bcphi1(2)
  !Phi(-1)= bcphi1(3)
  !Phi(ndx-2)= bcphi2(1)
  !Phi(ndx-1)= bcphi2(2)
  !Phi(ndx)= bcphi2(3)
  302 continue
  !-------Phidt------------

  !--------const------------
  !---------Phi-------------
!  Phicgp(:)=bcphi1(1)
!  Phicgm(:)=bcphi1(1)
  !---------Phi-------------

  !-------Phidt-----------
!  Phidt(:)=bcphi1(1)
  !-------Phdt-----------
  !-------Phi1step-----------
  !Phi1step(:)=bcphi1
!  Phi1step(:)= Phigrd(-1)
!  Phi2step(:)= Phigrd(-1)
  !-------Phi1step-----------
  !--------const------------
end subroutine INITIAL



subroutine BC(mode)
  use comvar
  use grvvar
  integer :: i,mode
  double precision , dimension(1:2) :: pl,pr

  if(mode==1) then
     !---------kotei-----------
     !------kotei & free-----------
     !goto 100
     !---------Phi-------------
     !Phicgp(1)= bcphi1(1)
     !Phicgp(0)= bcphi1(2)
     !Phicgp(-1)= bcphi1(3)
     !Phicgp(ndx-2)= bcphi2(1)
     !Phicgp(ndx-1)= bcphi2(2)
     !Phicgp(ndx)= bcphi2(3)
    ! Phicgp(ndx)=Phicgp(ndx-1)
    ! Phicgp(ndx-1)= Phicgp(ndx-2)
    ! Phicgp(ndx-2)= Phicgp(ndx-3)
    ! Phicgp(-1)= Phicgp(0)
    ! Phicgp(0)= Phicgp(1)
    ! Phicgp(1)= Phicgp(2)

    ! Phicgm(ndx)=Phicgm(ndx-1)
    ! Phicgm(ndx-1)= Phicgm(ndx-2)
    ! Phicgm(ndx-2)= Phicgm(ndx-3)
    ! Phicgm(-1)= Phicgm(0)
    ! Phicgm(0)= Phicgm(1)
     ! Phicgm(1)= Phicgm(2)
     !------xbc-----
     Phicgp(ndx,:)=Phiexa(ndx,:)
     Phicgp(ndx-1,:)= Phiexa(ndx-1,:)
     !Phicgp(ndx-2)= Phicgp(ndx-3)
     Phicgp(-1,:)= Phiexa(-1,:)
     Phicgp(0,:)= Phiexa(0,:)
!     Phicgp(ndx,:)=Phiexa(2,:)
!     Phicgp(ndx-1,:)= Phiexa(1,:)
     !Phicgp(ndx-2)= Phicgp(ndx-3)
!     Phicgp(-1,:)= Phiexa(ndx-3,:)
!     Phicgp(0,:)= Phiexa(ndx-2,:)
     !Phicgp(1)= Phicgp(2)

     Phicgm(ndx,:)=Phicgm(ndx-1,:)
     Phicgm(ndx-1,:)= Phicgm(ndx-2,:)
     !Phicgm(ndx-2)= Phicgm(ndx-3)
     Phicgm(-1,:)= Phicgm(0,:)
     Phicgm(0,:)= Phicgm(1,:)
!     Phicgm(ndx,:)=0.d0
!     Phicgm(ndx-1,:)=0.d0
     !Phicgm(ndx-2)= Phicgm(ndx-3)
!     Phicgm(-1,:)=0.d0
!     Phicgm(0,:)=0.d0
     !------xbc-----

     !------ybc-----
     Phicgp(:,ndy)=Phicgp(:,2)
     Phicgp(:,ndy-1)= Phicgp(:,1)
     !Phicgp(ndx-2)= Phicgp(ndx-3)
     Phicgp(:,-1)= Phicgp(:,ndy-3)
     Phicgp(:,0)= Phicgp(:,ndy-2)
     !Phicgp(1)= Phicgp(2)

     Phicgm(:,ndy)=Phicgm(:,2)
     Phicgm(:,ndy-1)= Phicgm(:,1)
     !Phicgm(ndx-2)= Phicgm(ndx-3)
     Phicgm(:,-1)= Phicgm(:,ndy-3)
     Phicgm(:,0)= Phicgm(:,ndy-2)
     !Phicgm(1)= Phicgm(2)
     !------ybc-----

!     Phicgm(1)= bcphi1(1)
!     Phicgm(0)= bcphi1(2)
!     Phicgm(-1)= bcphi1(3)
     !Phicgm(-1)= Phicgm(0)
     !Phicgm(0)= Phicgm(1)
     !Phicgm(1)= Phicgm(2)
!     Phicgm(ndx-2)= bcphi2(1)
!     Phicgm(ndx-1)= bcphi2(2)
!     Phicgm(ndx)= bcphi2(3)


     !Phicgp(1)= bcphi1(1)
     !Phicgp(0)= bcphi1(1)
     !Phicgp(-1)= bcphi1(1)
     !Phicgp(ndx-2)= bcphi2(1)
     !Phicgp(ndx-1)= bcphi2(1)
     !Phicgp(ndx)= bcphi2(1)
     !Phicgp(ndx)=Phicgp(ndx-1)
     !Phicgp(ndx-1)= Phicgp(ndx-2)
     !Phicgp(ndx-2)= Phicgp(ndx-3)

     !Phicgm(1)= bcphi1(1)
     !Phicgm(0)= bcphi1(1)
     !Phicgm(-1)= bcphi1(1)
     !Phicgm(-1)= Phicgm(0)
     !Phicgm(0)= Phicgm(1)
     !Phicgm(1)= Phicgm(2)
     !Phicgm(ndx-2)= bcphi2(1)
     !Phicgm(ndx-1)= bcphi2(1)
     !Phicgm(ndx)= bcphi2(1)


     !---------Phi-------------

     !Phi(1)= bcphi1(1)
     !Phi(0)= bcphi1(1)
     !Phi(-1)= bcphi1(1)
     !Phi(ndx-2)= bcphi2(1)
     !Phi(ndx-1)= bcphi2(1)
     !Phi(ndx)= bcphi2(1)
     !---------Phi-------------
  end if

  if(mode==2) then
     !-------Phi1step-----------
     !Phi1step(1)= bcphi1(1)
     !Phi1step(0)= bcphi1(2)
     !Phi1step(-1)= bcphi1(3)
     !Phi1step(ndx-2)= bcphi2(1)
     !Phi1step(ndx-1)= bcphi2(2)
     !Phi1step(ndx)= bcphi2(3)
     !-------Phi1step-----------
     !100 continue
     !---------kotei-----------

  end if

  if(mode==3)then
     !-------Phi1step+cg-----------
     !goto 700
     !Phi1step(1)= Phigrd(1)
     !Phi1step(0)= Phigrd(0)
     !Phi1step(-1)=Phigrd(-1)
    ! Phi1step(-1)= Phi1step(0)
    ! Phi1step(0)= Phi1step(1)
    ! Phi1step(1)=Phi1step(2)
!     Phi1step(-1,:)= Phi1step(ndx-3)
!     Phi1step(0,:)= Phi1step(ndx-2)
     !----xbcx---
     Phi1step(-1,:)= Phigrd(-1,:)
     Phi1step(0,:)= Phigrd(0,:)

     Phi1step(ndx,:)= Phigrd(ndx,:)
     Phi1step(ndx-1,:)= Phigrd(ndx-1,:)
     !----xbcx---

     !----ybcx---
     Phi1step(:,-1)= Phi1step(:,ndx-3)
     Phi1step(:,0)= Phi1step(:,ndx-2)

     Phi1step(:,ndx)= Phi1step(:,2)
     Phi1step(:,ndx-1)= Phi1step(:,1)
     !----ybcx---
     !Phi1step(ndx-2)= Phigrd(ndx-2)
     !Phi1step(ndx-1)= Phigrd(ndx-1)
     !Phi1step(ndx)= Phigrd(ndx)
   !  Phi1step(ndx)= Phi1step(ndx-1)
   !  Phi1step(ndx-1)= Phi1step(ndx-2)
   !  Phi1step(ndx-2)= Phi1step(ndx-3)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
!     Phi1step(ndx,:)= Phi1step(2)
!     Phi1step(ndx-1,:)= Phi1step(1)
     !700 continue


     !goto 700
     !Phi1step(1)= Phigrd(1)
     !Phi1step(0)= Phigrd(1)
     !Phi1step(-1)=Phigrd(1)
     !Phi1step(-1)= Phi1step(0)
     !Phi1step(0)= Phi1step(1)
     !Phi1step(1)=Phi1step(2)
     !Phi1step(ndx-2)= Phigrd(ndx-2)
     !Phi1step(ndx-1)= Phigrd(ndx-2)
     !Phi1step(ndx)= Phigrd(ndx-2)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
     !700 continue
     !-------Phi1step-----------
  end if

  if(mode==4) then
     !-------Phi1step-cg-----------
     !goto 701
     !Phi1step(1)= -Phigrd(1)
     !Phi1step(0)= -Phigrd(1)
     !Phi1step(-1)=-Phigrd(1)
     !Phi1step(ndx-2)= -Phigrd(ndx-2)
     !Phi1step(ndx-1)= -Phigrd(ndx-2)
     !!Phi1step(ndx)= -Phigrd(ndx-2)
     !Phi2step(1)= -Phigrd(1)
     !Phi2step(0)= -Phigrd(0)
     !Phi2step(-1)=-Phigrd(-1)
     !Phi2step(ndx-2)= -Phigrd(ndx-2)
     !Phi2step(ndx-1)= -Phigrd(ndx-1)
     !Phi2step(ndx)= -Phigrd(ndx)
     !Phi2step(ndx)= -Phi2step(ndx-1)
     !Phi2step(ndx-1)= -Phi2step(ndx-2)
     !Phi2step(ndx-2)= -Phi2step(ndx-3)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
     !701 continue

     !----xbcx---
     Phi2step(-1,:)= Phi2step(0,:)
     Phi2step(0,:)=  Phi2step(1,:)
     Phi2step(ndx,:)=Phi2step(ndx-1,:)
     Phi2step(ndx-1,:)=Phi2step(ndx-2,:)
     !----xbcx---

     !----ybcx---
     Phi2step(:,-1)= Phi2step(:,ndx-3)
     Phi2step(:,0)=  Phi2step(:,ndx-2)
     Phi2step(:,ndx)=Phi2step(:,2)
     Phi2step(:,ndx-1)=Phi2step(:,1)
     !----ybcx---

     !Phi1step(ndx)= -Phigrd(ndx-2)
     !Phi2step(1)= -Phigrd(1)
     !Phi2step(0)= -Phigrd(1)
     !Phi2step(-1)=-Phigrd(1)
     !Phi2step(ndx-2)= -Phigrd(ndx-2)
     !Phi2step(ndx-1)= -Phigrd(ndx-2)
     !Phi2step(ndx)= -Phigrd(ndx-2)
     !-------Phi1step-----------
  end if

!  end if
end subroutine BC


subroutine time(dt)
  use comvar
  use grvvar
  double precision :: dt
  dt = dx/cg * coeff
  write(*,*) 'time cg' , dt
end subroutine time



subroutine timesource(Phiv,source,dt,mode)
  use comvar
  !use grvver
  integer i,mode
  double precision :: dt,sdt,mindt,maxdt , epsl = 1.0d-4
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source

  !mindt=1000.0d0
  maxdt=0.0d0

  if(mode==1) then
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           sdt = 0.5d0*dabs(Phiv(i)) / dabs(cg * G4pi * source(i) )
           !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
           !mindt=dmin1(mindt,sdt)
           maxdt=dmax1(maxdt,sdt)
        end if
     end do
     if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
        dt = sdt
     end if
  end if


  if(mode==2) then
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           sdt = 0.5d0*dabs(Phiv(i)) / dabs( cg * source(i) )
           !sdt = 0.05d0*dabs(Phiv(i)) / ( cg * source(i) )
           !mindt=dmin1(mindt,sdt)
           maxdt=dmax1(maxdt,sdt)
        end if
     end do
     write(*,*) maxdt,'maxdt'
     if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
        dt = sdt
     end if
  end if


  write(*,*) 'time source' , dt
end subroutine timesource


subroutine timesource2(Phiv,Phidt,source,dt)
  use comvar
  !use grvver
  integer i,mode
  double precision :: dt,sdt,mindt,maxdt
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source,Phidt
  do i=1,ndx-2
     if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
        sdt = 0.5d0*dsqrt( dabs(2*Phiv(i)-Phidt(i)) / (cg * cg * G4pi * source(i) ) )
        !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
        !mindt=dmin1(mindt,sdt)
        maxdt=dmax1(maxdt,sdt)
     end if
  end do
  if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
     dt = sdt
  end if


  write(*,*) 'time source2' , dt
end subroutine timesource2

subroutine muslcslv1D(Phiv,source,dt,mode,hazi)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam  !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phigrad,Phipre,fluxphi&
       ,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie,idm
  !DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G

 if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
    if(iwy.eq.1) then; Ncell = ndy; Ncm = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if

  !----kyoukai-----
   if(hazi==1)then
      is = 2
      ie = Ncell-3
   end if
   if(hazi==2)then
      is = 1
      ie = Ncell-2
   end if
  !----kyoukai-----
  nu2 = cg * dt / dx
  Phipre(:,:) = Phiv(:,:)
  !write(name,'(i5.5)') cnt
  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
!     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1 !一次なので大丈夫
              Phi2dt(ix,jy) = Phipre(ix,jy)- 0.5d0 * nu2 * ( Phiu(ix,jy) - Phiu(ixm,jym))
           end do
        end DO
!     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
!      DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is,ie
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy) = Phipre(ix,jy) - nu2 * (Phiu(ix,jy) - Phiu(ixm,jym))
           end do
        end DO
!     end DO
  end if
  !------------ul.solver.+cg-------------



  !------------ul.solver.-cg-------------
  if(mode==2) then

     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,11,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,11)
     !------------calcurate dt/2------------
!     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i=ist-1,ndx-ien+1
              Phi2dt(ix,jy) = Phipre(ix,jy) + 0.5d0 * nu2 * ( Phiu(ixp,jyp) - Phiu(ix,jy))
           end do
        end DO
!     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
!     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
           do i = is,ie
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy) = Phipre(ix,jy) + nu2 * (Phiu(ixp,jyp) - Phiu(ix,jy))
           end do
        end DO
!     end DO

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do

  end if
  !------------ul.solver.-cg-------------


  !--------------source------------------
  if(mode==3) then
!     do k = 1,ndz-2
        do j = 1,ndy-2
           do i = 1,ndx-2
              !Phiv(i,j) =  cg * G4pi * source(i,j) * dt * phiratio + Phipre(i,j)
               Phiv(i,j) =  cg * source(i,j) * dt + Phipre(i,j)
           end do
        end do
!     end do
  end if

  if(mode==4) then
!     do k = 1,ndz-2
        do j = 1,ndy-2
           do i = 1,ndx-2
              !Phiv(i,j) = -cg * source(i,j) * dt + Phipre(i,j)
              Phiv(i,j) = -cg * source(i,j) * dt + Phipre(i,j)
           end do
        end do
!     end do
  end if
  !--------------source------------------

!  close(201)
!  close(202)
  cnt=cnt+2
end subroutine muslcslv1D

!subroutine vanalbada(fg,gradfg,iwx,iwy,iwz)
subroutine vanalbada(Mnum,Lnum,Phipre,Phigrad,i_sta,i_end,dmein)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  !integer :: i , ip , im , flmt ,eps=1.0d-10
  integer :: Mnum,Lnum,Ncell,i_sta,i_end,k,dmein
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm
  integer :: i , ip , im
  !DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy) :: Phipre
  DOUBLE PRECISION, dimension(-1:dmein) :: Phigrad


  !if(iwx.eq.1) Ncell = ndx
  !if(iwy.eq.1) Ncell = ndy
  !if(iwz.eq.1) Ncell = ndz

  do i = i_sta-1 , i_end+1
     ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
     jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!     kz  = iwx*Lnum + iwy*Mnum + iwz*i
     ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
     jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!     kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
     ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
     jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!     kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

     delp = Phipre(ixp,jyp)-Phipre(ix,jy)
     delm = Phipre(ix,jy)-Phipre(ixm,jym)
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !Phigrad(ix,jy,kz) = flmt
     Phigrad(i) = flmt
  end do

end subroutine vanalbada


subroutine fluxcal(preuse,pre,uin,ep,kappa,mode,is,ie)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy) :: ul,ur,pre,preuse,uin
  DOUBLE PRECISION , dimension(-1:ndx) :: slop  !------------- need allocation --------------
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  !DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  !uin(:)=0.0d0
!  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
!     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
!        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndx;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              ul(ix,jy) = preuse(ix,jy) + 0.25d0 * ep * slop(i) &
                   * ((1.0d0-slop(i)*kappa)*(pre(ix,jy)-pre(ixm,jym)) + &
                   (1.0d0+slop(i)*kappa)*(pre(ixp,jyp) - pre(ix,jy))) !i+1/2
              uin(ix,jy)=ul(ix,jy)
              end do
              end DO
!              end DO
              !write(*,*) slop(127),'127slop'
              !uin(:)=ul(:)
           end if


           if(mode==4) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-1,ndx-ien+1
              ur(ix,jy) = preuse(ix,jy) - 0.25d0 * ep * slop(i) &
                   * ((1.0d0+slop(i)*kappa)*(pre(ix,jy)-pre(ixm,jym)) + &
                   (1.0d0-slop(i)*kappa)*(pre(ixp,jyp) - pre(ix,jy))) !i-1/2
              uin(ix,jy)=ur(ix,jy)
              end do
              end DO
!              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !uin(:)=ur(:)
           end if

           if(mode==10) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum! + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-2,ndx-ien+2
              ul(ix,jy) = preuse(ix,jy)
              uin(ix,jy)=ul(ix,jy)
              end do
              end DO
!              end DO
           end if

           if(mode==11) then
!              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
              ix  = iwx*i    + iwy*Mnum! + iwz*Mnum
              jy  = iwx*Mnum + iwy*i   ! + iwz*Lnum
!              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Mnum! + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)!+ iwz*Lnum
!              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Mnum !+ iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)!+ iwz*Lnum
!              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = is,ie
              ur(ix,jy) = preuse(ix,jy)
              uin(ix,jy)=ur(ix,jy)
              end do
              end DO
!              end DO
           end if


end subroutine fluxcal

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscle1/phi'//name//'.dat')
  do j=1,ndy-2
     do i=1,ndx-2
        write(21,*) x(i),y(j), Phicgp(i,j),Phi1step(i,j) , Phicgm(i,j),Phi2step(i,j) ,rho(i,j),&
             Phiexa(i,j),Phigrd(i,j)
     end do
     write(21,*)
  end do
  close(21)
   !i=1
   !write(*,*) Phi1step(i),-Phi2step(i),Phigrd(i)
  !i=1
  !write(*,*)-((Phicgp(i-1)+Phicgp(i-1))*0.5d0 - (Phicgp(i+1)+Phicgp(i+1))*0.5d0)*0.5d0/dx,&
  !     -((Phicgm(i-1)+Phicgm(i-1))*0.5d0 - (Phicgm(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx,&
  !     Phigrd(i)
  write(*,*) 'save step : ',in1
  in1=in1+1
end subroutine saveu




subroutine dfi2(Phiv,Phidt,U,dt,mode)
  use comvar
  integer i,mode
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,Phidummy,Phidt,U
  double precision nu2,dt

  nu2=dt*dt*cg*cg/dx/dx
  Phidummy(:)=Phiv(:)
  if(mode==1) then
     do i=ist,ndx-ien
        Phiv(i)= 2.0d0*Phidummy(i) - Phidt(i) + nu2 * (Phidummy(i-1) + Phidummy(i+1)  &
             - w1 * Phidummy(i)) - U(i) * G4pi * dt * dt * cg * cg
     end do
  end if

  if(mode==2) then
     do i=ist,ndx-ien
        Phiv(i)= Phidummy(i+1) + Phidummy(i-1) - Phidt(i) + nu2 * (Phidummy(i-1) + Phidummy(i+1)  &
             - w1 * Phidummy(i)) - U(i) * G4pi * dt * dt * cg * cg
     end do
  end if

  if(mode==3) then
     do i=ist,ndx-ien
        Phiv(i)= Phidummy(i+1) + Phidummy(i-1) - 0.5d0 * Phidt(i+1) - 0.5d0 * Phidt(i-1) + nu2 * (Phidummy(i-1) + Phidummy(i+1)  &
             - w1 * Phidummy(i)) - U(i) * G4pi * dt * dt * cg * cg
     end do
  end if
  Phidt(:)=Phidummy(:)
end subroutine dfi2


subroutine split(Phiv,Phidt,Phi1step,dt,mode)
  use comvar
  integer :: i,mode,cnt=0
  double precision :: dt,dtpre !前の
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,Phi1step,Phidt
  character(5) name


  !sp1 = 0.5d0*(Phiv(i) + )
end subroutine split



subroutine pbstep(dt)
USE comvar
!USE mpivar
USE grvvar
!INCLUDE 'mpif.h'
integer :: i,j,k,c=0
double precision :: dxx,dt,prdt=0.0d0,ratio=dsqrt(2.0d0),dddt,rx,ry,rz,rx2,ry2,rz2
character(3) fn
character(5) ciii
!DOUBLE PRECISION , dimension(-1:ndx,1:Dim) :: slop
DOUBLE PRECISION flx1,fly1,flz1,flx2,fly2,flz2,flx3,fly3,flz3,bc11,bc12,bc21,bc22
DOUBLE PRECISION fllx1,flly1,fllz1,fllx2,flly2,fllz2,fllx3,flly3,fllz3

!dxx=dx(1)
!dxx=-dx(1)
!dxx=-deltalength
dxx=dx
!iwx = 0; iwy = 1; iwz = 1
!call BCgrv(101)
!call BCgrv(102)
!bphi1l(j,k,1)
write(fn,'(i3.3)') NRANK
write(ciii,'(i5.5)') c
!open(3,file='/work/maedarn/3DMHD/test/bcstep'//fn//'.DAT')
!open(4,file='/work/maedarn/3DMHD/test/bcstepscnd'//fn//'.DAT')
!open(5,file='/work/maedarn/3DMHD/test/source'//fn//ciii//'.DAT')


sourcedt2(:,:,:)=sourcedt(:,:,:)
sourcedt(:,:,1)=Phi1step(:,:)
sourcedt(:,:,2)=Phi2step(:,:)
source(:,:,:)=0.0d0

dddt=1.0d0/dt/ratio/cg
write(*,*) 'dddt',dddt

!do
!call vanalbada(Mnum,Lnum,Phipre,Phigrad,i_sta,i_end,dmein)

!do k = 1,ndz-2
do j = 1,ndy-2
do i = 1,ndx-2

rx = (sourcedt(i,j,1)-sourcedt(i-1,j,1)+1.d-20)/(sourcedt(i+1,j,1)-sourcedt(i,j,1)+1.d-20)
ry = (sourcedt(i,j,2)-sourcedt(i,j-1,2)+1.d-20)/(sourcedt(i,j+1,2)-sourcedt(i,j,2)+1.d-20)
!rz = (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3)+1.d-20)/(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3)+1.d-20)
rx2 = (sourcedt(i-1,j,1)-sourcedt(i-2,j,1)+1.d-20)/(sourcedt(i,j,1)-sourcedt(i-1,j,1)+1.d-20)
ry2 = (sourcedt(i,j-1,2)-sourcedt(i,j-2,2)+1.d-20)/(sourcedt(i,j,2)-sourcedt(i,j-1,2)+1.d-20)
!rz2 = (sourcedt(i,j,k-1,3)-sourcedt(i,j,k-2,3)+1.d-20)/(sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3)+1.d-20)

flx1=dmax1(0.0d0,dmin1(1.0d0,rx))
fly2=dmax1(0.0d0,dmin1(1.0d0,ry))
!flz3=dmax1(0.0d0,dmin1(1.0d0,rz))
fllx1=dmax1(0.0d0,dmin1(1.0d0,rx2))
flly2=dmax1(0.0d0,dmin1(1.0d0,ry2))
!fllz3=dmax1(0.0d0,dmin1(1.0d0,rz2))


source(i,j,1)= (sourcedt(i,j,2)-sourcedt2(i,j,2))*dddt! &
    ! +1.0d0/dx*0.5d0*((sourcedt(i,j+1,2)-sourcedt(i,j,2))*fly2+(sourcedt(i,j,2)-sourcedt(i,j-1,2))*flly2 &
    ! )
source(i,j,2)= (sourcedt(i,j,1)-sourcedt2(i,j,1))*dddt! &
    ! +1.0d0/dx*0.5d0*((sourcedt(i+1,j,1)-sourcedt(i,j,1))*flx1+(sourcedt(i,j,1)-sourcedt(i+1,j,1))*fllx1 &
    !)
!source(i,j,k,3)= (sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt+(sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2+(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
!     +(sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1 + (sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 )

end do
end do
!end do

end subroutine pbstep

