module comvar
  implicit none
  integer, parameter :: ndx=130,ndy=130,ndz=130,laststep=15000,ist=2,ien=3 , istp=1,ienp=2!preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 :: ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  DOUBLE PRECISION :: cg = 1.0d0 , dx ,dy,dz!= Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  integer iwx,iwy,iwz,ifEVO,ifEVO2
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.1d0 !,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3,-1:ndy,-1:ndz) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=130,ndy2=130,ndz2=130 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x,y,z
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,-1:ndz2) ::  Phicgm ,rho, Phi1step , Phi2step ,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2,-1:ndy2,-1:ndz2) :: Phidt,Phigrd,Phiexa
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2,ifEVO=1


  call INITIAL()
  !call BC(1)
  !call BC(8)
  !call BC(28)
  !call BC(9)
  !call BC(29)
  call cllsub(1,dt)
  !call muslcslv1D(Phi,Phi1step,dt,13)
  ifEVO=1
  ifEVO2=1
  do i=1,laststep
     write(*,*) i,'timestep'
     call cllsub(2,dt)
     !call time(dt)


     !call timesource(Phicgp,Phi1step,dt,2)
     !call timesource(Phi1step,rho,dt,1)

     !call timesource(Phicgm,Phi2step,dt,2)
     !call timesource(Phi2step,rho,dt,1)

     !call BC(1)
     !call cnbn(Phicgp,Phicgm)
     !call cnbn(Phi1step,Phi2step)
     !call BC(1)

     !call BC()
     !call muslcslv1D(Phi1step,rho,dt,1)
     !call BC(4)
     !call BC(3)
     !call BC(18)
     !call BC(38)
     !call BC(19)
     !call BC(39)
     call cllsub(3,dt)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,3)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,3)
     !call muslcslv1D(Phi1step,rho,0.5d0*dt,3)
     call cllsub(3,dt)
     !call BC(4)
     !call BC(3)
     !call BC(18)
     !call BC(38)
     !call BC(19)
     !call BC(39)
!     if(ifEVO.eq.1) then
!        iwx=1; iwy=0; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=1; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=0; iwz=1
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        ifEVO = 2
!        goto 1000
!     end if
!     if(ifEVO.eq.2) then
!        iwx=0; iwy=1; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=0; iwz=1
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=1; iwy=0; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        ifEVO = 3
!        goto 1000
!     end if
!     if(ifEVO.eq.3) then
!        iwx=0; iwy=0; iwz=1
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=1; iwy=0; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!       call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=1; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        ifEVO = 4
!        goto 1000
!     end if
!     if(ifEVO.eq.4) then
!        iwx=1; iwy=0; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=0; iwz=1
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=1; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        ifEVO = 5
!        goto 1000
!     end if
!     if(ifEVO.eq.5) then
!        iwx=0; iwy=1; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=1; iwy=0; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=0; iwz=1
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        ifEVO = 6
!        goto 1000
!     end if
!     if(ifEVO.eq.6) then
!        iwx=0; iwy=0; iwz=1
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=0; iwy=1; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        iwx=1; iwy=0; iwz=0
!        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
!        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
!        ifEVO = 1
!        goto 1000
!     end if
!     1000 continue
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
     !call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,1)
     call cllsub(4,dt)

!     call BC(4)
!     call BC(3)
!     call BC(18)
!     call BC(38)
!     call BC(19)
!     call BC(39)
!     call BC(1)
!     call BC(8)
!     call BC(28)
!    call BC(9)
!     call BC(29)
     !Phidt(:)=Phi(:)
     call cllsub(3,dt)
     call cllsub(1,dt)
     call muslcslv1D(Phicgp,Phi1step,dt,4)
     call muslcslv1D(Phicgm,Phi2step,dt,4)
     call cllsub(1,dt)
     !call BC(1)
     !call BC(8)
     !call BC(28)
     !call BC(9)
     !call BC(29)
     !call muslcslv1D(Phicgp,Phi1step,dt,1)
     !call muslcslv1D(Phicgm,Phi2step,dt,2)
     call cllsub(5)
     !call BC(1)

     !call BC(4)
     !call BC(3)
     !call BC(18)
     !call BC(38)
     !call BC(19)
     !call BC(39)
     call cllsub(3,dt)
     !call muslcslv1D(Phi1step,rho,dt,3)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,3)
     call muslcslv1D(Phi2step,rho,0.5d0*dt,3)
     !call BC(4)
     !call BC(3)
     !call BC(18)
     !call BC(38)
     !call BC(19)
     !call BC(39)
     !call muslcslv1D(Phi1step,rho,dt,1)
     call cllsub(4,dt)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
     !call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
     !call BC(3)
     !call BC(4)
     !call BC(18)
     !call BC(38)
     !call BC(19)
     !call BC(39)
     ifEVO=1+mod(i,6)
     ifEVO2=1+mod(i,6)
     call saveu(sv)
  end do

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
  integer :: i,j,k
  double precision :: amp,pi=3.1415926535d0,haba

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

  !----------z--------------
  dz = Lbox/dble(ndz-2)
  z(1) = dz/2.0d0
  z(0) = z(1) - dz
  z(-1) = z(0) - dz
  do i=2,ndz
     z(i) = z(i-1) + dz
  end do
  !----------z--------------


  !---------Phi-------------
  Phicgp(:,:,:)=0.0d0
  Phicgm(:,:,:)=0.0d0
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:,:,:)=0.0d0
  Phi2step(:,:,:)=0.0d0
  !-------Phi1step-----------

  !-------Phidt-----------
  !Phidt(:)=0.0d0
  !-------Phdt-----------




  !---------rho-------------
  do k=-1,ndz
     do j=-1,ndy
        do i = -1,ndx
           if( dabs(x(i) - hcen) .le. h) then
              !if( dabs(y(i) - hcen) .le. h) then
              !if( dabs(z(i) - hcen) .le. h) then
              rho(i,j,k) = dinit1
              !rho(i) = 0.0d0
           else
              rho(i,j,k) = 0.0d0
              !rho(i) = dinit1*1.d-2
           end if
        end do
     end do
  end do
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  open(142,file='/Users/maeda/Desktop/kaiseki/testcode4/phiexact.DAT')
  open(143,file='/Users/maeda/Desktop/kaiseki/testcode4/INIden.DAT')
  open(144,file='/Users/maeda/Desktop/kaiseki/testcode4/phigrd.DAT')

  do k=-1,ndz
     do j=-1,ndy
        do i= -1,ndx
           if( dabs(x(i) - hcen) .le. h ) then
              !if( dabs(y(i) - hcen) .le. h) then
              !if( dabs(z(i) - hcen) .le. h) then
              Phiexa(i,j,k) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
              write(142,*) sngl(x(i)),sngl(y(j)),sngl(z(k)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
           else
              Phiexa(i,j,k) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
              write(142,*) sngl(x(i)) ,sngl(y(j)),sngl(z(k)) ,&
                   sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
           end if
           write(143,*) sngl(rho(i,j,k))
        end do
     end do
  end do

  !--- x-exact -----
  do k=-1,ndz
     do j=-1,ndy
        do i=0,ndx-1
           Phigrd(i,j,k)=(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx
           !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
        end do
        Phigrd(-1,j,k)=(-Phiexa(0,j,k)+Phiexa(1,j,k))/dx
        Phigrd(ndx,j,k)=(Phiexa(ndx-1,j,k)-Phiexa(ndx-2,j,k))/dx
     end do
  end do
        !Phigrd(-1)=(-Phiexa(0)+Phiexa(1))/dx
        !Phigrd(ndx)=(Phiexa(ndx-1)-Phiexa(ndx-2))/dx

  do k=1,ndz-2
      do j=1,ndy-2
         do i=1,ndx-2
            write(144,*) sngl(x(i)),sngl(y(j)),sngl(z(k)) , sngl(Phigrd(i,j,k)) !, Phiexa(i-1,j,k),Phiexa(i+1)
         end do
      end do
   end do

   do k=-1,ndz
      do j=-1,ndy

         bcphi1(1,j,k) = G4pi * dinit1 * h * dabs(x(1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
         bcphi2(1,j,k) = G4pi * dinit1 * h * dabs(x(ndx-2) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

         bcphi1(2,j,k) = G4pi * dinit1 * h * dabs(x(0) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
         bcphi2(2,j,k) = G4pi * dinit1 * h * dabs(x(ndx-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

         bcphi1(3,j,k) = G4pi * dinit1 * h * dabs(x(-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
         bcphi2(3,j,k) = G4pi * dinit1 * h * dabs(x(ndx) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
      end do
   end do
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


  do i = -1, ndx
     amp = 1.d-3
     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     !Phi1step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
  end do
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
  Phicgp(:,:,:)=bcphi1(1,1,1)
  Phicgm(:,:,:)=bcphi1(1,1,1)
  !---------Phi-------------

  !-------Phidt-----------
  !Phidt(:)=bcphi1(1)
  !-------Phdt-----------
  !-------Phi1step-----------
  !Phi1step(:)=bcphi1
  !-------Phi1step-----------
  !--------const------------
end subroutine INITIAL



subroutine BC(mode)
  use comvar
  use grvvar
  integer :: i,mode,j,k
  double precision , dimension(1:2) :: pl,pr

  if(mode==1) then
     !---------kotei-x-----------
     !goto 100
     !---------Phi-------------
     Phicgp(1,:,:)= bcphi1(1,:,:)
     Phicgp(0,:,:)= bcphi1(2,:,:)
     Phicgp(-1,:,:)= bcphi1(3,:,:)
     Phicgp(ndx-2,:,:)= bcphi2(1,:,:)
     Phicgp(ndx-1,:,:)= bcphi2(2,:,:)
     Phicgp(ndx,:,:)= bcphi2(3,:,:)

     Phicgm(1,:,:)= bcphi1(1,:,:)
     Phicgm(0,:,:)= bcphi1(2,:,:)
     Phicgm(-1,:,:)= bcphi1(3,:,:)
     Phicgm(ndx-2,:,:)= bcphi2(1,:,:)
     Phicgm(ndx-1,:,:)= bcphi2(2,:,:)
     Phicgm(ndx,:,:)= bcphi2(3,:,:)

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
     !---------kotei-x-----------
     !-------Phi1step+cg-----------
     !goto 700
     Phi1step(1,:,:)= Phigrd(1,:,:)
     Phi1step(0,:,:)= Phigrd(0,:,:)
     Phi1step(-1,:,:)=Phigrd(-1,:,:)
     Phi1step(ndx-2,:,:)= Phigrd(ndx-2,:,:)
     Phi1step(ndx-1,:,:)= Phigrd(ndx-1,:,:)
     Phi1step(ndx,:,:)= Phigrd(ndx,:,:)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
     !700 continue
     !-------Phi1step-----------
     !---------kotei-x-----------
  end if

  if(mode==4) then
     !---------kotei-x-----------
     !-------Phi1step-cg-----------
     !goto 701
     !Phi1step(1)= -Phigrd(1)
     !Phi1step(0)= -Phigrd(1)
     !Phi1step(-1)=-Phigrd(1)
     !Phi1step(ndx-2)= -Phigrd(ndx-2)
     !Phi1step(ndx-1)= -Phigrd(ndx-2)
     !Phi1step(ndx)= -Phigrd(ndx-2)
     Phi2step(1,:,:)= -Phigrd(1,:,:)
     Phi2step(0,:,:)= -Phigrd(2,:,:)
     Phi2step(-1,:,:)=-Phigrd(3,:,:)
     Phi2step(ndx-2,:,:)= -Phigrd(ndx-2,:,:)
     Phi2step(ndx-1,:,:)= -Phigrd(ndx-1,:,:)
     Phi2step(ndx,:,:)= -Phigrd(ndx,:,:)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
     !701 continue
     !-------Phi1step-----------
     !---------kotei-x-----------
  end if


  !----- y-periodic ------
  if(mode==8) then
     !-------period2-----------
     !     goto 102
     do k=-1,ndz
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgp(i,ndx-2,k)
           pr(1)= Phicgp(i,ndx-3,k)
           pl(1)= Phicgp(i,1,k)
           pl(2)= Phicgp(i,2,k)
           Phicgp(i,1,k)=pr(1)
           Phicgp(i,2,k)=pr(2)
           Phicgp(i,ndx-2,k)=pl(1)
           Phicgp(i,ndx-3,k)=pr(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==18) then
     do k=-1,ndz
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi1step(i,ndx-2,k)
           pr(1)= Phi1step(i,ndx-3,k)
           pl(1)= Phi1step(i,1,k)
           pl(2)= Phi1step(i,2,k)
           Phi1step(i,1,k)=pr(1)
           Phi1step(i,2,k)=pr(2)
           Phi1step(i,ndx-2,k)=pl(1)
           Phi1step(i,ndx-3,k)=pr(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if




  if(mode==28) then
     !-------period2-----------
     !     goto 102
     do k=-1,ndz
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgm(i,ndx-2,k)
           pr(1)= Phicgm(i,ndx-3,k)
           pl(1)= Phicgm(i,1,k)
           pl(2)= Phicgm(i,2,k)
           Phicgm(i,1,k)=pr(1)
           Phicgm(i,2,k)=pr(2)
           Phicgm(i,ndx-2,k)=pl(1)
           Phicgm(i,ndx-3,k)=pr(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==38) then
     do k=-1,ndz
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi2step(i,ndx-2,k)
           pr(1)= Phi2step(i,ndx-3,k)
           pl(1)= Phi2step(i,1,k)
           pl(2)= Phi2step(i,2,k)
           Phi2step(i,1,k)=pr(1)
           Phi2step(i,2,k)=pr(2)
           Phi2step(i,ndx-2,k)=pl(1)
           Phi2step(i,ndx-3,k)=pr(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if


  !----- y-periodic ------


  !----- z-periodic ------
  if(mode==9) then
     !-------period2-----------
     !     goto 102
     do j=-1,ndy
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgp(i,j,ndx-2)
           pr(1)= Phicgp(i,j,ndx-3)
           pl(1)= Phicgp(i,j,1)
           pl(2)= Phicgp(i,j,2)
           Phicgp(i,j,1)=pr(1)
           Phicgp(i,j,2)=pr(2)
           Phicgp(i,j,ndx-2)=pl(1)
           Phicgp(i,j,ndx-3)=pr(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==19) then
     do j=-1,ndy
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi1step(i,j,ndx-2)
           pr(1)= Phi1step(i,j,ndx-3)
           pl(1)= Phi1step(i,j,1)
           pl(2)= Phi1step(i,j,2)
           Phi1step(i,j,1)=pr(1)
           Phi1step(i,j,2)=pr(2)
           Phi1step(i,j,ndx-2)=pl(1)
           Phi1step(i,j,ndx-3)=pr(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if




  if(mode==29) then
     !-------period2-----------
     !     goto 102
     do j=-1,ndy
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgm(i,j,ndx-2)
           pr(1)= Phicgm(i,j,ndx-3)
           pl(1)= Phicgm(i,j,1)
           pl(2)= Phicgm(i,j,2)
           Phicgm(i,j,1)=pr(1)
           Phicgm(i,j,2)=pr(2)
           Phicgm(i,j,ndx-2)=pl(1)
           Phicgm(i,j,ndx-3)=pr(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==39) then
     do j=-1,ndy
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi2step(i,j,ndx-2)
           pr(1)= Phi2step(i,j,ndx-3)
           pl(1)= Phi2step(i,j,1)
           pl(2)= Phi2step(i,j,2)
           Phi2step(i,j,1)=pr(1)
           Phi2step(i,j,2)=pr(2)
           Phi2step(i,j,ndx-2)=pl(1)
           Phi2step(i,j,ndx-3)=pr(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if


  !----- z-periodic ------




  if(mode==10) then
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

     !--------free--------------
     !goto 112
     !-------Phi1step-----------
!     Phi1step(1)= Phi1step(2)
!     Phi1step(0)= Phi1step(1)
!     Phi1step(-1)=Phi1step(0)
!     Phi1step(ndx-2)= bcphi2(1)
!     Phi1step(ndx-1)= bcphi2(2)
!     Phi1step(ndx)= bcphi2(3)
     !-------Phi1step-----------
!130  continue
     !--------free--------------
  end if
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
  integer i,mode,j,k
  double precision :: dt,sdt,mindt,maxdt , epsl = 1.0d-4
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phiv,source

  !mindt=1000.0d0
  maxdt=0.0d0

  if(mode==1) then
     do k=1,ndx-2
        do j=1,ndx-2
           do i=1,ndx-2
              if((source(i,j,k) .ne. 0.0d0) .and. (Phiv(i,j,k) .ne. 0.0d0))then
                 sdt = 0.5d0*dabs(Phiv(i,j,k)) / (cg * G4pi * source(i,j,k) )
                 !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
                 !mindt=dmin1(mindt,sdt)
                 maxdt=dmax1(maxdt,sdt)
              end if
           end do
        end do
     end do
     if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
        dt = sdt
     end if
  end if


  if(mode==2) then
     do k=1,ndx-2
        do j=1,ndx-2
           do i=1,ndx-2
              if((source(i,j,k) .ne. 0.0d0) .and. (Phiv(i,j,k) .ne. 0.0d0))then
                 sdt = 0.5d0*dabs(Phiv(i,j,k)) / ( cg * source(i,j,k) )
                 !sdt = 0.05d0*dabs(Phiv(i)) / ( cg * source(i) )
                 !mindt=dmin1(mindt,sdt)
                 maxdt=dmax1(maxdt,sdt)
              end if
           end do
        end do
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


subroutine muslcslv1D(Phiv,source,dt,mode)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0,iwx,iwy,iwz
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  !DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre,fluxphi&
       ,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum

  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy; endif! BT1 = 1; BT2 = 2; VN = 4; end if

  nu2 = cg * dt / dx
  Phipre(:,:,:) = Phiv(:,:,:)
  !write(name,'(i5.5)') cnt

  !ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  !jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  !kz  = iwx*Lnum + iwy*Mnum + iwz*i
  !ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  !jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  !kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)


  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1 !一次なので大丈夫
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz)- 0.5d0 * nu2 * ( Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
      DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) - nu2 * (Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
           end do
        end DO
     end DO
     !write(*,*) Phiv(127),'127-3'
     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do
  end if
  !------------ul.solver.+cg-------------



  !------------ul.solver.-cg-------------
  if(mode==2) then
     !write(name,'(i5.5)') cnt
     !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/2modepre'//name//'.dat')
     !write(name,'(i5.5)') cnt+1
     !open(202,file='/Users/maeda/Desktop/kaiseki/testcode2/2modepos'//name//'.dat')
     !do i=-1,ndx
     !   write(201,*) i, Phiv(i)
     !end do

     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,11)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,11)
     !------------calcurate dt/2------------
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz) + 0.5d0 * nu2 * ( Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) + nu2 * (Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
           end do
        end DO
     end DO

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do

  end if
  !------------ul.solver.-cg-------------


  !--------------source------------------
  if(mode==3) then
     !write(name,'(i5.5)') cnt
     !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/3modepre'//name//'.dat')
     !write(name,'(i5.5)') cnt+1
     !open(202,file='/Users/maeda/Desktop/kaiseki/testcode2/3modepos'//name//'.dat')
     !do i=-1,ndx
     !   write(201,*) i, Phiv(i)
     !end do
     !write(*,*) 'in1'
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i=ist,ndx-ien
              Phiv(i,j,k) =  -cg * G4pi * source(i,j,k) * dt + Phipre(i,j,k)
           end do
        end do
     end do
     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do
  end if

  if(mode==4) then
     !write(name,'(i5.5)') cnt
     !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/4modepre'//name//'.dat')
     !write(name,'(i5.5)') cnt+1
     !open(202,file='/Users/maeda/Desktop/kaiseki/testcode2/4modepos'//name//'.dat')
     !do i=-1,ndx
     !   write(201,*) i, Phiv(i)
     !end do
     !write(*,*) source(ndx-ien), dt , Phipre(ndx-ien),cg * source(ndx-ien) * dt + Phipre(ndx-ien),'mode4'
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i=ist,ndx-ien
              Phiv(i,j,k) = cg * source(i,j,k) * dt + Phipre(i,j,k)
           end do
        end do
     end do
     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do
  end if
  !--------------source------------------



  !--------------sourcetest------------------
  if(mode==5) then
     do i=ist,ndx-ien
        Phiv(i) = cg * 0.5d0* (source(i+1) + source(i-1)) * dt + Phipre(i)
     end do
  end if

  if(mode==6) then
     !write(*,*) source(ndx-ien), dt , Phipre(ndx-ien),cg * source(ndx-ien) * dt + Phipre(ndx-ien),'mode4'
     sourcepri(:)=source(:)
     do i=ist,ndx-ien
        Phiv(i) = cg * source(i) * dt + Phipre(i) - 0.5d0 * cg * cg * dt * (source(i) - sourcepre(i))
     end do
     sourcepre(:)=sourcepri(:)
  end if
  !--------------sourcetest------------------

  if(mode==13) then
     sourcepre(:)=0.0d0
  end if

  close(201)
  close(202)
  cnt=cnt+2
end subroutine muslcslv1D

!subroutine vanalbada(fg,gradfg,iwx,iwy,iwz)
subroutine vanalbada(Mnum,Lnum,Phipre,Phigrad,i_sta,i_end,dmin)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  !integer :: i , ip , im , flmt ,eps=1.0d-10
  integer :: Mnum,Lnum,Ncell,i_sta,i_end,k,dmin
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm
  integer :: i , ip , im
  !DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phipre
  DOUBLE PRECISION, dimension(-1:dmin) :: Phigrad

!  do i = ist-1 , ndx - ien + 1
!     ip=i+1
!     im=i-1

!     delp = Phipre(ip)-Phipre(i)
!     delm = Phipre(i)-Phipre(ixm,jym,kzm)
!     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !flmt = (2.d0*delp*delm+eps)/(delp**2+delm**2+eps)
     !if(i==58) then
     !   write(*,*) delp , delm ,Phipre(ip),Phipre(i),Phipre(im) , flmt
     !end if
     !flmt = (2.d0*delp*delm+eps)/(delp**2+delm**2+eps)
!     Phigrad(i) = flmt
     !grdU(i,k) = flmt*( wave(ixp,jyp,kzp)-wave(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )
     !Phigrad(i) = flmt*( Phipre(ip)-Phipre(im) )/( 2.0d0 * dx )
  end do

  if(iwx.eq.1) Ncell = ndx
  if(iwy.eq.1) Ncell = ndy
  if(iwz.eq.1) Ncell = ndz

  do i = i_sta, Ncell+i_end
     ix  = iwx*i    + iwy*Lnum + iwz*Mnum
     jy  = iwx*Mnum + iwy*i    + iwz*Lnum
     kz  = iwx*Lnum + iwy*Mnum + iwz*i
     ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
     jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
     kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
     ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
     jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
     kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

     delp = Phipre(ixp,jyp,kzp)-Phipre(ix,jy,kz)
     delm = Phipre(ix,jy,kz)-Phipre(ixm,jym,kzm)
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !Phigrad(ix,jy,kz) = flmt
     Phigrad(i) = flmt
     !grdU(i,k) = flmt*( U(ixp,jyp,kzp,k)-U(ixm,jym,kzm,k) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

     !T = 1.27d0*U(ix,jy,kz,5)/(kb*U(ix,jy,kz,1))
     !delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
     !gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
     !grdU(i,k) = grdU(i,k)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))

  end do

end subroutine vanalbada

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file='/Users/maeda/Desktop/kaiseki/testcode4/phi'//name//'.dat')
  do i=1,ndx-2
     write(21,*) x(i), Phicgp(i),Phi1step(i) , Phicgm(i),Phi2step(i) ,&
          (Phicgp(i)+Phicgm(i))*0.5d0,(Phi1step(i)+Phi2step(i))*0.5d0, (Phicgp(i)-Phicgm(i)),&
          (Phicgp(i)+Phicgm(i))*0.5d0+dabs(Phicgp(i)-Phicgm(i))
  end do
  close(21)
  in1=in1+1
end subroutine saveu



subroutine fluxcal(preuse,pre,u,ep,kappa,mode,is,ie)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,-1:ndz) :: ul,ur,pre,preuse,u
  DOUBLE PRECISION , dimension(-1:dmin) :: slop
  integer :: i,mode,Ncell,Ncl,Ncm,dmin,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  !u(:)=0.0d0
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(pre,slop)
              do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              ul(ix,jy,kz) = preuse(ix,jy,kz) + 0.25d0 * ep * slop(i) &
                   * ((1.0d0-slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + (1.0d0+slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i+1/2
              u(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !u(:)=ul(:)
           end if


           if(mode==4) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(pre,slop)
              do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-1,ndx-ien+1
              ur(ix,jy,kz) = preuse(ix,jy,kz) - 0.25d0 * ep * slop(i) &
                   * ((1.0d0+slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + (1.0d0-slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i-1/2
              u(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !u(:)=ur(:)
           end if

           if(mode==10) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(pre,slop)
              do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-2,ndx-ien+2
              ul(ix,jy,kz) = preuse(ix,jy,kz)
              u(ix,jy,kz)=ul(ix,jy,kz)
              end do
           end if

           if(mode==11) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(pre,slop)
              do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = is,ie
              ur(ix,jy,kz) = preuse(ix,jy,kz)
              u(ix,jy,kz)=ur(ix,jy,kz)
              end do
           end if

           !if(mode==12) then
           !   do i = ist-1,ndx-ien+1
           !      ur(i) = preuse(i+1)
           !      u(i)=ur(i)
           !   end do
           !end if



  !-----------test--------------
  goto 400
  if(mode==2) then
     do i = ist-1,ndx-ien+1
        ur(i) = preuse(i+1) - 0.25d0 * ep * slop(i) &
             * ((1.0d0+slop(i)*kappa)*(pre(i+1)-pre(i)) + (1.0d0-slop(i)*kappa)*(pre(i+2) - pre(i+1))) !i+1/2
     end do
     u(:)=ur(:)
  end if


  if(mode==3) then
     do i = ist-1,ndx-ien+1
        ur(i) = preuse(i+1) - 0.25d0 * ep * slop(i+1) &
             * ((1.0d0+slop(i+1)*kappa)*(pre(i+1)-pre(i)) + (1.0d0-slop(i+1)*kappa)*(pre(i+2) - pre(i+1))) !i+1/2
     end do
     u(:)=ur(:)
  end if
  400 continue
  !-----------test--------------
end subroutine fluxcal



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

  write(name,'(i5.5)') cnt
  !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/pcg'//name//'.dat')

  !----------- +cg -----------
  if(mode==1) then
     !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/pcg'//name//'.dat')
     !do i=ist-1,ndx-ien+1
     do i=ist,ndx-ien
!     do i=ist+2,ndx-ien-2
        Phi1step(i)= ((Phiv(i) - Phidt(i)) / cg / dt + (Phiv(i+1) - Phiv(i-1)) * 0.5d0 /dx) !wind up?
        !write(201,*) sngl(x(i)) , Phi1step(i)
     !   write(201,*) i , Phi1step(i),Phiv(i) , Phidt(i),Phiv(i+1) , Phiv(i-1),Phiv(i+1) - Phiv(i-1),Phiv(i) - Phidt(i)
     end do
!     Phi1step(ist+1)=Phi1step(ist)
!     Phi1step(ndx-ien-1)=Phi1step(ndx-ien)
     !close(201)
     !write(*,*) Phi1step(ndx-ien),Phidt(ndx-ien),Phiv(ndx-ien-1),'???'
  end if
  !----------- +cg -----------

  !----------- -cg -----------
  if(mode==2) then
     !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/mcg'//name//'.dat')
     !do i=ist-1,ndx-ien+1
     do i=ist,ndx-ien
        Phi1step(i)= ((Phiv(i) - Phidt(i)) / cg / dt - (Phiv(i+1) - Phiv(i-1)) * 0.5d0 /dx)

        !Phi1step(i)= ((Phiv(i) - Phidt(i) ) / cg / dt -&
        !     (3.0d0*Phiv(i) - 4.0d0*Phiv(i-1)+Phiv(i-2)) * 0.5d0 /dx)
        !cal(i)= ((3.0d0*Phidt(i) -4.0d0* Phi2dt(i) + Phi3dt(i) )*0.5d0 / cg / dt +&
        !     (3.0d0*Phidt(i) - 4.0d0*Phidt(i+1)+Phidt(i+2)) * 0.5d0 /dx)
        !cal(i)= ((3.0d0*Phidt(i) -4.0d0* Phi2dt(i) + Phi3dt(i) )*0.5d0 / cg / dt -&
        !(3.0d0*Phidt(i) - 4.0d0*Phidt(i-1)+Phidt(i-2)) * 0.5d0 /dx)
        !write(201,*) sngl(x(i)) , Phi1step(i)
     !   write(201,*) i , Phi1step(i),Phiv(i) , Phidt(i),Phiv(i+1) , Phiv(i-1),Phiv(i+1) - Phiv(i-1),&
     !        Phiv(i) - Phidt(i),(Phiv(i+1) - Phiv(i-1)) * 0.5d0 /dx,(Phiv(i) - Phidt(i)) / cg / dt
     end do
     !close(201)
  end if
  !----------- -cg -----------
  cnt=cnt+1
  !dtpre=dt
end subroutine split

subroutine cllsub(mode,dt)
  use comvar
  use grvvar
  integer :: mode!,ifEVO=1,ifEVO2=1
  double precision dt
  if(mode==1) then
     call BC(1)
     call BC(8)
     call BC(28)
     call BC(9)
     call BC(29)
  end if

  if(mode==2) then
     call time(dt)
     call timesource(Phicgp,Phi1step,dt,2)
     call timesource(Phi1step,rho,dt,1)
     call timesource(Phicgm,Phi2step,dt,2)
     call timesource(Phi2step,rho,dt,1)
  end if

  if(mode==3) then
     call BC(4)
     call BC(3)
     call BC(18)
     call BC(38)
     call BC(19)
     call BC(39)
  end if

  if(mode==4) then
     if(ifEVO.eq.1) then
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        !ifEVO = 2
        goto 1000
     end if
     if(ifEVO.eq.2) then
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        !ifEVO = 3
        goto 1000
     end if
     if(ifEVO.eq.3) then
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        !ifEVO = 4
        goto 1000
     end if
     if(ifEVO.eq.4) then
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        !ifEVO = 5
        goto 1000
     end if
     if(ifEVO.eq.5) then
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        !ifEVO = 6
        goto 1000
     end if
     if(ifEVO.eq.6) then
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
        !ifEVO = 1
        goto 1000
     end if
1000 continue
  end if


  if(mode==5) then
     if(ifEVO2.eq.1) then
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        !ifEVO2 = 2
        goto 1100
     end if
     if(ifEVO2.eq.2) then
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        !ifEVO2 = 3
        goto 1100
     end if
     if(ifEVO2.eq.3) then
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        !ifEVO2 = 4
        goto 1100
     end if
     if(ifEVO2.eq.4) then
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        !ifEVO2 = 5
        goto 1100
     end if
     if(ifEVO2.eq.5) then
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        !ifEVO2 = 6
        goto 1100
     end if
     if(ifEVO2.eq.6) then
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp,Phi1step,dt,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2)
        !ifEVO2 = 1
        goto 1100
     end if
1100 continue
  end if
end subroutine cllsub
