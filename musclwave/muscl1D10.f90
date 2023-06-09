module comvar
  implicit none
  integer, parameter :: ndx=514,laststep=28000!,ist=1,ien=2 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  integer :: ist,ien
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  DOUBLE PRECISION :: cg = 1.0d0 , dx != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.1d0 !,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=514 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x , Phicgm ,rho, Phi1step , Phi2step ,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phidt,Phigrd,Phiexa
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2


  call INITIAL()
  call BC(1)
  !call muslcslv1D(Phi,Phi1step,dt,13)

  do i=1,laststep
     call time(dt)
     write(*,*) i ,'step'


     call timesource(Phicgp,Phi1step,dt,2)
     call timesource(Phi1step,rho,dt,1)

     call timesource(Phicgm,Phi2step,dt,2)
     call timesource(Phi2step,rho,dt,1)

     !call BC(1)
     !call cnbn(Phicgp,Phicgm)
     !call cnbn(Phi1step,Phi2step)
     !call BC(1)

     !call BC()
     !call muslcslv1D(Phi1step,rho,dt,1)
     call BC(4)
     call BC(3)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,3,1)
     !call muslcslv1D(Phi1step,rho,0.5d0*dt,3)
     call BC(4)
     call BC(3)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,1)

     call BC(4)
     call BC(3)
     call BC(1)
     !Phidt(:)=Phi(:)
     call muslcslv1D(Phicgp,Phi1step,dt,4,1)
     call muslcslv1D(Phicgm,Phi2step,dt,4,2)
     call BC(1)
     call muslcslv1D(Phicgp,Phi1step,dt,1,1)
     call muslcslv1D(Phicgm,Phi2step,dt,2,2)
     !call BC(1)

     call BC(4)
     call BC(3)
     !call muslcslv1D(Phi1step,rho,dt,3)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
     call muslcslv1D(Phi2step,rho,0.5d0*dt,3,1)
     call BC(4)
     call BC(3)
     !call muslcslv1D(Phi1step,rho,dt,1)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
     !call BC(3)
     !call BC(4)
     !call BC(1)


     !call saveu(sv)
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
  integer :: i
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


  !---------Phi-------------
  Phicgp(:)=0.0d0
  Phicgm(:)=0.0d0
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:)=0.0d0
  Phi2step(:)=0.0d0
  !-------Phi1step-----------

  !-------Phidt-----------
  Phidt(:)=0.0d0
  !-------Phdt-----------




  !---------rho-------------
  do i = -1,ndx
     if( dabs(x(i) - hcen) .le. h) then
        rho(i) = dinit1
        !rho(i) = 0.0d0
     else
        rho(i) = 0.0d0
        !rho(i) = dinit1*1.d-2
     end if
  end do
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  open(142,file='/Users/maeda/Desktop/kaiseki/testcode6/phiexact.DAT')
  open(143,file='/Users/maeda/Desktop/kaiseki/testcode6/INIden.DAT')
  open(144,file='/Users/maeda/Desktop/kaiseki/testcode6/phigrd.DAT')
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
     else
        Phiexa(i) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
     end if
     write(143,*) sngl(rho(i))
  end do


  do i=0,ndx-1
     Phigrd(i)=(-Phiexa(i-1)+Phiexa(i+1))*0.5d0/dx
     !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
  end do
  Phigrd(-1)=(-Phiexa(0)+Phiexa(1))/dx
  Phigrd(ndx)=(Phiexa(ndx-1)-Phiexa(ndx-2))/dx

  do i=-1,ndx
     write(144,*) sngl(x(i)) , Phigrd(i) !, Phiexa(i-1),Phiexa(i+1)
  end do

  bcphi1(1) = G4pi * dinit1 * h * dabs(x(1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
  bcphi2(1) = G4pi * dinit1 * h * dabs(x(ndx-2) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

  bcphi1(2) = G4pi * dinit1 * h * dabs(x(0) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
  bcphi2(2) = G4pi * dinit1 * h * dabs(x(ndx-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2

  bcphi1(3) = G4pi * dinit1 * h * dabs(x(-1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
  bcphi2(3) = G4pi * dinit1 * h * dabs(x(ndx) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
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
     Phi1step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
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
  Phicgp(:)=bcphi1(1)
  Phicgm(:)=bcphi1(1)
  !---------Phi-------------

  !-------Phidt-----------
  Phidt(:)=bcphi1(1)
  !-------Phdt-----------
  !-------Phi1step-----------
  !Phi1step(:)=bcphi1
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
     Phicgp(0)= bcphi1(2)
     Phicgp(-1)= bcphi1(3)
     !Phicgp(ndx-2)= bcphi2(1)
     !Phicgp(ndx-1)= bcphi2(2)
     !Phicgp(ndx)= bcphi2(3)
     Phicgp(ndx)=Phicgp(ndx-1)
     Phicgp(ndx-1)= Phicgp(ndx-2)
     !Phicgp(ndx-2)= Phicgp(ndx-3)

     !Phicgm(1)= bcphi1(1)
     !Phicgm(0)= bcphi1(2)
     !Phicgm(-1)= bcphi1(3)
     Phicgm(-1)= Phicgm(0)
     Phicgm(0)= Phicgm(1)
     !Phicgm(1)= Phicgm(2)
     !Phicgm(ndx-2)= bcphi2(1)
     Phicgm(ndx-1)= bcphi2(2)
     Phicgm(ndx)= bcphi2(3)


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
     Phi1step(-1)= Phi1step(0)
     Phi1step(0)= Phi1step(1)
     !Phi1step(1)=Phi1step(2)
     !Phi1step(ndx-2)= Phigrd(ndx-2)
     Phi1step(ndx-1)= Phigrd(ndx-1)
     Phi1step(ndx)= Phigrd(ndx)
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
     !Phi1step(ndx)= -Phigrd(ndx-2)
     !Phi2step(1)= -Phigrd(1)
     Phi2step(0)= -Phigrd(2)
     Phi2step(-1)=-Phigrd(3)
     !Phi2step(ndx-2)= -Phigrd(ndx-2)
     !Phi2step(ndx-1)= -Phigrd(ndx-1)
     !Phi2step(ndx)= -Phigrd(ndx)
     Phi2step(ndx)= -Phi2step(ndx-1)
     Phi2step(ndx-1)= -Phi2step(ndx-2)
     !Phi2step(ndx-2)= -Phi2step(ndx-3)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
     !701 continue
     !-------Phi1step-----------
  end if

  if(mode==5) then
     !--------free--------------
     goto 112
     !-------Phi1step-----------
     Phi1step(1)= Phi1step(2)
     Phi1step(0)= Phi1step(1)
     Phi1step(-1)=Phi1step(0)
     Phi1step(ndx-2)= Phi1step(ndx-3)
     Phi1step(ndx-1)= Phi1step(ndx-2)
     Phi1step(ndx)= Phi1step(ndx-1)
     !-------Phi1step-----------
112  continue
     !--------free--------------
  end if


  if(mode==6) then
     !---------kotei2-----------
     goto 105
     !---------Phi-------------
!     Phi(0)= Phi(1)
!     Phi(-1)= Phi(-1)
!     Phi(ndx-1)=Phi(ndx-2)
!     Phi(ndx)= Phi(ndx-1)
     !---------Phi-------------

     !-------Phi1step-----------
     Phi1step(1)= bcphi1(1)
     Phi1step(0)= bcphi1(2)
     Phi1step(-1)= bcphi1(3)
     Phi1step(ndx-2)= bcphi2(1)
     Phi1step(ndx-1)= bcphi2(2)
     Phi1step(ndx)= bcphi2(3)
     !-------Phi1step-----------
105  continue
     !---------kotei2-----------
  end if


  if(mode==7) then
     !---------perio-----------
     goto 101
     !---------Phi-------------
!     Phi(0)= Phi(ndx-2)
!     Phi(-1)= Phi(ndx-3)
!     Phi(ndx-1)= Phi(1)
!     Phi(ndx)= Phi(2)
     !---------Phi-------------

     !-------Phi1step-----------
     Phi1step(0)= Phi1step(ndx-2)
     Phi1step(-1)= Phi1step(ndx-3)
     Phi1step(ndx-1)= Phi1step(1)
     Phi1step(ndx)= Phi1step(2)
     !-------Phi1step-----------
101  continue
     !---------perio-----------
  end if



  if(mode==8) then
     !-------period2-----------
     goto 102
     !---------Phi-------------
!     pr(2)= Phi(ndx-2)
!     pr(1)= Phi(ndx-3)
!     pl(1)= Phi(1)
!     pl(2)= Phi(2)
!     Phi(1)=pr(1)
!     Phi(2)=pr(2)
!     Phi(ndx-2)=pl(1)
!     Phi(ndx-3)=pr(2)
     !---------Phi-------------
     !-------Phi1step-----------
!     pr(2)= Phi1step(ndx-2)
!     pr(1)= Phi1step(ndx-3)
!     pl(1)= Phi1step(1)
!     pl(2)= Phi1step(2)
!     Phi1step(1)=pr(1)
!     Phi1step(2)=pr(2)
!     Phi1step(ndx-2)=pl(1)
!     Phi1step(ndx-3)=pr(2)
     !-------Phi1step-----------
102  continue
     !-------period2-----------
  end if





  if(mode==9) then
     !-------katagawa-----
     !goto 130
     !---------Phi-------------
!     Phi(1)= bcphi1(1)
!     Phi(0)= bcphi1(2)
!     Phi(-1)= bcphi1(3)
!     Phi(ndx-2)= Phi(ndx-3)
!     Phi(ndx-1)= Phi(ndx-2)
!     Phi(ndx)= Phi(ndx-1)
     !---------Phi-------------
  end if


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
     Phi1step(1)= Phi1step(2)
     Phi1step(0)= Phi1step(1)
     Phi1step(-1)=Phi1step(0)
     Phi1step(ndx-2)= bcphi2(1)
     Phi1step(ndx-1)= bcphi2(2)
     Phi1step(ndx)= bcphi2(3)
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
  integer i,mode
  double precision :: dt,sdt,mindt,maxdt , epsl = 1.0d-4
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source

  !mindt=1000.0d0
  maxdt=0.0d0

  if(mode==1) then
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           sdt = 0.5d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
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
           sdt = 0.5d0*dabs(Phiv(i)) / ( cg * source(i) )
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


subroutine muslcslv1D(Phiv,source,dt,mode,bcmode)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0,bcmode
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name
  if(bcmode==1)then
     ist=2
     ien=2
  end if
  if(bcmode==2)then
     ist=1
     ien=3
  end if

  nu2 = cg * dt / dx
  Phipre(:) = Phiv(:)
  !write(name,'(i5.5)') cnt



  !------------ul.solver.+cg-------------
  if(mode==1) then
     !write(name,'(i5.5)') cnt
     !open(201,file='/Users/maeda/Desktop/kaiseki/testcode2/1modepre'//name//'.dat')
     !write(name,'(i5.5)') cnt+1
     !open(202,file='/Users/maeda/Desktop/kaiseki/testcode2/1modepos'//name//'.dat')
     !do i=-1,ndx
     !   write(201,*) i, Phiv(i)
     !end do
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     do i=ist-1,ndx-ien+1 !一次なので大丈夫
        Phi2dt(i) = Phipre(i) - 0.5d0 * nu2 * ( Phiu(i) - Phiu(i-1))
     end do
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,-1.d0,1)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     do i = ist , ndx-ien
        Phiv(i) = Phipre(i) - nu2 * (Phiu(i) - Phiu(i-1))
     end do
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
     do i=ist-1,ndx-ien+1
        Phi2dt(i) = Phipre(i) + 0.5d0 * nu2 * ( Phiu(i+1) - Phiu(i))
     end do
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,-1.0d0,4)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     do i = ist , ndx-ien
        Phiv(i) = Phipre(i) + nu2 * (Phiu(i+1) - Phiu(i))
     end do

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
     do i=ist,ndx-ien
        Phiv(i) =  -cg * G4pi * source(i) * dt + Phipre(i)
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
     do i=ist,ndx-ien
        Phiv(i) = cg * source(i) * dt + Phipre(i)
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
subroutine vanalbada(Phipre,Phigrad)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  !integer :: i , ip , im , flmt ,eps=1.0d-10
  integer :: i , ip , im
  DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre

  do i = ist-1 , ndx - ien + 1
     ip=i+1
     im=i-1

     delp = Phipre(ip)-Phipre(i)
     delm = Phipre(i)-Phipre(im)
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !flmt = (2.d0*delp*delm+eps)/(delp**2+delm**2+eps)
     !if(i==58) then
     !   write(*,*) delp , delm ,Phipre(ip),Phipre(i),Phipre(im) , flmt
     !end if
     !flmt = (2.d0*delp*delm+eps)/(delp**2+delm**2+eps)
     Phigrad(i) = flmt
     !grdU(i,k) = flmt*( wave(ixp,jyp,kzp)-wave(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )
     !Phigrad(i) = flmt*( Phipre(ip)-Phipre(im) )/( 2.0d0 * dx )
  end do
end subroutine vanalbada

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file='/Users/maeda/Desktop/kaiseki/testcode6/phi'//name//'.dat')
  do i=1,ndx-2
     if(i==1 .or. i==ndx-2) then
        write(21,*) x(i), Phicgp(i),Phi1step(i) , Phicgm(i),Phi2step(i) ,&
             (Phicgp(i)+Phicgm(i))*0.5d0,(Phi1step(i)+Phi2step(i))*0.5d0, (Phicgp(i)-Phicgm(i)),&
             (Phicgp(i)+Phicgm(i))*0.5d0+dabs(Phicgp(i)-Phicgm(i)),&
                                !-((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx,&
                                !(Phicgp(i)+Phicgm(i))*0.5d0-Phiexa(i),&
                                !-((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx-Phigrd(i),&
                                !(-((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx-Phigrd(i))/&
                                !((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx,&
                                !-((Phicgp(i)+Phicgp(i))*0.5d0 - (Phicgp(i)+Phicgp(i))*0.5d0)*0.5d0/dx,&
                                !-((Phicgm(i)+Phicgm(i))*0.5d0 - (Phicgm(i)+Phicgm(i))*0.5d0)*0.5d0/dx
                                !Phigrd(i),Phigrd(i)
             Phi1step(i),-Phi2step(i)
        write(*,*) Phi1step(i),-Phi2step(i),i,Phigrd(i)
     else
        write(21,*) x(i), Phicgp(i),Phi1step(i) , Phicgm(i),Phi2step(i) ,&
             (Phicgp(i)+Phicgm(i))*0.5d0,(Phi1step(i)+Phi2step(i))*0.5d0, (Phicgp(i)-Phicgm(i)),&
             (Phicgp(i)+Phicgm(i))*0.5d0+dabs(Phicgp(i)-Phicgm(i)),&
                                !-((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx,&
                                !(Phicgp(i)+Phicgm(i))*0.5d0-Phiexa(i),&
                                !-((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx-Phigrd(i),&
                                !(-((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx-Phigrd(i))/&
                                !((Phicgp(i-1)+Phicgm(i-1))*0.5d0 - (Phicgp(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx,&
             -((Phicgp(i-1)+Phicgp(i-1))*0.5d0 - (Phicgp(i+1)+Phicgp(i+1))*0.5d0)*0.5d0/dx,&
             -((Phicgm(i-1)+Phicgm(i-1))*0.5d0 - (Phicgm(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx
     end if
  end do
  close(21)
   !i=1
   !write(*,*) Phi1step(i),-Phi2step(i),Phigrd(i)
  !i=1
  !write(*,*)-((Phicgp(i-1)+Phicgp(i-1))*0.5d0 - (Phicgp(i+1)+Phicgp(i+1))*0.5d0)*0.5d0/dx,&
  !     -((Phicgm(i-1)+Phicgm(i-1))*0.5d0 - (Phicgm(i+1)+Phicgm(i+1))*0.5d0)*0.5d0/dx,&
  !     Phigrd(i)
  in1=in1+1
end subroutine saveu



subroutine fluxcal(preuse,pre,u,ep,kappa,mode)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx) :: ul,ur,pre,slop,preuse,u
  integer :: i,mode
  !u(:)=0.0d0
  call vanalbada(pre,slop)
  if(mode==1) then
     do i = ist-1,ndx-ien+1
        ul(i) = preuse(i) + 0.25d0 * ep * slop(i) &
             * ((1.0d0-slop(i)*kappa)*(pre(i)-pre(i-1)) + (1.0d0+slop(i)*kappa)*(pre(i+1) - pre(i))) !i+1/2
        u(i)=ul(i)
     end do
     write(*,*) slop(127),'127slop'
     !u(:)=ul(:)
  end if


  if(mode==4) then
     do i = ist-1,ndx-ien+1
        ur(i) = preuse(i) - 0.25d0 * ep * slop(i) &
             * ((1.0d0+slop(i)*kappa)*(pre(i)-pre(i-1)) + (1.0d0-slop(i)*kappa)*(pre(i+1) - pre(i))) !i-1/2
        u(i)=ur(i)
     end do
     !write(*,*) slop(127),'127slop'
     !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
     !write(*,*) u(2)
     !u(:)=ur(:)
  end if

  if(mode==10) then
     do i = ist-2,ndx-ien+2
        ul(i) = preuse(i)
        u(i)=ul(i)
     end do
  end if

  if(mode==11) then
     do i = ist-2,ndx-ien+2
        ur(i) = preuse(i)
        u(i)=ur(i)
     end do
  end if

  if(mode==12) then
     do i = ist-1,ndx-ien+1
        ur(i) = preuse(i+1)
        u(i)=ur(i)
     end do
  end if



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

