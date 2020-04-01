module comvar
  implicit none
  integer, parameter :: ndx=130,laststep=1000,ist=1,ien=2,svnum=1 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  DOUBLE PRECISION :: cg = 1.0d0 , dx != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0 ,meanrho!, kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  DOUBLE PRECISION , dimension(-1:ndx) :: Phidt1,Phidt2
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=130 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phi,v!,w1,w2
  DOUBLE PRECISION , dimension(-1:ndx2) :: x , Phicgm ,rho, Phi1step , Phi2step ,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phidt,Phigrd,Phiexa,wp1,wp2,xi1,xi2
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2,j
   DOUBLE PRECISION , dimension(-1:ndx2) :: Phidummy,Phidummy1,Phidummy2!,wp1,wp2
   DOUBLE PRECISION , dimension(-1:ndx2) :: xi1dummy,xi2dummy

   call INITIAL()
   !Phidt(:)=Phi(:)
   call BC(1)
   Phidt(:)=Phi(:)
   Phidt2(:)=Phi(:)

   do i=1,laststep
      !if(mod(i,svnum)==1) then
      !  call saveu(sv)
     !end if
      call time(dt)
     ! call timesource2(Phi,Phidt,rho,dt)
     write(*,*) i,dt,'step'

     Phidummy(:)=Phi(:)
     Phidummy1(:)=Phi(:)
     Phidummy2(:)=Phi(:)
     Phidt1(:)=Phidt(:)
     xi1dummy(:)=xi1(:)
     xi2dummy(:)=xi2(:)
    !  write(*,*)'ok2'
!     call muslcslv1D(Phi,rho,dt,3)
     !call BC(1)
     !call split(Phi,Phicgm,Phicgp,Phidt,dt)
     !call BC(1)
     ! write(*,*)'ok3'
     do j=-1,ndx
        wp1(j)=Phi(j)-cg*v(j)
        wp2(j)=Phi(j)+cg*v(j)
     end do
     call muslcslv1D(wp1,rho,dt,2)
     call muslcslv1D(wp2,rho,dt,1)
     do j=1,ndx-2
        wp1(j) = wp1(j) -  dt * xi2(j)
        wp2(j) = wp2(j) -  dt * xi1(j)
     end do
      do j=1,ndx-2
         Phi(j)=0.5d0*(wp1(j)+wp2(j))
         v(j)=0.5d0*(-wp1(j)+wp2(j))/cg
      end do
     call BC(1)
     !do j=1,ndx-2
     !   Phi(j) = 2.d0 * Phidt(j) - Phidt2(j) - G4pi * cg * cg * dt * dt * rho(j)
     !end do
     !call BC(1)
     Phidt2(:)=Phidt(:)
     Phidt(:)=Phidummy(:)
     !Phidt(:)=Phi(:)

     call muslcslv1D(xi1,rho,dt,1)
     call muslcslv1D(xi2,rho,dt,2)
     call BC(2)
     do j=1,ndx-2
        xi1(j) = xi1dummy(j) - G4pi * cg * cg * dt * rho(j)
        xi2(j) = xi2dummy(j) - G4pi * cg * cg * dt * rho(j)
     end do
     call BC(2)

     if(mod(i,svnum)==0) then
        call saveu(sv)
     end if
  end do
 ! call BC(3)
 ! call BC(4)
  call BC(1)
  call saveu(sv)
end program muscl1D

subroutine INITIAL()
  use comvar
  use grvvar
  integer :: i
  double precision :: amp,pi=3.1415926535d0,haba

  dinit1 = 2.0d0/G4pi/90.d0*1.d-10

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
  Phi(:) = 0.d0
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:)=0.0d0
  Phi2step(:)=0.0d0
  !Phi1step(:)=+G4pi*meanrho*cg*Lbox
  !Phi2step(:)=0.0d0
  !-------Phi1step-----------

  !-------Phidt-----------
  Phidt(:)=0.0d0
  !-------Phdt-----------

  !----------v--------------
  v(:)=0.d0
  !----------v--------------

  !----------xi--------------
  xi1(:)=0.d0
  xi2(:)=0.d0
  !----------xi--------------



  !---------rho-------------
  do i = -1,ndx
     if( dabs(x(i) - hcen) .le. h) then
        rho(i) = dinit1
        !rho(i) = 0.0d0
     else
        rho(i) = 0.0d0
        !rho(i) = dinit1
        !rho(i) = dinit1*1.d-2
     end if
  end do

  meanrho=0.d0
  do i = 1,ndx-2
     meanrho=meanrho+rho(i)
  end do
  meanrho=meanrho/dble(ndx-2)

  do i = -1,ndx
     rho(i)=rho(i)-meanrho
  end do


  Phi1step(:)=0.d0 !+G4pi*meanrho*cg*Lbox
  Phi2step(:)=0.d0 !+G4pi*meanrho*cg*Lbox
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  open(142,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscledf/phiexact.DAT')
  open(143,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscledf/INIden.DAT')
  open(144,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscledf/phigrd.DAT')
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

  !do i=0,ndx-1
  !   Phigrd(i)=-(-Phiexa(i-1)+Phiexa(i+1))*0.5d0/dx
     !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
  !end do
  !Phigrd(-1)=-(-Phiexa(0)+Phiexa(1))/dx
  !Phigrd(ndx)=-(Phiexa(ndx-1)-Phiexa(ndx-2))/dx

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
  !goto 201
  !do i = -1, ndx
  !   amp = 1.d-3
  !   Phi(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  !   Phi1step(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  !end do


!  do i = -1, ndx
!     amp = 1.d-3
!     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!     Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!     v(i)= cg*amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
!  end do
  !201 continue
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

  Phi(   0)=Phi(ndx-2)
  Phi(  -1)=Phi(ndx-3)
  Phi(ndx-1)=Phi(   1)
  Phi(ndx  )=Phi(   2)

  v(   0)=v(ndx-2)
  v(  -1)=v(ndx-3)
  v(ndx-1)=v(   1)
  v(ndx  )=v(   2)


  xi1(   0)=xi1(ndx-2)
  xi1(  -1)=xi1(ndx-3)
  xi1(ndx-1)=xi1(   1)
  xi1(ndx  )=xi1(   2)

  xi2(   0)=xi2(ndx-2)
  xi2(  -1)=xi2(ndx-3)
  xi2(ndx-1)=xi2(   1)
  xi2(ndx  )=xi2(   2)

!  Phicgp(   0)=Phicgp(ndx-2)
!  Phicgp(  -1)=Phicgp(ndx-3)
!  Phicgp(ndx-1)=Phicgp(   1)
!  Phicgp(ndx  )=Phicgp(   2)

!  Phicgm(   0)=Phicgm(ndx-2)
!  Phicgm(  -1)=Phicgm(ndx-3)
!  Phicgm(ndx-1)=Phicgm(   1)
!  Phicgm(ndx  )=Phicgm(   2)

  !Phi(   0)=Phiexa(   0)
  !Phi(  -1)=Phiexa(  -1)
  !Phi(nx-1)=Phiexa(nx-1)
  !Phi(nx  )=Phiexa(nx  )
end subroutine BC


subroutine time(dt)
  use comvar
  use grvvar
  double precision :: dt
  dt = dx/cg * coeff
  write(*,*) 'time cg' , dt
end subroutine time



subroutine muslcslv1D(Phiv,source,dt,mode)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name

  nu2 = cg * dt / dx
  Phipre(:) = Phiv(:)
  !write(name,'(i5.5)') cnt



  !------------ul.solver.+cg-------------
  if(mode==1) then
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
        write(*,*) nu2 * (Phiu(i) - Phiu(i-1))
     end do
     !write(*,*) Phiv(127),'127-3'
     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do
  end if
  !------------ul.solver.+cg-------------



  !------------ul.solver.-cg-------------
  if(mode==2) then

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
     write(*,*) 'mode1'

  end if
  !------------ul.solver.-cg-------------


  !--------------source------------------
  if(mode==3) then
     do i=ist,ndx-ien
        Phiv(i) =  -(cg**2 )* G4pi * source(i) * (dt**2) + 2.d0*Phipre(i)-Phidt1(i)
     end do

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do
  end if

  if(mode==4) then
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
  open(21,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscledf/phi'//name//'.dat')
  do i=1,ndx-2
        write(21,*) x(i),Phi(i), v(i),wp1(i),wp2(i),xi1(i),xi2(i)!Phicgp(i) , Phicgm(i) ,rho(i)!,&
  end do
  close(21)
  write(*,*) 'save step : ',in1
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


subroutine split(Phi,Phicgm,Phicgp,Phidt,dt)
  use comvar
  integer :: i,mode,cnt=0
  double precision :: dt,dtpre,itgphi,itgphi2,mean  !前の
  DOUBLE PRECISION, dimension(-1:ndx) :: Phi,Phicgp,Phicgm,Phidt,Phitgrd
  character(5) name
  write(*,*) dx

  ! write(*,*)'ok1'
  do i = -1,ndx
     Phitgrd(i)=(Phi(i)-Phidt(i))/dt
  end do

  mean=0.d0
   do i=1,ndx-2
      mean = Phitgrd(i)*dx/cg + mean
   end do
   mean = mean/dble(ndx-2)

  ! write(*,*)'ok14'
  !do i=1,ndx-2
  !   itgphi = Phitgrd(i)*dx + itgphi
  !end do
  !itgphi=itgphi/cg
  !itgphi=0.d0
  itgphi = -Phitgrd(-1)*dx/cg-mean! + itgphi
  !write(*,*) itgphi
!  itgphi2= -Phitgrd(ndx)*dx/cg! + itgphi
  Phicgm(-1)=0.5d0*(Phi(-1) - itgphi )
  Phicgp(-1)=0.5d0*(Phi(-1) + itgphi )
  !write(*,*)'ok13'
  Phicgm(0)=0.5d0*Phi(0)
  Phicgp(0)=0.5d0*Phi(0)
  write(*,*) Phicgm(-1), Phicgm(0),Phicgp(-1), Phicgp(0),mean

  itgphi=0.d0
  do i=1,ndx
     itgphi = Phitgrd(i)*dx/cg + itgphi-mean
     !write(*,*)itgphi,Phitgrd(i),Phi(i),Phidt(i)
     Phicgm(i)=0.5d0*(Phi(i) - itgphi )
     Phicgp(i)=0.5d0*(Phi(i) + itgphi )
     !write(*,*) Phitgrd(1),dx,cg,itgphi,Phicgm(1),Phicgp(1)
  end do

!  itgphi2=0.d0
!  do i=-1,ndx-2,-1
!     itgphi2 = Phitgrd(i)*dx/cg + itgphi2
     !Phicgm(i)=0.5d0*(Phi(i) - itgphi )
!     Phicgp(i)=0.5d0*(Phi(i) + itgphi2 )
     !write(*,*) Phitgrd(i),dx,cg,itgphi,Phicgm(i),Phicgp(i)
!  end do
  ! write(*,*)'ok11'
end subroutine split

subroutine timesource2(Phiv,Phidt,source,dt)
  use comvar
  !use grvver
  integer i,mode
  double precision :: dt,sdt1,mindt,maxdt,std2,rt=0.5d0
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source,Phidt
  maxdt=0.d0
  mindt=1.d10
  do i=1,ndx-2
     if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
        !sdt = 0.5d0*dsqrt( dabs(2*Phiv(i)-Phidt(i)) / (cg * cg * G4pi * dabs(source(i) )) )

        sdt = dsqrt(dabs( (-Phiv(i)+rt*Phiv(i)+Phidt(i))/(-2.d0*cg*cg*G4pi*source(i)) ))
        sdt2 = dsqrt(dabs( (-Phiv(i)-rt*Phiv(i)+Phidt(i))/(-2.d0*cg*cg*G4pi*source(i)) ))
        !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
        mindt=dmin1(mindt,sdt)
        mindt=dmin1(mindt,sdt2)
        !maxdt=dmax1(maxdt,sdt)
     end if
  end do
  !if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
  !   dt = sdt
  !end if
  if( (mindt < dt) .and. (mindt .ne. 0.0d0)) then
     dt = sdt
  end if


  write(*,*) 'time source2' , dt
end subroutine timesource2

