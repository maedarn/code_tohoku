module comvar
  implicit none
  integer :: ndx=130,laststep=1000
  DOUBLE PRECISION :: cg = 1.0d0 , dx , bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.1d0 ,  kappa=1.0d0/3.0d0
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=130 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x,Phi,rho,Phiexa, Phi1step
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i


  call INITIAL()
  call BC()

  do i=1,laststep
     call time(dt)
     !call BC()
     !call timesource(Phi,Phi1step,dt,2)
     !call timesource(Phi1step,rho,dt,1)
     !call BC()
     !------source-------
     !call  muslcslv1D(Phi,Phi1step,0.5d0*dt,4)
     !call  muslcslv1D(Phi1step,rho,0.5d0*dt,3)
     !------source-------

     call BC()
     !call  muslcslv1D(Phi,Phi1step,dt,5)
     call  muslcslv1D(Phi1step,rho,dt,6)


     !call BC()
     !------source-------
     !call  muslcslv1D(Phi,Phi1step,0.5d0*dt,4)
     !call  muslcslv1D(Phi1step,rho,0.5d0*dt,3)
     !------source-------

     call saveu(i)
  end do

end program muscl1D

subroutine INITIAL()
  use comvar
  use grvvar
  integer :: i
  double precision :: amp,pi=3.1415926535d0

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
  Phi(:)=0.0d0
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:)=0.0d0
  !-------Phi1step-----------




  !---------rho-------------
  do i = -1,ndx
     if( dabs(x(i) - hcen) .le. h) then
        !rho(i) = dinit1
        rho(i) = 0.0d0
     else
        rho(i) = 0.0d0
     end if
  end do
  !---------rho-------------



  !--------Phiexa-----------
  goto 200
  open(142,file='/Users/maeda/Desktop/kaiseki/testcode1/phiexact.DAT')
  open(143,file='/Users/maeda/Desktop/kaiseki/testcode1/INIden.DAT')
  do i= 1,ndx-2
     if( dabs(x(i) - hcen) .le. h ) then
        !Phiexa(i,j,k) = G4pi/2.0d0 * dinit1 * (x(i) - censh )**2
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
     else
        !Phiexa(i,j,k) = G4pi * dinit1 * Hsheet * dabs(x(i) - censh)  - G4pi/2.0d0 * dinit1 * Hsheet**2
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
     end if
     write(143,*) sngl(rho(i))
  end do
  bcphi1 = G4pi * dinit1 * h * dabs(x(1) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
  bcphi2 = G4pi * dinit1 * h * dabs(x(ndx-2) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
  close(142)
  close(143)
  200 continue
  !--------Phiexa-----------


  !---------wave--------
  !goto 201
  do i = -1, ndx
     amp = 1.d-3
     Phi(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
     Phi1step(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
  end do
  !201 continue
  !---------wave--------
end subroutine INITIAL



subroutine BC()
  use comvar
  use grvvar
  integer :: i
  !double precision , dimension(1:2) :: pl,pr

  !---------kotei-----------
  goto 100
  !---------Phi-------------
  Phi(1)= bcphi1
  Phi(0)= bcphi1
  Phi(-1)= bcphi1
  Phi(ndx-2)= bcphi2
  Phi(ndx-1)= bcphi2
  Phi(ndx)= bcphi2
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(1)= bcphi1
  Phi1step(0)= bcphi1
  Phi1step(-1)= bcphi1
  Phi1step(ndx-2)= bcphi2
  Phi1step(ndx-1)= bcphi2
  Phi1step(ndx)= bcphi2
  !-------Phi1step-----------
  100 continue
  !---------kotei-----------


  !---------perio-----------
  ! goto 101
  !---------Phi-------------
  Phi(0)= Phi(ndx-2)
  Phi(-1)= Phi(ndx-3)
  Phi(ndx-1)= Phi(1)
  Phi(ndx)= Phi(2)
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(0)= Phi1step(ndx-2)
  Phi1step(-1)= Phi1step(ndx-3)
  Phi1step(ndx-1)= Phi1step(1)
  Phi1step(ndx)= Phi1step(2)
  !-------Phi1step-----------
  !101 continue
  !---------perio-----------
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
  double precision :: dt,sdt,mindt
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source

  mindt=1000.0d0

  if(mode==1) then
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           !sdt = 0.5d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
           sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
           mindt=dmin1(mindt,sdt)
        end if
     end do
     if(sdt < dt) then
        dt = sdt
     end if
  end if


  if(mode==2) then
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           !sdt = 0.5d0*dabs(Phiv(i)) / ( source(i) )
           sdt = 0.2d0*dabs(Phiv(i)) / ( source(i) )
           mindt=dmin1(mindt,sdt)
        end if
     end do
     if(sdt < dt) then
        dt = sdt
     end if
  end if

  write(*,*) 'time source' , dt
end subroutine timesource



subroutine muslcslv1D(Phiv,source,dt,mode)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam !kappa -> comver  better?
  integer direction , mode , invdt , loopmode , dloop
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt



  Phipre(:) = Phiv(:)

  !------------calcurate dt/2------------
  Phi2dt(i) = Phipre(i) - 0.5d0 * dt / dx  ( )
  !------------calcurate dt/2------------

  !-------------MUSCL solver-------------
  nu2 = cg * dt / dx
  call vanalbada(Phipre,Phigrad)

  !------2nd step-----
  !ul(i,j,k) = u(i,j,k) + 0.25d0 * s * ((1-kappa*s)*gradum + (1+kappa*s)*gradup) !j-1
  !ur(i,j,k) = u(i,j,k) - 0.25d0 * s * ((1-kappa*s)*gradup + (1+kappa*s)*gradum) !j+1
  if(mode==1) then
     do i = 1 , ndx-2
        deltap = Phipre(i+2) - Phipre(i+1)
        deltam = Phipre(i+1) - Phipre(i  )
        fluxphi(i) = Phiv(i+1) - Phigrad(i) * 0.25d0 *( (1.0d0 - kappa * Phigrad(i)) * deltap &
             + (1.0d0 + kappa * Phigrad(i)) * deltam)  !ur_{j+1/2}
        !fluxf(i,j,k) = - cg * fluxf(i,j,k)
        !fluxf = 0.5d0 *
        !f(i,j,k) = f(i,j,k) - nu2 * 0.5d0 * (fluxf(i,j,k) - fluxf(i-1,j,k)) !dt/2 timestep
        !fdt2(i,j,k) = f(i,j,k) - nu2 * 0.5d0 * (fluxf(i,j,k) - fluxf(i-1,j,k)) !dt/2 timestep
     end do

     do i = 2 , ndx-3
        Phiv(i) = Phipre(i) - nu2 * (fluxphi(i+1) - fluxphi(i  ))
        !g(i,j,k) = gpre(i,j,k) - nu2 * (g(i  ,j,k) - g(i-1,j,k))
     end do
  end if

  if (mode==2) then
     do i = 1 , ndx-2
        deltap = Phipre(i+1) - Phipre(i  )
        deltam = Phipre(i  ) - Phipre(i-1)
        fluxphi(i) = Phiv(i  ) + Phigrad(i) * 0.25d0 *( (1.0d0 - kappa * Phigrad(i)) * deltam &
             + (1.0d0 + kappa * Phigrad(i)) * deltap)  !ul_{j+1/2}
        !fluxg(i,j,k) =   cg * fluxg(i,j,k)
     end do
     do i = 2 , ndx-3
        !Phiv(i,j,k) = Phipre(i,j,k) - nu2 * (fluxphi(i+1,j,k) - fluxphi(i  ,j,k))
        Phiv(i) = Phipre(i) - nu2 * (fluxphi(i  ) - fluxphi(i-1))
     end do
  end if

  if(mode==3) then
     write(*,*) 'in1'
     do i=2,ndx-3
        Phiv(i) =  -cg * G4pi * source(i) * dt + Phipre(i)
     end do
  end if

  if(mode==4) then
     do i=2,ndx-3
        Phiv(i) =  source(i) * dt + Phipre(i)
     end do
  end if






  !----------------periodic-----------------
  if(mode==5) then
     do i = 0 , ndx-1
        deltap = Phipre(i+2) - Phipre(i+1)
        deltam = Phipre(i+1) - Phipre(i  )
        fluxphi(i) = Phiv(i+1) - Phigrad(i) * 0.25d0 *( (1.0d0 - kappa * Phigrad(i)) * deltap &
             + (1.0d0 + kappa * Phigrad(i)) * deltam)  !ur_{j+1/2}
        !fluxf(i,j,k) = - cg * fluxf(i,j,k)
        !fluxf = 0.5d0 *
        !f(i,j,k) = f(i,j,k) - nu2 * 0.5d0 * (fluxf(i,j,k) - fluxf(i-1,j,k)) !dt/2 timestep
        !fdt2(i,j,k) = f(i,j,k) - nu2 * 0.5d0 * (fluxf(i,j,k) - fluxf(i-1,j,k)) !dt/2 timestep
     end do

     do i = 1 , ndx-2
        Phiv(i) = Phipre(i) - nu2 * (fluxphi(i+1) - fluxphi(i  ))
        !g(i,j,k) = gpre(i,j,k) - nu2 * (g(i  ,j,k) - g(i-1,j,k))
     end do
  end if

  if (mode==6) then
     do i = 0 , ndx-1
        deltap = Phipre(i+1) - Phipre(i  )
        deltam = Phipre(i  ) - Phipre(i-1)
        fluxphi(i) = Phiv(i  ) + Phigrad(i) * 0.25d0 *( (1.0d0 - kappa * Phigrad(i)) * deltam &
             + (1.0d0 + kappa * Phigrad(i)) * deltap)  !ul_{j+1/2}
        !fluxg(i,j,k) =   cg * fluxg(i,j,k)
     end do
     do i = 1 , ndx-2
        !Phiv(i,j,k) = Phipre(i,j,k) - nu2 * (fluxphi(i+1,j,k) - fluxphi(i  ,j,k))
        Phiv(i) = Phipre(i) - nu2 * (fluxphi(i  ) - fluxphi(i-1))
     end do
  end if
  !----------------periodic-----------------

  !deltap = g(i+1,j,k) - g(i  ,j,k)
  !deltam = g(i  ,j,k) - g(i-1,j,k)
  !fluxg(i,j,k) = g(i  ,j,k) + gradg * 0.25d0 *( (1.0d0 - kappa * gradg) * deltam + (1.0d0 + kappa * gradg) * deltap)  !ul_{j+1/2}
  !fluxg(i,j,k) =   cg * fluxg(i,j,k)


  !fluxr = 0.5d0 * cg * (ul(i,j,k)+ur(i,j,k) + ul(i,j,k) - ur(i,j,k)) !向きによって片方で良い
  !fluxl = 0.5d0 * cg * (ul(i,j,k)+ur(i,j,k) - ul(i,j,k) + ur(i,j,k)) !向きによって片方で良い
  !u(i,j,k) = u(i.j,k) - dt/dx * (fluxr-fluxl)
  !-------------MUSCL solver-------------
end subroutine muslcslv1D

!subroutine vanalbada(fg,gradfg,iwx,iwy,iwz)
subroutine vanalbada(Phipre,Phigrad)
  use comvar
  double precision :: delp , delm
  integer :: i , ip , im , flmt ,eps=1.0d-20
  DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre

  write(*,*) 'vanalbada'
  do i = 0 , ndx-1
     ip=i+1
     im=i-1

     delp = Phipre(ip)-Phipre(i)
     delm = Phipre(i)-Phipre(im)
     !if(i==10) then
     !   write(*,*) delp , delm ,Phipre(ip),Phipre(i)Phipre(im)
     !end if
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !grdU(i,k) = flmt*( wave(ixp,jyp,kzp)-wave(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )
     Phigrad(i) = flmt*( Phipre(ip)-Phipre(im) )/( 2.0d0 * dx )
  end do
end subroutine vanalbada

subroutine saveu(in1)
  use comvar
  use grvvar
  integer :: i,in1
  character(5) name

  write(name,'(i5.5)') in1
  open(21,file='/Users/maeda/Desktop/kaiseki/testcode1/phi'//name//'.dat')
  do i=1,ndx-2
     write(21,*) i, Phi(i),Phi1step(i)
  end do
  close(21)
end subroutine saveu



subroutine fluxcal(preuse,pre,u,ep,kappa,mode)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx2) :: ul,ur,pre,slope,upreuse
  integer :: i,mode
  u(:)=0.0d0
  call vanalbada(pre,slop)
  if(mode==1) then
  do i = 0,ndx-1
  ul(i) = preuse(i) + 0.25 * ep * slop(i) * ((1.0d0-slop(i)*kappa)*(pre(i)-pre(i-1)) + (1+slop(i)*kappa)*(pre(i+1) - pre(i))) !i+1/2
  end do
  u(:)=ul(:)
  end if

  if(mode==2) then
  do i = 0,ndx-1
  ur(i) = preuse(i+1) - 0.25 * ep * slop(i) * ((1.0d0+slop(i)*kappa)*(pre(i+1)-pre(i)) + (1-slop(i)*kappa)*(pre(i+2) - pre(i+1))) !i+1/2
  end do
  u(:)=ur(:)
  end if

end subroutine fluxcal
