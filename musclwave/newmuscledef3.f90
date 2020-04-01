module comvar
  implicit none
  integer, parameter :: ndx=130,laststep=1000000,ist=1,ien=2,svnum=1000 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
  DOUBLE PRECISION :: cg = 1.0d3 , dx , ch = 1.d0 != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d2 , h=10.d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0!,h=10.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d-1 ,meanrho!, kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  DOUBLE PRECISION , dimension(-1:ndx) :: Phidt1,Phidt2
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=130 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phi,v!,w1,w2
  DOUBLE PRECISION , dimension(-1:ndx2) :: x , Phicgm ,rho, Phi1step , Phi2step ,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phidt,Phigrd,Phiexa,wp1,wp2,xi1,xi2
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phidtn , fx
  double precision , dimension(1:2) :: pl,pr
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0!,h
  integer :: i,sv=0,iws,ws=2,j
  DOUBLE PRECISION , dimension(-1:ndx2) :: xi1dummy,xi2dummy

  call INITIAL()
  !call BC(1,dt)
  call time(dt)
  !call BC(11,dt)
  call saveu(sv)

  do i=1,laststep
     call time(dt)
!     call timesource2(Phi,Phidt,rho,dt)
!     call timesource(xi1,rho,dt,1)
!     call timesource(xi2,rho,dt,1)
   !  do j=-1,ndx
   !     wp1(j)=Phi(j)-cg*v(j)
   !     wp2(j)=Phi(j)+cg*v(j)
   !  end do
     call timesource(Phidtn,rho,dt,3)
   !  h = ch*dt+h
   !  call exa(1)
   !  call BC(10,dt)
     !call exa(h,mode)
     do j=-1,ndx
        wp1(j)=(Phidtn(j)+cg*fx(j))/dsqrt(2.d0*cg)
        wp2(j)=(Phidtn(j)-cg*fx(j))/dsqrt(2.d0*cg)
     end do
!     call timesource(wp1,xi1,dt,2)
!     call timesource(wp2,xi2,dt,2)
     write(*,*) i,dt,'step'

     !---------xi part-------------
!     call muslcslv1D(xi1,rho,0.5d0*dt,1)
!     call muslcslv1D(xi2,rho,0.5d0*dt,2)
!     call BC(2)
!     do j=1,ndx-2
!        xi1(j) = xi1(j) - 0.5d0 * G4pi * cg * cg * dt * rho(j)
!        xi2(j) = xi2(j) - 0.5d0 * G4pi * cg * cg * dt * rho(j)
!     end do
!     call BC(2)
     !---------xi part-------------

     !do j=-1,ndx
     !   wp1(j)=Phi(j)-cg*v(j)
     !   wp2(j)=Phi(j)+cg*v(j)
     !end do
  !    do j=1,ndx-2
  !      wp1(j) = wp1(j) -  G4pi * cg * cg * dt * rho(j) / dsqrt(2.d0 * cg)
  !      wp2(j) = wp2(j) -  G4pi * cg * cg * dt * rho(j) / dsqrt(2.d0 * cg)
  !   end do

     call muslcslv1D(wp1,rho,dt,2)
     call muslcslv1D(wp2,rho,dt,1)
   !  call BC(1)
 !    do j=1,ndx-2
 !       wp1(j) = wp1(j) -  dt * xi2(j)
 !       wp2(j) = wp2(j) -  dt * xi1(j)
 !    end do
!     do j=1,ndx-2
!        Phi(j)=0.5d0*(wp1(j)+wp2(j))
!        v(j)=0.5d0*(-wp1(j)+wp2(j))/cg
!     end do

   !  do j=1,ndx-2
   !     wp1(j) = wp1(j) -  G4pi * cg * cg * dt * rho(j) / dsqrt(2.d0 * cg)
   !     wp2(j) = wp2(j) -  G4pi * cg * cg * dt * rho(j) / dsqrt(2.d0 * cg)
   !  end do
     do j=1,ndx-2
        Phidtn(j)=0.5d0*((wp1(j)+wp2(j))*dsqrt(2.d0*cg))
        fx(j)=0.5d0*((-wp1(j)+wp2(j))*dsqrt(2.d0*cg))/cg
     end do
 !    call BC(1)
 !    call BC(10,dt)
     do j=1,ndx-2
        Phidtn(j) = Phidtn(j) -  G4pi * cg * cg * dt * rho(j)! / Phidtn(i)
        !fx(j) = fx(j) +  G4pi  * cg * dt * rho(j)
     end do
 !    call BC(10,dt)
 !    call BC(20,dt)
     call BC(2,dt)
     !---------xi part-------------
!     call muslcslv1D(xi1,rho, dt,1)
!     call muslcslv1D(xi2,rho, dt,2)
!     call BC(2)
!     do j=1,ndx-2
!        xi1(j) = xi1(j) -  G4pi * cg * cg * dt * rho(j) / dsqrt(2.d0 * cg)
!        xi2(j) = xi2(j) -  G4pi * cg * cg * dt * rho(j) / dsqrt(2.d0 * cg)
!     end do
!    call BC(2)
     !---------xi part-------------

     if(mod(i,svnum)==0) then
        call saveu(sv)
     end if
  end do
  call BC(10,dt)
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
  Phi(:) = 0.d0
  Phiexa(:) = 0.d0
  !---------Phi-------------

  !---------Phidt-fx-------------
  Phidtn(:) = 0.d0
  fx(:)=0.d0
  !---------Phi-------------

  !----------v--------------
  v(:)=0.d0
  !----------v--------------

  !----------xi--------------
  xi1(:)=0.d0
  xi2(:)=0.d0
  !----------xi--------------

  !----------wp--------------
  wp1(:)=0.d0
  wp2(:)=0.d0
  !----------wp--------------



  !---------rho-------------
!  do i = -1,ndx
!     if( dabs(x(i) - hcen) .le. h) then
!        rho(i) = dinit1
        !rho(i) = 0.0d0
!     else
!        rho(i) = 0.0d0
        !rho(i) = dinit1
        !rho(i) = dinit1*1.d-2
!     end if
!  end do

  do i = -1, ndx
     amp =  dinit1
     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     rho(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
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

  Phiexa(:)=0.d0
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

end subroutine INITIAL



subroutine BC(mode,dt)
  use comvar
  use grvvar
  integer :: i,mode
  double precision  :: dt,dtt=0.d0
  double precision  :: fxm1=0.d0,fx0=0.d0,fxnxm1=0.d0,fxnx=0.d0
  !double precision , dimension(1:2) :: pl,pr

  if(mode==1) then

     goto 111
     Phi(   0)=Phi(ndx-2)
     Phi(  -1)=Phi(ndx-3)
     Phi(ndx-1)=Phi(   1)
     Phi(ndx  )=Phi(   2)

     v(   0)=v(ndx-2)
     v(  -1)=v(ndx-3)
     v(ndx-1)=v(   1)
     v(ndx  )=v(   2)
111  continue
     wp1(   0)=wp1(ndx-2)
     wp1(  -1)=wp1(ndx-3)
     wp1(ndx-1)=wp1(   1)
     wp1(ndx  )=wp1(   2)

     wp2(   0)=wp2(ndx-2)
     wp2(  -1)=wp2(ndx-3)
     wp2(ndx-1)=wp2(   1)
     wp2(ndx  )=wp2(   2)
  end if

  if(mode==2) then
     Phidtn(   0)=Phidtn(ndx-2)
     Phidtn(  -1)=Phidtn(ndx-3)
     Phidtn(ndx-1)=Phidtn(   1)
     Phidtn(ndx  )=Phidtn(   2)
     !Phidtn(   0)=0.d0
     !Phidtn(  -1)=0.d0
     !Phidtn(ndx-1)=0.d0
     !Phidtn(ndx  )=0.d0

     fx(   0)=fx(ndx-2)
     fx(  -1)=fx(ndx-3)
     fx(ndx-1)=fx(   1)
     fx(ndx  )=fx(   2)
     !fx(   0)=Phigrd(0)
     !fx(  -1)=Phigrd(-1)
     !fx(ndx-1)=Phigrd(ndx-1)
     !fx(ndx  )=Phigrd(ndx)
!     xi1(   0)=xi1(ndx-2)
!     xi1(  -1)=xi1(ndx-3)
!     xi1(ndx-1)=xi1(   1)
!     xi1(ndx  )=xi1(   2)

!     xi2(   0)=xi2(ndx-2)
!     xi2(  -1)=xi2(ndx-3)
!     xi2(ndx-1)=xi2(   1)
!     xi2(ndx  )=xi2(   2)
  end if

  if(mode==10) then
     !write(*,*) 'dtt',dtt,dt
     Phidtn(0) = (Phiexa(0)-pl(1))/dtt
     Phidtn(-1) = (Phiexa(-1)-pl(2))/dtt
     Phidtn(ndx-1) = (Phiexa(ndx-1)-pr(1))/dtt
     Phidtn(ndx) = (Phiexa(ndx)-pr(2))/dtt

     fx(0)    =(Phiexa(1)-Phiexa(0))/dx
     fx(-1)   =(Phiexa(0)-Phiexa(-1))/dx
     fx(ndx-1)=(Phiexa(ndx-1)-Phiexa(ndx-2))/dx
     fx(ndx)  =(Phiexa(ndx)-Phiexa(ndx-1))/dx

     !write(*,*)'BC',Phidtn(0),Phiexa(0),pl(1),dtt
     pl(1) = Phiexa(0)
     pl(2) = Phiexa(-1)
     pr(1) = Phiexa(ndx-1)
     pr(2) = Phiexa(ndx)

     !write(*,*),'BC',Phidtn(0)

     dtt=dt
  end if


  if(mode==11) then
     pl(1) =0.d0
     pl(2) =0.d0
     pr(1) =0.d0
     pr(2) =0.d0

     dtt=dt
    !write(*,*)'bcini' ,dtt,pl(1)
  end if

  if(mode==20) then
    ! fxm1=fx(1)-fx(0)
    ! fx0=fx(0)-fx(-1)
    ! fxnxm1=fx(ndx-1)-fx(ndx-2)
    ! fxnx=fx(ndx)-fx(ndx-1)

     fx(   0)=fx(ndx-2)
     fx(  -1)=fx(ndx-3)
     fx(ndx-1)=fx(   1)
     fx(ndx  )=fx(   2)

     !Phidtn(0) = Phidtn(0)+cg*cg*dt*(-G4pi*rho(0)+(fx(1)-fx(0))/dx)
     !Phidtn(-1) =  Phidtn(-1)+cg*cg*dt*(-G4pi*rho(-1)+(fx(0)-fx(-1))/dx)
     !Phidtn(ndx-1) =  Phidtn(ndx-1)+cg*cg*dt*(-G4pi*rho(ndx)+(fx(ndx-1)-fx(ndx-2))/dx)
     !Phidtn(ndx) =  Phidtn(ndx)+cg*cg*dt*(-G4pi*rho(ndx)+(fx(ndx)-fx(ndx-1))/dx)

     !Phidtn(0) = Phidtn(0)+cg*cg*dt*(-G4pi*rho(0)+0.5d0*(fx(1)-2.0d0*fx(0)+fx(-1))/dx)
     !Phidtn(-1) =  Phidtn(-1)+cg*cg*dt*(-G4pi*rho(-1)+0.5d0*(fx(0)-2.0d0*fx(-1)+fx(ndx))/dx)
     !Phidtn(ndx-1) =  Phidtn(ndx-1)+cg*cg*dt*(-G4pi*rho(ndx)+0.5d0*(fx(ndx)-2.0d0*fx(ndx-1)+fx(ndx-2))/dx)
     !Phidtn(ndx) =  Phidtn(ndx)+cg*cg*dt*(-G4pi*rho(ndx)+0.5d0*(fx(-1)-2.0d0*fx(ndx)+fx(ndx-1))/dx)

!     Phidtn(0) = Phidtn(0)+cg*cg*dt*(-G4pi*rho(0)+0.25d0*(fx(1)-2.0d0*fx(0)+fx(-1)+fxm1)/dx)
!     Phidtn(-1) =  Phidtn(-1)+cg*cg*dt*(-G4pi*rho(-1)+0.25d0*(fx(0)-2.0d0*fx(-1)+fx(ndx)+fx0)/dx)
!     Phidtn(ndx-1) =  Phidtn(ndx-1)+cg*cg*dt*(-G4pi*rho(ndx)+0.25d0*(fx(ndx)-2.0d0*fx(ndx-1)+fx(ndx-2)+fxnxm1)/dx)
!     Phidtn(ndx) =  Phidtn(ndx)+cg*cg*dt*(-G4pi*rho(ndx)+0.25d0*(fx(-1)-2.0d0*fx(ndx)+fx(ndx-1)+fxnx)/dx)

     Phidtn(0) = Phidtn(0)+cg*cg*dt*(-G4pi*rho(0)+0.25d0*(fx(1)-fx(-1)+fxm1)/dx)
     Phidtn(-1) =  Phidtn(-1)+cg*cg*dt*(-G4pi*rho(-1)+0.25d0*(fx(0)-fx(ndx)+fx0)/dx)
     Phidtn(ndx-1) =  Phidtn(ndx-1)+cg*cg*dt*(-G4pi*rho(ndx-1)+0.25d0*(fx(ndx)-fx(ndx-2)+fxnxm1)/dx)
     Phidtn(ndx) =  Phidtn(ndx)+cg*cg*dt*(-G4pi*rho(ndx)+0.25d0*(fx(-1)-fx(ndx-1)+fxnx)/dx)

!     fxm1=0.5d0*(fx(1)-2.0d0*fx(0)+fx(-1))/dx
!     fx0=0.5d0*(fx(0)-2.0d0*fx(-1)+fx(ndx))/dx
!     fxnxm1=0.5d0*(fx(ndx)-2.0d0*fx(ndx-1)+fx(ndx-2))/dx
!     fxnx=0.5d0*(fx(-1)-2.0d0*fx(ndx)+fx(ndx-1))/dx

     fxm1=(fx(1)-fx(-1))
     fx0=(fx(0)-fx(ndx))
     fxnxm1=(fx(ndx)-fx(ndx-2))
     fxnx=(fx(-1)-fx(ndx-1))

    ! fx(   0)=fx(ndx-2)
    ! fx(  -1)=fx(ndx-3)
    ! fx(ndx-1)=fx(   1)
    ! fx(ndx  )=fx(   2)
  end if
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
        !write(*,*) nu2 * (Phiu(i) - Phiu(i-1))
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
     !write(*,*) 'mode1'

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

  do i=1,ndx-2
     Phi(i)=0.5d0*(dsqrt(2.d0*cg)*wp1(i)+dsqrt(2.d0*cg)*wp2(i))
     v(i)=0.5d0*(-dsqrt(2.d0*cg)*wp1(i)+dsqrt(2.d0*cg)*wp2(i))/cg
  end do

  write(name,'(i5.5)') in1
  open(21,file='/Users/maeda/Desktop/Dropbox/kaiseki-desktpo/testmuscledf/phi'//name//'.dat')
  do i=1,ndx-2
        write(21,*) x(i),Phidtn(i),fx(i),wp1(i),wp2(i) ,rho(i),Phiexa(i)!Phi(i), v(i),wp1(i),wp2(i),xi1(i),xi2(i)!Phicgp(i) , Phicgm(i) ,rho(i)!,&
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
  end if


  if(mode==4) then
     do i = ist-1,ndx-ien+1
        ur(i) = preuse(i) - 0.25d0 * ep * slop(i) &
             * ((1.0d0+slop(i)*kappa)*(pre(i)-pre(i-1)) + (1.0d0-slop(i)*kappa)*(pre(i+1) - pre(i))) !i-1/2
        u(i)=ur(i)
     end do
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


subroutine timesource(Phiv,source,dt,mode)
  use comvar
  !use grvver
  integer i,mode
  double precision :: dt,sdt1,mindt,maxdt,std2,rt=0.5d0
  DOUBLE PRECISION, dimension(-1:ndx) :: Phiv,source,Phidt
  !maxdt=0.d0
  !mindt=1.d10

  if(mode==1)then
     maxdt=0.d0
     mindt=1.d10
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           sdt = dabs( 0.5d0 * Phiv(i) /( G4pi * cg * cg  * source(i)/ dsqrt(2.d0*cg) ))
           mindt=dmin1(mindt,sdt)
        end if
     end do
     !if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
     !   dt = sdt
     !end if
     if( (mindt < dt) .and. (mindt .ne. 0.0d0)) then
        dt=dmin1(mindt,dt)
        !dt = sdt
     end if
  end if

  if(mode==2)then
     maxdt=0.d0
     mindt=1.d10
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           sdt = dabs(0.5d0 * Phiv(i) /( source(i) ))
           mindt=dmin1(mindt,sdt)
        end if
     end do
     !if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
     !   dt = sdt
     !end if
     if( (mindt < dt) .and. (mindt .ne. 0.0d0)) then
        dt=dmin1(mindt,dt)
        !dt = sdt
     end if
  end if

   if(mode==3)then
     maxdt=0.d0
     mindt=1.d10
     do i=1,ndx-2
        if((source(i) .ne. 0.0d0) .and. (Phiv(i) .ne. 0.0d0))then
           sdt = dabs(0.5d0 * Phiv(i) /( G4pi * cg * cg  * source(i) ))
           mindt=dmin1(mindt,sdt)
        end if
     end do
     !if( (maxdt < dt) .and. (maxdt .ne. 0.0d0)) then
     !   dt = sdt
     !end if
     if( (mindt < dt) .and. (mindt .ne. 0.0d0)) then
        dt=dmin1(mindt,dt)
        !dt = sdt
     end if
  end if

  write(*,*) 'time source' , dt
end subroutine timesource


subroutine exa(mode)
  use comvar
  use grvvar
  integer i,mode
  !double precision h


  if(mode==1)then
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

     !---------rho-------------
  end if


  if(mode==3)then
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
           rho(i) = rho(i)-meanrho
     end do
     !---------rho-------------
  end if


  !--------Phiexa-----------
  !goto 200
!  open(142,file=dir//'phiexact.DAT')
!  open(143,file=dir//'INIden.DAT')
!  open(144,file=dir//'phigrd.DAT')
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        !write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
     else
        Phiexa(i) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        !write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
     end if
     !write(143,*) sngl(rho(i))
  end do
  !write(*,*)'exa',h,Phiexa(0)
  !--------Phiexa-----------
end subroutine exa
