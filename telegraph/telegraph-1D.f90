module comvar
  implicit none
  integer, parameter :: ndx=130,laststep=1000,ist=1,ien=2,svnum=10 !preiodic:ist=1,ien=2 , kotei:ist=2,ien=3 : ndx=130
  !double precision, parameter :: Lbox=1.0d2 , h=10.0d0 , hcen=50.0d0 , dinit1=1.29988444d0,w1=2.0d0
DOUBLE PRECISION :: cg = 1.d0 , dx, Tdiff=0.10d0 != Lbox/dble(ndx-2) !, bcphi1 , bcphi2
  double precision :: Lbox=1.0d0 , h=0.2d0 , hcen=0.5d0 , dinit1=1.29988444d0,w1=2.0d0
  !double precision :: G=1.11142d-4, G4pi=12.56637d0*G , coeff=0.90d0 ,  kappa=1.0d0/3.0d0
  double precision ::  G4pi=12.56637d0*1.11142d-4 , coeff=0.5d0 ,meanrho!,  kappa=1.0d0/3.0d0
  DOUBLE PRECISION , dimension(1:3) :: bcphi1 , bcphi2 ,bcphigrd1 , bcphigrd2
  double precision :: osromega=0.001d0*3.1415926535d0
  !character(63) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-1D-1/simu/'
  character(86) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-1D-1/cg1-T01-rho1-cen-128-1000st/'
end module comvar

module grvvar
  implicit none
  integer, parameter :: ndx2=130 !パラメータ属性必要
  DOUBLE PRECISION , dimension(-1:ndx2) :: x,Phicgm,rho,Phi1step,Phi2step,Phicgp
  DOUBLE PRECISION , dimension(-1:ndx2) :: Phidt,Phigrd,Phiexa,dmyphi,dmystep
end module grvvar

program muscl1D
  !implicit none まちがった位置
  use comvar
  use grvvar
  implicit none
  DOUBLE PRECISION :: dt=0.0d0
  integer :: i,sv=0,iws,ws=2,ii


  call INITIAL()
  call BC(1)
  !call muslcslv1D(Phi,Phi1step,dt,13)
  call saveu(sv)

  do i=1,laststep
     if(mod(i,svnum)==0) then
        call saveu(sv)
     end if

     call time(dt)

     !call BC(4)
     !call BC(3)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,3)

     !call osr(dt*0.5d0,2,1)
     !call BC(5)
     do ii=1,ndx-2
        !Phi1step(ii) = dt*0.5d0*(Phicgp(ii)-Phi1step(ii))*0.5d0/Tdiff &
        !+dt*0.5d0*cg*cg*2.d0*Tdiff*G4pi*rho(ii) +Phi1step(ii) !dt*cg^2*2T*4piG
        
        !non-iso
        Phi1step(ii) = Phi1step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        +(Phicgp(ii) - cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))
        
        !iso
        !Phi1step(ii) = Phi1step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        !+(Phicgp(ii) - cg*cg*2.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))
     enddo
     do ii=1,ndx-2
        !Phi2step(ii) = dt*0.5d0*(Phicgm(ii)-Phi2step(ii))*0.5d0/Tdiff &
        !+dt*0.5d0*cg*cg*2.d0*Tdiff*G4pi*rho(ii) +Phi2step(ii) !dt*-cg^2*2T*4piG

        !non-iso
        Phi2step(ii) = Phi2step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        +(Phicgm(ii) - cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))

        !iso
        !Phi2step(ii) = Phi2step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        !+(Phicgm(ii) - cg*cg*2.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))
     enddo
     !call muslcslv1D(Phi2step,rho,dt*0.5d0,3)
     !call muslcslv1D(Phi1step,rho,0.5d0*dt,3)
     !call BC(4)
     !call BC(3)
     call BC(155)
     !call BC(255)
     !call osr(dt*0.5d0,2,0)
     !call BC(5)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,1)

 
     !call BC(4)
     !call BC(3)
     !call BC(1)
     !call BC(55)
     !Phidt(:)=Phi(:)
!     call muslcslv1D(Phicgp,Phi2step,dt,4)
     !call muslcslv1D(Phicgm,Phi1step,dt,4)
     !call osr(dt*1.d0,1,1)
     
     do ii=1,ndx-2
     !Phicgp(ii) = dt*(-Phicgp(ii)+Phi1step(ii))*0.5d0/Tdiff +Phicgp(ii)

     !non-iso
     Phicgp(ii) = Phicgp(ii)*dexp(-0.5d0/Tdiff * dt) &
     + Phi1step(ii)*(1.d0-dexp(-0.5d0/Tdiff * dt))

     !iso
     !Phicgp(ii) = Phicgp(ii)*dexp(-0.5d0/Tdiff * dt) &
     !+(Phi1step(ii) - cg*cg*2.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt))
     enddo
     do ii=1,ndx-2
     !Phicgm(ii) = dt*(-Phicgm(ii)+Phi2step(ii))*0.5d0/Tdiff +Phicgm(ii)

     !non-iso
     Phicgm(ii) = Phicgm(ii)*dexp(-0.5d0/Tdiff * dt) &
     + Phi2step(ii)*(1.d0-dexp(-0.5d0/Tdiff * dt))

     !iso
     !Phicgm(ii) = Phicgm(ii)*dexp(-0.5d0/Tdiff * dt) &
     !+(Phi2step(ii) - cg*cg*2.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt))
     enddo

     !call osr(dt*1.d0,1,0)
     !call BC(1)
     call BC(156)
     !call BC(256)
     !call BC(55)
     call muslcslv1D(Phicgp,Phi1step,dt,1)
     call muslcslv1D(Phicgm,Phi2step,dt,2)
     !call muslcslv1D(Phicgp,Phi1step,0.5d0*dt,1)
     !call muslcslv1D(Phicgm,Phi2step,0.5d0*dt,2)
     
     !*****-4piGrho******
     !dmyphi(:)=Phicgp(:)
     !dmystep(:)=Phi1step(:)
     !do ii=1,ndx-2
     !Phicgp(ii)  =0.5d0*dmyphi(ii) *(1.d0+dexp(-dt/Tdiff))+0.5d0*dmystep(ii)*(1.d0-dexp(-dt/Tdiff)) &
     !+cg*cg*Tdiff*Tdiff*G4pi*rho(ii)*(1.d0-dexp(-dt/Tdiff))-cg*cg*Tdiff*G4pi*rho(ii)*dt
     !Phi1step(ii)=0.5d0*dmystep(ii)*(1.d0+dexp(-dt/Tdiff))+0.5d0*dmyphi(ii)* (1.d0-dexp(-dt/Tdiff)) &
     !-cg*cg*Tdiff*Tdiff*G4pi*rho(ii)*(1.d0-dexp(-dt/Tdiff))-cg*cg*Tdiff*G4pi*rho(ii)*dt
     !enddo
     !dmyphi(:)=Phicgm(:)
     !dmystep(:)=Phi2step(:)
     !do ii=1,ndx-2
     !Phicgm(ii)  =0.5d0*dmyphi(ii) *(1.d0+dexp(-dt/Tdiff))+0.5d0*dmystep(ii)*(1.d0-dexp(-dt/Tdiff)) &
     !-cg*cg*Tdiff*Tdiff*G4pi*rho(ii)*(1.d0-dexp(-dt/Tdiff))-cg*cg*Tdiff*G4pi*rho(ii)*dt
     !Phi2step(ii)=0.5d0*dmystep(ii)*(1.d0+dexp(-dt/Tdiff))+0.5d0*dmyphi(ii)* (1.d0-dexp(-dt/Tdiff)) &
     !+cg*cg*Tdiff*Tdiff*G4pi*rho(ii)*(1.d0-dexp(-dt/Tdiff))-cg*cg*Tdiff*G4pi*rho(ii)*dt
     !enddo

    !*****-4piGrho******
    !dmyphi(:)=Phicgp(:)
    !dmystep(:)=Phi1step(:)
    !do ii=1,ndx-2
    !Phicgp(ii)  =0.5d0*dmyphi(ii) *(1.d0+dexp(-dt/Tdiff))+0.5d0*dmystep(ii)*(1.d0-dexp(-dt/Tdiff)) &
    !-cg*cg*Tdiff*G4pi*rho(ii)*dt
    !Phi1step(ii)=0.5d0*dmystep(ii)*(1.d0+dexp(-dt/Tdiff))+0.5d0*dmyphi(ii)* (1.d0-dexp(-dt/Tdiff)) &
    !-cg*cg*Tdiff*G4pi*rho(ii)*dt
    !enddo
    !dmyphi(:)=Phicgm(:)
    !dmystep(:)=Phi2step(:)
    !do ii=1,ndx-2
    !Phicgm(ii)  =0.5d0*dmyphi(ii) *(1.d0+dexp(-dt/Tdiff))+0.5d0*dmystep(ii)*(1.d0-dexp(-dt/Tdiff)) &
    !-cg*cg*Tdiff*G4pi*rho(ii)*dt
    !Phi2step(ii)=0.5d0*dmystep(ii)*(1.d0+dexp(-dt/Tdiff))+0.5d0*dmyphi(ii)* (1.d0-dexp(-dt/Tdiff)) &
    !-cg*cg*Tdiff*G4pi*rho(ii)*dt
    !enddo

     !call BC(256)
     !call BC(55)
     !call BC(1)
     !call muslcslv1D(Phicgp,Phi1step,0.5d0*dt,1)
     !call muslcslv1D(Phicgm,Phi2step,0.5d0*dt,2)
     !call BC(1)

     !call BC(4)
     !call BC(3)
     !call osr(dt*0.5d0,2,1)
     !call BC(5)
     !call muslcslv1D(Phi1step,rho,dt,3)
     !call muslcslv1D(Phi1step,rho,dt*0.5d0,3)
!     call muslcslv1D(Phi2step,rho,0.5d0*dt,3)

     
    do ii=1,ndx-2
        !Phi1step(ii) = dt*0.5d0*(Phicgp(ii)-Phi1step(ii))*0.5d0/Tdiff &
        !+dt*0.5d0*cg*cg*2.d0*Tdiff*G4pi*rho(ii) +Phi1step(ii) !dt*cg^2*2T*4piG
        
        !non-iso
        Phi1step(ii) = Phi1step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        +(Phicgp(ii) - cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))

        !iso
        !Phi1step(ii) = Phi1step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        !+(Phicgp(ii) - cg*cg*2.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))
     enddo
     do ii=1,ndx-2
        !Phi2step(ii) = dt*0.5d0*(Phicgm(ii)-Phi2step(ii))*0.5d0/Tdiff &
        !+dt*0.5d0*cg*cg*2.d0*Tdiff*G4pi*rho(ii) +Phi2step(ii) !dt*-cg^2*2T*4piG

        !non-iso
        Phi2step(ii) = Phi2step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        +(Phicgm(ii) - cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))

        !iso
        !Phi2step(ii) = Phi2step(ii)*dexp(-0.5d0/Tdiff * dt*0.5d0) &
        !+(Phicgm(ii) - cg*cg*2.d0*Tdiff*Tdiff*G4pi*rho(ii))*(1.d0-dexp(-0.5d0/Tdiff * dt*0.5d0))
     enddo
     !call BC(4)
     !call BC(3)
     call BC(155)
     !call BC(255)
     !call osr(dt*0.5d0,2,0)
     !call BC(5)
     !call muslcslv1D(Phi1step,rho,dt,1)
     call muslcslv1D(Phi1step,rho,dt*0.5d0,2)
     call muslcslv1D(Phi2step,rho,dt*0.5d0,1)
     !call BC(3)
     !call BC(4)
     !call BC(1)

     !write(*,*)'end'
  end do
  !call BC(3)
  !call BC(4)
  !call BC(1)
  !call BC(5)
  !call BC(55)
  call BC(155)
  call BC(156)
  !call BC(255)
  !call BC(256)
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
  double precision :: amp,pi=3.1415926535d0,haba,k,meanphi

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


  !---------Phi-------------
  Phicgp(:)=0.0d0
  Phicgm(:)=0.0d0

  !amp=1.d-3
  !k=3.1415926535d0*2.d0/Lbox * 5.d0
  !do i=-1,ndx
  !Phicgp(i)=amp*dsin(k*x(i))
  !Phicgm(i)=amp*dsin(k*x(i))
  !enddo
  !---------Phi-------------

  !-------Phi1step-----------
  Phi1step(:)=0.0d0
  Phi2step(:)=0.0d0
  !Phi1step(:)=+G4pi*meanrho*cg*Lbox
  !Phi2step(:)=0.0d0
  !do i=-1,ndx
  !Phi1step(i)=amp*dsin(k*x(i))+2.d0*Tdiff*(amp*cg/k*dcos(k*x(i)) - amp*cg*k*dcos(k*x(i))) !b+2T(db/dt-c*db/cx)
  !Phi2step(i)=amp*dsin(k*x(i))+2.d0*Tdiff*(amp*cg/k*dcos(k*x(i)) + amp*cg*k*dcos(k*x(i)))
  !enddo
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
     rho(i)=rho(i)!-meanrho
  end do
  


  !Phi1step(:)=0.d0 !+G4pi*meanrho*cg*Lbox
  !Phi2step(:)=0.d0 !+G4pi*meanrho*cg*Lbox

  !rho(:)=0.d0
  !---------rho-------------



  !--------Phiexa-----------
  !goto 200
  open(142,file=dir//'phiexact.DAT')
  open(143,file=dir//'INIden.DAT')
  open(144,file=dir//'phigrd.DAT')
  meanphi=0.d0
  do i= -1,ndx
     if( dabs(x(i) - hcen) .le. h ) then
        Phiexa(i) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2
        !Phiexa(i) = G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2 - G4pi * dinit1 * (x(i) - hcen )
        write(142,*) sngl(x(i)) ,  sngl(G4pi/2.0d0 * dinit1 * (x(i) - hcen )**2)
     else
        Phiexa(i) = G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2
        !Phiexa(i) = G4pi * dinit1 * h * x(i)  - G4pi/2.0d0 * dinit1 * h**2 - G4pi * dinit1 * h * hcen
        !Phiexa(i) = G4pi * dinit1 * (h - hcen) * x(i)  + G4pi/2.0d0 * dinit1 * (-h**2+hcen**2)
        write(142,*) sngl(x(i)) , sngl(G4pi * dinit1 * h * dabs(x(i) - hcen)  - G4pi/2.0d0 * dinit1 * h**2)
     end if
     write(143,*) sngl(rho(i))
     meanphi=meanphi+Phiexa(i)
  end do

  !do i= -1,ndx
  !   Phiexa(i)=Phiexa(i)-meanphi
  !   rho(i)=rho(i)-meanrho
  !enddo


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
!  goto 201
!  do i = -1, ndx
!     amp = 1.d-3
!     !Phi(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
!     Phi1step(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
!     Phi2step(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
!     Phicgp(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
!     Phicgm(i) =  amp*dsin(2.d0*pi*x(i)/Lbox)
!  end do

  goto 201
  do i = -1, ndx
     amp = 1.d-3
     haba=10.0d0
     !Phi(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     Phi1step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     Phi2step(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     Phicgp(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
     Phicgm(i) =  amp*dexp(-(x(i) - 0.5d0*Lbox)**2 /(2.0d0 * haba**2))
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

write(*,*)'INITIAL'
end subroutine INITIAL

subroutine osr(dt,mode,tmode)
use comvar
use grvvar
double precision :: dt,t1=0.d0,t2=0.d0,amp1=1.d0
integer mode,tmode

if(mode==1)then
if(tmode==1)then
t1=dt+t1
!write(*,*)'aa',t1,dt
endif

Phicgp((ndx-1)/2)=amp1*dsin(osromega*t1)
Phicgm((ndx-1)/2)=amp1*dsin(osromega*t1)
!Phicgp(ndx/2)=amp1*dsin(osromega*t1)
!Phicgm(ndx/2)=amp1*dsin(osromega*t1)
!write(*,*)Phicgp((ndx-1)/2),Phi1step((ndx-1)/2)
!write(*,*)'aa'
endif

if(mode==2) then
if(tmode==1)then
t2=dt+t2
!write(*,*)'bb',t2,dt
endif

Phi1step((ndx-1)/2)=2.d0*Tdiff*(amp1*osromega*dcos(osromega*t2)- cg*0.5d0*(Phicgp(ndx/2+1)-Phicgp(ndx/2-1))/dx &
                +0.5d0/Tdiff*dsin(osromega*t2))
Phi2step((ndx-1)/2)=2.d0*Tdiff*(amp1*osromega*dcos(osromega*t2)+ cg*0.5d0*(Phicgm(ndx/2+1)-Phicgm(ndx/2-1))/dx &
                +0.5d0/Tdiff*dsin(osromega*t2))
!write(*,*)'bb'
endif
write(*,*)'osr',Phicgp((ndx-1)/2),Phi1step((ndx-1)/2),t1,t2
end subroutine osr




subroutine BC(mode)
  use comvar
  use grvvar
  integer :: i,mode
  double precision , dimension(1:2) :: pl,pr

  if(mode==1) then
     Phicgp(ndx)=Phicgp(2)
     Phicgp(ndx-1)= Phicgp(1)
     !Phicgp(ndx-2)= Phicgp(ndx-3)
     Phicgp(-1)= Phicgp(ndx-3)
     Phicgp(0)= Phicgp(ndx-2)
     !Phicgp(1)= Phicgp(2)

     Phicgm(ndx)=Phicgm(2)
     Phicgm(ndx-1)= Phicgm(1)
     !Phicgm(ndx-2)= Phicgm(ndx-3)
     Phicgm(-1)= Phicgm(ndx-3)
     Phicgm(0)= Phicgm(ndx-2)
     !Phicgm(1)= Phicgm(2)

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
     Phi1step(-1)= Phi1step(ndx-3)
     Phi1step(0)= Phi1step(ndx-2)
     !Phi1step(ndx-2)= Phigrd(ndx-2)
     !Phi1step(ndx-1)= Phigrd(ndx-1)
     !Phi1step(ndx)= Phigrd(ndx)
   !  Phi1step(ndx)= Phi1step(ndx-1)
   !  Phi1step(ndx-1)= Phi1step(ndx-2)
   !  Phi1step(ndx-2)= Phi1step(ndx-3)
     !Phi1step(ndx/2)=0.0d0
     !Phi1step(ndx/2-1)=0.0d0
     Phi1step(ndx)= Phi1step(2)
     Phi1step(ndx-1)= Phi1step(1)
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

     Phi2step(-1)= Phi2step(ndx-3)
     Phi2step(0)=  Phi2step(ndx-2)
     Phi2step(ndx)=Phi2step(2)
     Phi2step(ndx-1)=Phi2step(1)

     !Phi1step(ndx)= -Phigrd(ndx-2)
     !Phi2step(1)= -Phigrd(1)
     !Phi2step(0)= -Phigrd(1)
     !Phi2step(-1)=-Phigrd(1)
     !Phi2step(ndx-2)= -Phigrd(ndx-2)
     !Phi2step(ndx-1)= -Phigrd(ndx-2)
     !Phi2step(ndx)= -Phigrd(ndx-2)
     !-------Phi1step-----------
  end if

  if(mode==5) then
     !--------free--------------
!     goto 112
     !-------Phi1step-----------
!     Phi1step(1)= Phi1step(2)
     Phi1step(-1)=Phi1step(0)
     Phi1step(0)= Phi1step(1)
!     Phi1step(ndx-2)= Phi1step(ndx-3)
     Phi1step(ndx)= Phi1step(ndx-1)
     Phi1step(ndx-1)= Phi1step(ndx-2)

     Phi2step(-1)=Phi2step(0)
     Phi2step(0) =Phi2step(1)
     Phi2step(ndx)  = Phi2step(ndx-1)
     Phi2step(ndx-1)= Phi2step(ndx-2)
     
     !-------Phi1step-----------
!112  continue
     !--------free--------------
  end if
  if(mode==55) then
       Phicgp(-1)=Phicgp(0)
       Phicgp(0) =Phicgp(1)
       Phicgp(ndx)  = Phicgp(ndx-1)
       Phicgp(ndx-1)= Phicgp(ndx-2)

       Phicgm(-1)=Phicgm(0)
       Phicgm(0) =Phicgm(1)
       Phicgm(ndx)  = Phicgm(ndx-1)
       Phicgm(ndx-1)= Phicgm(ndx-2)
    end if

   if(mode==155) then
      Phi1step(-1)=+cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1) ! b=cg/kappa * da/dx + a
      Phi1step(0) =+cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
      Phi1step(ndx)  = +cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
      Phi1step(ndx-1)= +cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
      !Phi1step(-1)=-cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
      !Phi1step(0) =-cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
      !Phi1step(ndx)  = -cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1)
      !Phi1step(ndx-1)= -cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
      write(*,*) -cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx), -cg*2.d0*Tdiff*Phigrd(ndx), Phiexa(ndx), Phigrd(ndx)

      Phi2step(-1)=-cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1) ! b=-cg/kappa * db/dx + b
      Phi2step(0) =-cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
      Phi2step(ndx)  = -cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
      Phi2step(ndx-1)= -cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
      !Phi2step(-1)=+cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
      !Phi2step(0) =+cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
      !Phi2step(ndx)  = +cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1)
      !Phi2step(ndx-1)= +cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)

    !Phi1step(-1)=Phiexa(-1)
    !Phi1step(0) =Phiexa(0)
    !Phi1step(ndx)  =Phiexa(ndx)
    !Phi1step(ndx-1)=Phiexa(ndx-1)

    !Phi2step(-1)=Phiexa(-1)
    !Phi2step(0) =Phiexa(0)
    !Phi2step(ndx)  =Phiexa(ndx)
    !Phi2step(ndx-1)=Phiexa(ndx-1)

     !Phi1step(-1)=Phi1step(0)
     !Phi1step(0) =Phi1step(1)
     !Phi1step(ndx)  = -cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
     !Phi1step(ndx-1)= -cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)

     !Phi2step(-1)=+cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1)
     !Phi2step(0) =+cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
     !Phi2step(ndx)  = Phi2step(ndx-1)
     !Phi2step(ndx-1)= Phi2step(ndx-2)


        !Phi1step(-1)=(-cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1))*2.d0
        !Phi1step(0) =(-cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0))*2.d0
        !Phi1step(ndx)  = (-cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx))*2.d0
        !Phi1step(ndx-1)= (-cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1))*2.d0

        !Phi2step(-1)=(+cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1))*2.d0
        !Phi2step(0) =(+cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0))*2.d0
        !Phi2step(ndx)  =( +cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx))*2.d0
        !Phi2step(ndx-1)=( +cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1))*2.d0

      
    !Phi1step(ndx-1)=cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
    !Phi1step(ndx)  =cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
    !Phi1step(0)  = -cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
    !Phi1step(-1) = -cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1)

    !Phi2step(ndx-1)=-cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
    !Phi2step(ndx)  =-cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
    !Phi2step(0)  = cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
    !Phi2step(-1) = cg*2.d0*Tdiff*Ph                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   igrd(-1)+Phiexa(-1)
   end if

if(mode==156) then
Phicgp(-1)=Phiexa(-1)
Phicgp(0) =Phiexa(0)
Phicgp(ndx)  = Phiexa(ndx)
Phicgp(ndx-1)= Phiexa(ndx-1)

Phicgm(-1)=Phiexa(-1)
Phicgm(0) =Phiexa(0)
Phicgm(ndx)  = Phiexa(ndx)
Phicgm(ndx-1)= Phiexa(ndx-1)

   !Phicgp(-1)=-cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1)
   !Phicgp(0) =-cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
   !Phicgp(ndx)  = -cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
   !Phicgp(ndx-1)= -cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)

   !Phicgm(-1)=+cg*2.d0*Tdiff*Phigrd(-1)+Phiexa(-1)
   !Phicgm(0) =+cg*2.d0*Tdiff*Phigrd(0)+Phiexa(0)
   !Phicgm(ndx)  = +cg*2.d0*Tdiff*Phigrd(ndx)+Phiexa(ndx)
   !Phicgm(ndx-1)= +cg*2.d0*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)

   !Phicgp(-1)=Phiexa(-1)
   !Phicgp(0) =Phiexa(0)
   !Phicgp(ndx)  = Phicgp(ndx-1)
   !Phicgp(ndx-1)= Phicgp(ndx-2)

   !Phicgm(-1)=Phicgm(0)
   !Phicgm(0) =Phicgm(1)
   !Phicgm(ndx)  = Phiexa(ndx)
   !Phicgm(ndx-1)= Phiexa(ndx-1)

   !Phicgp(-1)=Phiexa(-1)*2.d0
   !Phicgp(0) =Phiexa(0)*2.d0
   !Phicgp(ndx)  = Phiexa(ndx)*2.d0
   !Phicgp(ndx-1)= Phiexa(ndx-1)*2.d0

   !Phicgm(-1)=Phiexa(-1)*2.d0
   !Phicgm(0) =Phiexa(0)*2.d0
   !Phicgm(ndx)  = Phiexa(ndx)*2.d0
   !Phicgm(ndx-1)= Phiexa(ndx-1)*2.d0
end if
  

   if(mode==255) then
      Phi2step(-1)=-cg*Tdiff*Phigrd(-1)+Phiexa(-1) ! b=cg/kappa * da/dx + a
      Phi2step(0) =-cg*Tdiff*Phigrd(0)+Phiexa(0)
      Phi2step(ndx)  = -cg*Tdiff*Phigrd(ndx)+Phiexa(ndx)
      Phi2step(ndx-1)= -cg*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
      

      Phi1step(-1)=+cg*Tdiff*Phigrd(-1)+Phiexa(-1) ! b=-cg/kappa * db/dx + b
      Phi1step(0) =+cg*Tdiff*Phigrd(0)+Phiexa(0)
      Phi1step(ndx)  = +cg*Tdiff*Phigrd(ndx)+Phiexa(ndx)
      Phi1step(ndx-1)= +cg*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
   end if

if(mode==256) then
Phicgp(-1)=-cg*Tdiff*Phigrd(-1)+Phiexa(-1)
Phicgp(0) =-cg*Tdiff*Phigrd(0)+Phiexa(0)
Phicgp(ndx)  = -cg*Tdiff*Phigrd(ndx)+Phiexa(ndx)
Phicgp(ndx-1)= -cg*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)

Phicgm(-1)=+cg*Tdiff*Phigrd(-1)+Phiexa(-1)
Phicgm(0) =+cg*Tdiff*Phigrd(0)+Phiexa(0)
Phicgm(ndx)  = +cg*Tdiff*Phigrd(ndx)+Phiexa(ndx)
Phicgm(ndx-1)= +cg*Tdiff*Phigrd(ndx-1)+Phiexa(ndx-1)
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
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,11)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,11)
     !------------calcurate dt/2------------
     do i=ist-1,ndx-ien+1
        Phi2dt(i) = Phipre(i) + 0.5d0 * nu2 * ( Phiu(i+1) - Phiu(i))
     end do
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4)
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
     !write(*,*) 'in1'
     do i=ist,ndx-ien
        Phiv(i) =  -cg * G4pi * source(i) * dt + Phipre(i)
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
  integer :: i,in1,j
  character(5) name
  double precision errphi1,errphi2

  write(name,'(i5.5)') in1
  open(20,file = dir//'phi'//name//'.dat')
  open(21,file = dir//'phi'//name//'.csv')
  open(22,file = dir//'phi'//name//'.txt')
  write(21,*) 'x,','Phicgp,','Phicgm,','Phi1step,','Phi2step,','rho','Phiexa','Phigrd'
  do i=-1,ndx
    write(21,*) x(i),',',Phicgp(i),',',Phicgm(i),',',Phi1step(i),',',Phi2step(i),',',rho(i),',',Phiexa(i),',',phigrd(i)
    write(20,*) x(i),',',Phicgp(i),',',Phicgm(i),',',Phi1step(i),',',Phi2step(i),',',rho(i),',',Phiexa(i),',',phigrd(i)
    write(22,*) x(i),',',Phicgp(i),',',Phicgm(i),',',Phi1step(i),',',Phi2step(i),',',rho(i),',',Phiexa(i),',',phigrd(i)
  end do
  close(21)
  close(20)
  close(22)
  write(*,*) 'save step : ',in1


errphi1=0.d0
errphi2=0.d0

do i=1,ndx-2
   !errphi1=(Phiexa(i,j)+meanphi1-meanphiexa)**2
   !errphi2=(Phiexa(i,j)+meanphi2-meanphiexa)**2
   errphi1=(Phiexa(i)-Phicgp(i))**2.d0+errphi1
   errphi2=(Phiexa(i)-(0.5d0*(Phicgp(i)+Phicgm(i))))**2.d0+errphi2
end do
errphi1=dsqrt(errphi1/dble(ndx-2))
errphi2=dsqrt(errphi2/dble(ndx-2))

open(850,file=dir//'err-phi-exa2.dat',FORM='FORMATTED',position='append')
write(850,*) in1*svnum,',', errphi1,',', errphi2
close(850)


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
     !write(*,*) slop(127),'127slop'
     !u(:)=ul(:)
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
