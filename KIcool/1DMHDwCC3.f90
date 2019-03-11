!*******************************************************!
!      1D MHD with Cooling, Chemistry & Conduction      !
!                              programed by T. Inoue    !
!*******************************************************!

MODULE comvar
INTEGER, parameter :: nd=1024+2
double precision, dimension(-1:nd) :: x,dx,Va
double precision, dimension(-1:nd,8) :: U

double precision  :: gamma,gammi1,gammi2,gammi3,gampl1,gampl2,gampl3
double precision  :: CFL,facdep,tfinal,time
double precision  :: pmin,pmax,rmin,rmax
INTEGER :: Ncell,maxstp,nitera
INTEGER :: iflag,ifchem,ifthrm,ifrad
INTEGER :: BCx1,BCx2

double precision, parameter :: kb=8.63359d0, Kcond=1.6384d-2
double precision, parameter :: mH=1.d0, mHe=4.d0, mH2=2.d0, mC=12.d0, mCO=28.d0
double precision, parameter :: G0=1.d0, xc=1.4d-4, xo=3.2d-4, dv=2.d0, Tgr=5.d-3
double precision, dimension(-1:nd) :: ndp,ndH,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
double precision, dimension(-1:nd,2) :: Ntot,NH2,NnC,NCO,tCII
double precision  :: ndpmin,ndHmin,ndH2min,ndHemin,ndHepmin,ndCmin,ndCpmin,ndCOmin

CHARACTER(33) :: dir='/Users/maeda/Desktop/code/KIcool/'

END MODULE comvar


!======================================================================*
!                                 MAIN                                 *
!======================================================================*

PROGRAM MAIN

call INITIA
call EVOLVE

END PROGRAM MAIN

!======================================================================*
!                     Prepare an Initial State                         *
!======================================================================*
SUBROUTINE INITIA
USE comvar
double precision ::  ql1x,ql2x,dinit1,dinit2,pinit1,pinit2, &
           vinitx1,vinitx2,vinity1,vinity2,vinitz1,vinitz2, &
           binitx1,binitx2,binity1,binity2,binitz1,binitz2, &
           Hini1,Hini2,pini1,pini2,H2ini1,H2ini2,Heini1,Heini2, &
           Hepini1,Hepini2,Cini1,Cini2,COini1,COini2,Cpini1,Cpini2
character(3) :: NPENUM


double precision :: pi,kwave,xp,kz,Bz,vzinp
integer :: bzstrt , bzfin

open(8,file=dir//'INPUT.DAT')
  read(8,*)  Np1x,Np2x
  read(8,*)  ql1x,ql2x
  read(8,*)  dinit1,dinit2
  read(8,*)  vinitx1,vinitx2
  read(8,*)  vinity1,vinity2
  read(8,*)  vinitz1,vinitz2
  read(8,*)  pinit1,pinit2
  read(8,*)  binitx1,binitx2
  read(8,*)  binity1,binity2
  read(8,*)  binitz1,binitz2
  read(8,*)  Hini1,Hini2
  read(8,*)  pini1,pini2
  read(8,*)  H2ini1,H2ini2
  read(8,*)  Heini1,Heini2
  read(8,*)  Hepini1,Hepini2
  read(8,*)  Cini1,Cini2
  read(8,*)  COini1,COini2
  read(8,*)  Cpini1,Cpini2
  read(8,*)  CFL,facdep
  read(8,*)  maxstp,nitera,tfinal
  read(8,*)  BCx1,BCx2
  read(8,*)  iflag,ifchem,ifthrm,ifrad
close(8)

  pinit1=8.810807d3*kb*1.d-3
  pinit2=8.810807d3*kb*1.d-3
!UEQ ntot = 5.059326
goto 10002
 pinit1=5.233297d3*kb*1.d-3; pinit2=pinit1 !セミコロンは文の区切り。
 Hini1=0.4632418d1; pini1=0.9158390d-2; H2ini1=0.6125440d-5; Heini1=0.4171380d0; Hepini1=0.6050313d-3
 Cini1=0.2951134d-7; COini1=0.1783309d-15; Cpini1=0.7082769d-3
 Hini2=Hini1; pini2=pini1; H2ini2=H2ini1; Heini2=Heini1; Hepini2=Hepini1; Cini2=Cini1; COini2=COini1; Cpini2=Cpini1
 dinit1=mH*Hini1+mH*pini1+mH2*H2ini1+mHe*Heini1+mHe*Hepini1; dinit2=dinit1
10002 continue

gamma  = ( 5.d0*(Hini1+pini1+Heini1+Hepini1)+7.d0*H2ini1 )/( 3.d0*(Hini1+pini1+Heini1+Hepini1)+5.d0*H2ini1 )
gammi1 = gamma - 1.0d0; gammi2 = gamma - 2.0d0; gammi3 = gamma - 3.0d0
gampl1 = gamma + 1.0d0; gampl2 = gamma + 2.0d0; gampl3 = gamma + 3.0d0

pmin = 1.d-20
pmax = 1.d20
rmin = 1.d-20
rmax = 1.d20

ndHmin  = rmin*0.91d0; ndpmin  = 1.d-20; ndH2min = 1.d-20; ndHemin = rmin*0.09d0
ndHepmin= 1.d-20; ndCpmin = 1.d-20; ndCmin = 1.d-20; ndCOmin = 1.d-20


!***** shock tube test *****!

Ncell = Np1x + Np2x; if(Ncell.ne.nd-2) write(*,*) 'err at INPUT'

do i = -1, Ncell+2
  if(i.le.Np1x) then
    U(i,1) = dinit1
    U(i,2) = vinitx1
    U(i,3) = vinity1
    U(i,4) = vinitz1
    U(i,5) = pinit1
    U(i,6) = binitx1
    U(i,7) = binity1
    U(i,8) = binitz1
    ndH(i)   = Hini1
    ndp(i)   = pini1
    ndH2(i)  = H2ini1
    ndHe(i)  = Heini1
    ndHep(i) = Hepini1
    ndC(i)   = Cini1
    ndCO(i)  = COini1
    ndCp(i)  = Cpini1
  end if
  if(i.gt.Np1x) then
    U(i,1) = dinit2
    U(i,2) = vinitx2
    U(i,3) = vinity2
    U(i,4) = vinitz2
    U(i,5) = pinit2
    U(i,6) = binitx2
    U(i,7) = binity2
    U(i,8) = binitz2
    ndH(i)   = Hini2
    ndp(i)   = pini2
    ndH2(i)  = H2ini2
    ndHe(i)  = Heini2
    ndHep(i) = Hepini2
    ndC(i)   = Cini2
    ndCO(i)  = COini2
    ndCp(i)  = Cpini2
  end if
end do


! TEST ALFVEN
!-------------------------------------------
goto 2003
!pi = 3.1415d0
!kz = 1.0d0
!Bz = 0.01d0
bzstrt = 100
bzfin = 200
vzinp = 0.001d0

!do i = bzstrt , bzfin
 !  U(i,8) = Bz*dcos(kz * (dble(i-bzstrt)/dble(bzfin-bzstrt)) * pi)
!end do
do i = bzstrt , bzfin
   U(i,4) = vzinp
end do
2003 continue
!TEST ALFVEN
!-------------------------------------------


!write(*,*) U(-1,1),U(-1,2),U(-1,3),U(-1,4),U(-1,5),U(-1,6),U(-1,7),U(-1,8)
!write(*,*) U(0,1),U(0,2),U(0,3),U(0,4),U(0,5),U(0,6),U(0,7),U(0,8)
!write(*,*) U(Ncell+1,1),U(Ncell+1,2),U(Ncell+1,3),U(Ncell+1,4),U(Ncell+1,5),U(Ncell+1,6),U(Ncell+1,7),U(Ncell+1,8)
!write(*,*) U(Ncell+2,1),U(Ncell+2,2),U(Ncell+2,3),U(Ncell+2,4),U(Ncell+2,5),U(Ncell+2,6),U(Ncell+2,7),U(Ncell+2,8)





do i = -1, Ncell+2
  dx(i) = ql1x/dble(Np1x)
end do

x(-1) = -dx(0) !負
do i = 0, Ncell+2
   x(i) = x(i-1) + dx(i) !位置座標の設定
end do

!***** TItest *****!
goto 101
pi = 3.1415d0
kwave = 1.d0
do i=1,Ncell
  xp = 0.5d0*(x(i)+x(i-1))
  ndH(i) = ndH(i) + 0.05d0*Hini1*dsin(2.d0*pi*kwave*xp)
  U(i,1) = mH*ndH(i)+mH*ndp(i)+mH2*ndH2(i)+mHe*ndHe(i)+mHe*ndHep(i)
end do
101 continue

!***** Read Initial Conditions *****!
open(2,file=dir//'tsave.DAT')
  read(2,'(d25.17)') amp
  read(2,'(i8)') nunit
close(2)

if(nunit.eq.1) goto 118
  open(unit=8,file=dir//'000.dat',FORM='UNFORMATTED') !CONVERT='LITTLE_ENDIAN'
    read(8) (x(i),U(i,1),U(i,2),U(i,3),U(i,4),U(i,5),U(i,6),U(i,7),U(i,8), &
             ndH(i),ndp(i),ndH2(i),ndHe(i),ndHep(i),ndC(i),ndCO(i),ndCp(i),i=-1,Ncell+2)
  close(8)
  do i=1,Ncell; dx(i) = x(i)-x(i-1); end do
  call BC(Ncell,dx(-1),BCx1,BCx2)
118  continue


do i=1,Ncell+1
  nde(i) = ndp(i)+ndHep(i)+ndCp(i)
  ndtot(i) = ndp(i)+ndH(i)+2.d0*ndH2(i)+ndHe(i)+ndHep(i)
  Ntot(i,1)=0.d0; NH2(i,1)=0.d0; NnC(i,1)=0.d0; NCO(i,1)=0.d0; tCII(i,1)=0.d0
  Ntot(i,2)=0.d0; NH2(i,2)=0.d0; NnC(i,2)=0.d0; NCO(i,2)=0.d0; tCII(i,2)=0.d0
end do

if(ifrad.eq.2) then; do l=1,20; call SHIELD(); end do; end if

END SUBROUTINE INITIA


!=====================================================================*
!                  Integration of The Evolution                       *
!=====================================================================*

SUBROUTINE EVOLVE
USE comvar
double precision  :: t(10000),dt, tLMT
character(3) fnunit


open(2,file=dir//'tsave.DAT')
  read(2,'(d25.17)') time
  read(2,'(i8)') nunit
close(2)
open(3,file=dir//'time.DAT')
  do i = 1, nunit
  read(3,'(d25.17)') t(i)
  end do
close(3)

open(100,file='cooling.dat')
write(100,*) time,U(510,1)/1.27d0,U(510,5)*1.52d0/1.38d0 * 1.d2,U(510,6)
do in10 = 1, maxstp
  if(time.ge.tfinal) goto 9000
  call SAVEU(nunit,dt,t,0)

  do in20 = 1, nitera
    if(time.ge.tfinal) goto 9000
!***** Determine time-step dt *****
    dt = tfinal
    call Couran(tLMT)
    dt = dmin1( dt, CFL * tLMT )
    
    if(ifthrm.eq.2) then
      call Stblty(tLMT)
      dt = dmin1( dt, 0.2d0 * tLMT )
    end if
    
    if(mod(in20,10).eq.1) write(*,*) in20,time,dt
    if(time+dt.gt.tfinal) dt = tfinal - time

!***** Source parts *****
    call SOURCE(dt*0.5d0)
!***** Godunov parts *****
    call MHD(dt)
!***** Source parts *****
    call SOURCE(dt*0.5d0)
!************************
    time = time + dt
    write(100,*) time,U(512,1)/1.27d0,U(512,5)*1.52d0/1.38d0 * 1.d2,U(510,7)
  end do
end do

9000 continue
call SAVEU(nunit,dt,t,1)
close(100)
END SUBROUTINE EVOLVE


!----------------------------------------------------------- SAVE VARIABLES ---!

SUBROUTINE SAVEU(nunit,dt,t,msig)
USE comvar

integer :: nunit,msig
double precision  :: dt,t(10000)
character(5) filenm

write(filenm,'("bw",I3.3)') nunit
100 format(E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3, &
         E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3,E19.10e3) !再利用のため100をつけておく
!100 format(E21.12,E21.12,E21.12,E21.12,E21.12,E21.12,E21.12,E21.12,E21.12, &
 !E21.12,E21.12,E21.12,E21.12,E21.12,E21.12,E21.12,E21.12)

open(10,file=dir//filenm//'.dat')
   write(10,100) ( 0.5d0*(x(i)+x(i-1)),U(i,1),U(i,2),U(i,3),U(i,4),U(i,5),U(i,6),U(i,7),U(i,8), &
                   ndH(i),ndp(i),ndH2(i),ndHe(i),ndHep(i),ndC(i),ndCO(i),ndCp(i),i=1,Ncell )
close(10)

!write(*,*) U(-1,1),U(-1,2),U(-1,3),U(-1,4),U(-1,5),U(-1,6),U(-1,7),U(-1,8)
!write(*,*) U(0,1),U(0,2),U(0,3),U(0,4),U(0,5),U(0,6),U(0,7),U(0,8)
!write(*,*) U(Ncell+1,1),U(Ncell+1,2),U(Ncell+1,3),U(Ncell+1,4),U(Ncell+1,5),U(Ncell+1,6),U(Ncell+1,7),U(Ncell+1,8)
!write(*,*) U(Ncell+2,1),U(Ncell+2,2),U(Ncell+2,3),U(Ncell+2,4),U(Ncell+2,5),U(Ncell+2,6),U(Ncell+2,7),U(Ncell+2,8)

!if(nunit >= 70) then

! write(*,*) '*************************************'
! do i= 420 , 450
!    write(*,*) '---------------------------'
!    write(*,*) i
! write(*,*) x(i),U(i,1),U(i,2),U(i,3),U(i,4),U(i,5),U(i,6),U(i,7),U(i,8)
! write(*,*) '---------------------------'
!enddo
! write(*,*) '*************************************'
!endif
!write(*,*) U(512,5), nunit

t(nunit) = time
open(3,file=dir//'time.DAT')
  do i = 1, nunit
    write(3,'(d25.17)') t(i)
  end do
close(3)

!write(*,*) '273'

if(msig.eq.1) then
!   write(*,*) 'time = ',time
  open(2,file=dir//'tsave.DAT')
    write(2,'(d25.17)') time
    write(2,'(i8)') nunit
  close(2)

  open(10,file=dir//'000.dat',FORM='UNFORMATTED') !CONVERT='LITTLE_ENDIAN'
    write(10) ( x(i),U(i,1),U(i,2),U(i,3),U(i,4),U(i,5),U(i,6),U(i,7),U(i,8), &
                ndH(i),ndp(i),ndH2(i),ndHe(i),ndHep(i),ndC(i),ndCO(i),ndCp(i),i=-1,Ncell+2 )
  close(10)
end if

nunit = nunit + 1

END SUBROUTINE SAVEU


!=====================================================================*
!           Boundary conditions                                       *
!=====================================================================*
SUBROUTINE BC(Ncell,Q,BC1,BC2)
integer :: BC1,BC2
double precision  :: Q(-1:Ncell+2)

if(BC1.eq.1) then !Free Boundary Condition
  Q(0 ) = Q(1)
  Q(-1) = Q(1)
end if
if(BC2.eq.1) then
  Q(Ncell+1) = Q(Ncell)
  Q(Ncell+2) = Q(Ncell)
end if

if(BC1.eq.3) then !Periodic Boundary Condition
  Q(0 ) = Q(Ncell)
  Q(-1) = Q(Ncell-1)
end if
if(BC2.eq.3) then
  Q(Ncell+1) = Q(1)
  Q(Ncell+2) = Q(2)
end if


END SUBROUTINE BC


  

!*--------------------------------------------------------------*-
  !f77 -*!
!*                 3D ideal MHD
  !*!
!*                 Lagrange and Eulerian remapping scheme
  !*!
!*
  !*!
!*                 E = P/(gamma-1) + rho*v^2/2 + B^2/2
  !*!
!*
  !*!
!*======================================================================*!
SUBROUTINE MHD(dt)
USE comvar
      
integer :: i
double precision  :: dt,xlag(-1:nd), ekin,emag, invd

do i = -1, Ncell+2
  xlag(i) = x(i)
end do

call RHS(dt,xlag)
call MOC(dt,xlag)

do i = 1, Ncell
  U(i,2) = U(i,2) * U(i,1); U(i,3) = U(i,3) * U(i,1)
  U(i,4) = U(i,4) * U(i,1); U(i,5) = U(i,5) * U(i,1)
end do

if(iflag.eq.2) then
  do i=0,Ncell
    x(i) = xlag(i)
  end do
  do i=1,Ncell
    dx(i) = x(i)-x(i-1)
  end do
  call BC(Ncell,dx(-1),BCx1,BCx2)
  x(-1)      = x(0) -dx(0)
 ! x(-2)      = x(-1)-dx(-1) 要確認
  x(Ncell+1) = x(Ncell)   + dx(Ncell+1)
  x(Ncell+2) = x(Ncell+1) + dx(Ncell+2)
else
  call REMAP(dt,xlag)
end if

do i = 1, Ncell
  U(i,1) = dmax1( 1.27d0*rmin, U(i,1) )
  invd = 1.d0/U(i,1)
  U(i,2) = U(i,2)*invd
  U(i,3) = U(i,3)*invd
  U(i,4) = U(i,4)*invd
  gammi1 =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
  gammi1 = ( 2.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+2.d0*ndH2(i) )/gammi1
  ekin = 0.5d0 * ( U(i,2)**2 + U(i,3)**2 + U(i,4)**2 ) * U(i,1)
  emag = 0.5d0 * ( U(i,6)**2 + U(i,7)**2 + U(i,8)**2 )
  U(i,5) = (  U(i,5)-ekin-emag  )*gammi1
  U(i,5) = dmax1( pmin, U(i,5) )
  U(i,5) = dmin1( pmax, U(i,5) )

    ndH(i) = dmax1( rmin,ndH(i)   )
    ndp(i) = dmax1( ndpmin,ndp(i)   )
   ndH2(i) = dmax1( ndH2min ,ndH2(i)  )
   ndHe(i) = dmax1( rmin*0.09d0,ndHe(i)  )
  ndHep(i) = dmax1( ndHepmin,ndHep(i) )
    ndC(i) = dmax1( ndCmin  ,ndC(i)   )
   ndCp(i) = dmax1( ndCpmin ,ndCp(i)  )
   ndCO(i) = dmax1( ndCOmin ,ndCO(i)  )
end do

END SUBROUTINE MHD


!======================================================================*
!           Right Hand Sides of Euler's Equations
!======================================================================*
SUBROUTINE RHS(dt,xlag)
USE comvar

double precision  :: dt,grdQ(-1:nd,8),xlag(-1:nd)
double precision  :: QL1,QL2,QL3,QL4,QR1,QR2,QR3,QR4
double precision  :: Pa(0:nd)
double precision  :: depend1,depend2,cm
double precision  :: dm(-1:nd)
double precision  :: ndHm,ndpm,ndHem,ndHepm,ndH2m

CALL BC(Ncell,U(-1,1),BCx1,BCx2); call VLIMIT(Ncell,U(-1,1),grdQ(-1,1),dx,0,1)
CALL BC(Ncell,U(-1,2),BCx1,BCx2); call VLIMIT(Ncell,U(-1,2),grdQ(-1,2),dx,0,1)
CALL BC(Ncell,U(-1,5),BCx1,BCx2); call VLIMIT(Ncell,U(-1,5),grdQ(-1,5),dx,0,1)
CALL BC(Ncell,U(-1,7),BCx1,BCx2); call VLIMIT(Ncell,U(-1,7),grdQ(-1,7),dx,0,1)
CALL BC(Ncell,U(-1,8),BCx1,BCx2); call VLIMIT(Ncell,U(-1,8),grdQ(-1,8),dx,0,1)
CALL BC(Ncell,ndH(-1),BCx1,BCx2)
CALL BC(Ncell,ndp(-1),BCx1,BCx2)
CALL BC(Ncell,ndH2(-1),BCx1,BCx2)
CALL BC(Ncell,ndHe(-1),BCx1,BCx2)
CALL BC(Ncell,ndHep(-1),BCx1,BCx2)


do i = 0, Ncell

  ix  = i
  ixp = i+1

  gamma =   3.d0*(ndH(ix )+ndp(ix )+ndHe(ix )+ndHep(ix ))+5.d0*ndH2(ix )
  gamma = ( 5.d0*(ndH(ix )+ndp(ix )+ndHe(ix )+ndHep(ix ))+7.d0*ndH2(ix ) )/gamma
  cm = dsqrt(  (gamma * U(ix ,5) + U(ix ,7)**2 + U(ix ,8)**2) / U(ix ,1 )  )
  depend1 = 0.5d0*(dx(i  ) - cm * dt)*facdep
  depend1 = dmax1(0.d0, depend1)
  gamma =   3.d0*(ndH(ixp)+ndp(ixp)+ndHe(ixp)+ndHep(ixp))+5.d0*ndH2(ixp)
  gamma = ( 5.d0*(ndH(ixp)+ndp(ixp)+ndHe(ixp)+ndHep(ixp))+7.d0*ndH2(ixp) )/gamma
  cm = dsqrt(  (gamma * U(ixp,5) + U(ixp,7)**2 + U(ixp,8)**2) / U(ixp,1)  )
  depend2 = 0.5d0*(dx(i+1) - cm * dt)*facdep
  depend2 = dmax1(0.d0, depend2)

  QL1 = U(ix ,1) + depend1 * grdQ(i  ,1)
  QL1 = (0.5d0-dsign(0.5d0,-QL1))*QL1 + (0.5d0+dsign(0.5d0,-QL1))*U(ix ,1)
  QR1 = U(ixp,1) - depend2 * grdQ(i+1,1)
  QR1 = (0.5d0-dsign(0.5d0,-QR1))*QR1 + (0.5d0+dsign(0.5d0,-QR1))*U(ixp,1)
  QL2 = U(ix ,2) + depend1 * grdQ(i  ,2)
  QR2 = U(ixp,2) - depend2 * grdQ(i+1,2)
  QL3 = U(ix ,5) + depend1 * grdQ(i  ,5)
  QL3 = (0.5d0-dsign(0.5d0,-QL3))*QL3 + (0.5d0+dsign(0.5d0,-QL3))*U(ix ,5)
  QR3 = U(ixp,5) - depend2 * grdQ(i+1,5)
  QR3 = (0.5d0-dsign(0.5d0,-QR3))*QR3 + (0.5d0+dsign(0.5d0,-QR3))*U(ixp,5)
  QL4 = (U(ix ,7)+depend1*grdQ(i  ,7))**2 + (U(ix ,8)+depend1*grdQ(i  ,8))**2
  QR4 = (U(ixp,7)-depend2*grdQ(i+1,7))**2 + (U(ixp,8)-depend2*grdQ(i+1,8))**2

  ndHm =0.5d0*( ndH(ix)+ ndH(ixp)); ndpm =0.5d0*(  ndp(ix)+  ndp(ixp))
  ndHem=0.5d0*(ndHe(ix)+ndHe(ixp));ndHepm=0.5d0*(ndHep(ix)+ndHep(ixp))
  ndH2m=0.5d0*(ndH2(ix)+ndH2(ixp))
  gamma = ( 5.d0*(ndHm+ndpm+ndHem+ndHepm)+7.d0*ndH2m )/( 3.d0*(ndHm+ndpm+ndHem+ndHepm)+5.d0*ndH2m )
  gammi1 = gamma - 1.0d0
  gammi2 = gamma - 2.0d0
  gammi3 = gamma - 3.0d0
  gampl1 = gamma + 1.0d0
  gampl2 = gamma + 2.0d0
  gampl3 = gamma + 3.0d0
!---------------------------------------------------------*
  call RIEMAN( QL1,QL2,QL3,QL4,QR1,QR2,QR3,QR4,Pa(i),Va(i), 5 )
!---------------------------------------------------------*

end do

!***** store lagrangian mass *****
do i = 0, Ncell
  dm(i)   = dx(i) * U(i,1)
  xlag(i) = xlag(i) + Va(i) * dt
end do
!***** Conservation laws *****

do i = 1, Ncell
  ix  = i
  ixm = i-1

  gammi1 =   3.d0*(ndH(ix)+ndp(ix)+ndHe(ix)+ndHep(ix))+5.d0*ndH2(ix)
  gammi1 = ( 2.d0*(ndH(ix)+ndp(ix)+ndHe(ix)+ndHep(ix))+2.d0*ndH2(ix) )/gammi1
  U(ix,5 )  = U(ix,5)/( U(ix,1)*gammi1 ) &
            + 0.5d0*( U(ix,2)**2.d0+U(ix,3)**2.d0+U(ix,4)**2.d0 ) &
            + 0.5d0*( U(ix,6)**2.d0+U(ix,7)**2.d0+U(ix,8)**2.d0 )/U(ix,1)
  U(ix,1)  = dm(i) / ( xlag(ix) - xlag(ixm) )
  U(ix,2)  = U(ix,2) - dt / dm(i) * ( Pa(i) - Pa(i-1) )
  U(ix,5)  = U(ix,5) - dt / dm(i) * ( Pa(i)*Va(i)-Pa(i-1)*Va(i-1) )
  U(ix,7) = U(ix,7)*dx(i) / ( xlag(ix) - xlag(ixm) )
  U(ix,8) = U(ix,8)*dx(i) / ( xlag(ix) - xlag(ixm) )

    ndH(ix) = dx(i)*  ndH(ix)/( xlag(ix) - xlag(ixm) )
    ndp(ix) = dx(i)*  ndp(ix)/( xlag(ix) - xlag(ixm) )
   ndH2(ix) = dx(i)* ndH2(ix)/( xlag(ix) - xlag(ixm) )
   ndHe(ix) = dx(i)* ndHe(ix)/( xlag(ix) - xlag(ixm) )
  ndHep(ix) = dx(i)*ndHep(ix)/( xlag(ix) - xlag(ixm) )
    ndC(ix) = dx(i)*  ndC(ix)/( xlag(ix) - xlag(ixm) )
   ndCp(ix) = dx(i)* ndCp(ix)/( xlag(ix) - xlag(ixm) )
   ndCO(ix) = dx(i)* ndCO(ix)/( xlag(ix) - xlag(ixm) )
end do

END SUBROUTINE RHS


!*======================================================================*
!*                 The Riemann Solver for Adiabatic Gases               *
!*                    with Tangential Magnetic Fields                   *
!*======================================================================*
SUBROUTINE RIEMAN( uLi1,uLi2,uLi3,uLi4,uRi1,uRi2,uRi3,uRi4,Pai,Vai, itrn )
USE comvar
integer :: itrn
double precision  :: Pai,Vai,uLi1,uLi2,uLi3,uLi4,uRi1,uRi2,uRi3,uRi4
double precision  :: D1i,D2i,V1i,V2i,P1i,P2i,B1i,B2i,GP1i,GP2i,GP1B1i,GP2B2i,pp1i,pp2i,qM1i,qM2i,qM1sqi,qM2sqi,PaOLDi

D1i = uLi1
V1i = uLi2
P1i = uLi3
B1i = 0.5d0 * uLi4
D2i = uRi1
V2i = uRi2
P2i = uRi3
B2i = 0.5d0 * uRi4

!*****    B1 = 0.5 * B1^2     B2 = 0.5 * B2^2

P1i    = P1i + B1i
P2i    = P2i + B2i
GP1i   = GAMmi1*P1i
GP2i   = GAMmi1*P2i
GP1B1i = GAMmi3*P1i - GAMmi2*B1i*2.d0
GP2B2i = GAMmi3*P2i - GAMmi2*B2i*2.d0

Pai  = ( P1i + P2i )*0.5d0

do loop = 1, itrn

  pp1i = GAMpl3*Pai + GP1B1i 
  qM1sqi = pp1i+dsqrt( pp1i*pp1i-8.d0*(Pai-P1i)*(GAMpl1*Pai+GP1i) )
  qM1i   = dsqrt( qM1sqi*D1i )*0.5d0
  pp2i = GAMpl3*Pai + GP2B2i 
  qM2sqi = pp2i+dsqrt( pp2i*pp2i-8.d0*(Pai-P2i)*(GAMpl1*Pai+GP2i) )
  qM2i   = dsqrt( qM2sqi*D2i )*0.5d0
  PaOLDi = Pai
  Pai    = qM2i*P1i + qM1i*P2i + qM2i*qM1i*(V1i-V2i)
  Pai    = Pai/( qM2i + qM1i )
  if(Pai.le.0.0) Pai = 0.1d0*PaOLDi

end do

Vai = qM1i*V1i + qM2i*V2i + P1i - P2i
Vai = Vai / ( qM2i + qM1i )

END SUBROUTINE RIEMAN

!=====================================================================*
!                          Eulerian Remaping
!=====================================================================*
SUBROUTINE REMAP(dt,xlag)
USE comvar

double precision  :: dt,xlag(-1:nd)
double precision  :: dxlag(-1:nd), F(0:nd,8), grdU(-1:nd,8)
double precision  :: depend, vadt, dxi
double precision  :: grdC(-1:nd,8),G(0:nd,8)


do i = 1,Ncell
  dxlag(i) = xlag(i)-xlag(i-1)
end do
CALL BC(Ncell,dxlag(-1),BCx1,BCx2)
  
CALL BC(Ncell,U(-1,1),BCx1,BCx2); call VLIMIT(Ncell,U(-1,1),grdU(-1,1),dxlag,0,1)
CALL BC(Ncell,U(-1,2),BCx1,BCx2); call VLIMIT(Ncell,U(-1,2),grdU(-1,2),dxlag,0,1)
CALL BC(Ncell,U(-1,3),BCx1,BCx2); call VLIMIT(Ncell,U(-1,3),grdU(-1,3),dxlag,0,1)
CALL BC(Ncell,U(-1,4),BCx1,BCx2); call VLIMIT(Ncell,U(-1,4),grdU(-1,4),dxlag,0,1)
CALL BC(Ncell,U(-1,5),BCx1,BCx2); call VLIMIT(Ncell,U(-1,5),grdU(-1,5),dxlag,0,1)
CALL BC(Ncell,U(-1,7),BCx1,BCx2); call VLIMIT(Ncell,U(-1,7),grdU(-1,7),dxlag,0,1)
CALL BC(Ncell,U(-1,8),BCx1,BCx2); call VLIMIT(Ncell,U(-1,8),grdU(-1,8),dxlag,0,1)
  
CALL BC(Ncell,  ndH(-1),BCx1,BCx2); call VLIMIT(Ncell,  ndH(-1),grdC(-1,1),dxlag,0,1)
CALL BC(Ncell,  ndp(-1),BCx1,BCx2); call VLIMIT(Ncell,  ndp(-1),grdC(-1,2),dxlag,0,1)
CALL BC(Ncell, ndH2(-1),BCx1,BCx2); call VLIMIT(Ncell, ndH2(-1),grdC(-1,3),dxlag,0,1)
CALL BC(Ncell, ndHe(-1),BCx1,BCx2); call VLIMIT(Ncell, ndHe(-1),grdC(-1,4),dxlag,0,1)
CALL BC(Ncell,ndHep(-1),BCx1,BCx2); call VLIMIT(Ncell,ndHep(-1),grdC(-1,5),dxlag,0,1)
CALL BC(Ncell,  ndC(-1),BCx1,BCx2); call VLIMIT(Ncell,  ndC(-1),grdC(-1,6),dxlag,0,1)
CALL BC(Ncell, ndCp(-1),BCx1,BCx2); call VLIMIT(Ncell, ndCp(-1),grdC(-1,7),dxlag,0,1)
CALL BC(Ncell, ndCO(-1),BCx1,BCx2); call VLIMIT(Ncell, ndCO(-1),grdC(-1,8),dxlag,0,1)


do i = 0, Ncell
  ix  = i
  ixp = i+1
  vadt= Va(ix)*dt
  if(vadt.le.0.d0) then
    depend = 0.5d0 * (dxlag(i+1) + vadt) * facdep
    F(i,1) = -vadt * ( U(ixp,1) - depend*grdU(i+1,1) )
    F(i,2) = -vadt * ( U(ixp,2) - depend*grdU(i+1,2) )
    F(i,3) = -vadt * ( U(ixp,3) - depend*grdU(i+1,3) )
    F(i,4) = -vadt * ( U(ixp,4) - depend*grdU(i+1,4) )
    F(i,5) = -vadt * ( U(ixp,5) - depend*grdU(i+1,5) )
    F(i,7) = -vadt * ( U(ixp,7) - depend*grdU(i+1,7) )
    F(i,8) = -vadt * ( U(ixp,8) - depend*grdU(i+1,8) )
    G(i,1) = -vadt * ( ndH(ixp) - depend*grdC(i+1,1) )
    G(i,2) = -vadt * ( ndp(ixp) - depend*grdC(i+1,2) )
    G(i,3) = -vadt * (ndH2(ixp) - depend*grdC(i+1,3) )
    G(i,4) = -vadt * (ndHe(ixp) - depend*grdC(i+1,4) )
    G(i,5) = -vadt * (ndHep(ixp)- depend*grdC(i+1,5) )
    G(i,6) = -vadt * ( ndC(ixp) - depend*grdC(i+1,6) )
    G(i,7) = -vadt * (ndCp(ixp) - depend*grdC(i+1,7) )
    G(i,8) = -vadt * (ndCO(ixp) - depend*grdC(i+1,8) )
  else  !*** va > 0
    depend = 0.5d0 * (dxlag(i  ) - vadt) * facdep
    F(i,1) = -vadt * ( U(ix ,1) + depend*grdU(i  ,1) )
    F(i,2) = -vadt * ( U(ix ,2) + depend*grdU(i  ,2) )
    F(i,3) = -vadt * ( U(ix ,3) + depend*grdU(i  ,3) )
    F(i,4) = -vadt * ( U(ix ,4) + depend*grdU(i  ,4) )
    F(i,5) = -vadt * ( U(ix ,5) + depend*grdU(i  ,5) )
    F(i,7) = -vadt * ( U(ix ,7) + depend*grdU(i  ,7) )
    F(i,8) = -vadt * ( U(ix ,8) + depend*grdU(i  ,8) )
    G(i,1) = -vadt * ( ndH(ix ) + depend*grdC(i  ,1) )
    G(i,2) = -vadt * ( ndp(ix ) + depend*grdC(i  ,2) )
    G(i,3) = -vadt * (ndH2(ix ) + depend*grdC(i  ,3) )
    G(i,4) = -vadt * (ndHe(ix ) + depend*grdC(i  ,4) )
    G(i,5) = -vadt * (ndHep(ix )+ depend*grdC(i  ,5) )
    G(i,6) = -vadt * ( ndC(ix ) + depend*grdC(i  ,6) )
    G(i,7) = -vadt * (ndCp(ix ) + depend*grdC(i  ,7) )
    G(i,8) = -vadt * (ndCO(ix ) + depend*grdC(i  ,8) )
  end if
end do

do i = 1, Ncell
  dxi    = 1.d0/(x(i)-x(i-1))
  U(i,1) = ( U(i,1)*dxlag(i) + F(i,1) - F(i-1,1) )*dxi
  U(i,2) = ( U(i,2)*dxlag(i) + F(i,2) - F(i-1,2) )*dxi
  U(i,3) = ( U(i,3)*dxlag(i) + F(i,3) - F(i-1,3) )*dxi
  U(i,4) = ( U(i,4)*dxlag(i) + F(i,4) - F(i-1,4) )*dxi
  U(i,5) = ( U(i,5)*dxlag(i) + F(i,5) - F(i-1,5) )*dxi
  U(i,7) = ( U(i,7)*dxlag(i) + F(i,7) - F(i-1,7) )*dxi
  U(i,8) = ( U(i,8)*dxlag(i) + F(i,8) - F(i-1,8) )*dxi
  ndH(i) = ( ndH(i)*dxlag(i) + G(i,1) - G(i-1,1) )*dxi
  ndp(i) = ( ndp(i)*dxlag(i) + G(i,2) - G(i-1,2) )*dxi
 ndH2(i) = (ndH2(i)*dxlag(i) + G(i,3) - G(i-1,3) )*dxi
 ndHe(i) = (ndHe(i)*dxlag(i) + G(i,4) - G(i-1,4) )*dxi
ndHep(i) = (ndHep(i)*dxlag(i)+ G(i,5) - G(i-1,5) )*dxi
  ndC(i) = ( ndC(i)*dxlag(i) + G(i,6) - G(i-1,6) )*dxi
 ndCp(i) = (ndCp(i)*dxlag(i) + G(i,7) - G(i-1,7) )*dxi
 ndCO(i) = (ndCO(i)*dxlag(i) + G(i,8) - G(i-1,8) )*dxi
end do
  
END SUBROUTINE REMAP


!======================================================================*
!                     TVD Limiting (Van Albada)
!======================================================================*
SUBROUTINE VLIMIT(Ncell,U,grdU,dx,i_sta,i_end)

integer :: Ncell
double precision, parameter :: eps = 1.d-10
double precision  :: U(-1:Ncell+2), grdU(-1:Ncell+2), dx(-1:Ncell+2)
double precision  :: delp,delm,flmt

do i = i_sta, Ncell+i_end
  ix  = i
  ixp = i+1
  ixm = i-1

  delp = U(ixp)-U(ix )
  delm = U(ix )-U(ixm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
  grdU(i) = flmt*( U(ixp)-U(ixm) )/( dx(i)+0.5d0*dx(i-1)+0.5d0*dx(i+1) )
end do

END SUBROUTINE


!=====================================================================*
!     Courant Condition for Euler code
!=====================================================================*
SUBROUTINE Couran(tCFL)
USE comvar
double precision :: tCFL,c2

tCFL = tfinal
do i = 1, Ncell
  gamma =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
  gamma = ( 5.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+7.d0*ndH2(i) )/gamma
  c2 = ( gamma * U(i,5) + U(i,6)**2.d0+U(i,7)**2.d0+U(i,8)**2.d0 ) / U(i,1)
  tCFL = dmin1(tCFL, dx(i)/(dsqrt(c2) + dabs(U(i,2))) )
end do
if(tCFL.lt.0.d0) write(*,*) time,NRANK,'err at Couran'

END SUBROUTINE Couran


!*-------------------------------------------------------------*- f90 -*!
!*          Solve Lorentz Force With Method of Characteristic          *!
!*                                   programed by T.Inoue              *!
!*                                        2006.11.                     *!
!*=====================================================================*!

SUBROUTINE MOC(dt,xlag)
USE comvar

double precision  :: dt, xlag(-1:nd), grdU(-1:nd,8), dxlag(-1:nd)
double precision  :: ca, depend1, depend2, deni, Bnap
double precision  :: QL11,QL12,QL13,QR11,QR12,QR13
double precision  :: QL21,QL22,QL23,QR21,QR22,QR23
double precision  :: Vay(-1:nd), Bay(-1:nd), Vaz(-1:nd), Baz(-1:nd)


do i = 1,Ncell
  dxlag(i) = xlag(i)-xlag(i-1)
end do
CALL BC(Ncell,dxlag(-1),BCx1,BCx2)
  
CALL BC(Ncell,U(-1,1),BCx1,BCx2); call VLIMIT(Ncell,U(-1,1),grdU(-1,1),dxlag,0,1)
CALL BC(Ncell,U(-1,3),BCx1,BCx2); call VLIMIT(Ncell,U(-1,3),grdU(-1,3),dxlag,0,1)
CALL BC(Ncell,U(-1,4),BCx1,BCx2); call VLIMIT(Ncell,U(-1,4),grdU(-1,4),dxlag,0,1)
CALL BC(Ncell,U(-1,7),BCx1,BCx2); call VLIMIT(Ncell,U(-1,7),grdU(-1,7),dxlag,0,1)
CALL BC(Ncell,U(-1,8),BCx1,BCx2); call VLIMIT(Ncell,U(-1,8),grdU(-1,8),dxlag,0,1)
CALL BC(Ncell,U(-1,6),BCx1,BCx2)  
!***** method of characteristic *****

do i = 0, Ncell

  ix  = i
  ixp = i+1

  Bnap = U(i,6)
  deni = dsqrt( (dxlag(i)+dxlag(i+1))/(dxlag(i)*U(ixp,1)+dxlag(i+1)*U(ix,1)) )
  ca = dabs( Bnap * deni )*0.5d0*dt
  ca = dmin1( dxlag(i  ), ca )
  depend1 = ( 0.5d0*dxlag(i  ) - ca )*facdep
  ca = dabs( Bnap * deni )*0.5d0*dt
  ca = dmin1( dxlag(i+1), ca )
  depend2 = ( 0.5d0*dxlag(i+1) - ca )*facdep
  QL11 = U(ix ,3) + depend1 * grdU(i  ,3)
  QR11 = U(ixp,3) - depend2 * grdU(i+1,3)
  QL12 = U(ix ,7) + depend1 * grdU(i  ,7)
  QR12 = U(ixp,7) - depend2 * grdU(i+1,7)
  QL13 = U(ix ,1) + depend1 * grdU(i  ,1)
  QL13 = (0.5d0-dsign(0.5d0,-QL13))*QL13 + (0.5d0+dsign(0.5d0,-QL13))*U(ix ,1)
  QR13 = U(ixp,1) - depend2 * grdU(i+1,1)
  QR13 = (0.5d0-dsign(0.5d0,-QR13))*QR13 + (0.5d0+dsign(0.5d0,-QR13))*U(ixp,1)
  QL21 = U(ix ,4) + depend1 * grdU(i  ,4)
  QR21 = U(ixp,4) - depend2 * grdU(i+1,4)
  QL22 = U(ix ,8) + depend1 * grdU(i  ,8) ! Bz cc
  QR22 = U(ixp,8) - depend2 * grdU(i+1,8)
  QL23 = QL13
  QR23 = QR13
!--------------------------------------------------------------*
  call CHAREQ( QL11,QL12,QL13,QR11,QR12,QR13,Bay(i),Vay(i) )
  call CHAREQ( QL21,QL22,QL23,QR21,QR22,QR23,Baz(i),Vaz(i) )
!--------------------------------------------------------------*
end do

do i = 1, Ncell
  ix  = i
  ixp = i+1

  U(ix,3) = U(ix,3) + dt / U(ix,1) / dxlag(i) * ( Bay(i) - Bay(i-1) ) * U(i,6)
  U(ix,4) = U(ix,4) + dt / U(ix,1) / dxlag(i) * ( Baz(i) - Baz(i-1) ) * U(i,6)
  U(ix,5) = U(ix,5) + dt / U(ix,1) / dxlag(i) * &
     (  (0.5d0*U(i,6)*Va(i  )+Bay(i  )*Vay(i  )+Baz(i  )*Vaz(i  )) * U(i,6)  &
       -(0.5d0*U(i,6)*Va(i-1)+Bay(i-1)*Vay(i-1)+Baz(i-1)*Vaz(i-1)) * U(i,6)  )
  U(ix,7) = U(ix,7) + dt / dxlag(i) * ( Vay(i) - Vay(i-1) ) * U(i,6)
  U(ix,8) = U(ix,8) + dt / dxlag(i) * ( Vaz(i) - Vaz(i-1) ) * U(i,6)
end do

END SUBROUTINE MOC


SUBROUTINE CHAREQ( QL1,QL2,QL3, QR1,QR2,QR3, Bay, Vay )
double precision :: QL1,QL2,QL3, QR1,QR2,QR3, Bay, Vay
double precision :: sqdyp, sqdym, vyp, vym, Byp, Bym

vyp   = QL1; vym   = QR1 ; Byp   = QL2; Bym   = QR2
sqdyp = dsqrt( QL3 ); sqdym = dsqrt( QR3 )

Vay = ( sqdyp*vyp + sqdym*vym + Bym - Byp )/(sqdyp+sqdym)
Bay = ( sqdyp*sqdym*(vym-vyp)+sqdyp*Bym+sqdym*Byp )/(sqdyp+sqdym)

END SUBROUTINE CHAREQ


SUBROUTINE SOURCE(dt)
USE comvar

double precision  dt
double precision :: ndpold,ndHold,ndH2old,ndHeold,ndHepold,ndCold,ndCpold,ndCOold,T
double precision :: zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec, &
          kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
double precision :: temp1,temp2,temp3,omeps,eps
double precision :: Tn(-1:nd),Pn(-1:nd),Qx(-1:nd)
double precision :: rtTx,tcd,CooL

do i = 1, Ncell
    nde(i) = ndp(i)+ndHep(i)+ndCp(i)
  ndtot(i) = ndp(i)+ndH(i)+2.d0*ndH2(i)+ndHe(i)+ndHep(i)
end do

if(ifrad.eq.2) then
  call SHIELD()
end if

!***** 2nd order explicit scheme ( with cooling function ) *****
if(ifthrm.eq.2) then

CALL BC(Ncell,U(-1,1)  ,BCx1,BCx2)
CALL BC(Ncell,U(-1,5)  ,BCx1,BCx2)
CALL BC(Ncell,ndH(-1)  ,BCx1,BCx2)
CALL BC(Ncell,ndp(-1)  ,BCx1,BCx2)
CALL BC(Ncell,ndH2(-1) ,BCx1,BCx2)
CALL BC(Ncell,ndHe(-1) ,BCx1,BCx2)
CALL BC(Ncell,ndHep(-1),BCx1,BCx2)

do i = 0, Ncell+1
  Tn(i) = U(i,5)/( kb*(ndp(i)+ndH(i)+ndH2(i)+ndHe(i)+ndHep(i)) )
  Pn(i) = U(i,5)
end do
do i = 0, Ncell
  rtTx  = dsqrt(  ( dx(i)*Tn(i+1)+dx(i+1)*Tn(i) ) /( dx(i) + dx(i+1) )  )
  Qx(i) = Kcond * rtTx * ( Tn(i+1)-Tn(i) )*2.d0 /( dx(i) + dx(i+1) )
end do
do i = 1, Ncell
!----- Cooling ---------------------------------------------------------
  Call Fcool( CooL,Tn(i),i)
  gammi1 =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
  gammi1 = ( 2.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+2.d0*ndH2(i) )/gammi1
  if( dt .le. 0.2d0*Pn(i)/(gammi1*dabs(CooL)) ) then
    U(i,5) = U(i,5) - gammi1*CooL*dt*0.5d0 !explicit
  else
    Call IMC( U(i,5),ndH(i)+ndp(i)+ndHe(i)+ndHep(i)+ndH2(i),dt*0.5d0,i ) !implicit
  end if
!----- Conduction ------------------------------------------------------
  U(i,5) = U(i,5) + gammi1*dt*0.5d0*(Qx(i)-Qx(i-1))/dx(i)
  U(i,5) = dmax1(U(i,5),pmin)
  U(i,5) = dmin1(U(i,5),pmax)
end do


CALL BC(Ncell,U(-1,5),BCx1,BCx2)

do i = 0, Ncell+1
  Tn(i) = U(i,5)/( kb*(ndp(i)+ndH(i)+ndH2(i)+ndHe(i)+ndHep(i)) )
end do
do i = 0, Ncell
  rtTx  = dsqrt(  ( dx(i)*Tn(i+1)+dx(i+1)*Tn(i) ) /( dx(i) + dx(i+1) )  )
  Qx(i) = Kcond * rtTx * ( Tn(i+1)-Tn(i) )*2.d0 /( dx(i) + dx(i+1) )
end do
do i = 1, Ncell
!----- Cooling ---------------------------------------------------------
  Call Fcool( CooL,Tn(i),i)
  gammi1 =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
  gammi1 = ( 2.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+2.d0*ndH2(i) )/gammi1
  if( dt .le. 0.2d0*Pn(i)/(gammi1*dabs(CooL)) ) then
    U(i,5) = Pn(i) - gammi1*CooL*dt !explicit
  else
    Call IMC( Pn(i),ndH(i)+ndp(i)+ndHe(i)+ndHep(i)+ndH2(i),dt,i ) !implicit
    U(i,5) = Pn(i)
  end if
!----- Conduction ------------------------------------------------------
  U(i,5) = U(i,5) + gammi1*dt*(Qx(i)-Qx(i-1))/dx(i)
  U(i,5) = dmax1(U(i,5),pmin)
  U(i,5) = dmin1(U(i,5),pmax)
end do

end if
!*******************************(C3)

!*******************************(A1)
if(ifchem.eq.2) then

do i = 1, Ncell
  ndpold=ndp(i); ndHold=ndH(i); ndH2old=ndH2(i); ndHeold=ndHe(i)
  ndHepold=ndHep(i); ndCold=ndC(i); ndCpold=ndCp(i); ndCOold=ndCO(i)
  T = U(i,5)/kb/( ndpold+ndHold+ndH2old+ndHeold+ndHepold )
  call RATES(i,T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec, &
             kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)
  
! H recombination & ionization by CR
  temp1 = kHrec*nde(i)
  temp2 = dexp(-dt*(zeta+temp1))
  temp3 = 1.d0/(zeta+temp1)
  call Omexp(omeps,dt*(zeta+temp1))
  ndH(i) = ( temp1*ndHold + zeta*ndHold*temp2 + temp1*ndpold*omeps )*temp3
  ndp(i) = (  zeta*ndpold +temp1*ndpold*temp2 +  zeta*ndHold*omeps )*temp3
  ndHold = ndH(i); ndpold = ndp(i); nde(i) = ndp(i)+ndHep(i)+ndCp(i)

!He recombination & ionization by CR
  temp1 = kHerec*nde(i)
  temp2 = dexp(-dt*(zeta+temp1))
  temp3 = 1.d0/(zeta+temp1)
  call Omexp(omeps,dt*(zeta+temp1))
  ndHe(i)  = ( temp1*ndHeold + zeta*ndHeold*temp2 + temp1*ndHepold*omeps )*temp3
  ndHep(i) = ( zeta*ndHepold +temp1*ndHepold*temp2+  zeta*ndHeold*omeps )*temp3
  ndHeold = ndHe(i); ndHepold = ndHep(i); nde(i) = ndp(i)+ndHep(i)+ndCp(i)

!H2 formation & dissociation by UV & electron collision
  temp1 = kH2ph + kH2de*nde(i)
  temp2 = temp1 + 2.d0*kH2*ndtot(i)
  temp3 = dexp(-dt*temp2)
  call Omexp(omeps,dt*temp2)
  ndH(i) = ( ndHold+2.d0*ndH2old*omeps )*temp1 + 2.d0*kH2*ndHold*ndtot(i)*temp3
  ndH(i) = ndH(i)/temp2
  ndH2(i)= ( 2.d0*ndH2old+ndHold*omeps )*kH2*ndtot(i) + temp1*ndH2old*temp3
  ndH2(i)= ndH2(i)/temp2
  ndHold = ndH(i); ndH2old = ndH2(i)
  
!H2 dissociation by H collision
  temp1 = ndHold + 2.d0*ndH2old
  temp2 = dexp(-dt*temp1*kH2dH)
  temp3 = 1.d0/(ndHold + 2.d0*ndH2old*temp2)
  ndH(i) = ndHold*temp1*temp3
  ndH2(i)= ndH2old*temp1*temp2*temp3
  ndHold = ndH(i); ndH2old = ndH2(i)
  
!CO formation
  temp1 = dexp(-dt*kCO*ndtot(i))
  eps = dt*kCO*ndtot(i); call Omexp(omeps,eps)
  ndCO(i)  = ndCOold + omeps*ndCpold
  ndCp(i)  = temp1*ndCpold
  ndCOold = ndCO(i); ndCpold = ndCp(i); nde(i) = ndp(i)+ndHep(i)+ndCp(i)
  
!CO dissociation by UV & electron collision & H collision
  temp1 = dexp(-dt*(kCOph+kCOde*nde(i)+kCOdH*ndHold))
  eps = dt*(kCOph+kCOde*nde(i)+kCOdH*ndHold); call Omexp(omeps,eps)
  ndC(i) = ndCold + omeps*ndCOold
  ndCO(i)= temp1*ndCOold
  ndCold = ndC(i); ndCOold = ndCO(i)
  
!C recombination & ionization by CR & UV
  temp1 = kCrec*nde(i)
  temp2 = dexp(-dt*(kCi+temp1))
  temp3 = 1.d0/(kCi+temp1)
  call Omexp(omeps,dt*(kCi+temp1))
  ndC(i) = ( temp1*ndCold +   kCi*ndCold*temp2 + temp1*ndCpold*omeps )*temp3
  ndCp(i)= (  kCi*ndCpold + temp1*ndCpold*temp2+   kCi*ndCold*omeps )*temp3
  ndCold = ndC(i); ndCpold = ndCp(i); nde(i) = ndp(i)+ndHep(i)+ndCp(i)

!H, He, C ionization by e collision
   ndH(i) = ndHold *dexp(-dt*kHie *nde(i))
  ndHe(i) = ndHeold*dexp(-dt*kHeie*nde(i))
   ndC(i) = ndCold *dexp(-dt*kCie *nde(i))
  call Omexp(omeps,dt*kHie *nde(i));  ndp(i) = omeps*ndHold + ndpold
  call Omexp(omeps,dt*kHeie*nde(i)); ndHep(i)= omeps*ndHeold+ ndHepold
  call Omexp(omeps,dt*kCie *nde(i)); ndCp(i) = omeps*ndCold + ndCpold
  ndHold = ndH(i); ndHeold  =  ndHe(i);  ndCold =  ndC(i)
  ndpold = ndp(i); ndHepold = ndHep(i); ndCpold = ndCp(i) 
  nde(i) = ndp(i)+ndHep(i)+ndCp(i)

!H, He, C ionization by H or H2 or p collision
   temp1 = 1.d0/( 1.d0+kHiH*dt*ndHold )
   ndH(i) = ndHold*temp1
   ndp(i) = ndpold+temp1*kHiH*dt*ndHold**2
   temp1 = ndH(i)+ndp(i)+ndH2old
  ndHe(i) = ndHeold*dexp(-dt*kHeiH*temp1)
   ndC(i) = ndCold *dexp(-dt*kCiH *temp1)
  call Omexp(omeps,dt*kHeiH*temp1); ndHep(i)= omeps*ndHeold+ ndHepold
  call Omexp(omeps,dt*kCie *temp1); ndCp(i) = omeps*ndCold + ndCpold
  ndHold = ndH(i); ndHeold  =  ndHe(i);  ndCold =  ndC(i)
  ndpold = ndp(i); ndHepold = ndHep(i); ndCpold = ndCp(i) 
  nde(i) = ndp(i)+ndHep(i)+ndCp(i)

!H2 dissiciation by Hep recombination
  temp1 = ndHepold-ndH2old
  if(temp1.ge.0.d0) then
    temp3 = dexp(-dt*kH2dHep*temp1)
    temp2 = 1.d0/( ndHepold-ndH2old*temp3 )
    ndHep(i) = ndHepold*temp1*temp2
     ndH2(i) =  ndH2old*temp1*temp2*temp3
    temp3 = ndHepold*( 1.d0-temp1*temp2 )
      ndp(i) =  ndpold + temp3
      ndH(i) =  ndHold + temp3
     ndHe(i) = ndHeold + temp3
  else
    temp3 = dexp(dt*kH2dHep*temp1)
    temp2 = 1.d0/( ndHepold*temp3-ndH2old )
    ndHep(i) = ndHepold*temp1*temp2*temp3
     ndH2(i) =  ndH2old*temp1*temp2
    temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
      ndp(i) =  ndpold + temp3
      ndH(i) =  ndHold + temp3
     ndHe(i) = ndHeold + temp3
  end if
   ndHold =  ndH(i);   ndpold =   ndp(i); ndH2old = ndH2(i)
  ndHeold = ndHe(i); ndHepold = ndHep(i) 

!CO dissiciation by Hep recombination
  temp1 = ndHepold-ndCOold
  if(temp1.ge.0.d0) then
    temp3 = dexp(-dt*kCOdHep*temp1)
    temp2 = 1.d0/( ndHepold-ndCOold*temp3 )
    ndHep(i) = ndHepold*temp1*temp2
     ndCO(i) =  ndCOold*temp1*temp2*temp3
    temp3 = ndHepold*( 1.d0-temp1*temp2 )
     ndCp(i) = ndCpold + temp3
     ndHe(i) = ndHeold + temp3
  else
    temp3 = dexp(dt*kCOdHep*temp1)
    temp2 = 1.d0/( ndHepold*temp3-ndCOold )
    ndHep(i) = ndHepold*temp1*temp2*temp3
     ndCO(i) =  ndCOold*temp1*temp2
    temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
     ndCp(i) = ndCpold + temp3
     ndHe(i) = ndHeold + temp3
  end if

    ndH(i) = dmax1( rmin,ndH(i)   )
    ndp(i) = dmax1( ndpmin,ndp(i)   )
   ndH2(i) = dmax1( ndH2min ,ndH2(i)  )
   ndHe(i) = dmax1( rmin*0.09d0,ndHe(i)  )
  ndHep(i) = dmax1( ndHepmin,ndHep(i) )
    ndC(i) = dmax1( ndCmin  ,ndC(i)   )
   ndCp(i) = dmax1( ndCpmin ,ndCp(i)  )
   ndCO(i) = dmax1( ndCOmin ,ndCO(i)  )

  nde(i) = ndp(i)+ndHep(i)+ndCp(i)
  ndtot(i) = ndp(i)+ndH(i)+2.d0*ndH2(i)+ndHe(i)+ndHep(i)
  U(i,1) = mH*ndp(i)+mH*ndH(i)+mH2*ndH2(i)+mHe*ndHe(i)+mHe*ndHep(i)
end do

end if

END SUBROUTINE SOURCE


SUBROUTINE Stblty(tLMT)
USE comvar

double precision  tLMT,alpha,tauC,Nn,Tn
double precision  CooL

tLMT = tfinal

do i = 1, Ncell
  Nn = ndH(i)+ndp(i)+ndH2(i)+ndHe(i)+ndHep(i)
  Tn = U(i,5)/(kb*Nn)
  gammi1 =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
  gammi1 = ( 2.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+2.d0*ndH2(i) )/gammi1
!*** Stability for Conduction ***!
  alpha = gammi1*Kcond*dsqrt(Tn)
  alpha = Nn*kb*dx(i)**2/alpha
!*** Avoid over cooling ***!
  Call Fcool(CooL,Tn,i)
  tauC = U(i,5)/gammi1/dabs(CooL)

  tLMT =  dmin1( tLMT, alpha )
  tLMT =  dmin1( tLMT, tauC  )
end do
if(tLMT.lt.0.d0) write(*,*) time,NRANK,'err at Stblty'

END SUBROUTINE Stblty


SUBROUTINE IMC( P,n,dt,i )
USE comvar
double precision P,n,T,dt
double precision Pu,Pd,Pm,fev,iud
double precision CooL,Pmold,nkbi
Pu    = 1.d10
Pd    = 1.d-10
Pm    = 1.d1
Pmold = -1.d0
nkbi  = 1.d0/(kb*n)
do kkk = 1,30
  T = Pm*nkbi
  call Fcool(CooL,T,i)
  fev = Pm + CooL*dt*gammi1 - P
  iud = dsign(1.d0,fev)
  Pd  = dmax1(-iud*Pm,Pd)
  Pu  = dmin1(2.d0*Pm/(iud+1.d0+1.d-10),Pu)
  Pm  = 1.d1**(0.5d0*(dlog10(Pu)+dlog10(Pd)))
  if(dabs(Pm-Pmold).lt.1.d-5) goto 835
  Pmold = Pm
end do
835 continue
P = Pm
END SUBROUTINE IMC


SUBROUTINE Omexp(omeps,eps)
double precision :: omeps,eps
if(eps.ge.1.d-4) omeps = 1.d0-dexp(-eps)
if(eps.lt.1.d-4) omeps = eps - 0.5d0*eps**2 + 0.16666666666666667d0*eps**3 - 4.16666666666666667d-2*eps**4
END SUBROUTINE Omexp



SUBROUTINE Fcool(CooL,T,i)
USE comvar
double precision :: CooL,T,Av1,Av2,x1,x2
double precision :: Laml,Lamc,Lamo,Lamd,LCOr,LCOH,LCOH2,pha,ncr
double precision :: Gampe,Gamcr,Gampd
double precision :: ATN1,ATN2,SHLD1,SHLD2
double precision :: tau1,tau2,ct1,ct2,ym1,ym2,fes
double precision :: tC1,tC2,fesC1,fesC2,tO1,tO2,fesO1,fesO2
double precision :: n1,n2,b21

double precision :: heat , KIcool,samm,sammm,sammmm,st,sss,check1,check2

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

double precision ::   NN,AbsL,T0
  NN = U(i,1) /1.27d0
AbsL = 3.94656d1
T0   = 1.0d3
!-------------------------( Cooling & Heating )
Laml = 1.0d7 * dexp( -1.184d5/(T*T0+1.0d3) )
Laml = dmin1(Laml,5.d3)
Lamc = 1.4d-2 * dsqrt(T*T0) * dexp( -9.2d1/(T*T0) )
CooL = ( NN**2.d0 * ( Laml + Lamc ) - NN ) * AbsL

!heat = 39.351d0
!heat = (0.01d0) * (1.2d0)
!KIcool = heat * (10.0d0**7 * dexp(-118400.0d0/(T+1000.0d0)) + 0.014d0 * dsqrt(T) * dexp(-92.0d0/T))
!sammmm = (1.0d7)*dexp(-1.184d5/(T*T0+1.0d3))

!samm= 1.4d-2* dsqrt(T*T0)*dexp(-9.2d1/(T*T0))

!st = dsqrt(T*T0)
!sss =  sammmm + 1.4d-2 * st * samm
!check2 = AbsL*(NN - NN**2.d0 * (samm+sammmm))
!check1 = (NN**2.d0 * (samm+sammmm)-NN)*AbsL
!CooL = check1
write(*,*) check1, check2
!CooL = 0.0d0

!sammm = CooL*dt*gammi1
!write(*,*) CooL

!Av1  = 1.63542d-3*Ntot(i,1); x1 = 6.1714d3*NH2(i,1)
!Av2  = 1.63542d-3*Ntot(i,2); x2 = 6.1714d3*NH2(i,2)
!ATN1 = ( dexp(-2.5d0*Av1)     + dexp(-2.5d0*Av2) ) * 0.5d0
!ATN2 = ( dexp(-3.77358d0*Av1) + dexp(-3.77358d0*Av2) ) * 0.5d0
!pha  = G0 * dsqrt(1.d3*T) * ATN1 / nde(i)
!------------------------- Lya Cooling
!Laml = ndH(i)*nde(i) * 1.44049d9*dexp( -1.184d2/T)
!------------------------- CII Cooling L
!call fesc(tCII(i,1),fesC1); call fesc(tCII(i,2),fesC2); b21 = fesC1+fesC2
!call LEVC2(T,b21,ndH(i),ndH2(i),nde(i),ndCp(i),n1,n2)
!Lamc = 6.0157d7*n2*b21
!***------------------------- OI Cooling
!tO1  = 11.6618d0*Ntot(i,1)*xo ; tO2 = 11.6618d0*Ntot(i,2)*xo
!call fesc(tO1,fesO1); call fesc(tO2,fesO2)
!Lamo = (dmax1(ndtot(i)*xo-ndCO(i),0.d0)/xo) * (ndH(i)+0.5d0*ndH2(i)) * 1.23916d1 * (T**0.4d0) * dexp( -0.228d0/T )
!Lamo = Lamo * (fesO1+fesO2)
!***------------------------- DustRec Cooling
!Lamd = nde(i)*ndtot(i)*6.06236d0 * (T**0.94d0) * ( pha**( 0.462628d0/(T**6.8d-2) ) )
!***------------------------- CO Cooling
!ncr  = 3.3d6*(T**0.75d0)/(ndH(i)+ndp(i)+ndH2(i)+ndHe(i)+ndHep(i))
!tau1 = 1.33194d1*NCO(i,1)/(T*dv); tau2 = 1.33194d1*NCO(i,2)/(T*dv)
!ct1  = tau1*dsqrt( 6.283d0*dlog(2.13d0+(tau1*0.36788d0)**2) ); ct2 = tau2*dsqrt( 6.283d0*dlog(2.13d0+(tau2*0.36788d0)**2) )
!ym1  = dlog( 1.d0+ct1/(1.d0+1.d1*ncr) ); ym2  = dlog( 1.d0+ct2/(1.d0+1.d1*ncr) )
!fes  = (2.d0+ym1+0.6d0*ym1**2)/(1.d0+ct1+ncr+1.5d0*dsqrt(ncr)) + (2.d0+ym2+0.6d0*ym2**2)/(1.d0+ct2+ncr+1.5d0*dsqrt(ncr))
!LCOr = ndCO(i) * 1.91505d10*(T**2)*fes
!LCOH = ndH(i) *ndCO(i) * 7.96086d4*dsqrt(T)*dexp(-(2.d0/T)**3.43d0)*dexp(-3.08d0/T)
!LCOH2= ndH2(i)*ndCO(i) * 3.60834d4*T*dexp(-(3.14d2/T)**0.333d0)*dexp(-3.08d0/T)
!***------------------------- Photo-electric Heating
!Gampe = ndtot(i) * 2.56526d3 * G0*ATN1 * &
!       ( 7.382d-3*(T**0.7d0)/(1.d0+2.d-4*pha) + 4.9d-2/(1.d0+4.d-3*(pha**0.73d0)) )
!------------------------- CR Heating
!Gamcr = (ndH(i)+ndHe(i)+ndH2(i)) * 1.89435d0
!***------------------------- Photo-destruction Heating
!SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
!SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))
!Gampd = ndH2(i) * 4.16362d4 * G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0

!CooL  = Laml + Lamc + Lamo + Lamd + LCOr + LCOH + LCOH2 - Gampe - Gamcr - Gampd

END SUBROUTINE Fcool


SUBROUTINE fesc(tau,fes)
double precision tau,fes
fes =   (0.5d0+dsign(0.5d0,0.1d0-1.d-16-tau))*(0.5d0-0.585d0*tau+0.4563d0*tau**2) &
      + (0.5d0+dsign(0.5d0,11.9025d0-(tau-3.55d0)**2))*(1.d0-dexp(-2.34d0*tau))/(4.68d0*tau+1.d-8) &
      + (0.5d0+dsign(0.5d0,tau-7.d0))/(4.d0*tau*dsqrt(dlog(dmax1(tau*0.56419d0,4.d0)))+1.d-8)
END SUBROUTINE fesc


SUBROUTINE LEVC2(T,b21,ndH,ndH2,nde,ndCp,n1,n2)
double precision :: T,b21,ndH,ndH2,nde,ndCp,n1,n2
double precision :: c21,c12,A21
c21 = 8.854d-8*nde/dsqrt(T) + 9.399d-10*(T**0.07d0)*(ndH+0.5d0*ndH2); c12 = c21*dexp(-0.092d0/T); A21 = 2.4d-6*b21
n2  = ndCp*c12/(c12+c21+A21); n1  = ndCp*(c21+A21)/(c12+c21+A21)
END SUBROUTINE LEVC2


SUBROUTINE SHIELD()
USE comvar
double precision, dimension(-1:nd) :: temp
double precision :: tNtot,tNH2,tNC,tNCO,ttau,dxh
double precision :: dndtot,dndH2,dndC,dndCO,dtemp,grad,dxhp,dxhm,dxi
double precision :: T,fesC1,fesC2,n1,n2,b21,dvin
  
dvin = 1.d0/dv
do i = 1, Ncell
  call fesc(tCII(i,1),fesC1); call fesc(tCII(i,2),fesC2); b21 = fesC1+fesC2
  T = U(i,5)/( kb*(ndp(i)+ndH(i)+ndH2(i)+ndHe(i)+ndHep(i)) )
  call LEVC2(T,b21,ndH(i),ndH2(i),nde(i),ndCp(i),n1,n2)
  temp(i) = 1.15563d1*(n1-n2)*dvin
end do

CALL BC(Ncell,ndtot(-1),BCx1,BCx2)
CALL BC(Ncell,ndH2(-1) ,BCx1,BCx2)
CALL BC(Ncell,ndC(-1)  ,BCx1,BCx2)
CALL BC(Ncell,ndCO(-1) ,BCx1,BCx2)
CALL BC(Ncell,temp(-1) ,BCx1,BCx2)
  
tNtot=0.d0; tNH2=0.d0; tNC=0.d0; tNCO=0.d0; ttau=0.d0
do i = 1, Ncell
  dxh  = 0.5d0*dx(i)
  dxhp = 0.5d0*dx(i+1)
  dxhm = 0.5d0*dx(i-1)
  dxi  = 1.d0/( dx(i)+dxhp+dxhm )
  grad   = ( ndtot(i+1)-ndtot(i-1) ) * dxi; dndtot = 0.5d0*dxh**2*grad
  grad   = ( ndH2(i+1) -ndH2(i-1)  ) * dxi; dndH2  = 0.5d0*dxh**2*grad
  grad   = ( ndC(i+1)  -ndC(i-1)   ) * dxi; dndC   = 0.5d0*dxh**2*grad
  grad   = ( ndCO(i+1) -ndCO(i-1)  ) * dxi; dndCO  = 0.5d0*dxh**2*grad
  grad   = ( temp(i+1) -temp(i-1)  ) * dxi; dtemp  = 0.5d0*dxh**2*grad

  tNtot = tNtot + dxh * ndtot(i) - dndtot
  tNH2  = tNH2  + dxh *  ndH2(i) - dndH2
  tNC   = tNC   + dxh *   ndC(i) - dndC
  tNCO  = tNCO  + dxh *  ndCO(i) - dndCO
  ttau  = ttau  + dxh *  temp(i) - dtemp
  Ntot(i,1) = tNtot
  NH2(i,1)  = tNH2
  NnC(i,1)  = tNC
  NCO(i,1)  = tNCO
  tCII(i,1) = ttau
  tNtot = tNtot + dxh * ndtot(i) + dndtot
  tNH2  = tNH2  + dxh *  ndH2(i) + dndH2
  tNC   = tNC   + dxh *   ndC(i) + dndC
  tNCO  = tNCO  + dxh *  ndCO(i) + dndCO
  ttau  = ttau  + dxh *  temp(i) + dtemp
end do
tNtot=0.d0; tNH2=0.d0; tNC=0.d0; tNCO=0.d0; ttau=0.d0
do i = Ncell, 1, -1
  dxh  = 0.5d0*dx(i)
  dxhp = 0.5d0*dx(i+1)
  dxhm = 0.5d0*dx(i-1)
  dxi  = 1.d0/( dx(i)+dxhp+dxhm )
  grad   = ( ndtot(i+1)-ndtot(i-1) ) * dxi; dndtot = 0.5d0*dxh**2*grad
  grad   = ( ndH2(i+1) -ndH2(i-1)  ) * dxi; dndH2  = 0.5d0*dxh**2*grad
  grad   = ( ndC(i+1)  -ndC(i-1)   ) * dxi; dndC   = 0.5d0*dxh**2*grad
  grad   = ( ndCO(i+1) -ndCO(i-1)  ) * dxi; dndCO  = 0.5d0*dxh**2*grad
  grad   = ( temp(i+1) -temp(i-1)  ) * dxi; dtemp  = 0.5d0*dxh**2*grad

  tNtot = tNtot + dxh * ndtot(i) + dndtot
  tNH2  = tNH2  + dxh *  ndH2(i) + dndH2
  tNC   = tNC   + dxh *   ndC(i) + dndC
  tNCO  = tNCO  + dxh *  ndCO(i) + dndCO
  ttau  = ttau  + dxh *  temp(i) + dtemp
  Ntot(i,2) = tNtot
  NH2(i,2)  = tNH2
  NnC(i,2)  = tNC
  NCO(i,2)  = tNCO
  tCII(i,2) = ttau
  tNtot = tNtot + dxh * ndtot(i) - dndtot
  tNH2  = tNH2  + dxh *  ndH2(i) - dndH2
  tNC   = tNC   + dxh *   ndC(i) - dndC
  tNCO  = tNCO  + dxh *  ndCO(i) - dndCO
  ttau  = ttau  + dxh *  temp(i) - dtemp
end do

END SUBROUTINE SHIELD


SUBROUTINE RATES(i,T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec, &
                 kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)
USE comvar
double precision :: T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi, &
          kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
double precision :: Av1,Av2,x1,x2,ATN2,ATN3,ATN4,SHLD1,SHLD2,SHLC1,SHLC2

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

Av1  = 1.63542d-3*Ntot(i,1); x1 = 6.1714d3*NH2(i,1)
Av2  = 1.63542d-3*Ntot(i,2); x2 = 6.1714d3*NH2(i,2)
ATN2 = ( dexp(-3.77358d0*Av1) +dexp(-3.77358d0*Av2) )*0.5d0
ATN3 = ( dexp(-2.3585d0*Av1)  +dexp(-2.3585d0*Av2) )*0.5d0
ATN4 = ( dexp(-3.2d0*Av1)     +dexp(-3.2d0*Av2) )*0.5d0

!H Ionization
zeta  = 9.4671d-4
!H Recombination
kHrec = (0.5d0+dsign(0.5d0,15.78d0-T))*0.45d0*dlog(1.578d2/T) &
       +(0.5d0-dsign(0.5d0,15.78d0-T))*0.4d0*dsqrt(1.578d2/T)
kHrec = 2.05572d1*(T**(-0.5d0))*kHrec
!He Recombination
kHerec= 6.37623d1*(T**(-0.672d0))
!H2 formation
kH2   = 3.4569d-3*dsqrt(T)/(1.d0+1.26491d0*dsqrt(T+Tgr)+2.d0*T+8.d0*T**2)
!***H2 Photo-dissociation
SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))
kH2ph = 1.04138d3 * G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0
!H2 destruction by H collision
kH2dH = 1.07294d5*dexp(-4.39d1/T)
!H2 destruction by e collision
kH2de = 1.133d6*(T**0.35d0)*dexp(-1.02d2/T)
!***CO formation
kCO   = 1.57785d-2/(1.d0+G0*ATN3/(ndH2(i)*xo))
!***CO Photo-dissociation
kCOph = 3.155d3*G0*ATN3
!***C ionization
SHLC1 = dexp(-2.6d0*Av1-3.0857d1*NnC(i,1)-4.33743413d-4*NH2(i,1))/(1.d0+4.33743413d-4*NH2(i,1))
SHLC2 = dexp(-2.6d0*Av2-3.0857d1*NnC(i,2)-4.33743413d-4*NH2(i,2))/(1.d0+4.33743413d-4*NH2(i,2))
kCi   = 3.97d0*zeta + 6.62697d3*G0 * ( SHLC1+SHLC2 )*0.5d0
!C recombination
kCrec = 2.69267d1*(T**(-0.62d0))*(1.d0+1.d-4*dsqrt(nde(i))/(T**2))
!CO destruction by e collision
kCOde = 1.4d4*(T**(-0.5d0))*dexp(-1.14d2/T)
!CO destruction by H collision
kCOdH = 1.91858d3*(T**(-0.5d0))*dexp(-7.77d1/T)

!H  ionization by e collision
kHie  = 5.76d4*dsqrt(T)*dexp(-1.58d2/T)
!He ionization by e collision
kHeie = 2.76d4*(T**0.43d0)*dexp(-2.85d2/T)
!C  ionization by e collision
kCie  = 1.74d5*(T**0.4d0)*dexp(-1.31d2/T)

!H  ionization by H collision
kHiH  = 9.79d0*dsqrt(T)*dexp(-1.58d2/T)
!He ionization by H or p or H2 collision
kHeiH = 4.57d0*(T**0.43d0)*dexp(-2.85d2/T)
!C  ionization by H or p or H2 collision
kCiH  = 3.01d1*(T**0.4d0)*dexp(-1.31d2/T)
!H2 destruction by Hep recombination
kH2dHep = 1.1676d0*dexp(-0.035d0/T)
!CO destruction by Hep recombination
kCOdHep = 5.0491d4

END SUBROUTINE RATES

