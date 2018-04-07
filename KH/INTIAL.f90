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
!goto 2003
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
!2003 continue
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
