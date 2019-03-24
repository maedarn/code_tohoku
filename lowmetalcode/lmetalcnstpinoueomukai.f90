module comvar
  DOUBLE PRECISION , dimension(:), allocatable:: U
  DOUBLE PRECISION, parameter :: kb=8.63359d0, Kcond=1.6384d-2
  DOUBLE PRECISION  :: gamma,gammi1,gammi2,gammi3,gampl1,gampl2,gampl3
  DOUBLE PRECISION  :: CFL,facdep,tfinal,time,phr(-1:400)
  DOUBLE PRECISION  :: pmin,pmax,rmin,rmax
  INTEGER :: Ncellx,Ncelly,Ncellz,iwx,iwy,iwz,maxstp,nitera
  INTEGER :: ifchem,ifthrm,ifrad,ifgrv,loopbc=2
  character(36) :: dir='/Users/maeda/Desktop/cooling/mt0.dat'
  DOUBLE PRECISION dx , dy,dz
end module comvar
MODULE chmvar
DOUBLE PRECISION, parameter :: mH=1.d0, mHe=4.d0, mH2=2.d0, mC=12.d0, mCO=28.d0!, TCMB=5.d-2
DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-4, xo=3.2d-4, dv=2.d0, Tgr=5.d-3, fgr=1.d0
!DOUBLE PRECISION, parameter :: G0=1.d0, xc=0.28d-4, xo=0.64d-4, dv=2.d0, Tgr=5.d-3, fgr=0.2d0 !1/5 solar metal
!POP0
integer :: md = 1 !savecount
!REAL*8, parameter :: G0=1.d0, xc=1.4d-5, xo=3.2d-5, dv=2.d0, Tgr=5.d-3, fgr=1.d-1!, Pen=1.d5 !POP1
!REAL*8, parameter :: G0=1.d0, xc=1.4d-6, xo=3.2d-6, dv=2.d0, Tgr=5.d-3, fgr=1.d-2!, Pen=1.d5 !POP2
!REAL*8, parameter :: G0=1.d0, xc=1.4d-7, xo=3.2d-7, dv=2.d0, Tgr=5.d-3, fgr=1.d-3!, Pen=1.d5 !POP3
!REAL*8, parameter :: G0=1.d0, xc=1.4d-8, xo=3.2d-8, dv=2.d0, Tgr=5.d-3, fgr=1.d-4!, Pen=1.d5 !POP4
!REAL*8, parameter :: G0=1.d0, xc=1.4d-5, xo=3.2d-5, dv=3.d0, Tgr=5.d-3, fgr=1.d-1!, Pen=1.d5 !POP1
!REAL*8, parameter :: G0=1.d0, xc=1.4d-6, xo=3.2d-6, dv=3.d0, Tgr=5.d-3, fgr=1.d-2!, Pen=1.d5 !POP2
!REAL*8, parameter :: G0=1.d0, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3!, Pen=1.d5 !POP3
!REAL*8, parameter :: G0=1.d0, xc=1.4d-8, xo=3.2d-8, dv=3.d0, Tgr=5.d-3, fgr=1.d-4!, Pen=1.d5 !POP4
!REAL*8, parameter :: G0=1.d-4, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POPA
!REAL*8, parameter :: G0=1.d-3, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POPB
!REAL*8, parameter :: G0=1.d-2, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POPC
DOUBLE PRECISION :: rhoold

DOUBLE PRECISION :: ndp,ndH,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
DOUBLE PRECISION, dimension(:), allocatable :: Ntot,NH2,NnC,NCO,tCII
DOUBLE PRECISION  :: ndpmin,ndHmin,ndH2min,ndHemin,ndHepmin,ndCmin,ndCpmin,ndCOmin
END MODULE chmvar

program main
  use comvar
  use chmvar

  write(*,*) 'ok1'
  allocate(U(8))
  allocate(Ntot(2), NH2(2),NnC(2),NCO(2),tCII(2) )
  call INIT
  write(*,*) 'ok2'
  call evlv
  write(*,*) 'ok3'
  deallocate(U)
  deallocate(Ntot, NH2,NnC,NCO,tCII)

end program main

subroutine init
  use comvar
  use chmvar
  integer :: Np1x, Np2x, Np1y, Np2y, Np1z, Np2z, nunit, ix, jy, kz,b,c
  double precision ::  ql1x,ql2x,ql1y,ql2y,ql1z,ql2z,dinit1,dinit2,pinit1,pinit2, &
       vinitx1,vinitx2,vinity1,vinity2,vinitz1,vinitz2,           &
       binitx1,binitx2,binity1,binity2,binitz1,binitz2
  double precision :: theta,pi,amp,xpi,ypi,zpi,phase1,phase2,phase3,kx,ky,kzz,kw
  double precision :: Hini,pini,H2ini,Heini,Hepini,Cini,COini,Cpini,dBC,totn1
  character*3 :: NPENUM
  integer i3,i4

  CFL=0.49d0
  maxstp=100000
  tfinal=1.d2
  ifchem=2
  ifthrm=0
  ifrad=0
  ifgrv=0

  pinit1=1.d5 * 1.38d0/1.52d0 * 1.d-2 ; pinit2=pinit1
  totn1=10.d0
  Hini=0.92d0*totn1; pini=0.3d-4*Hini; H2ini=0.d0*totn1; Heini=0.8d-1*totn1; Hepini=0.3d-4*Heini
  Cini=0.d0*totn1; COini=0.d0*totn1; Cpini=xc*totn1
  !dinit=1.24215093832927
  dinit1=mH*Hini+mH*pini+mH2*H2ini+mHe*Heini+mHe*Hepini; dinit2=dinit1

  write(*,*)dinit1
  !Hini=10.d0*0.92d0;pini=Hini*3.d-4;H2ini=Hini*1.d-50;Hmini=Hini*1.d-50
  !Heini=10.d0*0.08d0;Hepini=Heini*3.d-4;Cpini=10.d0*xc;COini=Cpini*1.d-50; Cini=Cpini*1.d-50

  gamma  = ( 5.d0*(Hini+pini+Heini+Hepini)+7.d0*H2ini )/( 3.d0*(Hini+pini+Heini+Hepini)+5.d0*H2ini )
  gammi1 = gamma - 1.0d0; gammi2 = gamma - 2.0d0; gammi3 = gamma - 3.0d0
  gampl1 = gamma + 1.0d0; gampl2 = gamma + 2.0d0; gampl3 = gamma + 3.0d0
  pi     = 3.14159265358979323846d0

  pmin = 1.829797d0 * 8.6336d0   !p/kb=1.d3
  pmax = 1.d10  !604.5288d0 !p/kb =7.d4
  rmin = 0.1949628d0
  rmax = 1.d10  !4168.669d0*1.27d0

  ndHmin  = rmin*0.91d0; ndpmin  = 1.d-20; ndH2min = 1.d-20; ndHemin = rmin*0.09d0
  ndHepmin= 1.d-20; ndCpmin = 1.d-20; ndCmin = 1.d-20; ndCOmin = 1.d-20

  dx = 1.d-1
  dy=dx
  dz=dx
  U(:)=0.d0
  U(1) = dinit1
  U(5) = pinit1
  ndH   = Hini
  ndp   = pini
  ndH2  = H2ini
  ndHe  = Heini
  ndHep = Hepini
  ndC   = Cini
  ndCO  = COini
  ndCp  = Cpini
  nde  = ndp+ndHep+ndCp
  ndtot = ndH+ndp+2.d0*ndH2+ndHe+ndHep
  Ntot(1)=0.d0; NH2(1)=0.d0; NnC(1)=0.d0; tCII(1)=0.d0
  Ntot(2)=0.d0; NH2(2)=0.d0; NnC(2)=0.d0; tCII(2)=0.d0

  time=0.d0
END subroutine init

subroutine evlv
  use comvar
  use chmvar
  integer i,j,k,svc
  double precision :: dt=0.d0,tLMT

  svc=0
  open(128,file=dir)
  do i =1,maxstp
  !do i =1,2
     svc=mod(i,md)
     call Couran(tLMT)
     call Stblty(tLMT)
     dt = tLMT
   !  write(*,*)U(1)
     call source(0.5d0*dt)
   !  write(*,*)U(1)
     call IandO(dt,svc)
   !  write(*,*)U(1)
     call source(0.5d0*dt)
   !  write(*,*)U(1)
     time = time + dt
     if(U(1)>2000d0*1.27)then
        write(*,*)'OKfor2000'
        goto 1902
     end if
  end do
1902 continue
  close(128)
end subroutine evlv

SUBROUTINE Couran(tCFL)
USE comvar
USE chmvar

double precision :: tCFL,c2

tCFL = tfinal
gamma =   3.d0*(ndH+ndp+ndHe+ndHep)+5.d0*ndH2
gamma = ( 5.d0*(ndH+ndp+ndHe+ndHep)+7.d0*ndH2)/gamma
c2 = ( gamma * U(5)) / U(1)
tCFL = dmin1(tCFL,dx/(dsqrt(c2) + dabs(U(2))) )
if(tCFL.lt.0.d0) write(5,*) 'err at Couran'

END SUBROUTINE Couran


SUBROUTINE SOURCE(dt)
USE comvar
USE chmvar

double precision  dt
DOUBLE PRECISION :: ndpold,ndHold,ndH2old,ndHeold,ndHepold,ndCold,ndCpold,ndCOold,T
DOUBLE PRECISION :: zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,&
     kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
DOUBLE PRECISION :: temp1,temp2,temp3,omeps,eps
DOUBLE PRECISION :: Tn,Pn,Qx,Qy,Qz
double precision  :: mmean,rtTx,rtTy,rtTz,tcd,CooL

nde = ndp+ndHep+ndCp
ndtot = ndp+ndH+2.d0*ndH2+ndHe+ndHep



!***** 2nd order explicit scheme ( with cooling function ) *****
if(ifthrm.eq.2) then

   Tn = U(5)/( kb*(ndp+ndH+ndH2+ndHe+ndHep) )
   Pn = U(5)


   Call Fcool( CooL,Tn)
   gammi1 =   3.d0*(ndH+ndp+ndHe+ndHep)+5.d0*ndH2
   gammi1 = ( 2.d0*(ndH+ndp+ndHe+ndHep)+2.d0*ndH2)/gammi1
   if( dt .le. 0.2d0*Pn/(gammi1*dabs(CooL)) ) then
      U(5) = U(5) - gammi1*CooL*dt*0.5d0 !explicit
   else
      Call IMC( U(5),ndH+ndp+ndHe+ndHep+ndH2,dt*0.5d0 ) !implicit
   end if

   Tn = U(5)/( kb*(ndp+ndH+ndH2+ndHe+ndHep) )

   !----- Cooling ---------------------------------------------------------
   Call Fcool( CooL,Tn)
   gammi1 =   3.d0*(ndH+ndp+ndHe+ndHep)+5.d0*ndH2
   gammi1 = ( 2.d0*(ndH+ndp+ndHe+ndHep)+2.d0*ndH2 )/gammi1
   if( dt .le. 0.2d0*Pn/(gammi1*dabs(CooL)) ) then
      U(5) = Pn - gammi1*CooL*dt !explicit
   else
      Call IMC( Pn,ndH+ndp+ndHe+ndHep+ndH2,dt ) !implicit
      U(5) = Pn
   end if

end if

if(ifchem.eq.2) then

  ndpold=ndp; ndHold=ndH; ndH2old=ndH2; ndHeold=ndHe
  ndHepold=ndHep; ndCold=ndC; ndCpold=ndCp; ndCOold=ndCO
  T = U(5)/kb/( ndpold+ndHold+ndH2old+ndHeold+ndHepold )
  call RATES(T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,&
       kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)

  ! H recombination & ionization by CR
  temp1 = kHrec*nde
  temp2 = dexp(-dt*(zeta+temp1))
  temp3 = 1.d0/(zeta+temp1)
  call Omexp(omeps,dt*(zeta+temp1))
  ndH = ( temp1*ndHold + zeta*ndHold*temp2 + temp1*ndpold*omeps )*temp3
  ndp = (  zeta*ndpold +temp1*ndpold*temp2 +  zeta*ndHold*omeps )*temp3
  ndHold = ndH; ndpold = ndp; nde = ndp+ndHep+ndCp


  !He recombination & ionization by CR
  temp1 = kHerec*nde
  temp2 = dexp(-dt*(zeta+temp1))
  temp3 = 1.d0/(zeta+temp1)
  call Omexp(omeps,dt*(zeta+temp1))
  ndHe  = ( temp1*ndHeold + zeta*ndHeold*temp2 + temp1*ndHepold*omeps )*temp3
  ndHep = ( zeta*ndHepold +temp1*ndHepold*temp2+  zeta*ndHeold*omeps )*temp3
  ndHeold = ndHe; ndHepold = ndHep; nde = ndp+ndHep+ndCp

  !H2 formation & dissociation by UV & electron collision
  temp1 = kH2ph + kH2de*nde
  temp2 = temp1 + 2.d0*kH2*ndtot
  temp3 = dexp(-dt*temp2)
  call Omexp(omeps,dt*temp2)
  ndH = ( ndHold+2.d0*ndH2old*omeps )*temp1 + 2.d0*kH2*ndHold*ndtot*temp3
  ndH = ndH/temp2
  ndH2= ( 2.d0*ndH2old+ndHold*omeps )*kH2*ndtot + temp1*ndH2old*temp3
  ndH2= ndH2/temp2
  ndHold = ndH; ndH2old = ndH2


  !H2 dissociation by H collision
  temp1 = ndHold + 2.d0*ndH2old
  temp2 = dexp(-dt*temp1*kH2dH)
  temp3 = 1.d0/(ndHold + 2.d0*ndH2old*temp2)
  ndH = ndHold*temp1*temp3
  ndH2= ndH2old*temp1*temp2*temp3
  ndHold = ndH; ndH2old = ndH2

!CO formation
  temp1 = dexp(-dt*kCO*ndtot)
  eps = dt*kCO*ndtot
  call Omexp(omeps,eps)
  ndCO  = ndCOold + omeps*ndCpold
  ndCp  = temp1*ndCpold
  ndCOold = ndCO; ndCpold = ndCp; nde = ndp+ndHep+ndCp

  !CO dissociation by UV & electron collision & H collision
  temp1 = dexp(-dt*(kCOph+kCOde*nde+kCOdH*ndHold))
  eps = dt*(kCOph+kCOde*nde+kCOdH*ndHold)
  call Omexp(omeps,eps)
  ndC = ndCold + omeps*ndCOold
  ndCO= temp1*ndCOold
  ndCold = ndC; ndCOold = ndCO

  !C recombination & ionization by CR & UV
  temp1 = kCrec*nde
  temp2 = dexp(-dt*(kCi+temp1))
  temp3 = 1.d0/(kCi+temp1)
  call Omexp(omeps,dt*(kCi+temp1))
  ndC = ( temp1*ndCold +   kCi*ndCold*temp2 + temp1*ndCpold*omeps )*temp3
  ndCp= (  kCi*ndCpold + temp1*ndCpold*temp2+   kCi*ndCold*omeps )*temp3
  ndCold = ndC; ndCpold = ndCp; nde = ndp+ndHep+ndCp


  !H, He, C ionization by e collision
  ndH = ndHold *dexp(-dt*kHie *nde)
  ndHe = ndHeold*dexp(-dt*kHeie*nde)
  ndC = ndCold *dexp(-dt*kCie *nde)
  call Omexp(omeps,dt*kHie *nde);  ndp = omeps*ndHold + ndpold
  call Omexp(omeps,dt*kHeie*nde); ndHep= omeps*ndHeold+ ndHepold
  call Omexp(omeps,dt*kCie *nde); ndCp = omeps*ndCold + ndCpold
  ndHold = ndH; ndHeold  =  ndHe;  ndCold =  ndC
  ndpold = ndp; ndHepold = ndHep; ndCpold = ndCp
  nde = ndp+ndHep+ndCp

  !H, He, C ionization by H or H2 or p collision
  temp1 = 1.d0/( 1.d0+kHiH*dt*ndHold )
  ndH = ndHold*temp1
  ndp = ndpold+temp1*kHiH*dt*ndHold**2
  temp1 = ndH+ndp+ndH2old
  ndHe = ndHeold*dexp(-dt*kHeiH*temp1)
  ndC = ndCold *dexp(-dt*kCiH *temp1)
  call Omexp(omeps,dt*kHeiH*temp1); ndHep= omeps*ndHeold+ ndHepold
  call Omexp(omeps,dt*kCiH *temp1); ndCp = omeps*ndCold + ndCpold
  ndHold = ndH; ndHeold  =  ndHe;  ndCold =  ndC
  ndpold = ndp; ndHepold = ndHep; ndCpold = ndCp
  nde = ndp+ndHep+ndCp

  !H2 dissiciation by Hep recombination
  temp1 = ndHepold-ndH2old
  if(temp1.ge.0.d0) then
     temp3 = dexp(-dt*kH2dHep*temp1)
     temp2 = 1.d0/( ndHepold-ndH2old*temp3 )
     ndHep = ndHepold*temp1*temp2
     ndH2 =  ndH2old*temp1*temp2*temp3
     temp3 = ndHepold*( 1.d0-temp1*temp2 )
     ndp =  ndpold + temp3
     ndH =  ndHold + temp3
     ndHe = ndHeold + temp3
  else
    temp3 = dexp(dt*kH2dHep*temp1)
    temp2 = 1.d0/( ndHepold*temp3-ndH2old )
    ndHep = ndHepold*temp1*temp2*temp3
     ndH2 =  ndH2old*temp1*temp2
    temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
      ndp =  ndpold + temp3
      ndH =  ndHold + temp3
     ndHe = ndHeold + temp3
  end if
   ndHold =  ndH;   ndpold =   ndp; ndH2old = ndH2
  ndHeold = ndHe; ndHepold = ndHep



!CO dissiciation by Hep recombination
  temp1 = ndHepold-ndCOold
  if(temp1.ge.0.d0) then
     temp3 = dexp(-dt*kCOdHep*temp1)
     temp2 = 1.d0/( ndHepold-ndCOold*temp3 )
     ndHep = ndHepold*temp1*temp2
     ndCO =  ndCOold*temp1*temp2*temp3
     temp3 = ndHepold*( 1.d0-temp1*temp2 )
     ndCp = ndCpold + temp3
     ndHe = ndHeold + temp3
  else
     temp3 = dexp(dt*kCOdHep*temp1)
     temp2 = 1.d0/( ndHepold*temp3-ndCOold )
     ndHep = ndHepold*temp1*temp2*temp3
     ndCO =  ndCOold*temp1*temp2
     temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
     ndCp = ndCpold + temp3
     ndHe = ndHeold + temp3
  end if

    ndH = dmax1( ndHmin  ,ndH   )
    ndp = dmax1( ndpmin  ,ndp   )
   ndH2 = dmax1( ndH2min ,ndH2  )
   ndHe = dmax1( ndHemin ,ndHe  )
  ndHep = dmax1( ndHepmin,ndHep )
    ndC = dmax1( ndCmin  ,ndC   )
   ndCp = dmax1( ndCpmin ,ndCp  )
   ndCO = dmax1( ndCOmin ,ndCO  )

  nde = ndp+ndHep+ndCp
  ndtot = ndp+ndH+2.d0*ndH2+ndHe+ndHep
  U(1) = mH*ndp+mH*ndH+mH2*ndH2+mHe*ndHe+mHe*ndHep

end if


END SUBROUTINE SOURCE



SUBROUTINE Stblty(tLMT)
  USE comvar
  USE chmvar

  double precision  tLMT,alpha,tauC,Nn,Tn,dl
  double precision  CooL

  !tLMT = tfinal

  Nn = ndp+ndH+ndH2+ndHe+ndHep
  Tn = U(5)/(kb*Nn)
  !-----------------------------( Courant condition )
  dl    = dmin1(dx,dy,dz)
  gammi1 =   3.d0*(ndH+ndp+ndHe+ndHep)+5.d0*ndH2
  gammi1 = ( 2.d0*(ndH+ndp+ndHe+ndHep)+2.d0*ndH2 )/gammi1
  !-----------------------------( Avoid over cooling )
  Call Fcool(CooL,Tn)
  tauC  =  U(5)/gammi1/dabs(CooL)
  tauC  =  dmax1( 0.2d0*tauC , 2.5d-4 )
  if(tauc<tLMT) then
     tLMT=tauC
  end if
  if(tLMT.lt.0.d0) write(5,*) 'err at Stblty'

END SUBROUTINE Stblty


SUBROUTINE Fcool(CooL,T)
USE comvar
USE chmvar
double precision :: CooL,T,Av1,Av2,x1,x2
double precision :: Laml,Lamc,Lamo,Lamd,LCOr,LCOH,LCOH2,pha,ncr
double precision :: Gampe,Gamcr,Gampd
double precision :: ATN1,ATN2,SHLD1,SHLD2
double precision :: tau1,tau2,ct1,ct2,ym1,ym2,fes
double precision :: tC1,tC2,fesC1,fesC2,tO1,tO2,fesO1,fesO2
double precision :: n1,n2,b21

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

Av1  = fgr*1.63542d-3*Ntot(1); x1 = 6.1714d3*NH2(1)
Av2  = fgr*1.63542d-3*Ntot(2); x2 = 6.1714d3*NH2(2)
ATN1 = ( dexp(-2.5d0*Av1)     + dexp(-2.5d0*Av2) ) * 0.5d0
ATN2 = ( dexp(-3.77358d0*Av1) + dexp(-3.77358d0*Av2) ) * 0.5d0
pha  = G0 * dsqrt(1.d3*T) * ATN1 / nde
!------------------------- Lya Cooling
Laml = ndH*nde * 1.44049d9*dexp( -1.184d2/T)
Laml = Laml + ndp*nde * 1.48d6*dexp( -8.d1/T) !ion
!_________________________ CII Cooling L
call fesc(tCII(1),fesC1); call fesc(tCII(2),fesC2); b21 = fesC1+fesC2
call LEVC2(T,b21,ndH,ndH2,nde,ndCp,n1,n2)
Lamc = 6.0157d7*n2*b21
!***------------------------- OI Cooling
tO1  = 3.40057d-6*Ntot(1)/xo ; tO2 = 3.40057d-6*Ntot(2)/xo
call fesc(tO1,fesO1); call fesc(tO2,fesO2)
Lamo = fgr*(dmax1(ndtot*xo-ndCO,0.d0)/xo) * (ndH+0.5d0*ndH2) * 1.23916d1 * (T**0.4d0) * dexp( -0.228d0/T )
Lamo = Lamo * (fesO1+fesO2)
!***------------------------- DustRec Cooling
Lamd = fgr*nde*ndtot*6.06236d0 * (T**0.94d0) * ( pha**( 0.462628d0/(T**6.8d-2) ) )
!***------------------------- CO Cooling
ncr  = 3.3d6*(T**0.75d0)/(ndH+ndp+ndH2+ndHe+ndHep)
tau1 = 1.33194d1*NCO(1)/(T*dv); tau2 = 1.33194d1*NCO(2)/(T*dv)
ct1  = tau1*dsqrt( 6.283d0*dlog(2.13d0+(tau1*0.36788d0)**2) ); ct2 = tau2*dsqrt( 6.283d0*dlog(2.13d0+(tau2*0.36788d0)**2) )
ym1  = dlog( 1.d0+ct1/(1.d0+1.d1*ncr) ); ym2  = dlog( 1.d0+ct2/(1.d0+1.d1*ncr) )
fes  = (2.d0+ym1+0.6d0*ym1**2)/(1.d0+ct1+ncr+1.5d0*dsqrt(ncr)) + (2.d0+ym2+0.6d0*ym2**2)/(1.d0+ct2+ncr+1.5d0*dsqrt(ncr))
LCOr = ndCO * 1.91505d10*(T**2)*fes
LCOH = ndH *ndCO * 7.96086d4*dsqrt(T)*dexp(-(2.d0/T)**3.43d0)*dexp(-3.08d0/T)
LCOH2= ndH2*ndCO * 3.60834d4*T*dexp(-(3.14d2/T)**0.333d0)*dexp(-3.08d0/T)
!***------------------------- Photo-electric Heating
Gampe = fgr*ndtot * 2.56526d3 * G0*ATN1 * &
       ( 7.382d-3*(T**0.7d0)/(1.d0+2.d-4*pha) + 4.9d-2/(1.d0+4.d-3*(pha**0.73d0)) )
!------------------------- CR Heating
Gamcr = (ndH+ndHe+ndH2) * 1.89435d0
!***------------------------- Photo-destruction Heating
SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))
Gampd = ndH2 * 4.16362d4 * G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0

CooL  = Laml + Lamc + Lamo + Lamd + LCOr + LCOH + LCOH2 - Gampe - Gamcr - Gampd

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
c21 = 8.854d-8*nde/dsqrt(T) + 9.399d-10*(T**0.07d0)*(ndH+ndH2); c12 = 2.d0*c21*dexp(-0.092d0/T); A21 = 2.4d-6*b21
n2  = ndCp*c12/(c12+c21+A21); n1  = ndCp*(c21+A21)/(c12+c21+A21)
END SUBROUTINE LEVC2


SUBROUTINE IMC( P,n,dt )
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
  call Fcool(CooL,T)
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


SUBROUTINE RATES(T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,&
     kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)
USE comvar
USE chmvar
DOUBLE PRECISION :: T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,&
     kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
DOUBLE PRECISION :: Av1,Av2,x1,x2,ATN2,ATN3,ATN4,ATN5,SHLD1,SHLD2,SHLC1,SHLC2

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

Av1  = 1.63542d-3*Ntot(1)*fgr; x1 = 6.1714d3*NH2(1)
Av2  = 1.63542d-3*Ntot(2)*fgr; x2 = 6.1714d3*NH2(2)
ATN2 = ( dexp(-3.77358d0*Av1) +dexp(-3.77358d0*Av2) )*0.5d0
ATN3 = ( dexp(-2.3585d0*Av1)  +dexp(-2.3585d0*Av2) )*0.5d0

!( 1 pc * 5.e-22  = 1.54285d-3 )
!( 1 pc / 9.98337d16 = 3.09084d1 )
ATN4 = dmin1(1.d0, (3.0857d1*NH2(1))**(-3.d-2)*dexp(-1.54285d-3*NH2(1)) ) &
      +dmin1(1.d0, (3.0857d1*NH2(2))**(-3.d-2)*dexp(-1.54285d-3*NH2(2)) )
ATN4 = ATN4*0.5d0

!( 1 pc / 3.e15   = 1.02857d3 )
!( 1 pc / 8.2189020e14 = 3.75439d3 )
!( 1 pc / 4.2e16  = 7.3469d1 )
!( 4.35068d15 / 1 pc = 1.40995d-3 )
!( 8.97293d18 / 1 pc = 2.90791d0  )
!if( k==33 .and. j==33 .and. i==33) write(*,*) NCO(33,33,33,1),NCO(33,33,33,2),'11'

ATN5 = &
 (0.5d0-dsign(0.5d0,NCO(1)-1.40995d-3))*dexp(-1.d0*(1.02857d3*NCO(1))**0.6) &
+(0.5d0+dsign(0.5d0,NCO(1)-1.40995d-3))*dmin1( (3.75439d3*NCO(1)+1.d-100)**(-0.75d0),(7.3469d1*NCO(1)+1.d-100)**(-1.3d0) )
ATN5 = ATN5 + &
 (0.5d0-dsign(0.5d0,NCO(2)-1.40995d-3))*dexp(-1.d0*(1.02857d3*NCO(2))**0.6) &
+(0.5d0+dsign(0.5d0,NCO(2)-1.40995d-3))*dmin1( (3.75439d3*NCO(2)+1.d-100)**(-0.75d0),(7.3469d1*NCO(2)+1.d-100)**(-1.3d0) )
ATN5 = ATN5*0.5d0
!H Ionization
zeta  = 9.4671d-4
!H Recombination
kHrec = (0.5d0+dsign(0.5d0,15.78d0-T))*0.45d0*dlog(1.578d2/T) &
       +(0.5d0-dsign(0.5d0,15.78d0-T))*0.4d0*dsqrt(1.578d2/T)
kHrec = 2.05572d1*(T**(-0.5d0))*kHrec
!He Recombination
kHerec= 6.37623d1*(T**(-0.672d0))
!H2 formation
kH2   = fgr*3.4569d-3*dsqrt(T)/(1.d0+1.26491d0*dsqrt(T+Tgr)+2.d0*T+8.d0*T**2)
!***H2 Photo-dissociation
SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))
kH2ph = 1.04138d3 * G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0
!H2 destruction by H collision
kH2dH = 1.07294d5*dexp(-4.39d1/T)
!H2 destruction by e collision
kH2de = 1.133d6*(T**0.35d0)*dexp(-1.02d2/T)
!***CO formation
kCO   = 1.57785d-2/(1.d0+G0*ATN3/(ndH2*xo))
!***CO Photo-dissociation
kCOph = 3.155d3*G0*ATN3*ATN4*ATN5
!***C ionization
SHLC1 = dexp(-2.6d0*Av1-3.0857d1*NnC(1)-4.33743413d-4*NH2(1))/(1.d0+4.33743413d-4*NH2(1))
SHLC2 = dexp(-2.6d0*Av2-3.0857d1*NnC(2)-4.33743413d-4*NH2(2))/(1.d0+4.33743413d-4*NH2(2))
kCi   = 3.97d0*zeta + 6.62697d3*G0 * ( SHLC1+SHLC2 )*0.5d0
!C recombination
kCrec = 2.69267d1*(T**(-0.62d0))*(1.d0+1.d-4*dsqrt(nde)/(T**2))
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


SUBROUTINE Omexp(omeps,eps)
DOUBLE PRECISION :: omeps,eps

omeps = ( 0.5d0+dsign(0.5d0,eps-1.d-4) )*( 1.d0-dexp(-eps) ) + &
( 0.5d0-dsign(0.5d0,eps-1.d-4) )*( eps-0.5d0*eps**2+0.16666666666666667d0*eps**3-4.16666666666666667d-2*eps**4 )

END SUBROUTINE Omexp

subroutine IandO(dt,svc)
USE comvar
USE chmvar
integer svc
double precision  dt,T30,G30,cs30
double precision :: timeIO=0.d0
DOUBLE PRECISION :: Tn
double precision  :: CooL,rhoold1,Coolold=0.d0,l1,l2,tcoo11,tcool2,lf1,lf2


Tn = U(5)/( kb*(ndp+ndH+ndH2+ndHe+ndHep) )

gammi1 =   3.d0*(ndH+ndp+ndHe+ndHep)+5.d0*ndH2
gammi1 = ( 2.d0*(ndH+ndp+ndHe+ndHep)+2.d0*ndH2)/gammi1
G30=gammi1+1.d0
cs30=dsqrt(U(5)/U(1)*G30)
!write(*,*)'cool1',U(5),abs(CooLold),gammi1
tcool1=U(5)/dabs(CooLold)/gammi1
tcool2=U(5)/dabs(CooLold)/gammi1/U(1)
l1=cs30*tcool1
l2=cs30*tcool2
lf1=dsqrt(Kcond*Tn/dabs(CooLold))
lf2=dsqrt(Kcond*Tn/U(1)/dabs(CooLold))
if(svc==0) then
write(128,*) U(1),U(5),U(1)/1.27d0,ndCp,Tn*1.d3,G30&
     ,cs30,l1,l2,lf1,lf2,tcool1,tcool2,coolold,timeIO
end if
!tauC  =  U(i,j,k,5)/gammi1/dabs(CooL)
!----- Cooling ---------------------------------------------------------
  Call Fcool( CooL,Tn)
  gammi1 =   3.d0*(ndH+ndp+ndHe+ndHep)+5.d0*ndH2
  gammi1 = ( 2.d0*(ndH+ndp+ndHe+ndHep)+2.d0*ndH2 )/gammi1
  rhoold1=U(1)
  !U(1)=rhoold+gammi1*U(i,j,k,1)/(gammi1+1.d0)/U(i,j,k,5)*CooL*dt
  !U(5) = U(5) - gammi1*CooL*dt*0.5d0
  !U(i,j,k,1)=rhoold1 - gammi1*U(i,j,k,1)/(gammi1+1.d0)/U(i,j,k,5)*CooL*dt
  U(1)=rhoold1+gammi1/(gammi1+1.d0)/U(5)*CooL*dt
  !U(1)=rhoold1-gammi1/(gammi1+1.d0)/U(5)*CooL*dt
  U(1)=dmax1(U(1),1.d-10)
  ndH=U(1)/rhoold1*ndH
  ndp=U(1)/rhoold1*ndp
  ndH2=U(1)/rhoold1*ndH2
  ndHe=U(1)/rhoold1*ndHe
  ndHep=U(1)/rhoold1*ndHep
  ndC=U(1)/rhoold1*ndC
  ndCp=U(1)/rhoold1*ndCp
  ndCO=U(1)/rhoold1*ndCO
  Coolold=CooL
  rhoold=rhoold1
!write(*,*) 'IandO2',U(30,30,30,1),U(30,30,30,5),U(30,30,30,1)/1.27d0,Tn(30,30,30)
timeIO=timeIO+dt
end subroutine IandO
