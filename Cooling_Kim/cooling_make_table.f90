MODULE comvar
  !INTEGER, parameter :: ndx=130, ndy=130, ndz=130, ndmax=130, Dim=3 !1024^3
INTEGER, parameter :: ndx=66, ndy=66, ndz=66, ndmax=66, Dim=3 !512^3
!INTEGER, parameter :: ndx=34, ndy=34, ndz=34, ndmax=34, Dim=3
DOUBLE PRECISION, dimension(-1:ndx) :: x,dx
DOUBLE PRECISION, dimension(-1:ndy) :: y,dy
DOUBLE PRECISION, dimension(-1:ndz) :: z,dz
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U, Bcc, Blg, Vfc, EMF,Bug
DOUBLE PRECISION, dimension(:,:,:),   allocatable :: dnc, xlag, dxlagM

DOUBLE PRECISION, parameter :: kb=8.63359d0, Kcond=1.6384d-2
DOUBLE PRECISION  :: gamma,gammi1,gammi2,gammi3,gampl1,gampl2,gampl3
DOUBLE PRECISION  :: CFL,facdep,tfinal,time,phr(-1:400)
DOUBLE PRECISION  :: pmin,pmax,rmin,rmax!,Tmin
INTEGER :: Ncellx,Ncelly,Ncellz,iwx,iwy,iwz,maxstp,nitera
INTEGER :: ifchem,ifthrm,ifrad,ifgrv,iffed,loopbc=2

DOUBLE PRECISION  :: dx1,dy1,dz1,prss1,Rst1
INTEGER :: idum1,idum2
DOUBLE PRECISION  :: nad
character(45) :: dir='/Users/maedarn/Dropbox/code/code/Cooling_Kim/' !samplecnv2

integer :: time_begin_c,time_end_c1,time_end_c2,time_end_c3,time_end_c4,time_end_c5,CountPerSec, CountMax
integer :: time_end_c6,time_end_c7,time_end_c8
!double precision :: Sigmo,Tth1=2.d4,Tth2=3.5d4,asigmo=1.d1

DOUBLE PRECISION  :: dt1=1.d-4
END MODULE comvar

MODULE chmvar
double precision :: Zsolar=1.d-2,xHe=0.1d0
DOUBLE PRECISION, parameter :: mH=1.d0, mHe=4.d0, mH2=2.d0, mC=12.d0, mCO=28.d0, TCMB=5.d-2
!DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-4, xo=3.2d-4, dv=2.d0, Tgr=5.d-3, fgr=1.d0
DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-6, xo=3.2d-6, dv=2.d0, Tgr=5.d-3, fgr=1.d-2
!DOUBLE PRECISION, parameter :: G0=1.d0, xc=0.28d-4, xo=0.64d-4, dv=2.d0, Tgr=5.d-3, fgr=0.2d0 !1/5 solar metal
!DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-7, xo=3.2d-7, dv=2.d0, Tgr=5.d-3, fgr=1.d-3 !1/1000 solar metal
!POP0
!REAL*8, parameter :: G0=1.d0, xc=1.4d-5, xo=3.2d-5, dv=3.d0, Tgr=5.d-3, fgr=1.d-1, Pen=1.d5 !POP1
!REAL*8, parameter :: G0=1.d0, xc=1.4d-6, xo=3.2d-6, dv=3.d0, Tgr=5.d-3, fgr=1.d-2, Pen=1.d5 !POP2
!REAL*8, parameter :: G0=1.d0, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POP3
!REAL*8, parameter :: G0=1.d0, xc=1.4d-8, xo=3.2d-8, dv=3.d0, Tgr=5.d-3, fgr=1.d-4, Pen=1.d5 !POP4
!REAL*8, parameter :: G0=1.d-4, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POPA
!REAL*8, parameter :: G0=1.d-3, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POPB
!REAL*8, parameter :: G0=1.d-2, xc=1.4d-7, xo=3.2d-7, dv=3.d0, Tgr=5.d-3, fgr=1.d-3, Pen=1.d5 !POPC

!DOUBLE PRECISION, allocatable :: ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
!DOUBLE PRECISION, dimension(2) :: Ntot,NH2,NnC,NCO,tCII
DOUBLE PRECISION  :: ndpmin,ndHmin,ndH2min,ndHmmin,ndHemin,ndHepmin,ndCmin,ndCpmin,ndCOmin

integer, parameter :: ifile = 30, imtrx=189,num_lin=30*(3+32)/2,itrchm=1000000,itrcool=1000
integer, dimension(1:ifile) :: in_mtl
double precision, dimension(1:ifile) :: Zmetals
real(4), dimension(1:ifile,1:ifile+3,1:imtrx):: Mtl
double precision :: Zi=1.d-2,alphaB,ni,ne !,kb=1.38d-16
END MODULE chmvar

program main
use comvar
use chmvar
!implicit none
integer, parameter :: ix = 200,cnt=10!,itrcool=1000
integer :: i,j,k,jm,j_lin

double precision, parameter :: Tmin=1.d4,Tmax=1.d8,Tmp1p=1.d-2,Tmp2p=1.d6
double precision :: Tmp1 = 1.0d-2,Tmp2 = 1.0d5
double precision, parameter :: Rhomin=1.d-2,Rhomax=1.d4,CooLth=1.d-5

double precision :: dlogT,dlogT1,dlogRho,Rho1,dlogRho1
double precision, dimension(:), allocatable :: Tmp,Rho
double precision, dimension(:), allocatable :: Lambda,Lambdaff,Lambdarr,Lambdameta,LambdaHe
double precision, dimension(:,:), allocatable :: Lambda_CIE
double precision :: neold,nech,neth=1.d-6,Tch,Tth=1.d-4
DOUBLE PRECISION:: ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
DOUBLE PRECISION, dimension(2) :: Ntot,NH2,NnC,NCO,tCII
DOUBLE PRECISION :: ndph,ndHh,ndH2h,ndHmh,ndHeh,ndHeph,ndCh,ndCph,ndCOh,ndeh,ndtoth
DOUBLE PRECISION, dimension(2) :: Ntoth,NH2h,NnCh,NCOh,tCIIh
DOUBLE PRECISION :: ndpmd,ndHmd,ndH2md,ndHmmd,ndHemd,ndHepmd,ndCmd,ndCpmd,ndCOmd,ndemd,ndtotmd
DOUBLE PRECISION, dimension(2) :: Ntotmd,NH2md,NnCmd,NCOmd,tCIImd
DOUBLE PRECISION :: ndpl,ndHl,ndH2l,ndHml,ndHel,ndHepl,ndCl,ndCpl,ndCOl,ndel,ndtotl
DOUBLE PRECISION, dimension(2) :: Ntotl,NH2l,NnCl,NCOl,tCIIl
DOUBLE PRECISION :: dt,Tmpmid,CooL1,CooL2,CooLmid,nechl,nechh,nechmd,neoldl,neoldh,neoldmd,dinit1,Rhopre
double precision :: Lamll,Lamcl,Lamol,Lamdl,LCOrl,LCOHl,LCOH2l,LamH2l,Lffl,Lrrl,LCIEl,LCIEHel
double precision :: Gampel,Gamcrl,Gampdl,Lnebl
double precision :: Lamlh,Lamch,Lamoh,Lamdh,LCOrh,LCOHh,LCOH2h,LamH2h,Lffh,Lrrh,LCIEh,LCIEHeh
double precision :: Gampeh,Gamcrh,Gampdh,Lnebh
double precision :: Lamlmd,Lamcmd,Lamomd,Lamdmd,LCOrmd,LCOHmd,LCOH2md,LamH2md,Lffmd,Lrrmd,LCIEmd,LCIEHemd
double precision :: Gampemd,Gamcrmd,Gampdmd,Lnebmd,LCIE
character(2) :: num

allocate(Tmp(0:ix),Rho(0:ix))
allocate(Lambda(0:ix),Lambdaff(0:ix),Lambdarr(0:ix),Lambdameta(0:ix),LambdaHe(0:ix))
allocate(Lambda_CIE(0:ix,1:num_lin))
!allocate(Mtl(1:ifile,1:ifile+3,1:imtrx))

do i=1,ifile
in_mtl(i)=i+3
enddo


dlogT = (-dlog10(Tmin) + dlog10(Tmax))/dble(ix)
dlogRho = (-dlog10(Rhomin) + dlog10(Rhomax))/dble(ix)
Tmp(0) = Tmin
Rho(0) = Rhomin
dlogT1=dlog10(Tmin)
dlogRho1=dlog10(Rhomin)

do i=1,ix
!dlogT1=dlog10(Tmin)
dlogT1=dlogT1+dlogT
!dlogRho1=dlog10(Rhomin)
dlogRho1=dlogRho1+dlogRho
!Tmp(i)=Tmp(i-1)+10.d0**dlogT
Tmp(i)=10.d0**dlogT1
Rho(i)=10.d0**dlogRho1
write(*,*) Rho(i),Tmp(i)
enddo

call metal_cool_corona()

do k=1,ifile
write(num,'(I2.2)') k
open(100,file=dir//'Mtln'//num//'.dat' ,access='stream',FORM='FORMATTED')!, position='append')

do i=0,ix
T=Tmp(i)
call linear(Mtl(k,1,:),Mtl(k,in_mtl(k),:),imtrx,sngl(T),LCIE)
write(100,*) T,LCIE
write(*,*) k,T,LCIE
enddo

close(100)

enddo




deallocate(Tmp,Rho)
deallocate(Lambda,Lambdaff,Lambdarr,Lambdameta,LambdaHe)
!deallocate(Mtl)
end program main

!http://wise-obs.tau.ac.il/~orlyg/ion_by_ion/

SUBROUTINE chem(T,ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII,dt)
USE comvar
USE chmvar
double precision  dt
DOUBLE PRECISION :: ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
DOUBLE PRECISION, dimension(2) :: Ntot,NH2,NnC,NCO,tCII
DOUBLE PRECISION :: ndpold,ndHold,ndH2old,ndHmold,ndHeold,ndHepold,ndCold,ndCpold,ndCOold,T
DOUBLE PRECISION :: zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,&
kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
DOUBLE PRECISION :: kHm,kH2m,kHmde
DOUBLE PRECISION :: temp1,temp2,temp3,omeps,eps
!DOUBLE PRECISION, dimension(:,:,:), allocatable :: Tn,Pn,Qx,Qy,Qz
double precision  :: mmean,rtTx,rtTy,rtTz,tcd,CooL

nde = ndp+ndHep+ndCp
ndtot = ndp+ndH+2.d0*ndH2+ndHe+ndHep

!if(ifrad.eq.2) then
!call SHIELD()
call SHIELD(T,ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)
!end if

ndpold=ndp; ndHold=ndH; ndH2old=ndH2; ndHmold=ndHm; ndHeold=ndHe
ndHepold=ndHep; ndCold=ndC; ndCpold=ndCp; ndCOold=ndCO
!T = U(i,j,k,5)/kb/( ndpold+ndHold+ndH2old+ndHeold+ndHepold )
!call RATES(T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH, &
!kCiH,kCOdHep,kH2dHep,kHm,kH2m,kHmde)
call RATES(T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie, &
kHiH,kHeiH,kCiH,kCOdHep,kH2dHep,kHm,kH2m,kHmde,ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)

!! H recombination & ionization by CR
temp1 = kHrec*nde
temp2 = dexp(-dt*(zeta+temp1))
temp3 = 1.d0/(zeta+temp1)
call Omexp(omeps,dt*(zeta+temp1))
ndH = ( temp1*ndHold + zeta*ndHold*temp2 + temp1*ndpold*omeps )*temp3
ndp = (  zeta*ndpold +temp1*ndpold*temp2 +  zeta*ndHold*omeps )*temp3
ndHold = ndH; ndpold = ndp; nde = ndp+ndHep+ndCp

!! He recombination & ionization by CR
temp1 = kHerec*nde
temp2 = dexp(-dt*(zeta+temp1))
temp3 = 1.d0/(zeta+temp1)
call Omexp(omeps,dt*(zeta+temp1))
ndHe  = ( temp1*ndHeold + zeta*ndHeold*temp2 + temp1*ndHepold*omeps )*temp3
ndHep = ( zeta*ndHepold +temp1*ndHepold*temp2+  zeta*ndHeold*omeps )*temp3
ndHeold = ndHe; ndHepold = ndHep; nde = ndp+ndHep+ndCp

!! H2 formation & dissociation by UV & electron collision
temp1 = kH2ph + kH2de*nde
temp2 = temp1 + 2.d0*kH2*ndtot
temp3 = dexp(-dt*temp2)
call Omexp(omeps,dt*temp2)
ndH = ( ndHold+2.d0*ndH2old*omeps )*temp1 + 2.d0*kH2*ndHold*ndtot*temp3
ndH = ndH/temp2
ndH2= ( 2.d0*ndH2old+ndHold*omeps )*kH2*ndtot + temp1*ndH2old*temp3
ndH2= ndH2/temp2
ndHold = ndH; ndH2old = ndH2

!! H2 dissociation by H collision
temp1 = ndHold + 2.d0*ndH2old
temp2 = dexp(-dt*temp1*kH2dH)
temp3 = 1.d0/(ndHold + 2.d0*ndH2old*temp2)
ndH = ndHold*temp1*temp3
ndH2= ndH2old*temp1*temp2*temp3
ndHold = ndH; ndH2old = ndH2

!! CO formation
temp1 = dexp(-dt*kCO*ndtot)
eps = dt*kCO*ndtot; call Omexp(omeps,eps)
ndCO  = ndCOold + omeps*ndCpold
ndCp  = temp1*ndCpold
ndCOold = ndCO; ndCpold = ndCp; nde = ndp+ndHep+ndCp


!! CO dissociation by UV & electron collision & H collision
temp1 = dexp(-dt*(kCOph+kCOde*nde+kCOdH*ndHold))
eps = dt*(kCOph+kCOde*nde+kCOdH*ndHold); call Omexp(omeps,eps)
ndC = ndCold + omeps*ndCOold
ndCO= temp1*ndCOold
ndCold = ndC; ndCOold = ndCO

!! C recombination & ionization by CR & UV
temp1 = kCrec*nde
temp2 = dexp(-dt*(kCi+temp1))
temp3 = 1.d0/(kCi+temp1)
call Omexp(omeps,dt*(kCi+temp1))
ndC = ( temp1*ndCold +   kCi*ndCold*temp2 + temp1*ndCpold*omeps )*temp3
ndCp= (  kCi*ndCpold + temp1*ndCpold*temp2+   kCi*ndCold*omeps )*temp3
ndCold = ndC; ndCpold = ndCp; nde = ndp+ndHep+ndCp

!! H, He, C ionization by e collision
 ndH = ndHold *dexp(-dt*kHie *nde)
ndHe = ndHeold*dexp(-dt*kHeie*nde)
 ndC = ndCold *dexp(-dt*kCie *nde)
call Omexp(omeps,dt*kHie *nde);  ndp = omeps*ndHold + ndpold
call Omexp(omeps,dt*kHeie*nde); ndHep= omeps*ndHeold+ ndHepold
call Omexp(omeps,dt*kCie *nde); ndCp = omeps*ndCold + ndCpold
ndHold = ndH; ndHeold  =  ndHe;  ndCold =  ndC
ndpold = ndp; ndHepold = ndHep; ndCpold = ndCp
nde = ndp+ndHep+ndCp


!! H, He, C ionization by H or H2 or p collision
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


!! H2 dissiciation by Hep recombination
temp1 = ndHepold-ndH2old
if(temp1.gt.1.d-50) then
  temp3 = dexp(-dt*kH2dHep*temp1)
  temp2 = 1.d0/( ndHepold-ndH2old*temp3 )
  ndHep = ndHepold*temp1*temp2
   ndH2 =  ndH2old*temp1*temp2*temp3
  temp3 = ndHepold*( 1.d0-temp1*temp2 )
    ndp =  ndpold + temp3
    ndH =  ndHold + temp3
   ndHe = ndHeold + temp3
endif
if(temp1.lt.-1.d-50) then
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
nde = ndp+ndHep+ndCp


!! CO dissiciation by Hep recombination
temp1 = ndHepold-ndCOold
if(temp1.gt.1.d-50) then
  temp3 = dexp(-dt*kCOdHep*temp1)
  temp2 = 1.d0/( ndHepold-ndCOold*temp3 )
  ndHep = ndHepold*temp1*temp2
   ndCO =  ndCOold*temp1*temp2*temp3
  temp3 = ndHepold*( 1.d0-temp1*temp2 )
   ndCp = ndCpold + temp3
   ndHe = ndHeold + temp3
endif
if(temp1.lt.-1.d-50) then
  temp3 = dexp(dt*kCOdHep*temp1)
  temp2 = 1.d0/( ndHepold*temp3-ndCOold )
  ndHep = ndHepold*temp1*temp2*temp3
   ndCO =  ndCOold*temp1*temp2
  temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
   ndCp = ndCpold + temp3
   ndHe = ndHeold + temp3
end if
ndHepold = ndHep; ndHeold = ndHe
 ndCOold =  ndCO; ndCpold = ndCp
nde = ndp+ndHep+ndCp

!equilibrium H- method
temp1 = kHmde+(kH2m+kHm)*ndHold
temp2 = 1.d0/(temp1+kH2m*kHm*nde*dt*ndHold)
ndH  = ndHold*temp1*temp2
ndH2 = ndH2old - 0.5d0*(ndH-ndHold)
ndHm = kHm*ndHold*nde/temp1
ndHold = ndH; ndHmold = ndHm; ndH2old = ndH2
nde = ndp+ndHep+ndCp

nde = ndp+ndHep+ndCp
ndtot = ndp+ndH+2.d0*ndH2+ndHe+ndHep
!U(i,j,k,1) = mH*ndp+mH*ndH+mH2*ndH2+mHe*ndHe+mHe*ndHep

!write(*,*) 'n',nde,ndtot,T
end SUBROUTINE chem




SUBROUTINE Fcool(CooL,ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII,T,&
Laml,Lamc,Lamo,Lamd,LCOr,LCOH,LCOH2,LamH2,Lff,Lrr,LCIE,LCIEHe,Gampe,Gamcr,Gampd,Lneb)
USE comvar
USE chmvar
double precision :: CooL,T,Av1,Av2,x1,x2,pha,ncr
double precision :: Laml,Lamc,Lamo,Lamd,LCOr,LCOH,LCOH2,LamH2,Lff,Lrr,LCIE,LCIEHe,Lneb
double precision :: LamH2H,LamH2H2,LamH2He,LamH2p,LamH2e
double precision :: alpha_B
double precision :: Sigmo,Tth1=2.d1,Tth2=3.5d1,asigmo=1.d1 !Tth1=2.d4,Tth2=3.5d4 K
double precision :: Gampe,Gamcr,Gampd
double precision :: ATN1,ATN2,SHLD1,SHLD2
double precision :: tau1,tau2,ct1,ct2,ym1,ym2,fes
double precision :: tC1,tC2,fesC1,fesC2,tO1,tO2,fesO1,fesO2
double precision :: n1,n2,b21,fneb
double precision, dimension(1:ifile) :: LCIE_each
DOUBLE PRECISION :: ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,ndHtot
DOUBLE PRECISION, dimension(2) :: Ntot,NH2,NnC,NCO,tCII

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

ndHtot = ndH+ndp+2.d0*ndH2+ndHm

Av1  = fgr*1.63542d-3*Ntot(1); x1 = 6.1714d3*NH2(1)
Av2  = fgr*1.63542d-3*Ntot(2); x2 = 6.1714d3*NH2(2)
ATN1 = ( dexp(-2.5d0*Av1)     + dexp(-2.5d0*Av2) ) * 0.5d0
ATN2 = ( dexp(-3.77358d0*Av1) + dexp(-3.77358d0*Av2) ) * 0.5d0
pha  = G0 * dsqrt(1.d3*T) * ATN1 / nde
!------------------------- Lya Cooling
!-spitzer
Laml = ndH*nde * 1.44049d9*dexp( -1.184d2/T)
Laml = Laml + ndH**2 * (1.54d8/dsqrt(T)+1.305d6*dsqrt(T)) * dexp( -1.184d2/T )
!Laml = Laml + ndp*nde * 1.48d6*dexp( -8.d1/T) !ion (imitating High-Temp region)
!-Kim
!Laml = ndH*nde*1.184d5*1.38d-16*5.31d-8*((T/1.d1)**0.15 / (1.d0+(T/5.d1)**0.65))*dexp(-11.84d0/(T/1.d1))*1.9733d27 !/ndtot
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
!------------------------- H2 cooling
LamH2H = 0.d0
if((T.gt.1.d-2).and.(T.le.1.d-1)) then
  LamH2H = ( -16.818d0+37.384d0*dlog10(T)+58.145d0*(dlog10(T))**2+48.656d0*(dlog10(T))**3 &
                     +20.1598d0*(dlog10(T))**4+3.848d0*(dlog10(T))**5 )
  LamH2H = ndH*ndH2*10.d0**(LamH2H)
endif
if((T.gt.1.d-1).and.(T.le.1.d0)) then
  LamH2H = ( -24.311d0+3.569d0*dlog10(T)-11.333d0*(dlog10(T))**2-27.85d0*(dlog10(T))**3 &
                    -21.328d0*(dlog10(T))**4-4.252d0*(dlog10(T))**5 )
  LamH2H = ndH*ndH2*10.d0**(LamH2H)
endif
if((T.gt.1.d0).and.(T.le.6.d0)) then
  LamH2H = ( -24.311d0+4.645d0*dlog10(T)-3.721d0*(dlog10(T))**2+5.937d0*(dlog10(T))**3 &
                     -5.5108d0*(dlog10(T))**4+1.5538d0*(dlog10(T))**5 )
  LamH2H = ndH*ndH2*10.d0**(LamH2H)
endif
LamH2H2 = 0.d0
if((T.gt.1.d-1).and.(T.le.6.d0)) then
  LamH2H2 = ( -23.962d0+2.094d0*dlog10(T)-0.7715d0*(dlog10(T))**2+0.4369d0*(dlog10(T))**3 &
                      -0.14913d0*(dlog10(T))**4-0.03364d0*(dlog10(T))**5 )
  LamH2H2 = ndH2*ndH2*10.d0**(LamH2H2)
endif
LamH2He = 0.d0
if((T.gt.1.d-2).and.(T.le.6.d0)) then
  LamH2He = ( -23.689d0+2.189d0*dlog10(T)-0.8152d0*(dlog10(T))**2+0.29d0*(dlog10(T))**3 &
                      -0.166d0*(dlog10(T))**4+0.1919d0*(dlog10(T))**5 )
  LamH2He = ndHe*ndH2*10.d0**(LamH2He)
endif
LamH2p = 0.d0
if((T.gt.1.d-2).and.(T.le.1.d1)) then
  LamH2p = ( -21.7167d0+1.3866d0*dlog10(T)-0.379d0*(dlog10(T))**2+0.1145d0*(dlog10(T))**3 &
                      -0.2321d0*(dlog10(T))**4+0.05854d0*(dlog10(T))**5 )
  LamH2p = ndp*ndH2*10.d0**(LamH2p)
endif
LamH2e = 0.d0
if((T.gt.1.d-2).and.(T.le.2.d-1)) then
  LamH2e = ( -34.286d0-48.537d0*dlog10(T)-77.121d0*(dlog10(T))**2-51.352d0*(dlog10(T))**3 &
                      -15.169d0*(dlog10(T))**4-0.9812d0*(dlog10(T))**5 )
  LamH2e = nde*ndH2*10.d0**(LamH2e)
endif
if((T.gt.2.d-1).and.(T.le.1.d1)) then
  LamH2e = ( -22.19d0+1.573d0*dlog10(T)-0.2134d0*(dlog10(T))**2+0.9615d0*(dlog10(T))**3 &
                      -0.910d0*(dlog10(T))**4+0.1375d0*(dlog10(T))**5 )
  LamH2e = nde*ndH2*10.d0**(LamH2e)
endif
LamH2 = (LamH2H+LamH2H2+LamH2He+LamH2p+LamH2e)*1.9733d27
!------------------------- free-free cooling
!Lff=1.422d-25*(1.d0+(0.44d0)/(1.d0+0.058d0*(dlog(T/10.d0**(5.4d0)/Zi**2))))*Zi**2.d0 * (T/1.d4)**0.5d0 * ndp * nde
!Lff=2.96d0*(1.d0+(0.44d0)/(1.d0+0.058d0*(dlog(T/10.d0**(2.4d0)/Zi**2))))*Zi**2.d0 * (T/1.d1)**0.5d0 * ndp * nde
Lff=(1.d0+(0.44d0)/(1.d0+0.058d0*(dlog(T/10.d0**(2.4d0)/Zi**2))))*Zi**2.d0 * (T/1.d1)**0.5d0 &
* ndp * nde * 1.422d-25 * 1.9733d27
!------------------------- recomb.-rad.(case B)
alpha_B=2.54d0*1.d-13*Zi*(T/1.d1/Zi/Zi)**(-0.8163d0-0.0208d0*dlog(T/1.d1/Zi/Zi))
Lrr=alpha_B * ndp * nde * (0.684d0-0.0416d0*dlog(T/1.d1/Zi/Zi))* 1.38d-16 * T * 1.d3 * 1.9733d27
!------------------------- CIE
!---He
call linear(Mtl(2,1,:),Mtl(2,in_mtl(2),:),imtrx,sngl(T*1.d3),LCIEHe)
!LCIEHe=LCIEHe*2.08d25
LCIEHe=LCIEHe*Zmetals(2)*1.9733d27 * ndHtot * ndHtot * ndHtot / Zmetals(1)
!---other metals
LCIE_each(:)=0.d0
LCIE=0.d0
do k=3,ifile
call linear(Mtl(k,1,:),Mtl(k,in_mtl(k),:),imtrx,sngl(T*1.d3),LCIE_each(k))
!LCIE=LCIE_each(k)*Zmetals(k)*fgr*2.08d25+LCIE
LCIE=LCIE_each(k)*Zmetals(k)*fgr*1.9733d27  * ndHtot * ndHtot * ndHtot / Zmetals(1) +LCIE
enddo
!------------------------- neb (Kim+32)
fneb = 6.92d-1*((dlog(T/1.d1))**0) - 5.86d-1*((dlog(T/1.d1))**1) + &
8.16d-1*((dlog(T/1.d1))**2) - 5.05d-1*(dlog((T/1.d1))**3) &
+ 1.18d-1*((dlog(T/1.d1))**4) + 7.66d-3*((dlog(T/1.d1))**5) - 5.08d-3*((dlog(T/1.d1))**6)
Lneb=Zsolar*ndtot*nde* (3.68d-23 * dexp(-3.86d1/T) * fneb) &
/ ( (T/1.d1)**(0.5d0) * (1.d0+0.12d0*(nde/100.d0)**(0.38d0-0.12d0*dlog(T/1.d1)))) * 1.9733d27

!write(*,*) fneb,Lneb,LCIE

!----Kim+23
Sigmo=1.d0/(1.d0+dexp(-asigmo*(T-0.5d0*(Tth1+Tth2))/(Tth2-Tth1)))

!CooL  = Laml + Lamc + Lamo + Lamd + LCOr + LCOH + LCOH2 - Gampe - Gamcr - Gampd + LamH2
!CooL  = (Lamc + Lamo + Lamd + LCOr + LCOH + LCOH2) * (1.d0-Sigmo) + Sigmo * (LCIEHe + LCIE) &
!- Gampe - Gamcr - Gampd + LamH2 + Laml + Lrr + Lff

Lamc = Lamc * (1.d0-Sigmo)
Lamo = Lamo * (1.d0-Sigmo)
Lamd = Lamd * (1.d0-Sigmo)
LCOr = LCOr * (1.d0-Sigmo)
LCOH = LCOH * (1.d0-Sigmo)
LCOH2 = LCOH2 * (1.d0-Sigmo)
LCIEHe = Sigmo * (LCIEHe)
LCIE = Sigmo * (LCIE)
Lneb = Lneb * (1.d0-Sigmo)

CooL  = (Lamc + Lamo + Lamd + LCOr + LCOH + LCOH2) + (LCIEHe + LCIE) &
- Gampe - Gamcr - Gampd + LamH2 + Laml + Lrr + Lff + Lneb

!write(*,*) 'FCooL', Lamc, Lamo, Lamd, LCOr, LCOH, LCOH2, LCIEHe,LCIE, Gampe, Gamcr, Gampd, LamH2, Laml, Lrr, Lff
END SUBROUTINE Fcool

SUBROUTINE linear(xa,ya,m,x,y)
integer m,i,ms
real(4) :: xa(1:m),ya(1:m)
real(4) :: x
double precision :: y1,y2,t,y

!ms=m
do i=1,m
   if(x-xa(i).le.0.d0) then
      ms=i
      !write(*,*) 'lin'
      go to 12
   endif
enddo
ms=m
12   continue

if(ms.eq.1) ms=2
y1=dble(ya(ms-1))
y2=dble(ya(ms))
!x=xa(ms)
t=(dble(x)-dble(xa(ms-1)))/(dble(xa(ms))-dble(xa(ms-1)))
y=(1.d0-t)*y1+t*y2
END SUBROUTINE linear


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


SUBROUTINE SHIELD(T,ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)
USE comvar
USE chmvar
double precision :: Nttot,NtH2,NtC,NtCO,tau,temp
double precision :: tNtot,tNH2,tNC,tNCO,ttau,dxh
double precision :: T,fesC1,fesC2,n1,n2,b21,dvin
DOUBLE PRECISION :: ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
DOUBLE PRECISION, dimension(2) :: Ntot,NH2,NnC,NCO,tCII

!ALLOCATE( Nttot(ndy,ndz,0:NSPLTx-1),NtH2(ndy,ndz,0:NSPLTx-1),NtC(ndy,ndz,0:NSPLTx-1),NtCO(ndy,ndz,0:NSPLTx-1) )
!ALLOCATE(   tau(ndy,ndz,0:NSPLTx-1),temp(ndx,ndy,ndz) )

dxh = 1.d0 ! 1 pc

dvin = 1.d0/dv
!do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
call fesc(tCII(1),fesC1); call fesc(tCII(2),fesC2); b21 = fesC1+fesC2
  !T = U(i,j,k,5)/( kb*(ndp+ndH+ndH2+ndHe+ndHep) )
call LEVC2(T,b21,ndH,ndH2,nde,ndCp,n1,n2)
  temp = 1.15563d1*(n1-n2)*dvin
!end do; end do;end do

!    tNtot =  dxh * ndtot
!    tNH2  =  dxh *  ndH2
!    tNC   =  dxh *   ndC
!    tNCO  =  dxh *  ndCO
!    ttau  =  dxh *  temp
!    Ntot(1) = tNtot
!    NH2(1)  = tNH2
!    NnC(1)  = tNC
!    NCO(1)  = tNCO
!    tCII(1) = ttau
!    tNtot = dxh * ndtot
!    tNH2  = dxh *  ndH2
!    tNC   = dxh *   ndC
!    tNCO  = dxh *  ndCO
!    ttau  = dxh *  temp

tNtot =  0.d0
tNH2  =  0.d0
tNC   =  0.d0
tNCO  =  0.d0
ttau  =  0.d0
Ntot(1) = 0.d0
NH2(1)  = 0.d0
NnC(1)  = 0.d0
NCO(1)  = 0.d0
tCII(1) = 0.d0
tNtot = 0.d0
tNH2  = 0.d0
tNC   = 0.d0
tNCO  = 0.d0
ttau  = 0.d0

    Nttot=tNtot; NtH2=tNH2; NtC=tNC; NtCO=tNCO; tau=ttau

    Ntot(2) = Ntot(1)
     NH2(2) =  NH2(1)
     NnC(2) =  NnC(1)
     NCO(2) =  NCO(1)
    tCII(2) = tCII(1)

!DEALLOCATE( Nttot,NtH2,NtC,NtCO )
!DEALLOCATE( tau,temp )
END SUBROUTINE SHIELD



!SUBROUTINE RATES(i,j,k,T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)
SUBROUTINE RATES(T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie, &
kHiH,kHeiH,kCiH,kCOdHep,kH2dHep,kHm,kH2m,kHmde,ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)
USE comvar
USE chmvar
DOUBLE PRECISION :: T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,&
kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
DOUBLE PRECISION :: kHm,kH2m,kHmde
DOUBLE PRECISION :: Av1,Av2,x1,x2,ATN2,ATN3,ATN4,ATN5,SHLD1,SHLD2,SHLC1,SHLC2
DOUBLE PRECISION :: ndp,ndH,ndH2,ndHm,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
DOUBLE PRECISION, dimension(2) :: Ntot,NH2,NnC,NCO,tCII

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

!H- formation
if(T.le.6.d0) kHm = -17.845d0+0.762d0*log10(T*1.d3)+0.1523d0*(log10(T*1.d3))**2-0.03274d0*(log10(T*1.d3))**3
if(T.gt.6.d0) kHm = -16.4199d0+0.1998d0*(log10(T*1.d3))**2-5.447d-3*(log10(T*1.d3))**4+4.0415d-5*(log10(T*1.d3))**6
kHm = 3.1557d13*(10.d0**kHm)
!H2 formation by Hm
kH2m = 3.1557d13*1.3d-9
!Hm destruction
kHmde = 3.1557d13*2.4d-7*(dexp(-0.5d0*Av1)+dexp(-0.5d0*Av2))*0.5d0*(G0*0.5848d0)
END SUBROUTINE RATES


SUBROUTINE Omexp(omeps,eps)
DOUBLE PRECISION :: omeps,eps
omeps = ( 0.5d0+dsign(0.5d0,eps-1.d-4) )*( 1.d0-dexp(-eps) ) + &
( 0.5d0-dsign(0.5d0,eps-1.d-4) )*( eps-0.5d0*eps**2+0.16666666666666667d0*eps**3-4.16666666666666667d-2*eps**4 )
END SUBROUTINE Omexp


subroutine metal_cool_corona()
!http://wise-obs.tau.ac.il/~orlyg/ion_by_ion/
use comvar
USE chmvar
integer :: i,j,k,imaxlp
!integer, dimension(1:ifile) :: in_mtl
!H,He,Li,Be,B,C,N,O,F,Ne
!Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca
!Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn
integer :: icrn1=4,icrn2=5,icrn3=6,icrn4=7,icrn5=8,icrn6=9,icrn7=10,icrn8=11,icrn9=12,icrn10=13, &
icrn11=14,icrn12=15,icrn13=16,icrn14=17,icrn15=18,icrn16=19,icrn17=20,icrn18=21,icrn19=22,icrn20=23, &
icrn21=24,icrn22=25,icrn23=26,icrn24=27,icrn25=28,icrn26=29,icrn27=30,icrn28=31,icrn29=32,icrn30=33
character(2) :: nm
double precision :: xcpre,xopre,xnHtotinitpre

Zmetals(1)=1.d0!; Zmetals(2)=8.33d-2
Zmetals(2)=1.d-1
Zmetals(3)=2.04d-9; Zmetals(4)=2.63d-11; Zmetals(5)=6.17d-10
Zmetals(6)=2.45d-4; Zmetals(7)=6.03d-5; Zmetals(8)=4.57d-4
Zmetals(9)=3.02d-8; Zmetals(10)=1.95d-4; Zmetals(11)=2.14d-6
Zmetals(12)=3.39d-5; Zmetals(13)=2.95d-6; Zmetals(14)=3.24d-5
Zmetals(15)=3.20d-7; Zmetals(16)=1.38d-5; Zmetals(17)=1.91d-7
Zmetals(18)=2.51d-6; Zmetals(19)=1.32d-7; Zmetals(20)=2.29d-6; Zmetals(21)=1.48d-9
Zmetals(22)=1.05d-7; Zmetals(23)=1.08d-8; Zmetals(24)=4.68d-7; Zmetals(25)=2.88d-7
Zmetals(26)=2.82d-5; Zmetals(27)=8.32d-8; Zmetals(28)=1.78d-6; Zmetals(29)=1.62d-8
Zmetals(30)=3.98d-8

Mtl(:,:,:)=0.e0

do k=1,ifile
write(nm,'(I2.2)') k
open(14,file=dir//nm//'.dat')
imaxlp=k+3
do j=1,imtrx
!read(14,'E8.2E2') (Mtl,i=1,imaxlp)
read(14,*) (Mtl(k,i,j),i=1,imaxlp)
enddo
close(14)
enddo

!xcpre=Zmetals(6)
!xopre=Zmetals(8)
!xnHtotinitpre=Zmetals(1)
!do k=1,ifile
!write(*,*)'Zmetal',Zmetals(k),Zmetals(k)/xnHtotinitpre*xnHtotinit
!Zmetals(k)=Zmetals(k)/Zmetals(1)*xnHtotinit
 !Zmetals(k)=Zmetals(k)/xcpre*xc
 !Zmetals(k)=Zmetals(k)/xopre*xo
!enddo
end subroutine metal_cool_corona
