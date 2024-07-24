MODULE COMVAR
integer, parameter :: N_sp=50,j_max=25
double precision :: G_0=1.d0,zeta=1.d-17,T_rad=10.d0
double precision :: A_v,xNc_H,xNc_H2,xNc_HD
double precision :: y(N_sp)
double precision ::xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
double precision :: xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,&
xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
double precision :: T_gr_K,tau_cont
double precision :: xm_p=1.67d-24,xk_B=1.38d-16,G=6.67d-8,&
pi=3.14159265358979d0,yHe=8.333d-2,yD=3.d-5,h_Pl=6.63d-27
double precision ::A0,DT0,sigma,eta_T,xmu_mol
double precision :: xLmbd_net
double precision :: Z_gas=1.0d0,Z_dust=1.0d0,Z_metal
double precision :: sigma_B=5.67d-5
END MODULE COMVAR
!new
PROGRAM ZR
!computes the evolutionary path of
!metal polluted protostellar clouds
!with radiation
!updated 2008 Jun
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer :: i_ev=1,itchem,i
!double precision :: Z_gas=1.0d0,Z_dust=1.0d0,Z_metal
double precision :: rho,y_Hp,y_H2,y_Dp,y_HD,y_H,y_D,yC,yO,y_C,y_Cp
double precision :: xmu,gamma,e,xLmbd_ch
double precision :: t,dt,t_chem,esc_cnt
double precision :: xlmbd_J,radius,xNcol,xNcolJ
double precision :: t_col,v_turb,v_D,xMJ
double precision :: c_H2

open(8,file='data/nT.dat',status='unknown')
open(11,file='data/1.dat',status='unknown')
open(12,file='data/2.dat',status='unknown')
!open(13,file='3.dat',status='unknown')
open(14,file='data/4.dat',status='unknown')
open(15,file='data/data_y.dat',status='unknown')
!open(16,file='data/ar.dat',status='unknown')

!initial condition
xnH=1.d-1
rho=(1.d0+4.d0*yHe)*xm_p*xnH
Z_metal_gas=Z_gas

do i=1,N_sp
y(i)=0.d0
enddo
!set abundance of non-zero species
y_Hp=1.d-4
y_H2=1.d-6
y_Dp=0.d0
y_HD=0.d0
y_H=1.d0-y_Hp-2.d0*y_H2-y_HD

y(1)=y_H
y(2)=y_H2
y(4)=y_Hp
y(8)=yHe

y_D=yD-y_Dp-y_HD
y(12)=y_D
y(13)=y_HD
y(14)=y_Dp

!for fiducial dust model
yC=0.927d-4*Z_metal_gas
yO=3.568d-4*Z_metal_gas

!no dust model (Anders & Grevesse 1989)
!yC=3.58d-4*Z_metal_gas
!yC=3.97d-4*Z_metal_gas
!yO=8.49d-4*Z_metal_gas

y_C=0.d0*yC
y_Cp=1.0d0*yC
y(17)=y_C
y(23)=y_Cp
y(30)=yO

y(3)=y_Hp+y_Dp+y_Cp

!dust Z
Z_metal=Z_dust

xmu=(1.d0+4.d0*yHe)/(y(1)+y(2)+y(3)+y(4)+y(8)+y(9)+y(10))

gamma=1.d0+(1.d0+4.d0*yHe)/(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))+c_H2(T_K)*y(2)))

e=xk_B*T_K/((gamma-1.d0)*xmu*xm_p)


T_gr_K=1.d0
xLmbd_ch=0.d0

t=0.d0
dt=1.d-1
itchem=0
t_chem=1.d-1
esc_cnt=1.d0

!temporal evolution
do i=1,1000000
!collapse timescale, radius, bulk velocity
!1) self-gravitating
xlmbd_J=dsqrt(pi*xk_B*T_K/(G*xmu*xm_p*rho))
radius=xlmbd_J
xNcolJ=xnH*xlmbd_J
!2) constant column density
!xNcol=1.d19
!radius=xNcol/xnH

t_col=dsqrt(3.d0*pi/(32.d0*G*rho)) !free-fall time
v_bulk=radius/3.d0/t_col !collapse velocity
v_turb=0.d0 !turb velocity
v_bulk=dsqrt(v_bulk**2+v_turb**2) !mean value

!continuum cooling
A_v=Z_metal*xnH*radius*5.3d-22 !Av by dust
call  rad_cool(T_K,T_gr_K,radius,esc_cnt,&
dt,xmu,gamma,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,i,i_ev)
tau_cont=tau_cnt

!PE heating/cooling
call  phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)

!CR heating
Gmm_CR=Gamma_CR(zeta,y_a,y_m,yHe)

v_D=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))
xNc_H2=dmin1(1.d0,v_D/v_bulk)*y_H2*xnH*radius
v_D=dsqrt(2.d0*xk_B*T_K/(1.d0*xm_p))
xNc_H=dmin1(1.d0,v_D/v_bulk)*y_H*xnH*radius
v_D=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p))
xNc_HD=dmin1(1.d0,v_D/v_bulk)*y_HD*xnH*radius

call chemcool(dt,t_chem,xmu,gamma,xLmbd_ch)

y_H=y(1)
y_H2=y(2)
y_e=y(3)
y_HD=y(13)

xLmbd_net=xLmbd_cnt+xLmbd_line+xLmbd_gr+xLmbd_ch &
-Gmm_pe-Gmm_CR

xMJ=rho*radius**3/2.d33

!update values of rho & e
call update(rho,e,p,t_col,dt,i_ev)

!data output
if((mod(i,10).eq.0).and.(i_ev==0))  then
write(*,101) t/t_col, T_K, y_H2, y_e
101 format(4e10.3)
endif
if((mod(i,10).eq.0).and.(i_ev==1))  then
write(*,102) xnH, T_K, T_gr_K, y_H2, y_e
102 format(5e10.3)
write(8,108) xnH, T_K
write(11,201) xnH,T_K,T_gr_K,p
!write(12,202) xnH,y
write(12,202) xnH,y(2),y(23),y(17),y(33)
!write(13,203) xnH,Gmm_cmp,xLmbd_cnt,xLmbd_line,xLmbd_gr,
!& xLmbd_ch,Gmm_pe,Gmm_CR
!write(*,*)'xLmbd_CI_1',xLmbd_CI
!xLmbd_CI=1.d-20!dsign(dmax1(dabs(xLmbd_CI),1.d-10),xLmbd_CI)
!write(*,*)'xLmbd_CI_2',xLmbd_CI
write(14,204) xnH, xLmbd_CII+xLmbd_CI,xLmbd_OI,&
xLmbd_CII,xLmbd_CI
!write(14,204) xnH, xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,
!& xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
!write(15,108) xMJ,T_K
!write(15,108) xNH,y(3),y(1),y(2),y(4)
!& ,y(13),y(33),y(32),y(34),y(23),y(17),y(30)
!write(16,*) xnH,T_K,T_gr_K,A_v,xNc_H2,y
endif
if((mod(i,75).eq.0).and.(i_ev==1))  then
write(15,*) xNH,y(3),y(1),y(2),y(4)&
,y(13),y(33),y(32),y(34),y(23),y(17),y(30)
endif
108     format(2e10.3)
201     format(4e10.3)
202     format(51e10.3)
203     format(8e11.3)
204     format(10e11.3)

!end of calculation
!         if(xnH.gt.1.d23) then
   if(xnH.gt.1.d14) then
      close(11)
      close(12)
      close(13)
      close(14)
      close(15)
!            close(16)
      stop
   endif

!dt for next step
   if(xLmbd_net.ne.0.d0) then
      t_cool=e/dabs(xLmbd_net)
   else
      t_cool=1.d50
   endif

   if(itchem.eq.0) then
      dt=dmax1(1.d0*t_chem,1.2d0*dt)
      if(dmin1(2.d-2*t_col,2.d-2*t_cool).le.dt) then; itchem=1; endif
   endif
   if((itchem.ne.0).and.(xnH < 1.d0)) then
      dt=dmin1(2.d-2*t_col,2.d-2*t_cool)
   endif
   if((itchem.ne.0).and.(xnH < 1.d14)) then
      dt=dmin1(2.d-2*t_col,2.d-2*t_cool)
   else
      dt=dmin1(5.d-2*t_col,5.d-2*t_cool)
   endif
!         print *, itchem, 'dt=', t_chem,
!     &        dmin1(2.d-2*t_col,2.d-2*t_cool)
!         pause
enddo
end PROGRAM ZR


SUBROUTINE update(rho,e,p,t_col,dt,i_ev)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer :: i_ev,it
double precision :: rho,e,p,rhoo,eo,po,t_col,dt
double precision :: t_ff,ft,gamma_eff
double precision :: de,drho,Gmm_cmp,dlgp,dlgrho

rhoo=rho
eo=e
po=p
t_ff=dsqrt(3.d0*pi/(32.d0*G*rho))
if(xnH < 1.d2) then
   ft=1.d0
else
   call coltime(gamma_eff,ft)
endif

do it=1,1000
   t_col=ft*t_ff
   drho=(rhoo/t_col)*dt
   Gmm_cmp=(gamma-1.d0)*eo/t_col
   de=(Gmm_cmp-xLmbd_net)*dt
   if(i_ev == 1) then
      rho=rhoo+drho
      xnH=rho/((1.d0+4.d0*yHe)*xm_p)
   endif
   e=eo+de

   T_K=e*((gamma-1.d0)*xmu*xm_p)/xk_B
   p=rho*xk_B*T_K/(xmu*xm_p)

   if(xnH < 1.d2) then
      ftn=1.d0
   else
      dlgp=(p-po)/p
      dlgrho=dt/t_col
      gamma_eff=dlgp/dlgrho
      call coltime(gamma_eff,ftn)
   endif
   err=(ftn-ft)/ftn
   if(dabs(err) < 1.d-4) exit

   if(it > 1) ftn=ft-err*(ft-fto)/(err-erro)
   erro=err
   fto=ft
   ft=ftn
enddo
t=t+dt
if((i_ev==0).and.(t/t_col > 1.d0)) then
   i_ev=1
   print *, '***********************'
!         pause
endif
!return
END SUBROUTINE update


SUBROUTINE coltime(gamma_eff,ft)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: gamma_eff,ft
double precision :: f

if(gamma_eff < 0.83d0) then
   f=0.d0
elseif(gamma_eff < 1.d0) then
   f=0.6d0+2.5d0*(gamma_eff-1.d0)-6.d0*(gamma_eff-1.d0)**2
elseif(gamma_eff < 1.33333d0) then
   f=1.d0+0.2d0*(gamma_eff-4.d0/3.d0)-2.9d0*(gamma_eff-4.d0/3.d0)**2
else
   f=1.d0
endif

if(f > 0.95d0) then
   ft=1.d0/dsqrt(0.05d0)
else
   ft=1.d0/dsqrt(1.d0-f)
endif
!return
END SUBROUTINE coltime

double precision FUNCTION Q_bg(T_nu)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision, parameter :: B_ex=5.d-4
double precision :: T_vis,delta,x,T_nu
double precision :: Q_bg_CMB,Q_bg_vis,Q_bg_gr

T_vis=6.d3
delta=4.41d3*(B_ex/T_vis**4) !stars
x=T_nu/T_rad
if(x.gt.1.d2) then
   Q_bg_CMB=0.d0
else
   Q_bg_CMB=1.d0/(dexp(x)-1.d0)
endif
Q_bg_vis=delta/(dexp(T_nu/T_vis)-1.d0)
x=T_nu/T_gr_K
if(x.gt.1.d2) then
   Q_bg_gr=0.d0
else
   Q_bg_gr=dmin1(tau_cont,1.d0)/(dexp(x)-1.d0)
endif
Q_bg=Q_bg_CMB+Q_bg_vis+Q_bg_gr
!      Q_bg=Q_bg_CMB
!      Q_bg=0.d0
return
END FUNCTION


FUNCTION func1(Tp,t,gamma,xmu,xm_p,xLmbd,dt)
IMPLICIT REAL*8(a-h,o-z)
double precision :: func1,Tp,t,gamma,xmu,xm_p,xLmbd,dt
double precision :: xk_B=1.380662d-16

func1=Tp-t-(gamma-1.d0)*xmu*xm_p*xLmbd*dt/xk_B
return
END FUNCTION func1


SUBROUTINE rad_cool(Tp,Tp_gr,radius,esc_cnt,dt,xmu,gamma,&
xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,i,i_ev)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer, parameter :: maxit=100
integer :: i,i_ev
double precision :: Tp,Tp_gr,radius,esc_cnt,dt,xmu,gamma,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr
double precision :: Tp0,Tp1,fL,f,Tpsec,TpL,swap,xLmbd_cnt1


func(t,xLmbd)=Tp-t-(gamma-1.d0)*xmu*xm_p*xLmbd*dt/xk_B !dif, function (T^n+1-T^n-T*dt/tcool)

Tp0=Tp
Tp1=Tp
if(tau_cnt < 10.d0) then
   T_K=Tp1
   call line_cool(radius,xLmbd_line,i)
else
   xLmbd_line=0.d0
endif

call cnt(Tp1,Tp_gr,radius,esc_cnt,xLmbd_cnt,xLmbd_gr)
if(dabs(Tp1/T_rad-1.d0) > 1.d0) return ! strong radiation
xLmbd=xLmbd_ch+xLmbd_line+xLmbd_cnt+xLmbd_gr
fL=func(Tp1,xLmbd)

Tp2=1.01d0*Tp1
call cnt(Tp2,Tp_gr,radius,esc_cnt,xLmbd_cnt,xLmbd_gr)
xLmbd=xLmbd_ch+xLmbd_line+xLmbd_cnt+xLmbd_gr
f=func(Tp2,xLmbd)

if(dabs(fL).lt.dabs(f)) then
   Tpsec=Tp1
   TpL=Tp2
   swap=fL
   fL=f
   f=swap
else
   TpL=Tp1
   Tpsec=Tp2
endif

!change T by cooling
do it=1,maxit
   dTp=(TpL-Tpsec)*f/(f-fL)
   TpL=Tpsec
   fL=f
   Tpsec=Tpsec+dTp
   call cnt(Tpsec,Tp_gr,radius,esc_cnt,xLmbd_cnt1,xLmbd_gr1)
   xLmbd1=xLmbd_ch+xLmbd_line+xLmbd_cnt1+xLmbd_gr1
   f=func(Tpsec,xLmbd1)
   err=dabs(f/Tpsec)
   if(err < 1.d-6) go to 12
end do
12   continue
xLmbd_cnt=xLmbd_cnt1
xLmbd_gr=xLmbd_gr1
T_K=Tp0
END SUBROUTINE rad_cool


SUBROUTINE cnt(T_K1,T_gr_K1,radius,esc_cnt,xLmbd_cnt,xLmbd_gr)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: T_K1,T_gr_K1,radius,esc_cnt,xLmbd_cnt,xLmbd_gr
!      PARAMETER (xk_vis=40.d0,B_ex=5.d-4)
!      PARAMETER (xk_vis=40.d0,B_ex=0.d0)


yHe=8.333d-2
!      xJ_ex=B_ex*dexp(-A_v/2.5d0)
xJ_ex=0.d0
call grtemp(xnH,esc_cnt,xJ_ex,T_gr_K1)
rho=(1.d0+4.d0*yHe)*xm_p*xnH
xk_gr=xkp_gr(rho,T_gr_K1)*Z_metal
xk_gas=xk_prm(T_K1,rho)
xk_cnt=xk_gas+xk_gr
tau_cnt=xk_cnt*rho*radius
if(tau_cnt.gt.1.d0) then
   esc_cnt=1.d0/(tau_cnt**2)
else
   esc_cnt=1.d0
endif
xLmbd_cnt=4.d0*sigma_B*(T_K1**4-T_rad**4)*xk_gas*esc_cnt
xLmbd_gr=4.d0*sigma_B*(T_gr_K1**4-T_rad**4)*xk_gr*esc_cnt
!     &     -xk_vis*Z_metal*xJ_ex
END SUBROUTINE cnt

SUBROUTINE line_cool(radius,xLmbd_line,i)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer, parameter :: N_p=100
double precision, SAVE :: func(N_p),esc_CI(N_p),esc_CII(N_p),esc_OI(N_p),&
esc_CImeta(N_p),esc_CIImeta(N_p),esc_OImeta(N_p)
double precision, EXTERNAL :: pop_CI,pop_OI,pop_CII,pop_CImeta,pop_OImeta,pop_CIImeta
double precision :: xLd_H2,xLd_HD,xLd
!double precision, SAVE :: esc_CI,esc_CII,esc_OI,esc_CImeta,esc_CIImeta,esc_OImeta

if(i.eq.1) then
   do l=1,N_p
      esc_CI(l)=1.d0
      esc_CII(l)=1.d0
      esc_OI(l)=1.d0
      esc_CImeta(l)=1.d0
      esc_CIImeta(l)=1.d0
      esc_OImeta(l)=1.d0
   enddo
endif

y_e=y(3)
y_a=y(1)
y_m=y(2)
y_Hp=y(4)
y_HD=y(13)
y_CO=y(33)
y_OH=y(32)
y_H2O=y(34)
y_CII=y(23)
y_CI=y(17)
y_OI=y(30)

! column density
xNc_H=xnH*radius
xNc_H2=y_m*xNc_H
xNc_CO=y_CO*xNc_H
xNc_OH=y_OH*xNc_H
xNc_H2O=y_H2O*xNc_H
xNc_CII=y_CII*xNc_H
xNc_CI=y_CI*xNc_H
xNc_OI=y_OI*xNc_H
xNc_HD=y_HD*xNc_H

v_th=dsqrt(2.d0*xk_B*T_K/xm_p)

!     H2 cooling
v_D=v_th/dsqrt(2.d0)
xNc_H2=dmin1(1.d0,v_D/v_bulk)*xNc_H2
xLmbd_H2=y_m*xLd_H2(xnH,T_K,y_a,y_m,y_e,y_Hp,xNc_H2,tau_cnt)*xnH/((1.d0+4.d0*yHe)*xm_p)

!     HD cooling
v_D=v_th/dsqrt(3.d0)
xNc_HD=dmin1(1.d0,v_D/v_bulk)*xNc_HD
xLmbd_HD=y_HD*xLd_HD(xnH,T_K,y_a,y_m,xNc_HD,tau_cnt)*xnH/((1.d0+4.d0*yHe)*xm_p)

!     CO cooling
xmu_mol=28.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CO
call COcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_CO)
call COcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_CO)
xLmbd_CO=y_CO*(xLd_CO-xLdr_CO)*xnH/((1.d0+4.d0*yHe)*xm_p)

!     OH cooling
xmu_mol=17.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OH
call OHcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_OH)
call OHcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_OH)
xLmbd_OH=y_OH*(xLd_OH-xLdr_OH)*xnH/((1.d0+4.d0*yHe)*xm_p)

!     H2O cooling
xmu_mol=18.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_H2O
call H2Ocool(xnH,T_K,y_m,xNc,tau_cnt,xLd_H2O)
call H2Ocool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_H2O)
xLmbd_H2O=y_H2O*(xLd_H2O-xLdr_H2O)*xnH/((1.d0+4.d0*yHe)*xm_p)

!     CII cooling
xmu_C=12.d0
xmu_O=16.d0
v_D=v_th/dsqrt(xmu_C)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CII
!     fine structure line
xLmbd_CII=y_CII*xLd(pop_CII,esc_CII,1)*xnH/((1.d0+4.d0*yHe)*xm_p)
!     metastable line
if(T_K > 2.d2) then
   xLmbd_CIImeta=y_CII*xLd(pop_CIImeta,esc_CIImeta,1)*xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_CII=xLmbd_CII+xLmbd_CIImeta
endif

!     CI cooling
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CI
!     fine structure line
xLmbd_CI=y_CI*xLd(pop_CI,esc_CI,3)*xnH/((1.d0+4.d0*yHe)*xm_p)

if(T_K > 3.d3) then
   xLmbd_CImeta=y_CI*xLd(pop_CImeta,esc_CImeta,3)*xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_CI=xLmbd_CI+xLmbd_CImeta
endif

!     OI cooling
v_D=v_th/dsqrt(xmu_O)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OI
!     fine structure line
if(T_K > 15.d0) then
   xLmbd_OI=y_OI*xLd(pop_OI,esc_OI,3)*xnH/((1.d0+4.d0*yHe)*xm_p)
else
   xLmbd_OI=0.d0
endif
!     metastable line
if(T_K > 3.5d3) then
   xLmbd_OImeta=y_OI*xLd(pop_OImeta,esc_OImeta,3)*xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_OI=xLmbd_OI+xLmbd_OImeta
endif

!     HI Lya cooling
if(xnH < 1.d5) then
   T_5=1.d-5*T_K
   xLmbd_Lya=y_e*y_a*7.50d-19/(1.d0+dsqrt(T_5))*dexp(-1.18348d5/T_K)*xnH/((1.d0+4.d0*yHe)*xm_p)
else
   xLmbd_Lya=0.d0
endif
!     total cooling rate by lines
xLmbd_line=xLmbd_H2+xLmbd_HD+xLmbd_CO+xLmbd_OH+xLmbd_H2O+xLmbd_CII+xLmbd_CI+xLmbd_OI+xLmbd_Lya
END SUBROUTINE line_cool


!include "subs/chemistry.f"
!include "subs/reaction.f"  ! full reactions
!include "subs/grain.f"
!include "subs/xk_prm.f"
!include "subs/H2.f"
!include "subs/HD.f"
!include "subs/lines.f"
!include "subs/CR.f"
!include "subs/math.f"
!include "subs/CO.f"
!include "subs/H2O.f"
!include "subs/OH.f"


SUBROUTINE chemcool(dt,tchem,xmu,gamma,xLmbdch)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer, parameter :: maxit=100
double precision :: ytmp(N_sp)
double precision :: Tmp_f,dTp
double precision :: tchem,tchem1,tchem2
double precision :: Tp,Tp_gr,dt,xmu,gamma,xLmbdch,Tp1,xmu1,gamma1,xLmbdch1&
,xmu2,gamma2,xLmbdch2
!     input    Tp,xnH,dt,y,xmu,gamma
!     output   y,xmu,gamma,xLmbdch

func(Tmp_f,xLmbdch)=Tp-Tmp_f-(gamma-1.d0)*xmu*xm_p*xLmbdch*dt/xk_B

Tp_gr=T_gr_K
Tp=T_K

Tp1=Tp
do j=1,N_sp
   ytmp(j)=y(j)
enddo
call chemreact(xnH,Z_metal,Tp1,Tp_gr,ytmp,dt,tchem1,xmu1,gamma1,xLmbdch1)
fL=func(Tp1,xLmbdch1)
if((xnH.le.1.d16).and.(Tp.le.1.65d3)) go to 12

Tp2=0.99d0*Tp1
do j=1,N_sp
   ytmp(j)=y(j)
enddo
call chemreact(xnH,Z_metal,Tp2,Tp_gr,ytmp,dt,tchem2,xmu2,gamma2,xLmbdch2)
f=func(Tp2,xLmbdch2)

if(dabs(fL).lt.dabs(f)) then
   Tpsec=Tp1
   TpL=Tp2
   swap=fL
   fL=f
   f=swap
else
   TpL=Tp1
   Tpsec=Tp2
endif

do i=1,maxit
   dTp=(TpL-Tpsec)*f/(f-fL)
   TpL=Tpsec
   fL=f
   Tpsec=Tpsec+dTp
   do j=1,N_sp
      ytmp(j)=y(j)
   enddo
   call chemreact(xnH,Z_metal,Tpsec,Tp_gr,ytmp,dt,tchem1,xmu1,gamma1,xLmbdch1)
   f=func(Tpsec,xLmbdch1)
   if(dabs(f/Tpsec).le.1.d-7) go to 12
enddo

12   continue
do j=1,N_sp
   y(j)=ytmp(j)
enddo

xmu=xmu1
gamma=gamma1
xLmbdch=xLmbdch1
tchem=tchem1
END SUBROUTINE chemcool


SUBROUTINE chemreact(xnH,Z_metal,T_K,T_gr_K,y   ,dt,tchem ,xmu ,gamma ,xLmbdch )
!******************************************
!*      Version 4+      (13 Sep 1998)     *
!*      no radiation                      *
!******************************************
!USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer, parameter :: N_react=675
integer, parameter :: N_sp=50
integer :: isp
double precision, parameter :: eps=1.d-4,eps_y=1.d-10,xnH_eq=3.d16
double precision :: y_init(N_sp),y_tmp(N_sp),dy(N_sp),ddy(N_sp), xk(N_react),&
r_f(N_sp),r_f_fw(N_sp),r_f_bw(N_sp), dr_fdy(N_sp,N_sp),A(N_sp,N_sp),y(N_sp)
double precision :: xnH,T_K,T_gr_K,Z_metal,dt,xmu,gamma,xLmbdch,tchem

do isp=1,N_sp
   y_init(isp)=y(isp)
   dy(isp)=0.d0
enddo

if(xnH.gt.xnH_eq) go to 50

!**************************************
!*    Non Equilibrium (Low Density)   *
!**************************************
xH=y(1)
xH2=2.d0*y(2)
xHe=y(8)

call react_coef(xnH,T_K,T_gr_K,Z_metal,xk)

do itr=1,20
   call react_rat(xk,xnH,y,r_f)
   do jsp=1,N_sp
      if(dabs(y(jsp)).le.1.d-100) then
         delta_y=eps_y
      else
         delta_y=eps*y(jsp)
      endif
      do isp=1,N_sp
         if(isp.eq.jsp) then
            y_tmp(isp)=y(isp)+5.d-1*delta_y
         else
            y_tmp(isp)=y(isp)
         endif
      enddo
      call react_rat(xk,xnH,y_tmp,r_f_fw)
      
      y_tmp(jsp)=y(jsp)-5.d-1*delta_y
      call react_rat(xk,xnH,y_tmp,r_f_bw)
      
      do isp=1,N_sp
          dr_fdy(isp,jsp)=(r_f_fw(isp)-r_f_bw(isp))/delta_y
          r_f_big=max(dabs(r_f_fw(isp)),dabs(r_f_bw(isp)))
          if(r_f_big.ne.0.d0) then
             dr_f=dabs(r_f_fw(isp)-r_f_bw(isp))
             if(dr_f/r_f_big.lt.1.d-15) then
                dr_fdy(isp,jsp)=0.d0
             endif
          endif
       enddo
   enddo

!C ------ set the matrix A
   do isp=1,N_sp
      do jsp=1,N_sp
         if(isp.ne.jsp) then
            A(isp,jsp)=-dt*dr_fdy(isp,jsp)
         else
            A(isp,jsp)=1.d0-dt*dr_fdy(isp,jsp)
         endif
      enddo
   enddo
   
!C ------- set the vector ddy
   do isp=1,N_sp
      ddy(isp)=r_f(isp)*dt-dy(isp)
   enddo

!C ------- solve linear equations
   call gaussj(A,N_sp,N_sp,ddy,1,1)
   
   do isp=1,N_sp
      dy(isp)=dy(isp)+ddy(isp)
      y(isp)=y(isp)+ddy(isp)
      if(y(isp).lt.0.d0) then; y(isp)=y_init(isp); endif
   enddo
   
   err_max=0.d0
   do isp=1,N_sp
      if(y(isp).ne.0.d0) then
         err=dabs(ddy(isp)/y(isp))
      else
         err=0.d0
      endif
      err_max=max(err,err_max)
   enddo

   if(err_max.lt.1.d-8) then
      tchem=1.d20
      do isp=1,N_sp
         if(y(isp)-y_init(isp).ne.0.d0) then
            tch=dabs((y(isp)+y_init(isp))/(2.d0*(y(isp)-y_init(isp))))*dt
            tchem=dmin1(tch,tchem)
         endif
      enddo
      go to 100
   endif
enddo

go to 100

50   continue
!***********************************************
!*        Equillibrium ( High Density )        *
!***********************************************
call equichem(xnH,T_K,y_e,y_HI,y_HII,y_H2,y_HeI,y_HeII,y_HeIII)
do i=1,N_sp
   y(i)=0.d0
enddo
y(3)=y_e
y(1)=y_HI
y(4)=y_HII
y(2)=y_H2
y(8)=y_HeI
y(9)=y_HeII
y(10)=y_HeIII

do i=1,N_sp
   dy(i)=y(i)-y_init(i)
enddo


100   continue

yHe=y(8)+y(9)+y(10)

xmu=(1.d0+4.d0*yHe)/(y(1)+y(2)+y(3)+y(4)+y(8)+y(9)+y(10))

gamma=1.d0+(1.d0+4.d0*yHe)&
/(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))+c_H2(T_K)*y(2)))

xm_p=1.67d-24
!C     H_2 formation/dissociation cooling
if(xnH.lt.1.d13) then
   xn_cr=1.d6/dsqrt(T_K)/(1.6d0*y(1)*dexp(-(4.d2/T_K)**2)&
   +1.4d0*y(2)*dexp(-1.2d4/(T_K+1.2d3)))
   crit=1.d0/(1.d0+xn_cr/xnH)
!C     2 H + grain -> H2
   rtgrain=xk(23)*y(1)
!C     H- + H -> H2 +  e
   rtHm=xk(8)*y(7)*y(1)
!C     H2+ + H  -> H2 + H+
   rtH2p=xk(10)*y(5)*y(1)
!C     3 H  -> H2 +  H
!C     2H  + H2  -> 2H2
!c     2 H + He  -> H2 + He :43
   rt3b=(xk(19)*(y(1)**3)+xk(20)*(y(1)**2)*y(2)+xk(43)*(y(1)**2)*y(8))*xnH
!C     H2 + H  -> 3 H
!C     2 H2  -> 2 H + H2
!c     H2 + He -> 2 H + He :37
   rtdis=xk(13)*y(2)*y(1)+xk(21)*(y(2)**2)+xk(37)*y(2)*y(8)
!C     H2 + ph. -> 2H
   rtphdis=xk(590)*y(2)/xnH

   xL_H2=-(rtgrain*(0.2d0+4.2d0*crit)+rtHm*3.53d0*crit &
   +rtH2p*1.83d0*crit+(rt3b*crit-rtdis)*4.48d0+rtphdis*0.4d0)&
   *xnH*1.60219d-12
   dyHp=dy(4)+(xk(2)*y(4)*y(3)*xnH-xk(544)*y(1)-xk(548)*y(2))*dt
   dyHep=dy(9)+(xk(4)*y(9)*y(3)*xnH-xk(545)*y(8))*dt
   dyHepp=dy(10)
else
   dyH2=dy(2)
   xL_H2=-7.18d-12*dyH2/dt
   dyHp=dy(4)
   dyHep=dy(9)
   dyHepp=dy(10)
endif
xLmbdch=((2.18d-11*dyHp+3.94d-11*dyHep+12.66d-11*dyHepp)/dt+xL_H2)/((1.d0+4.d0*yHe)*xm_p)

END SUBROUTINE chemreact


SUBROUTINE equichem(xnH,T_K,y_e,y_HI,y_HII,y_H2,y_HeI,y_HeII,y_HeIII)
!***********************************************
!*        Equillibrium  Chemistry    (Saha)    *
!***********************************************
!USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision, parameter :: yHe=8.333d-2,rmasse=9.109534d-28,h=6.626176d-27,rmassp=1.67d-24
double precision :: T_eV,phyc,chiH0,chiHe0,chiHep,chiH2,zH0,zHe0,zHep,zH2
double precision :: xnH,T_K,y_e,y_HI,y_HII,y_H2,y_HeI,y_HeII,y_HeIII &
,pi=3.14159265358979d0,yD=3.d-5

T_eV=8.61735d-5*T_K
phyc=2.d0*(2.d0*1.60219d-12*pi*rmasse/h**2)**1.5d0
chiH0=13.6d0
chiHe0=24.6d0
chiHep=79.0d0-chiHe0
chiH2=4.478d0
zH0=2.d0
zHe0=1.d0
zHep=2.d0
zH2=H2_pf(T_K)

fH0=phyc*(T_eV**1.5d0)*dexp(-chiH0/T_eV)/zH0/xnH
fHe0=phyc*(T_eV**1.5d0)*dexp(-chiHe0/T_eV)*(zHep/zHe0)/xnH
fHep=phyc*(T_eV**1.5d0)*dexp(-chiHep/T_eV)/zHep/xnH
fH2=(zH0**2/zH2)*(1.60219d-12*pi*rmassp/h**2)**1.5d0 &
*(T_eV**1.5d0)*dexp(-chiH2/T_eV)/xnH

if(fH2.eq.0.d0) then
   y_HI=0.d0
   y_HII=0.d0
   y_H2=5.d-1
   y_HeI=yHe
   y_HeII=0.d0
   y_HeIII=0.d0
   return
elseif(fH0.eq.0.d0) then
   y_HI=2.d0/(1.d0+dsqrt(1.d0+8.d0/fH2))
   y_HII=0.d0
   y_H2=5.d-1*(1.d0-y_HI)
   y_e=0.d0
   y_HeI=yHe
   y_HeII=0.d0
   y_HeIII=0.d0
   return
else
   y_e=1.d0+yHe
!     Newton-Raphson Loop
   do itr=1,1000
      b=5.d-1*fH2*(1.d0+fH0/y_e)
      d=2.d0*fH2/(b+dsqrt(b**2+2.d0*fH2))
      Fe=y_e**3+fHe0*(y_e**2+(fHep-yHe)*y_e-2.d0*fHep*yHe)&
      -5.d-1*fH0*(y_e+fHe0*(1.d0+fHep/y_e))*d
      if((dabs(Fe)**(1.d0/3.d0))/y_e.le.1.d-4) go to 12
      dFe=3.d0*y_e**2+fHe0*(2.d0*y_e+(fHep-yHe))-5.d-1*fH0*((1.d0-fHe0*fHep/y_e**2)*d &
      +(y_e+fHe0*(1.d0+fHep/y_e))*(-5.d-1*fH2*fH0/y_e**2)*(-1.d0+5.d-1*fH2*(1.d0+fH0/y_e) &
      /dsqrt(b**2+2.d0*fH2)))
      y_e=y_e-Fe/dFe
enddo
12      continue

   y_HI=5.d-1*d
   y_HII=fH0*y_HI/y_e
   y_H2=y_HI**2/fH2
   y_HeI=yHe/(1.d0+fHe0/y_e*(1.d0+fHep/y_e))
   y_HeII=fHe0*y_HeI/y_e
   y_HeIII=fHep*y_HeII/y_e
   return
endif
END SUBROUTINE equichem

SUBROUTINE COcool(xnH,T_K,y_H2,xNc_CO,tau_cnt,xLd_CO)
!     CO cooling function
IMPLICIT REAL*8(a-h,o-z)
double precision :: xnH,T_K,y_H2,xNc_CO,tau_cnt,xLd_CO
double precision :: xk_B=1.38d-16, xm_p=1.67d-24
double precision :: xlTa(1:13)
data xlTa/0.477d0, 0.778d0, 1.000d0, 1.301d0, 1.477d0,&
1.699d0, 1.903d0, 2.000d0, 2.477d0, 2.778d0,&
3.000d0, 3.176d0, 3.301d0/
double precision :: xlNa(1:11)
data xlNa/14.0d0, 14.5d0, 15.0d0, 15.5d0, 16.0d0,&
16.5d0, 17.0d0, 17.5d0, 18.0d0, 18.5d0, 19.0d0/
double precision :: aL0a(1:13)
data  aL0a/0.2585d2, 0.2517d2, 0.2477d2, 0.2438d2, 0.2421d2,&
0.2403d2,0.2389d2, 0.2382d2, 0.2342d2, 0.2313d2, 0.2291d2, 0.2263d2,&
0.2228d2/
double precision :: aLLTEa(1:13,1:11)
data aLLTEa/0.2251d2, 0.2165d2, 0.2108d2, 0.2035d2, 0.1994d2,&
0.1945d2,&
0.1901d2, 0.1880d2, 0.1781d2, 0.1723d2, 0.1686d2, 0.1666d2,&
0.1655d2, 0.2254d2, 0.2167d2, 0.2109d2, 0.2035d2, 0.1995d2,&
0.1945d2, 0.1901d2, 0.1880d2, 0.1781d2, 0.1723d2, 0.1686d2,&
0.1666d2, 0.1655d2, 0.2262d2, 0.2171d2, 0.2111d2, 0.2037d2,&
0.1996d2, 0.1946d2, 0.1901d2, 0.1880d2, 0.1781d2, 0.1723d2,&
0.1686d2, 0.1666d2, 0.1655d2, 0.2282d2, 0.2183d2, 0.2118d2,&
0.2040d2, 0.1998d2, 0.1947d2, 0.1902d2, 0.1881d2, 0.1782d2,&
0.1723d2, 0.1687d2, 0.1666d2, 0.1655d2, 0.2317d2, 0.2208d2,&
0.2137d2, 0.2051d2, 0.2005d2, 0.1952d2, 0.1905d2, 0.1883d2,&
0.1782d2, 0.1723d2, 0.1687d2, 0.1666d2, 0.1655d2, 0.2362d2,&
0.2244d2, 0.2167d2, 0.2073d2, 0.2023d2, 0.1964d2, 0.1913d2,&
0.1890d2, 0.1785d2, 0.1725d2, 0.1688d2, 0.1667d2, 0.1656d2,&
0.2405d2, 0.2284d2, 0.2204d2, 0.2105d2, 0.2052d2, 0.1987d2,&
0.1932d2, 0.1906d2, 0.1792d2, 0.1728d2, 0.1690d2, 0.1669d2,&
0.1658d2, 0.2450d2, 0.2326d2, 0.2244d2, 0.2142d2, 0.2086d2,&
0.2019d2, 0.1960d2, 0.1933d2, 0.1808d2, 0.1738d2, 0.1697d2,&
0.1675d2, 0.1663d2, 0.2500d2, 0.2373d2, 0.2287d2, 0.2182d2,&
0.2124d2, 0.2055d2, 0.1995d2, 0.1966d2, 0.1834d2, 0.1759d2,&
0.1715d2, 0.1691d2, 0.1678d2, 0.2549d2, 0.2417d2, 0.2330d2,&
0.2223d2, 0.2165d2, 0.2094d2, 0.2032d2, 0.2003d2, 0.1867d2,&
0.1789d2, 0.1748d2, 0.1726d2, 0.1712d2, 0.2597d2, 0.2463d2,&
0.2376d2, 0.2266d2, 0.2206d2, 0.2135d2, 0.2071d2, 0.2042d2,&
0.1903d2, 0.1826d2, 0.1793d2, 0.1774d2, 0.1761d2/
double precision :: xlnha(1:13,1:11)
data xlnha/0.323d1, 0.326d1, 0.329d1, 0.349d1, 0.367d1, 0.397d1,&
0.430d1, 0.446d1, 0.517d1, 0.547d1, 0.553d1, 0.530d1,&
0.470d1, 0.319d1, 0.324d1, 0.327d1, 0.348d1, 0.366d1,&
0.396d1, 0.430d1, 0.445d1, 0.516d1, 0.547d1, 0.553d1,&
0.530d1, 0.470d1, 0.305d1, 0.316d1, 0.322d1, 0.345d1,&
0.364d1, 0.394d1, 0.429d1, 0.445d1, 0.516d1, 0.547d1,&
0.553d1, 0.530d1, 0.470d1, 0.273d1, 0.297d1, 0.307d1,&
0.334d1, 0.356d1, 0.389d1, 0.426d1, 0.442d1, 0.515d1,&
0.546d1, 0.552d1, 0.530d1, 0.470d1, 0.226d1, 0.257d1,&
0.272d1, 0.309d1, 0.335d1, 0.374d1, 0.416d1, 0.434d1,&
0.513d1, 0.545d1, 0.551d1, 0.529d1, 0.468d1, 0.176d1,&
0.207d1, 0.224d1, 0.265d1, 0.295d1, 0.342d1, 0.392d1,&
0.414d1, 0.506d1, 0.541d1, 0.548d1, 0.526d1, 0.464d1,&
0.126d1, 0.157d1, 0.174d1, 0.215d1, 0.247d1, 0.295d1,&
0.349d1, 0.374d1, 0.486d1, 0.530d1, 0.539d1, 0.517d1,&
0.453d1, 0.760d0, 0.107d1, 0.124d1, 0.165d1, 0.197d1,&
0.245d1, 0.300d1, 0.325d1, 0.447d1, 0.502d1, 0.516d1,&
0.494d1, 0.427d1, 0.260d0, 0.573d0, 0.742d0, 0.115d1,&
0.147d1, 0.195d1, 0.250d1, 0.275d1, 0.398d1, 0.457d1,&
0.473d1, 0.452d1, 0.384d1,-0.240d0, 0.727d-1, 0.242d0,&
0.652d0, 0.966d0, 0.145d1, 0.200d1, 0.225d1, 0.348d1,&
0.407d1, 0.424d1, 0.403d1, 0.335d1,-0.740d0,-0.427d0,&
-0.258d0, 0.152d0, 0.466d0, 0.954d0, 0.150d1, 0.175d1,&
0.298d1, 0.357d1, 0.374d1, 0.353d1, 0.285d1/
double precision :: alphaa(1:13,1:11)
data alphaa/0.514d0, 0.465d0, 0.439d0, 0.409d0, 0.392d0, 0.370d0,&
0.361d0, 0.357d0, 0.385d0, 0.437d0, 0.428d0, 0.354d0,&
0.322d0, 0.505d0, 0.460d0, 0.436d0, 0.407d0, 0.391d0,&
0.368d0, 0.359d0, 0.356d0, 0.385d0, 0.437d0, 0.427d0,&
0.354d0, 0.322d0, 0.486d0, 0.448d0, 0.428d0, 0.401d0,&
0.385d0, 0.364d0, 0.356d0, 0.352d0, 0.383d0, 0.436d0,&
0.427d0, 0.352d0, 0.320d0, 0.465d0, 0.433d0, 0.416d0,&
0.388d0, 0.373d0, 0.353d0, 0.347d0, 0.345d0, 0.380d0,&
0.434d0, 0.425d0, 0.349d0, 0.316d0, 0.466d0, 0.448d0,&
0.416d0, 0.378d0, 0.360d0, 0.338d0, 0.332d0, 0.330d0,&
0.371d0, 0.429d0, 0.421d0, 0.341d0, 0.307d0, 0.503d0,&
0.487d0, 0.450d0, 0.396d0, 0.367d0, 0.334d0, 0.322d0,&
0.317d0, 0.355d0, 0.419d0, 0.414d0, 0.329d0, 0.292d0,&
0.566d0, 0.538d0, 0.492d0, 0.435d0, 0.403d0, 0.362d0,&
0.339d0, 0.329d0, 0.343d0, 0.406d0, 0.401d0, 0.317d0,&
0.276d0, 0.598d0, 0.574d0, 0.529d0, 0.473d0, 0.441d0,&
0.404d0, 0.381d0, 0.370d0, 0.362d0, 0.410d0, 0.392d0,&
0.316d0, 0.272d0, 0.603d0, 0.594d0, 0.555d0, 0.503d0,&
0.473d0, 0.440d0, 0.423d0, 0.414d0, 0.418d0, 0.446d0,&
0.404d0, 0.335d0, 0.289d0, 0.613d0, 0.623d0, 0.582d0,&
0.528d0, 0.499d0, 0.469d0, 0.457d0, 0.451d0, 0.470d0,&
0.487d0, 0.432d0, 0.364d0, 0.310d0, 0.634d0, 0.645d0,&
0.596d0, 0.546d0, 0.519d0, 0.492d0, 0.483d0, 0.479d0,&
0.510d0, 0.516d0, 0.448d0, 0.372d0, 0.313d0/


xn_c=(1.d0-y_H2)*xnH

xlT=dlog10(T_K)
cs=1.d-5*dsqrt(2.d0*xk_B*T_K/(28.d0*xm_p))
xNc_CO=xNc_CO+1.d-10
xlN=dlog10(xNc_CO/cs)

call linear(xlTa,aL0a,13,xlT,aL0)
call bilinear(xlTa,xlNa,aLLTEa,13,11,xlT,xlN,aLLTE)
call bilinear(xlTa,xlNa,xlnha,13,11,xlT,xlN,xlnh)
call bilinear(xlTa,xlNa,alphaa,13,11,xlT,xlN,alpha)

xL0inv=10.d0**aL0
xLLTEinv=10.d0**aLLTE
xn_h=10.d0**xlnh
xLinv=xL0inv+xn_c*xLLTEinv&
+xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
xL=1.d0/xLinv

xLd_CO=(1.d0-y_H2)*xL*dexp(-tau_cnt)
END SUBROUTINE COcool



FUNCTION Gamma_CR(zeta,y_a,y_m,y_He)
!     calculates heating rate (ergs/s/g) due to CR
IMPLICIT REAL*8(a-h,o-z)
double precision :: zeta,y_a,y_m,y_He
Gamma_CR=3.26d12*(0.46d0*y_a+0.50d0*y_He+0.94d0*y_m)*zeta/(1.d0+4.d0*y_He)
!return
END FUNCTION Gamma_CR


SUBROUTINE grtemp(xnH,T,esc_cnt,T_rad,xJ_ex,T_gr)
IMPLICIT REAL*8(a-h,o-z)
double precision :: xnH,T,esc_cnt,T_rad,xJ_ex,T_gr
double precision, PARAMETER :: yHe=8.333d-2,xk_vis=40.d0
double precision :: xm_p=1.67d-24

funcL(T_d)=4.9d-2*vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)&
*dsqrt(T*1.d-3)*(0.354+0.5d0*yHe)*(T-T_d)
!      funcL(T_d)=4.9d-2*vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
!     &     *dsqrt(T*1.d-3)*(1.d0-0.8d0*dexp(-75.d0/T))*(T-T_d)
!     for stellar radiation
!      funcD(T_d)=vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)/
!     &     vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,0.d0)
func(T_d)=xkp_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)&
*(T_d**4)*esc_cnt &
!     collision with gas particles
     -xnH*funcL(T_d) &
!     stellar radiation
!     &     -4.41d3*xk_vis*funcD(T_d)*xJ_ex
!     CMB
     -xkp_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d) &
     *(T_rad**4)*esc_cnt

if(esc_cnt.eq.0.d0) then
   T_gr=T
   go to 13
endif

x_l=T_gr
x_u=T_gr
if(func(T_gr).gt.0.d0) then
   do i=1,1000
      x_l=0.99d0*x_l
      if(func(x_l).lt.0.d0) go to 10
   enddo
elseif(func(T_gr).lt.0.d0) then
   do i=1,1000
      x_u=1.01d0*x_u
      if(func(x_u).gt.0.d0) go to 10
   enddo
else
   T_gr=T
   return
endif

10      continue

x1=x_u
x2=x_l

fmid=func(x2)
f=func(x1)
if(f.eq.0.d0) then
   T_gr=x1
   go to 13
elseif(f*fmid.ge.0.d0) then
   pause 'root must be bracketed in T_gr'
endif
if(f.lt.0.d0)then
  T_gr=x1
  dx=x2-x1
else
  T_gr=x2
  dx=x1-x2
endif
do j=1,100
  dx=dx*.5d0
  xmid=T_gr+dx
  fmid=func(xmid)
  if(fmid.le.0.d0)T_gr=xmid
  if(T_gr.ne.0.d0) then
     if(abs(dx/T_gr).lt.1.d-12 .or. fmid.eq.0.d0) then
        go to 13
     endif
  else
     if(fmid.eq.0.d0) then
        go to 13
     endif
  endif
enddo

pause 'too many bisections in T_gr'
13   continue
if(T_gr > 1.d3) then
   rho=(1.d0+4.d0*yHe)*xm_p*xnH
   call vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)
   if(T_gr > T_ol) then
      T_gr=T
      return
   endif
else
   return
endif
END SUBROUTINE grtemp



FUNCTION xkp_gr(rho,T)
IMPLICIT REAL*8(a-h,o-z)
double precision :: rho,T
!     Planck mean mass absorption coefficient ("opacity")
!     owing to Dust Grain (solar metalicity)
!     USES vaptemp,linear
double precision :: Ta(18),xka(18)
DATA (Ta(i),i=1,3)/0.d0,50.d0,100.d0/,&
(Ta(i),i=6,10)/150.d0,200.d0,250.d0,300.d0,350.d0/
DATA (xka(i),i=1,3)/0.d0,1.d0,2.9d0/,&
(xka(i),i=6,18)/3.8d0,4.7d0,5.d0,5.25d0,5.3d0,5.3d0,&
4.3d0,4.3d0,0.9d0,0.9d0,0.4d0,0.4d0,0.d0/

call vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)

Ta(4)=0.975d0*T_ice
Ta(5)=1.025d0*T_ice
Ta(11)=0.975d0*T_vo
Ta(12)=1.025d0*T_vo
Ta(13)=0.975d0*T_ro
Ta(14)=1.025d0*T_ro
Ta(15)=0.975d0*T_tr
Ta(16)=1.025d0*T_tr
Ta(17)=0.975d0*T_ir
Ta(18)=1.025d0*T_ir

xka(4)=1.9d0*(Ta(4)/50.d0)-0.9d0
xka(5)=1.6d0*(Ta(5)/50.d0)-1.d0

if(T.le.Ta(2)) then
   xkp_gr=xka(2)*(T/Ta(2))**2
elseif(T.le.Ta(18)) then
   call linear(Ta,xka,18,T,xkp_gr)
else
   xkp_gr=0.d0
endif
return
END FUNCTION xkp_gr



FUNCTION vol_gr(rho,T_K)
IMPLICIT REAL*8(a-h,o-z)
!     gives grain volume per unit mass of gas (solar metalicity)
double precision :: rho,T_K
!     USES vaptemp
x(t,t1)=dmin1(1.d0,ddim(1.025d0,t/t1)/0.05d0)

call vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)

vol_gr=(0.78d-4-0.46d-4*x(T_K,T_tr))*x(T_K,T_ir)&
+7.19d-4*x(T_K,T_ol)+2.16d-4*x(T_K,T_pyr)&
+1.18d-4*x(T_K,T_tr)+23.53d-4*x(T_K,T_ro)&
+6.02d-4*x(T_K,T_vo)+12.93d-4*x(T_K,T_ice)

return
END FUNCTION vol_gr



SUBROUTINE vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)
!     Vaporization Temperatures (K) of Grain Species
!     T_ice   for   Water ice
!     T_vo    for   Volatile organics
!     T_ro    for   Refractory organics
!     T_tr    for   Troilite
!     T_ir    for   Metallic iron
!     T_pyr   for   Orthopyroxene
!     T_ol    for   Olivine
IMPLICIT REAL*8(a-h,o-z)
double precision :: rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol
!     USES linear
double precision :: xlg_rhoa(13),T_icea(13),T_ira(13),&
T_pyra(13),T_ola(13)
DATA (xlg_rhoa(i),i=1,13)&
/-24.d0,-22.d0,-20.d0,-18.d0,-16.d0,-14.d0,&
-12.d0,-10.d0,-8.d0,-6.d0,-4.d0,-2.d0,0.d0/
DATA (T_icea(i),i=1,13)&
/85.d0,91.d0,98.d0,106.d0,115.d0,125.d0,&
138.d0,153.d0,172.d0,197.d0,230.d0,271.d0,320.d0/
DATA (T_ira(i),i=1,13)&
/694.d0,728.d0,775.d0,835.d0,908.d0,994.d0,&
1100.d0,1230.d0,1395.d0,1612.d0,1908.d0,2283.d0,2737.d0/
DATA (T_pyra(i),i=1,13)&
/794.d0,827.d0,867.d0,920.d0,980.d0,1049.d0,&
1129.d0,1222.d0,1331.d0,1462.d0,1621.d0,1808.d0,2023.d0/
DATA (T_ola(i),i=1,13)&
/791.d0,826.d0,872.d0,929.d0,997.d0,1076.d0,&
1168.d0,1277.d0,1408.d0,1570.d0,1774.d0,2020.d0,2308.d0/

xlg_rho=dlog10(rho)
call linear(xlg_rhoa,T_icea,13,xlg_rho,T_ice)
T_vo=375.d0
T_ro=575.d0
T_tr=680.d0
call linear(xlg_rhoa,T_ira,13,xlg_rho,T_ir)
call linear(xlg_rhoa,T_pyra,13,xlg_rho,T_pyr)
call linear(xlg_rhoa,T_ola,13,xlg_rho,T_ol)
!return

END SUBROUTINE vaptemp


SUBROUTINE phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)
IMPLICIT REAL*8(a-h,o-z)
double precision :: xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe
double precision :: xm_p=1.67d-24, T_ro=575.d0
if(T_gr_K > T_ro) then
   Gmm_pe=0.d0
   return
endif

yHe=8.333d-2
rho=(1.d0+4.d0*yHe)*xm_p*xnH
x=G_0*dsqrt(T_K)/(y_e*xnH)
eps=4.9d-2/(1.d0+4.0d-3*(x**0.73d0))&
+3.7d-2*(T_K*1.d-4)**0.7d0/(1.d0+2.0d-4*x)

Gam_pe=1.0d-24*eps*G_0*dexp(-1.8d0*A_v)

beta=0.74d0/(T_K**0.068d0)
xLd_pe=4.65d-30*(T_K**0.94d0)*(x**beta)*y_e

Gmm_pe=Z_metal*(Gam_pe*xnH-xLd_pe*xnH**2)/rho

!return
END SUBROUTINE phelectr


FUNCTION c_H2(T_K)
IMPLICIT REAL*8(a-h,o-z)
double precision :: T_K,c_H2,E_T,T_K_b,T_K_f,x
if(T_K.gt.1.d3) then
   c_rot=1.d0
   go to 10
endif

eps=1.d-3
T_K_b=(1.d0-eps)*T_K
T_K_f=(1.d0+eps)*T_K
Zp_b=0.d0
Zp_f=0.d0
Xp_b=0.d0
Xp_f=0.d0
Zo_b=0.d0
Zo_f=0.d0
Xo_b=0.d0
Xo_f=0.d0
do K=0,20
   E_T=85.4d0*dble(K*(K+1))
   if(mod(K,2).eq.0) then
      dZp_b=dble(2*K+1)*dexp(-E_T/T_K_b)
      dZp_f=dble(2*K+1)*dexp(-E_T/T_K_f)
      dXp_b=E_T*dZp_b
      dXp_f=E_T*dZp_f
      Zp_b=Zp_b+dZp_b
      Zp_f=Zp_f+dZp_f
      Xp_b=Xp_b+dXp_b
      Xp_f=Xp_f+dXp_f
   else
      dZo_b=dble(2*K+1)*dexp(-E_T/T_K_b)
      dZo_f=dble(2*K+1)*dexp(-E_T/T_K_f)
      dXo_b=E_T*dZo_b
      dXo_f=E_T*dZo_f
      Zo_b=Zo_b+dZo_b
      Zo_f=Zo_f+dZo_f
      Xo_b=Xo_b+dXo_b
      Xo_f=Xo_f+dXo_f
   endif
enddo
Ep_rot_f=Xp_f/Zp_f
Ep_rot_b=Xp_b/Zp_b
Eo_rot_f=Xo_f/Zo_f
Eo_rot_b=Xo_b/Zo_b
E_rot_f=0.25d0*Ep_rot_f+0.75d0*Eo_rot_f
E_rot_b=0.25d0*Ep_rot_b+0.75d0*Eo_rot_b
c_rot=(E_rot_f-E_rot_b)/(2.d0*eps*T_K)

10   continue

 x=6.1d3/T_K
 if(x.gt.1.d2) then
    c_vib=0.d0
 else
    c_vib=(x**2)*dexp(x)/(dexp(x)-1)**2
 endif

 c_H2=1.5d0+c_rot+c_vib

!return
END FUNCTION c_H2


!********************************************************************
!********************************************************************
!**********************  MOLECULAR HYDROGEN  ************************
!********************************************************************
!********************************************************************

function H2_pf(T_K1)
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: H2_pf,T_K1
integer, parameter :: iv_max=5!,j_max=25
double precision :: ET(0:iv_max,0:j_max),ET_H2_BFM(0:iv_max,0:j_max)
do iv=0,iv_max
   do j=0,j_max
      ET(iv,j)=E_H2_BFM(iv,j)
   enddo
enddo

ET_diss=5.196d4

z_p=0.d0
z_o=0.d0
do iv=0,iv_max
   do j=0,j_max
      if(ET(iv,j)-ET(0,0).le.ET_diss) then
         if(mod(j,2).eq.0) then
            z_p=z_p+dble(2*j+1)*dexp(-(ET(iv,j)-ET(0,0))/T_K1)
         else
            z_o=z_o+dble(2*j+1)*dexp(-(ET(iv,j)-ET(0,1))/T_K1)
         endif
      endif
   enddo
enddo

H2_pf=0.25d0*z_p+0.75d0*z_o
!return
end function H2_pf


function E_H2_BFM(iv,J)
IMPLICIT REAL*8(a-h,o-z)
!     Borysow, Frommhold, and Moraldi (1989) ApJ,336,495
!     Equation (A11) (in K)
double precision :: E_H2_BFM
integer :: iv,J

Ev0=0.38496d0
Ev1=-0.04609d0
Ev2=0.00178d0
Ev3=-7.d-5
Ev4=2.9511d-6

Bv0=54.438d0
Bv1=4.6063d0
Bv2=-2.0050d0
Bv3=0.19260d0
Bv4=-6.2953d-3

Dv0=-0.72593d0
Dv1=5.9990d0
Dv2=-1.5187d0
Dv3=0.12721d0
Dv4=-3.3391d-3

Fv0=-12.662d0
Fv1=20.047d0
Fv2=-4.7873d0
Fv3=0.36900d0
Fv4=-9.003d-3

Gv0=-24.006d0
Gv1=33.989d0
Gv2=-7.6841d0
Gv3=0.52413d0
Gv4=-1.1297d-2

Hv0=-22.384d0
Hv1=31.150d0
Hv2=-6.6139d0
Hv3=0.39226d0
Hv4=-9.496d-3

Ov0=-10.541d0
Ov1=14.746d0
Ov2=-2.9476d0
Ov3=0.16016d0
Ov4=-6.5005d-3

Pv0=-2.0021d0
Pv1=2.8400d0
Pv2=-0.54654d0
Pv3=0.031636d0
Pv4=-2.4398d-3

vv=dble(iv)+0.5d0
Ev=Ev0+Ev1*vv+Ev2*vv**2+Ev3*vv**3+Ev4*vv**4
Bv=Bv0+Bv1*vv+Bv2*vv**2+Bv3*vv**3+Bv4*vv**4
Dv=Dv0+Dv1*vv+Dv2*vv**2+Dv3*vv**3+Dv4*vv**4
Fv=Fv0+Fv1*vv+Fv2*vv**2+Fv3*vv**3+Fv4*vv**4
Gv=Gv0+Gv1*vv+Gv2*vv**2+Gv3*vv**3+Gv4*vv**4
Hv=Hv0+Hv1*vv+Hv2*vv**2+Hv3*vv**3+Hv4*vv**4
Ov=Ov0+Ov1*vv+Ov2*vv**2+Ov3*vv**3+Ov4*vv**4
Pv=Pv0+Pv1*vv+Pv2*vv**2+Pv3*vv**3+Pv4*vv**4

rrot=dble(J*(J+1))
E_H2_BFM=1.43879d0*&
(-Ev*1.d5+Bv*rrot-Dv*1.d-2*rrot**2+Fv*1.d-5*rrot**3&
-Gv*1.d-8*rrot**4+Hv*1.d-11*rrot**5-Ov*1.d-14*rrot**6&
+Pv*1.d-17*rrot**7)

!return
end function E_H2_BFM


FUNCTION xLd_H2(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt)
IMPLICIT REAL*8(a-h,o-z)
!     xnH*xn(H2)*xLd_H2 : cooling rate owing to H2 per unit volume
double precision :: xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt,xLd_H2,Q_bg
!     USES beta_esc, Q_bg
double precision :: A(0:2,0:2,0:22,0:22),ET(0:2,0:22),p(0:22),&
f(0:2,0:22),f_v(0:2)
double precision :: xk_B=1.380662d-16,h_P=6.626176d-27,pi=3.14159265358979d0,&
xm_p=1.67d-24,c_light=2.99792458d10,beta_esc,esc
data (((A(i,ii,j,j+2),ii=0,i-1),i=1,2),j=0,20)&
/8.54d-7 ,3.47d-7 ,1.29d-6&
,4.23d-7 ,1.61d-7 ,6.40d-7&
,2.90d-7 ,1.03d-7 ,4.41d-7&
,2.09d-7 ,6.98d-8 ,3.18d-7&
,1.50d-7 ,4.72d-8 ,2.28d-7&
,1.06d-7 ,3.15d-8 ,1.62d-7&
,7.38d-8 ,2.05d-8 ,1.12d-7&
,4.98d-8 ,1.31d-8 ,7.50d-8&
,3.27d-8 ,8.07d-9 ,4.88d-8&
,2.07d-8 ,4.83d-9 ,3.06d-8&
,1.27d-8 ,2.79d-9 ,1.85d-8&
,7.44d-9 ,1.55d-9 ,1.07d-8&
,4.16d-9 ,8.27d-10,5.84d-9&
,2.19d-9 ,4.21d-10,3.00d-9&
,1.08d-9 ,2.04d-10,1.43d-9&
,4.87d-10,9.45d-11,6.17d-10&
,1.96d-10,4.23d-11,2.33d-10&
,6.71d-11,1.91d-11,7.32d-11&
,1.82d-11,9.63d-12,1.73d-11&
,3.38d-12,6.26d-12,2.46d-12&
,2.96d-13,5.78d-12,1.16d-13/

data (((A(i,ii,j,j),ii=0,i-1),i=1,2),j=1,20)&
/4.29d-7 ,1.94d-7 ,6.37d-7&
,3.03d-7 ,1.38d-7 ,4.50d-7&
,2.78d-7 ,1.29d-7 ,4.12d-7&
,2.65d-7 ,1.25d-7 ,3.91d-7&
,2.55d-7 ,1.23d-7 ,3.74d-7&
,2.45d-7 ,1.21d-7 ,3.58d-7&
,2.34d-7 ,1.20d-7 ,3.40d-7&
,2.23d-7 ,1.18d-7 ,3.22d-7&
,2.12d-7 ,1.17d-7 ,3.03d-7&
,1.99d-7 ,1.15d-7 ,2.84d-7&
,1.87d-7 ,1.13d-7 ,2.63d-7&
,1.74d-7 ,1.11d-7 ,2.42d-7&
,1.61d-7 ,1.08d-7 ,2.21d-7&
,1.47d-7 ,1.05d-7 ,2.01d-7&
,1.34d-7 ,1.02d-7 ,1.80d-7&
,1.21d-7 ,9.88d-8 ,1.60d-7&
,1.08d-7 ,9.50d-8 ,1.41d-7&
,9.61d-8 ,9.07d-8 ,1.22d-7&
,8.43d-8 ,8.62d-8 ,1.05d-7&
,7.32d-8 ,8.14d-8 ,8.85d-8/

data (((A(i,ii,j,j-2),ii=0,i),i=0,2),j=2,20)&
/2.94d-11,2.53d-7 ,2.79d-11&
,1.27d-7 ,3.68d-7 ,2.56d-11&
,4.76d-10,3.47d-7,4.50d-10&
,1.90d-7 ,4.98d-7 ,4.12d-10&
,2.76d-9 ,3.98d-7 ,2.59d-9&
,2.38d-7 ,5.60d-7 ,2.37d-9&
,9.84d-9 ,4.21d-7 ,9.21d-9&
,2.77d-7 ,5.77d-7 ,8.37d-9&
,2.64d-8 ,4.19d-7 ,2.46d-8&
,3.07d-7 ,5.57d-7 ,2.22d-8&
,5.88d-8 ,3.96d-7 ,5.44d-8&
,3.28d-7 ,5.05d-7 ,4.88d-8&
,1.14d-7 ,3.54d-7 ,1.05d-7&
,3.39d-7 ,4.30d-7 ,9.33d-8&
,2.00d-7 ,2.98d-7 ,1.82d-7&
,3.40d-7 ,3.38d-7 ,1.61d-7&
,3.24d-7 ,2.34d-7 ,2.92d-7&
,3.30d-7 ,2.41d-7 ,2.55d-7&
,4.90d-7 ,1.68d-7 ,4.38d-7&
,3.12d-7 ,1.49d-7 ,3.79d-7&
,7.03d-7 ,1.05d-7 ,6.21d-7&
,2.85d-7 ,7.20d-8 ,5.32d-7&
,9.64d-7 ,5.30d-8 ,8.42d-7&
,2.53d-7 ,1.96d-8 ,7.13d-7&
,1.27d-6 ,1.65d-8 ,1.10d-6&
,2.16d-7 ,1.49d-11,9.18d-7&
,1.62d-6 ,4.26d-10,1.38d-6&
,1.78d-7 ,1.91d-8 ,1.14d-6&
,2.00d-6 ,8.38d-9 ,1.69d-6&
,1.39d-7 ,8.07d-8 ,1.37d-6&
,2.41d-6 ,4.27d-8 ,2.00d-6&
,1.01d-7 ,1.86d-7 ,1.61d-6&
,2.83d-6 ,1.04d-7 ,2.32d-6&
,6.80d-8 ,3.35d-7 ,1.84d-6&
,3.26d-6 ,1.93d-7 ,2.64d-6&
,3.98d-8 ,5.24d-7 ,2.05d-6&
,3.68d-6 ,3.08d-7 ,2.93d-6&
,1.84d-8 ,7.49d-7 ,2.23d-6/

g_0=1.d0
g_1=1.d0
g_2=1.d0

A_10=8.3d-7
A_20=4.1d-7
A_21=1.1d-6

!     Hollenbach & McKee (1979) corrected with (1989)
gamma_10H=1.0d-12*dsqrt(T_K)*dexp(-1.d3/T_K)
gamma_20H=1.6d-12*dsqrt(T_K)*dexp(-(4.d2/T_K)**2)
gamma_21H=4.5d-12*dsqrt(T_K)*dexp(-(5.d2/T_K)**2)

gamma_10H2=1.4d-12*dsqrt(T_K)*dexp(-1.81d4/(T_K+1.2d3))
gamma_20H2=0.d0
gamma_21H2=gamma_10H2

do iv=0,2
   do j=0,22
      ET(iv,j)=E_H2_BFM(iv,J)
   enddo
enddo

DT_10=ET(1,0)-ET(0,0)
DT_20=ET(2,0)-ET(0,0)
DT_21=ET(2,0)-ET(1,0)

!     Draine, Reberge, Dalgarno (1983)
gamma_10e=3.7d-11*(T_K**0.5d0)/(1.d0+0.5d0*DT_10/T_K)
gamma_20e=2.5d-12*(T_K**0.5d0)/(1.d0+0.5d0*DT_20/T_K)
gamma_21e=3.7d-11*(T_K**0.5d0)/(1.d0+0.5d0*DT_21/T_K)

!     Krstic (2002), my fit
gamma_10Hp=1.4d-4*T_K**(-1.344d0)&
*(1.d0+4.005d-9*T_K**2.066d0)&
*dexp((DT_10-9589d0)/T_K)
gamma_20Hp=4.585d-5*T_K**(-1.291d0)&
*(1.d0+1.378d-8*T_K**1.903d0)&
*dexp((DT_20-1.933d4)/T_K)
gamma_21Hp=5.593d-3*T_K**(-1.770d0)&
*(1.d0+3.505d-10*T_K**2.406d0)&
*dexp((DT_21-1.460d4)/T_K)

xn_a=y_H*xnH
xn_m=y_H2*xnH
xn_e=y_e*xnH
xn_Hp=y_Hp*xnH

C_10=gamma_10H*xn_a+gamma_10H2*xn_m&
+gamma_10e*xn_e+gamma_10Hp*xn_Hp
C_20=gamma_20H*xn_a+gamma_20H2*xn_m&
+gamma_20e*xn_e+gamma_20Hp*xn_Hp
C_21=gamma_21H*xn_a+gamma_21H2*xn_m&
+gamma_21e*xn_e+gamma_21Hp*xn_Hp

C_01=C_10*dexp(-DT_10/T_K)
C_02=C_20*dexp(-DT_20/T_K)
C_12=C_21*dexp(-DT_21/T_K)

Q_10=Q_bg(DT_10)
Q_20=Q_bg(DT_20)
Q_21=Q_bg(DT_21)

R_10=A_10*(1.d0+Q_10)+C_10
R_20=A_20*(1.d0+Q_20)+C_20
R_21=A_21*(1.d0+Q_21)+C_21
R_01=(g_1/g_0)*A_10*Q_10+C_01
R_02=(g_2/g_0)*A_20*Q_20+C_02
R_12=(g_2/g_1)*A_21*Q_21+C_12

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )&
/( (R_01+R_02+R_20)*(R_10+R_12+R_21)&
-(R_01-R_21)*(R_10-R_20) )
f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

f_v(0)=f_0
f_v(1)=f_1
f_v(2)=f_2


f_para=0.25d0
f_ortho=0.75d0
do iv=0,2
   z_para=0.d0
   z_ortho=0.d0
   p(0)=1.d0
   p(1)=1.d0
   do j=0,18
      if(mod(j,2).eq.0) then
         z_para=z_para+p(j)
      else
         z_ortho=z_ortho+p(j)
      endif

!     H2-H collision
!            Galli & Palla (1998)
      if(j == 0) then
         gamma_H=2.93d-14+1.21d-15*T_K+2.16d-19*T_K**2+1.32d-21*T_K**3
      elseif(j == 1) then
         gamma_H=8.34d-14+5.97d-16*T_K+7.76d-19*T_K**2+1.72d-21*T_K**3
      elseif(j == 2) then
         gamma_H=7.54d-14+2.65d-16*T_K+2.14d-19*T_K**2+2.65d-21*T_K**3
      elseif(j == 3) then
         gamma_H=2.95d-14+1.21d-16*T_K+5.53d-19*T_K**2+2.51d-21*T_K**3
      else
!     HM(1989)
         gamma_H=4.6d-12*(2.d0*dble(j+2)-3.d0)*dsqrt(T_K)&
         *dsqrt(1.d0+2.d0*85.25d0*(2.d0*dble(j+2)-1.d0)/T_K)&
         *dexp(-(1.d1*85.25d0*(2.d0*dble(j+2)-1.d0))&
         /(T_K+85.25d0*dble(j+2)*dble(j+3))&
         -0.1187d0*(4.d0*dble(j+2)-2.d0))
      endif

!     H2-H2 collision
      gamma_H2=(3.3d-12+6.6d-15*T_K)&
      *0.276d0*dble((j+2)**2)&
      *dexp(-(dble(j+2)/3.18d0)**1.7d0)

!     H2-Hp collision : Gerlich(1990)
xmeV=11.604d0
if(j == 0) then
!        2-0
   gamma_Hp=2.889d-10*dexp(0.171d0*xmeV/T_K)
elseif(j == 1) then
!        3-1
   gamma_Hp=8.202d-10*dexp(0.250d0*xmeV/T_K)
elseif(j == 2) then
!        4-2
   gamma_Hp=5.700d-10*dexp(0.159d0*xmeV/T_K)
elseif(j == 3) then
!        5-3
   gamma_Hp=8.883d-10*dexp(0.223d0*xmeV/T_K)
elseif(j == 4) then
!        6-4
   gamma_Hp=5.209d-10*dexp(0.124d0*xmeV/T_K)
elseif(j == 5) then
!        7-5
   gamma_Hp=7.977d-10*dexp(0.177d0*xmeV/T_K)
elseif(j == 6) then
!        8-6
   gamma_Hp=4.489d-10*dexp(0.101d0*xmeV/T_K)
else
!        9-7
   gamma_Hp=5.915d-10*dexp(0.094d0*xmeV/T_K)
endif

!     H2-e collisional rate from Draine(1990)
   xj=dble(j)
   x=(ET(iv,j+2)-ET(iv,j))/T_K
   gamma_e=1.d-10*(xj+2.d0)*(xj+1.d0)/(2.d0*xj+5.d0)*(1.d0+1.5d0/x)

      C_ul=gamma_H*xn_a+gamma_H2*xn_m+gamma_Hp*xn_Hp+gamma_e*xn_e
      DT=ET(iv,j+2)-ET(iv,j)
      Q_ul=Q_bg(DT)
      g_u=2.d0*dble(j)+5.d0
      g_l=2.d0*dble(j)+1.d0
      r=(g_u/g_l)*(A(iv,iv,j+2,j)*Q_ul+C_ul*dexp(-DT/T_K))/(A(iv,iv,j+2,j)*(1.d0+Q_ul)+C_ul)
      p(j+2)=p(j)*r
   enddo
   
   do j=0,20
      if(mod(j,2).eq.0) then
         f(iv,j)=(p(j)/z_para)*f_para*f_v(iv)
      else
         f(iv,j)=(p(j)/z_ortho)*f_ortho*f_v(iv)
      endif
   enddo
enddo

v_th=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))

xLd_H2=0.d0
do ji=0,18
   do ivi=1,2
      do ivf=0,ivi-1
         jf=ji+2
         DT=ET(ivi,ji)-ET(ivf,jf)
         DE=DT*xk_B
         xnu_Hz=DE/h_P
         xNc_l=xNc_H2*f(ivf,jf)
         xNc_u=xNc_H2*f(ivi,ji)
         g_l=dble(2*jf+1)
         g_u=dble(2*ji+1)
         Q_ul=Q_bg(DT)
         tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3&
         *(xNc_l*g_u/g_l-xNc_u)/v_th
         esc=beta_esc(tau_ul,tau_cnt)
         if(f(ivi,ji).ne.0.d0) then
            xLd_H2=xLd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc*(1.d0-Q_ul*((g_u/g_l)&
            *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
         endif
      enddo
   enddo
enddo

do ji=1,20
   do ivi=1,2
      do ivf=0,ivi-1
         jf=ji
         DT=ET(ivi,ji)-ET(ivf,jf)
         DE=DT*xk_B
         xnu_Hz=DE/h_P
         xNc_l=xNc_H2*f(ivf,jf)
         xNc_u=xNc_H2*f(ivi,ji)
         g_l=dble(2*jf+1)
         g_u=dble(2*ji+1)
         Q_ul=Q_bg(DT)
         tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3*(xNc_l*g_u/g_l-xNc_u)/v_th
         esc=beta_esc(tau_ul,tau_cnt)
         if(f(ivi,ji).ne.0.d0) then
            xLd_H2=xLd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc*(1.d0-Q_ul*((g_u/g_l)&
            *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
         endif
      enddo
   enddo
enddo

do ji=2,20
   do ivi=0,2
      do ivf=0,ivi
         jf=ji-2
         DT=ET(ivi,ji)-ET(ivf,jf)
         DE=DT*xk_B
         xnu_Hz=DE/h_P
         xNc_l=xNc_H2*f(ivf,jf)
         xNc_u=xNc_H2*f(ivi,ji)
         g_l=dble(2*jf+1)
         g_u=dble(2*ji+1)
         Q_ul=Q_bg(DT)
         tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3*(xNc_l*g_u/g_l-xNc_u)/v_th
         esc=beta_esc(tau_ul,tau_cnt)
         if(f(ivi,ji).ne.0.d0) then
            xLd_H2=xLd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc*(1.d0-Q_ul*((g_u/g_l)&
            *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
         endif
      enddo
   enddo
enddo
!return
END FUNCTION xLd_H2


SUBROUTINE H2Ocool(xnH,T_K,y_H2,xNc_H2O,tau_cnt,xLd_H2O)
IMPLICIT REAL*8(a-h,o-z)
!c     H2O cooling function by Neufeld et al.
double precision :: xnH,T_K,y_H2,xNc_H2O,tau_cnt,xLd_H2O
double precision :: xk_B=1.38d-16, xm_p=1.67d-24
!c     100, 200, 400, 1000, 2000, 4000K
double precision ::  xlTa(1:6)
data xlTa/2.000d0, 2.301d0, 2.602d0,&
3.000d0, 3.301d0, 3.602d0/
double precision :: xlNa(1:10)
data xlNa/10.0d0, 11.0d0, 12.0d0, 13.0d0, 14.0d0,&
15.0d0, 16.0d0, 17.0d0, 18.0d0, 19.0d0/
double precision :: aL0a(1:6)
data aL0a/24.35d0, 23.87d0, 23.42d0,&
22.88d0, 22.50d0, 22.14d0/
double precision :: aLLTEa(1:6,1:10)
data aLLTEa/14.59d0, 13.85d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,&
14.59d0, 13.86d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,&
14.60d0, 13.86d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,&
14.68d0, 13.88d0, 13.17d0, 12.32d0, 11.86d0, 11.64d0,&
14.98d0, 14.05d0, 13.25d0, 12.34d0, 11.87d0, 11.65d0,&
15.53d0, 14.46d0, 13.53d0, 12.49d0, 11.97d0, 11.72d0,&
16.22d0, 15.05d0, 14.02d0, 12.87d0, 12.35d0, 12.06d0,&
17.00d0, 15.74d0, 14.63d0, 13.46d0, 12.97d0, 12.66d0,&
17.83d0, 16.50d0, 15.32d0, 14.16d0, 13.69d0, 13.36d0,&
18.70d0, 17.31d0, 16.07d0, 14.94d0, 14.46d0, 14.13d0/
double precision :: xlnha(1:6,1:10)
data xlnha/9.00d0, 9.04d0, 9.19d0, 9.50d0, 9.67d0, 9.60d0,&
8.99d0, 9.04d0, 9.19d0, 9.50d0, 9.67d0, 9.60d0,&
8.96d0, 9.03d0, 9.19d0, 9.50d0, 9.66d0, 9.59d0,&
8.74d0, 8.89d0, 9.11d0, 9.47d0, 9.65d0, 9.59d0,&
8.11d0, 8.37d0, 8.73d0, 9.31d0, 9.56d0, 9.53d0,&
7.20d0, 7.51d0, 7.95d0, 8.74d0, 9.15d0, 9.20d0,&
6.22d0, 6.53d0, 6.99d0, 7.87d0, 8.38d0, 8.50d0,&
5.22d0, 5.57d0, 6.03d0, 6.94d0, 7.48d0, 7.64d0,&
4.24d0, 4.59d0, 5.09d0, 6.02d0, 6.59d0, 6.78d0,&
3.21d0, 3.58d0, 4.10d0, 5.08d0, 5.69d0, 5.89d0/
double precision :: alphaa(1:6,1:10)
data alphaa/0.43d0, 0.42d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0,&
0.43d0, 0.42d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0,&
0.42d0, 0.41d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0,&
0.41d0, 0.39d0, 0.37d0, 0.35d0, 0.33d0, 0.33d0,&
0.42d0, 0.38d0, 0.34d0, 0.33d0, 0.32d0, 0.32d0,&
0.45d0, 0.38d0, 0.34d0, 0.32d0, 0.30d0, 0.30d0,&
0.47d0, 0.40d0, 0.35d0, 0.32d0, 0.29d0, 0.30d0,&
0.50d0, 0.42d0, 0.36d0, 0.32d0, 0.28d0, 0.29d0,&
0.52d0, 0.44d0, 0.37d0, 0.31d0, 0.27d0, 0.28d0,&
0.53d0, 0.45d0, 0.39d0, 0.31d0, 0.27d0, 0.27d0/

!     10, 20, 30, 50, 80, 100K
double precision :: xlTb(1:6)
data xlTb/1.000d0, 1.301d0, 1.477d0,&
1.699d0, 1.903d0, 2.000d0/
double precision :: aL0bo(1:6)
data aL0bo/26.81d0, 25.88d0, 25.43d0,&
24.96d0, 24.58d0, 24.41d0/
double precision :: aLLTEbo(1:6,1:10)
data aLLTEbo/17.94d0, 16.71d0, 16.08d0, 15.41d0, 14.85d0, 14.60d0,&
17.96d0, 16.72d0, 16.09d0, 15.42d0, 14.86d0, 14.60d0,&
18.14d0, 16.86d0, 16.19d0, 15.47d0, 14.88d0, 14.62d0,&
18.77d0, 17.36d0, 16.58d0, 15.72d0, 15.02d0, 14.73d0,&
19.70d0, 18.11d0, 17.25d0, 16.27d0, 15.47d0, 15.11d0,&
20.67d0, 18.93d0, 18.05d0, 17.01d0, 16.12d0, 15.72d0,&
21.58d0, 19.81d0, 18.91d0, 17.82d0, 16.88d0, 16.45d0,&
22.53d0, 20.71d0, 19.80d0, 18.69d0, 17.70d0, 17.25d0,&
23.50d0, 21.64d0, 20.72d0, 19.58d0, 18.56d0, 18.10d0,&
24.41d0, 22.58d0, 21.65d0, 20.49d0, 19.46d0, 18.99d0/
double precision :: xlnhbo(1:6,1:10)
data  xlnhbo/8.81d0, 8.90d0, 9.03d0, 9.14d0, 9.19d0, 9.20d0,&
8.79d0, 8.88d0, 9.01d0, 9.13d0, 9.20d0, 9.20d0,&
8.60d0, 8.73d0, 8.88d0, 9.03d0, 9.13d0, 9.15d0,&
7.96d0, 8.14d0, 8.31d0, 8.55d0, 8.78d0, 8.85d0,&
7.01d0, 7.21d0, 7.40d0, 7.69d0, 8.00d0, 8.11d0,&
6.02d0, 6.20d0, 6.41d0, 6.71d0, 7.04d0, 7.17d0,&
5.03d0, 5.21d0, 5.42d0, 5.70d0, 6.06d0, 6.19d0,&
4.02d0, 4.21d0, 4.41d0, 4.71d0, 5.05d0, 5.19d0,&
3.02d0, 3.22d0, 3.41d0, 3.71d0, 4.06d0, 4.20d0,&
2.03d0, 2.21d0, 2.42d0, 2.70d0, 3.06d0, 3.20/
double precision :: alphabo(1:6,1:10)
data alphabo/0.71d0, 0.49d0, 0.48d0, 0.46d0, 0.45d0, 0.45d0,&
0.64d0, 0.49d0, 0.48d0, 0.45d0, 0.45d0, 0.45d0,&
0.57d0, 0.50d0, 0.47d0, 0.45d0, 0.44d0, 0.44d0,&
0.56d0, 0.53d0, 0.48d0, 0.44d0, 0.42d0, 0.41d0,&
0.59d0, 0.60d0, 0.53d0, 0.47d0, 0.45d0, 0.43d0,&
0.72d0, 0.64d0, 0.58d0, 0.52d0, 0.49d0, 0.47d0,&
0.85d0, 0.68d0, 0.61d0, 0.56d0, 0.52d0, 0.50d0,&
0.86d0, 0.70d0, 0.63d0, 0.58d0, 0.55d0, 0.53d0,&
0.93d0, 0.72d0, 0.66d0, 0.61d0, 0.57d0, 0.55d0,&
0.87d0, 0.73d0, 0.67d0, 0.62d0, 0.59d0, 0.56/

double precision :: aL0bp(1:6)
data aL0bp/27.01d0, 25.73d0, 25.24d0,&
24.75d0, 24.38d0, 24.22d0/
double precision :: aLLTEbp(1:6,1:10)
data aLLTEbp/17.72d0, 16.60d0, 16.12d0, 15.43d0, 14.86d0, 14.60d0,&
17.76d0, 16.63d0, 16.13d0, 15.43d0, 14.86d0, 14.60d0,&
18.07d0, 16.85d0, 16.24d0, 15.48d0, 14.88d0, 14.62d0,&
18.83d0, 17.41d0, 16.61d0, 15.72d0, 15.03d0, 14.73d0,&
19.68d0, 18.11d0, 17.26d0, 16.28d0, 15.47d0, 15.11d0,&
20.50d0, 18.94d0, 18.06d0, 17.01d0, 16.12d0, 15.72d0,&
21.37d0, 19.83d0, 18.93d0, 17.82d0, 16.87d0, 16.45d0,&
22.28d0, 20.75d0, 19.81d0, 18.69d0, 17.70d0, 17.25d0,&
23.22d0, 21.69d0, 20.73d0, 19.58d0, 18.56d0, 18.10d0,&
24.19d0, 22.60d0, 21.65d0, 20.49d0, 19.45d0, 18.99d0/
double precision :: xlnhbp(1:6,1:10)
data xlnhbp/9.30d0, 9.11d0, 8.94d0, 8.70d0, 8.55d0, 8.50d0,&
9.25d0, 9.06d0, 8.91d0, 8.69d0, 8.54d0, 8.49d0,&
8.94d0, 8.82d0, 8.71d0, 8.56d0, 8.46d0, 8.43d0,&
8.17d0, 8.16d0, 8.16d0, 8.12d0, 8.10d0, 8.11d0,&
7.21d0, 7.24d0, 7.28d0, 7.30d0, 7.32d0, 7.35d0,&
6.20d0, 6.24d0, 6.31d0, 6.33d0, 6.35d0, 6.38d0,&
5.21d0, 5.25d0, 5.30d0, 5.32d0, 5.35d0, 5.39d0,&
4.22d0, 4.26d0, 4.31d0, 4.33d0, 4.36d0, 4.40d0,&
3.23d0, 3.24d0, 3.31d0, 3.33d0, 3.37d0, 3.41d0,&
2.21d0, 2.25d0, 2.30d0, 2.32d0, 2.37d0, 2.41d0/
double precision :: alphabp(1:6,1:10)
data alphabp/0.49d0, 0.72d0, 0.69d0, 0.53d0, 0.46d0, 0.44d0,&
0.52d0, 0.65d0, 0.66d0, 0.53d0, 0.46d0, 0.44d0,&
0.49d0, 0.65d0, 0.63d0, 0.51d0, 0.45d0, 0.43d0,&
0.57d0, 0.73d0, 0.65d0, 0.49d0, 0.43d0, 0.41d0,&
0.75d0, 0.68d0, 0.64d0, 0.52d0, 0.45d0, 0.42d0,&
0.76d0, 0.70d0, 0.68d0, 0.55d0, 0.47d0, 0.43d0,&
0.79d0, 0.73d0, 0.69d0, 0.58d0, 0.48d0, 0.44d0,&
0.80d0, 0.75d0, 0.71d0, 0.58d0, 0.50d0, 0.46d0,&
0.82d0, 0.77d0, 0.73d0, 0.60d0, 0.52d0, 0.48d0,&
0.91d0, 0.80d0, 0.74d0, 0.62d0, 0.53d0, 0.50d0/

xn_c=(1.d0-y_H2)*xnH

xlT=dlog10(T_K)
cs=1.d-5*dsqrt(2.d0*xk_B*T_K/(18.d0*xm_p))
xNc_H2O=xNc_H2O+1.d-10
if(xlT > 2.d0) then
   xlN=dlog10(xNc_H2O/cs)
   call linear(xlTa,aL0a,6,xlT,aL0)
   call bilinear(xlTa,xlNa,aLLTEa,6,10,xlT,xlN,aLLTE)
   call bilinear(xlTa,xlNa,xlnha,6,10,xlT,xlN,xlnh)
   call bilinear(xlTa,xlNa,alphaa,6,10,xlT,xlN,alpha)
   xL0inv=10.d0**aL0
   xLLTEinv=10.d0**aLLTE
   xn_h=10.d0**xlnh
   xLinv=xL0inv+xn_c*xLLTEinv&
   +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
   xL=1.d0/xLinv
else
   xlNo=dlog10(0.75d0*xNc_H2O/cs)
   call linear(xlTb,aL0bo,6,xlT,aL0)
   call bilinear(xlTb,xlNa,aLLTEbo,6,10,xlT,xlNo,aLLTE)
   call bilinear(xlTb,xlNa,xlnhbo,6,10,xlT,xlNo,xlnh)
   call bilinear(xlTb,xlNa,alphabo,6,10,xlT,xlNo,alpha)
   xL0inv=10.d0**aL0
   xLLTEinv=10.d0**aLLTE
   xn_h=10.d0**xlnh
   xLinv=xL0inv+xn_c*xLLTEinv&
   +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
   xLo=1.d0/xLinv

   xlNp=dlog10(0.25d0*xNc_H2O/cs)
   call linear(xlTb,aL0bp,6,xlT,aL0p)
   call bilinear(xlTb,xlNa,aLLTEbp,6,10,xlT,xlNp,aLLTE)
   call bilinear(xlTb,xlNa,xlnhbp,6,10,xlT,xlNp,xlnh)
   call bilinear(xlTb,xlNa,alphabp,6,10,xlT,xlNp,alpha)
   xL0inv=10.d0**aL0
   xLLTEinv=10.d0**aLLTE
   xn_h=10.d0**xlnh
   xLinv=xL0inv+xn_c*xLLTEinv&
   +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
   xLp=1.d0/xLinv

   xL=0.75d0*xLo+0.25d0*xLp
   xL=1.148d0*xL
endif

xLd_H2O=(1.d0-y_H2)*xL*dexp(-tau_cnt)
END SUBROUTINE H2Ocool


FUNCTION xLd_HD(xnH,T_K,y_H,y_H2,xNc_HD,tau_cnt)
IMPLICIT REAL*8(a-h,o-z)
! HD cooling from Galli and Palla (1998), A&A, 335, 403
! cooling is for T < 3000 K, in erg cm^3 s-1
double precision :: xLd_HD,xnH,T_K,y_H,y_H2,xNc_HD,tau_cnt,Q_bg,beta_esc
double precision ::  xk_B=1.38066d-16,h_Pl=6.62618d-27,pi=3.14159265358979d0,xm_p=1.67d-24
!     USES Q_bg, beta_esc
!------------------------------------------------
!     Various parameters for HD
!-----------------------------------------------

g0=1.d0
g1=3.d0
g2=5.d0
g3=7.d0
!     radiative transitions
a10d = 5.12d-8
a21d = 4.86d-7
a32d = 1.72d-6
a20d = 7.05d-12
a31d = 1.15d-10
!     frequencies in Hz
xnu10d = 2.67d12
xnu21d = 5.32d12
xnu32d = 7.93d12
xnu20d = xnu21d+xnu10d
xnu31d = xnu32d+xnu21d
!     in K
DT10d = h_Pl*xnu10d/xk_B
DT21d = h_Pl*xnu21d/xk_B
DT32d = h_Pl*xnu32d/xk_B
DT20d = h_Pl*xnu20d/xk_B
DT31d = h_Pl*xnu31d/xk_B

!     ---------------------------------------------------------------
!     Collisional transitions for HD, from Shaefer (1990) [HD,He]
!     with a corrective factor for [HD,H] from Wright & Morton (1979)
!     ---------------------------------------------------------------
c10d = (4.4d-12+3.6d-13*T_K**0.77d0)*y_H*xnH/1.27d0
c21d = (4.1d-12+2.0d-13*T_K**0.92d0)*y_H*xnH/1.27d0
c32d = (2.4d-12+8.7d-14*T_K**1.03d0)*y_H*xnH/1.27d0
c20d = (3.4d-13+1.1d-14*T_K**1.12d0)*y_H*xnH/1.27d0
c31d = (3.2d-13+1.3d-15*T_K**1.47d0)*y_H*xnH/1.27d0

c01d = (g1/g0)*c10d*dexp(-h_Pl*xnu10d/(xk_B*T_K))
c12d = (g2/g1)*c21d*dexp(-h_Pl*xnu21d/(xk_B*T_K))
c23d = (g3/g2)*c32d*dexp(-h_Pl*xnu32d/(xk_B*T_K))
c02d = (g2/g0)*c20d*dexp(-h_Pl*xnu20d/(xk_B*T_K))
c13d = (g3/g1)*c31d*dexp(-h_Pl*xnu31d/(xk_B*T_K))

q10d = Q_bg(DT10d)
q21d = Q_bg(DT21d)
q32d = Q_bg(DT32d)
q20d = Q_bg(DT20d)
q31d = Q_bg(DT31d)
!---------------------------------------------------------------
! Level populations for HD
!---------------------------------------------------------------

w10d=a10d*(1.d0+q10d)+c10d
w21d=a21d*(1.d0+q21d)+c21d
w32d=a32d*(1.d0+q32d)+c32d
w20d=a20d*(1.d0+q20d)+c20d
w31d=a31d*(1.d0+q31d)+c31d
w30d=0.d0

w01d=(g1/g0)*a10d*q10d+c01d
w12d=(g2/g1)*a21d*q21d+c12d
w23d=(g3/g2)*a32d*q32d+c23d
w02d=(g2/g0)*a20d*q20d+c02d
w13d=(g3/g1)*a31d*q31d+c13d
w03d=0.d0

w0d = w01d+w02d+w03d
w1d = w10d+w12d+w13d
w2d = w20d+w21d+w23d

an0d = -(w10d*w21d*w32d+w30d*w2d*w1d+w20d*w31d*w12d-&
w12d*w21d*w30d+w2d*w31d*w10d+w32d*w20d*w1d)
an1d = -w0d*w21d*w32d+w20d*w31d*w02d-w30d*w2d*w01d-&
w02d*w21d*w30d-w2d*w31d*w0d-w32d*w20d*w01d
an2d = -(w0d*w1d*w32d+w10d*w31d*w02d+w30d*w12d*w01d+&
w02d*w1d*w30d+w12d*w31d*w0d-w32d*w10d*w01d)
an3d = -w0d*w1d*w2d+w10d*w21d*w02d+w20d*w12d*w01d+&
w02d*w1d*w20d+w12d*w21d*w0d+w2d*w10d*w01d
antotd = an0d+an1d+an2d+an3d
f0d = an0d/antotd
f1d = an1d/antotd
f2d = an2d/antotd
f3d = an3d/antotd

xNc0d = f0d*xNc_HD
xNc1d = f1d*xNc_HD
xNc2d = f2d*xNc_HD
xNc3d = f3d*xNc_HD

v_D=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p))

tau10d = (a10d/8.d0/pi)*(3.d10/xnu10d)**3&
*(xNc0d*g1/g0-xNc1d)/v_D
tau21d = (a21d/8.d0/pi)*(3.d10/xnu21d)**3&
*(xNc1d*g2/g1-xNc2d)/v_D
tau32d = (a32d/8.d0/pi)*(3.d10/xnu32d)**3&
*(xNc2d*g3/g2-xNc3d)/v_D
tau20d = (a20d/8.d0/pi)*(3.d10/xnu20d)**3&
*(xNc0d*g2/g0-xNc2d)/v_D
tau31d = (a31d/8.d0/pi)*(3.d10/xnu31d)**3&
*(xNc1d*g3/g1-xNc3d)/v_D

esc10d=beta_esc(tau10d,tau_cnt)
esc21d=beta_esc(tau21d,tau_cnt)
esc32d=beta_esc(tau32d,tau_cnt)
esc20d=beta_esc(tau20d,tau_cnt)
esc31d=beta_esc(tau31d,tau_cnt)

if(f1d /= 0.d0) then
   s10d=1.d0/(g1*f0d/(g0*f1d)-1.d0)
   x10d=esc10d*(1.d0-q10d/s10d)
else
   x10d=0.d0
endif
if(f2d /= 0.d0) then
   s21d=1.d0/(g2*f1d/(g1*f2d)-1.d0)
   s20d=1.d0/(g2*f0d/(g0*f2d)-1.d0)
   x21d=esc21d*(1.d0-q21d/s21d)
   x20d=esc20d*(1.d0-q20d/s20d)
else
   x21d=0.d0
   x20d=0.d0
endif
if(f3d /= 0.d0) then
   s32d=1.d0/(g3*f2d/(g2*f3d)-1.d0)
   s31d=1.d0/(g3*f1d/(g1*f3d)-1.d0)
   x32d=esc32d*(1.d0-q32d/s32d)
   x31d=esc31d*(1.d0-q31d/s31d)
else
   x32d=0.d0
   x31d=0.d0
endif

!---------------------------------------------------------------
! Computes cooling rate for HD
!---------------------------------------------------------------
xLd10d = x10d*f1d*a10d*h_Pl*xnu10d/xnH
xLd21d = x21d*f2d*a21d*h_Pl*xnu21d/xnH
xLd32d = x32d*f3d*a32d*h_Pl*xnu32d/xnH
xLd20d = x20d*f2d*a20d*h_Pl*xnu20d/xnH
xLd31d = x31d*f3d*a31d*h_Pl*xnu31d/xnH
xLd_HD = xLd10d+xLd21d+xLd32d+xLd20d+xLd31d
return
END FUNCTION xLd_HD

FUNCTION xLd(pop,esc,N_line)
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd,pop
integer :: N_line
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p),desc(N_p),&
esc_f(N_p),func_f(N_p),A(N_p,N_p),&
esc_min(N_p)

err_max_min=1.d10
do 10 itr=1,1000
   call pop(esc,func,xLd)
   do j=1,N_line
      if(esc(j).eq.0.d0) then
         delta_esc=1.d-10
      else
         delta_esc=1.d-5*esc(j)
      endif
      do jj=1,N_p
         if(jj.eq.j) then
            esc_f(jj)=esc(jj)+delta_esc
         else
            esc_f(jj)=esc(jj)
         endif
      enddo
      
      call pop(esc_f,func_f,xLd)
! ------ set the matrix A
      do i=1,N_line
         A(i,j)=(func_f(i)-func(i))/(esc_f(j)-esc(j))
      enddo
   enddo
! ------- set the vector desc
   do i=1,N_line
      desc(i)=-func(i)
   enddo

! ------- solve linear equations
   call gaussj(A,N_line,N_p,desc,1,1)
   if(itr.eq.1000) go to 20

   fact=1.d0
   if(itr.gt.20) then
      if(err_max.gt.1.d0) fact=1.d-2
   endif

   do i=1,N_line
      if(esc(i)*desc(i).ne.0.d0) then
         fact=dmin1(fact,4.d-1*dabs(esc(i)/desc(i)))
      endif
   enddo

   do i=1,N_line
      esc(i)=esc(i)+fact*desc(i)
   enddo

   err_max=0.d0
   do i=1,N_line
      if(esc(i).ne.0.d0) then
         err=dabs(desc(i)/esc(i))
      else
         err=0.d0
      endif
      err_max=max(err,err_max)
   enddo
   if(err_max.lt.err_max_min) then
      err_max_min=err_max
      do i=1,N_line
         esc_min(i)=esc(i)
      enddo
   endif

   if(itr.le.10) then
      err_th=1.d-12
   elseif(itr.le.40) then
      err_th=1.d-8
   else
      err_th=1.d-4
   endif
   if(err_max.lt.err_th) go to 20

10   continue

20   continue

do i=1,N_line
   esc(i)=esc_min(i)
enddo

call pop(esc,func,xLd)
return
END FUNCTION xLd

SUBROUTINE pop_rot(esc,func,xLd_rot)
!     xnH*rn(molecule)*xLd_rot : cooling rate owing to
!     rotational transitions per unit volume
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd_rot
!     USES beta_esc
integer, PARAMETER :: N_p=100
double precision :: f_LTE(0:N_p),A(N_p),f(0:N_p),aa(N_p),&
DDD(0:N_p),xnn(0:N_p),esc(N_p),func(N_p),Q(N_p),S(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

Z=0.d0
do j=0,J_max
   Z=Z+dble(2*j+1)*dexp(-DT0*dble(j*(j+1))/T_K)
enddo

v_T=dsqrt(8.d0*xk_B*T_K*(1.d0-y_m)/(pi*xm_p))
do j=0,J_max
   f_LTE(j)=dble(2*j+1)*dexp(-DT0*dble(j*(j+1))/T_K)/Z
enddo

do j=1,J_max
   A(j)=2.d0*A0*dble(j**4)/dble(2*j+1)
   DT=2.d0*dble(j)*DT0
   Q(j)=Q_bg(DT)
   aa(j)=A(j)*esc(j)/(xnH*sigma*v_T)/(1.d0-y_m)
enddo

do j=0,J_max
   if(j == 0) then
      AAA=0.d0
      BBB=(dble(2*j+3)/dble(2*j+1))*aa(j+1)*Q(j+1)+1.d0
      CCC=aa(j+1)*(1.d0+Q(j+1))
   elseif(j == J_max) then
      AAA=(dble(2*j+1)/dble(2*j-1))*aa(j)*Q(j)
      BBB=aa(j)*(1.d0+Q(j))+1.d0
      CCC=0.d0
   else
      AAA=(dble(2*j+1)/dble(2*j-1))*aa(j)*Q(j)
      BBB=aa(j)*(1.d0+Q(j))+(dble(2*j+3)/dble(2*j+1))*aa(j+1)*Q(j+1)+1.d0
      CCC=aa(j+1)*(1.d0+Q(j+1))
   endif

   if(j == 0) then
      DDD(j)=CCC/BBB
      xnn(j)=f_LTE(0)/BBB
   else
      DDD(j)=CCC/(BBB-AAA*DDD(j-1))
      xnn(j)=(f_LTE(j)+AAA*xnn(j-1))/(BBB-AAA*DDD(j-1))
   endif
enddo
!     back substitution
f(J_max)=xnn(J_max)
do j=J_max-1,0,-1
   f(j)=DDD(j)*f(j+1)+xnn(j)
enddo

xLd_rot=0.d0
xm_mol=xmu_mol*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/xm_mol)
do j=1,J_max
   DE=DT0*dble(2*j)*xk_B
   xnu=DE/h_Pl
   xNc_l=xNc*f(j-1)
   xNc_u=xNc*f(j)
   g_l=dble(2*j-1)
   g_u=dble(2*j+1)
   if(f(j).eq.0.d0) then
      S(j)=1.d0/(dexp(DT0*2.d0*dble(j)/T_K)-1.d0)
   elseif((g_u*f(j-1)/(g_l*f(j))-1.d0).eq.0.d0) then
      S(j)=1.d50
   else
      S(j)=1.d0/(g_u*f(j-1)/(g_l*f(j))-1.d0)
   endif
   tau=(A(j)/8.d0/pi)*(3.d10/xnu)**3*(xNc_l*g_u/g_l-xNc_u)/v_th/eta_T
   func(j)=esc(j)-beta_esc(tau,tau_cnt)
   xLd_rot=xLd_rot+f(j)*A(j)*DE*beta_esc(tau,tau_cnt)*(1.d0-Q(j)/S(j))/xnH
enddo
END SUBROUTINE pop_rot


SUBROUTINE pop_CII(esc,func,xLd_CII)
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd_CII
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

g_0=2.d0 !Degeneracy
g_1=4.d0 !Degeneracy

!esc_10=esc(1)
esc_10=1.d0
DT_10=9.2d1 !T_line
DE_10=DT_10*1.38d-16 !T_line*k_B
Q_10=0.d0 !Q_bg(DT_10)

A_10=2.4d-6 !Einstein A coefficient
gamma_e=2.8d-7/dsqrt(T_K*1.d-2) !Sigma v
gamma_H=8.0d-10*(T_K*1.d-2)**7.d-2 !Sigma v
gamma_H2=5.d-1*gamma_H !Sigma v
C_10=xnH*(y_e*gamma_e+y_a*gamma_H+y_m*gamma_H2) !Einstein C coefficient
C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K) !Einstein C coefficient

R_01=(g_1/g_0)*A_10*esc_10*Q_10+C_01
R_10=esc_10*A_10*(1.d0+Q_10)+C_10
f_1=R_01/(R_10+R_01)
f_0=R_10/(R_10+R_01)

DE_10=DT_10*xk_B
xnu=DE_10/h_Pl
xm_C=12.d0*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/xm_C)
xNc_1=xNc*f_1
xNc_0=xNc*f_0
!tau_10=(A_10/8.d0/pi)*(3.d10/xnu)**3*(xNc_0*g_1/g_0-xNc_1)
!&     /v_th
tau_10=0.d0

if(g_1*f_0/(g_0*f_1)-1.d0.eq.0.d0) then
   S_10=1.d50
else
   S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
endif

func(1)=esc_10-beta_esc(tau_10,tau_cnt)

xLd_CII=DE_10*A_10*f_1*esc_10*(1.d0-Q_10/S_10)/xnH
return

END SUBROUTINE pop_CII

SUBROUTINE pop_CI(esc,func,xLd_CI)
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd_CI
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

g_0=1.d0
g_1=3.d0
g_2=5.d0

!esc_10=esc(1)
!esc_20=esc(2)
!esc_21=esc(3)
esc_10=1.d0
esc_20=1.d0
esc_21=1.d0

DT_10=2.4d1
DT_20=6.3d1
DT_21=3.9d1

DE_10=DT_10*1.38d-16
DE_20=DT_20*1.38d-16
DE_21=DT_21*1.38d-16

!Q_10=Q_bg(DT_10)
!Q_20=Q_bg(DT_20)
!Q_21=Q_bg(DT_21)

Q_10=0.d0
Q_20=0.d0
Q_21=0.d0

A_10=7.9d-8
A_20=2.0d-14
A_21=2.7d-7

gamma10_e=3.0d-9
gamma20_e=5.0d-9
gamma21_e=1.5d-8

gamma10_H=1.6d-10*(T_K*1.d-2)**0.14d0
gamma20_H=9.2d-11*(T_K*1.d-2)**0.26d0
gamma21_H=2.9d-10*(T_K*1.d-2)**0.26d0

gamma10_H2=5.d-2*gamma10_H
gamma20_H2=5.d-1*gamma20_H
gamma21_H2=5.d-1*gamma21_H

C_10=xnH*(y_e*gamma10_e+y_a*gamma10_H+y_m*gamma10_H2)
C_20=xnH*(y_e*gamma20_e+y_a*gamma20_H+y_m*gamma20_H2)
C_21=xnH*(y_e*gamma21_e+y_a*gamma21_H+y_m*gamma21_H2)
C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)
C_02=(g_2/g_0)*C_20*dexp(-DT_20/T_K)
C_12=(g_2/g_1)*C_21*dexp(-DT_21/T_K)

R_10=esc_10*A_10*(1.d0+Q_10)+C_10
R_20=esc_20*A_20*(1.d0+Q_20)+C_20
R_21=esc_21*A_21*(1.d0+Q_21)+C_21
R_01=(g_1/g_0)*esc_10*A_10*Q_10+C_01
R_02=(g_2/g_0)*esc_20*A_20*Q_20+C_02
R_12=(g_2/g_1)*esc_21*A_21*Q_21+C_12

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )&
/( (R_01+R_02+R_20)*(R_10+R_12+R_21)&
-(R_01-R_21)*(R_10-R_20) )
f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

xnu_10=DE_10/h_Pl
xnu_20=DE_20/h_Pl
xnu_21=DE_21/h_Pl

rm_C=12.d0*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/rm_C)

xNc_0=xNc*f_0
xNc_1=xNc*f_1
xNc_2=xNc*f_2

!tau_10=(A_10/8.d0/pi)*(3.d10/xnu_10)**3
!&     *(xNc_0*g_1/g_0-xNc_1)/v_th
!tau_20=(A_20/8.d0/pi)*(3.d10/xnu_20)**3
!&     *(xNc_0*g_2/g_0-xNc_2)/v_th
!tau_21=(A_21/8.d0/pi)*(3.d10/xnu_21)**3
!&     *(xNc_1*g_2/g_1-xNc_2)/v_th

tau_10=0.d0
tau_20=0.d0
tau_21=0.d0

if(g_1*f_0/(g_0*f_1)-1.d0.eq.0.d0) then
   S_10=1.d50
else
   S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
endif
if(g_2*f_0/(g_0*f_2)-1.d0.eq.0.d0) then
   S_20=1.d50
else
   S_20=1.d0/(g_2*f_0/(g_0*f_2)-1.d0)
endif
if(g_2*f_1/(g_1*f_2)-1.d0.eq.0.d0) then
   S_21=1.d50
else
   S_21=1.d0/(g_2*f_1/(g_1*f_2)-1.d0)
endif

func(1)=esc_10-beta_esc(tau_10,tau_cnt)
func(2)=esc_20-beta_esc(tau_20,tau_cnt)
func(3)=esc_21-beta_esc(tau_21,tau_cnt)

x_10=esc_10*(1.d0-Q_10/S_10)
x_20=esc_20*(1.d0-Q_20/S_20)
x_21=esc_21*(1.d0-Q_21/S_21)

xLd_CI=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20&
+DE_21*A_21*f_2*x_21)/xnH
return

END SUBROUTINE pop_CI

SUBROUTINE pop_OI(esc,func,xLd_OI)
!     xnH*rn(OI)*xLd_OI : cooling rate owing to OI per unit volume
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd_OI
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

g_0=5.d0
g_1=3.d0
g_2=1.d0

esc_10=1.0d0!esc(1)
esc_20=1.0d0!esc(2)
esc_21=1.0d0!esc(3)

DT_10=2.3d2
DT_20=3.28d2
DT_21=9.8d1

DE_10=DT_10*1.38d-16 !level energy
DE_20=DT_20*1.38d-16
DE_21=DT_21*1.38d-16

!Q_10=Q_bg(DT_10)
!Q_20=Q_bg(DT_20)
!Q_21=Q_bg(DT_21)
Q_10=0.d0
Q_20=0.d0
Q_21=0.d0

A_10=9.0d-5 !Einstein A
A_20=1.0d-10
A_21=1.7d-5

gamma10_e=1.4d-8 !Einstein C by e-
gamma20_e=1.4d-8
gamma21_e=5.0d-9

gamma10_H=9.2d-11*(T_K*1.d-2)**0.67d0 !Einstein C by H
gamma20_H=4.3d-11*(T_K*1.d-2)**0.80d0
gamma21_H=1.1d-10*(T_K*1.d-2)**0.44d0

gamma10_H2=5.d-1*gamma10_H !Einstein C by H2
gamma20_H2=5.d-1*gamma20_H
gamma21_H2=5.d-2*gamma21_H

C_10=xnH*(y_e*gamma10_e+y_a*gamma10_H+y_m*gamma10_H2)
C_20=xnH*(y_e*gamma20_e+y_a*gamma20_H+y_m*gamma20_H2)
C_21=xnH*(y_e*gamma21_e+y_a*gamma21_H+y_m*gamma21_H2)
C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)
C_02=(g_2/g_0)*C_20*dexp(-DT_20/T_K)
C_12=(g_2/g_1)*C_21*dexp(-DT_21/T_K)

!esc_10=esc(1)
!esc_20=esc(2)
!esc_21=esc(3)
esc_10=1.0d0!esc(1)
esc_20=1.0d0!esc(2)
esc_21=1.0d0!esc(3)

R_10=esc_10*A_10*(1.d0+Q_10)+C_10 !Reaction rate
R_20=esc_20*A_20*(1.d0+Q_20)+C_20
R_21=esc_21*A_21*(1.d0+Q_21)+C_21
R_01=(g_1/g_0)*esc_10*A_10*Q_10+C_01
R_02=(g_2/g_0)*esc_20*A_20*Q_20+C_02
R_12=(g_2/g_1)*esc_21*A_21*Q_21+C_12

!population
f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )&
/( (R_01+R_02+R_20)*(R_10+R_12+R_21)&
-(R_01-R_21)*(R_10-R_20) )
f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

xnu_10=DE_10/h_Pl !frequency
xnu_20=DE_20/h_Pl
xnu_21=DE_21/h_Pl

rm_O=16.d0*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/rm_O) !thermal velocity of O

xNc_0=xNc*f_0
xNc_1=xNc*f_1
xNc_2=xNc*f_2

!tau_10=(A_10/8.d0/pi)*(3.d10/xnu_10)**3
!&     *(xNc_0*g_1/g_0-xNc_1)/v_th
!tau_20=(A_20/8.d0/pi)*(3.d10/xnu_20)**3
!&     *(xNc_0*g_2/g_0-xNc_2)/v_th
!tau_21=(A_21/8.d0/pi)*(3.d10/xnu_21)**3
!&     *(xNc_1*g_2/g_1-xNc_2)/v_th
tau_10=0.d0
tau_20=0.d0
tau_21=0.d0

if(g_1*f_0/(g_0*f_1)-1.d0.eq.0.d0) then
   S_10=1.d50
else
   S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
endif
if(g_2*f_0/(g_0*f_2)-1.d0.eq.0.d0) then
   S_20=1.d50
else
   S_20=1.d0/(g_2*f_0/(g_0*f_2)-1.d0)
endif
if(g_2*f_1/(g_1*f_2)-1.d0.eq.0.d0) then
   S_21=1.d50
else
   S_21=1.d0/(g_2*f_1/(g_1*f_2)-1.d0)
endif

func(1)=esc_10-beta_esc(tau_10,tau_cnt)
func(2)=esc_20-beta_esc(tau_20,tau_cnt)
func(3)=esc_21-beta_esc(tau_21,tau_cnt)

x_10=esc_10*(1.d0-Q_10/S_10)
x_20=esc_20*(1.d0-Q_20/S_20)
x_21=esc_21*(1.d0-Q_21/S_21)

xLd_OI=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20&
+DE_21*A_21*f_2*x_21)/xnH

END SUBROUTINE pop_OI


SUBROUTINE pop_CIImeta(esc,func,xLd_CIImeta)
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd_CIImeta
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

g_0=6.d0
g_1=12.d0

esc_10=1.d0!esc(1)
DT_10=6.2d4
DE_10=DT_10*1.38d-16
Q_10=0.d0!Q_bg(DT_10)

A_10=3.6d0
gamma_e=2.3d-8/dsqrt(T_K*1.d-4)
C_10=xnH*y_e*gamma_e
C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)

R_01=(g_1/g_0)*A_10*esc_10*Q_10+C_01
R_10=esc_10*A_10*(1.d0+Q_10)+C_10
f_1=R_01/(R_10+R_01)
f_0=R_10/(R_10+R_01)

DE_10=DT_10*xk_B
xnu=DE_10/h_Pl
xm_C=12.d0*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/xm_C)
xNc_1=xNc*f_1
xNc_0=xNc*f_0
!tau_10=(A_10/8.d0/pi)*(3.d10/xnu)**3*(xNc_0*g_1/g_0-xNc_1)
!&     /v_th
tau_10=0.d0

if(g_1*f_0/(g_0*f_1)-1.d0 .eq.0.d0) then
   S_10=1.d50
else
   S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
endif

func(1)=esc_10-beta_esc(tau_10,tau_cnt)

xLd_CIImeta=DE_10*A_10*f_1*esc_10*(1.d0-Q_10/S_10)/xnH

!return

END SUBROUTINE pop_CIImeta

SUBROUTINE pop_CImeta(esc,func,xLd_CImeta)
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: Ld_CImeta
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

g_0=9.d0
g_1=5.d0
g_2=1.d0

esc_10=1.d0!esc(1)
esc_20=1.d0!esc(2)
esc_21=1.d0!esc(3)

DT_10=1.5d4
DT_20=3.2d4
DT_21=1.7d4

DE_10=DT_10*1.38d-16
DE_20=DT_20*1.38d-16
DE_21=DT_21*1.38d-16

Q_10=0.d0!Q_bg(DT_10)
Q_20=0.d0!Q_bg(DT_20)
Q_21=0.d0!Q_bg(DT_21)

A_10=3.4d-4
A_20=2.6d-3
A_21=0.5d0

T4=1.d-4*T_K
if(T4 < 1.d0) then
   gamma10_e=2.7d-8*T4**0.57
   gamma20_e=1.3d-8*T4**0.57
   gamma21_e=1.2d-8*T4**0.57
else
   gamma10_e=2.7d-8*T4**(-0.13)
   gamma20_e=1.3d-8
   gamma21_e=1.2d-8
endif

C_10=xnH*y_e*gamma10_e
C_20=xnH*y_e*gamma20_e
C_21=xnH*y_e*gamma21_e
C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)
C_02=(g_2/g_0)*C_20*dexp(-DT_20/T_K)
C_12=(g_2/g_1)*C_21*dexp(-DT_21/T_K)

R_10=esc_10*A_10*(1.d0+Q_10)+C_10
R_20=esc_20*A_20*(1.d0+Q_20)+C_20
R_21=esc_21*A_21*(1.d0+Q_21)+C_21
R_01=(g_1/g_0)*esc_10*A_10*Q_10+C_01
R_02=(g_2/g_0)*esc_20*A_20*Q_20+C_02
R_12=(g_2/g_1)*esc_21*A_21*Q_21+C_12

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )&
/( (R_01+R_02+R_20)*(R_10+R_12+R_21)&
-(R_01-R_21)*(R_10-R_20) )
f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

xnu_10=DE_10/h_Pl
xnu_20=DE_20/h_Pl
xnu_21=DE_21/h_Pl

rm_C=12.d0*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/rm_C)

xNc_0=xNc*f_0
xNc_1=xNc*f_1
xNc_2=xNc*f_2

!tau_10=(A_10/8.d0/pi)*(3.d10/xnu_10)**3
!&     *(xNc_0*g_1/g_0-xNc_1)/v_th
!tau_20=(A_20/8.d0/pi)*(3.d10/xnu_20)**3
!&     *(xNc_0*g_2/g_0-xNc_2)/v_th
!tau_21=(A_21/8.d0/pi)*(3.d10/xnu_21)**3
!&     *(xNc_1*g_2/g_1-xNc_2)/v_th
tau_10=0.d0
tau_20=0.d0
tau_21=0.d0


if(g_1*f_0/(g_0*f_1)-1.d0.eq.0.d0) then
   S_10=1.d50
else
   S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
endif
if(g_2*f_0/(g_0*f_2)-1.d0.eq.0.d0) then
   S_20=1.d50
else
   S_20=1.d0/(g_2*f_0/(g_0*f_2)-1.d0)
endif
if(g_2*f_1/(g_1*f_2)-1.d0.eq.0.d0) then
   S_21=1.d50
else
   S_21=1.d0/(g_2*f_1/(g_1*f_2)-1.d0)
endif

func(1)=esc_10-beta_esc(tau_10,tau_cnt)
func(2)=esc_20-beta_esc(tau_20,tau_cnt)
func(3)=esc_21-beta_esc(tau_21,tau_cnt)

x_10=esc_10*(1.d0-Q_10/S_10)
x_20=esc_20*(1.d0-Q_20/S_20)
x_21=esc_21*(1.d0-Q_21/S_21)

xLd_CImeta=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20&
+DE_21*A_21*f_2*x_21)/xnH

!return
END SUBROUTINE pop_CImeta

SUBROUTINE pop_OImeta(esc,func,xLd_OImeta)
!     xnH*rn(OI)*xLd_OI : cooling rate owing to OI per unit volume
use COMVAR
IMPLICIT REAL*8(a-h,o-z)
double precision :: xLd_OImeta
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p)
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

g_0=9.d0
g_1=5.d0
g_2=1.d0

esc_10=1.d0!esc(1)
esc_20=1.d0!esc(2)
esc_21=1.d0!esc(3)

DT_10=2.3d4
DT_20=4.9d4
DT_21=2.6d4

DE_10=DT_10*1.38d-16
DE_20=DT_20*1.38d-16
DE_21=DT_21*1.38d-16

Q_10=0.d0!Q_bg(DT_10)
Q_20=0.d0!Q_bg(DT_20)
Q_21=0.d0!Q_bg(DT_21)

A_10=6.3d-3
A_20=6.7d-2
A_21=1.3d0

T4=1.d-4*T_K
if(T4 < 1.d0) then
   gamma10_e=5.1d-9*T4**0.57d0
   gamma20_e=2.5d-9*T4**0.57d0
   gamma21_e=5.2d-9*T4**0.57d0
else
   gamma10_e=5.1d-9*T4**0.17d0
   gamma20_e=2.5d-9*T4**0.13d0
   gamma21_e=5.2d-9*T4**0.15d0
endif

C_10=xnH*y_e*gamma10_e
C_20=xnH*y_e*gamma20_e
C_21=xnH*y_e*gamma21_e
C_01=(g_1/g_0)*C_10*dexp(-DT_10/T_K)
C_02=(g_2/g_0)*C_20*dexp(-DT_20/T_K)
C_12=(g_2/g_1)*C_21*dexp(-DT_21/T_K)

esc_10=1.d0!esc(1)
esc_20=1.d0!esc(2)
esc_21=1.d0!esc(3)

R_10=esc_10*A_10*(1.d0+Q_10)+C_10
R_20=esc_20*A_20*(1.d0+Q_20)+C_20
R_21=esc_21*A_21*(1.d0+Q_21)+C_21
R_01=(g_1/g_0)*esc_10*A_10*Q_10+C_01
R_02=(g_2/g_0)*esc_20*A_20*Q_20+C_02
R_12=(g_2/g_1)*esc_21*A_21*Q_21+C_12

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )&
/( (R_01+R_02+R_20)*(R_10+R_12+R_21)&
-(R_01-R_21)*(R_10-R_20) )
f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

if(f_1 == 0.d0) then
   xLd_OImeta=0.d0
   return
endif


xnu_10=DE_10/h_Pl
xnu_20=DE_20/h_Pl
xnu_21=DE_21/h_Pl

rm_O=16.d0*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/rm_O)

xNc_0=xNc*f_0
xNc_1=xNc*f_1
xNc_2=xNc*f_2

!tau_10=(A_10/8.d0/pi)*(3.d10/xnu_10)**3
!&     *(xNc_0*g_1/g_0-xNc_1)/v_th
!tau_20=(A_20/8.d0/pi)*(3.d10/xnu_20)**3
!&     *(xNc_0*g_2/g_0-xNc_2)/v_th
!tau_21=(A_21/8.d0/pi)*(3.d10/xnu_21)**3
!&     *(xNc_1*g_2/g_1-xNc_2)/v_th
tau_10=0.d0
tau_20=0.d0
tau_21=0.d0

if(g_1*f_0/(g_0*f_1)-1.d0.eq.0.d0) then
   S_10=1.d50
else
   S_10=1.d0/(g_1*f_0/(g_0*f_1)-1.d0)
endif
if(g_2*f_0/(g_0*f_2)-1.d0.eq.0.d0) then
   S_20=1.d50
else
   S_20=1.d0/(g_2*f_0/(g_0*f_2)-1.d0)
endif
if(g_2*f_1/(g_1*f_2)-1.d0.eq.0.d0) then
   S_21=1.d50
else
   S_21=1.d0/(g_2*f_1/(g_1*f_2)-1.d0)
endif

func(1)=esc_10-beta_esc(tau_10,tau_cnt)
func(2)=esc_20-beta_esc(tau_20,tau_cnt)
func(3)=esc_21-beta_esc(tau_21,tau_cnt)

x_10=esc_10*(1.d0-Q_10/S_10)
x_20=esc_20*(1.d0-Q_20/S_20)
x_21=esc_21*(1.d0-Q_21/S_21)

xLd_OImeta=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20&
+DE_21*A_21*f_2*x_21)/xnH
!return
END SUBROUTINE pop_OImeta

SUBROUTINE thick_lev(J_thick)
!double precision :: thick_lev
USE COMVAR
IMPLICIT REAL*8(a-h,o-z)
integer J_thick,j
integer, PARAMETER :: N_p=100
double precision :: f_LTE(0:N_p),A(N_p),f(0:N_p),aa(N_p),P_J(N_p),esc(N_p)
!COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
!COMMON /mol_lines/ A0,DT0,sigma,eta_T,xmu_mol,J_max
!double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
!xm_p=1.67d-24

J_thick=0

do j=1,J_max
   esc(j)=1.d0
enddo

Z=0.d0
do j=0,J_max
   Z=Z+dble(2*j+1)*dexp(-DT0*dble(j*(j+1))/T_K)
enddo

v_T=dsqrt(8.d0*xk_B*T_K*(1.d0-y_m)/(pi*xm_p))
do j=0,J_max
   f_LTE(j)=dble(2*j+1)*dexp(-DT0*dble(j*(j+1))/T_K)/Z
enddo

do j=1,J_max
   A(j)=2.d0*A0*dble(j**4)/dble(2*j+1)
   aa(j)=A(j)*esc(j)/(xnH*sigma*v_T)/(1.d0-y_m)
enddo

do J=1,J_max
   P_J(J)=aa(J)/(1.d0+aa(J))
   f(J)=f_LTE(J)*P_J(J)/aa(J)
   do jj=J+1,J_max
      P_J(jj)=P_J(jj-1)*(aa(jj)/(1.d0+aa(jj)))
      f(J)=f(J)+f_LTE(jj)*P_J(jj)/aa(J)
   enddo
enddo
f(0)=f_LTE(0)+aa(1)*f(1)

xm_mol=xmu_mol*xm_p
v_th=dsqrt(2.d0*xk_B*T_K/xm_mol)
do j=1,J_max
   DE=DT0*2.d0*dble(j)*xk_B
   xnu=DE/h_Pl
   xNc_l=xNc*f(j-1)
   xNc_u=xNc*f(j)
   g_l=dble(2*j-1)
   g_u=dble(2*j+1)
   tau_line=(A(j)/8.d0/pi)*(3.d10/xnu)**3&
   *(xNc_l*g_u/g_l-xNc_u)/v_th/eta_T
   if(tau_line.gt.1.d-2) then
      J_thick=j
   endif
enddo

!return
END SUBROUTINE thick_lev



FUNCTION beta_esc(tau_L,tau_C)
IMPLICIT REAL*8(a-h,o-z)
double precision :: beta_esc,tau_L,tau_C

if(tau_L.lt.0.d0) then
   beta_esc=1.d0
elseif(tau_L.lt.1.d-5) then
   beta_esc=dexp(-tau_C)
else
   beta_esc=dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
endif

return
END FUNCTION beta_esc


SUBROUTINE OHcool(xnH,T_K,y_H2,xNc_OH,tau_cnt,xLd_OH)
IMPLICIT REAL*8(a-h,o-z)
double precision :: xnH,T_K,y_H2,xNc_OH,tau_cnt,xLd_OH
!     OH cooling function by Leiden data
double precision :: xk_B=1.38d-16, xm_p=1.67d-24
double precision :: xlTa(1:6)
data xlTa/1.477d0, 1.699d0, 1.903d0,&
2.000d0, 2.477d0, 2.778d0/
double precision :: xlNa(1:9)
data xlNa/10.0d0, 11.d0, 12.0d0, 13.d0, 14.d0,&
15.d0, 16.d0, 17.d0, 18.d0/
double precision :: aL0a(1:6)
data aL0a/0.2531d2, 0.2453d2, 0.2402d2, 0.2382d2, 0.2316d2,&
0.2290d2/
double precision :: aLLTEa(1:6,1:9)
data aLLTEa/0.1620d2, 0.1544d2, 0.1490d2, 0.1467d2, 0.1381d2,&
0.1358d2,&
0.1620d2, 0.1545d2, 0.1490d2, 0.1467d2, 0.1381d2, 0.1358d2,&
0.1622d2, 0.1546d2, 0.1491d2, 0.1467d2, 0.1381d2, 0.1359d2,&
0.1636d2, 0.1555d2, 0.1496d2, 0.1471d2, 0.1383d2, 0.1359d2,&
0.1702d2, 0.1600d2, 0.1525d2, 0.1495d2, 0.1394d2, 0.1367d2,&
0.1774d2, 0.1660d2, 0.1577d2, 0.1545d2, 0.1449d2, 0.1415d2,&
0.1858d2, 0.1738d2, 0.1650d2, 0.1616d2, 0.1510d2, 0.1477d2,&
0.1942d2, 0.1828d2, 0.1740d2, 0.1705d2, 0.1584d2, 0.1540d2,&
0.2038d2, 0.1921d2, 0.1836d2, 0.1802d2, 0.1683d2, 0.1635d2/
double precision :: xlnha(1:6,1:9)
data xlnha/0.905d1,  0.898d1,  0.890d1,  0.888d1,  0.877d1,&
0.872d1,&
0.905d1,  0.898d1,  0.890d1,  0.887d1,  0.877d1,  0.871d1,&
0.903d1,  0.897d1,  0.889d1,  0.887d1,  0.876d1,  0.871d1,&
0.888d1,  0.885d1,  0.881d1,  0.880d1,  0.874d1,  0.869d1,&
0.820d1,  0.827d1,  0.836d1,  0.840d1,  0.850d1,  0.850d1,&
0.725d1,  0.745d1,  0.771d1,  0.779d1,  0.799d1,  0.801d1,&
0.626d1,  0.646d1,  0.675d1,  0.685d1,  0.722d1,  0.729d1,&
0.526d1,  0.546d1,  0.575d1,  0.585d1,  0.623d1,  0.633d1,&
0.426d1,  0.446d1,  0.475d1,  0.485d1,  0.523d1,  0.533d1/
double precision :: alphaa(1:6,1:9)
data alphaa/0.353d0,  0.543d0,  0.536d0,  0.530d0,  0.466d0,&
0.434d0,&
0.354d0,  0.543d0,  0.536d0,  0.530d0,  0.466d0,  0.434d0,&
0.364d0,  0.544d0,  0.536d0,  0.529d0,  0.466d0,  0.434d0,&
0.399d0,  0.556d0,  0.534d0,  0.525d0,  0.463d0,  0.432d0,&
0.655d0,  0.586d0,  0.532d0,  0.521d0,  0.447d0,  0.424d0,&
0.641d0,  0.538d0,  0.538d0,  0.538d0,  0.448d0,  0.443d0,&
0.649d0,  0.572d0,  0.552d0,  0.516d0,  0.381d0,  0.370d0,&
0.566d0,  0.618d0,  0.580d0,  0.538d0,  0.381d0,  0.363d0,&
0.701d0,  0.645d0,  0.598d0,  0.550d0,  0.383d0,  0.369d0/

xn_c=(1.d0-y_H2)*xnH

xlT=dlog10(T_K)
cs=1.d-5*dsqrt(2.d0*xk_B*T_K/(17.d0*xm_p))
xNc_OH=xNc_OH+1.d-10
xlN=dlog10(xNc_OH/cs)
call linear(xlTa,aL0a,6,xlT,aL0)
call bilinear(xlTa,xlNa,aLLTEa,6,9,xlT,xlN,aLLTE)
call bilinear(xlTa,xlNa,xlnha,6,9,xlT,xlN,xlnh)
call bilinear(xlTa,xlNa,alphaa,6,9,xlT,xlN,alpha)

xL0inv=10.d0**aL0
xLLTEinv=10.d0**aLLTE
xn_h=10.d0**xlnh
xLinv=xL0inv&
+xn_c*xLLTEinv&
+xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
xL=1.d0/xLinv
xLd_OH=(1.d0-y_H2)*xL*dexp(-tau_cnt)
END


SUBROUTINE react_coef(xnH,T_K,T_gr_K,Z_metal,xk)
IMPLICIT REAL*8(a-h,o-z)
!USE COMVAR
integer, PARAMETER :: N_react=675
double precision :: xnH,T_K,T_gr_K,Z_metal
!common /UVCR/ zeta, G_0
!common /coldens/A_v,xNc_H,xNc_H2,xNc_HD
!common /xcrit/xH,xH2,xHe
double precision :: xk(N_react)

T_eV=8.61735d-5*T_K
xlnT_eV=dlog(T_eV)
T300=T_K/300.d0
xlnT=dlog(T_K)
xlgT=dlog10(T_K)
xlgT2=xlgT**2
xlgT3=xlgT**3
xlgT4=xlgT**4
xlgT5=xlgT**5
xlgT6=xlgT**6
xlgT7=xlgT**7

xlT4=dlog10(T_K/1.d4)
xncr_H=10.d0**(3.d0-0.416d0*xlT4-0.327d0*xlT4**2)
xncr_H2=10.d0**(4.845d0-1.3d0*xlT4+1.62d0*xlT4**2)
xncr_He=10.d0**(5.0792d0*(1.d0-1.23d-5*(T_K-2000.d0)))

xncr=(xH/xncr_H+xH2/xncr_H2+xHe/xncr_He)**(-1.d0)
xcr=xnH/xncr


!  1)   H     +   e     ->   H+    + 2 e
!  none(UMIST)
!     GA08 12
xk(1)=dexp(-32.71396786d0+(13.536556d0&
+(-5.73932875d0+(1.56315498d0+(-0.2877056d0&
+(3.48255977d-2+(-2.63197617d-3&
+(1.11954395d-4-2.03914985d-6*xlnT_eV)*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
!  2)   H+    +   e     ->   H     +   ph.
!  4009(UMIST)
!  GA08 13
!   case A
!      xk(2)=1.269d-13*(315614.d0/T_K)**1.503d0
!     &     *(1.d0+(604625.d0/T_K)**0.470d0)**(-1.923d0)
!   case B
xk(2)=2.753d-14*(315614.d0/T_K)**1.500d0&
*(1.d0+(115188.d0/T_K)**0.407d0)**(-2.242d0)
!  3)   He    +   e     ->   He+   + 2 e
!  UMIST none
!  GA08 17
xk(3)=dexp(-44.09864886d0+(23.91596563d0&
+(-10.7532302d0+(3.05803875d0+(-0.56851189d0&
+(6.79539123d-2+(-5.00905610d-3+(2.06723616d-4&
-3.64916141d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
!  4)   He+   +   e     ->   He    +   ph.
!  4010
!  GA08 19
!     case A
xkrrA=1.d-11*T_K**(-0.5d0)*(12.72d0-1.615d0*xlgT&
-0.3162d0*xlgT2+0.0493d0*xlgT3)
!     case B
xkrrB=1.d-11*T_K**(-0.5d0)*(11.19d0-1.676d0*xlgT&
-0.2852d0*xlgT2+0.04433d0*xlgT3)
!     dielectric
xkdi=1.9d-3*T_K**(-1.5d0)*dexp(-473421.d0/T_K)&
*(1.d0+0.3d0*dexp(-94684.d0/T_K))
!     optically thick
xkthick=0.68d0*xkrrA+0.32d0*xkrrB+xkdi
xk(4)=xkthick
!  5)   He+   +   e     ->   He++  + 2 e
!  none
!  GA08 18
xk(5)=dexp(-68.71040990d0+(43.93347633d0+(-18.4806699d0&
+(4.70162649d0+(-0.76924663d0+(8.113042d-2&
+(-5.32402063d-3+(1.97570531d-4-3.16558106d-6*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)
!  6)   He++  +   e     ->   He+   +   ph.
!  none
!  GA08 20
!     case A
xk6A=2.538d-13*(1262456d0/T_K)**1.503d0&
*(1.d0+(2418500.d0/T_K)**0.470d0)**(-1.923d0)
!     case B
xk6B=5.506d-14*(1262456d0/T_K)**1.500d0&
*(1.d0+(460752.d0/T_K)**0.407d0)**(-2.242d0)
xk(6)=xk6B
!  7)   H     +   e     ->   H-    +   ph.
!  4006,4007,4008
!  GA08 1
if(T_K < 6000.d0) then
   xk(7)=10.d0**(-17.845d0+0.762d0*xlgT&
+0.1523d0*xlgT2-0.03274d0*xlgT3)
else
   xk(7)=10.d0**(-16.4199d0+0.1998d0*xlgT2&
-5.447d-3*xlgT4+4.0415d-5*xlgT6)
endif
!  8)   H-    +   H     ->   H2    +   e
!  4031
!  GA08 2
xk(8)=1.3d-9
!  9)   H     +   H+    ->   H2+   +   ph.
!  4078, 4079
!  GA08 3
xk(9)=1.d1**(-19.38d0-1.523d0*xlgT&
+1.118d0*xlgT2-0.1269d0*xlgT3)
!  10)   H2+   +   H     ->   H2    +   H+
!  2937
!  GA08 4
xk(10)=6.4d-10
!  11)   H2    +   H+    ->   H2+   +   H
!  none
!  GA08 7
if(T_K < 100.d0) then
   xk(11)=0.d0
else
   xk(11)=(-3.3232183d-7+3.3735382d-7*xlnT&
-1.4491368d-7*(xlnT**2)&
+3.4172805d-8*(xlnT**3)-4.7813720d-9*(xlnT**4)&
+3.9731542d-10*(xlnT**5)&
-1.8171411d-11*(xlnT**6)+3.5311932d-13*(xlnT**7))&
*dexp(-21237.15d0/T_K)
endif
!  12)   H2    +   e     -> 2 H     +   e
!  4592
!  GA08 8
!     v=0
xk12v0=4.49d-9*(T_K**0.11d0)*dexp(-101858.d0/T_K)
!     LTE
xk12LTE=1.91d-9*(T_K**0.136d0)*dexp(-53407.1d0/T_K)
xlgk12=(xcr/(1.d0+xcr))*dlog10(xk12LTE)&
+(1.d0/(1.d0+xcr))*dlog10(xk12v0)
xk(12)=10.d0**xlgk12
!  13)   H2    +   H     -> 3 H
!  4584
!     GA08 9
!     v=0
xk13v0=6.67d-12*(T_K**0.5d0)*dexp(-(1.d0+63593.d0/T_K))
!     LTE
xk13LTE=3.52d-9*dexp(-43900.d0/T_K)
xlgk13=(xcr/(1.d0+xcr))*dlog10(xk13LTE)&
+(1.d0/(1.d0+xcr))*dlog10(xk13v0)
xk(13)=10.d0**xlgk13
!  14)   H-    +   e     ->   H     + 2 e
!   none
!     GA08 (14)
xk(14)=dexp(-18.01849334d0+(2.3608522d0&
+(-0.28274430d0+(1.62331664d-2+(-3.36501203d-2&
+(1.17832978d-2+(-1.65619470d-3+(1.06827520d-4&
-2.63128581d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
!  15)   H-    +   H+    -> 2 H
!   3489
!   GA08 (5)
xk(15)=2.4d-6*(T_K**(-0.5d0))*(1.d0+T_K/20000.d0)
!  16)   H-    +   H+    ->   H2+   +   e
!   none
!   GA08 (16)
if(T_K.le.8.d3) then
   xk(16)=6.9d-9*(T_K**(-0.35d0))
else
   xk(16)=9.6d-7*(T_K**(-0.90d0))
endif
!  17)   H2+   +   e     -> 2 H
!   3520
!   GA08 (6)
if(T_K < 617.d0) then
   xk(17)=1.d-8
else
   xk(17)=1.32d-6*T_K**(-0.76d0)
endif
!  18)   H2+   +   H-    ->   H2    +   H
!   3490
!   GA08 (21)
xk(18)=1.4d-7*T300**(-0.5d0)
!  19) 3 H               ->   H2    +   H
!   none
!   Glover 08 (GA08;30)
xk(19)=7.7d-31*T_K**(-0.464d0)
!  20) 2 H     +   H2    -> 2 H2
!  none
!   GA08 (31)
xk(20)=xk(19)/8.d0
!  21) 2 H2              -> 2 H     +   H2
!  4593
!   GA08 (10)
!     v=0
xk21v0=(5.996d-30*T_K**4.1881d0/(1.d0+6.761d-6*T_K)**5.6881d0)&
*dexp(-54657.4d0/T_K)
!     LTE
xk21LTE=1.3d-9*dexp(-53300.d0/T_K)
xlgk21=(xcr/(1.d0+xcr))*dlog10(xk21LTE)&
+(1.d0/(1.d0+xcr))*dlog10(xk21v0)
xk(21)=10.d0**xlgk21
!  22)
!   none
xk(22)=0.d0
!  23) 2 H     +   grain ->   H2
!   none
!     Tielens & Hollenbach (1985)
if(T_gr_K.eq.0.d0) then
   f_a=1.d0
else
   f_a=1.d0/(1.d0+dexp(7.5d2*(1.d0/7.5d1-1.d0/T_gr_K)))
endif
xk(23)=6.d-17*dsqrt(T_K/3.d2)*Z_metal*f_a&
/(1.d0+4.d-2*dsqrt(T_K+T_gr_K)&
+2.d-3*T_K+8.d-6*T_K**2)
!  24)   He+   +   H2    ->   H+    +   H     +   He
!   642
!   GA08 (24)
xk(24)=3.70d-14*dexp(-35.d0/T_K)
!  25)   H2+   +   He    ->   HeH+  +   H
!   641
xk(25)=1.30d-10
!  26)   H2+   +   H2    ->   H3+   +   H
!   640
xk(26)=2.08d-9
!  27)   H3+   +   H-    -> 2 H2
!   4601
xk(27)=2.30d-7*(T_K/300.d0)**(-0.50d0)
!  28)   He+   +   H     ->   H+    +   He
!   2938, 2939
!   GA08 (26)
xk(28)=1.2d-15*dexp(T_K/300.d0)**0.25d0
!  29)   He+   +   H-    ->   H     +   He
!   3491
!   GA08 (28)
 xk(29)=2.32d-7*(T_K/300.d0)**(-0.52d0)&
*dexp(T_K/22400.d0)
!  30)   He+   +   H2    ->   H2+   +   He
!   3060
!   GA08 (25)
xk(30)=7.20d-15
!  31)   HeH+  +   H     ->   H2+   +   He
!   550
xk(31)=9.10d-10
!  32)   HeH+  +   H2    ->   H3+   +   He
!   643
xk(32)=1.50d-9
!  33)
!   none
xk(33)=0.d0
!  34)   H3+   +   e     ->   H2    +   H
!   3522
xk(34)=2.34d-8*(T_K/300.d0)**(-0.52d0)
!  35)   H3+   +   e     -> 3 H
!   3521
xk(35)=4.36d-8*(T_K/300.d0)**(-0.52d0)
!  36)   HeH+  +   e     ->   H     +   He
!   3523
xk(36)=1.00d-8*(T_K/300.d0)**(-0.60d0)
!  37)   H2    +   He    -> 2 H     +   He
!  GA08 (11)
!     v=0
xk37v0=10.d0**(-27.029d0+3.801d0*xlgT-29487.d0/T_K)
!     LTE
xk37LTE=10.d0**(-2.729d0-1.75d0*xlgT-23474.d0/T_K)
!
xlgk37=(xcr/(1.d0+xcr))*dlog10(xk37LTE)&
+(1.d0/(1.d0+xcr))*dlog10(xk37v0)
!
xk(37)=10.d0**xlgk37
!
!  38)   H-    +   H     -> 2 H     +   e
!  GA08 (15)
if(T_eV < 0.1d0) then
   xk(38)=2.5634d-9*T_eV**1.78186d0
else
   xk(38)=dexp(-20.372609d0+(1.13944933d0&
+(-1.4210135d-1+(8.4644554d-3+(-1.4327641d-3&
+(2.0122503d-4+(8.6639632d-5+(-2.5850097d-5&
+(2.4555012d-6-8.0683825d-8*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)&
*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
endif
!  39)   H-    +   H2+   -> 3 H
!  GA08 (22)
xk(39)=1.4d-7*T300**(-0.5d0)
!  40)   H2    +   e     ->   H-    +   H
!  GA08 (23)
xk(40)=2.7d-8*(T_K**(-1.27d0))*dexp(-43000.d0/T_K)
!  41)   He    +   H+    ->   He+   +   H
!  GA08 (27)
if(T_K < 1.d4) then
   xk(41)=1.26d-9*(T_K**(-0.75d0))*dexp(-127500.d0/T_K)
else
   xk(41)=4.0d-37*T_K**4.74d0
endif
!  42)   He    +   H-    ->   He    +   H   +   e
!  GA08 (29)
xk(42)=4.1d-17*(T_K**2)*dexp(-19870.d0/T_K)
!  43) 2 H     +   He    ->   H2    +   He
!  GA08 (32)
xk(43)=6.9d-32*T_K**(-0.4d0)
xk(44:46)=0.d0
!     47-50: used for additional D reactions
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     D    reactions                           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!  51)   D+  +   e    ->   D    +   ph.
!  GA08 (33)
xk(51)=xk(2)
!  52)   D   +   H+   ->   D+   +   H
!  GA08 (34)
if(T_K < 2.d5) then
   xk(52)=2.0d-10*(T_K**0.402d0)*dexp(-37.1d0/T_K)&
-3.31d-17*(T_K**1.48d0)
else
   xk(52)=3.44d-10*(T_K**0.35d0)
endif
!  53)   D+  +   H    ->   D    +   H+
!  GA08 (35)
xk(53)=2.06d-10*(T_K**0.396d0)*dexp(-33.d0/T_K)&
+2.03d-9*(T_K**(-0.332d0))
!  54)   D   +   H    ->   HD   +   ph.
!  GA08 (36)
if(T_K < 200.d0) then
   xk(54)=1.d-25*(2.80202d0-6.63697d0*xlnT&
+4.75619d0*xlnT**2-1.39325d0*xlnT**3&
+0.178259d0*xlnT**4-0.00817097d0*xlnT**5)
else
   xk(54)=1.d-25*dexp(507.207d0-370.889d0*xlnT&
+104.854d0*xlnT**2-14.4192d0*xlnT**3&
+0.971469d0*xlnT**4-0.0258076d0*xlnT**5)
endif
!  55)   D   +   H2   ->   H    +   HD
!  GA08 (37)
if(T_K < 2000.d0) then
   xk(55)=10.d0**(-56.4737d0+5.88886d0*xlgT&
+7.19692d0*xlgT**2+2.25069d0*xlgT**3&
-2.16903d0*xlgT**4+0.317887d0*xlgT**5)
else
   xk(55)=3.17d-10*dexp(-5207.d0/T_K)
endif
!  56)   HD+ +   H    ->   H+   +   HD
!  GA08 (38)
xk(56)=xk(10)
!  57)   D+  +   H2   ->   H+   +   HD
!  GA08 (39)
xk(57)=(0.417d0+0.846d0*xlgT-0.137d0*(xlgT**2))*1.d-9
!  58)   HD  +   H    ->   H2   +   D
!  GA08 (40)
if(T_K < 200.d0) then
   xk(58)=5.25d-11*dexp(-4430.d0/T_K)
else
   xk(58)=5.25d-11*dexp(-4430.d0/T_K+173900.d0/T_K**2)
endif
!  59)   HD  +   H+   ->   H2   +   D+
!  GA08 (41)
xk(59)=1.1d-9*dexp(-488.d0/T_K)
!  60)   D   +   H+   ->   HD+  +   ph.
!  GA08 (42)
xk(60)=3.9d-19*(T300**1.8d0)*dexp(20.d0/T_K)
!  61)   D+  +   H    ->   HD+  +   ph.
!  GA08 (43)
xk(61)=3.9d-19*(T300**1.8d0)*dexp(20.d0/T_K)
!  62)   HD+ +   e    ->   H    +   D
!  GA08 (44)
xk(62)=7.2d-8*T_K**(-0.5d0)
!  63)   D   +   e    ->   D-   +   ph.
!  GA08 (51)
xk(63)=xk(7)
!  64)   D+  +   D-   ->   2D
!  GA08 (68)
xk(64)=xk(15)
!  65)   H+  +   D-   ->   D    +   H
!  GA08 (67)
xk(65)=xk(15)
!  66)   H-  +   D    ->   H    +   D-
!  GA08 (53)
xk(66)=6.4d-9*(T300**0.41d0)
!  67)   D-  +   H    ->   D    +   H-
!  GA08 (52)
xk(67)=6.4d-9*(T300**0.41d0)
!  68)   D-  +   H    ->   HD   +   e
!  GA08 (55)
xk(68)=0.5d0*xk(8)
!  69)   D   +   e    ->   D+   +  2 e
!  GA08 (45)
xk(69)=xk(1)
!  70)   He+ +   D    ->   D+   +   He
!  GA08 (46)
xk(70)=1.1d-15*T300**0.25d0
!  71)   He  +   D+   ->   D    +   He+
!  GA08 (47)
if(T_K < 10000.d0) then
   xk(71)=1.85d-9*(T_K**(-0.75d0))*dexp(-127500.d0/T_K)
else
   xk(71)=5.9d-37*T_K**4.74d0
endif
!  72)   H2+ +   D    ->  HD+   +   H
!  GA08 (48)
xk(72)=1.07d-9*(T300**0.062d0)*dexp(-T_K/41400.d0)
!  73)   HD+ +   D    ->  HD    +   D+
!  GA08 (49)
xk(73)=xk(10)
!  74)   HD+ +   H    ->  H2+   +   D
!  GA08 (50)
xk(74)=1.0d-9*dexp(-154.d0/T_K)
!  75)   HD  +   e    ->  H     +   D-
!  GA08 (57)
xk(75)=1.35d-9*(T_K**(-1.27d0))*dexp(-43000.d0/T_K)
!  76)   HD  +   e    ->  D     +   H-
!  GA08 (58)
xk(76)=1.35d-9*(T_K**(-1.27d0))*dexp(-43000.d0/T_K)
!  77)   H+  +   D-   ->  HD+   +   e
!  GA08 (60)
xk(77)=1.1d-9*T300**(-0.4d0)
!  78)   D+  +   H-   ->  HD+   +   e
!  GA08 (61)
xk(78)=1.1d-9*T300**(-0.4d0)
!  79)   D-  +   e    ->  D     + 2 e
!  GA08 (63)
xk(79)=xk(14)
!  80)   D-  +   H    ->  D     +   H    +   e
!  GA08 (64)
xk(80)=xk(38)
!  81)   D-  +   He   ->  D     +   He   +   e
!  GA08 (65)
xk(81)=1.5d-17*(T_K**2)*dexp(-19870.d0/T_K)
!  82)   D+  +   H-   ->  D     +   H
!  GA08 (66)
xk(82)=xk(15)
!  83)   H2+ +   D-   ->  H2    +   D
!  GA08 (69)
xk(83)=1.7d-7*T300**(-0.5d0)
!  84)   H2+ +   D-   -> 2 H    +   D
!  GA08 (70)
xk(84)=1.7d-7*T300**(-0.5d0)
!  85)   HD+ +   H-   ->   HD   +   H
!  GA08 (71)
xk(85)=1.5d-7*T300**(-0.5d0)
!  86)   HD+ +   H-   ->   D    + 2 H
!  GA08 (72)
xk(86)=1.5d-7*T300**(-0.5d0)
!  87)   HD+ +   D-   ->  HD    +   D
!  GA08 (73)
xk(87)=1.9d-7*T300**(-0.5d0)
!  88)   HD+ +   D-   ->  2 D   +   H
!  GA08 (74)
xk(88)=1.9d-7*T300**(-0.5d0)
!  89)   He+ +   D-   ->  He    +   D
!  GA08 (79)
xk(89)=3.03d-7*(T300**(-0.52d0))*dexp(T_K/22400.d0)
!  90)   D   +   H2+  ->  H2    +   D+
!  GA08 (81)
xk(90)=xk(10)
!  91)  H2+  +   D    ->  HD    +   H+
!  GA08 (82)
xk(91)=1.0d-9
!  92)  HD+  +   H    ->  H2    +   D+
!  GA08 (83)
xk(92)=1.0d-9
!  93)  H2   +   D+   ->  H2+   +   D
!  GA08 (90)
xk(93)=xk(11)
!  94)  H2   +   D+   ->  HD+   +   H
!  GA08 (91)
xk(94)=(1.04d-9+9.52d-9*(T_K/10000.d0)&
-1.81d-9*(T_K/10000.d0)**2)*dexp(-21000.d0/T_K)
!  95)  HD   +   H+   ->  HD+   +   H
!  GA08 (92)
xk(95)=xk(11)
!  96)  HD   +   H+   ->  H2+   +   D
!  GA08 (93)
xk(96)=1.0d-9*dexp(-21600.d0/T_K)
! 97)  HD   +   D+   ->  HD+   +   D
! GA08 (94)
xk(97)=xk(11)
! 98)  HD   +   He+  ->  HD+   +   He
! GA08 (101)
xk(98)=xk(30)
! 99)  HD   +   He+  ->  He    +   H+   +   D
! GA08 (102)
xk(99)=1.85d-14*dexp(-35.d0/T_K)
!100)  HD   +   He+  ->  He    +   H    +   D+
! GA08 (103)
xk(100)=1.85d-14*dexp(-35.d0/T_K)
!D50)  HD   +   H    -> 2 H    +   D
! GA08 (108)
xk(50)=xk(13)
!D49)  HD   +   H2   ->  H     +   D    +   H2
! GA08 (109)
xk(49)=xk(21)
!D48)  HD   +   He   ->  H     +   D    +   He
! GA08 (110)
xk(48)=xk(37)
!D47)  HD   +   e    ->  H     +   D    +   e
! GA08 (111)
!    v=0
xk47v0=5.09d-9*(T_K**0.128d0)*dexp(-103258.d0/T_K)
!    LTE
xk47LTE=1.04d-9*(T_K**0.218d0)*dexp(-53070.7d0/T_K)
xlgk47=(xcr/(1.d0+xcr))*dlog10(xk47LTE)&
+(1.d0/(1.d0+xcr))*dlog10(xk47v0)
xk(47)=10.d0**xlgk47

!    metal reactions
! 101)   H     +   CH    ->   C     +   H2
!  1
xk(101)= 1.31d-10*dexp(-80.d0/T_K)
! 102)   H     +   CH    ->   C     + 2 H
!  4585
xk(102)=6.00d-09*dexp(-40200.d0/T_K)

! 103)   H     +   CH2   ->   CH    +   H2
!  2
xk(103)=6.64d-11
! 104)   H     +   CH3   ->   CH2   +   H2
!  4
xk(104)=1.00d-10*dexp(-7600.d0/T_K)
! 105)   H     +   CH4   ->   H2    +   CH3
!  7
xk(105)=5.94d-13*(T_K/300.d0)**3.d0*dexp(-4045.d0/T_K)
! 106)   H     +   OH    ->   H2    +   O
!  8
xk(106)=6.99d-14*(T_K/300.d0)**2.80d0*dexp(-1950.d0/T_K)
! 107)   H     +   OH    ->   O     + 2 H
!  4586
xk(107)=6.00d-9*dexp(-50900.d0/T_K)
! 108)   H     +   H2O   ->   OH    +   H2
!  10
xk(108)=1.59d-11*(T_K/300.d0)**1.2d0*dexp(-9610.d0/T_K)
! 109)   H     +   H2O   ->   OH    + 2 H
!  4587
xk(109)=5.80d-9*dexp(-52900.d0/T_K)
! 110)   H     +   C2    ->   CH    +   C
!  11
xk(110)=4.67d-10*(T_K/300.d0)**0.50d0*dexp(-30450.d0/T_K)
! 111)   H     +   CO    ->   C    +   OH
!  15
xk(111)=1.10d-10*(T_K/300.d0)**0.50d0*dexp(-77700.d0/T_K)
! 112)   H     +   H2CO  ->   HCO   +   H2
!  21
xk(112)=4.85d-12*(T_K/300.d0)**1.9d0*dexp(-1379.d0/T_K)
! 113)   H     +   O2    ->   OH    +   O
!  25
xk(113)=2.61d-10*dexp(-8156.d0/T_K)
! 114)   H     +   O2    -> 2 O     +   H
!  4591
xk(114)=6.00d-9*dexp(-52300.d0/T_K)
! 115)   H     +   O2H   ->   H2O   +   O
!  27
xk(115)=5.00d-11*dexp(-866.d0/T_K)
! 116)   H     +   O2H   ->   H2    +   O2
!  28
xk(116)=2.06d-11*(T_K/300.d0)**0.84d0*dexp(-277.d0/T_K)
! 117)   H     +   O2H   -> 2 OH
!  26
xk(117)=1.66d-10*dexp(-413.d0/T_K)
! 118)   H     +   H2O2  ->   O2H   +   H2
!  31
xk(118)=8.00d-11*dexp(-4000.d0/T_K)
! 119)   H     +   CO2   ->   OH    +   CO
!  36
xk(119)=3.38d-10*dexp(-13163.d0/T_K)
! 120)   C     +   H2    ->   CH    +   H
!  47
xk(120)=6.64d-10*dexp(-11700.d0/T_K)
! 121)   C     +   OH    ->   CH    +   O
!  72
xk(121)=2.25d-11*(T_K/300.d0)**0.50d0*dexp(-14800.d0/T_K)
! 122)   C     +   CO    ->   C2    +   O
!  79
xk(122)=2.94d-11*(T_K/300.d0)**0.50d0*dexp(-58025.d0/T_K)
! 123)   O     +   H2    ->   OH    +   H
!  53
xk(123)=3.14d-13*(T_K/300.d0)**2.7d0*dexp(-3150.d0/T_K)
! 124)   O     +   CH    ->   OH    +   C
! 138
xk(124)=2.52d-11*dexp(-2381.d0/T_K)
! 125)   O     +   CH2   ->   CH    +   OH
! 240
xk(125)=4.98d-10*dexp(-6000.d0/T_K)
! 126)   O     +   CH2   ->   HCO   +   H
! 243
xk(126)=5.01d-11
! 127)   O     +   CH4   ->   CH3   +   OH
! 311
xk(127)=2.29d-12*(T_K/300.d0)**2.2d0*dexp(-3820.d0/T_K)
! 128)   O     +   H2O   -> 2 OH
! 314
xk(128)=1.85d-11*(T_K/300.d0)**0.95d0*dexp(-8571.d0/T_K)
! 129)   O     +   H2CO  ->   HCO   +   OH
! 339
xk(129)=1.07d-11*(T_K/300.d0)**1.17d0*dexp(-1242.d0/T_K)
! 130)   O     +   H2O2  ->   O2H   +   OH
! 354
xk(130)=8.54d-14*(T_K/300.d0)**3.25d0*dexp(-1200.d0/T_K)
! 131)   O     +   CO2   ->   CO    +   O2
! 370
xk(131)=2.46d-11*dexp(-26567.d0/T_K)
! 132)   H+    +   O     ->   O+    +   H
! 2945, 2946
if(T_K <1.d4) then
   xk(132)=7.31d-10*(T_K/300.d0)**0.23d0*dexp(-225.9d0/T_K)
else
   xk(132)=3.04d-10*(T_K/300.d0)**0.47d0*dexp(11.5d0/T_K)
endif
! 133)   H2    +   CH    ->   CH2   +   H
!  48
xk(133)=5.46d-10*dexp(-1943.d0/T_K)
! 134)   H2    +   CH    ->   C    +   H   +   H2
!  4594
xk(134)=6.00d-9*dexp(-40200.d0/T_K)
! 135)   H2    +   CH2   ->   CH3   +   H
!  50
xk(135)=5.18d-11*(T_K/300.d0)**0.17d0*dexp(-6400.d0/T_K)
! 136)   H2    +   CH3   ->   CH4   +   H
!  52
xk(136)=6.86d-14*(T_K/300.d0)**2.74d0*dexp(-4740.d0/T_K)
! 137)   H2    +   OH    ->   H2O   +   H
!  55
xk(137)=2.05d-12*(T_K/300.d0)**1.52d0*dexp(-1736.d0/T_K)
! 138)   H2    +   OH    ->   O     +   H     +   H2
!  4595
xk(138)=6.00d-9*dexp(-50900.d0/T_K)
! 139)   H2    +   H2O   ->   OH    +   H     +   H2
!  4596
xk(139)=5.80d-9*dexp(-52900.d0/T_K)
! 140)   H2    +   O2    ->   O2H   +   H
!  59
xk(140)=2.40d-10*dexp(-28500.d0/T_K)
! 141)   H2    +   O2    -> 2 OH
!  58
xk(141)=3.16d-10*dexp(-21890.d0/T_K)
! 142)   H2    +   O2    -> 2 O     +   H2
!  4598
xk(142)=6.00d-9*dexp(-52300.d0/T_K)
! 143)   H2    +   O2H   ->   H2O2  +   H
!  61
xk(143)=4.38d-12*dexp(-10751.d0/T_K)
! 144)
!  none
xk(144)=0.d0
! 145)   H3+   +   O2    ->   O2H+  +   H2
!  788
xk(145)=9.30d-10*dexp(-100.d0/T_K)
! 146)   C+    +   H2    ->   CH+   +   H
!  645
xk(146)=1.00d-10*dexp(-4640.d0/T_K)
! 147)   CH    +   CH4   ->   CH2   +   CH3
!  141
xk(147)=2.28d-11*(T_K/300.d0)**0.70d0*dexp(-3000.d0/T_K)
! 148)   CH    +   OH    ->   HCO   +   H
!  143
xk(148)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-5000.d0/T_K)
! 149)   CH    +   HCO   ->   CH2   +   CO
!  147
xk(149)=2.87d-12*(T_K/300.d0)**0.70d0*dexp(-500.d0/T_K)
! 150)   CH    +   H2CO  ->   CH2   +   HCO
!  152
xk(150)=9.21d-12*(T_K/300.d0)**0.70d0*dexp(-2000.d0/T_K)
! 151)   CH    +   O2    ->   HCO   +   O
!  155
xk(151)=1.44d-11*(T_K/300.d0)**0.70d0*dexp(-3000.d0/T_K)
! 152)   CH    +   O2H   ->   CH2   +   O2
!  159
xk(152)=2.94d-13*(T_K/300.d0)**0.50d0*dexp(-7550.d0/T_K)
! 153)   CH    +   O2H   ->   HCO   +   OH
!  160
xk(153)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
! 154)   CH    +   CO2   ->   HCO   +   CO
!  161
xk(154)=2.94d-13*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
! 155) 2 CH2             ->   CH    +   CH3
!  236
xk(155)=4.00d-10*dexp(-5000.d0/T_K)
! 156)   CH2   +   CH4   -> 2 CH3
!  244
xk(156)=7.13d-12*dexp(-5050.d0/T_K)
! 157)   CH2   +   HCO   ->   CH3   +   CO
!  251
xk(157)=3.00d-11
! 158)   CH2   +   H2CO  ->   CH3   +   HCO
!  254
xk(158)=3.30d-13*dexp(-3270.d0/T_K)
! 159)   CH2+  +   H     ->   CH+   +   H2
!  553
xk(159)=1.00d-9*dexp(-7080.d0/T_K)
! 160)   CH3   +   H2CO  ->   CH4   +   HCO
!  299
xk(160)=1.34d-15*(T_K/300.d0)**5.05d0*dexp(-1636.d0/T_K)
! 161)   CH3+  +   H     ->   CH2+  +   H2
!  554
xk(161)=7.00d-10*dexp(-10560.d0/T_K)
! 162)   OH    +   CH2   ->   CH3   +   O
!  245
xk(162)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
! 163)   OH    +   CH2   ->   H2O   +   CH
!  246
xk(163)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
!164)   OH    +   CH2   ->   H2CO  +   H
! 247
xk(164)=3.00d-11
!165)   OH    +   CH3   ->   H2O   +   CH2
! 290
xk(165)=1.20d-10*dexp(-1400.d0/T_K)
!166)   OH    +   CH4   ->   H2O   +   CH3
! 422
xk(166)=3.77d-13*(T_K/300.d0)**2.42d0*dexp(-1162.d0/T_K)
!167) 2 OH              ->   H2O   +   O
! 426
xk(167)=1.65d-12*(T_K/300.d0)**1.14d0*dexp(-50.d0/T_K)
!168)   OH    +   CO    ->   CO2   +   H
! 437
xk(168)=2.81d-13*dexp(-176.d0/T_K)
!169)   OH    +   H2O2  ->   H2O   +   O2H
! 449
xk(169)=5.26d-12*dexp(-307.d0/T_K)
!170)   H2O   +   CH3   ->   CH4   +   OH
! 294
xk(170)=2.30d-15*(T_K/300.d0)**3.47d0*dexp(-6681.d0/T_K)
!171)   CO    +   O2    ->   CO2   +   O
! 504
xk(171)=5.99d-12*dexp(-24075.d0/T_K)
!172)   CO    +   O2H   ->   CO2   +   OH
! 505
xk(172)=5.60d-10*dexp(-12160.d0/T_K)
!173)   O2    +   CH2   ->   HCO   +   OH
! 258
xk(173)=4.10d-11*dexp(-750.d0/T_K)
!174)   O2    +   CH2   ->   H2CO  +   O
! 259
xk(174)=3.65d-11*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
!175)   O2    +   CH3   ->   O2H   +   CH2
! 303
xk(175)=5.30d-12*dexp(-34975.d0/T_K)
!176)   O2    +   CH3   ->   H2CO  +   OH
! 302
xk(176)=5.64d-13*dexp(-4500.d0/T_K)
!177)   O2    +   CH4   ->   CH3   +   O2H
! 424
xk(177)=6.70d-11*dexp(-28640.d0/T_K)
!178)   O2    +   HCO   ->   O2H   +   CO
! 520
if(T_K > 200.d0) then
   xk(178)=1.58d-12*(T_K/300.d0)**1.24d0*dexp(353.d0/T_K)
else
   xk(178)=5.58d-12
endif
!179)   O2H   +   CH3   ->   CH4   +   O2
! 306
xk(179)=6.00d-12
!180)   O2H   +   H2O   ->   H2O2  +   OH
! 460
xk(180)=4.65d-11*dexp(-16500.d0/T_K)
!181)   O2H   +   HCO   ->   H2CO  +   O2
! 521
xk(181)=5.00d-11
!182)   O2H   +   H2CO  ->   H2O2  +   HCO
! 531
xk(182)=3.30d-12*dexp(-5870.d0/T_K)
!183) 2 O2H             ->   H2O2  +   O2
! 540
if(T_K > 200.d0) then
   xk(183)=3.81d-15*(T_K/300.d0)**3.71d0*dexp(1761.d0/T_K)
else
   xk(183)=5.64d-12
endif
!184)   H     +   HCO   ->   CO    +   H2
! 18
xk(184)=2.00d-10
!185)   C     +   H     ->   CH    +   ph.
! 4081
xk(185)=1.00d-17
!186) 2 C               ->   C2    +   ph.
! 4106
xk(186)=4.36d-18*(T_K/300.d0)**0.35d0*dexp(-161.3d0/T_K)
!187)   C     +   O     ->   CO    +   ph.
! 4112, 4113
if(T_K > 300.d0) then
   xk(187)=3.09d-17*(T_K/300.d0)**0.33d0*dexp(-1629.d0/T_K)
else
   xk(187)=2.10d-19
endif
!188)   C     +   CH    ->   C2    +   H
! 63
xk(188)=6.59d-11
!189)   C     +   OH    ->   CO    +   H
! 73
xk(189)=1.00d-10
!190)   C     +   HCO   ->   CH    +   CO
! 83
xk(190)=1.00d-10
!191)   C     +   O2    ->   CO    +   O
! 92, 93
if(T_K < 295.d0) then
   xk(191)=4.70d-11*(T_K/300.d0)**(-0.34d0)
else
   xk(191)=2.48d-12*(T_K/300.d0)**1.54d0*dexp(613.d0/T_K)
endif
!192)   CH3   +   HCO   ->   CH4   +   CO
! 297
xk(192)=2.00d-10
!193)   O     +   H     ->   OH    +   ph.
! 4083
xk(193)=9.90d-19*(T_K/300.d0)**(-0.38d0)
!194) 2 O               ->   O2    +   ph.
! 4136
xk(194)=4.90d-20*(T_K/300.d0)**1.58d0
!195)   O     +   CH    ->   CO    +   H
! 139, 140
if(T_K < 2000.d0) then
   xk(195)=6.60d-11
else
   xk(195)=1.02d-10*dexp(-914.d0/T_K)
endif
!196)   O     +   CH    ->   HCO+  +   e
! 4600
xk(196)=2.00d-11*(T_K/300.d0)**0.44d0
!197)   O     +   CH2   ->   CO    + 2 H
! 241
xk(197)=1.33d-10
!198)   O     +   CH3   ->   H2CO  +   H
! 286
xk(198)=1.30d-10
!199)   O     +   OH    ->   O2    +   H
! 312
if(T_K > 158.d0) then
   xk(199)=1.77d-11*dexp(178.d0/T_K)
else
   xk(199)=5.46d-11
endif
!200)   O     +   C2    ->   CO    +   C
! 315, 316
if(T_K < 300.d0) then
   xk(200)=1.d-10
else
   xk(200)=6.d-10
endif
!201)   O     +   HCO   ->   CO2   +   H
! 332
xk(201)=5.00d-11
!202)   O     +   HCO   ->   OH    +   CO
! 333
xk(202)=5.00d-11
!203)   O     +   O2H   ->   OH    +   O2
! 351
if(T_K > 200.d0) then
   xk(203)=3.17d-11*dexp(174.d0/T_K)
else
   xk(203)=7.57d-11
endif
!204)   OH    +   HCO   ->   H2O   +   CO
! 439
xk(204)=1.70d-10
!205)   OH    +   H2CO  ->   H2O   +   HCO
! 443
if(T_K > 200.d0) then
   xk(205)=2.22d-12*(T_K/300.d0)**1.42d0*dexp(416.d0/T_K)
else
   xk(205)=9.99d-12
endif
!206)   OH    +   O2H   ->   H2O   +   O2
! 448
if(T_K > 200.d0) then
   xk(206)=3.66d-11*(T_K/300.d0)**(-0.13d0)*dexp(244.d0/T_K)
else
   xk(206)=1.31d-10
endif
!207) 2 HCO             ->   H2CO  +   CO
! 516
xk(207)=3.00d-11
!208)   H+    +   CH    ->   CH+   +   H
! 2940
xk(208)=1.90d-9*(T_K/300.d0)**(-0.5d0)
!209)   H+    +   CH2   ->   CH+   +   H2
! 552
xk(209)=1.40d-9
!210)   H+    +   CH2   ->   CH2+  +   H
! 2941
xk(210)=1.40d-9
!211)   H+    +   CH3   ->   CH3+  +   H
! 2943
xk(211)=3.40d-9
!212)   H+    +   CH4   ->   CH3+  +   H2
! 555
xk(212)=2.30d-9
!213)   H+    +   CH4   ->   CH4+  +   H
! 2948
xk(213)=1.50d-9
!214)   H+    +   OH    ->   OH+   +   H
! 2949
xk(214)=2.10d-9*(T_K/300.d0)**(-0.5d0)
!215)   H+    +   H2O   ->   H2O+  +   H
! 2951
xk(215)=6.90d-9*(T_K/300.d0)**(-0.5d0)
!216)   H+    +   C2    ->   C2+   +   H
! 2952
xk(216)=3.10d-9
!217)   H+    +   HCO   ->   CO+   +   H2
! 566
xk(217)=9.40d-10*(T_K/300.d0)**(-0.5d0)
!218)   H+    +   HCO   ->   H2+   +   CO
! 567
xk(218)=9.40d-10*(T_K/300.d0)**(-0.5d0)
!219)   H+    +   HCO   ->   HCO+  +   H
! 2963
xk(219)=9.40d-10*(T_K/300.d0)**(-0.5d0)
!220)   H+    +   H2CO  ->   H2CO+ +   H
! 2967
xk(220)=2.96d-9*(T_K/300.d0)**(-0.5d0)
!221)   H+    +   H2CO  ->   HCO+  +   H2
! 575
xk(221)=3.57d-9*(T_K/300.d0)**(-0.5d0)
!222)   H+    +   O2    ->   O2+   +   H
! 2972
xk(222)=2.00d-9
!223)   H+    +   CO2   ->   HCO+  +   O
! 606
xk(223)=3.50d-9
!224)   H-    +   C     ->   CH    +   e
! 4032
xk(224)=1.00d-9
!225)   H-    +   O     ->   OH    +   e
! 4039
xk(225)=1.00d-9
!226)   H-    +   CH    ->   CH2   +   e
! 4034
xk(226)=1.00d-10
!227)   H-    +   CH2   ->   CH3   +   e
! 4036
xk(227)=1.00d-9
!228)   H-    +   CH3   ->   CH4   +   e
! 4038
xk(228)=1.00d-9
!229)   H-    +   OH    ->   H2O   +   e
! 4042
xk(229)=1.00d-10
!230)   H-    +   CO    ->   HCO   +   e
! 4048
xk(230)=5.00d-11
!231)   H-    +   HCO   ->   H2CO  +   e
! 4049
xk(231)=1.00d-9
!232)   H2+   +   C     ->   CH+   +   H
! 644
xk(232)=2.40d-9
!233)   H2+   +   O     ->   OH+   +   H
! 655
xk(233)=1.50d-9
!234)   H2+   +   CH    ->   CH+   +   H2
! 3061
xk(234)=7.10d-10*(T_K/300.d0)**(-0.5d0)
!235)   H2+   +   CH    ->   CH2+  +   H
! 646
xk(235)=7.10d-10*(T_K/300.d0)**(-0.5d0)
!236)   H2+   +   CH2   ->   CH3+  +   H
! 650
xk(236)=1.00d-9
!237)   H2+   +   CH2   ->   CH2+  +   H2
! 3062
xk(237)=1.00d-9
!238)   H2+   +   CH4   ->   CH4+  +   H2
! 3065
xk(238)=1.40d-9
!239)   H2+   +   CH4   ->   CH5+  +   H
! 659
xk(239)=1.14d-10
!240)   H2+   +   CH4   ->   CH3+  +   H     +  H2
! 662
xk(240)=2.30d-9
!241)   H2+   +   OH    ->   OH+   +   H2
! 3066
xk(241)=7.60d-10*(T_K/300.d0)**(-0.5d0)
!242)   H2+   +   OH    ->   H2O+  +   H
! 663
xk(242)=7.60d-10*(T_K/300.d0)**(-0.5d0)
!243)   H2+   +   H2O   ->   H2O+  +   H2
! 3068
xk(243)=3.90d-9*(T_K/300.d0)**(-0.5d0)
!244)   H2+   +   H2O   ->   H3O+  +   H
! 668
xk(244)=3.40d-9*(T_K/300.d0)**(-0.5d0)
!245)   H2+   +   C2    ->   C2+   +   H2
! 3069
xk(245)=1.10d-9
!246)   H2+   +   CO    ->   CO+   +   H2
! 3074
xk(246)=6.44d-10
!247)   H2+   +   CO    ->   HCO+  +   H
! 682
xk(247)=2.16d-9
!248)   H2+   +   HCO   ->   H3+   +   CO
! 689
xk(248)=1.00d-9*(T_K/300.d0)**(-0.5d0)
!249)   H2+   +   HCO   ->   HCO+  +   H2
! 3076
xk(249)=1.00d-9*(T_K/300.d0)**(-0.5d0)
!250)   H2+   +   H2CO  ->   H2CO+ +   H2
! 3078
xk(250)=1.40d-9*(T_K/300.d0)**(-0.5d0)
!251)   H2+   +   H2CO  ->   HCO+  +   H     +   H2
! 691
xk(251)=1.40d-9*(T_K/300.d0)**(-0.5d0)
!252)   H2+   +   O2    ->   O2H+  +   H
! 694
xk(252)=1.90d-9
!253)   H2+   +   O2    ->   O2+   +   H2
! 3080
xk(253)=8.00d-10
!254)   H2+   +   CO2   ->   HCO2+ +   H
! 713
xk(254)=2.35d-9
!255)   H3+   +   C     ->   CH+   +   H2
! 749
xk(255)=2.00d-9
!256)   H3+   +   O     ->   OH+   +   H2
! 754
xk(256)=8.40d-10
!257)   H3+   +   CH    ->   CH2+  +   H2
! 750
xk(257)=1.20d-9*(T_K/300.d0)**(-0.5d0)
!258)   H3+   +   CH2   ->   CH3+  +   H2
! 751
xk(258)=1.70d-9
!259)   H3+   +   CH3   ->   CH4+  +   H2
! 753
xk(259)=2.10d-9
!260)   H3+   +   CH4   ->   CH5+  +   H2
! 757
xk(260)=2.40d-9
!261)   H3+   +   OH    ->   H2O+  +   H2
! 758
xk(261)=1.30d-9*(T_K/300.d0)**(-0.5d0)
!262)   H3+   +   H2O   ->   H3O+  +   H2
! 760
xk(262)=5.90d-9*(T_K/300.d0)**(-0.5d0)
!263)   H3+   +   CO    ->   HCO+  +   H2
! 771
xk(263)=1.70d-9
!264)   H3+   +   HCO   ->   H2CO+ +   H2
! 777
xk(264)=1.70d-9*(T_K/300.d0)**(-0.5d0)
!265)   H3+   +   H2CO  ->   H3CO+ +   H2
! 781
xk(265)=6.30d-9*(T_K/300.d0)**(-0.5d0)
!266)   H3+   +   CO2   ->   HCO2+ +   H2
! 818
xk(266)=2.00d-9
!267)   He+   +   CH    ->   C+    +   H     +   He
! 913
xk(267)=1.10d-9*(T_K/300.d0)**(-0.5d0)
!268)   He+   +   CH    ->   CH+   +   He
! 3083
xk(268)=5.00d-10*(T_K/300.d0)**(-0.5d0)
!269)   He+   +   CH2   ->   CH+   +   H     +   Hec
! 915
xk(269)=7.50d-10
!270)   He+   +   CH2   ->   C+    +   H2    +   He
! 914
xk(270)=7.50d-10
!271)   He+   +   CH3   ->   CH+   +   H2    +   He
! 917
xk(271)=1.80d-9
!272)   He+   +   CH4   ->   CH3   +   H+    +   He
! 922
xk(272)=4.80d-10
!273)   He+   +   CH4   ->   CH+   +   H2    +   He    +   H
! 920
xk(273)=2.40d-10
!274)   He+   +   CH4   ->   CH2+  +   H2    +   He
! 921
xk(274)=9.50d-10
!275)   He+   +   CH4   ->   CH3+  +   He    +   H
! 923
xk(275)=8.50d-11
!276)   He+   +   CH4   ->   CH4+  +   He
! 3084
xk(276)=5.10d-11
!277)   He+   +   OH    ->   O+    +   H     +   He
! 924
xk(277)=1.10d-9*(T_K/300.d0)**(-0.5d0)
!278)   He+   +   H2O   ->   H+    +   OH    +   He
! 927
xk(278)=2.04d-10*(T_K/300.d0)**(-0.5d0)
!279)   He+   +   H2O   ->   OH+   +   H     +   He
! 928
xk(279)=2.86d-10*(T_K/300.d0)**(-0.5d0)
!280)   He+   +   H2O   ->   H2O+  +   He
! 3086
xk(280)=6.05d-11*(T_K/300.d0)**(-0.5d0)
!281)   He+   +   C2    ->   C+    +   C    +   He
! 930
xk(281)=1.60d-9
!282)   He+   +   C2    ->   C2+   +   He
! 3087
xk(282)=5.00d-10
!283)   He+   +   CO    ->   C+    +   O     +   He
! 948, 949
xk(283)=1.60d-9
!     xk(283)=1.40d-9*(T_K/300.d0)**(-0.5d0)
!284)   He+   +   HCO   ->   CO+   +   H     +   He
! 957
xk(284)=4.90d-10*(T_K/300.d0)**(-0.5d0)
!285)   He+   +   HCO   ->   CH+   +   O     +   He
! 955
xk(285)=4.90d-10*(T_K/300.d0)**(-0.5d0)
!286)   He+   +   HCO   ->   HeH+  +   CO
! 956
xk(286)=3.00d-10*(T_K/300.d0)**(-0.5d0)
!287)   He+   +   H2CO  ->   CO+   +   H2    +   He
! 965
xk(287)=1.88d-9*(T_K/300.d0)**(-0.5d0)
!288)   He+   +   H2CO  ->   HCO+  +   H     +   He
! 966
xk(288)=1.14d-9*(T_K/300.d0)**(-0.5d0)
!289)   He+   +   O2    ->   O+    +   O     +   He
! 977
xk(289)=1.10d-9
!290)   He+   +   O2    ->   O2+   +   He
! 3095
xk(290)=3.30d-11
!291)   He+   +   CO2   ->   O2+   +   C    +   He
! 1025
xk(291)=1.10d-11
!292)   He+   +   CO2   ->   O+    +   CO    +   He
! 1023
xk(292)=1.00d-10
!293)   He+   +   CO2   ->   CO+   +   O     +   He
! 1024
xk(293)=8.70d-10
!294)   He+   +   CO2   ->   C+    +   O2    +   He
! 1026
xk(294)=4.00d-11
!295)   C+    +   H     ->   CH+   +   ph.
! 4082
xk(295)=1.70d-17
!296)   C+    +   O     ->   CO+   +   ph.
! 4114, 4115
if(T_K < 300.d0) then
   xk(296)=2.50d-18
else
   xk(296)=3.14d-18*(T_K/300.d0)**(-0.15d0)*dexp(-68.d0/T_K)
endif
!297)   C+    +   H-    ->   H     +   C
! 3493
xk(297)=2.30d-7*(T_K/300.d0)**(-0.50d0)
!298)   C+    +   H2    ->   CH2+  +   ph.
! 4089
xk(298)=4.00d-16*(T_K/300.d0)**(-0.20d0)
!299)   C+    +   CH    ->   C2+   +   H
! 1173
xk(299)=3.80d-10*(T_K/300.d0)**(-0.50d0)
!300)   C+    +   CH    ->   CH+   +   C
! 3100
xk(300)=3.80d-10*(T_K/300.d0)**(-0.50d0)
!301)   C+    +   CH2   ->   CH2+  +   C
! 3101
xk(301)=5.20d-10
!302)   C+    +   OH    ->   CO+   +   H
! 1186
xk(302)=7.70d-10*(T_K/300.d0)**(-0.50d0)
!303)   C+    +   H2O   ->   HCO+  +   H
! 1191
xk(303)=9.00d-10*(T_K/300.d0)**(-0.50d0)
!304)   C+    +   HCO   ->   HCO+  +   C
! 3111
xk(304)=4.80d-10*(T_K/300.d0)**(-0.50d0)
!305)   C+    +   HCO   ->   CH+   +   CO
! 1212
xk(305)=4.80d-10*(T_K/300.d0)**(-0.50d0)
!306)   C+    +   H2CO  ->   CH2+  +   CO
! 1225
xk(306)=2.34d-9*(T_K/300.d0)**(-0.50d0)
!307)   C+    +   H2CO  ->   HCO+  +   CH
! 1226
xk(307)=7.80d-10*(T_K/300.d0)**(-0.50d0)
!308)   C+    +   H2CO  ->   H2CO+ +   C
! 3114
xk(308)=7.80d-10*(T_K/300.d0)**(-0.50d0)
!309)   C+    +   O2    ->   CO+   +   O
! 1234
xk(309)=3.80d-10
!310)   C+    +   O2    ->   O+    +   CO
! 1235
xk(310)=6.20d-10
!311)   C+    +   CO2   ->   CO+   +   CO
! 1287
xk(311)=1.10d-9
!312)   CH+   +   H     ->   C+    +   H2
! 551
xk(312)=7.50d-10
!313)   CH+   +   C     ->   C2+   +   H
! 1174
xk(313)=1.20d-9
!314)   CH+   +   O     ->   CO+   +   H
! 1404
xk(314)=3.50d-10
!315)   CH+   +   H2    ->   CH2+  +   H
! 647
xk(315)=1.20d-9
!316)   CH+   +   CH    ->   C2+   +   H2
! 1397
xk(316)=7.40d-10*(T_K/300.d0)**(-0.50d0)
!317)   CH+   +   OH    ->   CO+   +   H2
! 1412
xk(317)=7.50d-10*(T_K/300.d0)**(-0.50d0)
!318)   CH+   +   H2O   ->   HCO+  +   H2
! 1420
xk(318)=2.90d-9*(T_K/300.d0)**(-0.50d0)
!319)   CH+   +   H2O   ->   H3O+  +   C
! 1417
xk(319)=5.80d-10*(T_K/300.d0)**(-0.50d0)
!320)   CH+   +   H2O   ->   H2CO+ +   H
! 1419
xk(320)=5.80d-10*(T_K/300.d0)**(-0.50d0)
!321)
! none
xk(321)=0.d0
!322)   CH+   +   HCO   ->   CH2+  +   CO
! 1440
xk(322)=4.60d-10*(T_K/300.d0)**(-0.50d0)
!323)   CH+   +   HCO   ->   HCO+  +   CH
! 3174
xk(323)=4.60d-10*(T_K/300.d0)**(-0.50d0)
!324)   CH+   +   H2CO  ->   CH3+  +   CO
! 1444
xk(324)=9.60d-10*(T_K/300.d0)**(-0.50d0)
!325)   CH+   +   H2CO  ->   HCO+  +   CH2
! 1446
xk(325)=9.60d-10*(T_K/300.d0)**(-0.50d0)
!326)   CH+   +   H2CO  ->   H3CO+ +   C
! 1445
xk(326)=9.60d-10*(T_K/300.d0)**(-0.50d0)
!327)   CH+   +   O2    ->   HCO+  +   O
! 1452
xk(327)=9.70d-10
!328)   CH+   +   O2    ->   HCO   +   O+
! 1453
xk(328)=1.00d-11
!329)   CH+   +   O2    ->   CO+   +   OH
! 1451
xk(329)=1.00d-11
!330)   CH+   +   CO2   ->   HCO+  +   CO
! 1465
xk(330)=1.60d-9
!331)   CH2+  +   O     ->   HCO+  +   H
! 1554
xk(331)=7.50d-10
!332)   CH2+  +   H2    ->   CH3+  +   H
! 651
xk(332)=1.60d-9
!333)   CH2+  +   H2O   ->   H3CO+ +   H
! 1563
xk(333)=1.20d-9*(T_K/300.d0)**(-0.50d0)
!334)   CH2+  +   HCO   ->   CH3+  +   CO
! 1577
xk(334)=4.50d-10*(T_K/300.d0)**(-0.50d0)
!335)   CH2+  +   H2CO  ->   HCO+  +   CH3
! 1581
xk(335)=2.81d-9*(T_K/300.d0)**(-0.50d0)
!336)   CH2+  +   O2    ->   HCO+  +   OH
! 1586
xk(336)=9.10d-10
!337)   CH2+  +   CO2   ->   H2CO+ +   CO
! 1593
xk(337)=1.60d-9
!338)   CH3+  +   O     ->   HCO+  +   H2
! 1649
xk(338)=4.00d-10
!339)   CH3+  +   O     ->   H2CO+ +   H
! 1648
xk(339)=4.00d-11
!340)   CH3+  +   H2    ->   CH5+  +   ph.
! 4091
xk(340)=1.30d-14*(T_K/300.d0)**(-1.00d0)
!341)   CH3+  +   OH    ->   H2CO+ +   H2
! 1652
xk(341)=7.20d-10*(T_K/300.d0)**(-0.50d0)
!342)   CH3+  +   HCO   ->   CH4+  +   CO
! 1669
xk(342)=4.40d-10*(T_K/300.d0)**(-0.50d0)
!343)   CH3+  +   HCO   ->   HCO+  +   CH3
! 3227
xk(343)=4.40d-10*(T_K/300.d0)**(-0.50d0)
!344)   CH3+  +   H2CO  ->   HCO+  +   CH4
! 1674
xk(344)=1.60d-9*(T_K/300.d0)**(-0.50d0)
!345)   CH3+  +   O2    ->   H3CO+ +   O
! 1676
xk(345)=5.00d-12
!346)   O+    +   H     ->   H+    +   O
! 2944
xk(346)=5.66d-10*(T_K/300.d0)**0.36d0*dexp(8.6d0/T_K)
!347)   O+    +   H-    ->   H     +   O
! 3495
xk(347)=2.30d-7*(T_K/300.d0)**(-0.50d0)
!348)   O+    +   H2    ->   OH+   +   H
! 656
xk(348)=1.70d-9
!349)   O+    +   CH    ->   CH+   +   O
! 3162
xk(349)=3.50d-10*(T_K/300.d0)**(-0.50d0)
!350)   O+    +   CH    ->   CO+   +   H
! 1405
xk(350)=3.50d-10*(T_K/300.d0)**(-0.50d0)
!351)   O+    +   CH2   ->   CH2+  +   O
! 3203
xk(351)=9.70d-10
!352)   O+    +   CH4   ->   CH3+  +   OH
! 1711
xk(352)=1.10d-10
!353)   O+    +   CH4   ->   CH4+  +   O
! 3232
xk(353)=8.90d-10
!354)   O+    +   OH    ->   OH+   +   O
! 3233
xk(354)=3.60d-10*(T_K/300.d0)**(-0.50d0)
!355)   O+    +   OH    ->   O2+   +   H
! 1714
xk(355)=3.60d-10*(T_K/300.d0)**(-0.50d0)
!356)   O+    +   H2O   ->   H2O+  +   O
! 3235
xk(356)=3.20d-9*(T_K/300.d0)**(-0.50d0)
!357)   O+    +   C2    ->   C2+   +   O
! 3236
xk(357)=4.80d-10
!358)   O+    +   C2    ->   CO+   +   C
! 1720
xk(358)=4.80d-10
!359)   O+    +   HCO   ->   OH+   +   CO
! 1738
xk(359)=4.30d-10*(T_K/300.d0)**(-0.50d0)
!360)   O+    +   HCO   ->   HCO+  +   O
! 3245
xk(360)=4.30d-10*(T_K/300.d0)**(-0.50d0)
!361)   O+    +   H2CO  ->   HCO+  +   OH
! 1741
xk(361)=1.40d-9*(T_K/300.d0)**(-0.50d0)
!362)   O+    +   H2CO  ->   H2CO+ +   O
! 3246
xk(362)=2.10d-9*(T_K/300.d0)**(-0.50d0)
!363)   O+    +   O2    ->   O2+   +   O
! 3247
xk(363)=1.90d-11
!364)   O+    +   CO2   ->   O2+   +   CO
! 1768
xk(364)=9.40d-10
!365)   CH4+  +   H     ->   CH3+  +   H2
! 556
xk(365)=1.00d-11
!366)   CH4+  +   O     ->   CH3+  +   OH
! 1713
xk(366)=1.00d-9
!367)   CH4+  +   H2    ->   CH5+  +   H
! 660, 661
if(T_K > 300.d0) then
   xk(367)=3.30d-11
else
   xk(367)=3.40d-11*(T_K/300.d0)**(-1.35d0)*dexp(-23.d0/T_K)
endif
!368)   CH4+  +   CH4   ->   CH5+  +   CH3
! 1840
xk(368)=1.50d-9
!369)   CH4+  +   H2O   ->   H3O+  +   CH3
! 1845
xk(369)=2.60d-9*(T_K/300.d0)**(-0.50d0)
!370)   CH4+  +   CO    ->   HCO+  +   CH3
! 1863
xk(370)=1.40d-9
!371)   CH4+  +   H2CO  ->   H2CO+ +   CH4
! 3272
xk(371)=1.62d-9*(T_K/300.d0)**(-0.50d0)
!372)   CH4+  +   H2CO  ->   H3CO+ +   CH3
! 1869
xk(372)=1.98d-9*(T_K/300.d0)**(-0.50d0)
!373)   CH4+  +   O2    ->   O2+   +   CH4
! 3273
xk(373)=3.90d-10
!374)   CH4+  +   CO2   ->   HCO2+ +   CH3
! 1889
xk(374)=1.20d-9
!375)   OH+   +   C     ->   CH+   +   O
! 1185
xk(375)=1.20d-9
!376)   OH+   +   O     ->   O2+   +   H
! 1715
xk(376)=7.10d-10
!377)   OH+   +   H2    ->   H2O+  +   H
! 664
xk(377)=1.01d-9
!378)   OH+   +   CH    ->   CH+   +   OH
! 3164
xk(378)=3.50d-10*(T_K/300.d0)**(-0.50d0)
!379)   OH+   +   CH    ->   CH2+  +   O
! 1411
xk(379)=3.50d-10*(T_K/300.d0)**(-0.50d0)
!380)   OH+   +   CH2   ->   CH2+  +   OH
! 3205
xk(380)=4.80d-10
!381)   OH+   +   CH2   ->   CH3+  +   O
! 1558
xk(381)=4.80d-10
!382)   OH+   +   CH4   ->   H3O+  +   CH2
! 1842
xk(382)=1.31d-9
!383)   OH+   +   CH4   ->   CH5+  +   O
! 1841
xk(383)=1.95d-10
!384)   OH+   +   OH    ->   H2O+  +   O
! 1923
xk(384)=7.00d-10*(T_K/300.d0)**(-0.50d0)
!385)   OH+   +   H2O   ->   H2O+  +   OH
! 3279
xk(385)=1.59d-9*(T_K/300.d0)**(-0.50d0)
!386)   OH+   +   H2O   ->   H3O+  +   O
! 1927
xk(386)=1.30d-9*(T_K/300.d0)**(-0.50d0)
!387)   OH+   +   C2    ->   C2+   +   OH
! 3280
xk(387)=4.80d-10
!388)   OH+   +   CO    ->   HCO+  +   O
! 1937
xk(388)=1.05d-9
!389)   OH+   +   HCO   ->   H2O+  +   CO
! 1943
xk(389)=2.80d-10*(T_K/300.d0)**(-0.50d0)
!390)   OH+   +   HCO   ->   HCO+  +   OH
! 3287
xk(390)=2.80d-10*(T_K/300.d0)**(-0.50d0)
!391)   OH+   +   HCO   ->   H2CO+ +   O
! 1944
xk(391)=2.80d-10*(T_K/300.d0)**(-0.50d0)
!392)   OH+   +   H2CO  ->   H3CO+ +   O
! 1949
xk(392)=1.12d-9*(T_K/300.d0)**(-0.50d0)
!393)   OH+   +   H2CO  ->   H2CO+ +   OH
! 3289
xk(393)=7.44d-10*(T_K/300.d0)**(-0.50d0)
!394)   OH+   +   O2    ->   O2+   +   OH
! 3291
xk(394)=5.90d-10
!395)   OH+   +   CO2   ->   HCO2+ +   O
! 1962
xk(395)=1.44d-9
!396)   CH5+  +   H     ->   CH4+  +   H2
! 557
xk(396)=1.50d-10
!397)   CH5+  +   C     ->   CH+   +   CH4
! 1189
xk(397)=1.20d-9
!398)   CH5+  +   O     ->   H3O+  +   CH2
! 1717
xk(398)=2.20d-10
!399)   CH5+  +   O     ->   H3CO+ +   H2
! 1718
xk(399)=4.40d-12
!400)   CH5+  +   CH    ->   CH2+  +   CH4
! 1416
xk(400)=6.90d-10*(T_K/300.d0)**(-0.50d0)
!401)   CH5+  +   CH2   ->   CH3+  +   CH4
! 1562
xk(401)=9.60d-10
!402)   CH5+  +   OH    ->   H2O+  +   CH4
! 1926
xk(402)=7.00d-10*(T_K/300.d0)**(-0.50d0)
!403)   CH5+  +   H2O   ->   H3O+  +   CH4
! 2021
xk(403)=3.70d-9*(T_K/300.d0)**(-0.50d0)
!404)   CH5+  +   CO    ->   HCO+  +   CH4
! 2028
xk(404)=1.00d-9
!405)   CH5+  +   HCO   ->   H2CO+ +   CH4
! 2032
xk(405)=8.50d-10*(T_K/300.d0)**(-0.50d0)
!406)   CH5+  +   H2CO  ->   H3CO+ +   CH4
! 2033
xk(406)=4.50d-9*(T_K/300.d0)**(-0.50d0)
!407)   CH5+  +   CO2   ->   HCO2+ +   CH4
! 2039
xk(407)=3.20d-11
!408)   H2O+  +   C    ->   CH+   +   OH
! 1190
xk(408)=1.10d-9
!409)   H2O+  +   O     ->   O2+   +   H2
! 1719
xk(409)=4.00d-11
!410)   H2O+  +   H2    ->   H3O+  +   H
! 669
xk(410)=6.40d-10
!411)   H2O+  +   CH    ->   CH+   +   H2O
! 3166
xk(411)=3.40d-10*(T_K/300.d0)**(-0.50d0)
!412)   H2O+  +   CH    ->   CH2+  +   OH
! 1418
xk(412)=3.40d-10*(T_K/300.d0)**(-0.50d0)
!413)   H2O+  +   CH2   ->   CH3+  +   OH
! 1564
xk(413)=4.70d-10
!414)   H2O+  +   CH2   ->   CH2+  +   H2O
! 3206
xk(414)=4.70d-10
!415)   H2O+  +   CH4   ->   H3O+  +   CH3
! 1846
xk(415)=1.40d-9
!416)   H2O+  +   OH    ->   H3O+  +   O
! 1928
xk(416)=6.90d-10*(T_K/300.d0)**(-0.50d0)
!417)   H2O+  +   H2O   ->   H3O+  +   OH
! 2043
xk(417)=2.10d-9*(T_K/300.d0)**(-0.50d0)
!418)   H2O+  +   C2    ->   C2+   +   H2O
! 3319
xk(418)=4.70d-10
!419)   H2O+  +   CO    ->   HCO+  +   OH
! 2058
xk(419)=5.00d-10
!420)   H2O+  +   HCO   ->   H3O+  +   CO
! 2062
xk(420)=2.80d-10*(T_K/300.d0)**(-0.50d0)
!421)   H2O+  +   HCO   ->   HCO+  +   H2O
! 3328
xk(421)=2.80d-10*(T_K/300.d0)**(-0.50d0)
!422)   H2O+  +   HCO   ->   H2CO+ +   OH
! 2063
xk(422)=2.80d-10*(T_K/300.d0)**(-0.50d0)
!423)   H2O+  +   H2CO  ->   H2CO+ +   H2O
! 3330
xk(423)=1.41d-9*(T_K/300.d0)**(-0.50d0)
!424)   H2O+  +   H2CO  ->   H3CO+ +   OH
! 2068
xk(424)=6.62d-10*(T_K/300.d0)**(-0.50d0)
!425)   H2O+  +   O2    ->   O2+   +   H2O
! 3332
xk(425)=4.60d-10
!426)   H3O+  +   C    ->   HCO+  +   H2
! 1194
xk(426)=1.00d-11
!427)   H3O+  +   H-    ->   OH    +   H2    +   H
! 4603
xk(427)=2.30d-7*(T_K/300.d0)**(-0.50d0)
!428)   H3O+  +   H-    ->   H2O   +   H2
! 4604
xk(428)=2.30d-7*(T_K/300.d0)**(-0.50d0)
!429)   H3O+  +   CH    ->   CH2+  +   H2O
! 1421
xk(429)=6.80d-10*(T_K/300.d0)**(-0.50d0)
!430)   H3O+  +   CH2   ->   CH3+  +   H2O
! 1565
xk(430)=9.40d-10
!431)   H3O+  +   H2CO  ->   H3CO+ +   H2O
! 2122
xk(431)=3.40d-9*(T_K/300.d0)**(-0.50d0)
!432)   C2+   +   C     ->   C+    +   C2
! 3104
xk(432)=1.10d-10
!433)   C2+   +   O     ->   CO+   +   C
! 1721
xk(433)=3.10d-10
!434)   C2+   +   CH    ->   CH+   +   C2
! 3168
xk(434)=3.20d-10*(T_K/300.d0)**(-0.50d0)
!435)   C2+   +   CH2   ->   CH2+  +   C2
! 3207
xk(435)=4.50d-10
!436)   C2+   +   OH    ->   OH+   +   C2
! 3281
xk(436)=6.50d-10*(T_K/300.d0)**(-0.50d0)
!437)   C2+   +   HCO   ->   HCO+  +   C2
! 3356
xk(437)=3.80d-10*(T_K/300.d0)**(-0.50d0)
!438)   C2+   +   O2    ->   CO+   +   CO
! 2187
xk(438)=8.00d-10
!439)   CO+   +   H     ->   H+    +   CO
! 2960
xk(439)=7.50d-10
!440)   CO+   +   C     ->   C+    +   CO
! 3107
xk(440)=1.10d-10
!441)   CO+   +   O     ->   O+    +   CO
! 3241
xk(441)=1.40d-10
!442)   CO+   +   H2    ->   HCO+  +   H
! 683
xk(442)=7.50d-10
!443)   CO+   +   CH    ->   CH+   +   CO
! 3171
xk(443)=3.20d-10*(T_K/300.d0)**(-0.50d0)
!444)   CO+   +   CH    ->   HCO+  +   C
! 1436
xk(444)=3.20d-10*(T_K/300.d0)**(-0.50d0)
!445)   CO+   +   CH2   ->   CH2+  +   CO
! 3209
xk(445)=4.30d-10
!446)   CO+   +   CH2   ->   HCO+  +   CH
! 1573
xk(446)=4.30d-10
!447)   CO+   +   CH4   ->   CH4+  +   CO
! 3270
xk(447)=7.93d-10
!448)   CO+   +   CH4   ->   HCO+  +   CH3
! 1864
xk(448)=4.55d-10
!449)   CO+   +   OH    ->   OH+   +   CO
! 3285
xk(449)=3.10d-10*(T_K/300.d0)**(-0.50d0)
!450)   CO+   +   OH    ->   HCO+  +   O
! 1938
xk(450)=3.10d-10*(T_K/300.d0)**(-0.50d0)
!451)   CO+   +   H2O   ->   H2O+  +   CO
! 3324
xk(451)=1.72d-9*(T_K/300.d0)**(-0.50d0)
!452)   CO+   +   H2O   ->   HCO+  +   OH
! 2059
xk(452)=8.84d-10*(T_K/300.d0)**(-0.50d0)
!453)   CO+   +   C2    ->   C2+   +   CO
! 3354
xk(453)=8.40d-10
!454)   CO+   +   HCO   ->   HCO+  +   CO
! 3410
xk(454)=7.40d-10*(T_K/300.d0)**(-0.50d0)
!455)   CO+   +   H2CO  ->   HCO+  +   HCO
! 2450
xk(455)=1.65d-9*(T_K/300.d0)**(-0.50d0)
!456)   CO+   +   H2CO  ->   H2CO+ +   CO
! 3412
xk(456)=1.35d-9*(T_K/300.d0)**(-0.50d0)
!457)   CO+   +   O2    ->   O2+   +   CO
! 3413
xk(457)=1.20d-10
!458)   HCO+  +   C     ->   CH+   +   CO
! 1213
xk(458)=1.10d-9
!459)   HCO+  +   H-    ->   CO    +   H2
! 4605
xk(459)=2.30d-7*(T_K/300.d0)**(-0.50d0)
!460)   HCO+  +   CH    ->   CH2+  +   CO
! 1441
xk(460)=6.30d-10*(T_K/300.d0)**(-0.50d0)
!461)   HCO+  +   CH2   ->   CH3+  +   CO
! 1578
xk(461)=8.60d-10
!462)   HCO+  +   OH    ->   H2O+  +   CO
! 1945
xk(462)=6.20d-10*(T_K/300.d0)**(-0.50d0)
!463)   HCO+  +   OH    ->   HCO2+ +   H
! 1942
xk(463)=1.00d-9*(T_K/300.d0)**(-0.50d0)
!464)   HCO+  +   H2O   ->   H3O+  +   CO
! 2064
xk(464)=2.50d-9*(T_K/300.d0)**(-0.50d0)
!465)   HCO+  +   HCO   ->   H2CO+ +   CO
! 2587
xk(465)=7.30d-10*(T_K/300.d0)**(-0.50d0)
!466)   HCO+  +   H2CO  ->   H3CO+ +   CO
! 2592
xk(466)=3.30d-9*(T_K/300.d0)**(-0.50d0)
!467)   H2CO+ +   CH    ->   CH+   +   H2CO
! 3176
xk(467)=3.10d-10*(T_K/300.d0)**(-0.50d0)
!468)   H2CO+ +   CH    ->   CH2+  +   HCO
! 1448
xk(468)=3.10d-10*(T_K/300.d0)**(-0.50d0)
!469)   H2CO+ +   CH2   ->   CH3+  +   HCO
! 1582
xk(469)=4.30d-10
!470)   H2CO+ +   CH2   ->   CH2+  +   H2CO
! 3212
xk(470)=4.30d-10
!471)   H2CO+ +   CH4   ->   H3CO+ +   CH3
! 1870
xk(471)=9.35d-11
!472)   H2CO+ +   H2O   ->   H3O+  +   HCO
! 2069
xk(472)=2.60d-9*(T_K/300.d0)**(-0.50d0)
!473)   H2CO+ +   HCO   ->   HCO+  +   H2CO
! 3442
xk(473)=3.60d-10*(T_K/300.d0)**(-0.50d0)
!474)   H2CO+ +   HCO   ->   H3CO+ +   CO
! 2593
xk(474)=3.60d-10*(T_K/300.d0)**(-0.50d0)
!475)   H2CO+ +   H2CO  ->   H3CO+ +   HCO
! 2712
xk(475)=3.20d-9*(T_K/300.d0)**(-0.50d0)
!476)   H2CO+ +   O2    ->   HCO+  +   O2H
! 2715
xk(476)=7.70d-11
!477)   H3CO+ +   CH    ->   CH2+  +   H2CO
! 1450
xk(477)=6.20d-10*(T_K/300.d0)**(-0.50d0)
!478)   H3CO+ +   H2O   ->   H3O+  +   H2CO
! 2076
xk(478)=2.30d-10*(T_K/300.d0)**(-0.50d0)
!479)   O2+   +   C    ->   CO+   +   O
! 1237
xk(479)=5.20d-11
!480)   O2+   +   C     ->   C+    +   O2
! 3119
xk(480)=5.20d-11
!481)   O2+   +   CH    ->   CH+   +   O2
! 3177
xk(481)=3.10d-10*(T_K/300.d0)**(-0.50d0)
!482)   O2+   +   CH    ->   HCO+  +   O
! 1454
xk(482)=3.10d-10*(T_K/300.d0)**(-0.50d0)
!483)   O2+   +   CH2   ->   H2CO+ +   O
! 1585
xk(483)=4.30d-10
!484)   O2+   +   CH2   ->   CH2+  +   O2
! 3213
xk(484)=4.30d-10
!485)   O2+   +   C2    ->   CO+   +   CO
! 2188
xk(485)=4.10d-10
!486)   O2+   +   C2    ->   C2+   +   O2
! 3358
xk(486)=4.10d-10
!487)   O2+   +   HCO   ->   O2H+  +   CO
! 2599
xk(487)=3.60d-10*(T_K/300.d0)**(-0.50d0)
!488)   O2+   +   HCO   ->   HCO+  +   O2
! 3443
xk(488)=3.60d-10*(T_K/300.d0)**(-0.50d0)
!489)   O2+   +   H2CO  ->   HCO+  +   O2    +   H
! 2714
xk(489)=2.30d-10*(T_K/300.d0)**(-0.50d0)
!490)   O2+   +   H2CO  ->   H2CO+ +   O2
! 3461
xk(490)=2.07d-9*(T_K/300.d0)**(-0.50d0)
!491)   O2H+  +   C     ->   CH+   +   O2
! 1243
xk(491)=1.00d-9
!492)   O2H+  +   O     ->   OH+   +   O2
! 1749
xk(492)=6.20d-10
!493)   O2H+  +   H2    ->   H3+   +   O2
! 697
xk(493)=6.40d-10
!494)   O2H+  +   CH    ->   CH2+  +   O2
! 1461
xk(494)=6.20d-10*(T_K/300.d0)**(-0.50d0)
!495)   O2H+  +   CH2   ->   CH3+  +   O2
! 1589
xk(495)=8.50d-10
!496)   O2H+  +   OH    ->   H2O+  +   O2
! 1959
xk(496)=6.10d-10*(T_K/300.d0)**(-0.50d0)
!497)   O2H+  +   H2O   ->   H3O+  +   O2
! 2085
xk(497)=8.20d-10*(T_K/300.d0)**(-0.50d0)
!498)   O2H+  +   CO    ->   HCO+  +   O2
! 2453
xk(498)=8.40d-10
!499)   O2H+  +   HCO   ->   H2CO+ +   O2
! 2605
xk(499)=7.10d-10*(T_K/300.d0)**(-0.50d0)
!500)   O2H+  +   H2CO  ->   H3CO+ +   O2
! 2720
xk(500)=9.80d-10*(T_K/300.d0)**(-0.50d0)
!501)   O2H+  +   CO2   ->   HCO2+ +   O2
! 2788
xk(501)=1.10d-9
!502)   HCO2+ +   !    ->   CH+   +   CO2
! 1296
xk(502)=1.00d-9
!503)   HCO2+ +   O     ->   HCO+  +   O2
! 1773
xk(503)=1.00d-9
!504)   HCO2+ +   CH4   ->   CH5+  +   CO2
! 1892
xk(504)=7.80d-10
!505)   HCO2+ +   H2O   ->   H3O+  +   CO2
! 2102
xk(505)=2.30d-9*(T_K/300.d0)**(-0.50d0)
!506)   HCO2+ +   CO    ->   HCO+  +   CO2
! 2459
xk(506)=7.80d-10
!507)   C+    +   e     ->   C     +   ph.
! 4012, 4013, 4014
if(T_K < 7950.d0) then
   xk(507)=4.67d-12*(T_K/300.d0)**(-0.6d0)
elseif(T_K < 21140.d0) then
   xk(507)=1.23d-17*(T_K/300.d0)**2.49d0*dexp(21845.6/T_K)
else
   xk(507)=9.62d-8*(T_K/300.d0)**(-1.37d0)*dexp(-115786.2/T_K)
endif
!508)   CH+   +   e     ->   C    +   H
! 3524
xk(508)=1.50d-7*(T_K/300.d0)**(-0.42d0)
!509)   CH2+  +   e     ->   C    +   H2
! 3527
xk(509)=7.68d-8*(T_K/300.d0)**(-0.60d0)
!510)   CH2+  +   e     ->   CH    +   H
! 3525
xk(510)=1.60d-7*(T_K/300.d0)**(-0.60d0)
!511)   CH3+  +   e     ->   CH2   +   H
! 3529
xk(511)=7.75d-8*(T_K/300.d0)**(-0.50d0)
!512)   CH3+  +   e     ->   CH    +   H2
! 3530
xk(512)=1.95d-7*(T_K/300.d0)**(-0.50d0)
!513)   CH3+  +   e     ->   CH    + 2 H
! 3531
xk(513)=2.00d-7*(T_K/300.d0)**(-0.40d0)
!514)   CH3+  +   e     ->   CH3   +   ph.
! 4017
xk(514)=1.10d-10*(T_K/300.d0)**(-0.50d0)
!515)   O+    +   e     ->   O     +   ph.
! 4019
xk(515)=3.24d-12*(T_K/300.d0)**(-0.66d0)
!516)   CH4+  +   e     ->   CH3   +   H
! 3538
xk(516)=1.75d-7*(T_K/300.d0)**(-0.50d0)
!517)   CH4+  +   e     ->   CH2   + 2 H
! 3539
xk(517)=1.75d-7*(T_K/300.d0)**(-0.50d0)
!518)   OH+   +   e     ->   O     +   H
! 3540
xk(518)=3.75d-8*(T_K/300.d0)**(-0.50d0)
!519)   CH5+  +   e     ->   CH4   +   H
! 3544
xk(519)=1.40d-8*(T_K/300.d0)**(-0.52d0)
!520)   CH5+  +   e     ->   CH3   +   H2
! 3543
xk(520)=1.40d-8*(T_K/300.d0)**(-0.52d0)
!521)   H2O+  +   e     ->   OH    +   H
! 3550
xk(521)=8.60d-8*(T_K/300.d0)**(-0.50d0)
!522)   H2O+  +   e     ->   O     +   H2
! 3549
xk(522)=3.90d-8*(T_K/300.d0)**(-0.50d0)
!523)   H3O+  +   e     ->   H2O   +   H
! 3563
xk(523)=1.08d-7*(T_K/300.d0)**(-0.50d0)
!524)   H3O+  +   e     ->   OH    + 2 H
! 3561
xk(524)=2.58d-7*(T_K/300.d0)**(-0.50d0)
!525)   C2+   +   e     -> 2 C
! 3567
xk(525)=3.00d-7*(T_K/300.d0)**(-0.50d0)
!526)   CO+   +   e     ->   O     +   C
! 3587
xk(526)=2.00d-7*(T_K/300.d0)**(-0.48d0)
!527)   HCO+  +   e     ->   CO    +   H
! 3601
xk(527)=2.40d-7*(T_K/300.d0)**(-0.69d0)
!528)   H2CO+ +   e     ->   HCO   +   H
! 3613
xk(528)=1.00d-7*(T_K/300.d0)**(-0.50d0)
!529)   H2CO+ +   e     ->   CO    + 2 H
! 3612
xk(529)=5.00d-7*(T_K/300.d0)**(-0.50d0)
!530)   H2CO+ +   e     ->   H2CO  +   ph.
! 4023
xk(530)=1.10d-10*(T_K/300.d0)**(-0.70d0)
!531)   H3CO+ +   e     ->   CO    +   H     +   H2
! 3624
xk(531)=2.00d-7*(T_K/300.d0)**(-0.50d0)
!532)   H3CO+ +   e     ->   HCO   + 2 H
! 3625
xk(532)=2.00d-7*(T_K/300.d0)**(-0.50d0)
!533)   H3CO+ +   e     ->   H2CO  +   H
! 3626
xk(533)=2.00d-7*(T_K/300.d0)**(-0.50d0)
!534)   O2+   +   e     -> 2 O
! 3633
xk(534)=1.95d-7*(T_K/300.d0)**(-0.70d0)
!535)   O2H+  +   e     ->   O2    +   H
! 3641
xk(535)=3.00d-7*(T_K/300.d0)**(-0.50d0)
!536)   HCO2+ +   e     ->   CO2   +   H
! 3731
xk(536)=6.00d-8*(T_K/300.d0)**(-0.64d0)
!537)   HCO2+ +   e     ->   CO    +   OH
! 3730
xk(537)=3.20d-7*(T_K/300.d0)**(-0.64d0)
!538)   H2    +   !    ->   CH2   +   ph.
! 4088
xk(538)=1.00d-17
!539)   OH    +   CH3   ->   CH4   +   O
! 289
xk(539)=3.27d-14*(T_K/300.d0)**2.20d0*dexp(-2240.d0/T_K)
!540)   H+    +   H2CO  ->   CO+   +   H2   +   H
! 574
xk(540)=1.06d-9*(T_K/300.d0)**(-0.5d0)
!541)   He+   +   !    ->   C+    +   He
! 3082
xk(541)=6.30d-15*(T_K/300.d0)**0.75d0
!542)   He+   +   H2CO  ->   H2CO+ +   He
! 3093
xk(542)=9.69d-10*(T_K/300.d0)**(-0.5d0)
!543)   He+   +   H2CO  ->   CH2+  +   O    +   He
! 964
xk(543)=1.71d-9*(T_K/300.d0)**(-0.5d0)
!544)   H     +   CR    ->   H+    +   e
! 4385
!     xk(544)=4.60d-1*zeta
xk(544)=5.98d-18*(zeta/1.36d-17)
!545)   He    +   CR    ->   He+   +   e
! 4390
!     xk(545)=5.00d-1*zeta
xk(545)=6.50d-18*(zeta/1.36d-17)
!546)   C    +   CR    ->   C+    +   e
! 4391
!     xk(546)=1.77d0*zeta
xk(546)=2.30d-17*(zeta/1.36d-17)
!547)   O     +   CR    ->   O+    +   e
! 4393
!     xk(547)=2.62d0*zeta
xk(547)=3.40d-17*(zeta/1.36d-17)
!548)   H2    +   CR    ->   H+    +   H     +   e
! 4386
!     xk(548)=1.69d-2*zeta
xk(548)=2.20d-19*(zeta/1.36d-17)
!549)   H2    +   CR    ->   H2+   +   e
! 4389
!     xk(549)=9.23d-1*zeta
xk(549)=1.20d-17*(zeta/1.36d-17)
!550)   H2    +   CR    -> 2 H
! 4387
!     xk(550)=1.00d-1*zeta
 xk(550)=1.30d-18*(zeta/1.36d-17)
!551)   H2    +   CR    ->   H+    +   H-
! 4388
!     xk(551)=3.00d-4*zeta
 xk(551)=3.90d-21*(zeta/1.36d-17)
!552)   CO    +   CR    ->   CO+   +   e
! 4394
!     xk(552)=3.00d0*zeta
 xk(552)=3.90d-17*(zeta/1.36d-17)
!553)   C     +   ph.   ->   C+    +   e
! 4173
xk(553)=3.00d-10*dexp(-3.d0*A_v)*(G_0/1.71d0)
!554)   H-    +   ph.   ->   H     +   e
! 4169
xk(554)=2.40d-7*dexp(-0.5d0*A_v)*(G_0/1.71d0)
!555)   H2+   +   ph.   ->   H+    +   H
! 4170
xk(555)=5.70d-10*dexp(-1.9d0*A_v)*(G_0/1.71d0)
!556)   H3+   +   ph.   ->   H2+   +   H
! 4171
xk(556)=5.00d-15*dexp(-2.3d0*A_v)*(G_0/1.71d0)
!557)   H3+   +   ph.   ->   H+    +   H2
! 4172
xk(557)=5.00d-15*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!558)   CH    +   ph.   ->   CH+   +   e
! 4175
xk(558)=7.60d-10*dexp(-2.8d0*A_v)*(G_0/1.71d0)
!559)   CH    +   ph.   ->   C    +   H
! 4176
xk(559)=8.60d-10*dexp(-1.2d0*A_v)*(G_0/1.71d0)
!560)   CH+   +   ph.   ->   C+    +   H
! 4177
xk(560)=2.50d-10*dexp(-2.5d0*A_v)*(G_0/1.71d0)
!561)   CH2   +   ph.   ->   CH2+  +   e
! 4178
xk(561)=1.00d-9*dexp(-2.3d0*A_v)*(G_0/1.71d0)
!562)   CH2   +   ph.   ->   CH    +   H
! 4179
xk(562)=7.20d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!563)   CH2+  +   ph.   ->   CH+   +   H
! 4180
xk(563)=4.60d-11*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!564)   CH3   +   ph.   ->   CH    +   H2
! 4188
xk(564)=2.50d-10*dexp(-1.9d0*A_v)*(G_0/1.71d0)
!565)   CH3   +   ph.   ->   CH2   +   H
! 4186
xk(565)=2.50d-10*dexp(-1.9d0*A_v)*(G_0/1.71d0)
!566)   CH3   +   ph.   ->   CH3+  +   e
! 4185
xk(566)=1.00d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
!567)   CH3+  +   ph.   ->   CH2+  +   H
! 4187
xk(567)=1.00d-9*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!568)   CH3+  +   ph.   ->   CH+   +   H2
! 4189
xk(568)=1.00d-9*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!569)   CH4   +   ph.   ->   CH3   +   H
! 4193
xk(569)=2.20d-10*dexp(-2.2d0*A_v)*(G_0/1.71d0)
!570)   CH4   +   ph.   ->   CH2   +   H2
! 4194
xk(570)=9.80d-10*dexp(-2.2d0*A_v)*(G_0/1.71d0)
!571)   CH4   +   ph.   ->   CH    +   H     +   H2
! 4195
xk(571)=2.20d-10*dexp(-2.2d0*A_v)*(G_0/1.71d0)
!572)   OH    +   ph.   ->   O     +   H
! 4200
xk(572)=3.50d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!573)   OH    +   ph.   ->   OH+   +   e
! 4198
xk(573)=1.60d-12*dexp(-3.1d0*A_v)*(G_0/1.71d0)
!574)   OH+   +   ph.   ->   H+    +   O
! none
xk(574)=0.d0
!575)   H2O   +   ph.   ->   OH    +   H
! 4205
xk(575)=5.90d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!576)   H2O   +   ph.   ->   H2O+  +   e
! 4206
xk(576)=3.30d-11*dexp(-3.9d0*A_v)*(G_0/1.71d0)
!577)   C2    +   ph.   -> 2 C
! 4210
xk(577)=1.50d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
!578)   C2    +   ph.   ->   C2+   +   e
! 4211
xk(578)=4.10d-10*dexp(-3.5d0*A_v)*(G_0/1.71d0)
!579)   C2+   +   ph.   ->   C+    +   C
! 4212
xk(579)=1.00d-11*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!580)   CO    +   ph.   ->   C     +   O
! 4225
shield=3.5d0*A_v+xNc_H2/1.6d21
xk(580)=2.00d-10*dexp(-shield) *(G_0/1.71d0)
!581)   HCO   +   ph.   ->   H     +   CO
! 4231
xk(581)=1.10d-9*dexp(-0.8d0*A_v)*(G_0/1.71d0)
!582)   HCO   +   ph.   ->   HCO+  +   e
! 4232
xk(582)=5.60d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
!583)   H2CO  +   ph.   ->   CO    +   H2
! 4243
xk(583)=7.00d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!584)   H2CO  +   ph.   ->   CO    + 2 H
! 4242
xk(584)=7.00d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
!585)   H2CO  +   ph.   ->   H2CO+ +   e
! 4241
xk(585)=4.70d-10*dexp(-2.8d0*A_v)*(G_0/1.71d0)
!586)   H2CO  +   ph.   ->   HCO+  +   e     +   H
! 4244
xk(586)=1.40d-11*dexp(-3.1d0*A_v)*(G_0/1.71d0)
!587)   O2    +   ph.   -> 2 O
! 4247
xk(587)=6.90d-10*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!588)   O2    +   ph.   ->   O2+   +   e
! 4248
xk(588)=5.60d-11*dexp(-3.7d0*A_v)*(G_0/1.71d0)
!589)   CO2   +   ph.   ->   CO    +   O
! 4284
xk(589)=1.40d-9*dexp(-2.5d0*A_v)*(G_0/1.71d0)
!590)   H2   +   ph.   ->   2H
! none
xm_p=1.67d-24
xk_B=1.38d-16
D5=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))*1.d-5
x=xNc_H2/8.465d13
fsh=0.9379d0/(1.d0+x/D5)**1.879d0&
+(0.03465d0/(1.d0+x)**0.473d0)&
*dexp(-2.293d-4*(1.d0+x)**0.5d0)
!     if(xNc_H2.le.1.d14) then
!        fsh=1.d0
!     else
!        fsh=(1.d14/xNc_H2)**0.75d0
!     endif
D_0=G_0/1.71d0
xk(590)=7.70d-11*fsh*dexp(-2.5d0*A_v)*D_0

!591)   HD   +   ph.   ->    H     +   D
! none
D5=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p))*1.d-5
x=xNc_HD/8.465d13
fsh=0.9379d0/(1.d0+x/D5)**1.879d0&
+(0.03465d0/(1.d0+x)**0.473d0)&
*dexp(-2.293d-4*(1.d0+x)**0.5d0)
x=xNc_H/2.848d23
fH=dexp(-0.149d0*x)/(1.d0+x)**1.620d0
x=xNc_H/2.848d23
fH2=dexp(-5.2d-3*x)/(1.d0+x)**0.238d0
fsh=fsh*fH*fH2
xk(591)=7.70d-11*fsh*dexp(-2.5d0*A_v)*D_0
!supplemented reactions
!601)   H    +   HCO   ->   O   +   CH2
! 17
xk(601)=6.61d-11*dexp(-51598.d0/T_K)
!602)   H    +   H2O2  ->  H2O   +   OH
! 30
xk(602)=1.70d-11*dexp(-1800.d0/T_K)
!603)   C    +   CH2   ->   2 CH
! 64
xk(603)=2.69d-12*dexp(-23550.d0/T_K)
!604)   CH   +    O2   ->   CO   +   OH
! 154
xk(604)=2.60d-11
!605)   H    +   HCO   ->   O   +   CH2
! 17
xk(605)=8.00d-11
!606)   CH2  +   O2   ->   CO2    +   2 H
! 256
xk(606)=3.65d-11*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
!607)   CH2  +   O2   ->   CO2    +   H2
! 257
xk(607)=2.92d-11*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
!608)   CH2  +   O2   ->   CO     +   H2O
! 260
xk(608)=2.48d-10*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
!609)   2 CH3         ->   CH4    +   CH2
! 283
xk(609)=7.13d-12*dexp(-5052.d0/T_K)
!610)   CH3  +   O    ->   CO     +   H2    +   H
! 287
xk(610)=3.60d-11*dexp(-202.d0/T_K)
!611)   CH3  +   OH   ->   H2CO   +   H2
! 291
xk(611)=1.70d-12
!612)   CH3  +   O2   ->   HCO    +   H2O
! 301
xk(612)=1.66d-12
!613)   O    +  H2CO  ->   CO     +   OH  +   H
! 338
xk(613)=1.00d-10
!614)   C2   +   O2   ->   2 CO
! 464
xk(614)=1.50d-11*dexp(-4300.d0/T_K)
!615)   2 HCO         ->   2 CO   +   H2
! 515
xk(615)=3.63d-11
!616)   HCO   +   O2   ->   CO2   +   OH
! 519
xk(616)=7.60d-13
!617)    O    +   OH   ->   O2    +   H
! 544
xk(617)=3.50d-11
xk(618)=0.d0
xk(619)=0.d0
xk(620)=0.d0
!621)    H3+  +    O    ->  H2O+  +   H
! 755
xk(621)=3.60d-10
xk(622)=0.d0
xk(623)=0.d0
xk(624)=0.d0
xk(625)=0.d0
!626)   O+    +  CO     ->  CO+   +   O
! 3242
xk(626)=4.90d-12*(T_K/300.d0)**0.5d0*dexp(-4580.d0/T_K)
xk(627)=0.d0
xk(628)=0.d0
xk(629)=0.d0
xk(630)=0.d0
xk(631)=0.d0
!632)   CH2+  +  e      ->  C     +  2 H
! 3526
xk(632)=4.03d-7*(T_K/300.d0)**(-0.6d0)
!633)   CH5+  +  e      ->  CH3   +  2 H
! 3545
xk(633)=1.96d-7*(T_K/300.d0)**(-0.52d0)
!634)   CH5+  +  e      ->  CH2   +    H2    +   H
! 3546
xk(634)=4.76d-8*(T_K/300.d0)**(-0.52d0)
!635)   CH5+  +  e      ->  CH    +  2 H2
! 3547
xk(635)=8.40d-9*(T_K/300.d0)**(-0.52d0)
!636)   H2O+  +  e      ->  O     +  2 H
! 3548
xk(636)=3.05d-7*(T_K/300.d0)**(-0.5d0)
!637)   H3O+  +  e      ->  O     +    H2   +    H
! 3560
xk(637)=5.60d-9*(T_K/300.d0)**(-0.5d0)
!638)   H3O+  +  e      ->  OH    +    H2
! 3562
xk(638)=6.02d-8*(T_K/300.d0)**(-0.5d0)
xk(639)=0.d0
!640)   HCO2+ +  e      ->  CO    +    O   +    H
! 3732
xk(640)=8.10d-7*(T_K/300.d0)**(-0.64d0)
!641)   H+    +  He     ->  HeH+  +    ph.
! 4080
xk(641)=5.26d-20*(T_K/300.d0)**(-0.51d0)
!642)   H     +  OH     ->  H2O   +    ph.
! 4084
xk(642)=5.26d-18*(T_K/300.d0)**(-5.22d0)*dexp(-90.d0/T_K)
!643)   H2    +  CH     ->  CH3   +    ph.
! 4090
xk(643)=5.09d-18*(T_K/300.d0)**(-0.71d0)*dexp(-11.6d0/T_K)
!644)   C+    +  C     ->  C2+   +    ph.
! 4107
xk(644)=4.01d-18*(T_K/300.d0)**0.17d0*dexp(-101.5d0/T_K)
!645)   C     +  O+     ->  CO+   +    ph.
! 4111
if(T_K > 2000.d0) then
   xk(645)=4.69d-11*(T_K/300.d0)**(-3.08d0)&
*dexp(2114.d0/T_K)
else
   xk(645)=3.91d-11
endif
!646)   CH2+  +  ph.    ->  CH    +    H+
! 4181
xk(646)=4.60d-11*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!647)   CH2+  +  ph.    ->  C+    +    H2
! 4182
xk(647)=4.60d-11*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!648)   CH4+  +  ph.    ->  CH2+  +    H2
! 4196
xk(648)=3.40d-11*dexp(-2.d0*A_v)*(G_0/1.71d0)
!649)   CH4+  +  ph.    ->  CH3+  +    H
! 4197
xk(649)=8.00d-12*dexp(-2.d0*A_v)*(G_0/1.71d0)
!650)   OH+   +  ph.    ->  O+    +    H
! 4201
xk(650)=1.00d-12*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!651)   H2O+   +  ph.   ->  OH+   +    H
! 4207
xk(651)=1.00d-12*dexp(-2.d0*A_v)*(G_0/1.71d0)
!652)   CO+   +  ph.   ->   C+    +    O
! 4226
xk(652)=1.30d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
!653)   HCO+  +  ph.   ->   CO+   +    H
! 4233
xk(653)=5.40d-12*dexp(-2.d0*A_v)*(G_0/1.71d0)
!654)   O2+   +  ph.   ->   O+    +    O
! 4249
xk(654)=2.40d-11*dexp(-2.d0*A_v)*(G_0/1.71d0)
!655)   H2O2  +  ph.   -> 2 OH
! 4258
xk(655)=8.30d-10*dexp(-1.8d0*A_v)*(G_0/1.71d0)
!albedo
omega=0.6d0
!656)   C     + CR ph. ->   C+    +    e
! 4396
xk(656)=zeta*255.d0/(1.d0-omega)
!657)   CH    + CR ph. ->   C     +    H
! 4397
xk(657)=zeta*365.d0/(1.d0-omega)
!658)   CH+   + CR ph. ->   C+    +    H
! 4398
xk(658)=zeta*88.d0/(1.d0-omega)
!659)   CH2   + CR ph. ->   CH2+  +    e
! 4399
xk(659)=zeta*250.d0/(1.d0-omega)
!660)   CH2   + CR ph. ->   CH    +    H
! 4400
xk(660)=zeta*250.d0/(1.d0-omega)
!661)   CH3   + CR ph. ->   CH3+  +    e
! 4403
xk(661)=zeta*250.d0/(1.d0-omega)
!662)   CH3   + CR ph. ->   CH2   +    H
! 4404
xk(662)=zeta*250.d0/(1.d0-omega)
!663)   CH3   + CR ph. ->   CH    +    H2
! 4405
xk(663)=zeta*250.d0/(1.d0-omega)
!664)   CH4   + CR ph. ->   CH2   +    H2
! 4408
xk(664)=zeta*1169.5d0/(1.d0-omega)
!665)   OH    + CR ph. ->   O     +    H
! 4409
xk(665)=zeta*254.5d0/(1.d0-omega)
!666)   H2O   + CR ph. ->   OH    +    H
! 4413
xk(666)=zeta*485.5d0/(1.d0-omega)
!667)   C2    + CR ph. -> 2 C
! 4416
xk(667)=zeta*119.5d0/(1.d0-omega)
!668)   CO    + CR ph. ->   O     +    C
! 4427
xk(668)=zeta*(T_K/300.d0)**1.17d0*105.d0/(1.d0-omega)
!669)  HCO    + CR ph. ->   CO    +    H
! 4432
xk(669)=zeta*210.5d0/(1.d0-omega)
!670)  HCO    + CR ph. ->   HCO+  +    e
! 4433
xk(670)=zeta*584.5d0/(1.d0-omega)
!671)  H2CO   + CR ph. ->   CO    +    H2
! 4440
xk(671)=zeta*1329.5d0/(1.d0-omega)
!672)  O2     + CR ph. -> 2 O
! 4445
xk(672)=zeta*375.5d0/(1.d0-omega)
!673)  O2     + CR ph. ->   O2+   +    e
! 4446
xk(673)=zeta*58.5d0/(1.d0-omega)
!674)  H2O2   + CR ph. -> 2 OH
! 4452
xk(674)=zeta*750.d0/(1.d0-omega)
!675)  CO2    + CR ph. ->   CO    +    O
! 4474
xk(675)=zeta*854.d0/(1.d0-omega)

return
END

SUBROUTINE react_rat(xk,xnH,y,r_f_tot)
IMPLICIT REAL*8(a-h,o-z)
!double precision :: xk,xnH,y,r_f_tot
!****************************************************************
!*     dy(i)/dt=r_f_tot(i)                                      *
!*     This subroutine returns reaction rate for each spieces,  *
!*     r_f_tot(i).                                              *
!****************************************************************
!    N_sp = number of spiecies
!    N_react = number of reactions
integer, PARAMETER :: N_sp=50,N_react=675
double precision ::  y(N_sp),r_f(N_react,N_sp),r_f_tot(N_sp),xk(N_react)
!****************************************************************
!*     SPECIES                                                  *
!*     1 : H      2 : H2     3 : e      4 : H+     5 : H2+      *
!*     6 : H3+    7 : H-                                        *
!*     8 : He     9 : He+    10: He++   11: HeH+                *
!*     12: D      13: HD     14: D+     15: HD+    16: D-      *
!*     17: C      18: C2     19: CH     20: CH2    21: CH3      *
!*     22: CH4    23: C+     24: C2+    25: CH+    26: CH2+     *
!*     27: CH3+   28: CH4+   29: CH5+                           *
!*     30: O      31: O2     32: OH     33: CO     34: H2O      *
!*     35: HCO    36: O2H    37: CO2    38: H2CO   39: H2O2     *
!*     40: O+     41: O2+    42: OH+    43: CO+    44: H2O+     *
!*     45: HCO+   46: O2H+   47: H3O+   48: H2CO+  49: HCO2+    *
!*     50: H3CO+                                                *
!****************************************************************

do isp=1,N_sp
   do ire=1,N_react
      r_f(ire,isp)=0.d0
   enddo
enddo

!********* primordial gas reactions *********

!  1)   H     +   e     ->   H+    + 2 e
!        (1         3          4       2*3)
rate=xk(1)*y(1)*y(3)*xnH
r_f(1,1)=-rate
r_f(1,3)=rate
r_f(1,4)=rate
!  2)   H+    +   e     ->   H     +   ph.
!        (4         3          1)
rate=xk(2)*y(4)*y(3)*xnH
r_f(2,4)=-rate
r_f(2,3)=-rate
r_f(2,1)=rate
!  3)   He    +   e     ->   He+   + 2 e
!        (8         3          9       2*3)
rate=xk(3)*y(8)*y(3)*xnH
r_f(3,8)=-rate
r_f(3,9)=rate
r_f(3,3)=rate
!  4)   He+   +   e     ->   He    +   ph.
!        (9         3          8)
rate=xk(4)*y(9)*y(3)*xnH
r_f(4,9)=-rate
r_f(4,3)=-rate
r_f(4,8)=rate
!  5)   He+   +   e     ->   He++  + 2 e
!        (9         3          10      2*3)
rate=xk(5)*y(9)*y(3)*xnH
r_f(5,9)=-rate
r_f(5,3)=rate
r_f(5,10)=rate
!  6)   He++  +   e     ->   He+   +   ph.
!        (10        3          9)
rate=xk(6)*y(10)*y(3)*xnH
r_f(6,10)=-rate
r_f(6,3)=-rate
r_f(6,9)=rate
!  7)   H     +   e     ->   H-    +   ph.
!        (1         3          7)
rate=xk(7)*y(1)*y(3)*xnH
r_f(7,1)=-rate
r_f(7,3)=-rate
r_f(7,7)=rate
!  8)   H-    +   H     ->   H2    +   e
!        (7         1          2         3)
rate=xk(8)*y(7)*y(1)*xnH
r_f(8,7)=-rate
r_f(8,1)=-rate
r_f(8,2)=rate
r_f(8,3)=rate
!  9)   H     +   H+    ->   H2+   +   ph.
!        (1         4          5)
rate=xk(9)*y(1)*y(4)*xnH
r_f(9,1)=-rate
r_f(9,4)=-rate
r_f(9,5)=rate
! 10)   H2+   +   H     ->   H2    +   H+
!        (5         1          2         4)
rate=xk(10)*y(5)*y(1)*xnH
r_f(10,5)=-rate
r_f(10,1)=-rate
r_f(10,2)=rate
r_f(10,4)=rate
! 11)   H2    +   H+    ->   H2+   +   H
!        (2         4          5         1)
rate=xk(11)*y(2)*y(4)*xnH
r_f(11,2)=-rate
r_f(11,4)=-rate
r_f(11,5)=rate
r_f(11,1)=rate
! 12)   H2    +   e     -> 2 H     +   e
!        (2         3        2*1         3)
rate=xk(12)*y(2)*y(3)*xnH
r_f(12,2)=-rate
r_f(12,1)=2.d0*rate
! 13)   H2    +   H     -> 3 H
!        (2         1        3*1)
rate=xk(13)*y(2)*y(1)*xnH
r_f(13,2)=-rate
r_f(13,1)=2.d0*rate
! 14)   H-    +   e     ->   H     + 2 e
!        (7         3          1       2*3)
rate=xk(14)*y(7)*y(3)*xnH
r_f(14,7)=-rate
r_f(14,3)=rate
r_f(14,1)=rate
! 15)   H-    +   H+    -> 2 H
!        (7         4        2*1)
rate=xk(15)*y(7)*y(4)*xnH
r_f(15,7)=-rate
r_f(15,4)=-rate
r_f(15,1)=2.d0*rate
! 16)   H-    +   H+    ->   H2+   +   e
!        (7         4          5         3)
rate=xk(16)*y(7)*y(4)*xnH
r_f(16,7)=-rate
r_f(16,4)=-rate
r_f(16,5)=rate
r_f(16,3)=rate
! 17)   H2+   +   e     -> 2 H
!        (5         3        2*1)
rate=xk(17)*y(5)*y(3)*xnH
r_f(17,5)=-rate
r_f(17,3)=-rate
r_f(17,1)=2.d0*rate
! 18)   H2+   +   H-    ->   H2    +   H
!        (5         7          2         1)
rate=xk(18)*y(5)*y(7)*xnH
r_f(18,5)=-rate
r_f(18,7)=-rate
r_f(18,2)=rate
r_f(18,1)=rate
! 19) 3 H               ->   H2    +   H
!      (3*1                    2         1)
rate=xk(19)*(y(1)**3)*(xnH**2)
r_f(19,1)=-2.d0*rate
r_f(19,2)=rate
! 20) 2 H     +   H2    -> 2 H2
!      (2*1         2        2*2)
rate=xk(20)*(y(1)**2)*y(2)*(xnH**2)
r_f(20,1)=-2.d0*rate
r_f(20,2)=rate
! 21) 2 H2              -> 2 H     +   H2
!      (2*2                  2*1         2)
rate=xk(21)*(y(2)**2)*xnH
r_f(21,2)=-rate
r_f(21,1)=2.d0*rate
! 22) 2 H               ->   H+    +   e     +   H
!      (2*1                    4         3         1)
rate=xk(22)*(y(1)**2)*xnH
r_f(22,1)=-rate
r_f(22,4)=rate
r_f(22,3)=rate
! 23) 2 H     +   grain ->   H2
!      (2*1                    2) special
rate=xk(23)*y(1)*xnH
r_f(23,1)=-2.d0*rate
r_f(23,2)=rate
! 24)   He+   +   H2    ->   H+    +   H     +   He
!        (9         2          4         1         8)
rate=xk(24)*y(9)*y(2)*xnH
r_f(24,9)=-rate
r_f(24,2)=-rate
r_f(24,4)=rate
r_f(24,1)=rate
r_f(24,8)=rate
! 25)   H2+   +   He    ->   HeH+  +   H
!        (5         8          11        1)
rate=xk(25)*y(5)*y(8)*xnH
r_f(25,5)=-rate
r_f(25,8)=-rate
r_f(25,11)=rate
r_f(25,1)=rate
! 26)   H2+   +   H2    ->   H3+   +   H
!        (5         2          6         1)
rate=xk(26)*y(5)*y(2)*xnH
r_f(26,5)=-rate
r_f(26,2)=-rate
r_f(26,6)=rate
r_f(26,1)=rate
! 27)   H3+   +   H-    -> 2 H2
!        (6         7        2*2)
rate=xk(27)*y(6)*y(7)*xnH
r_f(27,6)=-rate
r_f(27,7)=-rate
r_f(27,2)=2.d0*rate
! 28)   He+   +   H     ->   H+    +   He
!        (9         1          4         8)
rate=xk(28)*y(9)*y(1)*xnH
r_f(28,9)=-rate
r_f(28,1)=-rate
r_f(28,4)=rate
r_f(28,8)=rate
! 29)   He+   +   H-    ->   H     +   He
!        (9         7          1         8)
rate=xk(29)*y(9)*y(7)*xnH
r_f(29,9)=-rate
r_f(29,7)=-rate
r_f(29,1)=rate
r_f(29,8)=rate
! 30)   He+   +   H2    ->   H2+   +   He
!        (9         2          5         8)
rate=xk(30)*y(9)*y(2)*xnH
r_f(30,9)=-rate
r_f(30,2)=-rate
r_f(30,5)=rate
r_f(30,8)=rate
! 31)   HeH+  +   H     ->   H2+   +   He
!        (11        1          5         8)
rate=xk(31)*y(11)*y(1)*xnH
r_f(31,11)=-rate
r_f(31,1)=-rate
r_f(31,5)=rate
r_f(31,8)=rate
! 32)   HeH+  +   H2    ->   H3+   +   He
!        (11        2          6         8)
rate=xk(32)*y(11)*y(2)*xnH
r_f(32,11)=-rate
r_f(32,2)=-rate
r_f(32,6)=rate
r_f(32,8)=rate
! 33)   H2+   +   e     ->   H2    +   ph.
!        (5         3          2)
rate=xk(33)*y(5)*y(3)*xnH
r_f(33,5)=-rate
r_f(33,3)=-rate
r_f(33,2)=rate
! 34)   H3+   +   e     ->   H2    +   H
!        (6         3          2         1)
rate=xk(34)*y(6)*y(3)*xnH
r_f(34,6)=-rate
r_f(34,3)=-rate
r_f(34,2)=rate
r_f(34,1)=rate
! 35)   H3+   +   e     -> 3 H
!        (6         3        3*1)
rate=xk(35)*y(6)*y(3)*xnH
r_f(35,6)=-rate
r_f(35,3)=-rate
r_f(35,1)=3.d0*rate
! 36)   HeH+  +   e     ->   H     +   He
!        (11        3          1         8)
rate=xk(36)*y(11)*y(3)*xnH
r_f(36,11)=-rate
r_f(36,3)=-rate
r_f(36,1)=rate
r_f(36,8)=rate
! 37)   H2    +   He    -> 2 H     +   He
!      (2         8        2x1         8)
rate=xk(37)*y(2)*y(8)*xnH
r_f(37,2)=-rate
r_f(37,1)=2.d0*rate
! 38)   H-    +   H     -> 2 H     +   e
!      (7         1        2x1         3)
rate=xk(38)*y(7)*y(1)*xnH
r_f(38,7)=-rate
r_f(38,1)=rate
r_f(38,3)=rate
! 39)   H-    +   H2+   -> 3 H
!      (7         5        3x1     )
rate=xk(39)*y(7)*y(5)*xnH
r_f(39,7)=-rate
r_f(39,5)=-rate
r_f(39,1)=3.d0*rate
! 40)   H2    +   e     ->   H-    +   H
!      (2         3          7         1)
rate=xk(40)*y(2)*y(3)*xnH
r_f(40,2)=-rate
r_f(40,3)=-rate
r_f(40,7)=rate
r_f(40,1)=rate
! 41)   He    +   H+    ->   He+   +   H
!      (8         4          9         1)
rate=xk(41)*y(8)*y(4)*xnH
r_f(41,8)=-rate
r_f(41,4)=-rate
r_f(41,9)=rate
r_f(41,1)=rate
! 42)   He    +   H-    ->   He    +   H   +   e
!      (8         7          8         1       3)
rate=xk(42)*y(8)*y(7)*xnH
r_f(42,7)=-rate
r_f(42,1)=rate
r_f(42,3)=rate
! 43) 2 H     +   He    ->   H2    +   He
!    (2x1         8          2         8)
rate=xk(43)*(y(1)**2)*y(8)*(xnH**2)
r_f(43,1)=-2.d0*rate
r_f(43,2)=rate

!    D reactions
! 51)   D+  +   e    ->   D    +   ph.
!    (  14      3         12   )
rate=xk(51)*y(14)*y(3)*xnH
r_f(51,14)=-rate
r_f(51,3)=-rate
r_f(51,12)=rate
! 52)   D   +   H+   ->   D+   +   H
!    (  12      4         14       1 )
rate=xk(52)*y(12)*y(4)*xnH
r_f(52,12)=-rate
r_f(52,4)=-rate
r_f(52,14)=rate
r_f(52,1)=rate
! 53)   D+  +   H    ->   D    +   H+
!    (  14      1         12       4 )
rate=xk(53)*y(14)*y(1)*xnH
r_f(53,14)=-rate
r_f(53,1)=-rate
r_f(53,12)=rate
r_f(53,4)=rate
! 54)   D   +   H    ->   HD   +   ph.
!    (  12      1         13  )
rate=xk(54)*y(12)*y(1)*xnH
r_f(54,12)=-rate
r_f(54,1)=-rate
r_f(54,13)=rate
! 55)   D   +   H2   ->   H    +   HD
!    (  12      2         1        13  )
rate=xk(55)*y(12)*y(2)*xnH
r_f(55,12)=-rate
r_f(55,2)=-rate
r_f(55,1)=rate
r_f(55,13)=rate
! 56)   HD+ +   H    ->   H+   +   HD
!    (  15      1         4        13  )
rate=xk(56)*y(15)*y(1)*xnH
r_f(56,15)=-rate
r_f(56,1)=-rate
r_f(56,4)=rate
r_f(56,13)=rate
! 57)   D+  +   H2   ->   H+   +   HD
!    (  14      2         4        13  )
rate=xk(57)*y(14)*y(2)*xnH
r_f(57,14)=-rate
r_f(57,2)=-rate
r_f(57,4)=rate
r_f(57,13)=rate
! 58)   HD  +   H    ->   H2   +   D
!    (  13      1         2        12  )
rate=xk(58)*y(13)*y(1)*xnH
r_f(58,13)=-rate
r_f(58,1)=-rate
r_f(58,2)=rate
r_f(58,12)=rate
! 59)   HD  +   H+   ->   H2   +   D+
!    (  13      4         2        14  )
rate=xk(59)*y(13)*y(4)*xnH
r_f(59,13)=-rate
r_f(59,4)=-rate
r_f(59,2)=rate
r_f(59,14)=rate
! 60)   D   +   H+   ->   HD+  +   ph.
!    (  12      4         15   )
rate=xk(60)*y(12)*y(4)*xnH
r_f(60,12)=-rate
r_f(60,4)=-rate
r_f(60,15)=rate
! 61)   D+  +   H    ->   HD+  +   ph.
!    (  14      1         15   )
rate=xk(61)*y(14)*y(1)*xnH
r_f(61,14)=-rate
r_f(61,1)=-rate
r_f(61,15)=rate
! 62)   HD+ +   e    ->   H    +   D
!    (  15      3         1        12  )
rate=xk(62)*y(15)*y(3)*xnH
r_f(62,15)=-rate
r_f(62,3)=-rate
r_f(62,1)=rate
r_f(62,12)=rate
! 63)   D   +   e    ->   D-   +   ph.
!    (  12      3         16  )
rate=xk(63)*y(12)*y(3)*xnH
r_f(63,12)=-rate
r_f(63,3)=-rate
r_f(63,16)=rate
! 64)   D+  +   D-   ->   2D
!    (  14      16        2*12)
rate=xk(64)*y(14)*y(16)*xnH
r_f(64,14)=-rate
r_f(64,16)=-rate
r_f(64,12)=2.d0*rate
! 65)   H+  +   D-   ->   D    +   H
!    (  4       16        12       1  )
rate=xk(65)*y(4)*y(16)*xnH
r_f(65,4)=-rate
r_f(65,16)=-rate
r_f(65,12)=rate
r_f(65,1)=rate
! 66)   H-  +   D    ->   H    +   D-
!    (  7       12        1        16 )
rate=xk(66)*y(7)*y(12)*xnH
r_f(66,7)=-rate
r_f(66,12)=-rate
r_f(66,1)=rate
r_f(66,16)=rate
! 67)   D-  +   H    ->   D    +   H-
!    (  16      1         12       7  )
rate=xk(67)*y(16)*y(1)*xnH
r_f(67,16)=-rate
r_f(67,1)=-rate
r_f(67,12)=rate
r_f(67,7)=rate
! 68)   D-  +   H    ->   HD   +   e
!    (  16      1         13       3  )
rate=xk(68)*y(16)*y(1)*xnH
r_f(68,16)=-rate
r_f(68,1)=-rate
r_f(68,13)=rate
r_f(68,3)=rate
! 69)   D   +   e    ->   D+   +  2 e
!    ( 12       3         14      2x3    )
rate=xk(69)*y(12)*y(3)*xnH
r_f(69,12)=-rate
r_f(69,14)=rate
r_f(69,3)=rate
! 70)   He+ +   D    ->   D+   +   He
!    (  9      12         14       8   )
rate=xk(70)*y(9)*y(12)*xnH
r_f(70,9)=-rate
r_f(70,12)=-rate
r_f(70,14)=rate
r_f(70,8)=rate
! 71)   He  +   D+   ->   D    +   He+
!    (  8       14        12       9   )
rate=xk(71)*y(8)*y(14)*xnH
r_f(71,8)=-rate
r_f(71,14)=-rate
r_f(71,12)=rate
r_f(71,9)=rate
! 72)   H2+ +   D    ->  HD+   +   H
!    (  5      12        15        1   )
rate=xk(72)*y(5)*y(12)*xnH
r_f(72,5)=-rate
r_f(72,12)=-rate
r_f(72,15)=rate
r_f(72,1)=rate
! 73)   HD+ +   D    ->  HD    +   D+
!    (  15     12        13        14   )
rate=xk(73)*y(15)*y(12)*xnH
r_f(73,15)=-rate
r_f(73,12)=-rate
r_f(73,13)=rate
r_f(73,14)=rate
! 74)   HD+ +   H    ->  H2+   +   D
!    (  15      1        5        12   )
rate=xk(74)*y(15)*y(1)*xnH
r_f(74,15)=-rate
r_f(74,1)=-rate
r_f(74,5)=rate
r_f(74,12)=rate
! 75)   HD  +   e    ->  H     +   D-
!    (  13      3        1         16  )
rate=xk(75)*y(13)*y(3)*xnH
r_f(75,13)=-rate
r_f(75,3)=-rate
r_f(75,1)=rate
r_f(75,16)=rate
! 76)   HD  +   e    ->  D     +   H-
!    (  13      3       12         7   )
rate=xk(76)*y(13)*y(3)*xnH
r_f(76,13)=-rate
r_f(76,3)=-rate
r_f(76,12)=rate
r_f(76,7)=rate
! 77)   H+  +   D-   ->  HD+   +   e
!    (  4       16       15        3   )
rate=xk(77)*y(4)*y(16)*xnH
r_f(77,4)=-rate
r_f(77,16)=-rate
r_f(77,15)=rate
r_f(77,3)=rate
! 78)   D+  +   H-   ->  HD+   +   e
!    (  14      7        15        3  )
rate=xk(78)*y(14)*y(7)*xnH
r_f(78,14)=-rate
r_f(78,7)=-rate
r_f(78,15)=rate
r_f(78,3)=rate
! 79)   D-  +   e    ->  D     + 2 e
!    (  16      3       12       2x3     )
rate=xk(79)*y(16)*y(3)*xnH
r_f(79,16)=-rate
r_f(79,12)=rate
r_f(79,3)=rate
! 80)   D-  +   H    ->  D     +   H    +   e
!    (  16      1       12         1        3 )
rate=xk(80)*y(16)*y(1)*xnH
r_f(80,16)=-rate
r_f(80,12)=rate
r_f(80,3)=rate
! 81)   D-  +   He   ->  D     +   He   +   e
!    (  16      8       12         8        3 )
rate=xk(81)*y(16)*y(8)*xnH
r_f(81,16)=-rate
r_f(81,12)=rate
r_f(81,3)=rate
! 82)   D+  +   H-   ->  D     +   H
!    (  14      7       12         1   )
rate=xk(82)*y(14)*y(7)*xnH
r_f(82,14)=-rate
r_f(82,7)=-rate
r_f(82,12)=rate
r_f(82,1)=rate
! 83)   H2+ +   D-   ->  H2    +   D
!    (  5       16       2        12    )
rate=xk(83)*y(5)*y(16)*xnH
r_f(83,5)=-rate
r_f(83,16)=-rate
r_f(83,2)=rate
r_f(83,12)=rate
! 84)   H2+ +   D-   -> 2 H    +   D
!    (  5       16      2x1       12   )
rate=xk(84)*y(5)*y(16)*xnH
r_f(84,5)=-rate
r_f(84,16)=-rate
r_f(84,1)=2.d0*rate
r_f(84,12)=rate
! 85)   HD+ +   H-   ->   HD   +   H
!    (  15      7         13       1   )
rate=xk(85)*y(15)*y(7)*xnH
r_f(85,15)=-rate
r_f(85,7)=-rate
r_f(85,13)=rate
r_f(85,1)=rate
! 86)   HD+ +   H-   ->   D    + 2 H
!    (  15      7        12      2x1     )
rate=xk(86)*y(15)*y(7)*xnH
r_f(86,15)=-rate
r_f(86,7)=-rate
r_f(86,12)=rate
r_f(86,1)=2.d0*rate
! 87)   HD+ +   D-   ->  HD    +   D
!    (  15      16       13       12   )
rate=xk(87)*y(15)*y(16)*xnH
r_f(87,15)=-rate
r_f(87,16)=-rate
r_f(87,13)=rate
r_f(87,12)=rate
! 88)   HD+ +   D-   ->  2 D   +   H
!    (  15      16       2x12      1    )
rate=xk(88)*y(15)*y(16)*xnH
r_f(88,15)=-rate
r_f(88,16)=-rate
r_f(88,12)=2.d0*rate
r_f(88,1)=rate
! 89)   He+ +   D-   ->  He    +   D
!    (  9       16       8        12   )
rate=xk(89)*y(9)*y(16)*xnH
r_f(89,9)=-rate
r_f(89,16)=-rate
r_f(89,8)=rate
r_f(89,12)=rate
! 90)   D   +   H2+  ->  H2    +   D+
!    ( 12       5        2         14  )
rate=xk(90)*y(12)*y(5)*xnH
r_f(90,12)=-rate
r_f(90,5)=-rate
r_f(90,2)=rate
r_f(90,14)=rate
! 91)  H2+  +   D    ->  HD    +   H+
!    ( 5       12        13        4    )
rate=xk(91)*y(5)*y(12)*xnH
r_f(91,5)=-rate
r_f(91,12)=-rate
r_f(91,13)=rate
r_f(91,4)=rate
! 92)  HD+  +   H    ->  H2    +   D+
!    ( 15       1        2         14  )
rate=xk(92)*y(15)*y(1)*xnH
r_f(92,15)=-rate
r_f(92,1)=-rate
r_f(92,2)=rate
r_f(92,14)=rate
! 93)  H2   +   D+   ->  H2+   +   D
!    ( 2        14       5        12   )
rate=xk(93)*y(2)*y(14)*xnH
r_f(93,2)=-rate
r_f(93,14)=-rate
r_f(93,5)=rate
r_f(93,12)=rate
! 94)  H2   +   D+   ->  HD+   +   H
!    ( 2        14       15        1   )
rate=xk(94)*y(2)*y(14)*xnH
r_f(94,2)=-rate
r_f(94,14)=-rate
r_f(94,15)=rate
r_f(94,1)=rate
! 95)  HD   +   H+   ->  HD+   +   H
!    ( 13       4        15        1   )
rate=xk(95)*y(13)*y(4)*xnH
r_f(95,13)=-rate
r_f(95,4)=-rate
r_f(95,15)=rate
r_f(95,1)=rate
! 96)  HD   +   H+   ->  H2+   +   D
!    ( 13       4        5        12   )
rate=xk(96)*y(13)*y(4)*xnH
r_f(96,13)=-rate
r_f(96,4)=-rate
r_f(96,5)=rate
r_f(96,12)=rate
! 97)  HD   +   D+   ->  HD+   +   D
!    ( 13       14       15       12   )
rate=xk(97)*y(13)*y(14)*xnH
r_f(97,13)=-rate
r_f(97,14)=-rate
r_f(97,15)=rate
r_f(97,12)=rate
! 98)  HD   +   He+  ->  HD+   +   He
!    ( 13       9        15        8   )
rate=xk(98)*y(13)*y(9)*xnH
r_f(98,13)=-rate
r_f(98,9)=-rate
r_f(98,15)=rate
r_f(98,8)=rate
! 99)  HD   +   He+  ->  He    +   H+   +   D
!    ( 13       9        8         4       12 )
rate=xk(99)*y(13)*y(9)*xnH
r_f(99,13)=-rate
r_f(99,9)=-rate
r_f(99,8)=rate
r_f(99,4)=rate
r_f(99,12)=rate
!100)  HD   +   He+  ->  He    +   H    +   D+
!    ( 13       9        8         1        14)
rate=xk(100)*y(13)*y(9)*xnH
r_f(100,13)=-rate
r_f(100,9)=-rate
r_f(100,8)=rate
r_f(100,1)=rate
r_f(100,14)=rate
!D50)  HD   +   H    -> 2 H    +   D
!    ( 13       1       2x1       12 )
rate=xk(50)*y(13)*y(1)*xnH
r_f(50,13)=-rate
r_f(50,1)=rate
r_f(50,12)=rate
!D49)  HD   +   H2   ->  H     +   D    +   H2
!    ( 13       2        1        12        2)
rate=xk(49)*y(13)*y(2)*xnH
r_f(49,13)=-rate
r_f(49,1)=rate
r_f(49,12)=rate
!D48)  HD   +   He   ->  H     +   D    +   He
!    ( 13       8        1        12        8 )
rate=xk(48)*y(13)*y(8)*xnH
r_f(48,13)=-rate
r_f(48,1)=rate
r_f(48,12)=rate
!D47)  HD   +   e    ->  H     +   D    +   e
!    ( 13       3        1        12        3 )
rate=xk(47)*y(13)*y(3)*xnH
r_f(47,13)=-rate
r_f(47,1)=rate
r_f(47,12)=rate


!    C,O reactions
! 101)   H     +   CH    ->   C     +   H2
!        (1         19         17        2)
rate=xk(101)*y(1)*y(19)*xnH
r_f(101,1)=-rate
r_f(101,19)=-rate
r_f(101,17)=rate
r_f(101,2)=rate
! 102)   H     +   CH    ->   C     + 2 H
!        (1         19         17      2*1)
rate=xk(102)*y(1)*y(19)*xnH
r_f(102,19)=-rate
r_f(102,17)=rate
r_f(102,1)=rate
! 103)   H     +   CH2   ->   CH    +   H2
!        (1         20         19        2)
rate=xk(103)*y(1)*y(20)*xnH
r_f(103,1)=-rate
r_f(103,20)=-rate
r_f(103,19)=rate
r_f(103,2)=rate
! 104)   H     +   CH3   ->   CH2   +   H2
!        (1         21         20        2)
rate=xk(104)*y(1)*y(21)*xnH
r_f(104,1)=-rate
r_f(104,21)=-rate
r_f(104,20)=rate
r_f(104,2)=rate
! 105)   H     +   CH4   ->   H2    +   CH3
!        (1         22         2         21)
rate=xk(105)*y(1)*y(22)*xnH
r_f(105,1)=-rate
r_f(105,22)=-rate
r_f(105,2)=rate
r_f(105,21)=rate
! 106)   H     +   OH    ->   H2    +   O
!        (1         32         2         30)
rate=xk(106)*y(1)*y(32)*xnH
r_f(106,1)=-rate
r_f(106,32)=-rate
r_f(106,2)=rate
r_f(106,30)=rate
! 107)   H     +   OH    ->   O     + 2 H
!        (1         32         30      2*1)
rate=xk(107)*y(1)*y(32)*xnH
r_f(107,32)=-rate
r_f(107,30)=rate
r_f(107,1)=rate
! 108)   H     +   H2O   ->   OH    +   H2
!        (1         34         32        2)
rate=xk(108)*y(1)*y(34)*xnH
r_f(108,1)=-rate
r_f(108,34)=-rate
r_f(108,32)=rate
r_f(108,2)=rate
! 109)   H     +   H2O   ->   OH    + 2 H
!        (1         34         32      2*1)
rate=xk(109)*y(1)*y(34)*xnH
r_f(109,34)=-rate
r_f(109,32)=rate
r_f(109,1)=rate
! 110)   H     +   C2    ->   CH    +   C
!        (1         18         19        17)
rate=xk(110)*y(1)*y(18)*xnH
r_f(110,1)=-rate
r_f(110,18)=-rate
r_f(110,19)=rate
r_f(110,17)=rate
! 111)   H     +   CO    ->   C     +   OH
!        (1         33         17        32)
rate=xk(111)*y(1)*y(33)*xnH
r_f(111,1)=-rate
r_f(111,33)=-rate
r_f(111,17)=rate
r_f(111,32)=rate
! 112)   H     +   H2CO  ->   HCO   +   H2
!        (1         38         35        2)
rate=xk(112)*y(1)*y(38)*xnH
r_f(112,1)=-rate
r_f(112,38)=-rate
r_f(112,35)=rate
r_f(112,2)=rate
! 113)   H     +   O2    ->   OH    +   O
!        (1         31         32        30)
rate=xk(113)*y(1)*y(31)*xnH
r_f(113,1)=-rate
r_f(113,31)=-rate
r_f(113,32)=rate
r_f(113,30)=rate
! 114)   H     +   O2    -> 2 O     +   H
!        (1         31       2*30        1)
rate=xk(114)*y(1)*y(31)*xnH
r_f(114,31)=-rate
r_f(114,30)=2.d0*rate
! 115)   H     +   O2H   ->   H2O   +   O
!        (1         36         34        30)
rate=xk(115)*y(1)*y(36)*xnH
r_f(115,1)=-rate
r_f(115,36)=-rate
r_f(115,34)=rate
r_f(115,30)=rate
! 116)   H     +   O2H   ->   H2    +   O2
!        (1         36         2         31)
rate=xk(116)*y(1)*y(36)*xnH
r_f(116,1)=-rate
r_f(116,36)=-rate
r_f(116,2)=rate
r_f(116,31)=rate
! 117)   H     +   O2H   -> 2 OH
!        (1         36       2*32)
rate=xk(117)*y(1)*y(36)*xnH
r_f(117,1)=-rate
r_f(117,36)=-rate
r_f(117,32)=2.d0*rate
! 118)   H     +   H2O2  ->   O2H   +   H2
!        (1         39         36        2)
rate=xk(118)*y(1)*y(39)*xnH
r_f(118,1)=-rate
r_f(118,39)=-rate
r_f(118,36)=rate
r_f(118,2)=rate
! 119)   H     +   CO2   ->   OH    +   CO
!        (1         37         32        33)
rate=xk(119)*y(1)*y(37)*xnH
r_f(119,1)=-rate
r_f(119,37)=-rate
r_f(119,32)=rate
r_f(119,33)=rate
! 120)   C    +   H2    ->   CH    +   H
!        (17        2          19        1)
rate=xk(120)*y(17)*y(2)*xnH
r_f(120,17)=-rate
r_f(120,2)=-rate
r_f(120,19)=rate
r_f(120,1)=rate
! 121)   C    +   OH    ->   CH    +   O
!        (17        32         19        30)
rate=xk(121)*y(17)*y(32)*xnH
r_f(121,17)=-rate
r_f(121,32)=-rate
r_f(121,19)=rate
r_f(121,30)=rate
! 122)   C     +   CO    ->   C2    +   O
!        (17        33         18        30)
rate=xk(122)*y(17)*y(33)*xnH
r_f(122,17)=-rate
r_f(122,33)=-rate
r_f(122,18)=rate
r_f(122,30)=rate
! 123)   O     +   H2    ->   OH    +   H
!        (30        2          32        1)
rate=xk(123)*y(30)*y(2)*xnH
r_f(123,30)=-rate
r_f(123,2)=-rate
r_f(123,32)=rate
r_f(123,1)=rate
! 124)   O     +   CH    ->   OH    +   C
!        (30        19         32        17)
rate=xk(124)*y(30)*y(19)*xnH
r_f(124,30)=-rate
r_f(124,19)=-rate
r_f(124,32)=rate
r_f(124,17)=rate
! 125)   O     +   CH2   ->   CH    +   OH
!        (30        20         19        32)
rate=xk(125)*y(30)*y(20)*xnH
r_f(125,30)=-rate
r_f(125,20)=-rate
r_f(125,19)=rate
r_f(125,32)=rate
! 126)   O     +   CH2   ->   HCO   +   H
!        (30        20         35        1)
rate=xk(126)*y(30)*y(20)*xnH
r_f(126,30)=-rate
r_f(126,20)=-rate
r_f(126,35)=rate
r_f(126,1)=rate
! 127)   O     +   CH4   ->   CH3   +   OH
!        (30        22         21        32)
rate=xk(127)*y(30)*y(22)*xnH
r_f(127,30)=-rate
r_f(127,22)=-rate
r_f(127,21)=rate
r_f(127,32)=rate
!  128)   O     +   H2O   -> 2 OH
!        (30        34       2*32)
rate=xk(128)*y(30)*y(34)*xnH
r_f(128,30)=-rate
r_f(128,34)=-rate
r_f(128,32)=2.d0*rate
!  129)   O     +   H2CO  ->   HCO   +   OH
!        (30        38         35        32)
rate=xk(129)*y(30)*y(38)*xnH
r_f(129,30)=-rate
r_f(129,38)=-rate
r_f(129,35)=rate
r_f(129,32)=rate
!  130)   O     +   H2O2  ->   O2H   +   OH
!        (30        39         36        32)
rate=xk(130)*y(30)*y(39)*xnH
r_f(130,30)=-rate
r_f(130,39)=-rate
r_f(130,36)=rate
r_f(130,32)=rate
!  131)   O     +   CO2   ->   CO    +   O2
!        (30        37         33        31)
rate=xk(131)*y(30)*y(37)*xnH
r_f(131,30)=-rate
r_f(131,37)=-rate
r_f(131,33)=rate
r_f(131,31)=rate
!  132)   H+    +   O     ->   O+    +   H
!        (4         30         40        1)
rate=xk(132)*y(4)*y(30)*xnH
r_f(132,4)=-rate
r_f(132,30)=-rate
r_f(132,40)=rate
r_f(132,1)=rate
!  133)   H2    +   CH    ->   CH2   +   H
!        (2         19         20        1)
rate=xk(133)*y(2)*y(19)*xnH
r_f(133,2)=-rate
r_f(133,19)=-rate
r_f(133,20)=rate
r_f(133,1)=rate
!  134)   H2    +   CH    ->   C    +   H     +   H2
!        (2         19         17        1         2)
rate=xk(134)*y(2)*y(19)*xnH
r_f(134,19)=-rate
r_f(134,17)=rate
r_f(134,1)=rate
!  135)   H2    +   CH2   ->   CH3   +   H
!        (2         20         21        1)
rate=xk(135)*y(2)*y(20)*xnH
r_f(135,2)=-rate
r_f(135,20)=-rate
r_f(135,21)=rate
r_f(135,1)=rate
!  136)   H2    +   CH3   ->   CH4   +   H
!        (2         21         22        1)
rate=xk(136)*y(2)*y(21)*xnH
r_f(136,2)=-rate
r_f(136,21)=-rate
r_f(136,22)=rate
r_f(136,1)=rate
!  137)   H2    +   OH    ->   H2O   +   H
!        (2         32         34        1)
rate=xk(137)*y(2)*y(32)*xnH
r_f(137,2)=-rate
r_f(137,32)=-rate
r_f(137,34)=rate
r_f(137,1)=rate
!  138)   H2    +   OH    ->   O     +   H     +   H2
!        (2         32         30        1         2)
rate=xk(138)*y(2)*y(32)*xnH
r_f(138,32)=-rate
r_f(138,30)=rate
r_f(138,1)=rate
!  139)   H2    +   H2O   ->   OH    +   H     +   H2
!        (2         34         32        1         2)
rate=xk(139)*y(2)*y(34)*xnH
r_f(139,34)=-rate
r_f(139,32)=rate
r_f(139,1)=rate
!  140)   H2    +   O2    ->   O2H   +   H
!        (2         31         36        1)
rate=xk(140)*y(2)*y(31)*xnH
r_f(140,2)=-rate
r_f(140,31)=-rate
r_f(140,36)=rate
r_f(140,1)=rate
!  141)   H2    +   O2    -> 2 OH
!        (2         31       2*32)
rate=xk(141)*y(2)*y(31)*xnH
r_f(141,2)=-rate
r_f(141,31)=-rate
r_f(141,32)=2.d0*rate
!  142)   H2    +   O2    -> 2 O     +   H2
!        (2         31       2*30        2)
rate=xk(142)*y(2)*y(31)*xnH
r_f(142,31)=-rate
r_f(142,30)=2.d0*rate
!  143)   H2    +   O2H   ->   H2O2  +   H
!        (2         36         39        1)
rate=xk(143)*y(2)*y(36)*xnH
r_f(143,2)=-rate
r_f(143,36)=-rate
r_f(143,39)=rate
r_f(143,1)=rate
!  144)   H2    +   CO2   ->   H2O   +   CO
!        (2         37         34        33)
rate=xk(144)*y(2)*y(37)*xnH
r_f(144,2)=-rate
r_f(144,37)=-rate
r_f(144,34)=rate
r_f(144,33)=rate
!  145)   H3+   +   O2    ->   O2H+  +   H2
!        (6         31         46        2)
rate=xk(145)*y(6)*y(31)*xnH
r_f(145,6)=-rate
r_f(145,31)=-rate
r_f(145,46)=rate
r_f(145,2)=rate
!  146)   C+    +   H2    ->   CH+   +   H
!        (23        2          25        1)
rate=xk(146)*y(23)*y(2)*xnH
r_f(146,23)=-rate
r_f(146,2)=-rate
r_f(146,25)=rate
r_f(146,1)=rate
!  147)   CH    +   CH4   ->   CH2   +   CH3
!        (19        22         20        21)
rate=xk(147)*y(19)*y(22)*xnH
r_f(147,19)=-rate
r_f(147,22)=-rate
r_f(147,20)=rate
r_f(147,21)=rate
!  148)   CH    +   OH    ->   HCO   +   H
!        (19        32         35        1)
rate=xk(148)*y(19)*y(32)*xnH
r_f(148,19)=-rate
r_f(148,32)=-rate
r_f(148,35)=rate
r_f(148,1)=rate
!  149)   CH    +   HCO   ->   CH2   +   CO
!        (19        35         20        33)
rate=xk(149)*y(19)*y(35)*xnH
r_f(149,19)=-rate
r_f(149,35)=-rate
r_f(149,20)=rate
r_f(149,33)=rate
!  150)   CH    +   H2CO  ->   CH2   +   HCO
!        (19        38         20        35)
rate=xk(150)*y(19)*y(38)*xnH
r_f(150,19)=-rate
r_f(150,38)=-rate
r_f(150,20)=rate
r_f(150,35)=rate
!  151)   CH    +   O2    ->   HCO   +   O
!        (19        31         35        30)
rate=xk(151)*y(19)*y(31)*xnH
r_f(151,19)=-rate
r_f(151,31)=-rate
r_f(151,35)=rate
r_f(151,30)=rate
!  152)   CH    +   O2H   ->   CH2   +   O2
!        (19        36         20        31)
rate=xk(152)*y(19)*y(36)*xnH
r_f(152,19)=-rate
r_f(152,36)=-rate
r_f(152,20)=rate
r_f(152,31)=rate
!  153)   CH    +   O2H   ->   HCO   +   OH
!        (19        36         35        32)
rate=xk(153)*y(19)*y(36)*xnH
r_f(153,19)=-rate
r_f(153,36)=-rate
r_f(153,35)=rate
r_f(153,32)=rate
!  154)   CH    +   CO2   ->   HCO   +   CO
!        (19        37         35        33)
rate=xk(154)*y(19)*y(37)*xnH
r_f(154,19)=-rate
r_f(154,37)=-rate
r_f(154,35)=rate
r_f(154,33)=rate
!  155) 2 CH2             ->   CH    +   CH3
!      (2*20                   19        21)
rate=xk(155)*(y(20)**2)*xnH
r_f(155,20)=-2.d0*rate
r_f(155,19)=rate
r_f(155,21)=rate
!  156)   CH2   +   CH4   -> 2 CH3
!        (20        22       2*21)
rate=xk(156)*y(20)*y(22)*xnH
r_f(156,20)=-rate
r_f(156,22)=-rate
r_f(156,21)=2.d0*rate
!  157)   CH2   +   HCO   ->   CH3   +   CO
!        (20        35         21        33)
rate=xk(157)*y(20)*y(35)*xnH
r_f(157,20)=-rate
r_f(157,35)=-rate
r_f(157,21)=rate
r_f(157,33)=rate
!  158)   CH2   +   H2CO  ->   CH3   +   HCO
!        (20        38         21        35)
rate=xk(158)*y(20)*y(38)*xnH
r_f(158,20)=-rate
r_f(158,38)=-rate
r_f(158,21)=rate
r_f(158,35)=rate
!  159)   CH2+  +   H     ->   CH+   +   H2
!        (26        1          25        2)
rate=xk(159)*y(26)*y(1)*xnH
r_f(159,26)=-rate
r_f(159,1)=-rate
r_f(159,25)=rate
r_f(159,2)=rate
!  160)   CH3   +   H2CO  ->   CH4   +   HCO
!        (21        38         22        35)
rate=xk(160)*y(21)*y(38)*xnH
r_f(160,21)=-rate
r_f(160,38)=-rate
r_f(160,22)=rate
r_f(160,35)=rate
!  161)   CH3+  +   H     ->   CH2+  +   H2
!        (27        1          26        2)
rate=xk(161)*y(27)*y(1)*xnH
r_f(161,27)=-rate
r_f(161,1)=-rate
r_f(161,26)=rate
r_f(161,2)=rate
!  162)   OH    +   CH2   ->   CH3   +   O
!        (32        20         21        30)
rate=xk(162)*y(32)*y(20)*xnH
r_f(162,32)=-rate
r_f(162,20)=-rate
r_f(162,21)=rate
r_f(162,30)=rate
!  163)   OH    +   CH2   ->   H2O   +   CH
!        (32        20         34        19)
rate=xk(163)*y(32)*y(20)*xnH
r_f(163,32)=-rate
r_f(163,20)=-rate
r_f(163,34)=rate
r_f(163,19)=rate
! 164)   OH    +   CH2   ->   H2CO  +   H
!        (32        20         38        1)
rate=xk(164)*y(32)*y(20)*xnH
r_f(164,32)=-rate
r_f(164,20)=-rate
r_f(164,38)=rate
r_f(164,1)=rate
! 165)   OH    +   CH3   ->   H2O   +   CH2
!        (32        21         34        20)
rate=xk(165)*y(32)*y(21)*xnH
r_f(165,32)=-rate
r_f(165,21)=-rate
r_f(165,34)=rate
r_f(165,20)=rate
! 166)   OH    +   CH4   ->   H2O   +   CH3
!        (32        22         34        21)
rate=xk(166)*y(32)*y(22)*xnH
r_f(166,32)=-rate
r_f(166,22)=-rate
r_f(166,34)=rate
r_f(166,21)=rate
! 167) 2 OH              ->   H2O   +   O
!      (2*32                   34        30)
rate=xk(167)*(y(32)**2)*xnH
r_f(167,32)=-2.d0*rate
r_f(167,34)=rate
r_f(167,30)=rate
! 168)   OH    +   CO    ->   CO2   +   H
!        (32        33         37        1)
rate=xk(168)*y(32)*y(33)*xnH
r_f(168,32)=-rate
r_f(168,33)=-rate
r_f(168,37)=rate
r_f(168,1)=rate
! 169)   OH    +   H2O2  ->   H2O   +   O2H
!        (32        39         34        36)
rate=xk(169)*y(32)*y(39)*xnH
r_f(169,32)=-rate
r_f(169,39)=-rate
r_f(169,34)=rate
r_f(169,36)=rate
! 170)   H2O   +   CH3   ->   CH4   +   OH
!        (34        21         22        32)
rate=xk(170)*y(34)*y(21)*xnH
r_f(170,34)=-rate
r_f(170,21)=-rate
r_f(170,22)=rate
r_f(170,32)=rate
! 171)   CO    +   O2    ->   CO2   +   O
!        (33        31         37        30)
rate=xk(171)*y(33)*y(31)*xnH
r_f(171,33)=-rate
r_f(171,31)=-rate
r_f(171,37)=rate
r_f(171,30)=rate
! 172)   CO    +   O2H   ->   CO2   +   OH
!        (33        36         37        32)
rate=xk(172)*y(33)*y(36)*xnH
r_f(172,33)=-rate
r_f(172,36)=-rate
r_f(172,37)=rate
r_f(172,32)=rate
! 173)   O2    +   CH2   ->   HCO   +   OH
!        (31        20         35        32)
rate=xk(173)*y(31)*y(20)*xnH
r_f(173,31)=-rate
r_f(173,20)=-rate
r_f(173,35)=rate
r_f(173,32)=rate
! 174)   O2    +   CH2   ->   H2CO  +   O
!        (31        20         38        30)
rate=xk(174)*y(31)*y(20)*xnH
r_f(174,31)=-rate
r_f(174,20)=-rate
r_f(174,38)=rate
r_f(174,30)=rate
! 175)   O2    +   CH3   ->   O2H   +   CH2
!        (31        21         36        20)
rate=xk(175)*y(31)*y(21)*xnH
r_f(175,31)=-rate
r_f(175,21)=-rate
r_f(175,36)=rate
r_f(175,20)=rate
! 176)   O2    +   CH3   ->   H2CO  +   OH
!        (31        21         38        32)
rate=xk(176)*y(31)*y(21)*xnH
r_f(176,31)=-rate
r_f(176,21)=-rate
r_f(176,38)=rate
r_f(176,32)=rate
! 177)   O2    +   CH4   ->   CH3   +   O2H
!        (31        22         21        36)
rate=xk(177)*y(31)*y(22)*xnH
r_f(177,31)=-rate
r_f(177,22)=-rate
r_f(177,21)=rate
r_f(177,36)=rate
! 178)   O2    +   HCO   ->   O2H   +   CO
!        (31        35         36        33)
rate=xk(178)*y(31)*y(35)*xnH
r_f(178,31)=-rate
r_f(178,35)=-rate
r_f(178,36)=rate
r_f(178,33)=rate
! 179)   O2H   +   CH3   ->   CH4   +   O2
!        (36        21         22        31)
rate=xk(179)*y(36)*y(21)*xnH
r_f(179,36)=-rate
r_f(179,21)=-rate
r_f(179,22)=rate
r_f(179,31)=rate
! 180)   O2H   +   H2O   ->   H2O2  +   OH
!        (36        34         39        32)
rate=xk(180)*y(36)*y(34)*xnH
r_f(180,36)=-rate
r_f(180,34)=-rate
r_f(180,39)=rate
r_f(180,32)=rate
! 181)   O2H   +   HCO   ->   H2CO  +   O2
!        (36        35         38        31)
rate=xk(181)*y(36)*y(35)*xnH
r_f(181,36)=-rate
r_f(181,35)=-rate
r_f(181,38)=rate
r_f(181,31)=rate
! 182)   O2H   +   H2CO  ->   H2O2  +   HCO
!        (36        38         39        35)
rate=xk(182)*y(36)*y(38)*xnH
r_f(182,36)=-rate
r_f(182,38)=-rate
r_f(182,39)=rate
r_f(182,35)=rate
! 183) 2 O2H             ->   H2O2  +   O2
!      (2*36                   39        31)
rate=xk(183)*(y(36)**2)*xnH
r_f(183,36)=-2.d0*rate
r_f(183,39)=rate
r_f(183,31)=rate
! 184)   H     +   HCO   ->   CO    +   H2
!        (1         35         33        2)
rate=xk(184)*y(1)*y(35)*xnH
r_f(184,1)=-rate
r_f(184,35)=-rate
r_f(184,33)=rate
r_f(184,2)=rate
! 185)   !    +   H     ->   CH    +   ph.
!        (17        1          19)
rate=xk(185)*y(17)*y(1)*xnH
r_f(185,17)=-rate
r_f(185,1)=-rate
r_f(185,19)=rate
! 186) 2 C               ->   C2    +   ph.
!      (2*17                   18)
rate=xk(186)*(y(17)**2)*xnH
r_f(186,17)=-2.d0*rate
r_f(186,18)=rate
! 187)   C     +   O     ->   CO    +   ph.
!        (17        30         33)
rate=xk(187)*y(17)*y(30)*xnH
r_f(187,17)=-rate
r_f(187,30)=-rate
r_f(187,33)=rate
! 188)   C     +   CH    ->   C2    +   H
!        (17        19         18        1)
rate=xk(188)*y(17)*y(19)*xnH
r_f(188,17)=-rate
r_f(188,19)=-rate
r_f(188,18)=rate
r_f(188,1)=rate
! 189)   C     +   OH    ->   CO    +   H
!        (17        32         33        1)
rate=xk(189)*y(17)*y(32)*xnH
r_f(189,17)=-rate
r_f(189,32)=-rate
r_f(189,33)=rate
r_f(189,1)=rate
! 190)   C     +   HCO   ->   CH    +   CO
!        (17        35         19        33)
rate=xk(190)*y(17)*y(35)*xnH
r_f(190,17)=-rate
r_f(190,35)=-rate
r_f(190,19)=rate
r_f(190,33)=rate
! 191)   C     +   O2    ->   CO    +   O
!        (17        31         33        30)
rate=xk(191)*y(17)*y(31)*xnH
r_f(191,17)=-rate
r_f(191,31)=-rate
r_f(191,33)=rate
r_f(191,30)=rate
! 192)   CH3   +   HCO   ->   CH4   +   CO
!        (21        35         22        33)
rate=xk(192)*y(21)*y(35)*xnH
r_f(192,21)=-rate
r_f(192,35)=-rate
r_f(192,22)=rate
r_f(192,33)=rate
! 193)   O     +   H     ->   OH    +   ph.
!        (30        1          32)
rate=xk(193)*y(30)*y(1)*xnH
r_f(193,30)=-rate
r_f(193,1)=-rate
r_f(193,32)=rate
! 194) 2 O               ->   O2    +   ph.
!      (2*30                   31)
rate=xk(194)*(y(30)**2)*xnH
r_f(194,30)=-2.d0*rate
r_f(194,31)=rate
! 195)   O     +   CH    ->   CO    +   H
!        (30        19         33        1)
rate=xk(195)*y(30)*y(19)*xnH
r_f(195,30)=-rate
r_f(195,19)=-rate
r_f(195,33)=rate
r_f(195,1)=rate
! 196)   O     +   CH    ->   HCO+  +   e
!        (30        19         45        3)
rate=xk(196)*y(30)*y(19)*xnH
r_f(196,30)=-rate
r_f(196,19)=-rate
r_f(196,45)=rate
r_f(196,3)=rate
! 197)   O     +   CH2   ->   CO    + 2 H
!        (30        20         33      2*1)
rate=xk(197)*y(30)*y(20)*xnH
r_f(197,30)=-rate
r_f(197,20)=-rate
r_f(197,33)=rate
r_f(197,1)=2.d0*rate
! 198)   O     +   CH3   ->   H2CO  +   H
!        (30        21         38        1)
rate=xk(198)*y(30)*y(21)*xnH
r_f(198,30)=-rate
r_f(198,21)=-rate
r_f(198,38)=rate
r_f(198,1)=rate
! 199)   O     +   OH    ->   O2    +   H
!        (30        32         31        1)
rate=xk(199)*y(30)*y(32)*xnH
r_f(199,30)=-rate
r_f(199,32)=-rate
r_f(199,31)=rate
r_f(199,1)=rate
! 200)   O     +   C2    ->   CO    +   C
!        (30        18         33        17)
rate=xk(200)*y(30)*y(18)*xnH
r_f(200,30)=-rate
r_f(200,18)=-rate
r_f(200,33)=rate
r_f(200,17)=rate
! 201)   O     +   HCO   ->   CO2   +   H
!        (30        35         37        1)
rate=xk(201)*y(30)*y(35)*xnH
r_f(201,30)=-rate
r_f(201,35)=-rate
r_f(201,37)=rate
r_f(201,1)=rate
! 202)   O     +   HCO   ->   OH    +   CO
!        (30        35         32        33)
rate=xk(202)*y(30)*y(35)*xnH
r_f(202,30)=-rate
r_f(202,35)=-rate
r_f(202,32)=rate
r_f(202,33)=rate
! 203)   O     +   O2H   ->   OH    +   O2
!        (30        36         32        31)
rate=xk(203)*y(30)*y(36)*xnH
r_f(203,30)=-rate
r_f(203,36)=-rate
r_f(203,32)=rate
r_f(203,31)=rate
! 204)   OH    +   HCO   ->   H2O   +   CO
!        (32        35         34        33)
rate=xk(204)*y(32)*y(35)*xnH
r_f(204,32)=-rate
r_f(204,35)=-rate
r_f(204,34)=rate
r_f(204,33)=rate
! 205)   OH    +   H2CO  ->   H2O   +   HCO
!        (32        38         34        35)
rate=xk(205)*y(32)*y(38)*xnH
r_f(205,32)=-rate
r_f(205,38)=-rate
r_f(205,34)=rate
r_f(205,35)=rate
! 206)   OH    +   O2H   ->   H2O   +   O2
!        (32        36         34        31)
rate=xk(206)*y(32)*y(36)*xnH
r_f(206,32)=-rate
r_f(206,36)=-rate
r_f(206,34)=rate
r_f(206,31)=rate
! 207) 2 HCO             ->   H2CO  +   CO
!      (2*35                   38        33)
rate=xk(207)*(y(35)**2)*xnH
r_f(207,35)=-2.d0*rate
r_f(207,38)=rate
r_f(207,33)=rate
!208)   H+    +   CH    ->   CH+   +   H
!        (4         19         25        1)
rate=xk(208)*y(4)*y(19)*xnH
r_f(208,4)=-rate
r_f(208,19)=-rate
r_f(208,25)=rate
r_f(208,1)=rate
!209)   H+    +   CH2   ->   CH+   +   H2
!        (4         20         25        2)
rate=xk(209)*y(4)*y(20)*xnH
r_f(209,4)=-rate
r_f(209,20)=-rate
r_f(209,25)=rate
r_f(209,2)=rate
!210)   H+    +   CH2   ->   CH2+  +   H
!        (4         20         26        1)
rate=xk(210)*y(4)*y(20)*xnH
r_f(210,4)=-rate
r_f(210,20)=-rate
r_f(210,26)=rate
r_f(210,1)=rate
!211)   H+    +   CH3   ->   CH3+  +   H
!        (4         21         27        1)
rate=xk(211)*y(4)*y(21)*xnH
r_f(211,4)=-rate
r_f(211,21)=-rate
r_f(211,27)=rate
r_f(211,1)=rate
!212)   H+    +   CH4   ->   CH3+  +   H2
!        (4         22         27        2)
rate=xk(212)*y(4)*y(22)*xnH
r_f(212,4)=-rate
r_f(212,22)=-rate
r_f(212,27)=rate
r_f(212,2)=rate
!213)   H+    +   CH4   ->   CH4+  +   H
!        (4         22         28        1)
rate=xk(213)*y(4)*y(22)*xnH
r_f(213,4)=-rate
r_f(213,22)=-rate
r_f(213,28)=rate
r_f(213,1)=rate
!214)   H+    +   OH    ->   OH+   +   H
!        (4         32         42        1)
rate=xk(214)*y(4)*y(32)*xnH
r_f(214,4)=-rate
r_f(214,32)=-rate
r_f(214,42)=rate
r_f(214,1)=rate
!215)   H+    +   H2O   ->   H2O+  +   H
!        (4         34         44        1)
rate=xk(215)*y(4)*y(34)*xnH
r_f(215,4)=-rate
r_f(215,34)=-rate
r_f(215,44)=rate
r_f(215,1)=rate
!216)   H+    +   C2    ->   C2+   +   H
!        (4         18         24        1)
rate=xk(216)*y(4)*y(18)*xnH
r_f(216,4)=-rate
r_f(216,18)=-rate
r_f(216,24)=rate
r_f(216,1)=rate
!217)   H+    +   HCO   ->   CO+   +   H2
!        (4         35         43        2)
rate=xk(217)*y(4)*y(35)*xnH
r_f(217,4)=-rate
r_f(217,35)=-rate
r_f(217,43)=rate
r_f(217,2)=rate
!218)   H+    +   HCO   ->   H2+   +   CO
!        (4         35         5         33)
rate=xk(218)*y(4)*y(35)*xnH
r_f(218,4)=-rate
r_f(218,35)=-rate
r_f(218,5)=rate
r_f(218,33)=rate
!219)   H+    +   HCO   ->   HCO+  +   H
!        (4         35         45        1)
rate=xk(219)*y(4)*y(35)*xnH
r_f(219,4)=-rate
r_f(219,35)=-rate
r_f(219,45)=rate
r_f(219,1)=rate
!220)   H+    +   H2CO  ->   H2CO+ +   H
!        (4         38         48        1)
rate=xk(220)*y(4)*y(38)*xnH
r_f(220,4)=-rate
r_f(220,38)=-rate
r_f(220,48)=rate
r_f(220,1)=rate
!221)   H+    +   H2CO  ->   HCO+  +   H2
!        (4         38         45        2)
rate=xk(221)*y(4)*y(38)*xnH
r_f(221,4)=-rate
r_f(221,38)=-rate
r_f(221,45)=rate
r_f(221,2)=rate
!222)   H+    +   O2    ->   O2+   +   H
!        (4         31         41        1)
rate=xk(222)*y(4)*y(31)*xnH
r_f(222,4)=-rate
r_f(222,31)=-rate
r_f(222,41)=rate
r_f(222,1)=rate
!223)   H+    +   CO2   ->   HCO+  +   O
!        (4         37         45        30)
rate=xk(223)*y(4)*y(37)*xnH
r_f(223,4)=-rate
r_f(223,37)=-rate
r_f(223,45)=rate
r_f(223,30)=rate
!224)   H-    +   C     ->   CH    +   e
!        (7         17         19        3)
rate=xk(224)*y(7)*y(17)*xnH
r_f(224,7)=-rate
r_f(224,17)=-rate
r_f(224,19)=rate
r_f(224,3)=rate
!225)   H-    +   O     ->   OH    +   e
!        (7         30         32        3)
rate=xk(225)*y(7)*y(30)*xnH
r_f(225,7)=-rate
r_f(225,30)=-rate
r_f(225,32)=rate
r_f(225,3)=rate
!226)   H-    +   CH    ->   CH2   +   e
!        (7         19         20        3)
rate=xk(226)*y(7)*y(19)*xnH
r_f(226,7)=-rate
r_f(226,19)=-rate
r_f(226,20)=rate
r_f(226,3)=rate
!227)   H-    +   CH2   ->   CH3   +   e
!        (7         20         21        3)
rate=xk(227)*y(7)*y(20)*xnH
r_f(227,7)=-rate
r_f(227,20)=-rate
r_f(227,21)=rate
r_f(227,3)=rate
!228)   H-    +   CH3   ->   CH4   +   e
!        (7         21         22        3)
rate=xk(228)*y(7)*y(21)*xnH
r_f(228,7)=-rate
r_f(228,21)=-rate
r_f(228,22)=rate
r_f(228,3)=rate
!229)   H-    +   OH    ->   H2O   +   e
!        (7         32         34        3)
rate=xk(229)*y(7)*y(32)*xnH
r_f(229,7)=-rate
r_f(229,32)=-rate
r_f(229,34)=rate
r_f(229,3)=rate
!230)   H-    +   CO    ->   HCO   +   e
!        (7         33         35        3)
rate=xk(230)*y(7)*y(33)*xnH
r_f(230,7)=-rate
r_f(230,33)=-rate
r_f(230,35)=rate
r_f(230,3)=rate
!231)   H-    +   HCO   ->   H2CO  +   e
!        (7         35         38        3)
rate=xk(231)*y(7)*y(35)*xnH
r_f(231,7)=-rate
r_f(231,35)=-rate
r_f(231,38)=rate
r_f(231,3)=rate
!232)   H2+   +   C     ->   CH+   +   H
!        (5         17         25        1)
rate=xk(232)*y(5)*y(17)*xnH
r_f(232,5)=-rate
r_f(232,17)=-rate
r_f(232,25)=rate
r_f(232,1)=rate
!233)   H2+   +   O     ->   OH+   +   H
!        (5         30         42        1)
rate=xk(233)*y(5)*y(30)*xnH
r_f(233,5)=-rate
r_f(233,30)=-rate
r_f(233,42)=rate
r_f(233,1)=rate
!234)   H2+   +   CH    ->   CH+   +   H2
!        (5         19         25        2)
rate=xk(234)*y(5)*y(19)*xnH
r_f(234,5)=-rate
r_f(234,19)=-rate
r_f(234,25)=rate
r_f(234,2)=rate
!235)   H2+   +   CH    ->   CH2+  +   H
!        (5         19         26        1)
rate=xk(235)*y(5)*y(19)*xnH
r_f(235,5)=-rate
r_f(235,19)=-rate
r_f(235,26)=rate
r_f(235,1)=rate
!236)   H2+   +   CH2   ->   CH3+  +   H
!        (5         20         27        1)
rate=xk(236)*y(5)*y(20)*xnH
r_f(236,5)=-rate
r_f(236,20)=-rate
r_f(236,27)=rate
r_f(236,1)=rate
!237)   H2+   +   CH2   ->   CH2+  +   H2
!        (5         20         26        2)
rate=xk(237)*y(5)*y(20)*xnH
r_f(237,5)=-rate
r_f(237,20)=-rate
r_f(237,26)=rate
r_f(237,2)=rate
!238)   H2+   +   CH4   ->   CH4+  +   H2
!        (5         22         28        2)
rate=xk(238)*y(5)*y(22)*xnH
r_f(238,5)=-rate
r_f(238,22)=-rate
r_f(238,28)=rate
r_f(238,2)=rate
!239)   H2+   +   CH4   ->   CH5+  +   H
!        (5         22         29        1)
rate=xk(239)*y(5)*y(22)*xnH
r_f(239,5)=-rate
r_f(239,22)=-rate
r_f(239,29)=rate
r_f(239,1)=rate
!240)   H2+   +   CH4   ->   CH3+  +   H     +  H2
!        (5         22         27        1        2)
rate=xk(240)*y(5)*y(22)*xnH
r_f(240,5)=-rate
r_f(240,22)=-rate
r_f(240,27)=rate
r_f(240,1)=rate
r_f(240,2)=rate
!241)   H2+   +   OH    ->   OH+   +   H2
!        (5         32         42        2)
rate=xk(241)*y(5)*y(32)*xnH
r_f(241,5)=-rate
r_f(241,32)=-rate
r_f(241,42)=rate
r_f(241,2)=rate
!242)   H2+   +   OH    ->   H2O+  +   H
!        (5         32         44        1)
rate=xk(242)*y(5)*y(32)*xnH
r_f(242,5)=-rate
r_f(242,32)=-rate
r_f(242,44)=rate
r_f(242,1)=rate
!243)   H2+   +   H2O   ->   H2O+  +   H2
!        (5         34         44        2)
rate=xk(243)*y(5)*y(34)*xnH
r_f(243,5)=-rate
r_f(243,34)=-rate
r_f(243,44)=rate
r_f(243,2)=rate
!244)   H2+   +   H2O   ->   H3O+  +   H
!        (5         34         47        1)
rate=xk(244)*y(5)*y(34)*xnH
r_f(244,5)=-rate
r_f(244,34)=-rate
r_f(244,47)=rate
r_f(244,1)=rate
!245)   H2+   +   C2    ->   C2+   +   H2
!        (5         18         24        2)
rate=xk(245)*y(5)*y(18)*xnH
r_f(245,5)=-rate
r_f(245,18)=-rate
r_f(245,24)=rate
r_f(245,2)=rate
!246)   H2+   +   CO    ->   CO+   +   H2
!        (5         33         43        2)
rate=xk(246)*y(5)*y(33)*xnH
r_f(246,5)=-rate
r_f(246,33)=-rate
r_f(246,43)=rate
r_f(246,2)=rate
!247)   H2+   +   CO    ->   HCO+  +   H
!        (5         33         45        1)
rate=xk(247)*y(5)*y(33)*xnH
r_f(247,5)=-rate
r_f(247,33)=-rate
r_f(247,45)=rate
r_f(247,1)=rate
!248)   H2+   +   HCO   ->   H3+   +   CO
!        (5         35         6         33)
rate=xk(248)*y(5)*y(35)*xnH
r_f(248,5)=-rate
r_f(248,35)=-rate
r_f(248,6)=rate
r_f(248,33)=rate
!249)   H2+   +   HCO   ->   HCO+  +   H2
!        (5         35         45        2)
rate=xk(249)*y(5)*y(35)*xnH
r_f(249,5)=-rate
r_f(249,35)=-rate
r_f(249,45)=rate
r_f(249,2)=rate
!250)   H2+   +   H2CO  ->   H2CO+ +   H2
!        (5         38         48        2)
rate=xk(250)*y(5)*y(38)*xnH
r_f(250,5)=-rate
r_f(250,38)=-rate
r_f(250,48)=rate
r_f(250,2)=rate
!251)   H2+   +   H2CO  ->   HCO+  +   H     +   H2
!        (5         38         45        1         2)
rate=xk(251)*y(5)*y(38)*xnH
r_f(251,5)=-rate
r_f(251,38)=-rate
r_f(251,45)=rate
r_f(251,1)=rate
r_f(251,2)=rate
!252)   H2+   +   O2    ->   O2H+  +   H
!        (5         31         46        1)
rate=xk(252)*y(5)*y(31)*xnH
r_f(252,5)=-rate
r_f(252,31)=-rate
r_f(252,46)=rate
r_f(252,1)=rate
!253)   H2+   +   O2    ->   O2+   +   H2
!        (5         31         41        2)
rate=xk(253)*y(5)*y(31)*xnH
r_f(253,5)=-rate
r_f(253,31)=-rate
r_f(253,41)=rate
r_f(253,2)=rate
!254)   H2+   +   CO2   ->   HCO2+ +   H
!        (5         37         49        1)
rate=xk(254)*y(5)*y(37)*xnH
r_f(254,5)=-rate
r_f(254,37)=-rate
r_f(254,49)=rate
r_f(254,1)=rate
!255)   H3+   +   C    ->   CH+   +   H2
!        (6         17         25        2)
rate=xk(255)*y(6)*y(17)*xnH
r_f(255,6)=-rate
r_f(255,17)=-rate
r_f(255,25)=rate
r_f(255,2)=rate
! 256)   H3+   +   O     ->   OH+   +   H2
!        (6         30         42        2)
rate=xk(256)*y(6)*y(30)*xnH
r_f(256,6)=-rate
r_f(256,30)=-rate
r_f(256,42)=rate
r_f(256,2)=rate
! 257)   H3+   +   CH    ->   CH2+  +   H2
!        (6         19         26        2)
rate=xk(257)*y(6)*y(19)*xnH
r_f(257,6)=-rate
r_f(257,19)=-rate
r_f(257,26)=rate
r_f(257,2)=rate
! 258)   H3+   +   CH2   ->   CH3+  +   H2
!        (6         20         27        2)
rate=xk(258)*y(6)*y(20)*xnH
r_f(258,6)=-rate
r_f(258,20)=-rate
r_f(258,27)=rate
r_f(258,2)=rate
! 259)   H3+   +   CH3   ->   CH4+  +   H2
!        (6         21         28        2)
rate=xk(259)*y(6)*y(21)*xnH
r_f(259,6)=-rate
r_f(259,21)=-rate
r_f(259,28)=rate
r_f(259,2)=rate
! 260)   H3+   +   CH4   ->   CH5+  +   H2
!        (6         22         29        2)
rate=xk(260)*y(6)*y(22)*xnH
r_f(260,6)=-rate
r_f(260,22)=-rate
r_f(260,29)=rate
r_f(260,2)=rate
! 261)   H3+   +   OH    ->   H2O+  +   H2
!        (6         32         44        2)
rate=xk(261)*y(6)*y(32)*xnH
r_f(261,6)=-rate
r_f(261,32)=-rate
r_f(261,44)=rate
r_f(261,2)=rate
! 262)   H3+   +   H2O   ->   H3O+  +   H2
!        (6         34         47        2)
rate=xk(262)*y(6)*y(34)*xnH
r_f(262,6)=-rate
r_f(262,34)=-rate
r_f(262,47)=rate
r_f(262,2)=rate
! 263)   H3+   +   CO    ->   HCO+  +   H2
!        (6         33         45        2)
rate=xk(263)*y(6)*y(33)*xnH
r_f(263,6)=-rate
r_f(263,33)=-rate
r_f(263,45)=rate
r_f(263,2)=rate
! 264)   H3+   +   HCO   ->   H2CO+ +   H2
!        (6         35         48        2)
rate=xk(264)*y(6)*y(35)*xnH
r_f(264,6)=-rate
r_f(264,35)=-rate
r_f(264,48)=rate
r_f(264,2)=rate
! 265)   H3+   +   H2CO  ->   H3CO+ +   H2
!        (6         38         50        2)
rate=xk(265)*y(6)*y(38)*xnH
r_f(265,6)=-rate
r_f(265,38)=-rate
r_f(265,50)=rate
r_f(265,2)=rate
! 266)   H3+   +   CO2   ->   HCO2+ +   H2
!        (6         37         49        2)
rate=xk(266)*y(6)*y(37)*xnH
r_f(266,6)=-rate
r_f(266,37)=-rate
r_f(266,49)=rate
r_f(266,2)=rate
! 267)   He+   +   CH    ->   C+    +   H     +   He
!        (9         19         23        1         8)
rate=xk(267)*y(9)*y(19)*xnH
r_f(267,9)=-rate
r_f(267,19)=-rate
r_f(267,23)=rate
r_f(267,1)=rate
r_f(267,8)=rate
! 268)   He+   +   CH    ->   CH+   +   He
!        (9         19         25        8)
rate=xk(268)*y(9)*y(19)*xnH
r_f(268,9)=-rate
r_f(268,19)=-rate
r_f(268,25)=rate
r_f(268,8)=rate
! 269)   He+   +   CH2   ->   CH+   +   H     +   He
!        (9         20         25        1         8)
rate=xk(269)*y(9)*y(20)*xnH
r_f(269,9)=-rate
r_f(269,20)=-rate
r_f(269,25)=rate
r_f(269,1)=rate
r_f(269,8)=rate
! 270)   He+   +   CH2   ->   C+    +   H2    +   He
!        (9         20         23        2         8)
rate=xk(270)*y(9)*y(20)*xnH
r_f(270,9)=-rate
r_f(270,20)=-rate
r_f(270,23)=rate
r_f(270,2)=rate
r_f(270,8)=rate
! 271)   He+   +   CH3   ->   CH+   +   H2    +   He
!        (9         21         25        2         8)
rate=xk(271)*y(9)*y(21)*xnH
r_f(271,9)=-rate
r_f(271,21)=-rate
r_f(271,25)=rate
r_f(271,2)=rate
r_f(271,8)=rate
! 272)   He+   +   CH4   ->   CH3   +   H+    +   He
!        (9         22         21        4         8)
rate=xk(272)*y(9)*y(22)*xnH
r_f(272,9)=-rate
r_f(272,22)=-rate
r_f(272,21)=rate
r_f(272,4)=rate
r_f(272,8)=rate
! 273)   He+   +   CH4   ->   CH+   +   H2    +   He    +   H
!        (9         22         25        2         8         1)
rate=xk(273)*y(9)*y(22)*xnH
r_f(273,9)=-rate
r_f(273,22)=-rate
r_f(273,25)=rate
r_f(273,2)=rate
r_f(273,8)=rate
r_f(273,1)=rate
! 274)   He+   +   CH4   ->   CH2+  +   H2    +   He
!        (9         22         26        2         8)
rate=xk(274)*y(9)*y(22)*xnH
r_f(274,9)=-rate
r_f(274,22)=-rate
r_f(274,26)=rate
r_f(274,2)=rate
r_f(274,8)=rate
! 275)   He+   +   CH4   ->   CH3+  +   He    +   H
!        (9         22         27        8         1)
rate=xk(275)*y(9)*y(22)*xnH
r_f(275,9)=-rate
r_f(275,22)=-rate
r_f(275,27)=rate
r_f(275,8)=rate
r_f(275,1)=rate
! 276)   He+   +   CH4   ->   CH4+  +   He
!        (9         22         28        8)
rate=xk(276)*y(9)*y(22)*xnH
r_f(276,9)=-rate
r_f(276,22)=-rate
r_f(276,28)=rate
r_f(276,8)=rate
! 277)   He+   +   OH    ->   O+    +   H     +   He
!        (9         32         40        1         8)
rate=xk(277)*y(9)*y(32)*xnH
r_f(277,9)=-rate
r_f(277,32)=-rate
r_f(277,40)=rate
r_f(277,1)=rate
r_f(277,8)=rate
! 278)   He+   +   H2O   ->   H+    +   OH    +   He
!        (9         34         4         32        8)
rate=xk(278)*y(9)*y(34)*xnH
r_f(278,9)=-rate
r_f(278,34)=-rate
r_f(278,4)=rate
r_f(278,32)=rate
r_f(278,8)=rate
! 279)   He+   +   H2O   ->   OH+   +   H     +   He
!        (9         34         42        1         8)
rate=xk(279)*y(9)*y(34)*xnH
r_f(279,9)=-rate
r_f(279,34)=-rate
r_f(279,42)=rate
r_f(279,1)=rate
r_f(279,8)=rate
! 280)   He+   +   H2O   ->   H2O+  +   He
!        (9         34         44        8)
rate=xk(280)*y(9)*y(34)*xnH
r_f(280,9)=-rate
r_f(280,34)=-rate
r_f(280,44)=rate
r_f(280,8)=rate
! 281)   He+   +   C2    ->   C+    +   C    +   He
!        (9         18         23        17        8)
rate=xk(281)*y(9)*y(18)*xnH
r_f(281,9)=-rate
r_f(281,18)=-rate
r_f(281,23)=rate
r_f(281,17)=rate
r_f(281,8)=rate
! 282)   He+   +   C2    ->   C2+   +   He
!        (9         18         24        8)
rate=xk(282)*y(9)*y(18)*xnH
r_f(282,9)=-rate
r_f(282,18)=-rate
r_f(282,24)=rate
r_f(282,8)=rate
! 283)   He+   +   CO    ->   C+    +   O     +   He
!        (9         33         23        30        8)
rate=xk(283)*y(9)*y(33)*xnH
r_f(283,9)=-rate
r_f(283,33)=-rate
r_f(283,23)=rate
r_f(283,30)=rate
r_f(283,8)=rate
! 284)   He+   +   HCO   ->   CO+   +   H     +   He
!        (9         35         43        1         8)
rate=xk(284)*y(9)*y(35)*xnH
r_f(284,9)=-rate
r_f(284,35)=-rate
r_f(284,43)=rate
r_f(284,1)=rate
r_f(284,8)=rate
! 285)   He+   +   HCO   ->   CH+   +   O     +   He
!        (9         35         25        30        8)
rate=xk(285)*y(9)*y(35)*xnH
r_f(285,9)=-rate
r_f(285,35)=-rate
r_f(285,25)=rate
r_f(285,30)=rate
r_f(285,8)=rate
! 286)   He+   +   HCO   ->   HeH+  +   CO
!        (9         35         11        33)
rate=xk(286)*y(9)*y(35)*xnH
r_f(286,9)=-rate
r_f(286,35)=-rate
r_f(286,11)=rate
r_f(286,33)=rate
! 287)   He+   +   H2CO  ->   CO+   +   H2    +   He
!        (9         38         43        2         8)
rate=xk(287)*y(9)*y(38)*xnH
r_f(287,9)=-rate
r_f(287,38)=-rate
r_f(287,43)=rate
r_f(287,2)=rate
r_f(287,8)=rate
! 288)   He+   +   H2CO  ->   HCO+  +   H     +   He
!        (9         38         45        1         8)
rate=xk(288)*y(9)*y(38)*xnH
r_f(288,9)=-rate
r_f(288,38)=-rate
r_f(288,45)=rate
r_f(288,1)=rate
r_f(288,8)=rate
! 289)   He+   +   O2    ->   O+    +   O     +   He
!        (9         31         40        30        8)
rate=xk(289)*y(9)*y(31)*xnH
r_f(289,9)=-rate
r_f(289,31)=-rate
r_f(289,40)=rate
r_f(289,30)=rate
r_f(289,8)=rate
! 290)   He+   +   O2    ->   O2+   +   He
!        (9         31         41        8)
rate=xk(290)*y(9)*y(31)*xnH
r_f(290,9)=-rate
r_f(290,31)=-rate
r_f(290,41)=rate
r_f(290,8)=rate
! 291)   He+   +   CO2   ->   O2+   +   C    +   He
!        (9         37         41        17        8)
rate=xk(291)*y(9)*y(37)*xnH
r_f(291,9)=-rate
r_f(291,37)=-rate
r_f(291,41)=rate
r_f(291,17)=rate
r_f(291,8)=rate
! 292)   He+   +   CO2   ->   O+    +   CO    +   He
!        (9         37         40        33        8)
rate=xk(292)*y(9)*y(37)*xnH
r_f(292,9)=-rate
r_f(292,37)=-rate
r_f(292,40)=rate
r_f(292,33)=rate
r_f(292,8)=rate
! 293)   He+   +   CO2   ->   CO+   +   O     +   He
!        (9         37         43        30        8)
rate=xk(293)*y(9)*y(37)*xnH
r_f(293,9)=-rate
r_f(293,37)=-rate
r_f(293,43)=rate
r_f(293,30)=rate
r_f(293,8)=rate
! 294)   He+   +   CO2   ->   C+    +   O2    +   He
!        (9         37         23        31        8)
rate=xk(294)*y(9)*y(37)*xnH
r_f(294,9)=-rate
r_f(294,37)=-rate
r_f(294,23)=rate
r_f(294,31)=rate
r_f(294,8)=rate
! 295)   C+    +   H     ->   CH+   +   ph.
!        (23        1          25)
rate=xk(295)*y(23)*y(1)*xnH
r_f(295,23)=-rate
r_f(295,1)=-rate
r_f(295,25)=rate
! 296)   C+    +   O     ->   CO+   +   ph.
!        (23        30         43)
rate=xk(296)*y(23)*y(30)*xnH
r_f(296,23)=-rate
r_f(296,30)=-rate
r_f(296,43)=rate
! 297)   C+    +   H-    ->   H     +   C
!        (23        7          1         17)
rate=xk(297)*y(23)*y(7)*xnH
r_f(297,23)=-rate
r_f(297,7)=-rate
r_f(297,1)=rate
r_f(297,17)=rate
! 298)   C+    +   H2    ->   CH2+  +   ph.
!        (23        2          26)
rate=xk(298)*y(23)*y(2)*xnH
r_f(298,23)=-rate
r_f(298,2)=-rate
r_f(298,26)=rate
! 299)   C+    +   CH    ->   C2+   +   H
!        (23        19         24        1)
rate=xk(299)*y(23)*y(19)*xnH
r_f(299,23)=-rate
r_f(299,19)=-rate
r_f(299,24)=rate
r_f(299,1)=rate
! 300)   C+    +   CH    ->   CH+   +   C
!        (23        19         25        17)
rate=xk(300)*y(23)*y(19)*xnH
r_f(300,23)=-rate
r_f(300,19)=-rate
r_f(300,25)=rate
r_f(300,17)=rate
! 301)   C+    +   CH2   ->   CH2+  +   C
!        (23        20         26        17)
rate=xk(301)*y(23)*y(20)*xnH
r_f(301,23)=-rate
r_f(301,20)=-rate
r_f(301,26)=rate
r_f(301,17)=rate
! 302)   C+    +   OH    ->   CO+   +   H
!        (23        32         43        1)
rate=xk(302)*y(23)*y(32)*xnH
r_f(302,23)=-rate
r_f(302,32)=-rate
r_f(302,43)=rate
r_f(302,1)=rate
! 303)   C+    +   H2O   ->   HCO+  +   H
!        (23        34         45        1)
rate=xk(303)*y(23)*y(34)*xnH
r_f(303,23)=-rate
r_f(303,34)=-rate
r_f(303,45)=rate
r_f(303,1)=rate
! 304)   C+    +   HCO   ->   HCO+  +   C
!        (23        35         45        17)
rate=xk(304)*y(23)*y(35)*xnH
r_f(304,23)=-rate
r_f(304,35)=-rate
r_f(304,45)=rate
r_f(304,17)=rate
! 305)   C+    +   HCO   ->   CH+   +   CO
!        (23        35         25        33)
rate=xk(305)*y(23)*y(35)*xnH
r_f(305,23)=-rate
r_f(305,35)=-rate
r_f(305,25)=rate
r_f(305,33)=rate
! 306)   C+    +   H2CO  ->   CH2+  +   CO
!        (23        38         26        33)
rate=xk(306)*y(23)*y(38)*xnH
r_f(306,23)=-rate
r_f(306,38)=-rate
r_f(306,26)=rate
r_f(306,33)=rate
! 307)   C+    +   H2CO  ->   HCO+  +   CH
!        (23        38         45        19)
rate=xk(307)*y(23)*y(38)*xnH
r_f(307,23)=-rate
r_f(307,38)=-rate
r_f(307,45)=rate
r_f(307,19)=rate
! 308)   C+    +   H2CO  ->   H2CO+ +   C
!        (23        38         48        17)
rate=xk(308)*y(23)*y(38)*xnH
r_f(308,23)=-rate
r_f(308,38)=-rate
r_f(308,48)=rate
r_f(308,17)=rate
! 309)   C+    +   O2    ->   CO+   +   O
!        (23        31         43        30)
rate=xk(309)*y(23)*y(31)*xnH
r_f(309,23)=-rate
r_f(309,31)=-rate
r_f(309,43)=rate
r_f(309,30)=rate
! 310)   C+    +   O2    ->   O+    +   CO
!        (23        31         40        33)
rate=xk(310)*y(23)*y(31)*xnH
r_f(310,23)=-rate
r_f(310,31)=-rate
r_f(310,40)=rate
r_f(310,33)=rate
! 311)   C+    +   CO2   ->   CO+   +   CO
!        (23        37         43        33)
rate=xk(311)*y(23)*y(37)*xnH
r_f(311,23)=-rate
r_f(311,37)=-rate
r_f(311,43)=rate
r_f(311,33)=rate
! 312)   CH+   +   H     ->   C+    +   H2
!        (25        1          23        2)
rate=xk(312)*y(25)*y(1)*xnH
r_f(312,25)=-rate
r_f(312,1)=-rate
r_f(312,23)=rate
r_f(312,2)=rate
! 313)   CH+   +   C    ->   C2+   +   H
!        (25        17         24        1)
rate=xk(313)*y(25)*y(17)*xnH
r_f(313,25)=-rate
r_f(313,17)=-rate
r_f(313,24)=rate
r_f(313,1)=rate
! 314)   CH+   +   O     ->   CO+   +   H
!        (25        30         43        1)
rate=xk(314)*y(25)*y(30)*xnH
r_f(314,25)=-rate
r_f(314,30)=-rate
r_f(314,43)=rate
r_f(314,1)=rate
! 315)   CH+   +   H2    ->   CH2+  +   H
!        (25        2          26        1)
rate=xk(315)*y(25)*y(2)*xnH
r_f(315,25)=-rate
r_f(315,2)=-rate
r_f(315,26)=rate
r_f(315,1)=rate
! 316)   CH+   +   CH    ->   C2+   +   H2
!        (25        19         24        2)
rate=xk(316)*y(25)*y(19)*xnH
r_f(316,25)=-rate
r_f(316,19)=-rate
r_f(316,24)=rate
r_f(316,2)=rate
! 317)   CH+   +   OH    ->   CO+   +   H2
!        (25        32         43        2)
rate=xk(317)*y(25)*y(32)*xnH
r_f(317,25)=-rate
r_f(317,32)=-rate
r_f(317,43)=rate
r_f(317,2)=rate
! 318)   CH+   +   H2O   ->   HCO+  +   H2
!        (25        34         45        2)
rate=xk(318)*y(25)*y(34)*xnH
r_f(318,25)=-rate
r_f(318,34)=-rate
r_f(318,45)=rate
r_f(318,2)=rate
! 319)   CH+   +   H2O   ->   H3O+  +   C
!        (25        34         47        17)
rate=xk(319)*y(25)*y(34)*xnH
r_f(319,25)=-rate
r_f(319,34)=-rate
r_f(319,47)=rate
r_f(319,17)=rate
! 320)   CH+   +   H2O   ->   H2CO+ +   H
!        (25        34         48        1)
rate=xk(320)*y(25)*y(34)*xnH
r_f(320,25)=-rate
r_f(320,34)=-rate
r_f(320,48)=rate
r_f(320,1)=rate
! 321)   CH+   +   CO    ->   HCO+  +   C
!        (25        33         45        17)
rate=xk(321)*y(25)*y(33)*xnH
r_f(321,25)=-rate
r_f(321,33)=-rate
r_f(321,45)=rate
r_f(321,17)=rate
! 322)   CH+   +   HCO   ->   CH2+  +   CO
!        (25        35         26        33)
rate=xk(322)*y(25)*y(35)*xnH
r_f(322,25)=-rate
r_f(322,35)=-rate
r_f(322,26)=rate
r_f(322,33)=rate
! 323)   CH+   +   HCO   ->   HCO+  +   CH
!        (25        35         45        19)
rate=xk(323)*y(25)*y(35)*xnH
r_f(323,25)=-rate
r_f(323,35)=-rate
r_f(323,45)=rate
r_f(323,19)=rate
! 324)   CH+   +   H2CO  ->   CH3+  +   CO
!        (25        38         27        33)
rate=xk(324)*y(25)*y(38)*xnH
r_f(324,25)=-rate
r_f(324,38)=-rate
r_f(324,27)=rate
r_f(324,33)=rate
! 325)   CH+   +   H2CO  ->   HCO+  +   CH2
!        (25        38         45        20)
rate=xk(325)*y(25)*y(38)*xnH
r_f(325,25)=-rate
r_f(325,38)=-rate
r_f(325,45)=rate
r_f(325,20)=rate
! 326)   CH+   +   H2CO  ->   H3CO+ +   C
!        (25        38         50        17)
rate=xk(326)*y(25)*y(38)*xnH
r_f(326,25)=-rate
r_f(326,38)=-rate
r_f(326,50)=rate
r_f(326,17)=rate
! 327)   CH+   +   O2    ->   HCO+  +   O
!        (25        31         45        30)
rate=xk(327)*y(25)*y(31)*xnH
r_f(327,25)=-rate
r_f(327,31)=-rate
r_f(327,45)=rate
r_f(327,30)=rate
! 328)   CH+   +   O2    ->   HCO   +   O+
!        (25        31         35        40)
rate=xk(328)*y(25)*y(31)*xnH
r_f(328,25)=-rate
r_f(328,31)=-rate
r_f(328,35)=rate
r_f(328,40)=rate
! 329)   CH+   +   O2    ->   CO+   +   OH
!        (25        31         43        32)
rate=xk(329)*y(25)*y(31)*xnH
r_f(329,25)=-rate
r_f(329,31)=-rate
r_f(329,43)=rate
r_f(329,32)=rate
! 330)   CH+   +   CO2   ->   HCO+  +   CO
!        (25        37         45        33)
rate=xk(330)*y(25)*y(37)*xnH
r_f(330,25)=-rate
r_f(330,37)=-rate
r_f(330,45)=rate
r_f(330,33)=rate
! 331)   CH2+  +   O     ->   HCO+  +   H
!        (26        30         45        1)
rate=xk(331)*y(26)*y(30)*xnH
r_f(331,26)=-rate
r_f(331,30)=-rate
r_f(331,45)=rate
r_f(331,1)=rate
! 332)   CH2+  +   H2    ->   CH3+  +   H
!        (26        2          27        1)
rate=xk(332)*y(26)*y(2)*xnH
r_f(332,26)=-rate
r_f(332,2)=-rate
r_f(332,27)=rate
r_f(332,1)=rate
! 333)   CH2+  +   H2O   ->   H3CO+ +   H
!        (26        34         50        1)
rate=xk(333)*y(26)*y(34)*xnH
r_f(333,26)=-rate
r_f(333,34)=-rate
r_f(333,50)=rate
r_f(333,1)=rate
! 334)   CH2+  +   HCO   ->   CH3+  +   CO
!        (26        35         27        33)
rate=xk(334)*y(26)*y(35)*xnH
r_f(334,26)=-rate
r_f(334,35)=-rate
r_f(334,27)=rate
r_f(334,33)=rate
! 335)   CH2+  +   H2CO  ->   HCO+  +   CH3
!        (26        38         45        21)
rate=xk(335)*y(26)*y(38)*xnH
r_f(335,26)=-rate
r_f(335,38)=-rate
r_f(335,45)=rate
r_f(335,21)=rate
! 336)   CH2+  +   O2    ->   HCO+  +   OH
!        (26        31         45        32)
rate=xk(336)*y(26)*y(31)*xnH
r_f(336,26)=-rate
r_f(336,31)=-rate
r_f(336,45)=rate
r_f(336,32)=rate
! 337)   CH2+  +   CO2   ->   H2CO+ +   CO
!        (26        37         48        33)
rate=xk(337)*y(26)*y(37)*xnH
r_f(337,26)=-rate
r_f(337,37)=-rate
r_f(337,48)=rate
r_f(337,33)=rate
! 338)   CH3+  +   O     ->   HCO+  +   H2
!        (27        30         45        2)
rate=xk(338)*y(27)*y(30)*xnH
r_f(338,27)=-rate
r_f(338,30)=-rate
r_f(338,45)=rate
r_f(338,2)=rate
! 339)   CH3+  +   O     ->   H2CO+ +   H
!        (27        30         48        1)
rate=xk(339)*y(27)*y(30)*xnH
r_f(339,27)=-rate
r_f(339,30)=-rate
r_f(339,48)=rate
r_f(339,1)=rate
! 340)   CH3+  +   H2    ->   CH5+  +   ph.
!        (27        2          29)
rate=xk(340)*y(27)*y(2)*xnH
r_f(340,27)=-rate
r_f(340,2)=-rate
r_f(340,29)=rate
! 341)   CH3+  +   OH    ->   H2CO+ +   H2
!        (27        32         48        2)
rate=xk(341)*y(27)*y(32)*xnH
r_f(341,27)=-rate
r_f(341,32)=-rate
r_f(341,48)=rate
r_f(341,2)=rate
! 342)   CH3+  +   HCO   ->   CH4+  +   CO
!        (27        35         28        33)
rate=xk(342)*y(27)*y(35)*xnH
r_f(342,27)=-rate
r_f(342,35)=-rate
r_f(342,28)=rate
r_f(342,33)=rate
! 343)   CH3+  +   HCO   ->   HCO+  +   CH3
!        (27        35         45        21)
rate=xk(343)*y(27)*y(35)*xnH
r_f(343,27)=-rate
r_f(343,35)=-rate
r_f(343,45)=rate
r_f(343,21)=rate
! 344)   CH3+  +   H2CO  ->   HCO+  +   CH4
!        (27        38         45        22)
rate=xk(344)*y(27)*y(38)*xnH
r_f(344,27)=-rate
r_f(344,38)=-rate
r_f(344,45)=rate
r_f(344,22)=rate
! 345)   CH3+  +   O2    ->   H3CO+ +   O
!        (27        31         50        30)
rate=xk(345)*y(27)*y(31)*xnH
r_f(345,27)=-rate
r_f(345,31)=-rate
r_f(345,50)=rate
r_f(345,30)=rate
! 346)   O+    +   H     ->   H+    +   O
!        (40        1          4         30)
rate=xk(346)*y(40)*y(1)*xnH
r_f(346,40)=-rate
r_f(346,1)=-rate
r_f(346,4)=rate
r_f(346,30)=rate
! 347)   O+    +   H-    ->   H     +   O
!        (40        7          1         30)
rate=xk(347)*y(40)*y(7)*xnH
r_f(347,40)=-rate
r_f(347,7)=-rate
r_f(347,1)=rate
r_f(347,30)=rate
! 348)   O+    +   H2    ->   OH+   +   H
!        (40        2          42        1)
rate=xk(348)*y(40)*y(2)*xnH
r_f(348,40)=-rate
r_f(348,2)=-rate
r_f(348,42)=rate
r_f(348,1)=rate
! 349)   O+    +   CH    ->   CH+   +   O
!        (40        19         25        30)
rate=xk(349)*y(40)*y(19)*xnH
r_f(349,40)=-rate
r_f(349,19)=-rate
r_f(349,25)=rate
r_f(349,30)=rate
! 350)   O+    +   CH    ->   CO+   +   H
!        (40        19         43        1)
rate=xk(350)*y(40)*y(19)*xnH
r_f(350,40)=-rate
r_f(350,19)=-rate
r_f(350,43)=rate
r_f(350,1)=rate
! 351)   O+    +   CH2   ->   CH2+  +   O
!        (40        20         26        30)
rate=xk(351)*y(40)*y(20)*xnH
r_f(351,40)=-rate
r_f(351,20)=-rate
r_f(351,26)=rate
r_f(351,30)=rate
! 352)   O+    +   CH4   ->   CH3+  +   OH
!        (40        22         27        32)
rate=xk(352)*y(40)*y(22)*xnH
r_f(352,40)=-rate
r_f(352,22)=-rate
r_f(352,27)=rate
r_f(352,32)=rate
! 353)   O+    +   CH4   ->   CH4+  +   O
!        (40        22         28        30)
rate=xk(353)*y(40)*y(22)*xnH
r_f(353,40)=-rate
r_f(353,22)=-rate
r_f(353,28)=rate
r_f(353,30)=rate
! 354)   O+    +   OH    ->   OH+   +   O
!        (40        32         42        30)
rate=xk(354)*y(40)*y(32)*xnH
r_f(354,40)=-rate
r_f(354,32)=-rate
r_f(354,42)=rate
r_f(354,30)=rate
! 355)   O+    +   OH    ->   O2+   +   H
!        (40        32         41        1)
rate=xk(355)*y(40)*y(32)*xnH
r_f(355,40)=-rate
r_f(355,32)=-rate
r_f(355,41)=rate
r_f(355,1)=rate
! 356)   O+    +   H2O   ->   H2O+  +   O
!        (40        34         44        30)
rate=xk(356)*y(40)*y(34)*xnH
r_f(356,40)=-rate
r_f(356,34)=-rate
r_f(356,44)=rate
r_f(356,30)=rate
! 357)   O+    +   C2    ->   C2+   +   O
!        (40        18         24        30)
rate=xk(357)*y(40)*y(18)*xnH
r_f(357,40)=-rate
r_f(357,18)=-rate
r_f(357,24)=rate
r_f(357,30)=rate
! 358)   O+    +   C2    ->   CO+   +   C
!        (40        18         43        17)
rate=xk(358)*y(40)*y(18)*xnH
r_f(358,40)=-rate
r_f(358,18)=-rate
r_f(358,43)=rate
r_f(358,17)=rate
! 359)   O+    +   HCO   ->   OH+   +   CO
!        (40        35         42        33)
rate=xk(359)*y(40)*y(35)*xnH
r_f(359,40)=-rate
r_f(359,35)=-rate
r_f(359,42)=rate
r_f(359,33)=rate
! 360)   O+    +   HCO   ->   HCO+  +   O
!        (40        35         45        30)
rate=xk(360)*y(40)*y(35)*xnH
r_f(360,40)=-rate
r_f(360,35)=-rate
r_f(360,45)=rate
r_f(360,30)=rate
! 361)   O+    +   H2CO  ->   HCO+  +   OH
!        (40        38         45        32)
rate=xk(361)*y(40)*y(38)*xnH
r_f(361,40)=-rate
r_f(361,38)=-rate
r_f(361,45)=rate
r_f(361,32)=rate
! 362)   O+    +   H2CO  ->   H2CO+ +   O
!        (40        38         48        30)
rate=xk(362)*y(40)*y(38)*xnH
r_f(362,40)=-rate
r_f(362,38)=-rate
r_f(362,48)=rate
r_f(362,30)=rate
! 363)   O+    +   O2    ->   O2+   +   O
!        (40        31         41        30)
rate=xk(363)*y(40)*y(31)*xnH
r_f(363,40)=-rate
r_f(363,31)=-rate
r_f(363,41)=rate
r_f(363,30)=rate
! 364)   O+    +   CO2   ->   O2+   +   CO
!        (40        37         41        33)
rate=xk(364)*y(40)*y(37)*xnH
r_f(364,40)=-rate
r_f(364,37)=-rate
r_f(364,41)=rate
r_f(364,33)=rate
! 365)   CH4+  +   H     ->   CH3+  +   H2
!        (28        1          27        2)
rate=xk(365)*y(28)*y(1)*xnH
r_f(365,28)=-rate
r_f(365,1)=-rate
r_f(365,27)=rate
r_f(365,2)=rate
! 366)   CH4+  +   O     ->   CH3+  +   OH
!        (28        30         27        32)
rate=xk(366)*y(28)*y(30)*xnH
r_f(366,28)=-rate
r_f(366,30)=-rate
r_f(366,27)=rate
r_f(366,32)=rate
! 367)   CH4+  +   H2    ->   CH5+  +   H
!        (28        2          29        1)
rate=xk(367)*y(28)*y(2)*xnH
r_f(367,28)=-rate
r_f(367,2)=-rate
r_f(367,29)=rate
r_f(367,1)=rate
! 368)   CH4+  +   CH4   ->   CH5+  +   CH3
!        (28        22         29        21)
rate=xk(368)*y(28)*y(22)*xnH
r_f(368,28)=-rate
r_f(368,22)=-rate
r_f(368,29)=rate
r_f(368,21)=rate
! 369)   CH4+  +   H2O   ->   H3O+  +   CH3
!        (28        34         47        21)
rate=xk(369)*y(28)*y(34)*xnH
r_f(369,28)=-rate
r_f(369,34)=-rate
r_f(369,47)=rate
r_f(369,21)=rate
! 370)   CH4+  +   CO    ->   HCO+  +   CH3
!        (28        33         45        21)
rate=xk(370)*y(28)*y(33)*xnH
r_f(370,28)=-rate
r_f(370,33)=-rate
r_f(370,45)=rate
r_f(370,21)=rate
! 371)   CH4+  +   H2CO  ->   H2CO+ +   CH4
!        (28        38         48        22)
rate=xk(371)*y(28)*y(38)*xnH
r_f(371,28)=-rate
r_f(371,38)=-rate
r_f(371,48)=rate
r_f(371,22)=rate
! 372)   CH4+  +   H2CO  ->   H3CO+ +   CH3
!        (28        38         50        21)
rate=xk(372)*y(28)*y(38)*xnH
r_f(372,28)=-rate
r_f(372,38)=-rate
r_f(372,50)=rate
r_f(372,21)=rate
! 373)   CH4+  +   O2    ->   O2+   +   CH4
!        (28        31         41        22)
rate=xk(373)*y(28)*y(31)*xnH
r_f(373,28)=-rate
r_f(373,31)=-rate
r_f(373,41)=rate
r_f(373,22)=rate
! 374)   CH4+  +   CO2   ->   HCO2+ +   CH3
!        (28        37         49        21)
rate=xk(374)*y(28)*y(37)*xnH
r_f(374,28)=-rate
r_f(374,37)=-rate
r_f(374,49)=rate
r_f(374,21)=rate
! 375)   OH+   +   C    ->   CH+   +   O
!        (42        17         25        30)
rate=xk(375)*y(42)*y(17)*xnH
r_f(375,42)=-rate
r_f(375,17)=-rate
r_f(375,25)=rate
r_f(375,30)=rate
! 376)   OH+   +   O     ->   O2+   +   H
!        (42        30         41        1)
rate=xk(376)*y(42)*y(30)*xnH
r_f(376,42)=-rate
r_f(376,30)=-rate
r_f(376,41)=rate
r_f(376,1)=rate
! 377)   OH+   +   H2    ->   H2O+  +   H
!        (42        2          44        1)
rate=xk(377)*y(42)*y(2)*xnH
r_f(377,42)=-rate
r_f(377,2)=-rate
r_f(377,44)=rate
r_f(377,1)=rate
! 378)   OH+   +   CH    ->   CH+   +   OH
!        (42        19         25        32)
rate=xk(378)*y(42)*y(19)*xnH
r_f(378,42)=-rate
r_f(378,19)=-rate
r_f(378,25)=rate
r_f(378,32)=rate
! 379)   OH+   +   CH    ->   CH2+  +   O
!        (42        19         26        30)
rate=xk(379)*y(42)*y(19)*xnH
r_f(379,42)=-rate
r_f(379,19)=-rate
r_f(379,26)=rate
r_f(379,30)=rate
! 380)   OH+   +   CH2   ->   CH2+  +   OH
!        (42        20         26        32)
rate=xk(380)*y(42)*y(20)*xnH
r_f(380,42)=-rate
r_f(380,20)=-rate
r_f(380,26)=rate
r_f(380,32)=rate
! 381)   OH+   +   CH2   ->   CH3+  +   O
!        (42        20         27        30)
rate=xk(381)*y(42)*y(20)*xnH
r_f(381,42)=-rate
r_f(381,20)=-rate
r_f(381,27)=rate
r_f(381,30)=rate
! 382)   OH+   +   CH4   ->   H3O+  +   CH2
!        (42        22         47        20)
rate=xk(382)*y(42)*y(22)*xnH
r_f(382,42)=-rate
r_f(382,22)=-rate
r_f(382,47)=rate
r_f(382,20)=rate
! 383)   OH+   +   CH4   ->   CH5+  +   O
!        (42        22         29        30)
rate=xk(383)*y(42)*y(22)*xnH
r_f(383,42)=-rate
r_f(383,22)=-rate
r_f(383,29)=rate
r_f(383,30)=rate
! 384)   OH+   +   OH    ->   H2O+  +   O
!        (42        32         44        30)
rate=xk(384)*y(42)*y(32)*xnH
r_f(384,42)=-rate
r_f(384,32)=-rate
r_f(384,44)=rate
r_f(384,30)=rate
! 385)   OH+   +   H2O   ->   H2O+  +   OH
!        (42        34         44        32)
rate=xk(385)*y(42)*y(34)*xnH
r_f(385,42)=-rate
r_f(385,34)=-rate
r_f(385,44)=rate
r_f(385,32)=rate
! 386)   OH+   +   H2O   ->   H3O+  +   O
!        (42        34         47        30)
rate=xk(386)*y(42)*y(34)*xnH
r_f(386,42)=-rate
r_f(386,34)=-rate
r_f(386,47)=rate
r_f(386,30)=rate
! 387)   OH+   +   C2    ->   C2+   +   OH
!        (42        18         24        32)
rate=xk(387)*y(42)*y(18)*xnH
r_f(387,42)=-rate
r_f(387,18)=-rate
r_f(387,24)=rate
r_f(387,32)=rate
! 388)   OH+   +   CO    ->   HCO+  +   O
!        (42        33         45        30)
rate=xk(388)*y(42)*y(33)*xnH
r_f(388,42)=-rate
r_f(388,33)=-rate
r_f(388,45)=rate
r_f(388,30)=rate
! 389)   OH+   +   HCO   ->   H2O+  +   CO
!        (42        35         44        33)
rate=xk(389)*y(42)*y(35)*xnH
r_f(389,42)=-rate
r_f(389,35)=-rate
r_f(389,44)=rate
r_f(389,33)=rate
! 390)   OH+   +   HCO   ->   HCO+  +   OH
!        (42        35         45        32)
rate=xk(390)*y(42)*y(35)*xnH
r_f(390,42)=-rate
r_f(390,35)=-rate
r_f(390,45)=rate
r_f(390,32)=rate
! 391)   OH+   +   HCO   ->   H2CO+ +   O
!        (42        35         48        30)
rate=xk(391)*y(42)*y(35)*xnH
r_f(391,42)=-rate
r_f(391,35)=-rate
r_f(391,48)=rate
r_f(391,30)=rate
! 392)   OH+   +   H2CO  ->   H3CO+ +   O
!        (42        38         50        30)
rate=xk(392)*y(42)*y(38)*xnH
r_f(392,42)=-rate
r_f(392,38)=-rate
r_f(392,50)=rate
r_f(392,30)=rate
! 393)   OH+   +   H2CO  ->   H2CO+ +   OH
!        (42        38         48        32)
rate=xk(393)*y(42)*y(38)*xnH
r_f(393,42)=-rate
r_f(393,38)=-rate
r_f(393,48)=rate
r_f(393,32)=rate
! 394)   OH+   +   O2    ->   O2+   +   OH
!        (42        31         41        32)
rate=xk(394)*y(42)*y(31)*xnH
r_f(394,42)=-rate
r_f(394,31)=-rate
r_f(394,41)=rate
r_f(394,32)=rate
! 395)   OH+   +   CO2   ->   HCO2+ +   O
!        (42        37         49        30)
rate=xk(395)*y(42)*y(37)*xnH
r_f(395,42)=-rate
r_f(395,37)=-rate
r_f(395,49)=rate
r_f(395,30)=rate
! 396)   CH5+  +   H     ->   CH4+  +   H2
!        (29        1          28        2)
rate=xk(396)*y(29)*y(1)*xnH
r_f(396,29)=-rate
r_f(396,1)=-rate
r_f(396,28)=rate
r_f(396,2)=rate
! 397)   CH5+  +   C     ->   CH+   +   CH4
!        (29        17         25        22)
rate=xk(397)*y(29)*y(17)*xnH
r_f(397,29)=-rate
r_f(397,17)=-rate
r_f(397,25)=rate
r_f(397,22)=rate
! 398)   CH5+  +   O     ->   H3O+  +   CH2
!        (29        30         47        20)
rate=xk(398)*y(29)*y(30)*xnH
r_f(398,29)=-rate
r_f(398,30)=-rate
r_f(398,47)=rate
r_f(398,20)=rate
! 399)   CH5+  +   O     ->   H3CO+ +   H2
!        (29        30         50        2)
rate=xk(399)*y(29)*y(30)*xnH
r_f(399,29)=-rate
r_f(399,30)=-rate
r_f(399,50)=rate
r_f(399,2)=rate
! 400)   CH5+  +   CH    ->   CH2+  +   CH4
!        (29        19         26        22)
rate=xk(400)*y(29)*y(19)*xnH
r_f(400,29)=-rate
r_f(400,19)=-rate
r_f(400,26)=rate
r_f(400,22)=rate
! 401)   CH5+  +   CH2   ->   CH3+  +   CH4
!        (29        20         27        22)
rate=xk(401)*y(29)*y(20)*xnH
r_f(401,29)=-rate
r_f(401,20)=-rate
r_f(401,27)=rate
r_f(401,22)=rate
! 402)   CH5+  +   OH    ->   H2O+  +   CH4
!        (29        32         44        22)
rate=xk(402)*y(29)*y(32)*xnH
r_f(402,29)=-rate
r_f(402,32)=-rate
r_f(402,44)=rate
r_f(402,22)=rate
! 403)   CH5+  +   H2O   ->   H3O+  +   CH4
!        (29        34         47        22)
rate=xk(403)*y(29)*y(34)*xnH
r_f(403,29)=-rate
r_f(403,34)=-rate
r_f(403,47)=rate
r_f(403,22)=rate
! 404)   CH5+  +   CO    ->   HCO+  +   CH4
!        (29        33         45        22)
rate=xk(404)*y(29)*y(33)*xnH
r_f(404,29)=-rate
r_f(404,33)=-rate
r_f(404,45)=rate
r_f(404,22)=rate
! 405)   CH5+  +   HCO   ->   H2CO+ +   CH4
!        (29        35         48        22)
rate=xk(405)*y(29)*y(35)*xnH
r_f(405,29)=-rate
r_f(405,35)=-rate
r_f(405,48)=rate
r_f(405,22)=rate
! 406)   CH5+  +   H2CO  ->   H3CO+ +   CH4
!        (29        38         50        22)
rate=xk(406)*y(29)*y(38)*xnH
r_f(406,29)=-rate
r_f(406,38)=-rate
r_f(406,50)=rate
r_f(406,22)=rate
! 407)   CH5+  +   CO2   ->   HCO2+ +   CH4
!        (29        37         49        22)
rate=xk(407)*y(29)*y(37)*xnH
r_f(407,29)=-rate
r_f(407,37)=-rate
r_f(407,49)=rate
r_f(407,22)=rate
! 408)   H2O+  +   C    ->   CH+   +   OH
!        (44        17         25        32)
rate=xk(408)*y(44)*y(17)*xnH
r_f(408,44)=-rate
r_f(408,17)=-rate
r_f(408,25)=rate
r_f(408,32)=rate
! 409)   H2O+  +   O     ->   O2+   +   H2
!        (44        30         41        2)
rate=xk(409)*y(44)*y(30)*xnH
r_f(409,44)=-rate
r_f(409,30)=-rate
r_f(409,41)=rate
r_f(409,2)=rate
! 410)   H2O+  +   H2    ->   H3O+  +   H
!        (44        2          47        1)
rate=xk(410)*y(44)*y(2)*xnH
r_f(410,44)=-rate
r_f(410,2)=-rate
r_f(410,47)=rate
r_f(410,1)=rate
! 411)   H2O+  +   CH    ->   CH+   +   H2O
!        (44        19         25        34)
rate=xk(411)*y(44)*y(19)*xnH
r_f(411,44)=-rate
r_f(411,19)=-rate
r_f(411,25)=rate
r_f(411,34)=rate
! 412)   H2O+  +   CH    ->   CH2+  +   OH
!        (44        19         26        32)
rate=xk(412)*y(44)*y(19)*xnH
r_f(412,44)=-rate
r_f(412,19)=-rate
r_f(412,26)=rate
r_f(412,32)=rate
! 413)   H2O+  +   CH2   ->   CH3+  +   OH
!        (44        20         27        32)
rate=xk(413)*y(44)*y(20)*xnH
r_f(413,44)=-rate
r_f(413,20)=-rate
r_f(413,27)=rate
r_f(413,32)=rate
! 414)   H2O+  +   CH2   ->   CH2+  +   H2O
!        (44        20         26        34)
rate=xk(414)*y(44)*y(20)*xnH
r_f(414,44)=-rate
r_f(414,20)=-rate
r_f(414,26)=rate
r_f(414,34)=rate
! 415)   H2O+  +   CH4   ->   H3O+  +   CH3
!        (44        22         47        21)
rate=xk(415)*y(44)*y(22)*xnH
r_f(415,44)=-rate
r_f(415,22)=-rate
r_f(415,47)=rate
r_f(415,21)=rate
! 416)   H2O+  +   OH    ->   H3O+  +   O
!        (44        32         47        30)
rate=xk(416)*y(44)*y(32)*xnH
r_f(416,44)=-rate
r_f(416,32)=-rate
r_f(416,47)=rate
r_f(416,30)=rate
! 417)   H2O+  +   H2O   ->   H3O+  +   OH
!        (44        34         47        32)
rate=xk(417)*y(44)*y(34)*xnH
r_f(417,44)=-rate
r_f(417,34)=-rate
r_f(417,47)=rate
r_f(417,32)=rate
! 418)   H2O+  +   C2    ->   C2+   +   H2O
!        (44        18         24        34)
rate=xk(418)*y(44)*y(18)*xnH
r_f(418,44)=-rate
r_f(418,18)=-rate
r_f(418,24)=rate
r_f(418,34)=rate
! 419)   H2O+  +   CO    ->   HCO+  +   OH
!        (44        33         45        32)
rate=xk(419)*y(44)*y(33)*xnH
r_f(419,44)=-rate
r_f(419,33)=-rate
r_f(419,45)=rate
r_f(419,32)=rate
! 420)   H2O+  +   HCO   ->   H3O+  +   CO
!        (44        35         47        33)
rate=xk(420)*y(44)*y(35)*xnH
r_f(420,44)=-rate
r_f(420,35)=-rate
r_f(420,47)=rate
r_f(420,33)=rate
! 421)   H2O+  +   HCO   ->   HCO+  +   H2O
!        (44        35         45        34)
rate=xk(421)*y(44)*y(35)*xnH
r_f(421,44)=-rate
r_f(421,35)=-rate
r_f(421,45)=rate
r_f(421,34)=rate
! 422)   H2O+  +   HCO   ->   H2CO+ +   OH
!        (44        35         48        32)
rate=xk(422)*y(44)*y(35)*xnH
r_f(422,44)=-rate
r_f(422,35)=-rate
r_f(422,48)=rate
r_f(422,32)=rate
! 423)   H2O+  +   H2CO  ->   H2CO+ +   H2O
!        (44        38         48        34)
rate=xk(423)*y(44)*y(38)*xnH
r_f(423,44)=-rate
r_f(423,38)=-rate
r_f(423,48)=rate
r_f(423,34)=rate
! 424)   H2O+  +   H2CO  ->   H3CO+ +   OH
!        (44        38         50        32)
rate=xk(424)*y(44)*y(38)*xnH
r_f(424,44)=-rate
r_f(424,38)=-rate
r_f(424,50)=rate
r_f(424,32)=rate
! 425)   H2O+  +   O2    ->   O2+   +   H2O
!        (44        31         41        34)
rate=xk(425)*y(44)*y(31)*xnH
r_f(425,44)=-rate
r_f(425,31)=-rate
r_f(425,41)=rate
r_f(425,34)=rate
! 426)   H3O+  +   C     ->   HCO+  +   H2
!        (47        17         45        2)
rate=xk(426)*y(47)*y(17)*xnH
r_f(426,47)=-rate
r_f(426,17)=-rate
r_f(426,45)=rate
r_f(426,2)=rate
! 427)   H3O+  +   H-    ->   OH    +   H2    +   H
!        (47        7          32        2         1)
rate=xk(427)*y(47)*y(7)*xnH
r_f(427,47)=-rate
r_f(427,7)=-rate
r_f(427,32)=rate
r_f(427,2)=rate
r_f(427,1)=rate
! 428)   H3O+  +   H-    ->   H2O   +   H2
!        (47        7          34        2)
rate=xk(428)*y(47)*y(7)*xnH
r_f(428,47)=-rate
r_f(428,7)=-rate
r_f(428,34)=rate
r_f(428,2)=rate
! 429)   H3O+  +   CH    ->   CH2+  +   H2O
!        (47        19         26        34)
rate=xk(429)*y(47)*y(19)*xnH
r_f(429,47)=-rate
r_f(429,19)=-rate
r_f(429,26)=rate
r_f(429,34)=rate
! 430)   H3O+  +   CH2   ->   CH3+  +   H2O
!        (47        20         27        34)
rate=xk(430)*y(47)*y(20)*xnH
r_f(430,47)=-rate
r_f(430,20)=-rate
r_f(430,27)=rate
r_f(430,34)=rate
! 431)   H3O+  +   H2CO  ->   H3CO+ +   H2O
!        (47        38         50        34)
rate=xk(431)*y(47)*y(38)*xnH
r_f(431,47)=-rate
r_f(431,38)=-rate
r_f(431,50)=rate
r_f(431,34)=rate
! 432)   C2+   +   C     ->   C+    +   C2
!        (24        17         23        18)
rate=xk(432)*y(24)*y(17)*xnH
r_f(432,24)=-rate
r_f(432,17)=-rate
r_f(432,23)=rate
r_f(432,18)=rate
! 433)   C2+   +   O     ->   CO+   +   C
!        (24        30         43        17)
rate=xk(433)*y(24)*y(30)*xnH
r_f(433,24)=-rate
r_f(433,30)=-rate
r_f(433,43)=rate
r_f(433,17)=rate
! 434)   C2+   +   CH    ->   CH+   +   C2
!        (24        19         25        18)
rate=xk(434)*y(24)*y(19)*xnH
r_f(434,24)=-rate
r_f(434,19)=-rate
r_f(434,25)=rate
r_f(434,18)=rate
! 435)   C2+   +   CH2   ->   CH2+  +   C2
!        (24        20         26        18)
rate=xk(435)*y(24)*y(20)*xnH
r_f(435,24)=-rate
r_f(435,20)=-rate
r_f(435,26)=rate
r_f(435,18)=rate
! 436)   C2+   +   OH    ->   OH+   +   C2
!        (24        32         42        18)
rate=xk(436)*y(24)*y(32)*xnH
r_f(436,24)=-rate
r_f(436,32)=-rate
r_f(436,42)=rate
r_f(436,18)=rate
! 437)   C2+   +   HCO   ->   HCO+  +   C2
!        (24        35         45        18)
rate=xk(437)*y(24)*y(35)*xnH
r_f(437,24)=-rate
r_f(437,35)=-rate
r_f(437,45)=rate
r_f(437,18)=rate
! 438)   C2+   +   O2    ->   CO+   +   CO
!        (24        31         43        33)
rate=xk(438)*y(24)*y(31)*xnH
r_f(438,24)=-rate
r_f(438,31)=-rate
r_f(438,43)=rate
r_f(438,33)=rate
! 439)   CO+   +   H     ->   H+    +   CO
!        (43        1          4         33)
rate=xk(439)*y(43)*y(1)*xnH
r_f(439,43)=-rate
r_f(439,1)=-rate
r_f(439,4)=rate
r_f(439,33)=rate
! 440)   CO+   +   C    ->   C+    +   CO
!        (43        17         23        33)
rate=xk(440)*y(43)*y(17)*xnH
r_f(440,43)=-rate
r_f(440,17)=-rate
r_f(440,23)=rate
r_f(440,33)=rate
! 441)   CO+   +   O     ->   O+    +   CO
!        (43        30         40        33)
rate=xk(441)*y(43)*y(30)*xnH
r_f(441,43)=-rate
r_f(441,30)=-rate
r_f(441,40)=rate
r_f(441,33)=rate
! 442)   CO+   +   H2    ->   HCO+  +   H
!        (43        2          45        1)
rate=xk(442)*y(43)*y(2)*xnH
r_f(442,43)=-rate
r_f(442,2)=-rate
r_f(442,45)=rate
r_f(442,1)=rate
! 443)   CO+   +   CH    ->   CH+   +   CO
!        (43        19         25        33)
rate=xk(443)*y(43)*y(19)*xnH
r_f(443,43)=-rate
r_f(443,19)=-rate
r_f(443,25)=rate
r_f(443,33)=rate
! 444)   CO+   +   CH    ->   HCO+  +   C
!        (43        19         45        17)
rate=xk(444)*y(43)*y(19)*xnH
r_f(444,43)=-rate
r_f(444,19)=-rate
r_f(444,45)=rate
r_f(444,17)=rate
! 445)   CO+   +   CH2   ->   CH2+  +   CO
!        (43        20         26        33)
rate=xk(445)*y(43)*y(20)*xnH
r_f(445,43)=-rate
r_f(445,20)=-rate
r_f(445,26)=rate
r_f(445,33)=rate
! 446)   CO+   +   CH2   ->   HCO+  +   CH
!        (43        20         45        19)
rate=xk(446)*y(43)*y(20)*xnH
r_f(446,43)=-rate
r_f(446,20)=-rate
r_f(446,45)=rate
r_f(446,19)=rate
! 447)   CO+   +   CH4   ->   CH4+  +   CO
!        (43        22         28        33)
rate=xk(447)*y(43)*y(22)*xnH
r_f(447,43)=-rate
r_f(447,22)=-rate
r_f(447,28)=rate
r_f(447,33)=rate
! 448)   CO+   +   CH4   ->   HCO+  +   CH3
!        (43        22         45        21)
rate=xk(448)*y(43)*y(22)*xnH
r_f(448,43)=-rate
r_f(448,22)=-rate
r_f(448,45)=rate
r_f(448,21)=rate
! 449)   CO+   +   OH    ->   OH+   +   CO
!        (43        32         42        33)
rate=xk(449)*y(43)*y(32)*xnH
r_f(449,43)=-rate
r_f(449,32)=-rate
r_f(449,42)=rate
r_f(449,33)=rate
! 450)   CO+   +   OH    ->   HCO+  +   O
!        (43        32         45        30)
rate=xk(450)*y(43)*y(32)*xnH
r_f(450,43)=-rate
r_f(450,32)=-rate
r_f(450,45)=rate
r_f(450,30)=rate
! 451)   CO+   +   H2O   ->   H2O+  +   CO
!        (43        34         44        33)
rate=xk(451)*y(43)*y(34)*xnH
r_f(451,43)=-rate
r_f(451,34)=-rate
r_f(451,44)=rate
r_f(451,33)=rate
! 452)   CO+   +   H2O   ->   HCO+  +   OH
!        (43        34         45        32)
rate=xk(452)*y(43)*y(34)*xnH
r_f(452,43)=-rate
r_f(452,34)=-rate
r_f(452,45)=rate
r_f(452,32)=rate
! 453)   CO+   +   C2    ->   C2+   +   CO
!        (43        18         24        33)
rate=xk(453)*y(43)*y(18)*xnH
r_f(453,43)=-rate
r_f(453,18)=-rate
r_f(453,24)=rate
r_f(453,33)=rate
! 454)   CO+   +   HCO   ->   HCO+  +   CO
!        (43        35         45        33)
rate=xk(454)*y(43)*y(35)*xnH
r_f(454,43)=-rate
r_f(454,35)=-rate
r_f(454,45)=rate
r_f(454,33)=rate
! 455)   CO+   +   H2CO  ->   HCO+  +   HCO
!        (43        38         45        35)
rate=xk(455)*y(43)*y(38)*xnH
r_f(455,43)=-rate
r_f(455,38)=-rate
r_f(455,45)=rate
r_f(455,35)=rate
! 456)   CO+   +   H2CO  ->   H2CO+ +   CO
!        (43        38         48        33)
rate=xk(456)*y(43)*y(38)*xnH
r_f(456,43)=-rate
r_f(456,38)=-rate
r_f(456,48)=rate
r_f(456,33)=rate
! 457)   CO+   +   O2    ->   O2+   +   CO
!        (43        31         41        33)
rate=xk(457)*y(43)*y(31)*xnH
r_f(457,43)=-rate
r_f(457,31)=-rate
r_f(457,41)=rate
r_f(457,33)=rate
! 458)   HCO+  +   C     ->   CH+   +   CO
!        (45        17         25        33)
rate=xk(458)*y(45)*y(17)*xnH
r_f(458,45)=-rate
r_f(458,17)=-rate
r_f(458,25)=rate
r_f(458,33)=rate
! 459)   HCO+  +   H-    ->   CO    +   H2
!        (45        7          33        2)
rate=xk(459)*y(45)*y(7)*xnH
r_f(459,45)=-rate
r_f(459,7)=-rate
r_f(459,33)=rate
r_f(459,2)=rate
! 460)   HCO+  +   CH    ->   CH2+  +   CO
!        (45        19         26        33)
rate=xk(460)*y(45)*y(19)*xnH
r_f(460,45)=-rate
r_f(460,19)=-rate
r_f(460,26)=rate
r_f(460,33)=rate
! 461)   HCO+  +   CH2   ->   CH3+  +   CO
!        (45        20         27        33)
rate=xk(461)*y(45)*y(20)*xnH
r_f(461,45)=-rate
r_f(461,20)=-rate
r_f(461,27)=rate
r_f(461,33)=rate
! 462)   HCO+  +   OH    ->   H2O+  +   CO
!        (45        32         44        33)
rate=xk(462)*y(45)*y(32)*xnH
r_f(462,45)=-rate
r_f(462,32)=-rate
r_f(462,44)=rate
r_f(462,33)=rate
! 463)   HCO+  +   OH    ->   HCO2+ +   H
!        (45        32         49        1)
rate=xk(463)*y(45)*y(32)*xnH
r_f(463,45)=-rate
r_f(463,32)=-rate
r_f(463,49)=rate
r_f(463,1)=rate
! 464)   HCO+  +   H2O   ->   H3O+  +   CO
!        (45        34         47        33)
rate=xk(464)*y(45)*y(34)*xnH
r_f(464,45)=-rate
r_f(464,34)=-rate
r_f(464,47)=rate
r_f(464,33)=rate
! 465)   HCO+  +   HCO   ->   H2CO+ +   CO
!        (45        35         48        33)
rate=xk(465)*y(45)*y(35)*xnH
r_f(465,45)=-rate
r_f(465,35)=-rate
r_f(465,48)=rate
r_f(465,33)=rate
! 466)   HCO+  +   H2CO  ->   H3CO+ +   CO
!        (45        38         50        33)
rate=xk(466)*y(45)*y(38)*xnH
r_f(466,45)=-rate
r_f(466,38)=-rate
r_f(466,50)=rate
r_f(466,33)=rate
! 467)   H2CO+ +   CH    ->   CH+   +   H2CO
!        (48        19         25        38)
rate=xk(467)*y(48)*y(19)*xnH
r_f(467,48)=-rate
r_f(467,19)=-rate
r_f(467,25)=rate
r_f(467,38)=rate
! 468)   H2CO+ +   CH    ->   CH2+  +   HCO
!        (48        19         26        35)
rate=xk(468)*y(48)*y(19)*xnH
r_f(468,48)=-rate
r_f(468,19)=-rate
r_f(468,26)=rate
r_f(468,35)=rate
! 469)   H2CO+ +   CH2   ->   CH3+  +   HCO
!        (48        20         27        35)
rate=xk(469)*y(48)*y(20)*xnH
r_f(469,48)=-rate
r_f(469,20)=-rate
r_f(469,27)=rate
r_f(469,35)=rate
! 470)   H2CO+ +   CH2   ->   CH2+  +   H2CO
!        (48        20         26        38)
rate=xk(470)*y(48)*y(20)*xnH
r_f(470,48)=-rate
r_f(470,20)=-rate
r_f(470,26)=rate
r_f(470,38)=rate
! 471)   H2CO+ +   CH4   ->   H3CO+ +   CH3
!        (48        22         50        21)
rate=xk(471)*y(48)*y(22)*xnH
r_f(471,48)=-rate
r_f(471,22)=-rate
r_f(471,50)=rate
r_f(471,21)=rate
! 472)   H2CO+ +   H2O   ->   H3O+  +   HCO
!        (48        34         47        35)
rate=xk(472)*y(48)*y(34)*xnH
r_f(472,48)=-rate
r_f(472,34)=-rate
r_f(472,47)=rate
r_f(472,35)=rate
! 473)   H2CO+ +   HCO   ->   HCO+  +   H2CO
!        (48        35         45        38)
rate=xk(473)*y(48)*y(35)*xnH
r_f(473,48)=-rate
r_f(473,35)=-rate
r_f(473,45)=rate
r_f(473,38)=rate
! 474)   H2CO+ +   HCO   ->   H3CO+ +   CO
!        (48        35         50        33)
rate=xk(474)*y(48)*y(35)*xnH
r_f(474,48)=-rate
r_f(474,35)=-rate
r_f(474,50)=rate
r_f(474,33)=rate
! 475)   H2CO+ +   H2CO  ->   H3CO+ +   HCO
!        (48        38         50        35)
rate=xk(475)*y(48)*y(38)*xnH
r_f(475,48)=-rate
r_f(475,38)=-rate
r_f(475,50)=rate
r_f(475,35)=rate
! 476)   H2CO+ +   O2    ->   HCO+  +   O2H
!        (48        31         45        36)
rate=xk(476)*y(48)*y(31)*xnH
r_f(476,48)=-rate
r_f(476,31)=-rate
r_f(476,45)=rate
r_f(476,36)=rate
! 477)   H3CO+ +   CH    ->   CH2+  +   H2CO
!        (50        19         26        38)
rate=xk(477)*y(50)*y(19)*xnH
r_f(477,50)=-rate
r_f(477,19)=-rate
r_f(477,26)=rate
r_f(477,38)=rate
! 478)   H3CO+ +   H2O   ->   H3O+  +   H2CO
!        (50        34         47        38)
rate=xk(478)*y(50)*y(34)*xnH
r_f(478,50)=-rate
r_f(478,34)=-rate
r_f(478,47)=rate
r_f(478,38)=rate
! 479)   O2+   +   C    ->   CO+   +   O
!        (41        17         43        30)
rate=xk(479)*y(41)*y(17)*xnH
r_f(479,41)=-rate
r_f(479,17)=-rate
r_f(479,43)=rate
r_f(479,30)=rate
! 480)   O2+   +   C     ->   C+    +   O2
!        (41        17         23        31)
rate=xk(480)*y(41)*y(17)*xnH
r_f(480,41)=-rate
r_f(480,17)=-rate
r_f(480,23)=rate
r_f(480,31)=rate
! 481)   O2+   +   CH    ->   CH+   +   O2
!        (41        19         25        31)
rate=xk(481)*y(41)*y(19)*xnH
r_f(481,41)=-rate
r_f(481,19)=-rate
r_f(481,25)=rate
r_f(481,31)=rate
! 482)   O2+   +   CH    ->   HCO+  +   O
!        (41        19         45        30)
rate=xk(482)*y(41)*y(19)*xnH
r_f(482,41)=-rate
r_f(482,19)=-rate
r_f(482,45)=rate
r_f(482,30)=rate
! 483)   O2+   +   CH2   ->   H2CO+ +   O
!        (41        20         48        30)
rate=xk(483)*y(41)*y(20)*xnH
r_f(483,41)=-rate
r_f(483,20)=-rate
r_f(483,48)=rate
r_f(483,30)=rate
! 484)   O2+   +   CH2   ->   CH2+  +   O2
!        (41        20         26        31)
rate=xk(484)*y(41)*y(20)*xnH
r_f(484,41)=-rate
r_f(484,20)=-rate
r_f(484,26)=rate
r_f(484,31)=rate
! 485)   O2+   +   C2    ->   CO+   +   CO
!        (41        18         43        33)
rate=xk(485)*y(41)*y(18)*xnH
r_f(485,41)=-rate
r_f(485,18)=-rate
r_f(485,43)=rate
r_f(485,33)=rate
! 486)   O2+   +   C2    ->   C2+   +   O2
!        (41        18         24        31)
rate=xk(486)*y(41)*y(18)*xnH
r_f(486,41)=-rate
r_f(486,18)=-rate
r_f(486,24)=rate
r_f(486,31)=rate
! 487)   O2+   +   HCO   ->   O2H+  +   CO
!        (41        35         46        33)
rate=xk(487)*y(41)*y(35)*xnH
r_f(487,41)=-rate
r_f(487,35)=-rate
r_f(487,46)=rate
r_f(487,33)=rate
! 488)   O2+   +   HCO   ->   HCO+  +   O2
!        (41        35         45        31)
rate=xk(488)*y(41)*y(35)*xnH
r_f(488,41)=-rate
r_f(488,35)=-rate
r_f(488,45)=rate
r_f(488,31)=rate
! 489)   O2+   +   H2CO  ->   HCO+  +   O2    +   H
!        (41        38         45        31        1)
rate=xk(489)*y(41)*y(38)*xnH
r_f(489,41)=-rate
r_f(489,38)=-rate
r_f(489,45)=rate
r_f(489,31)=rate
r_f(489,1)=rate
! 490)   O2+   +   H2CO  ->   H2CO+ +   O2
!        (41        38         48        31)
rate=xk(490)*y(41)*y(38)*xnH
r_f(490,41)=-rate
r_f(490,38)=-rate
r_f(490,48)=rate
r_f(490,31)=rate
! 491)   O2H+  +   C     ->   CH+   +   O2
!        (46        17         25        31)
rate=xk(491)*y(46)*y(17)*xnH
r_f(491,46)=-rate
r_f(491,17)=-rate
r_f(491,25)=rate
r_f(491,31)=rate
!  492)   O2H+  +   O     ->   OH+   +   O2
!        (46        30         42        31)
rate=xk(492)*y(46)*y(30)*xnH
r_f(492,46)=-rate
r_f(492,30)=-rate
r_f(492,42)=rate
r_f(492,31)=rate
! 493)   O2H+  +   H2    ->   H3+   +   O2
!        (46        2          6         31)
rate=xk(493)*y(46)*y(2)*xnH
r_f(493,46)=-rate
r_f(493,2)=-rate
r_f(493,6)=rate
r_f(493,31)=rate
! 494)   O2H+  +   CH    ->   CH2+  +   O2
!        (46        19         26        31)
rate=xk(494)*y(46)*y(19)*xnH
r_f(494,46)=-rate
r_f(494,19)=-rate
r_f(494,26)=rate
r_f(494,31)=rate
! 495)   O2H+  +   CH2   ->   CH3+  +   O2
!        (46        20         27        31)
rate=xk(495)*y(46)*y(20)*xnH
r_f(495,46)=-rate
r_f(495,20)=-rate
r_f(495,27)=rate
r_f(495,31)=rate
! 496)   O2H+  +   OH    ->   H2O+  +   O2
!        (46        32         44        31)
rate=xk(496)*y(46)*y(32)*xnH
r_f(496,46)=-rate
r_f(496,32)=-rate
r_f(496,44)=rate
r_f(496,31)=rate
! 497)   O2H+  +   H2O   ->   H3O+  +   O2
!        (46        34         47        31)
rate=xk(497)*y(46)*y(34)*xnH
r_f(497,46)=-rate
r_f(497,34)=-rate
r_f(497,47)=rate
r_f(497,31)=rate
! 498)   O2H+  +   CO    ->   HCO+  +   O2
!        (46        33         45        31)
rate=xk(498)*y(46)*y(33)*xnH
r_f(498,46)=-rate
r_f(498,33)=-rate
r_f(498,45)=rate
r_f(498,31)=rate
! 499)   O2H+  +   HCO   ->   H2CO+ +   O2
!        (46        35         48        31)
rate=xk(499)*y(46)*y(35)*xnH
r_f(499,46)=-rate
r_f(499,35)=-rate
r_f(499,48)=rate
r_f(499,31)=rate
! 500)   O2H+  +   H2CO  ->   H3CO+ +   O2
!        (46        38         50        31)
rate=xk(500)*y(46)*y(38)*xnH
r_f(500,46)=-rate
r_f(500,38)=-rate
r_f(500,50)=rate
r_f(500,31)=rate
! 501)   O2H+  +   CO2   ->   HCO2+ +   O2
!        (46        37         49        31)
rate=xk(501)*y(46)*y(37)*xnH
r_f(501,46)=-rate
r_f(501,37)=-rate
r_f(501,49)=rate
r_f(501,31)=rate
! 502)   HCO2+ +   C    ->   CH+   +   CO2
!        (49        17         25        37)
rate=xk(502)*y(49)*y(17)*xnH
r_f(502,49)=-rate
r_f(502,17)=-rate
r_f(502,25)=rate
r_f(502,37)=rate
! 503)   HCO2+ +   O     ->   HCO+  +   O2
!        (49        30         45        31)
rate=xk(503)*y(49)*y(30)*xnH
r_f(503,49)=-rate
r_f(503,30)=-rate
r_f(503,45)=rate
r_f(503,31)=rate
! 504)   HCO2+ +   CH4   ->   CH5+  +   CO2
!        (49        22         29        37)
rate=xk(504)*y(49)*y(22)*xnH
r_f(504,49)=-rate
r_f(504,22)=-rate
r_f(504,29)=rate
r_f(504,37)=rate
! 505)   HCO2+ +   H2O   ->   H3O+  +   CO2
!        (49        34         47        37)
rate=xk(505)*y(49)*y(34)*xnH
r_f(505,49)=-rate
r_f(505,34)=-rate
r_f(505,47)=rate
r_f(505,37)=rate
! 506)   HCO2+ +   CO    ->   HCO+  +   CO2
!        (49        33         45        37)
rate=xk(506)*y(49)*y(33)*xnH
r_f(506,49)=-rate
r_f(506,33)=-rate
r_f(506,45)=rate
r_f(506,37)=rate
! 507)   C+    +   e     ->   C     +   ph.
!        (23        3          17)
rate=xk(507)*y(23)*y(3)*xnH
r_f(507,23)=-rate
r_f(507,3)=-rate
r_f(507,17)=rate
! 508)   CH+   +   e     ->   C     +   H
!        (25        3          17        1)
rate=xk(508)*y(25)*y(3)*xnH
r_f(508,25)=-rate
r_f(508,3)=-rate
r_f(508,17)=rate
r_f(508,1)=rate
! 509)   CH2+  +   e     ->   C     +   H2
!        (26        3          17        2)
rate=xk(509)*y(26)*y(3)*xnH
r_f(509,26)=-rate
r_f(509,3)=-rate
r_f(509,17)=rate
r_f(509,2)=rate
! 510)   CH2+  +   e     ->   CH    +   H
!        (26        3          19        1)
rate=xk(510)*y(26)*y(3)*xnH
r_f(510,26)=-rate
r_f(510,3)=-rate
r_f(510,19)=rate
r_f(510,1)=rate
! 511)   CH3+  +   e     ->   CH2   +   H
!        (27        3          20        1)
rate=xk(511)*y(27)*y(3)*xnH
r_f(511,27)=-rate
r_f(511,3)=-rate
r_f(511,20)=rate
r_f(511,1)=rate
! 512)   CH3+  +   e     ->   CH    +   H2
!        (27        3          19        2)
rate=xk(512)*y(27)*y(3)*xnH
r_f(512,27)=-rate
r_f(512,3)=-rate
r_f(512,19)=rate
r_f(512,2)=rate
! 513)   CH3+  +   e     ->   CH    + 2 H
!        (27        3          19      2*1)
rate=xk(513)*y(27)*y(3)*xnH
r_f(513,27)=-rate
r_f(513,3)=-rate
r_f(513,19)=rate
r_f(513,1)=2.d0*rate
! 514)   CH3+  +   e     ->   CH3   +   ph.
!        (27        3          21)
rate=xk(514)*y(27)*y(3)*xnH
r_f(514,27)=-rate
r_f(514,3)=-rate
r_f(514,21)=rate
! 515)   O+    +   e     ->   O     +   ph.
!        (40        3          30)
rate=xk(515)*y(40)*y(3)*xnH
r_f(515,40)=-rate
r_f(515,3)=-rate
r_f(515,30)=rate
! 516)   CH4+  +   e     ->   CH3   +   H
!        (28        3          21        1)
rate=xk(516)*y(28)*y(3)*xnH
r_f(516,28)=-rate
r_f(516,3)=-rate
r_f(516,21)=rate
r_f(516,1)=rate
! 517)   CH4+  +   e     ->   CH2   + 2 H
!        (28        3          20      2*1)
rate=xk(517)*y(28)*y(3)*xnH
r_f(517,28)=-rate
r_f(517,3)=-rate
r_f(517,20)=rate
r_f(517,1)=2.d0*rate
! 518)   OH+   +   e     ->   O     +   H
!        (42        3          30        1)
rate=xk(518)*y(42)*y(3)*xnH
r_f(518,42)=-rate
r_f(518,3)=-rate
r_f(518,30)=rate
r_f(518,1)=rate
! 519)   CH5+  +   e     ->   CH4   +   H
!        (29        3          22        1)
rate=xk(519)*y(29)*y(3)*xnH
r_f(519,29)=-rate
r_f(519,3)=-rate
r_f(519,22)=rate
r_f(519,1)=rate
! 520)   CH5+  +   e     ->   CH3   +   H2
!        (29        3          21        2)
rate=xk(520)*y(29)*y(3)*xnH
r_f(520,29)=-rate
r_f(520,3)=-rate
r_f(520,21)=rate
r_f(520,2)=rate
! 521)   H2O+  +   e     ->   OH    +   H
!        (44        3          32        1)
rate=xk(521)*y(44)*y(3)*xnH
r_f(521,44)=-rate
r_f(521,3)=-rate
r_f(521,32)=rate
r_f(521,1)=rate
! 522)   H2O+  +   e     ->   O     +   H2
!        (44        3          30        2)
rate=xk(522)*y(44)*y(3)*xnH
r_f(522,44)=-rate
r_f(522,3)=-rate
r_f(522,30)=rate
r_f(522,2)=rate
! 523)   H3O+  +   e     ->   H2O   +   H
!        (47        3          34        1)
rate=xk(523)*y(47)*y(3)*xnH
r_f(523,47)=-rate
r_f(523,3)=-rate
r_f(523,34)=rate
r_f(523,1)=rate
! 524)   H3O+  +   e     ->   OH    + 2 H
!        (47        3          32      2*1)
rate=xk(524)*y(47)*y(3)*xnH
r_f(524,47)=-rate
r_f(524,3)=-rate
r_f(524,32)=rate
r_f(524,1)=2.d0*rate
! 525)   C2+   +   e     -> 2 C
!        (24        3        2*17)
rate=xk(525)*y(24)*y(3)*xnH
r_f(525,24)=-rate
r_f(525,3)=-rate
r_f(525,17)=2.d0*rate
! 526)   CO+   +   e     ->   O     +   C
!        (43        3          30        17)
rate=xk(526)*y(43)*y(3)*xnH
r_f(526,43)=-rate
r_f(526,3)=-rate
r_f(526,30)=rate
r_f(526,17)=rate
! 527)   HCO+  +   e     ->   CO    +   H
!        (45        3          33        1)
rate=xk(527)*y(45)*y(3)*xnH
r_f(527,45)=-rate
r_f(527,3)=-rate
r_f(527,33)=rate
r_f(527,1)=rate
! 528)   H2CO+ +   e     ->   HCO   +   H
!        (48        3          35        1)
rate=xk(528)*y(48)*y(3)*xnH
r_f(528,48)=-rate
r_f(528,3)=-rate
r_f(528,35)=rate
r_f(528,1)=rate
! 529)   H2CO+ +   e     ->   CO    + 2 H
!        (48        3          33      2*1)
rate=xk(529)*y(48)*y(3)*xnH
r_f(529,48)=-rate
r_f(529,3)=-rate
r_f(529,33)=rate
r_f(529,1)=2.d0*rate
! 530)   H2CO+ +   e     ->   H2CO  +   ph.
!        (48        3          38)
rate=xk(530)*y(48)*y(3)*xnH
r_f(530,48)=-rate
r_f(530,3)=-rate
r_f(530,38)=rate
! 531)   H3CO+ +   e     ->   CO    +   H     +   H2
!        (50        3          33        1         2)
rate=xk(531)*y(50)*y(3)*xnH
r_f(531,50)=-rate
r_f(531,3)=-rate
r_f(531,33)=rate
r_f(531,1)=rate
r_f(531,2)=rate
! 532)   H3CO+ +   e     ->   HCO   + 2 H
!        (50        3          35      2*1)
rate=xk(532)*y(50)*y(3)*xnH
r_f(532,50)=-rate
r_f(532,3)=-rate
r_f(532,35)=rate
r_f(532,1)=2.d0*rate
! 533)   H3CO+ +   e     ->   H2CO  +   H
!        (50        3          38        1)
rate=xk(533)*y(50)*y(3)*xnH
r_f(533,50)=-rate
r_f(533,3)=-rate
r_f(533,38)=rate
r_f(533,1)=rate
! 534)   O2+   +   e     -> 2 O
!        (41        3        2*30)
rate=xk(534)*y(41)*y(3)*xnH
r_f(534,41)=-rate
r_f(534,3)=-rate
r_f(534,30)=2.d0*rate
! 535)   O2H+  +   e     ->   O2    +   H
!        (46        3          31        1)
rate=xk(535)*y(46)*y(3)*xnH
r_f(535,46)=-rate
r_f(535,3)=-rate
r_f(535,31)=rate
r_f(535,1)=rate
! 536)   HCO2+ +   e     ->   CO2   +   H
!        (49        3          37        1)
rate=xk(536)*y(49)*y(3)*xnH
r_f(536,49)=-rate
r_f(536,3)=-rate
r_f(536,37)=rate
r_f(536,1)=rate
! 537)   HCO2+ +   e     ->   CO    +   OH
!        (49        3          33        32)
rate=xk(537)*y(49)*y(3)*xnH
r_f(537,49)=-rate
r_f(537,3)=-rate
r_f(537,33)=rate
r_f(537,32)=rate
! 538)    H2    +   C     ->   CH2   +   ph.
!        (2        17         20)
rate=xk(538)*y(2)*y(17)*xnH
r_f(538,2)=-rate
r_f(538,17)=-rate
r_f(538,20)=rate
! 539)    OH    +   CH3   ->   CH4   +   O
!        (32       21         22        30)
rate=xk(539)*y(32)*y(21)*xnH
r_f(539,32)=-rate
r_f(539,21)=-rate
r_f(539,22)=rate
r_f(539,30)=rate
! 540)    H+    +   H2CO  ->   CO+   +   H2   +   H
!       (4         38         43        2        1)
rate=xk(540)*y(4)*y(38)*xnH
r_f(540,4)=-rate
r_f(540,38)=-rate
r_f(540,43)=rate
r_f(540,2)=rate
r_f(540,1)=rate
! 541)    He+   +   C     ->   C+    +   He
!       (9         17         23        8)
rate=xk(541)*y(9)*y(17)*xnH
r_f(541,9)=-rate
r_f(541,17)=-rate
r_f(541,23)=rate
r_f(541,8)=rate
! 542)    He+   +   H2CO  ->   H2CO+ +   He
!       (9         38         48        8)
rate=xk(542)*y(9)*y(38)*xnH
r_f(542,9)=-rate
r_f(542,38)=-rate
r_f(542,48)=rate
r_f(542,8)=rate
! 543)    He+   +   H2CO  ->   CH2+  +   O    +   He
!       (9         38         26        30       8)
rate=xk(543)*y(9)*y(38)*xnH
r_f(543,9)=-rate
r_f(543,38)=-rate
r_f(543,26)=rate
r_f(543,30)=rate
r_f(543,8)=rate
! 544)   H     +   CR    ->   H+    +   e
!        (1                    4         3)
rate=xk(544)*y(1)
r_f(544,1)=-rate
r_f(544,4)=rate
r_f(544,3)=rate
! 545)   He    +   CR    ->   He+   +   e
!        (8                    9         3)
rate=xk(545)*y(8)
r_f(545,8)=-rate
r_f(545,9)=rate
r_f(545,3)=rate
! 546)   C     +   CR    ->   C+    +   e
!        (17                   23        3)
rate=xk(546)*y(17)
r_f(546,17)=-rate
r_f(546,23)=rate
r_f(546,3)=rate
! 547)   O     +   CR    ->   O+    +   e
!        (30                   40        3)
rate=xk(547)*y(30)
r_f(547,30)=-rate
r_f(547,40)=rate
r_f(547,3)=rate
! 548)   H2    +   CR    ->   H+    +   H     +   e
!        (2                    4         1         3)
rate=xk(548)*y(2)
r_f(548,2)=-rate
r_f(548,4)=rate
r_f(548,1)=rate
r_f(548,3)=rate
! 549)   H2    +   CR    ->   H2+   +   e
!        (2                    5         3)
rate=xk(549)*y(2)
r_f(549,2)=-rate
r_f(549,5)=rate
r_f(549,3)=rate
! 550)   H2    +   CR    -> 2 H
!        (2                  2*1)
rate=xk(550)*y(2)
r_f(550,2)=-rate
r_f(550,1)=2.d0*rate
! 551)   H2    +   CR    ->   H+    +   H-
!        (2                    4         7)
rate=xk(551)*y(2)
r_f(551,2)=-rate
r_f(551,4)=rate
r_f(551,7)=rate
! 552)   CO    +   CR    ->   CO+   +   e
!        (33                   43        3)
rate=xk(552)*y(33)
r_f(552,33)=-rate
r_f(552,43)=rate
r_f(552,3)=rate
! 553)   C     +   ph.   ->   C+    +   e
!        (17                   23        3)
rate=xk(553)*y(17)
r_f(553,17)=-rate
r_f(553,23)=rate
r_f(553,3)=rate
! 554)   H-    +   ph.   ->   H     +   e
!        (7                    1         3)
rate=xk(554)*y(7)
r_f(554,7)=-rate
r_f(554,1)=rate
r_f(554,3)=rate
! 555)   H2+   +   ph.   ->   H+    +   H
!        (5                    4         1)
rate=xk(555)*y(5)
r_f(555,5)=-rate
r_f(555,4)=rate
r_f(555,1)=rate
! 556)   H3+   +   ph.   ->   H2+   +   H
!        (6                    5         1)
rate=xk(556)*y(6)
r_f(556,6)=-rate
r_f(556,5)=rate
r_f(556,1)=rate
! 557)   H3+   +   ph.   ->   H+    +   H2
!        (6                    4         2)
rate=xk(557)*y(6)
r_f(557,6)=-rate
r_f(557,4)=rate
r_f(557,2)=rate
! 558)   CH    +   ph.   ->   CH+   +   e
!        (19                   25        3)
rate=xk(558)*y(19)
r_f(558,19)=-rate
r_f(558,25)=rate
r_f(558,3)=rate
! 559)   CH    +   ph.   ->   C     +   H
!        (19                   17        1)
rate=xk(559)*y(19)
r_f(559,19)=-rate
r_f(559,17)=rate
r_f(559,1)=rate
! 560)   CH+   +   ph.   ->   C+    +   H
!        (25                   23        1)
rate=xk(560)*y(25)
r_f(560,25)=-rate
r_f(560,23)=rate
r_f(560,1)=rate
! 561)   CH2   +   ph.   ->   CH2+  +   e
!        (20                   26        3)
rate=xk(561)*y(20)
r_f(561,20)=-rate
r_f(561,26)=rate
r_f(561,3)=rate
! 562)   CH2   +   ph.   ->   CH    +   H
!        (20                   19        1)
rate=xk(562)*y(20)
r_f(562,20)=-rate
r_f(562,19)=rate
r_f(562,1)=rate
! 563)   CH2+  +   ph.   ->   CH+   +   H
!        (26                   25        1)
rate=xk(563)*y(26)
r_f(563,26)=-rate
r_f(563,25)=rate
r_f(563,1)=rate
! 564)   CH3   +   ph.   ->   CH    +   H2
!        (21                   19        2)
rate=xk(564)*y(21)
r_f(564,21)=-rate
r_f(564,19)=rate
r_f(564,2)=rate
! 565)   CH3   +   ph.   ->   CH2   +   H
!        (21                   20        1)
rate=xk(565)*y(21)
r_f(565,21)=-rate
r_f(565,20)=rate
r_f(565,1)=rate
! 566)   CH3   +   ph.   ->   CH3+  +   e
!        (21                   27        3)
rate=xk(566)*y(21)
r_f(566,21)=-rate
r_f(566,27)=rate
r_f(566,3)=rate
! 567)   CH3+  +   ph.   ->   CH2+  +   H
!        (27                   26        1)
rate=xk(567)*y(27)
r_f(567,27)=-rate
r_f(567,26)=rate
r_f(567,1)=rate
! 568)   CH3+  +   ph.   ->   CH+   +   H2
!        (27                   25        2)
rate=xk(568)*y(27)
r_f(568,27)=-rate
r_f(568,25)=rate
r_f(568,2)=rate
! 569)   CH4   +   ph.   ->   CH3   +   H
!        (22                   21        1)
rate=xk(569)*y(22)
r_f(569,22)=-rate
r_f(569,21)=rate
r_f(569,1)=rate
! 570)   CH4   +   ph.   ->   CH2   +   H2
!        (22                   20        2)
rate=xk(570)*y(22)
r_f(570,22)=-rate
r_f(570,20)=rate
r_f(570,2)=rate
! 571)   CH4   +   ph.   ->   CH    +   H     +   H2
!        (22                   19        1         2)
rate=xk(571)*y(22)
r_f(571,22)=-rate
r_f(571,19)=rate
r_f(571,1)=rate
r_f(571,2)=rate
! 572)   OH    +   ph.   ->   O     +   H
!        (32                   30        1)
rate=xk(572)*y(32)
r_f(572,32)=-rate
r_f(572,30)=rate
r_f(572,1)=rate
! 573)   OH    +   ph.   ->   OH+   +   e
!        (32                   42        3)
rate=xk(573)*y(32)
r_f(573,32)=-rate
r_f(573,42)=rate
r_f(573,3)=rate
! 574)   OH+   +   ph.   ->   H+    +   O
!        (42                   4         30)
rate=xk(574)*y(42)
r_f(574,42)=-rate
r_f(574,4)=rate
r_f(574,30)=rate
! 575)   H2O   +   ph.   ->   OH    +   H
!        (34                   32        1)
rate=xk(575)*y(34)
r_f(575,34)=-rate
r_f(575,32)=rate
r_f(575,1)=rate
! 576)   H2O   +   ph.   ->   H2O+  +   e
!        (34                   44        3)
rate=xk(576)*y(34)
r_f(576,34)=-rate
r_f(576,44)=rate
r_f(576,3)=rate
! 577)   C2    +   ph.   -> 2 C
!        (18                 2*17)
rate=xk(577)*y(18)
r_f(577,18)=-rate
r_f(577,17)=2.d0*rate
! 578)   C2    +   ph.   ->   C2+   +   e
!        (18                   24        3)
rate=xk(578)*y(18)
r_f(578,18)=-rate
r_f(578,24)=rate
r_f(578,3)=rate
! 579)   C2+   +   ph.   ->   C+    +   C
!        (24                   23        17)
rate=xk(579)*y(24)
r_f(579,24)=-rate
r_f(579,23)=rate
r_f(579,17)=rate
! 580)   CO    +   ph.   ->   C    +   O
!        (33                   17        30)
rate=xk(580)*y(33)
r_f(580,33)=-rate
r_f(580,17)=rate
r_f(580,30)=rate
! 581)   HCO   +   ph.   ->   H     +   CO
!        (35                   1         33)
rate=xk(581)*y(35)
r_f(581,35)=-rate
r_f(581,1)=rate
r_f(581,33)=rate
! 582)   HCO   +   ph.   ->   HCO+  +   e
!        (35                   45        3)
rate=xk(582)*y(35)
r_f(582,35)=-rate
r_f(582,45)=rate
r_f(582,3)=rate
! 583)   H2CO  +   ph.   ->   CO    +   H2
!        (38                   33        2)
rate=xk(583)*y(38)
r_f(583,38)=-rate
r_f(583,33)=rate
r_f(583,2)=rate
! 584)   H2CO  +   ph.   ->   CO    + 2 H
!        (38                   33      2*1)
rate=xk(584)*y(38)
r_f(584,38)=-rate
r_f(584,33)=rate
r_f(584,1)=2.d0*rate
! 585)   H2CO  +   ph.   ->   H2CO+ +   e
!        (38                   48        3)
rate=xk(585)*y(38)
r_f(585,38)=-rate
r_f(585,48)=rate
r_f(585,3)=rate
! 586)   H2CO  +   ph.   ->   HCO+  +   e     +   H
!        (38                   45        3         1)
rate=xk(586)*y(38)
r_f(586,38)=-rate
r_f(586,45)=rate
r_f(586,3)=rate
r_f(586,1)=rate
! 587)   O2    +   ph.   -> 2 O
!        (31                 2*30)
rate=xk(587)*y(31)
r_f(587,31)=-rate
r_f(587,30)=2.d0*rate
! 588)   O2    +   ph.   ->   O2+   +   e
!        (31                   41        3)
rate=xk(588)*y(31)
r_f(588,31)=-rate
r_f(588,41)=rate
r_f(588,3)=rate
! 589)   CO2   +   ph.   ->   CO    +   O
!        (37                   33        30)
rate=xk(589)*y(37)
r_f(589,37)=-rate
r_f(589,33)=rate
r_f(589,30)=rate
! 590)   H2   +   ph.   -> 2 H
!       (2                   2*1)
rate=xk(590)*y(2)
r_f(590,2)=-rate
r_f(590,1)=2.d0*rate
!591)   HD   +   ph.   ->    H     +   D
!       (13                  1         12)
rate=xk(591)*y(13)
r_f(591,13)=-rate
r_f(591,1)=rate
r_f(591,12)=rate
!supplemented reactions
!601)   H    +   HCO   ->   O   +   CH2
!      (1        35         30      20  )
rate=xk(601)*y(1)*y(35)*xnH
r_f(601,1)=-rate
r_f(601,35)=-rate
r_f(601,30)=rate
r_f(601,20)=rate
!602)   H    +   H2O2  ->  H2O   +   OH
!      (1        39        34        32 )
rate=xk(602)*y(1)*y(39)*xnH
r_f(602,1)=-rate
r_f(602,39)=-rate
r_f(602,34)=rate
r_f(602,32)=rate
!603)   C    +   CH2   ->   2 CH
!      (17       20         2x19         )
rate=xk(603)*y(17)*y(20)*xnH
r_f(603,17)=-rate
r_f(603,20)=-rate
r_f(603,19)=2.d0*rate
!604)   CH   +    O2   ->   CO   +   OH
!      (19        31        33       32  )
rate=xk(604)*y(19)*y(31)*xnH
r_f(604,19)=-rate
r_f(604,31)=-rate
r_f(604,33)=rate
r_f(604,32)=rate
!605)   H    +   HCO   ->   O   +   CH2
!      (1        35         30      20  )
rate=xk(605)*y(1)*y(35)*xnH
r_f(605,1)=-rate
r_f(605,35)=-rate
r_f(605,30)=rate
r_f(605,20)=rate
!606)   CH2  +   O2   ->   CO2    +   2 H
!      (20       31        37         2x1 )
rate=xk(606)*y(20)*y(31)*xnH
r_f(606,20)=-rate
r_f(606,31)=-rate
r_f(606,37)=rate
r_f(606,1)=2.d0*rate
!607)   CH2  +   O2   ->   CO2    +   H2
!      (20       31        37         2  )
rate=xk(607)*y(20)*y(31)*xnH
r_f(607,20)=-rate
r_f(607,31)=-rate
r_f(607,37)=rate
r_f(607,2)=rate
!608)   CH2  +   O2   ->   CO     +   H2O
!      (20       31        33         34 )
rate=xk(608)*y(20)*y(31)*xnH
r_f(608,20)=-rate
r_f(608,31)=-rate
r_f(608,33)=rate
r_f(608,34)=rate
!609)   2 CH3         ->   CH4    +   CH2
!      (2x21               22         20 )
rate=xk(609)*(y(21)**2)*xnH
r_f(609,21)=-2.d0*rate
r_f(609,22)=rate
r_f(609,20)=rate
!610)   CH3  +   O    ->   CO     +   H2    +   H
!      (21       30        33         2         1 )
rate=xk(610)*y(21)*y(30)*xnH
r_f(610,21)=-rate
r_f(610,30)=-rate
r_f(610,33)=rate
r_f(610,2)=rate
r_f(610,1)=rate
!611)   CH3  +   OH   ->   H2CO   +   H2
!      (21       32        38         2 )
rate=xk(611)*y(21)*y(32)*xnH
r_f(611,21)=-rate
r_f(611,32)=-rate
r_f(611,38)=rate
r_f(611,2)=rate
!612)   CH3  +   O2   ->   HCO    +   H2O
!      (21       31        35         34 )
rate=xk(612)*y(21)*y(31)*xnH
r_f(612,21)=-rate
r_f(612,31)=-rate
r_f(612,35)=rate
r_f(612,34)=rate
!613)   O    +  H2CO  ->   CO     +   OH  +   H
!      (30      38         33         32      1 )
rate=xk(613)*y(30)*y(38)*xnH
r_f(613,30)=-rate
r_f(613,38)=-rate
r_f(613,33)=rate
r_f(613,32)=rate
r_f(613,1)=rate
!614)   C2   +   O2   ->   2 CO
!      (18       31        2x33          )
rate=xk(614)*y(18)*y(31)*xnH
r_f(614,18)=-rate
r_f(614,31)=-rate
r_f(614,33)=2.d0*rate
!615)   2 HCO         ->   2 CO   +   H2
!      (2x35               2x33       2   )
rate=xk(615)*(y(35)**2)*xnH
r_f(615,35)=-2.d0*rate
r_f(615,33)=2.d0*rate
r_f(615,2)=rate
!616)   HCO   +   O2   ->   CO2   +   OH
!      (35        31        37        32 )
rate=xk(616)*y(35)*y(31)*xnH
r_f(616,35)=-rate
r_f(616,31)=-rate
r_f(616,37)=rate
r_f(616,32)=rate
!617)    O    +   OH   ->   O2    +   H
!       (30       32        31        1  )
rate=xk(617)*y(30)*y(32)*xnH
r_f(617,30)=-rate
r_f(617,32)=-rate
r_f(617,31)=rate
r_f(617,1)=rate
!621)    H3+  +    O    ->  H2O+  +   H
!       (6         30       44        1 )
rate=xk(621)*y(6)*y(30)*xnH
r_f(621,6)=-rate
r_f(621,30)=-rate
r_f(621,44)=rate
r_f(621,1)=rate
!626)   O+    +  CO     ->  CO+   +   O
!      (40       33         43        30  )
rate=xk(626)*y(40)*y(33)*xnH
r_f(626,40)=-rate
r_f(626,33)=-rate
r_f(626,43)=rate
r_f(626,30)=rate
!632)   CH2+  +  e      ->  C    +  2 H
!      (26       3          17       2x1 )
rate=xk(632)*y(26)*y(3)*xnH
r_f(632,26)=-rate
r_f(632,3)=-rate
r_f(632,17)=rate
r_f(632,1)=2.d0*rate
!633)   CH5+  +  e      ->  CH3   +  2 H
!      (29       3          21       2x1 )
rate=xk(633)*y(29)*y(3)*xnH
r_f(633,29)=-rate
r_f(633,3)=-rate
r_f(633,21)=rate
r_f(633,1)=2.d0*rate
!634)   CH5+  +  e      ->  CH2   +    H2    +   H
!      (29       3          20         2         1)
rate=xk(634)*y(29)*y(3)*xnH
r_f(634,29)=-rate
r_f(634,3)=-rate
r_f(634,20)=rate
r_f(634,2)=rate
r_f(634,1)=rate
!635)   CH5+  +  e      ->  CH    +  2 H2
!      (29       3          19       2x2 )
rate=xk(635)*y(29)*y(3)*xnH
r_f(635,29)=-rate
r_f(635,3)=-rate
r_f(635,19)=rate
r_f(635,2)=2.d0*rate
!636)   H2O+  +  e      ->  O     +  2 H
!      (44       3          30       2x1 )
rate=xk(636)*y(44)*y(3)*xnH
r_f(636,44)=-rate
r_f(636,3)=-rate
r_f(636,30)=rate
r_f(636,1)=2.d0*rate
!637)   H3O+  +  e      ->  O     +    H2   +    H
!      (47       3          30         2         1)
rate=xk(637)*y(47)*y(3)*xnH
r_f(637,47)=-rate
r_f(637,3)=-rate
r_f(637,30)=rate
r_f(637,2)=rate
r_f(637,1)=rate
!638)   H3O+  +  e      ->  OH    +    H2
!      (47       3          32         2)
rate=xk(638)*y(47)*y(3)*xnH
r_f(638,47)=-rate
r_f(638,3)=-rate
r_f(638,32)=rate
r_f(638,2)=rate
!640)   HCO2+ +  e      ->  CO    +    O   +    H
!      (49       3          33         30       1)
rate=xk(640)*y(49)*y(3)*xnH
r_f(640,49)=-rate
r_f(640,3)=-rate
r_f(640,33)=rate
r_f(640,30)=rate
r_f(640,1)=rate
!641)   H+    +  He     ->  HeH+  +    ph.
!      (4        8          11       )
rate=xk(641)*y(4)*y(8)*xnH
r_f(641,4)=-rate
r_f(641,8)=-rate
r_f(641,11)=rate
!642)   H     +  OH     ->  H2O   +    ph.
!      (1        32         34        )
rate=xk(642)*y(1)*y(32)*xnH
r_f(642,1)=-rate
r_f(642,32)=-rate
r_f(642,34)=rate
!643)   H2    +  CH     ->  CH3   +    ph.
!      (2        19         21        )
rate=xk(643)*y(2)*y(19)*xnH
r_f(643,2)=-rate
r_f(643,19)=-rate
r_f(643,21)=rate
!644)   C+    +  C     ->  C2+   +    ph.
!      (23       17         24         )
rate=xk(644)*y(23)*y(17)*xnH
r_f(644,23)=-rate
r_f(644,17)=-rate
r_f(644,24)=rate
!645)   C     +  O+     ->  CO+   +    ph.
!      (17       40         43         )
rate=xk(645)*y(17)*y(40)*xnH
r_f(645,17)=-rate
r_f(645,40)=-rate
r_f(645,43)=rate
!646)   CH2+  +  ph.    ->  CH    +    H+
!      (26                  19         4)
rate=xk(646)*y(26)
r_f(646,26)=-rate
r_f(646,19)=rate
r_f(646,4)=rate
!647)   CH2+  +  ph.    ->  C+    +    H2
!      (26                  23         2)
rate=xk(647)*y(26)
r_f(647,26)=-rate
r_f(647,23)=rate
r_f(647,2)=rate
!648)   CH4+  +  ph.    ->  CH2+  +    H2
!      (28                  26         2)
rate=xk(648)*y(28)
r_f(648,28)=-rate
r_f(648,26)=rate
r_f(648,2)=rate
!649)   CH4+  +  ph.    ->  CH3+  +    H
!      (28                  27         1)
rate=xk(649)*y(28)
r_f(649,28)=-rate
r_f(649,27)=rate
r_f(649,1)=rate
!650)   OH+   +  ph.    ->  O+    +    H
!      (42                  40         1)
rate=xk(650)*y(42)
r_f(650,42)=-rate
r_f(650,40)=rate
r_f(650,1)=rate
!651)   H2O+   +  ph.   ->  OH+   +    H
!      (44                  42         1)
rate=xk(651)*y(44)
r_f(651,44)=-rate
r_f(651,42)=rate
r_f(651,1)=rate
!652)   CO+   +  ph.   ->   C+    +    O
!      (43                  23         30)
rate=xk(652)*y(43)
r_f(652,43)=-rate
r_f(652,23)=rate
r_f(652,30)=rate
!653)   HCO+  +  ph.   ->   CO+   +    H
!      (45                  43         1)
rate=xk(653)*y(45)
r_f(653,45)=-rate
r_f(653,43)=rate
r_f(653,1)=rate
!654)   O2+   +  ph.   ->   O+    +    O
!      (41                  40         30)
rate=xk(654)*y(41)
r_f(654,41)=-rate
r_f(654,40)=rate
r_f(654,30)=rate
!655)   H2O2  +  ph.   -> 2 OH
!      (39                2x32         )
rate=xk(655)*y(39)
r_f(655,39)=-rate
r_f(655,32)=2.d0*rate
!656)   C     + CR ph. ->   C+    +    e
!      (17                  23         3)
rate=xk(656)*y(17)
r_f(656,17)=-rate
r_f(656,23)=rate
r_f(656,3)=rate
!657)   CH    + CR ph. ->   C    +    H
!      (19                  17         1)
rate=xk(657)*y(19)
r_f(657,19)=-rate
r_f(657,17)=rate
r_f(657,1)=rate
!658)   CH+   + CR ph. ->   C+    +    H
!      (25                  23         1)
rate=xk(658)*y(25)
r_f(658,25)=-rate
r_f(658,23)=rate
r_f(658,1)=rate
!659)   CH2   + CR ph. ->   CH2+  +    e
!      (20                  26         3)
rate=xk(659)*y(20)
r_f(659,20)=-rate
r_f(659,26)=rate
r_f(659,3)=rate
!660)   CH2   + CR ph. ->   CH    +    H
!      (20                  19         1)
rate=xk(660)*y(20)
r_f(660,20)=-rate
r_f(660,19)=rate
r_f(660,1)=rate
!661)   CH3   + CR ph. ->   CH3+  +    e
!      (21                  27         3)
rate=xk(661)*y(21)
r_f(661,21)=-rate
r_f(661,27)=rate
r_f(661,3)=rate
!662)   CH3   + CR ph. ->   CH2   +    H
!      (21                  20         1)
rate=xk(662)*y(21)
r_f(662,21)=-rate
r_f(662,20)=rate
r_f(662,1)=rate
!663)   CH3   + CR ph. ->   CH    +    H2
!      (21                  19         2)
rate=xk(663)*y(21)
r_f(663,21)=-rate
r_f(663,19)=rate
r_f(663,2)=rate
!664)   CH4   + CR ph. ->   CH2   +    H2
!      (22                  20         2)
rate=xk(664)*y(22)
r_f(664,22)=-rate
r_f(664,20)=rate
r_f(664,2)=rate
!665)   OH    + CR ph. ->   O     +    H
!      (32                  30         1)
rate=xk(665)*y(32)
r_f(665,32)=-rate
r_f(665,30)=rate
r_f(665,1)=rate
!666)   H2O   + CR ph. ->   OH    +    H
!      (34                  32         1)
rate=xk(666)*y(34)
r_f(666,34)=-rate
r_f(666,32)=rate
r_f(666,1)=rate
!667)   C2    + CR ph. -> 2 C
!      (18                2x17         )
rate=xk(667)*y(18)
r_f(667,18)=-rate
r_f(667,17)=2.d0*rate
!668)   CO    + CR ph. ->   O     +    C
!      (33                  30         17)
rate=xk(668)*y(33)
r_f(668,33)=-rate
r_f(668,30)=rate
r_f(668,17)=rate
!669)  HCO    + CR ph. ->   CO    +    H
!     (35                   33         1)
rate=xk(669)*y(35)
r_f(669,35)=-rate
r_f(669,33)=rate
r_f(669,1)=rate
!670)  HCO    + CR ph. ->   HCO+  +    e
!     (35                   45         3)
rate=xk(670)*y(35)
r_f(670,35)=-rate
r_f(670,45)=rate
r_f(670,3)=rate
!671)  H2CO   + CR ph. ->   CO    +    H2
!     (38                   33         2)
rate=xk(671)*y(38)
r_f(671,38)=-rate
r_f(671,33)=rate
r_f(671,2)=rate
!672)  O2     + CR ph. -> 2 O
!     (31                 2x30        )
rate=xk(672)*y(31)
r_f(672,31)=-rate
r_f(672,30)=2.d0*rate
!673)  O2     + CR ph. ->   O2+   +    e
!     (31                   41         3)
rate=xk(673)*y(31)
r_f(673,31)=-rate
r_f(673,41)=rate
r_f(673,3)=rate
!674)  H2O2   + CR ph. -> 2 OH
!     (39                 2x32        )
rate=xk(674)*y(39)
r_f(674,39)=-rate
r_f(674,32)=2.d0*rate
!675)  CO2    + CR ph. ->   CO    +    O
!     (37                   33         30)
rate=xk(675)*y(37)
r_f(675,37)=-rate
r_f(675,33)=rate
r_f(675,30)=rate

!************************************************************

do isp=1,N_sp
   r_f_tot(isp)=0.d0
   do ire=1,N_react
      r_f_tot(isp)=r_f_tot(isp)+r_f(ire,isp)
   enddo
enddo

return
END

FUNCTION xk_prm(T,rho)
!    continuum opacity of primordial gas
!    input   T    temperature (K)
!            rho  density (g/cm^3)
!    output  xk   mean opacity (cm^2/g)
IMPLICIT REAL*8(a-h,o-z)
!    USES bilinear
DIMENSION rlgk(57,19),T6(57),rlgR(19)

DATA ((rlgk(i,j),j=1,19),i=1,8)&
!    T6=0.00075
/-12.43d0,-11.91d0,-11.39d0,-10.86d0,-10.36d0,&
-9.804d0,-9.350d0,-8.847d0,-8.346d0,-7.845d0,&
-7.345d0,-6.845d0,-6.344d0,-5.843d0,-5.343d0,&
-4.843d0,-4.343d0,-3.843d0,-3.343d0, &

!    T6=0.00100
-11.85d0,-11.35d0,-10.85d0,-10.35d0,-9.854d0,&
-9.354d0,-8.854d0,-8.354d0,-7.854d0,-7.354d0,&
-6.854d0,-6.354d0,-5.854d0,-5.354d0,-4.854d0,&
-4.354d0,-3.854d0,-3.354d0,-2.854d0,&

!    T6=0.00125
-11.46d0,-10.96d0,-10.46d0,-9.964d0,-9.464d0,&
-8.964d0,-8.464d0,-7.964d0,-7.464d0,-6.964d0,&
-6.464d0,-5.964d0,-5.464d0,-4.964d0,-4.464d0,&
-3.964d0,-3.464d0,-2.964d0,-2.464d0,&

!    T6=0.00150
-11.26d0,-10.74d0,-10.22d0,-9.704d0,-9.194d0,&
-8.634d0,-8.180d0,-7.677d0,-7.176d0,-6.675d0,&
-6.175d0,-5.675d0,-5.174d0,-4.673d0,-4.173d0,&
-3.673d0,-3.173d0,-2.673d0,-2.173d0,&

!    T6=0.00175
-11.86d0,-11.12d0,-10.42d0,-9.757d0,-9.139d0,&
-8.557d0,-8.018d0,-7.491d0,-6.979d0,-6.471d0,&
-5.967d0,-5.465d0,-4.962d0,-4.461d0,-3.961d0,&
-3.461d0,-2.961d0,-2.461d0,-1.961d0,&

!    T6=0.00200
-13.67d0,-12.51d0,-11.39d0,-10.43d0,-9.514d0,&
-8.757d0,-8.039d0,-7.433d0,-6.854d0,-6.318d0,&
-5.792d0,-5.281d0,-4.772d0,-4.269d0,-3.767d0,&
-3.266d0,-2.765d0,-2.264d0,-1.764d0,&

!    T6=0.00225
-12.58d0,-12.14d0,-11.67d0,-10.99d0,-10.27d0,&
-9.348d0,-8.440d0,-7.672d0,-6.923d0,-6.322d0,&
-5.727d0,-5.193d0,-4.661d0,-4.150d0,-3.640d0,&
-3.137d0,-2.634d0,-2.133d0,-1.632d0,&

!    T6=0.00250
-10.62d0,-10.50d0,-10.32d0,-10.07d0,-9.761d0,&
-9.375d0,-8.827d0,-8.025d0,-7.238d0,-6.478d0,&
-5.780d0,-5.178d0,-4.602d0,-4.067d0,-3.542d0,&
-3.032d0,-2.524d0,-2.020d0,-1.518d0/

DATA ((rlgk(i,j),j=1,19),i=9,16)&
!    T6=0.00275
/-9.597d0,-9.352d0,-9.106d0,-8.857d0,-8.606d0,&
-8.354d0,-8.075d0,-7.781d0,-7.270d0,-6.635d0,&
-5.954d0,-5.245d0,-4.617d0,-4.036d0,-3.489d0,&
-2.961d0,-2.446d0,-1.937d0,-1.432d0,&

!    T6=0.00300
-8.520d0,-8.331d0,-8.108d0,-7.879d0,-7.615d0,&
-7.345d0,-7.093d0,-6.844d0,-6.582d0,-6.318d0,&
-5.827d0,-5.300d0,-4.677d0,-4.038d0,-3.472d0,&
-2.918d0,-2.394d0,-1.875d0,-1.367d0,&

!    T6=0.00325
-7.730d0,-7.481d0,-7.232d0,-6.982d0,-6.732d0,&
-6.483d0,-6.233d0,-5.984d0,-5.741d0,-5.498d0,&
-5.264d0,-5.016d0,-4.573d0,-4.121d0,-3.559d0,&
-3.000d0,-2.473d0,-1.947d0,-1.436d0,&

!    T6=0.00350
-7.011d0,-6.747d0,-6.498d0,-6.255d0,-6.027d0,&
-5.788d0,-5.516d0,-5.249d0,-5.001d0,-4.758d0,&
-4.531d0,-4.309d0,-4.101d0,-3.832d0,-3.395d0,&
-2.937d0,-2.414d0,-1.894d0,-1.383d0,&

!    T6=0.00375
-6.359d0,-6.111d0,-5.862d0,-5.612d0,-5.362d0,&
-5.112d0,-4.862d0,-4.613d0,-4.364d0,-4.119d0,&
-3.879d0,-3.661d0,-3.474d0,-3.272d0,-3.051d0,&
-2.732d0,-2.291d0,-1.825d0,-1.328d0,&

!    T6=0.00400
-5.798d0,-5.552d0,-5.305d0,-5.057d0,-4.809d0,&
-4.560d0,-4.312d0,-4.061d0,-3.808d0,-3.560d0,&
-3.314d0,-3.089d0,-2.876d0,-2.704d0,-2.559d0,&
-2.329d0,-2.046d0,-1.667d0,-1.228d0,&

!    T6=0.00425
-5.288d0,-5.048d0,-4.807d0,-4.563d0,-4.317d0,&
-4.068d0,-3.819d0,-3.569d0,-3.319d0,-3.070d0,&
-2.822d0,-2.587d0,-2.357d0,-2.183d0,-2.027d0,&
-1.858d0,-1.684d0,-1.385d0,-1.048d0,&

!    T6=0.00450
-4.810d0,-4.585d0,-4.359d0,-4.119d0,-3.877d0,&
-3.630d0,-3.382d0,-3.133d0,-2.883d0,-2.634d0,&
-2.385d0,-2.145d0,-1.905d0,-1.718d0,-1.536d0,&
-1.400d0,-1.267d0,-1.034d0,-0.792d0/

DATA ((rlgk(i,j),j=1,19),i=17,24)&
!    T6=0.00475
/-4.344d0,-4.143d0,-3.941d0,-3.710d0,-3.479d0,&
-3.235d0,-2.991d0,-2.742d0,-2.494d0,-2.244d0,&
-1.995d0,-1.751d0,-1.509d0,-1.304d0,-1.104d0,&
-0.975d0,-0.843d0,-0.661d0,-0.472d0,&

!    T6=0.00500
-3.875d0,-3.713d0,-3.541d0,-3.328d0,-3.111d0,&
-2.875d0,-2.636d0,-2.390d0,-2.143d0,-1.893d0,&
-1.645d0,-1.399d0,-1.157d0,-0.937d0,-0.732d0,&
-0.590d0,-0.448d0,-0.300d0,-0.133d0,&

!    T6=0.00525
-3.410d0,-3.288d0,-3.147d0,-2.961d0,-2.764d0,&
-2.539d0,-2.309d0,-2.067d0,-1.824d0,-1.576d0,&
-1.329d0,-1.081d0,-0.839d0,-0.608d0,-0.400d0,&
-0.240d0,-0.091d0, 0.037d0, 0.191d0,&

!    T6=0.00550
-2.958d0,-2.872d0,-2.757d0,-2.606d0,-2.430d0,&
-2.224d0,-2.005d0,-1.771d0,-1.532d0,-1.286d0,&
-1.039d0,-0.792d0,-0.549d0,-0.312d0,-0.101d0,&
0.077d0, 0.229d0, 0.347d0, 0.496d0,&

!    T6=0.00575
-2.529d0,-2.467d0,-2.374d0,-2.257d0,-2.103d0,&
-1.922d0,-1.717d0,-1.493d0,-1.259d0,-1.017d0,&
-0.773d0,-0.527d0,-0.283d0,-0.042d0, 0.173d0,&
0.366d0, 0.516d0, 0.632d0, 0.780d0,&

!    T6=0.00600
-2.130d0,-2.080d0,-2.004d0,-1.916d0,-1.784d0,&
-1.631d0,-1.441d0,-1.234d0,-1.007d0,-0.771d0,&
-0.529d0,-0.284d0,-0.041d0, 0.201d0, 0.420d0,&
0.628d0, 0.776d0, 0.896d0, 1.042d0,&

!    T6=0.00625
-1.766d0,-1.717d0,-1.653d0,-1.585d0,-1.473d0,&
-1.349d0,-1.174d0,-0.986d0,-0.768d0,-0.542d0,&
-0.303d0,-0.060d0, 0.183d0, 0.425d0, 0.649d0,&
0.866d0, 1.014d0, 1.142d0, 1.285d0,&

!    T6=0.00650
-1.437d0,-1.380d0,-1.324d0,-1.269d0,-1.174d0,&
-1.074d0,-0.916d0,-0.751d0,-0.542d0,-0.327d0,&
-0.091d0, 0.147d0, 0.390d0, 0.633d0, 0.859d0,&
1.083d0, 1.230d0, 1.367d0, 1.507d0/

DATA ((rlgk(i,j),j=1,19),i=25,32)&
!    T6=0.00675
/-1.149d0,-1.070d0,-1.019d0,-0.968d0,-0.887d0,&
-0.806d0,-0.665d0,-0.523d0,-0.324d0,-0.123d0,&
0.108d0, 0.341d0, 0.583d0, 0.825d0, 1.055d0,&
1.284d0, 1.433d0, 1.580d0, 1.715d0,&

!    T6=0.00700
-0.900d0,-0.792d0,-0.737d0,-0.681d0,-0.615d0,&
-0.544d0,-0.424d0,-0.299d0,-0.115d0, 0.073d0,&
0.298d0, 0.523d0, 0.763d0, 1.002d0, 1.235d0,&
1.463d0, 1.620d0, 1.775d0, 1.907d0,&

!    T6=0.00800
-0.224d0,-0.120d0,-0.033d0, 0.037d0, 0.098d0,&
0.189d0, 0.286d0, 0.402d0, 0.538d0, 0.684d0,&
0.842d0, 1.019d0, 1.213d0, 1.405d0, 1.626d0,&
1.850d0, 2.077d0, 2.307d0, 2.540d0,&

!    T6=0.00900
-0.096d0, 0.191d0, 0.434d0, 0.632d0, 0.785d0,&
0.893d0, 0.987d0, 1.084d0, 1.181d0, 1.290d0,&
1.410d0, 1.556d0, 1.721d0, 1.875d0, 2.080d0,&
2.301d0, 2.510d0, 2.679d0, 2.780d0,&

!    T6=0.01000
-0.107d0, 0.234d0, 0.579d0, 0.929d0, 1.196d0,&
1.392d0, 1.527d0, 1.624d0, 1.708d0, 1.799d0,&
1.901d0, 2.020d0, 2.156d0, 2.270d0, 2.458d0,&
2.664d0, 2.861d0, 3.022d0, 3.120d0,&

!    T6=0.01100
-0.145d0, 0.193d0, 0.572d0, 0.991d0, 1.388d0,&
1.681d0, 1.885d0, 2.026d0, 2.118d0, 2.202d0,&
2.294d0, 2.391d0, 2.521d0, 2.606d0, 2.778d0,&
2.973d0, 3.160d0, 3.308d0, 3.386d0,&

!    T6=0.01200
-0.180d0, 0.139d0, 0.519d0, 0.960d0, 1.415d0,&
1.810d0, 2.106d0, 2.297d0, 2.437d0, 2.537d0,&
2.632d0, 2.733d0, 2.848d0, 2.905d0, 3.061d0,&
3.237d0, 3.422d0, 3.605d0, 3.775d0,&

!    T6=0.01400
-0.205d0, 0.048d0, 0.391d0, 0.825d0, 1.314d0,&
1.806d0, 2.241d0, 2.585d0, 2.837d0, 3.012d0,&
3.139d0, 3.259d0, 3.378d0, 3.386d0, 3.527d0,&
3.684d0, 3.853d0, 4.030d0, 4.211d0/

DATA ((rlgk(i,j),j=1,19),i=33,40)&
!    T6=0.01600
/-0.229d0, 0.018d0, 0.341d0, 0.741d0, 1.205d0,&
1.726d0, 2.218d0, 2.671d0, 3.041d0, 3.307d0,&
3.499d0, 3.651d0, 3.786d0, 3.762d0, 3.896d0,&
4.045d0, 4.197d0, 4.340d0, 4.462d0,&

!    T6=0.01800
-0.251d0,-0.024d0, 0.296d0, 0.709d0, 1.186d0,&
1.679d0, 2.175d0, 2.681d0, 3.135d0, 3.490d0,&
3.755d0, 3.949d0, 4.110d0, 4.076d0, 4.204d0,&
4.354d0, 4.495d0, 4.596d0, 4.626d0,&

!    T6=0.02000
-0.258d0,-0.054d0, 0.260d0, 0.684d0, 1.160d0,&
1.657d0, 2.170d0, 2.685d0, 3.182d0, 3.612d0,&
3.939d0, 4.180d0, 4.375d0, 4.351d0, 4.465d0,&
4.590d0, 4.727d0, 4.877d0, 5.041d0,&

!    T6=0.02500
-0.229d0,-0.046d0, 0.255d0, 0.676d0, 1.166d0,&
1.695d0, 2.231d0, 2.765d0, 3.300d0, 3.808d0,&
4.238d0, 4.582d0, 4.836d0, 4.831d0, 4.958d0,&
5.046d0, 5.097d0, 5.113d0, 5.096d0,&

!    T6=0.03000
-0.168d0, 0.014d0, 0.318d0, 0.743d0, 1.245d0,&
1.772d0, 2.317d0, 2.869d0, 3.440d0, 3.971d0,&
4.442d0, 4.832d0, 5.138d0, 5.101d0, 5.206d0,&
5.206d0, 5.206d0, 5.206d0, 5.206d0,&

!    T6=0.03500
-0.078d0, 0.124d0, 0.432d0, 0.845d0, 1.341d0,&
1.874d0, 2.420d0, 2.977d0, 3.552d0, 4.100d0,&
4.589d0, 4.995d0, 5.312d0, 5.250d0, 5.426d0,&
5.426d0, 5.426d0, 5.426d0, 5.426d0,&

!    T6=0.04000
-0.064d0, 0.221d0, 0.571d0, 0.985d0, 1.461d0,&
1.980d0, 2.516d0, 3.068d0, 3.631d0, 4.175d0,&
4.666d0, 5.073d0, 5.384d0, 5.392d0, 5.555d0,&
5.555d0, 5.555d0, 5.555d0, 5.555d0,&

!    T6=0.04500
-0.101d0, 0.211d0, 0.606d0, 1.084d0, 1.565d0,&
2.065d0, 2.583d0, 3.120d0, 3.664d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0/

DATA ((rlgk(i,j),j=1,19),i=41,48)&
!    T6=0.05000
/-0.138d0, 0.149d0, 0.541d0, 1.038d0, 1.568d0,&
2.091d0, 2.604d0, 3.125d0, 3.649d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.05500
-0.170d0, 0.099d0, 0.469d0, 0.941d0, 1.476d0,&
2.030d0, 2.570d0, 3.091d0, 3.603d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.06000
-0.198d0, 0.059d0, 0.409d0, 0.854d0, 1.376d0,&
1.925d0, 2.489d0, 3.032d0, 3.540d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.07000
-0.230d0,-0.011d0, 0.311d0, 0.736d0, 1.242d0,&
1.761d0, 2.329d0, 2.895d0, 3.431d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.08000
-0.242d0,-0.044d0, 0.260d0, 0.670d0, 1.169d0,&
1.685d0, 2.239d0, 2.799d0, 3.354d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.09000
-0.240d0,-0.055d0, 0.241d0, 0.649d0, 1.140d0,&
1.655d0, 2.202d0, 2.758d0, 3.316d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.10000
-0.240d0,-0.050d0, 0.243d0, 0.639d0, 1.128d0,&
1.652d0, 2.196d0, 2.749d0, 3.305d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.12000
-0.246d0,-0.045d0, 0.258d0, 0.661d0, 1.145d0,&
1.672d0, 2.210d0, 2.754d0, 3.292d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0/

DATA ((rlgk(i,j),j=1,19),i=49,57)&
!    T6=0.15000
/-0.286d0,-0.078d0, 0.229d0, 0.634d0, 1.112d0,&
1.621d0, 2.135d0, 2.646d0, 3.141d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.20000
-0.350d0,-0.196d0, 0.054d0, 0.400d0, 0.818d0,&
1.280d0, 1.758d0, 2.236d0, 2.700d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.25000
-0.387d0,-0.281d0,-0.090d0, 0.184d0, 0.539d0,&
0.957d0, 1.413d0, 1.881d0, 2.340d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.30000
-0.409d0,-0.332d0,-0.182d0, 0.041d0, 0.347d0,&
0.729d0, 1.163d0, 1.620d0, 2.077d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.40000
-0.431d0,-0.384d0,-0.280d0,-0.118d0, 0.123d0,&
0.446d0, 0.839d0, 1.275d0, 1.724d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.50000
-0.442d0,-0.409d0,-0.329d0,-0.200d0, 0.000d0,&
0.282d0, 0.641d0, 1.055d0, 1.493d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.60000
-0.449d0,-0.424d0,-0.357d0,-0.249d0,-0.074d0,&
0.179d0, 0.511d0, 0.905d0, 1.333d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=0.80000
-0.456d0,-0.439d0,-0.389d0,-0.305d0,-0.164d0,&
0.049d0, 0.341d0, 0.704d0, 1.115d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0,&

!    T6=1.00000
-0.460d0,-0.447d0,-0.405d0,-0.335d0,-0.213d0,&
-0.026d0, 0.238d0, 0.577d0, 0.972d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0, 0.000d0,&
0.000d0, 0.000d0, 0.000d0, 0.000d0/

DATA (T6(i),i=1,57)&
     /0.75d-3,1.00d-3,1.25d-3,1.50d-3,1.75d-3,2.00d-3,&
2.25d-3,2.50d-3,2.75d-3,3.00d-3,3.25d-3,&
3.50d-3,3.75d-3,4.00d-3,4.25d-3,4.50d-3,&
4.75d-3,5.00d-3,5.25d-3,5.50d-3,5.75d-3,&
6.00d-3,6.25d-3,6.50d-3,6.75d-3,7.00d-3,&
8.00d-3,9.00d-3,1.00d-2,1.10d-2,1.20d-2,&
1.40d-2,1.60d-2,1.80d-2,2.00d-2,2.50d-2,&
3.00d-2,3.50d-2,4.00d-2,4.50d-2,5.00d-2,&
5.50d-2,6.00d-2,7.00d-2,8.00d-2,9.00d-2,&
1.00d-1,1.20d-1,1.50d-1,2.00d-1,2.50d-1,&
3.00d-1,4.00d-1,5.00d-1,6.00d-1,8.00d-1,1.00d+0/

DATA (rlgR(i),i=1,19)&
/-5.0d+0,-4.5d+0,-4.0d+0,-3.5d+0,-3.0d+0,&
-2.5d+0,-2.0d+0,-1.5d+0,-1.0d+0,-0.5d+0,&
0.0d+0, 0.5d+0, 1.0d+0, 1.5d+0,2.0d+0,&
2.5d+0, 3.0d+0, 3.5d+0, 4.0d+0/

temp6=T/1.d6
rlogR=log10(rho/temp6**3)
 
if((temp6.ge.0.75d-3).and.(temp6.le.2.5d-2)&
.and.(rlogR.ge.-5.d0).and.(rlogR.le.4.d0)) then
   call bilinear(T6,rlgR,rlgk,57,19,temp6,rlogR,rlogk)
   xk_prm=1.d1**rlogk
elseif((temp6.ge.2.5d-2).and.(temp6.le.4.d-2)&
.and.(rlogR.ge.-5.d0).and.(rlogR.le.3.d0)) then
   call bilinear(T6,rlgR,rlgk,57,19,temp6,rlogR,rlogk)
   xk_prm=1.d1**rlogk
elseif((temp6.gt.4.d-2)&
.and.(rlogR.ge.-5.d0).and.(rlogR.le.-1.d0)) then
   call bilinear(T6,rlgR,rlgk,57,19,temp6,rlogR,rlogk)
   xk_prm=1.d1**rlogk
elseif((temp6.lt.0.75d-3).or.(rlogR.lt.-5.d0)) then
   xk_prm=0.d0
else
   xk_prm=1.d5
endif

return
END



SUBROUTINE linear(xa,ya,m,x,y)
IMPLICIT REAL*8(a-h,o-z)
integer :: m
double precision :: xa(m),ya(m),x,y
do 11 i=1,m
   if(x-xa(i).le.0.d0) then
      ms=i
      go to 12
   endif
11   continue
ms=m
12   continue
if(ms.eq.1) ms=2
y1=ya(ms-1)
y2=ya(ms)
t=(x-xa(ms-1))/(xa(ms)-xa(ms-1))
y=(1.d0-t)*y1+t*y2
return
END


SUBROUTINE gaussj(a,n,np,b,m,mp)
IMPLICIT REAL*8(a-h,o-z)
INTEGER m,mp,n,np
!     PARAMETER (NMAX=50,eps=1.d-13)
integer, PARAMETER :: NMAX=50
double precision, PARAMETER :: eps=0.d0
DOUBLE PRECISION a(np,np),b(np,mp),b_max(np,mp)
INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
DOUBLE PRECISION big,dum,pivinv
do 11 j=1,n
  ipiv(j)=0
11    continue

do ll=1,np
   do l=1,mp
      b_max(ll,l)=0.d0
   enddo
enddo

do 22 i=1,n
  big=0.d0
  do 13 j=1,n
    if(ipiv(j).ne.1)then
      do 12 k=1,n
        if (ipiv(k).eq.0) then
          if (abs(a(j,k)).ge.big)then
            big=abs(a(j,k))
            irow=j
            icol=k
          endif
        else if (ipiv(k).gt.1) then
          print *,'singular matrix in gaussj'
          go to 25
        endif
12          continue
    endif
13      continue
  ipiv(icol)=ipiv(icol)+1
  if (irow.ne.icol) then
    do 14 l=1,n
      dum=a(irow,l)
      a(irow,l)=a(icol,l)
      a(icol,l)=dum
14        continue
    do 15 l=1,m
      dum=b(irow,l)
      b(irow,l)=b(icol,l)
      b(icol,l)=dum
15        continue
  endif
  indxr(i)=irow
  indxc(i)=icol
  if (a(icol,icol).eq.0.d0) then
     !print *, 'singular matrix in gaussj'
     go to 25
  endif
  pivinv=1.d0/a(icol,icol)
  a(icol,icol)=1.d0
  do 16 l=1,n
    a(icol,l)=a(icol,l)*pivinv
16      continue
  do 17 l=1,m
    b(icol,l)=b(icol,l)*pivinv
17      continue

  do 21 ll=1,n
     if(ll.ne.icol)then
        dum=a(ll,icol)
        a(ll,icol)=0.d0
        do 18 l=1,n
           a(ll,l)=a(ll,l)-a(icol,l)*dum
           if(a(icol,l)*dum+a(ll,l).ne.a(ll,l)) then
              if(dabs(a(ll,l)/(a(icol,l)*dum)).lt.eps) then
                 a(ll,l)=0.d0
              endif
           endif
18           continue
        do 19 l=1,m
           b(ll,l)=b(ll,l)-b(icol,l)*dum
           if(b(icol,l)*dum+b(ll,l).ne.b(ll,l)) then
              if(dabs(b(ll,l)/(b(icol,l)*dum)).lt.eps) then
                 b(ll,l)=0.d0
              endif
           endif
           b_max(ll,l)=max(b_max(ll,l),dabs(b(icol,l)*dum))
19           continue
     endif
21     continue
22   continue

do ll=1,np
   do l=1,mp
      if(b_max(ll,l).ne.0.d0) then
         if(dabs(b(ll,l)/b_max(ll,l)).le.eps) then
            b(ll,l)=0.d0
         endif
      endif
   enddo
enddo


do 24 l=n,1,-1
  if(indxr(l).ne.indxc(l))then
    do 23 k=1,n
      dum=a(k,indxr(l))
      a(k,indxr(l))=a(k,indxc(l))
      a(k,indxc(l))=dum
23        continue
  endif
24    continue
25    continue
return
END
! (C) Copr. 1986-92 Numerical Recipes Software v%1jw#<?4210(93Y"+91.d0

SUBROUTINE bilinear(x1a,x2a,ya,m,n,x1,x2,y)
IMPLICIT REAL*8(a-h,o-z)
integer m,n
double precision :: x1a(m),x2a(n),ya(m,n),x1,x2,y
do 11 i=1,m
   if(x1-x1a(i).le.0.d0) then
      ms=i
      go to 12
   endif
11   continue
ms=m
12   continue
do 13 i=1,n
   if(x2-x2a(i).le.0.d0) then
      ns=i
      go to 14
   endif
13   continue
ns=n
14   continue
if(ms.eq.1) ms=2
if(ns.eq.1) ns=2
y1=ya(ms-1,ns-1)
y2=ya(ms,ns-1)
y3=ya(ms,ns)
y4=ya(ms-1,ns)
t=(x1-x1a(ms-1))/(x1a(ms)-x1a(ms-1))
u=(x2-x2a(ns-1))/(x2a(ns)-x2a(ns-1))
y=(1.d0-t)*(1.d0-u)*y1+t*(1.d0-u)*y2+t*u*y3+(1.d0-t)*u*y4
return
END




