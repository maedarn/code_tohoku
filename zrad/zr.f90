MODULE COMVAR
integer, parameter :: N_sp=50
double precision :: zeta, G_0
double precision :: T_rad
double precision :: A_v,xNc_H,xNc_H2,xNc_HD
double precision :: y(N_sp)
double precision ::xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
double precision :: xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,&
xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
double precision :: T_gr_K,tau_cont
double precision :: xm_p=1.67d-24,xk_B=1.38d-16,G=6.67d-8,&
pi=3.14159265358979d0,yHe=8.333d-2
double precision ::A0,DT0,sigma,eta_T,xmu_mol,J_max
double precision :: xLmbd_net
END MODULE COMVAR

PROGRAM ZR
!computes the evolutionary path of
!metal polluted protostellar clouds
!with radiation
!updated 2008 Jun
USE COMVAR
integer :: i_ev=1,itchem,i
double precision :: Z_gas=1.0d0,Z_dust=1.0d0,G_0=1.d0,zeta=1.d-17,T_rad=10.d0,Z_metal
double precision :: yHe=8.333d-2,yD=3.d-5
double precision :: xnH,rho,T_K,y_Hp,y_H2,y_Dp,y_HD,y_H,y_D,yC,yO,y_C,y_Cp
double precision :: xmu,gamma,e,T_gr_K,xLmbd_ch
double precision :: t,dt,t_chem,esc_cnt
double precision :: xlmbd_J,radius,xNcol,xNcolJ
double precision :: t_col,v_bulk,v_turb,v_bulk,v_D,xNc_H2,xNc_H,xNc_HD,xMJ

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
T_K=10.d0
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

xmu=(1.d0+4.d0*yHe)
&     /(y(1)+y(2)+y(3)+y(4)+y(8)+y(9)+y(10))

gamma=1.d0+(1.d0+4.d0*yHe)
&     /(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))
&     +c_H2(T_K)*y(2)))

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
!xlmbd_J=dsqrt(pi*xk_B*T_K/(G*xmu*xm_p*rho))
!radius=xlmbd_J
!xNcolJ=xnH*xlmbd_J
!2) constant column density
xNcol=1.d19
radius=xNcol/xnH

t_col=dsqrt(3.d0*pi/(32.d0*G*rho))
v_bulk=radius/3.d0/t_col
v_turb=0.d0
v_bulk=dsqrt(v_bulk**2+v_turb**2)
!continuum cooling
A_v=Z_metal*xnH*radius*5.3d-22
call  rad_cool(Z_metal,T_K,T_gr_K,radius,A_v,esc_cnt,&
dt,xmu,gamma,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,&
i,i_ev)
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

call chemcool(xnH,T_K,T_gr_K,Z_metal,&
y,dt,t_chem,xmu,gamma,xLmbd_ch)

y_H=y(1)
y_H2=y(2)
y_e=y(3)
y_HD=y(13)

xLmbd_net=xLmbd_cnt+xLmbd_line+xLmbd_gr+xLmbd_ch &
-Gmm_pe-Gmm_CR

xMJ=rho*radius**3/2.d33

!update values of rho & e
call update(i_ev)

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
      if(dmin1(2.d-2*t_col,2.d-2*t_cool).le.dt) itchem=1
   elseif(xnH < 1.d0) then
      dt=dmin1(2.d-2*t_col,2.d-2*t_cool)
   elseif(xnH < 1.d14) then
      dt=dmin1(2.d-2*t_col,2.d-2*t_cool)
   else
      dt=dmin1(5.d-2*t_col,5.d-2*t_cool)
   endif
!         print *, itchem, 'dt=', t_chem,
!     &        dmin1(2.d-2*t_col,2.d-2*t_cool)
!         pause
enddo
end PROGRAM ZR


SUBROUTINE update(rho,e,p,T_K,xnH,t_col,dt,i_ev)
USE COMVAR
integer :: i_ev,it
double precision :: rho,e,p,rhoo,eo,po,xnH,t_col,dt,T_K
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
double precision :: gamma_eff,ft
double precision :: f

if(gamma_eff < 0.83d0) then
   f=0.d0
elseif(gamma_eff < 1.d0) then
   f=0.6d0+2.5d0*(gamma_eff-1.d0)-6.d0*(gamma_eff-1.d0)**2
elseif(gamma_eff < 1.33333d0) then
   f=1.d0+0.2d0*(gamma_eff-4.d0/3.d0)
&        -2.9d0*(gamma_eff-4.d0/3.d0)**2
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
double precision, parameter :: B_ex=5.d-4
double precision :: T_vis,delta,x,T_nu
double precision :: Q_bg_CMB,Q_bg_vis,Q_bg_gr

T_vis=6.d3
delta=4.41d3*(B_ex/T_vis**4)
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


SUBROUTINE rad_cool(Z_metal,Tp,Tp_gr,
&     radius,A_v,esc_cnt,dt,xmu,gamma,
&     xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,i,i_ev)
USE COMVAR
integer, parameter :: maxit=100
integer :: i,i_ev
double precision :: Z_metal,Tp,Tp_gr,radius,A_v,esc_cnt,dt,xmu,gamma,y,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr
double precision :: Tp0,Tp1,fL,f,Tpsec,TpL,swap,xLmbd_cnt1,xLmbd_gr


func(t,xLmbd)=Tp-t-(gamma-1.d0)*xmu*xm_p*xLmbd*dt/xk_B !dif, function

Tp0=Tp
Tp1=Tp
if(tau_cnt < 10.d0) then
   T_K=Tp1
   call line_cool(radius,xLmbd_line,i)
else
   xLmbd_line=0.d0
endif
call cnt(Z_metal,xnH,Tp1,Tp_gr,
&     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
if(dabs(Tp1/T_rad-1.d0) > 1.d0) return
xLmbd=xLmbd_ch+xLmbd_line+xLmbd_cnt+xLmbd_gr
fL=func(Tp1,xLmbd)

Tp2=1.01d0*Tp1
call cnt(Z_metal,xnH,Tp2,Tp_gr,
&     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
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

do it=1,maxit
   dTp=(TpL-Tpsec)*f/(f-fL)
   TpL=Tpsec
   fL=f
   Tpsec=Tpsec+dTp
   call cnt(Z_metal,xnH,Tpsec,Tp_gr,
&     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt1,xLmbd_gr1)
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


SUBROUTINE cnt(Z_metal,xnH,T_K,T_gr_K,
&     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
USE COMVAR
!      PARAMETER (xk_vis=40.d0,B_ex=5.d-4)
!      PARAMETER (xk_vis=40.d0,B_ex=0.d0)

yHe=8.333d-2
!      xJ_ex=B_ex*dexp(-A_v/2.5d0)
xJ_ex=0.d0
call grtemp(xnH,T_K,esc_cnt,T_rad,xJ_ex,T_gr_K)
rho=(1.d0+4.d0*yHe)*xm_p*xnH
xk_gr=xkp_gr(rho,T_gr_K)*Z_metal
xk_gas=xk_prm(T_K,rho)
xk_cnt=xk_gas+xk_gr
tau_cnt=xk_cnt*rho*radius
if(tau_cnt.gt.1.d0) then
   esc_cnt=1.d0/(tau_cnt**2)
else
   esc_cnt=1.d0
endif
xLmbd_cnt=4.d0*sigma_B*(T_K**4-T_rad**4)*xk_gas*esc_cnt
xLmbd_gr=4.d0*sigma_B*(T_gr_K**4-T_rad**4)*xk_gr*esc_cnt
!     &     -xk_vis*Z_metal*xJ_ex
END SUBROUTINE cnt

SUBROUTINE line_cool(radius,xLmbd_line,i)
USE COMVAR
integer, parameter :: N_p=100
double precision :: func(N_p),esc_CI(N_p),esc_CII(N_p),esc_OI(N_p),&
esc_CImeta(N_p),esc_CIImeta(N_p),esc_OImeta(N_p)
EXTERNAL pop_CI,pop_OI,pop_CII,
&     pop_CImeta,pop_OImeta,pop_CIImeta
DATA xm_p/1.67d-24/,xk_B/1.38d-16/,pi/3.14159265358979d0/
SAVE esc_CI,esc_CII,esc_OI,
&     esc_CImeta,esc_CIImeta,esc_OImeta

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
xLmbd_H2=y_m
&     *xLd_H2(xnH,T_K,y_a,y_m,y_e,y_Hp,xNc_H2,tau_cnt)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)

!     HD cooling
v_D=v_th/dsqrt(3.d0)
xNc_HD=dmin1(1.d0,v_D/v_bulk)*xNc_HD
xLmbd_HD=y_HD*xLd_HD(xnH,T_K,y_a,y_m,xNc_HD,tau_cnt)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)

!     CO cooling
xmu_mol=28.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CO
call COcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_CO)
call COcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_CO)
xLmbd_CO=y_CO
&     *(xLd_CO-xLdr_CO)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)

!     OH cooling
xmu_mol=17.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OH
call OHcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_OH)
call OHcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_OH)
xLmbd_OH=y_OH
&     *(xLd_OH-xLdr_OH)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)

!     H2O cooling
xmu_mol=18.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_H2O
call H2Ocool(xnH,T_K,y_m,xNc,tau_cnt,xLd_H2O)
call H2Ocool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_H2O)
xLmbd_H2O=y_H2O
&     *(xLd_H2O-xLdr_H2O)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)

!     CII cooling
xmu_C=12.d0
xmu_O=16.d0
v_D=v_th/dsqrt(xmu_C)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CII
!     fine structure line
xLmbd_CII=y_CII
&     *xLd(pop_CII,esc_CII,1)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)
!     metastable line
if(T_K > 2.d2) then
   xLmbd_CIImeta=y_CII
&        *xLd(pop_CIImeta,esc_CIImeta,1)
&        *xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_CII=xLmbd_CII+xLmbd_CIImeta
endif

!     CI cooling
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CI
!     fine structure line
xLmbd_CI=y_CI
&     *xLd(pop_CI,esc_CI,3)
&     *xnH/((1.d0+4.d0*yHe)*xm_p)

if(T_K > 3.d3) then
   xLmbd_CImeta=y_CI
&        *xLd(pop_CImeta,esc_CImeta,3)
&        *xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_CI=xLmbd_CI+xLmbd_CImeta
endif

!     OI cooling
v_D=v_th/dsqrt(xmu_O)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OI
!     fine structure line
if(T_K > 15.d0) then
   xLmbd_OI=y_OI
&        *xLd(pop_OI,esc_OI,3)
&        *xnH/((1.d0+4.d0*yHe)*xm_p)
else
   xLmbd_OI=0.d0
endif
!     metastable line
if(T_K > 3.5d3) then
   xLmbd_OImeta=y_OI
&        *xLd(pop_OImeta,esc_OImeta,3)
&        *xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_OI=xLmbd_OI+xLmbd_OImeta
endif

!     HI Lya cooling
if(xnH < 1.d5) then
   T_5=1.d-5*T_K
   xLmbd_Lya=y_e*y_a*7.50d-19/(1.d0+dsqrt(T_5))
&        *dexp(-1.18348d5/T_K)
&        *xnH/((1.d0+4.d0*yHe)*xm_p)
else
   xLmbd_Lya=0.d0
endif
!     total cooling rate by lines
xLmbd_line=xLmbd_H2+xLmbd_HD+xLmbd_CO+xLmbd_OH+xLmbd_H2O
&     +xLmbd_CII+xLmbd_CI+xLmbd_OI+xLmbd_Lya
END SUBROUTINE line_cool


include "subs/chemistry.f"
include "subs/reaction.f"  ! full reactions
include "subs/grain.f"
include "subs/xk_prm.f"
include "subs/H2.f"
include "subs/HD.f"
include "subs/lines.f"
include "subs/CR.f"
include "subs/math.f"
include "subs/CO.f"
include "subs/H2O.f"
include "subs/OH.f"



SUBROUTINE chemcool(xnH,Tp,Tp_gr,Z_metal,y,dt,tchem,xmu,gamma,xLmbdch)
USE COMVAR
integer, parameter :: maxit=100
double precision :: ytmp(N_sp)
double precision :: Tmp_f,dTp
!     input    Tp,xnH,dt,y,xmu,gamma
!     output   y,xmu,gamma,xLmbdch

func(Tmp_f,xLmbdch)=Tp-Tmp_f-(gamma-1.d0)*xmu*xm_p*xLmbdch*dt/xk_B

Tp1=Tp
do j=1,N_sp
   ytmp(j)=y(j)
enddo
call chemreact(xnH,Tp1,Tp_gr,Z_metal,
&     ytmp,dt,tchem1,xmu1,gamma1,xLmbdch1)
fL=func(Tp1,xLmbdch1)
if((xnH.le.1.d16).and.(Tp.le.1.65d3)) go to 12

Tp2=0.99d0*Tp1
do j=1,N_sp
   ytmp(j)=y(j)
enddo
call chemreact(xnH,Tp2,Tp_gr,Z_metal,
&     ytmp,dt,tchem2,xmu2,gamma2,xLmbdch2)
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
   call chemreact(xnH,Tpsec,Tp_gr,Z_metal,
&        ytmp,dt,tchem1,xmu1,gamma1,xLmbdch1)
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


SUBROUTINE chemreact(xnH,T_K,T_gr_K,Z_metal,y,dt,tchem,xmu,gamma,xLmbdch)
!******************************************
!*      Version 4+      (13 Sep 1998)     *
!*      no radiation                      *
!******************************************
USE COMVAR
integer, parameter :: N_react=675
integer :: isp
double precision, parameter :: eps=1.d-4,eps_y=1.d-10,xnH_eq=3.d16
double precision :: y_init(N_sp),y_tmp(N_sp),dy(N_sp),ddy(N_sp), xk(N_react), r_f(N_sp),r_f_fw(N_sp),r_f_bw(N_sp), dr_fdy(N_sp,N_sp),A(N_sp,N_sp)

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
          dr_fdy(isp,jsp)=(r_f_fw(isp)-r_f_bw(isp))
&              /delta_y
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
      if(y(isp).lt.0.d0) y(isp)=y_init(isp)
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
            tch=dabs((y(isp)+y_init(isp))
&                 /(2.d0*(y(isp)-y_init(isp))))*dt
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

xmu=(1.d0+4.d0*yHe)
&     /(y(1)+y(2)+y(3)+y(4)+y(8)+y(9)+y(10))

gamma=1.d0+(1.d0+4.d0*yHe)
&     /(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))
&     +c_H2(T_K)*y(2)))

xm_p=1.67d-24
!C     H_2 formation/dissociation cooling
if(xnH.lt.1.d13) then
   xn_cr=1.d6/dsqrt(T_K)
&        /(1.6d0*y(1)*dexp(-(4.d2/T_K)**2)
&        +1.4d0*y(2)*dexp(-1.2d4/(T_K+1.2d3)))
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
   rt3b=(xk(19)*(y(1)**3)+xk(20)*(y(1)**2)*y(2)
&        +xk(43)*(y(1)**2)*y(8))*xnH
!C     H2 + H  -> 3 H
!C     2 H2  -> 2 H + H2
!c     H2 + He -> 2 H + He :37
   rtdis=xk(13)*y(2)*y(1)+xk(21)*(y(2)**2)
&        +xk(37)*y(2)*y(8)
!C     H2 + ph. -> 2H
   rtphdis=xk(590)*y(2)/xnH

   xL_H2=-(rtgrain*(0.2d0+4.2d0*crit)
&        +rtHm*3.53d0*crit
&        +rtH2p*1.83d0*crit
&        +(rt3b*crit-rtdis)*4.48d0
&        +rtphdis*0.4d0)
&        *xnH*1.60219d-12
   dyHp=dy(4)
&        +(xk(2)*y(4)*y(3)*xnH-xk(544)*y(1)-xk(548)*y(2))*dt
   dyHep=dy(9)+(xk(4)*y(9)*y(3)*xnH-xk(545)*y(8))*dt
   dyHepp=dy(10)
else
   dyH2=dy(2)
   xL_H2=-7.18d-12*dyH2/dt
   dyHp=dy(4)
   dyHep=dy(9)
   dyHepp=dy(10)
endif
xLmbdch=((2.18d-11*dyHp+3.94d-11*dyHep+12.66d-11*dyHepp)/dt
&     +xL_H2)/((1.d0+4.d0*yHe)*xm_p)

END SUBROUTINE chemreact


SUBROUTINE equichem(xnH,T_K,y_e,y_HI,y_HII,y_H2,y_HeI,y_HeII,y_HeIII)
!***********************************************
!*        Equillibrium  Chemistry    (Saha)    *
!***********************************************
USE COMVAR
double precision, parameter :: yHe=8.333d-2,rmasse=9.109534d-28,h=6.626176d-27,rmassp=1.67d-24
double precision :: T_eV,phyc,chiH0,chiHe0,chiHep,chiH2,zH0,zHe0,zHep,zH2

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
      Fe=y_e**3+fHe0*(y_e**2+(fHep-yHe)*y_e-2.d0*fHep*yHe)
&           -5.d-1*fH0*(y_e+fHe0*(1.d0+fHep/y_e))*d
      if((dabs(Fe)**(1.d0/3.d0))/y_e.le.1.d-4) go to 12
      dFe=3.d0*y_e**2+fHe0*(2.d0*y_e+(fHep-yHe))
&           -5.d-1*fH0*((1.d0-fHe0*fHep/y_e**2)*d
&           +(y_e+fHe0*(1.d0+fHep/y_e))*(-5.d-1*fH2*fH0/y_e**2)
&           *(-1.d0+5.d-1*fH2*(1.d0+fH0/y_e)
&           /dsqrt(b**2+2.d0*fH2)))
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
implicit real*8(a-h,o-z)
data xk_B/1.38d-16/, xm_p/1.67d-24/
real*8 :: xlTa(1:13)
data xlTa/0.477d0, 0.778d0, 1.000d0, 1.301d0, 1.477d0,
&     1.699d0, 1.903d0, 2.000d0, 2.477d0, 2.778d0,
&     3.000d0, 3.176d0, 3.301d0/
real*8 :: xlNa(1:11)
data xlNa/14.0d0, 14.5d0, 15.0d0, 15.5d0, 16.0d0,
&     16.5d0, 17.0d0, 17.5d0, 18.0d0, 18.5d0, 19.0d0/
real*8 :: aL0a(1:13)
data  aL0a/0.2585d2, 0.2517d2, 0.2477d2, 0.2438d2, 0.2421d2,
&     0.2403d2,
&     0.2389d2, 0.2382d2, 0.2342d2, 0.2313d2, 0.2291d2, 0.2263d2,
&     0.2228d2/
real*8 :: aLLTEa(1:13,1:11)
data aLLTEa/0.2251d2, 0.2165d2, 0.2108d2, 0.2035d2, 0.1994d2,
&     0.1945d2,
&     0.1901d2, 0.1880d2, 0.1781d2, 0.1723d2, 0.1686d2, 0.1666d2,
&     0.1655d2, 0.2254d2, 0.2167d2, 0.2109d2, 0.2035d2, 0.1995d2,
&     0.1945d2, 0.1901d2, 0.1880d2, 0.1781d2, 0.1723d2, 0.1686d2,
&     0.1666d2, 0.1655d2, 0.2262d2, 0.2171d2, 0.2111d2, 0.2037d2,
&     0.1996d2, 0.1946d2, 0.1901d2, 0.1880d2, 0.1781d2, 0.1723d2,
&     0.1686d2, 0.1666d2, 0.1655d2, 0.2282d2, 0.2183d2, 0.2118d2,
&     0.2040d2, 0.1998d2, 0.1947d2, 0.1902d2, 0.1881d2, 0.1782d2,
&     0.1723d2, 0.1687d2, 0.1666d2, 0.1655d2, 0.2317d2, 0.2208d2,
&     0.2137d2, 0.2051d2, 0.2005d2, 0.1952d2, 0.1905d2, 0.1883d2,
&     0.1782d2, 0.1723d2, 0.1687d2, 0.1666d2, 0.1655d2, 0.2362d2,
&     0.2244d2, 0.2167d2, 0.2073d2, 0.2023d2, 0.1964d2, 0.1913d2,
&     0.1890d2, 0.1785d2, 0.1725d2, 0.1688d2, 0.1667d2, 0.1656d2,
&     0.2405d2, 0.2284d2, 0.2204d2, 0.2105d2, 0.2052d2, 0.1987d2,
&     0.1932d2, 0.1906d2, 0.1792d2, 0.1728d2, 0.1690d2, 0.1669d2,
&     0.1658d2, 0.2450d2, 0.2326d2, 0.2244d2, 0.2142d2, 0.2086d2,
&     0.2019d2, 0.1960d2, 0.1933d2, 0.1808d2, 0.1738d2, 0.1697d2,
&     0.1675d2, 0.1663d2, 0.2500d2, 0.2373d2, 0.2287d2, 0.2182d2,
&     0.2124d2, 0.2055d2, 0.1995d2, 0.1966d2, 0.1834d2, 0.1759d2,
&     0.1715d2, 0.1691d2, 0.1678d2, 0.2549d2, 0.2417d2, 0.2330d2,
&     0.2223d2, 0.2165d2, 0.2094d2, 0.2032d2, 0.2003d2, 0.1867d2,
&     0.1789d2, 0.1748d2, 0.1726d2, 0.1712d2, 0.2597d2, 0.2463d2,
&     0.2376d2, 0.2266d2, 0.2206d2, 0.2135d2, 0.2071d2, 0.2042d2,
&     0.1903d2, 0.1826d2, 0.1793d2, 0.1774d2, 0.1761d2/
real*8 :: xlnha(1:13,1:11)
data xlnha/0.323d1, 0.326d1, 0.329d1, 0.349d1, 0.367d1, 0.397d1,
&     0.430d1, 0.446d1, 0.517d1, 0.547d1, 0.553d1, 0.530d1,
&     0.470d1, 0.319d1, 0.324d1, 0.327d1, 0.348d1, 0.366d1,
&     0.396d1, 0.430d1, 0.445d1, 0.516d1, 0.547d1, 0.553d1,
&     0.530d1, 0.470d1, 0.305d1, 0.316d1, 0.322d1, 0.345d1,
&     0.364d1, 0.394d1, 0.429d1, 0.445d1, 0.516d1, 0.547d1,
&     0.553d1, 0.530d1, 0.470d1, 0.273d1, 0.297d1, 0.307d1,
&     0.334d1, 0.356d1, 0.389d1, 0.426d1, 0.442d1, 0.515d1,
&     0.546d1, 0.552d1, 0.530d1, 0.470d1, 0.226d1, 0.257d1,
&     0.272d1, 0.309d1, 0.335d1, 0.374d1, 0.416d1, 0.434d1,
&     0.513d1, 0.545d1, 0.551d1, 0.529d1, 0.468d1, 0.176d1,
&     0.207d1, 0.224d1, 0.265d1, 0.295d1, 0.342d1, 0.392d1,
&     0.414d1, 0.506d1, 0.541d1, 0.548d1, 0.526d1, 0.464d1,
&     0.126d1, 0.157d1, 0.174d1, 0.215d1, 0.247d1, 0.295d1,
&     0.349d1, 0.374d1, 0.486d1, 0.530d1, 0.539d1, 0.517d1,
&     0.453d1, 0.760d0, 0.107d1, 0.124d1, 0.165d1, 0.197d1,
&     0.245d1, 0.300d1, 0.325d1, 0.447d1, 0.502d1, 0.516d1,
&     0.494d1, 0.427d1, 0.260d0, 0.573d0, 0.742d0, 0.115d1,
&     0.147d1, 0.195d1, 0.250d1, 0.275d1, 0.398d1, 0.457d1,
&     0.473d1, 0.452d1, 0.384d1,-0.240d0, 0.727d-1, 0.242d0,
&     0.652d0, 0.966d0, 0.145d1, 0.200d1, 0.225d1, 0.348d1,
&     0.407d1, 0.424d1, 0.403d1, 0.335d1,-0.740d0,-0.427d0,
&    -0.258d0, 0.152d0, 0.466d0, 0.954d0, 0.150d1, 0.175d1,
&     0.298d1, 0.357d1, 0.374d1, 0.353d1, 0.285d1/
real*8 :: alphaa(1:13,1:11)
data alphaa/0.514d0, 0.465d0, 0.439d0, 0.409d0, 0.392d0, 0.370d0,
&     0.361d0, 0.357d0, 0.385d0, 0.437d0, 0.428d0, 0.354d0,
&     0.322d0, 0.505d0, 0.460d0, 0.436d0, 0.407d0, 0.391d0,
&     0.368d0, 0.359d0, 0.356d0, 0.385d0, 0.437d0, 0.427d0,
&     0.354d0, 0.322d0, 0.486d0, 0.448d0, 0.428d0, 0.401d0,
&     0.385d0, 0.364d0, 0.356d0, 0.352d0, 0.383d0, 0.436d0,
&     0.427d0, 0.352d0, 0.320d0, 0.465d0, 0.433d0, 0.416d0,
&     0.388d0, 0.373d0, 0.353d0, 0.347d0, 0.345d0, 0.380d0,
&     0.434d0, 0.425d0, 0.349d0, 0.316d0, 0.466d0, 0.448d0,
&     0.416d0, 0.378d0, 0.360d0, 0.338d0, 0.332d0, 0.330d0,
&     0.371d0, 0.429d0, 0.421d0, 0.341d0, 0.307d0, 0.503d0,
&     0.487d0, 0.450d0, 0.396d0, 0.367d0, 0.334d0, 0.322d0,
&     0.317d0, 0.355d0, 0.419d0, 0.414d0, 0.329d0, 0.292d0,
&     0.566d0, 0.538d0, 0.492d0, 0.435d0, 0.403d0, 0.362d0,
&     0.339d0, 0.329d0, 0.343d0, 0.406d0, 0.401d0, 0.317d0,
&     0.276d0, 0.598d0, 0.574d0, 0.529d0, 0.473d0, 0.441d0,
&     0.404d0, 0.381d0, 0.370d0, 0.362d0, 0.410d0, 0.392d0,
&     0.316d0, 0.272d0, 0.603d0, 0.594d0, 0.555d0, 0.503d0,
&     0.473d0, 0.440d0, 0.423d0, 0.414d0, 0.418d0, 0.446d0,
&     0.404d0, 0.335d0, 0.289d0, 0.613d0, 0.623d0, 0.582d0,
&     0.528d0, 0.499d0, 0.469d0, 0.457d0, 0.451d0, 0.470d0,
&     0.487d0, 0.432d0, 0.364d0, 0.310d0, 0.634d0, 0.645d0,
&     0.596d0, 0.546d0, 0.519d0, 0.492d0, 0.483d0, 0.479d0,
&     0.510d0, 0.516d0, 0.448d0, 0.372d0, 0.313d0/


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
xLinv=xL0inv
&     +xn_c*xLLTEinv
&     +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
xL=1.d0/xLinv

xLd_CO=(1.d0-y_H2)*xL*dexp(-tau_cnt)
END SUBROUTINE COcool



FUNCTION Gamma_CR(zeta,y_a,y_m,y_He)
!     calculates heating rate (ergs/s/g) due to CR
IMPLICIT REAL*8(a-h,o-z)
Gamma_CR=3.26d12*(0.46d0*y_a+0.50d0*y_He+0.94d0*y_m)*zeta
&     /(1.d0+4.d0*y_He)
return
END FUNCTION Gamma_CR


SUBROUTINE grtemp(xnH,T,esc_cnt,T_rad,xJ_ex,T_gr)
IMPLICIT REAL*8(a-h,o-z)
PARAMETER (yHe=8.333d-2,xk_vis=40.d0)
DATA xm_p/1.67d-24/
funcL(T_d)=4.9d-2*vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
&     *dsqrt(T*1.d-3)*(0.354+0.5d0*yHe)*(T-T_d)
!      funcL(T_d)=4.9d-2*vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
!     &     *dsqrt(T*1.d-3)*(1.d0-0.8d0*dexp(-75.d0/T))*(T-T_d)
!     for stellar radiation
!      funcD(T_d)=vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)/
!     &     vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,0.d0)
func(T_d)=xkp_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
&     *(T_d**4)*esc_cnt
!     collision with gas particles
&     -xnH*funcL(T_d)
!     stellar radiation
!     &     -4.41d3*xk_vis*funcD(T_d)*xJ_ex
!     CMB
&     -xkp_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
&     *(T_rad**4)*esc_cnt

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
!     Planck mean mass absorption coefficient ("opacity")
!     owing to Dust Grain (solar metalicity)
!     USES vaptemp,linear
DIMENSION Ta(18),xka(18)
DATA (Ta(i),i=1,3)/0.d0,50.d0,100.d0/,
&     (Ta(i),i=6,10)/150.d0,200.d0,250.d0,300.d0,350.d0/
DATA (xka(i),i=1,3)/0.d0,1.d0,2.9d0/,
&     (xka(i),i=6,18)/3.8d0,4.7d0,5.d0,5.25d0,5.3d0,5.3d0,
&     4.3d0,4.3d0,0.9d0,0.9d0,0.4d0,0.4d0,0.d0/

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
!     gives grain volume per unit mass of gas (solar metalicity)
IMPLICIT REAL*8(a-h,o-z)
!     USES vaptemp
x(t,t1)=dmin1(1.d0,ddim(1.025d0,t/t1)/0.05d0)

call vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)

vol_gr=(0.78d-4-0.46d-4*x(T_K,T_tr))*x(T_K,T_ir)
&     +7.19d-4*x(T_K,T_ol)+2.16d-4*x(T_K,T_pyr)
&     +1.18d-4*x(T_K,T_tr)+23.53d-4*x(T_K,T_ro)
&     +6.02d-4*x(T_K,T_vo)+12.93d-4*x(T_K,T_ice)

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
!     USES linear
DIMENSION xlg_rhoa(13),T_icea(13),T_ira(13),
&     T_pyra(13),T_ola(13)
DATA (xlg_rhoa(i),i=1,13)
&     /-24.d0,-22.d0,-20.d0,-18.d0,-16.d0,-14.d0,
&     -12.d0,-10.d0,-8.d0,-6.d0,-4.d0,-2.d0,0.d0/
DATA (T_icea(i),i=1,13)
&     /85.d0,91.d0,98.d0,106.d0,115.d0,125.d0,
&     138.d0,153.d0,172.d0,197.d0,230.d0,271.d0,320.d0/
DATA (T_ira(i),i=1,13)
&     /694.d0,728.d0,775.d0,835.d0,908.d0,994.d0,
&     1100.d0,1230.d0,1395.d0,1612.d0,1908.d0,2283.d0,2737.d0/
DATA (T_pyra(i),i=1,13)
&     /794.d0,827.d0,867.d0,920.d0,980.d0,1049.d0,
&     1129.d0,1222.d0,1331.d0,1462.d0,1621.d0,1808.d0,2023.d0/
DATA (T_ola(i),i=1,13)
&     /791.d0,826.d0,872.d0,929.d0,997.d0,1076.d0,
&     1168.d0,1277.d0,1408.d0,1570.d0,1774.d0,2020.d0,2308.d0/

xlg_rho=dlog10(rho)
call linear(xlg_rhoa,T_icea,13,xlg_rho,T_ice)
T_vo=375.d0
T_ro=575.d0
T_tr=680.d0
call linear(xlg_rhoa,T_ira,13,xlg_rho,T_ir)
call linear(xlg_rhoa,T_pyra,13,xlg_rho,T_pyr)
call linear(xlg_rhoa,T_ola,13,xlg_rho,T_ol)
return

END SUBROUTINE vaptemp


SUBROUTINE phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)
IMPLICIT REAL*8(a-h,o-z)
DATA xm_p/1.67d-24/, T_ro/575.d0/
if(T_gr_K > T_ro) then
   Gmm_pe=0.d0
   return
endif

yHe=8.333d-2
rho=(1.d0+4.d0*yHe)*xm_p*xnH
x=G_0*dsqrt(T_K)/(y_e*xnH)
eps=4.9d-2/(1.d0+4.0d-3*(x**0.73d0))
&     +3.7d-2*(T_K*1.d-4)**0.7d0/(1.d0+2.0d-4*x)

Gam_pe=1.0d-24*eps*G_0*dexp(-1.8d0*A_v)

beta=0.74d0/(T_K**0.068d0)
xLd_pe=4.65d-30*(T_K**0.94d0)*(x**beta)*y_e

Gmm_pe=Z_metal*(Gam_pe*xnH-xLd_pe*xnH**2)/rho

return
END SUBROUTINE phelectr


FUNCTION c_H2(T_K)
IMPLICIT REAL*8(a-h,o-z)
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

return
END SUBROUTINE phelectr


!********************************************************************
!********************************************************************
!**********************  MOLECULAR HYDROGEN  ************************
!********************************************************************
!********************************************************************

function H2_pf(T_K)
implicit real*8(a-h,o-z)
parameter(iv_max=5,j_max=25)
dimension ET(0:iv_max,0:j_max)
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
            z_p=z_p+dble(2*j+1)*dexp(-(ET(iv,j)-ET(0,0))/T_K)
         else
            z_o=z_o+dble(2*j+1)*dexp(-(ET(iv,j)-ET(0,1))/T_K)
         endif
      endif
   enddo
enddo

H2_pf=0.25d0*z_p+0.75d0*z_o

return
end function H2_pf


function E_H2_BFM(iv,J)
!     Borysow, Frommhold, and Moraldi (1989) ApJ,336,495
!     Equation (A11) (in K)
implicit real*8(a-h,o-z)

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
E_H2_BFM=1.43879d0*
&     (-Ev*1.d5+Bv*rrot-Dv*1.d-2*rrot**2+Fv*1.d-5*rrot**3
&     -Gv*1.d-8*rrot**4+Hv*1.d-11*rrot**5-Ov*1.d-14*rrot**6
&     +Pv*1.d-17*rrot**7)

return
end function E_H2_BFM


FUNCTION xLd_H2(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt)
!     xnH*xn(H2)*xLd_H2 : cooling rate owing to H2 per unit volume
implicit real*8(a-h,o-z)
!     USES beta_esc, Q_bg
dimension A(0:2,0:2,0:22,0:22),ET(0:2,0:22),p(0:22),
&     f(0:2,0:22),f_v(0:2)
data xk_B/1.380662d-16/,h_P/6.626176d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/,c_light/2.99792458d10/
data (((A(i,ii,j,j+2),ii=0,i-1),i=1,2),j=0,20)
&     /8.54d-7 ,3.47d-7 ,1.29d-6
&     ,4.23d-7 ,1.61d-7 ,6.40d-7
&     ,2.90d-7 ,1.03d-7 ,4.41d-7
&     ,2.09d-7 ,6.98d-8 ,3.18d-7
&     ,1.50d-7 ,4.72d-8 ,2.28d-7
&     ,1.06d-7 ,3.15d-8 ,1.62d-7
&     ,7.38d-8 ,2.05d-8 ,1.12d-7
&     ,4.98d-8 ,1.31d-8 ,7.50d-8
&     ,3.27d-8 ,8.07d-9 ,4.88d-8
&     ,2.07d-8 ,4.83d-9 ,3.06d-8
&     ,1.27d-8 ,2.79d-9 ,1.85d-8
&     ,7.44d-9 ,1.55d-9 ,1.07d-8
&     ,4.16d-9 ,8.27d-10,5.84d-9
&     ,2.19d-9 ,4.21d-10,3.00d-9
&     ,1.08d-9 ,2.04d-10,1.43d-9
&     ,4.87d-10,9.45d-11,6.17d-10
&     ,1.96d-10,4.23d-11,2.33d-10
&     ,6.71d-11,1.91d-11,7.32d-11
&     ,1.82d-11,9.63d-12,1.73d-11
&     ,3.38d-12,6.26d-12,2.46d-12
&     ,2.96d-13,5.78d-12,1.16d-13/

data (((A(i,ii,j,j),ii=0,i-1),i=1,2),j=1,20)
&     /4.29d-7 ,1.94d-7 ,6.37d-7
&     ,3.03d-7 ,1.38d-7 ,4.50d-7
&     ,2.78d-7 ,1.29d-7 ,4.12d-7
&     ,2.65d-7 ,1.25d-7 ,3.91d-7
&     ,2.55d-7 ,1.23d-7 ,3.74d-7
&     ,2.45d-7 ,1.21d-7 ,3.58d-7
&     ,2.34d-7 ,1.20d-7 ,3.40d-7
&     ,2.23d-7 ,1.18d-7 ,3.22d-7
&     ,2.12d-7 ,1.17d-7 ,3.03d-7
&     ,1.99d-7 ,1.15d-7 ,2.84d-7
&     ,1.87d-7 ,1.13d-7 ,2.63d-7
&     ,1.74d-7 ,1.11d-7 ,2.42d-7
&     ,1.61d-7 ,1.08d-7 ,2.21d-7
&     ,1.47d-7 ,1.05d-7 ,2.01d-7
&     ,1.34d-7 ,1.02d-7 ,1.80d-7
&     ,1.21d-7 ,9.88d-8 ,1.60d-7
&     ,1.08d-7 ,9.50d-8 ,1.41d-7
&     ,9.61d-8 ,9.07d-8 ,1.22d-7
&     ,8.43d-8 ,8.62d-8 ,1.05d-7
&     ,7.32d-8 ,8.14d-8 ,8.85d-8/

data (((A(i,ii,j,j-2),ii=0,i),i=0,2),j=2,20)
&     /2.94d-11,2.53d-7 ,2.79d-11
&     ,1.27d-7 ,3.68d-7 ,2.56d-11
&     ,4.76d-10,3.47d-7,4.50d-10
&     ,1.90d-7 ,4.98d-7 ,4.12d-10
&     ,2.76d-9 ,3.98d-7 ,2.59d-9
&     ,2.38d-7 ,5.60d-7 ,2.37d-9
&     ,9.84d-9 ,4.21d-7 ,9.21d-9
&     ,2.77d-7 ,5.77d-7 ,8.37d-9
&     ,2.64d-8 ,4.19d-7 ,2.46d-8
&     ,3.07d-7 ,5.57d-7 ,2.22d-8
&     ,5.88d-8 ,3.96d-7 ,5.44d-8
&     ,3.28d-7 ,5.05d-7 ,4.88d-8
&     ,1.14d-7 ,3.54d-7 ,1.05d-7
&     ,3.39d-7 ,4.30d-7 ,9.33d-8
&     ,2.00d-7 ,2.98d-7 ,1.82d-7
&     ,3.40d-7 ,3.38d-7 ,1.61d-7
&     ,3.24d-7 ,2.34d-7 ,2.92d-7
&     ,3.30d-7 ,2.41d-7 ,2.55d-7
&     ,4.90d-7 ,1.68d-7 ,4.38d-7
&     ,3.12d-7 ,1.49d-7 ,3.79d-7
&     ,7.03d-7 ,1.05d-7 ,6.21d-7
&     ,2.85d-7 ,7.20d-8 ,5.32d-7
&     ,9.64d-7 ,5.30d-8 ,8.42d-7
&     ,2.53d-7 ,1.96d-8 ,7.13d-7
&     ,1.27d-6 ,1.65d-8 ,1.10d-6
&     ,2.16d-7 ,1.49d-11,9.18d-7
&     ,1.62d-6 ,4.26d-10,1.38d-6
&     ,1.78d-7 ,1.91d-8 ,1.14d-6
&     ,2.00d-6 ,8.38d-9 ,1.69d-6
&     ,1.39d-7 ,8.07d-8 ,1.37d-6
&     ,2.41d-6 ,4.27d-8 ,2.00d-6
&     ,1.01d-7 ,1.86d-7 ,1.61d-6
&     ,2.83d-6 ,1.04d-7 ,2.32d-6
&     ,6.80d-8 ,3.35d-7 ,1.84d-6
&     ,3.26d-6 ,1.93d-7 ,2.64d-6
&     ,3.98d-8 ,5.24d-7 ,2.05d-6
&     ,3.68d-6 ,3.08d-7 ,2.93d-6
&     ,1.84d-8 ,7.49d-7 ,2.23d-6/

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
gamma_10Hp=1.4d-4*T_K**(-1.344d0)
&     *(1.d0+4.005d-9*T_K**2.066d0)
&     *dexp((DT_10-9589d0)/T_K)
gamma_20Hp=4.585d-5*T_K**(-1.291d0)
&     *(1.d0+1.378d-8*T_K**1.903d0)
&     *dexp((DT_20-1.933d4)/T_K)
gamma_21Hp=5.593d-3*T_K**(-1.770d0)
&     *(1.d0+3.505d-10*T_K**2.406d0)
&     *dexp((DT_21-1.460d4)/T_K)

xn_a=y_H*xnH
xn_m=y_H2*xnH
xn_e=y_e*xnH
xn_Hp=y_Hp*xnH

C_10=gamma_10H*xn_a+gamma_10H2*xn_m
&     +gamma_10e*xn_e+gamma_10Hp*xn_Hp
C_20=gamma_20H*xn_a+gamma_20H2*xn_m
&     +gamma_20e*xn_e+gamma_20Hp*xn_Hp
C_21=gamma_21H*xn_a+gamma_21H2*xn_m
&     +gamma_21e*xn_e+gamma_21Hp*xn_Hp

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

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
&     /( (R_01+R_02+R_20)*(R_10+R_12+R_21)
&     -(R_01-R_21)*(R_10-R_20) )
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
         gamma_H=2.93d-14+1.21d-15*T_K
&              +2.16d-19*T_K**2+1.32d-21*T_K**3
      elseif(j == 1) then
         gamma_H=8.34d-14+5.97d-16*T_K
&              +7.76d-19*T_K**2+1.72d-21*T_K**3
      elseif(j == 2) then
         gamma_H=7.54d-14+2.65d-16*T_K
&              +2.14d-19*T_K**2+2.65d-21*T_K**3
      elseif(j == 3) then
         gamma_H=2.95d-14+1.21d-16*T_K
&              +5.53d-19*T_K**2+2.51d-21*T_K**3
      else
!     HM(1989)
         gamma_H=4.6d-12*(2.d0*dble(j+2)-3.d0)*dsqrt(T_K)
&              *dsqrt(1.d0+2.d0*85.25d0*(2.d0*dble(j+2)-1.d0)/T_K)
&              *dexp(-(1.d1*85.25d0*(2.d0*dble(j+2)-1.d0))
&              /(T_K+85.25d0*dble(j+2)*dble(j+3))
&              -0.1187d0*(4.d0*dble(j+2)-2.d0))
      endif

!     H2-H2 collision
      gamma_H2=(3.3d-12+6.6d-15*T_K)
&           *0.276d0*dble((j+2)**2)
&           *dexp(-(dble(j+2)/3.18d0)**1.7d0)

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
   gamma_e=1.d-10*(xj+2.d0)*(xj+1.d0)/(2.d0*xj+5.d0)
&      *(1.d0+1.5d0/x)

      C_ul=gamma_H*xn_a+gamma_H2*xn_m
&        +gamma_Hp*xn_Hp+gamma_e*xn_e
      DT=ET(iv,j+2)-ET(iv,j)
      Q_ul=Q_bg(DT)
      g_u=2.d0*dble(j)+5.d0
      g_l=2.d0*dble(j)+1.d0
      r=(g_u/g_l)*(A(iv,iv,j+2,j)*Q_ul+C_ul*dexp(-DT/T_K))
&           /(A(iv,iv,j+2,j)*(1.d0+Q_ul)+C_ul)
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
         tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3
&              *(xNc_l*g_u/g_l-xNc_u)/v_th
         esc=beta_esc(tau_ul,tau_cnt)
         if(f(ivi,ji).ne.0.d0) then
            xLd_H2=xLd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc
&                 *(1.d0-Q_ul*((g_u/g_l)
&                 *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
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
         tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3
&              *(xNc_l*g_u/g_l-xNc_u)/v_th
         esc=beta_esc(tau_ul,tau_cnt)
         if(f(ivi,ji).ne.0.d0) then
            xLd_H2=xLd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc
&                 *(1.d0-Q_ul*((g_u/g_l)
&                 *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
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
         tau_ul=(A(ivi,ivf,ji,jf)/8.d0/pi)*(c_light/xnu_Hz)**3
&              *(xNc_l*g_u/g_l-xNc_u)/v_th
         esc=beta_esc(tau_ul,tau_cnt)
         if(f(ivi,ji).ne.0.d0) then
            xLd_H2=xLd_H2+f(ivi,ji)*A(ivi,ivf,ji,jf)*DE*esc
&                 *(1.d0-Q_ul*((g_u/g_l)
&                 *(f(ivf,jf)/f(ivi,ji))-1.d0))/xnH
         endif
      enddo
   enddo
enddo
return
END FUNCTION xLd_H2


SUBROUTINE H2Ocool(xnH,T_K,y_H2,xNc_H2O,tau_cnt,xLd_H2O)
!c     H2O cooling function by Neufeld et al.
implicit real*8(a-h,o-z)
data xk_B/1.38d-16/, xm_p/1.67d-24/
!c     100, 200, 400, 1000, 2000, 4000K
real*8 :: xlTa(1:6)
data xlTa/2.000d0, 2.301d0, 2.602d0,
&     3.000d0, 3.301d0, 3.602d0/
real*8 :: xlNa(1:10)
data xlNa/10.0d0, 11.0d0, 12.0d0, 13.0d0, 14.0d0,
&     15.0d0, 16.0d0, 17.0d0, 18.0d0, 19.0d0/
real*8 :: aL0a(1:6)
data aL0a/24.35d0, 23.87d0, 23.42d0,
&     22.88d0, 22.50d0, 22.14d0/
real*8 :: aLLTEa(1:6,1:10)
data aLLTEa/14.59d0, 13.85d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,
&     14.59d0, 13.86d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,
&     14.60d0, 13.86d0, 13.16d0, 12.32d0, 11.86d0, 11.64d0,
&     14.68d0, 13.88d0, 13.17d0, 12.32d0, 11.86d0, 11.64d0,
&     14.98d0, 14.05d0, 13.25d0, 12.34d0, 11.87d0, 11.65d0,
&     15.53d0, 14.46d0, 13.53d0, 12.49d0, 11.97d0, 11.72d0,
&     16.22d0, 15.05d0, 14.02d0, 12.87d0, 12.35d0, 12.06d0,
&     17.00d0, 15.74d0, 14.63d0, 13.46d0, 12.97d0, 12.66d0,
&     17.83d0, 16.50d0, 15.32d0, 14.16d0, 13.69d0, 13.36d0,
&     18.70d0, 17.31d0, 16.07d0, 14.94d0, 14.46d0, 14.13d0/
real*8 :: xlnha(1:6,1:10)
data xlnha/9.00d0, 9.04d0, 9.19d0, 9.50d0, 9.67d0, 9.60d0,
&     8.99d0, 9.04d0, 9.19d0, 9.50d0, 9.67d0, 9.60d0,
&     8.96d0, 9.03d0, 9.19d0, 9.50d0, 9.66d0, 9.59d0,
&     8.74d0, 8.89d0, 9.11d0, 9.47d0, 9.65d0, 9.59d0,
&     8.11d0, 8.37d0, 8.73d0, 9.31d0, 9.56d0, 9.53d0,
&     7.20d0, 7.51d0, 7.95d0, 8.74d0, 9.15d0, 9.20d0,
&     6.22d0, 6.53d0, 6.99d0, 7.87d0, 8.38d0, 8.50d0,
&     5.22d0, 5.57d0, 6.03d0, 6.94d0, 7.48d0, 7.64d0,
&     4.24d0, 4.59d0, 5.09d0, 6.02d0, 6.59d0, 6.78d0,
&     3.21d0, 3.58d0, 4.10d0, 5.08d0, 5.69d0, 5.89d0/
real*8 :: alphaa(1:6,1:10)
data alphaa/0.43d0, 0.42d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0,
&     0.43d0, 0.42d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0,
&     0.42d0, 0.41d0, 0.39d0, 0.36d0, 0.34d0, 0.34d0,
&     0.41d0, 0.39d0, 0.37d0, 0.35d0, 0.33d0, 0.33d0,
&     0.42d0, 0.38d0, 0.34d0, 0.33d0, 0.32d0, 0.32d0,
&     0.45d0, 0.38d0, 0.34d0, 0.32d0, 0.30d0, 0.30d0,
&     0.47d0, 0.40d0, 0.35d0, 0.32d0, 0.29d0, 0.30d0,
&     0.50d0, 0.42d0, 0.36d0, 0.32d0, 0.28d0, 0.29d0,
&     0.52d0, 0.44d0, 0.37d0, 0.31d0, 0.27d0, 0.28d0,
&     0.53d0, 0.45d0, 0.39d0, 0.31d0, 0.27d0, 0.27d0/

!     10, 20, 30, 50, 80, 100K
real*8 :: xlTb(1:6)
data xlTb/1.000d0, 1.301d0, 1.477d0,
&     1.699d0, 1.903d0, 2.000d0/
real*8 :: aL0bo(1:6)
data aL0bo/26.81d0, 25.88d0, 25.43d0,
&     24.96d0, 24.58d0, 24.41d0/
real*8 :: aLLTEbo(1:6,1:10)
data aLLTEbo/17.94d0, 16.71d0, 16.08d0, 15.41d0, 14.85d0, 14.60d0,
&     17.96d0, 16.72d0, 16.09d0, 15.42d0, 14.86d0, 14.60d0,
&     18.14d0, 16.86d0, 16.19d0, 15.47d0, 14.88d0, 14.62d0,
&     18.77d0, 17.36d0, 16.58d0, 15.72d0, 15.02d0, 14.73d0,
&     19.70d0, 18.11d0, 17.25d0, 16.27d0, 15.47d0, 15.11d0,
&     20.67d0, 18.93d0, 18.05d0, 17.01d0, 16.12d0, 15.72d0,
&     21.58d0, 19.81d0, 18.91d0, 17.82d0, 16.88d0, 16.45d0,
&     22.53d0, 20.71d0, 19.80d0, 18.69d0, 17.70d0, 17.25d0,
&     23.50d0, 21.64d0, 20.72d0, 19.58d0, 18.56d0, 18.10d0,
&     24.41d0, 22.58d0, 21.65d0, 20.49d0, 19.46d0, 18.99d0/
real*8 :: xlnhbo(1:6,1:10)
data  xlnhbo/8.81d0, 8.90d0, 9.03d0, 9.14d0, 9.19d0, 9.20d0,
&      8.79d0, 8.88d0, 9.01d0, 9.13d0, 9.20d0, 9.20d0,
&      8.60d0, 8.73d0, 8.88d0, 9.03d0, 9.13d0, 9.15d0,
&      7.96d0, 8.14d0, 8.31d0, 8.55d0, 8.78d0, 8.85d0,
&      7.01d0, 7.21d0, 7.40d0, 7.69d0, 8.00d0, 8.11d0,
&      6.02d0, 6.20d0, 6.41d0, 6.71d0, 7.04d0, 7.17d0,
&      5.03d0, 5.21d0, 5.42d0, 5.70d0, 6.06d0, 6.19d0,
&      4.02d0, 4.21d0, 4.41d0, 4.71d0, 5.05d0, 5.19d0,
&      3.02d0, 3.22d0, 3.41d0, 3.71d0, 4.06d0, 4.20d0,
&      2.03d0, 2.21d0, 2.42d0, 2.70d0, 3.06d0, 3.20/
real*8 :: alphabo(1:6,1:10)
data alphabo/0.71d0, 0.49d0, 0.48d0, 0.46d0, 0.45d0, 0.45d0,
&     0.64d0, 0.49d0, 0.48d0, 0.45d0, 0.45d0, 0.45d0,
&     0.57d0, 0.50d0, 0.47d0, 0.45d0, 0.44d0, 0.44d0,
&     0.56d0, 0.53d0, 0.48d0, 0.44d0, 0.42d0, 0.41d0,
&     0.59d0, 0.60d0, 0.53d0, 0.47d0, 0.45d0, 0.43d0,
&     0.72d0, 0.64d0, 0.58d0, 0.52d0, 0.49d0, 0.47d0,
&     0.85d0, 0.68d0, 0.61d0, 0.56d0, 0.52d0, 0.50d0,
&     0.86d0, 0.70d0, 0.63d0, 0.58d0, 0.55d0, 0.53d0,
&     0.93d0, 0.72d0, 0.66d0, 0.61d0, 0.57d0, 0.55d0,
&     0.87d0, 0.73d0, 0.67d0, 0.62d0, 0.59d0, 0.56/
c
real*8 :: aL0bp(1:6)
data aL0bp/27.01d0, 25.73d0, 25.24d0,
&     24.75d0, 24.38d0, 24.22d0/
real*8 :: aLLTEbp(1:6,1:10)
data aLLTEbp/17.72d0, 16.60d0, 16.12d0, 15.43d0, 14.86d0, 14.60d0,
&     17.76d0, 16.63d0, 16.13d0, 15.43d0, 14.86d0, 14.60d0,
&     18.07d0, 16.85d0, 16.24d0, 15.48d0, 14.88d0, 14.62d0,
&     18.83d0, 17.41d0, 16.61d0, 15.72d0, 15.03d0, 14.73d0,
&     19.68d0, 18.11d0, 17.26d0, 16.28d0, 15.47d0, 15.11d0,
&     20.50d0, 18.94d0, 18.06d0, 17.01d0, 16.12d0, 15.72d0,
&     21.37d0, 19.83d0, 18.93d0, 17.82d0, 16.87d0, 16.45d0,
&     22.28d0, 20.75d0, 19.81d0, 18.69d0, 17.70d0, 17.25d0,
&     23.22d0, 21.69d0, 20.73d0, 19.58d0, 18.56d0, 18.10d0,
&     24.19d0, 22.60d0, 21.65d0, 20.49d0, 19.45d0, 18.99d0/
real*8 :: xlnhbp(1:6,1:10)
data xlnhbp/9.30d0, 9.11d0, 8.94d0, 8.70d0, 8.55d0, 8.50d0,
&     9.25d0, 9.06d0, 8.91d0, 8.69d0, 8.54d0, 8.49d0,
&     8.94d0, 8.82d0, 8.71d0, 8.56d0, 8.46d0, 8.43d0,
&     8.17d0, 8.16d0, 8.16d0, 8.12d0, 8.10d0, 8.11d0,
&     7.21d0, 7.24d0, 7.28d0, 7.30d0, 7.32d0, 7.35d0,
&     6.20d0, 6.24d0, 6.31d0, 6.33d0, 6.35d0, 6.38d0,
&     5.21d0, 5.25d0, 5.30d0, 5.32d0, 5.35d0, 5.39d0,
&     4.22d0, 4.26d0, 4.31d0, 4.33d0, 4.36d0, 4.40d0,
&     3.23d0, 3.24d0, 3.31d0, 3.33d0, 3.37d0, 3.41d0,
&     2.21d0, 2.25d0, 2.30d0, 2.32d0, 2.37d0, 2.41d0/
real*8 :: alphabp(1:6,1:10)
data alphabp/0.49d0, 0.72d0, 0.69d0, 0.53d0, 0.46d0, 0.44d0,
&     0.52d0, 0.65d0, 0.66d0, 0.53d0, 0.46d0, 0.44d0,
&     0.49d0, 0.65d0, 0.63d0, 0.51d0, 0.45d0, 0.43d0,
&     0.57d0, 0.73d0, 0.65d0, 0.49d0, 0.43d0, 0.41d0,
&     0.75d0, 0.68d0, 0.64d0, 0.52d0, 0.45d0, 0.42d0,
&     0.76d0, 0.70d0, 0.68d0, 0.55d0, 0.47d0, 0.43d0,
&     0.79d0, 0.73d0, 0.69d0, 0.58d0, 0.48d0, 0.44d0,
&     0.80d0, 0.75d0, 0.71d0, 0.58d0, 0.50d0, 0.46d0,
&     0.82d0, 0.77d0, 0.73d0, 0.60d0, 0.52d0, 0.48d0,
&     0.91d0, 0.80d0, 0.74d0, 0.62d0, 0.53d0, 0.50d0/

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
   xLinv=xL0inv+xn_c*xLLTEinv
&        +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
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
   xLinv=xL0inv+xn_c*xLLTEinv
&        +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
   xLo=1.d0/xLinv

   xlNp=dlog10(0.25d0*xNc_H2O/cs)
   call linear(xlTb,aL0bp,6,xlT,aL0p)
   call bilinear(xlTb,xlNa,aLLTEbp,6,10,xlT,xlNp,aLLTE)
   call bilinear(xlTb,xlNa,xlnhbp,6,10,xlT,xlNp,xlnh)
   call bilinear(xlTb,xlNa,alphabp,6,10,xlT,xlNp,alpha)
   xL0inv=10.d0**aL0
   xLLTEinv=10.d0**aLLTE
   xn_h=10.d0**xlnh
   xLinv=xL0inv+xn_c*xLLTEinv
&        +xL0inv*((xn_c/xn_h)**alpha)*(1.d0-xn_h*xLLTEinv/xL0inv)
   xLp=1.d0/xLinv

   xL=0.75d0*xLo+0.25d0*xLp
   xL=1.148d0*xL
endif

xLd_H2O=(1.d0-y_H2)*xL*dexp(-tau_cnt)
END SUBROUTINE H2Ocool


FUNCTION xLd_HD(xnH,T_K,y_H,y_H2,xNc_HD,tau_cnt)
! HD cooling from Galli and Palla (1998), A&A, 335, 403
! cooling is for T < 3000 K, in erg cm^3 s-1
implicit real*8 (a-h,o-z)
data xk_B/1.38066d-16/,h_Pl/6.62618d-27/,
&     pi/3.14159265358979d0/,xm_p/1.67d-24/
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

an0d = -(w10d*w21d*w32d+w30d*w2d*w1d+w20d*w31d*w12d-
&         w12d*w21d*w30d+w2d*w31d*w10d+w32d*w20d*w1d)
an1d = -w0d*w21d*w32d+w20d*w31d*w02d-w30d*w2d*w01d-
&         w02d*w21d*w30d-w2d*w31d*w0d-w32d*w20d*w01d
an2d = -(w0d*w1d*w32d+w10d*w31d*w02d+w30d*w12d*w01d+
&         w02d*w1d*w30d+w12d*w31d*w0d-w32d*w10d*w01d)
an3d = -w0d*w1d*w2d+w10d*w21d*w02d+w20d*w12d*w01d+
&         w02d*w1d*w20d+w12d*w21d*w0d+w2d*w10d*w01d
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

tau10d = (a10d/8.d0/pi)*(3.d10/xnu10d)**3
&         *(xNc0d*g1/g0-xNc1d)/v_D
tau21d = (a21d/8.d0/pi)*(3.d10/xnu21d)**3
&         *(xNc1d*g2/g1-xNc2d)/v_D
tau32d = (a32d/8.d0/pi)*(3.d10/xnu32d)**3
&         *(xNc2d*g3/g2-xNc3d)/v_D
tau20d = (a20d/8.d0/pi)*(3.d10/xnu20d)**3
&         *(xNc0d*g2/g0-xNc2d)/v_D
tau31d = (a31d/8.d0/pi)*(3.d10/xnu31d)**3
&         *(xNc1d*g3/g1-xNc3d)/v_D

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
PARAMETER (N_p=100)
DIMENSION esc(N_p),func(N_p),desc(N_p),
&     esc_f(N_p),func_f(N_p),A(N_p,N_p),
&     esc_min(N_p)

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
IMPLICIT REAL*8(a-h,o-z)
!     USES beta_esc
PARAMETER (N_p=100)
DIMENSION f_LTE(0:N_p),A(N_p),f(0:N_p),aa(N_p),
&     DDD(0:N_p),xnn(0:N_p),esc(N_p),func(N_p),Q(N_p),S(N_p)
COMMON /mol_lines/ A0,DT0,sigma,eta_T,xmu_mol,J_max
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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
   aa(j)=A(j)*esc(j)
&        /(xnH*sigma*v_T)/(1.d0-y_m)
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
      BBB=aa(j)*(1.d0+Q(j))
&           +(dble(2*j+3)/dble(2*j+1))*aa(j+1)*Q(j+1)+1.d0
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
   tau=(A(j)/8.d0/pi)*(3.d10/xnu)**3
&        *(xNc_l*g_u/g_l-xNc_u)/v_th/eta_T
   func(j)=esc(j)-beta_esc(tau,tau_cnt)
   xLd_rot=xLd_rot
&    +f(j)*A(j)*DE*beta_esc(tau,tau_cnt)*(1.d0-Q(j)/S(j))/xnH
enddo

return
END pop_rot


SUBROUTINE pop_CII(esc,func,xLd_CII)
IMPLICIT REAL*8(a-h,o-z)
PARAMETER(N_p=100)
DIMENSION esc(N_p),func(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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
IMPLICIT REAL*8(a-h,o-z)
PARAMETER(N_p=100)
DIMENSION esc(N_p),func(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
&     /( (R_01+R_02+R_20)*(R_10+R_12+R_21)
&     -(R_01-R_21)*(R_10-R_20) )
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

xLd_CI=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20
&     +DE_21*A_21*f_2*x_21)/xnH
return

END SUBROUTINE pop_CI

SUBROUTINE pop_OI(esc,func,xLd_OI)
!     xnH*rn(OI)*xLd_OI : cooling rate owing to OI per unit volume
IMPLICIT REAL*8(a-h,o-z)
PARAMETER(N_p=100)
DIMENSION esc(N_p),func(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

g_0=5.d0
g_1=3.d0
g_2=1.d0

esc_10=1.0d0!esc(1)
esc_20=1.0d0!esc(2)
esc_21=1.0d0!esc(3)

DT_10=2.3d2
DT_20=3.28d2
DT_21=9.8d1

DE_10=DT_10*1.38d-16
DE_20=DT_20*1.38d-16
DE_21=DT_21*1.38d-16

!Q_10=Q_bg(DT_10)
!Q_20=Q_bg(DT_20)
!Q_21=Q_bg(DT_21)
Q_10=0.d0
Q_20=0.d0
Q_21=0.d0

A_10=9.0d-5
A_20=1.0d-10
A_21=1.7d-5

gamma10_e=1.4d-8
gamma20_e=1.4d-8
gamma21_e=5.0d-9

gamma10_H=9.2d-11*(T_K*1.d-2)**0.67d0
gamma20_H=4.3d-11*(T_K*1.d-2)**0.80d0
gamma21_H=1.1d-10*(T_K*1.d-2)**0.44d0

gamma10_H2=5.d-1*gamma10_H
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

R_10=esc_10*A_10*(1.d0+Q_10)+C_10
R_20=esc_20*A_20*(1.d0+Q_20)+C_20
R_21=esc_21*A_21*(1.d0+Q_21)+C_21
R_01=(g_1/g_0)*esc_10*A_10*Q_10+C_01
R_02=(g_2/g_0)*esc_20*A_20*Q_20+C_02
R_12=(g_2/g_1)*esc_21*A_21*Q_21+C_12

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
&     /( (R_01+R_02+R_20)*(R_10+R_12+R_21)
&     -(R_01-R_21)*(R_10-R_20) )
f_1=(f_0*(R_01-R_21)+R_21)/(R_10+R_12+R_21)
f_2=(f_0*R_02+f_1*R_12)/(R_21+R_20)

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

xLd_OI=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20
&     +DE_21*A_21*f_2*x_21)/xnH

return
END SUBROUTINE pop_OI


SUBROUTINE pop_CIImeta(esc,func,xLd_CIImeta)
IMPLICIT REAL*8(a-h,o-z)
PARAMETER(N_p=100)
DIMENSION esc(N_p),func(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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

return

END SUBROUTINE pop_CIImeta

SUBROUTINE pop_CImeta(esc,func,xLd_CImeta)
IMPLICIT REAL*8(a-h,o-z)
PARAMETER(N_p=100)
DIMENSION esc(N_p),func(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
&     /( (R_01+R_02+R_20)*(R_10+R_12+R_21)
&     -(R_01-R_21)*(R_10-R_20) )
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

xLd_CImeta=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20
&     +DE_21*A_21*f_2*x_21)/xnH

return
END SUBROUTINE pop_CImeta

SUBROUTINE pop_OImeta(esc,func,xLd_OImeta)
!     xnH*rn(OI)*xLd_OI : cooling rate owing to OI per unit volume
IMPLICIT REAL*8(a-h,o-z)
PARAMETER(N_p=100)
DIMENSION esc(N_p),func(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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

f_0=( R_21*(R_10-R_20)+R_20*(R_10+R_12+R_21) )
&     /( (R_01+R_02+R_20)*(R_10+R_12+R_21)
&     -(R_01-R_21)*(R_10-R_20) )
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

xLd_OImeta=(DE_10*A_10*f_1*x_10+DE_20*A_20*f_2*x_20
&     +DE_21*A_21*f_2*x_21)/xnH
return
END SUBROUTINE pop_OImeta

SUBROUTINE thick_lev(J_thick)
IMPLICIT REAL*8(a-h,o-z)
PARAMETER (N_p=100)
DIMENSION f_LTE(0:N_p),A(N_p),f(0:N_p),aa(N_p),P_J(N_p),esc(N_p)
COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc
COMMON /mol_lines/ A0,DT0,sigma,eta_T,xmu_mol,J_max
DATA xk_B/1.38d-16/,h_Pl/6.63d-27/,pi/3.14159265358979d0/,
&     xm_p/1.67d-24/

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
   tau_line=(A(j)/8.d0/pi)*(3.d10/xnu)**3
&        *(xNc_l*g_u/g_l-xNc_u)/v_th/eta_T
   if(tau_line.gt.1.d-2) then
      J_thick=j
   endif
enddo

return
END SUBROUTINE thick_lev



FUNCTION beta_esc(tau_L,tau_C)
IMPLICIT REAL*8(a-h,o-z)

if(tau_L.lt.0.d0) then
   beta_esc=1.d0
elseif(tau_L.lt.1.d-5) then
   beta_esc=dexp(-tau_C)
else
   beta_esc=dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
endif

return
END FUNCTION beta_esc


