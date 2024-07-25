!-------------------
!------MODURES------
!-------------------
MODULE UVCR
double precision :: zeta, G_0
END MODULE UVCR
   
MODULE radtemp
double precision :: T_rad
END MODULE radtemp
   
MODULE coldens
double precision :: A_v,xNc_H,xNc_H2,xNc_HD
END MODULE coldens
   
MODULE lines
double precision :: xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc  
END MODULE lines
   
MODULE indcool
double precision ::  xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,&
Lmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
END MODULE indcool
   
MODULE radgr
double precision :: T_gr_K,tau_cont
END MODULE radgr
   
MODULE mol_lines 
double precision :: A0,DT0,sigma,eta_T,xmu_mol
integer :: J_max=25
END MODULE mol_lines 
   
MODULE xcrit
double precision ::xH,xH2,xHe
END MODULE xcrit

!-------------------
!-----INCLUDES------
!-------------------
include "subs3/chemistry.f90"
include "subs3/reaction.f90"  !full reactions
include "subs3/grain.f90"
include "subs3/xk_prm.f90"
include "subs3/H2.f90"
include "subs3/HD.f90"
include "subs3/lines.f90"
include "subs3/CR.f90"
include "subs3/math.f90"
include "subs3/CO.f90"
include "subs3/H2O.f90"
include "subs3/OH.f90"

!-------------------
!---MAIN PROGRAM----
!-------------------
PROGRAM ZR
!     computes the evolutionary path of 
!     metal polluted protostellar clouds 
!     with radiation
!     updated 2008 Jun
!     f90 ver. (2024.07.25), MAEDA
USE UVCR
USE radtemp
USE coldens
USE lines 
USE indcool
USE radgr
IMPLICIT REAL*8(a-h,o-z)
integer, PARAMETER :: N_sp=50 
double precision :: y(N_sp), t_ratio
double precision :: xm_p=1.67d-24,xk_B=1.38d-16,G=6.67d-8,&
                 pi=3.14159265358979d0

!---parameters---
Z_gas=1.0d-6
Z_dust=1.0d-6
!G_0=1.d0
!zeta=1.d-17
G_0=0.d0
zeta=0.d0
T_rad=10.d0
yHe=8.333d-2
yD=3.d-5
!t_ratio=4.d-3 !for Z=1.d0
t_ratio=2.d-2

!initial condition 
xnH=1.d-1
!xnH=1.d5
rho=(1.d0+4.d0*yHe)*xm_p*xnH    
T_K=300.d0
!T_K=10.d0
!T_K=8000.d0
Z_metal_gas=Z_gas

!evolution
i_ev=1
!---parameters---
            
open(8,file='data2/nT_z6_nm1_GZ0t.dat',status='unknown')
open(11,file='data2/1_z6_nm1_GZ0t.dat',status='unknown')
open(12,file='data2/2_z6_nm1_GZ0t.dat',status='unknown')
!open(13,file='3.dat',status='unknown')    
open(14,file='data2/4_z6_nm1_GZ0t.dat',status='unknown')
open(15,file='data2/data_y_z6_nm1_GZ0t.dat',status='unknown')
!open(16,file='data/ar.dat',status='unknown')

!************************************************************
do i=1,N_sp
      y(i)=0.d0
enddo

!---set abundance of non-zero species
!y_Hp=3.d-4
!y_H2=2.d-5
!y_Hp=0.5d0
!y_H2=0.d0
y_Hp=1.d-4
y_H2=1.d-6
y_Dp=0.d0
y_HD=0.d0
!y_HD=1.d-9

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
gamma=1.d0+(1.d0+4.d0*yHe)&
/(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))&
+c_H2(T_K)*y(2)))
e=xk_B*T_K/((gamma-1.d0)*xmu*xm_p)
!************************************************************


T_gr_K=1.d0
xLmbd_ch=0.d0
t=0.d0
dt=1.d-1
itchem=0
t_chem=1.d-1
esc_cnt=1.d0

!---temporal evolution

do i=1,1000000    
      !collapse timescale, radius, bulk velocity
      !1) self-gravitating
      xlmbd_J=dsqrt(pi*xk_B*T_K/(G*xmu*xm_p*rho))
      radius=xlmbd_J
      xNcolJ=xnH*xlmbd_J
      !2) constant column density
      !Ncol=1.d19
      !radius=xNcol/xnH

      t_col=dsqrt(3.d0*pi/(32.d0*G*rho))
      v_bulk=radius/3.d0/t_col
      v_turb=0.d0
      v_bulk=dsqrt(v_bulk**2+v_turb**2)
      !continuum cooling
      A_v=Z_metal*xnH*radius*5.3d-22

      !write(*,*) '-tst1-', i
      call  rad_cool(Z_metal,T_K,T_gr_K,radius,A_v,esc_cnt,&
      dt,xmu,gamma,y,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,i,i_ev)
      tau_cont=tau_cnt

      !write(*,*) '-tst2-', i
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
      !write(*,*) '-tst3-', i
      call chemcool(xnH,T_K,T_gr_K,Z_metal,&
           y,dt,t_chem,xmu,gamma,xLmbd_ch)

      y_H=y(1)
      y_H2=y(2)
      y_e=y(3)
      y_HD=y(13)

      xLmbd_net=xLmbd_cnt+xLmbd_line+xLmbd_gr+xLmbd_ch&
            -Gmm_pe-Gmm_CR         

      xMJ=rho*radius**3/2.d33
      
      !write(*,*) '-tst4-', i

      !update values of rho & e
      call update(i_ev,rho,e,p,xnH,t_col,xLmbd_net,xmu,gamma,dt,yHe,T_K,t)

      !---data output 
      if((mod(i,10).eq.0).and.(i_ev==0))  then 
            write(*,101) t/t_col, T_K, y_H2, y_e
101         format(4e10.3)
      endif
      if((mod(i,10).eq.0).and.(i_ev==1))  then             
            write(*,102) xnH, T_K, T_gr_K, y_H2, y_e
102         format(5e10.3)
            write(8,108) xnH, T_K
            write(11,201) xnH,T_K,T_gr_K,p
            !write(12,202) xnH,y
            write(12,202) xnH,y(2),y(23),y(17),y(33)             
            !write(13,203) xnH,Gmm_cmp,xLmbd_cnt,xLmbd_line,xLmbd_gr,&
            !           xLmbd_ch,Gmm_pe,Gmm_CR
            write(*,*)'xLmbd_CI_1',xLmbd_CI
            !xLmbd_CI=1.d-20!dsign(dmax1(dabs(xLmbd_CI),1.d-10),xLmbd_CI)
            !write(*,*)'xLmbd_CI_2',xLmbd_CI
            write(14,204) xnH, xLmbd_CII+xLmbd_CI,xLmbd_OI,&
                        xLmbd_CII,xLmbd_CI
            !write(14,204) xnH, xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,&
            !xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
            !write(15,108) xMJ,T_K
            !write(15,108) xNH,y(3),y(1),y(2),y(4)&
            !,y(13),y(33),y(32),y(34),y(23),y(17),y(30)
            !write(16,*) xnH,T_K,T_gr_K,A_v,xNc_H2,y
      endif
      if((mod(i,75).eq.0).and.(i_ev==1))  then
            write(15,*) xNH,y(3),y(1),y(2),y(4)&
            ,y(13),y(33),y(32),y(34),y(23),y(17),y(30)
      endif
108   format(2e10.3)
201   format(4e10.3)
202   format(51e10.3)
203   format(8e11.3)
204   format(10e11.3)
 
      !end of calculation
      !if(xnH.gt.1.d23) then
      if(xnH.gt.1.d14) then
            close(11)
            close(12)
            close(13)
            close(14)
            close(15)
            !close(16)
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
            if(dmin1(t_ratio*t_col,t_ratio*t_cool).le.dt) itchem=1
      elseif(xnH < 1.d0) then
            dt=dmin1(t_ratio*t_col,t_ratio*t_cool)
      elseif(xnH < 1.d14) then
            dt=dmin1(t_ratio*t_col,t_ratio*t_cool)
      else
            dt=dmin1(5.d-2*t_col,5.d-2*t_cool)
      endif
      !print *, itchem, 'dt=', t_chem,&
      !dmin1(2.d-2*t_col,2.d-2*t_cool)
      !pause

enddo
END PROGRAM ZR


!-------------------
!----SUBROUTINES----
!-------------------

SUBROUTINE update(i_ev,rho,e,p,xnH,t_col,xLmbd_net,xmu,gamma,dt,yHe,T_K,t)
implicit real*8(a-h,o-z)
double precision :: xm_p=1.67d-24,xk_B=1.38d-16,G=6.67d-8,&
                 pi=3.14159265358979d0

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
      !pause
endif

return
END SUBROUTINE update

SUBROUTINE coltime(gamma_eff,ft)
implicit real*8(a-h,o-z)      
if(gamma_eff < 0.83d0) then
   f=0.d0
elseif(gamma_eff < 1.d0) then
   f=0.6d0+2.5d0*(gamma_eff-1.d0)-6.d0*(gamma_eff-1.d0)**2
elseif(gamma_eff < 1.33333d0) then
   f=1.d0+0.2d0*(gamma_eff-4.d0/3.d0)&
           -2.9d0*(gamma_eff-4.d0/3.d0)**2
else
   f=1.d0
endif

if(f > 0.95d0) then
   ft=1.d0/dsqrt(0.05d0)
else
   ft=1.d0/dsqrt(1.d0-f)
endif

return
END SUBROUTINE coltime


FUNCTION Q_bg(T_nu)
   USE radgr
   USE radtemp
IMPLICIT REAL*8(a-h,o-z)
double precision, PARAMETER :: B_ex=5.d-4
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
END FUNCTION Q_bg


SUBROUTINE rad_cool(Z_metal,Tp,Tp_gr,&
        radius,A_v,esc_cnt,dt,xmu,gamma,y,&
             xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,i,i_ev)
             USE lines
             USE radtemp
IMPLICIT REAL*8(a-h,o-z)
integer, PARAMETER :: maxit=100,N_sp=50
double precision :: y(N_sp)
double precision :: xm_p=1.67d-24,xk_B=1.38d-16,G=6.67d-8,&
     pi=3.14159265358979d0

     func(t,xLmbd)=Tp-t-(gamma-1.d0)*xmu*xm_p*xLmbd*dt/xk_B

Tp0=Tp
Tp1=Tp
if(tau_cnt < 10.d0) then
   T_K=Tp1
   call line_cool(y,radius,xLmbd_line,i)
else
   xLmbd_line=0.d0
endif
call cnt(Z_metal,xnH,Tp1,Tp_gr,&
     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
if(dabs(Tp1/T_rad-1.d0) > 1.d0) return   
xLmbd=xLmbd_ch+xLmbd_line+xLmbd_cnt+xLmbd_gr
fL=func(Tp1,xLmbd)

Tp2=1.01d0*Tp1
call cnt(Z_metal,xnH,Tp2,Tp_gr,&
     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
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
   call cnt(Z_metal,xnH,Tpsec,Tp_gr,&
        radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt1,xLmbd_gr1)
   xLmbd1=xLmbd_ch+xLmbd_line+xLmbd_cnt1+xLmbd_gr1
   f=func(Tpsec,xLmbd1)
   err=dabs(f/Tpsec)
   if(err < 1.d-6) go to 12
enddo
12   continue 
xLmbd_cnt=xLmbd_cnt1
xLmbd_gr=xLmbd_gr1
T_K=Tp0

return
END SUBROUTINE rad_cool
!
!
SUBROUTINE cnt(Z_metal,xnH,T_K,T_gr_K,&
        radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
        USE radtemp
IMPLICIT REAL*8(a-h,o-z)
double precision :: sigma_B=5.67d-5,xm_p=1.67d-24 
!PARAMETER (xk_vis=40.d0,B_ex=5.d-4)
!PARAMETER (xk_vis=40.d0,B_ex=0.d0)

yHe=8.333d-2
!xJ_ex=B_ex*dexp(-A_v/2.5d0)
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
return
ENDSUBROUTINE cnt
!      
!     
SUBROUTINE line_cool(y,radius,xLmbd_line,i)
   USE radtemp
   USE lines
   USE mol_lines
   USE indcool 
IMPLICIT REAL*8(a-h,o-z)
integer, PARAMETER :: N_sp=50,N_p=100
double precision, parameter :: yHe=8.333d-2
DIMENSION y(N_sp),func(N_p),&
     esc_CI(N_p),esc_CII(N_p),esc_OI(N_p),&
     esc_CImeta(N_p),esc_CIImeta(N_p),esc_OImeta(N_p)
EXTERNAL pop_CI,pop_OI,pop_CII,&
     pop_CImeta,pop_OImeta,pop_CIImeta
double precision :: xm_p=1.67d-24,xk_B=1.38d-16,G=6.67d-8,&
pi=3.14159265358979d0
SAVE esc_CI,esc_CII,esc_OI,&
     esc_CImeta,esc_CIImeta,esc_OImeta

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
   
!H2 cooling
v_D=v_th/dsqrt(2.d0)
xNc_H2=dmin1(1.d0,v_D/v_bulk)*xNc_H2
xLmbd_H2=y_m&
     *xLd_H2(xnH,T_K,y_a,y_m,y_e,y_Hp,xNc_H2,tau_cnt)&
          *xnH/((1.d0+4.d0*yHe)*xm_p)

!HD cooling
v_D=v_th/dsqrt(3.d0)
xNc_HD=dmin1(1.d0,v_D/v_bulk)*xNc_HD      
xLmbd_HD=y_HD*xLd_HD(xnH,T_K,y_a,y_m,xNc_HD,tau_cnt)&
     *xnH/((1.d0+4.d0*yHe)*xm_p)

!CO cooling
xmu_mol=28.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CO
call COcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_CO)
call COcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_CO)
xLmbd_CO=y_CO&
     *(xLd_CO-xLdr_CO)&
          *xnH/((1.d0+4.d0*yHe)*xm_p)

!OH cooling
xmu_mol=17.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OH
call OHcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_OH)
call OHcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_OH)
xLmbd_OH=y_OH&
     *(xLd_OH-xLdr_OH)&
          *xnH/((1.d0+4.d0*yHe)*xm_p)

!H2O cooling
xmu_mol=18.d0
v_D=v_th/dsqrt(xmu_mol)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_H2O
call H2Ocool(xnH,T_K,y_m,xNc,tau_cnt,xLd_H2O)
call H2Ocool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_H2O)
xLmbd_H2O=y_H2O&
     *(xLd_H2O-xLdr_H2O)&
          *xnH/((1.d0+4.d0*yHe)*xm_p)

!CII cooling
xmu_C=12.d0
xmu_O=16.d0
v_D=v_th/dsqrt(xmu_C)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CII
!fine structure line
xLmbd_CII=y_CII&
     *xLd(pop_CII,esc_CII,1)&
          *xnH/((1.d0+4.d0*yHe)*xm_p)
!metastable line
if(T_K > 2.d2) then
   xLmbd_CIImeta=y_CII&
           *xLd(pop_CIImeta,esc_CIImeta,1)&
                   *xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_CII=xLmbd_CII+xLmbd_CIImeta
endif

!CI cooling
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CI
!fine structure line
xLmbd_CI=y_CI&
     *xLd(pop_CI,esc_CI,3)&
          *xnH/((1.d0+4.d0*yHe)*xm_p)
!metastable line
if(T_K > 3.d3) then
   xLmbd_CImeta=y_CI&
           *xLd(pop_CImeta,esc_CImeta,3)&
                   *xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_CI=xLmbd_CI+xLmbd_CImeta
endif
  
!OI cooling
v_D=v_th/dsqrt(xmu_O)
xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OI
!fine structure line
if(T_K > 15.d0) then
   xLmbd_OI=y_OI&
           *xLd(pop_OI,esc_OI,3)&
                   *xnH/((1.d0+4.d0*yHe)*xm_p)
else
   xLmbd_OI=0.d0
endif
!metastable line
if(T_K > 3.5d3) then
   xLmbd_OImeta=y_OI&
           *xLd(pop_OImeta,esc_OImeta,3)&
                   *xnH/((1.d0+4.d0*yHe)*xm_p)
   xLmbd_OI=xLmbd_OI+xLmbd_OImeta
endif

!HI Lya cooling 
if(xnH < 1.d5) then
   T_5=1.d-5*T_K
   xLmbd_Lya=y_e*y_a*7.50d-19/(1.d0+dsqrt(T_5))&
           *dexp(-1.18348d5/T_K)&
                   *xnH/((1.d0+4.d0*yHe)*xm_p)
else
   xLmbd_Lya=0.d0
endif

!total cooling rate by lines
xLmbd_line=xLmbd_H2+xLmbd_HD+xLmbd_CO+xLmbd_OH+xLmbd_H2O&
     +xLmbd_CII+xLmbd_CI+xLmbd_OI+xLmbd_Lya     
    
return
END  SUBROUTINE line_cool





















