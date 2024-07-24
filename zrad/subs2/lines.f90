FUNCTION xLd(pop,esc,N_line)
IMPLICIT REAL*8(a-h,o-z)
!double precision :: xLd
!integer :: N_line
integer, PARAMETER :: N_p=100
double precision :: esc(N_p),func(N_p),desc(N_p),&
esc_f(N_p),func_f(N_p),A(N_p,N_p),&
esc_min(N_p)

err_max_min=1.d10
do itr=1,1000
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
enddo
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
   !use COMVAR
   USE lines
   USE mol_lines
   IMPLICIT REAL*8(a-h,o-z)
   !double precision :: xLd_rot
   !     USES beta_esc
   integer, PARAMETER :: N_p=100
   double precision :: f_LTE(0:N_p),A(N_p),f(0:N_p),aa(N_p),&
   DDD(0:N_p),xnn(0:N_p),esc(N_p),func(N_p),Q(N_p),S(N_p)

   DOUBLE PRECISION :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
        xm_p=1.67d-24


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
   !use COMVAR
   USE lines
   IMPLICIT REAL*8(a-h,o-z)
   !double precision :: xLd_CII
   integer, PARAMETER :: N_p=100
   double precision :: esc(N_p),func(N_p)
   DOUBLE PRECISION :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
        xm_p=1.67d-24

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
   !use COMVAR
   USE lines  
   IMPLICIT REAL*8(a-h,o-z)
   !double precision :: xLd_CI
   integer, PARAMETER :: N_p=100
   double precision :: esc(N_p),func(N_p)
   DOUBLE PRECISION :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
        xm_p=1.67d-24


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
   !use COMVAR
   USE lines
   IMPLICIT REAL*8(a-h,o-z)
   !double precision :: xLd_OI
   integer, PARAMETER :: N_p=100
   double precision :: esc(N_p),func(N_p)
   double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
    xm_p=1.67d-24


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
   !use COMVAR
   USE lines
   IMPLICIT REAL*8(a-h,o-z)
   !double precision :: xLd_CIImeta
   integer, PARAMETER :: N_p=100
   double precision :: esc(N_p),func(N_p)
   double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
   xm_p=1.67d-24


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
    !use COMVAR
    USE lines
    IMPLICIT REAL*8(a-h,o-z)
    !double precision :: Ld_CImeta
    integer, PARAMETER :: N_p=100
    double precision :: esc(N_p),func(N_p)

    double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
    xm_p=1.67d-24


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

return
END SUBROUTINE pop_CImeta

SUBROUTINE pop_OImeta(esc,func,xLd_OImeta)
    !     xnH*rn(OI)*xLd_OI : cooling rate owing to OI per unit volume
    !use COMVAR
    USE lines
    IMPLICIT REAL*8(a-h,o-z)
    !double precision :: xLd_OImeta
    integer, PARAMETER :: N_p=100
    double precision :: esc(N_p),func(N_p)
    double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
    xm_p=1.67d-24

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
return
END SUBROUTINE pop_OImeta

SUBROUTINE thick_lev(J_thick)
   !double precision :: thick_lev
   !USE COMVAR
   USE lines 
   USE mol_lines 
   IMPLICIT REAL*8(a-h,o-z)
   integer, PARAMETER :: N_p=100
   double precision :: f_LTE(0:N_p),A(N_p),f(0:N_p),aa(N_p),P_J(N_p),esc(N_p)
   double precision :: xk_B=1.38d-16,h_Pl=6.63d-27,pi=3.14159265358979d0,&
   xm_p=1.67d-24

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
   !double precision :: beta_esc,tau_L,tau_C
   
   if(tau_L.lt.0.d0) then
      beta_esc=1.d0
   elseif(tau_L.lt.1.d-5) then
      beta_esc=dexp(-tau_C)
   else
      beta_esc=dexp(-tau_C)*(1.d0-dexp(-tau_L))/tau_L
   endif
   
   return
END FUNCTION beta_esc

