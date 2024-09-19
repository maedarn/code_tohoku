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
   
   !return
END FUNCTION c_H2



!********************************************************************
!********************************************************************
!**********************  MOLECULAR HYDROGEN  ************************
!********************************************************************
!********************************************************************

function H2_pf(T_K)
   USE mol_lines
   IMPLICIT REAL*8(a-h,o-z)
   double precision :: ET(0:iv_max,0:J_max)
   do iv=0,iv_max
      do j=0,J_max
         ET(iv,j)=E_H2_BFM(iv,j)
      enddo
   enddo
   
   ET_diss=5.196d4
   
   z_p=0.d0
   z_o=0.d0
   do iv=0,iv_max
      do j=0,J_max
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
   IMPLICIT REAL*8(a-h,o-z)
   !     Borysow, Frommhold, and Moraldi (1989) ApJ,336,495
   !     Equation (A11) (in K)
   
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
   
   return
end function E_H2_BFM

        

FUNCTION xLd_H2(xnH,T_K,y_H,y_H2,y_e,y_Hp,xNc_H2,tau_cnt)
   USE PYSCONST
   IMPLICIT REAL*8(a-h,o-z)
   !     xnH*xn(H2)*xLd_H2 : cooling rate owing to H2 per unit volume
   double precision :: A(0:2,0:2,0:22,0:22),ET(0:2,0:22),p(0:22),&
   f(0:2,0:22),f_v(0:2)
   double precision :: beta_esc,esc
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
            xnu_Hz=DE/h_Pl
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
            xnu_Hz=DE/h_Pl
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
            xnu_Hz=DE/h_Pl
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
   return
END FUNCTION xLd_H2