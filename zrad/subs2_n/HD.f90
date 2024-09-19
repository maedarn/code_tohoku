
FUNCTION xLd_HD(xnH,T_K,y_H,y_H2,xNc_HD,tau_cnt)
   USE PYSCONST
   IMPLICIT REAL*8(a-h,o-z)
   ! HD cooling from Galli and Palla (1998), A&A, 335, 403
   ! cooling is for T < 3000 K, in erg cm^3 s-1
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