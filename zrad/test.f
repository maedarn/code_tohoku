      PROGRAM ZR
C     computes the evolutionary path of 
C     metal polluted protostellar clouds 
C     with radiation
C     updated 2008 Jun
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (N_sp=50) 
      common /UVCR/ zeta, G_0
      common /radtemp/ T_rad
      common /coldens/A_v,xNc_H,xNc_H2,xNc_HD
      DIMENSION y(N_sp)
      COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc  
      common /indcool/ xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,
     &     xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
      common /indcool_1/ xLmbd_H2_1, xLmbd_HD_1, xLmbd_CO_1, xLmbd_OH_1,
     &     xLmbd_H2O_1, xLmbd_CII_1, xLmbd_CI_1, xLmbd_OI_1, xLmbd_Lya_1
      COMMON /radgr/ T_gr_K,tau_cont
      DATA xm_p/1.67d-24/,xk_B/1.38d-16/,G/6.67d-8/,
     &     pi/3.14159265358979d0/
c
c
      func(t,xLmbd)=t-xLmbd
      write(*,*) func(1.d0,1.d0),func(2.d0,1.d0)

      END PROGRAM 

