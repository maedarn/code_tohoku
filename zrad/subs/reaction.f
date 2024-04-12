
      SUBROUTINE react_coef
     &     (xnH,T_K,T_gr_K,Z_metal,xk)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (N_react=675)        
      common /UVCR/ zeta, G_0
      common /coldens/A_v,xNc_H,xNc_H2,xNc_HD
      common /xcrit/xH,xH2,xHe
      DIMENSION xk(N_react)

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
c
      xlT4=dlog10(T_K/1.d4)
      xncr_H=10.d0**(3.d0-0.416d0*xlT4-0.327d0*xlT4**2)
      xncr_H2=10.d0**(4.845d0-1.3d0*xlT4+1.62d0*xlT4**2)
      xncr_He=10.d0**(5.0792d0*(1.d0-1.23d-5*(T_K-2000.d0)))
c     
      xncr=(xH/xncr_H+xH2/xncr_H2+xHe/xncr_He)**(-1.d0)
      xcr=xnH/xncr
c

C  1)   H     +   e     ->   H+    + 2 e   
c  none(UMIST)
c     GA08 12
      xk(1)=dexp(-32.71396786d0+(13.536556d0
     &     +(-5.73932875d0+(1.56315498d0+(-0.2877056d0
     &     +(3.48255977d-2+(-2.63197617d-3
     &     +(1.11954395d-4-2.03914985d-6*xlnT_eV)*xlnT_eV)
     &     *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
C  2)   H+    +   e     ->   H     +   ph.
c  4009(UMIST) 
c  GA08 13
c   case A
c      xk(2)=1.269d-13*(315614.d0/T_K)**1.503d0
c     &     *(1.d0+(604625.d0/T_K)**0.470d0)**(-1.923d0)
c   case B
      xk(2)=2.753d-14*(315614.d0/T_K)**1.500d0
     &     *(1.d0+(115188.d0/T_K)**0.407d0)**(-2.242d0)
c
C  3)   He    +   e     ->   He+   + 2 e    
c  UMIST none
c  GA08 17
      xk(3)=dexp(-44.09864886d0+(23.91596563d0
     &      +(-10.7532302d0+(3.05803875d0+(-0.56851189d0
     &      +(6.79539123d-2+(-5.00905610d-3+(2.06723616d-4
     &      -3.64916141d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
     &      *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
C  4)   He+   +   e     ->   He    +   ph.        
c  4010
c  GA08 19
c     case A
      xkrrA=1.d-11*T_K**(-0.5d0)*(12.72d0-1.615d0*xlgT
     &     -0.3162d0*xlgT2+0.0493d0*xlgT3)
c     case B
      xkrrB=1.d-11*T_K**(-0.5d0)*(11.19d0-1.676d0*xlgT
     &     -0.2852d0*xlgT2+0.04433d0*xlgT3)
c     dielectric
      xkdi=1.9d-3*T_K**(-1.5d0)*dexp(-473421.d0/T_K)
     &     *(1.d0+0.3d0*dexp(-94684.d0/T_K))
c     optically thick
      xkthick=0.68d0*xkrrA+0.32d0*xkrrB+xkdi
      xk(4)=xkthick
C  5)   He+   +   e     ->   He++  + 2 e
c  none
c  GA08 18
      xk(5)=dexp(-68.71040990d0+(43.93347633d0+(-18.4806699d0
     &     +(4.70162649d0+(-0.76924663d0+(8.113042d-2
     &     +(-5.32402063d-3+(1.97570531d-4-3.16558106d-6*xlnT_eV)
     &     *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
     &     *xlnT_eV)*xlnT_eV)*xlnT_eV)
C  6)   He++  +   e     ->   He+   +   ph.        
c  none
c  GA08 20
c     case A
      xk6A=2.538d-13*(1262456d0/T_K)**1.503d0
     &     *(1.d0+(2418500.d0/T_K)**0.470d0)**(-1.923d0)
c     case B
      xk6B=5.506d-14*(1262456d0/T_K)**1.500d0
     &     *(1.d0+(460752.d0/T_K)**0.407d0)**(-2.242d0)
c
      xk(6)=xk6B
C  7)   H     +   e     ->   H-    +   ph.           
c  4006,4007,4008
c  GA08 1
      if(T_K < 6000.d0) then 
         xk(7)=10.d0**(-17.845d0+0.762d0*xlgT
     &     +0.1523d0*xlgT2-0.03274d0*xlgT3)
      else
         xk(7)=10.d0**(-16.4199d0+0.1998d0*xlgT2
     &        -5.447d-3*xlgT4+4.0415d-5*xlgT6)
      endif
C  8)   H-    +   H     ->   H2    +   e              
c  4031
c  GA08 2
      xk(8)=1.3d-9
C  9)   H     +   H+    ->   H2+   +   ph.
c  4078, 4079       
c  GA08 3
      xk(9)=1.d1**(-19.38d0-1.523d0*xlgT
     &     +1.118d0*xlgT2-0.1269d0*xlgT3)
C  10)   H2+   +   H     ->   H2    +   H+  
c  2937       
c  GA08 4
      xk(10)=6.4d-10
C  11)   H2    +   H+    ->   H2+   +   H        
c  none
c  GA08 7
      if(T_K < 100.d0) then
         xk(11)=0.d0
      else
         xk(11)=(-3.3232183d-7+3.3735382d-7*xlnT
     &        -1.4491368d-7*(xlnT**2)
     &        +3.4172805d-8*(xlnT**3)-4.7813720d-9*(xlnT**4)
     &        +3.9731542d-10*(xlnT**5)
     &        -1.8171411d-11*(xlnT**6)+3.5311932d-13*(xlnT**7))
     &        *dexp(-21237.15d0/T_K)
      endif
C  12)   H2    +   e     -> 2 H     +   e             
c  4592
c  GA08 8
c     v=0
      xk12v0=4.49d-9*(T_K**0.11d0)*dexp(-101858.d0/T_K)
c     LTE
      xk12LTE=1.91d-9*(T_K**0.136d0)*dexp(-53407.1d0/T_K)
c
      xlgk12=(xcr/(1.d0+xcr))*dlog10(xk12LTE)
     &     +(1.d0/(1.d0+xcr))*dlog10(xk12v0)
c      
      xk(12)=10.d0**xlgk12
C  13)   H2    +   H     -> 3 H              
c  4584
c     GA08 9
c     v=0
      xk13v0=6.67d-12*(T_K**0.5d0)*dexp(-(1.d0+63593.d0/T_K))
c     LTE
      xk13LTE=3.52d-9*dexp(-43900.d0/T_K)
c
      xlgk13=(xcr/(1.d0+xcr))*dlog10(xk13LTE)
     &     +(1.d0/(1.d0+xcr))*dlog10(xk13v0)
c      
      xk(13)=10.d0**xlgk13
c      
C  14)   H-    +   e     ->   H     + 2 e             
c   none
c     GA08 (14)
      xk(14)=dexp(-18.01849334d0+(2.3608522d0
     &      +(-0.28274430d0+(1.62331664d-2+(-3.36501203d-2
     &      +(1.17832978d-2+(-1.65619470d-3+(1.06827520d-4
     &      -2.63128581d-6*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
     &      *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
C  15)   H-    +   H+    -> 2 H                
c   3489
c   GA08 (5)
      xk(15)=2.4d-6*(T_K**(-0.5d0))*(1.d0+T_K/20000.d0)
C  16)   H-    +   H+    ->   H2+   +   e          
c   none
c   GA08 (16)
      if(T_K.le.8.d3) then
         xk(16)=6.9d-9*(T_K**(-0.35d0))
      else
         xk(16)=9.6d-7*(T_K**(-0.90d0))
      endif
C  17)   H2+   +   e     -> 2 H                
c   3520
c   GA08 (6)
      if(T_K < 617.d0) then 
         xk(17)=1.d-8
      else
         xk(17)=1.32d-6*T_K**(-0.76d0)
      endif
C  18)   H2+   +   H-    ->   H2    +   H         
c   3490
c   GA08 (21)
      xk(18)=1.4d-7*T300**(-0.5d0)
C  19) 3 H               ->   H2    +   H
c   none               
c   Glover 08 (GA08;30) 
      xk(19)=7.7d-31*T_K**(-0.464d0)
C  20) 2 H     +   H2    -> 2 H2                
c  none
c   GA08 (31)
      xk(20)=xk(19)/8.d0
C  21) 2 H2              -> 2 H     +   H2            
c  4593
c   GA08 (10)
c     v=0
      xk21v0=(5.996d-30*T_K**4.1881d0/(1.d0+6.761d-6*T_K)**5.6881d0)
     &     *dexp(-54657.4d0/T_K)
c     LTE
      xk21LTE=1.3d-9*dexp(-53300.d0/T_K)
c
      xlgk21=(xcr/(1.d0+xcr))*dlog10(xk21LTE)
     &     +(1.d0/(1.d0+xcr))*dlog10(xk21v0)
c      
      xk(21)=10.d0**xlgk21
c
C  22)
c   none
      xk(22)=0.d0
C  23) 2 H     +   grain ->   H2               
c   none
C     Tielens & Hollenbach (1985) 
      if(T_gr_K.eq.0.d0) then
         f_a=1.d0
      else
         f_a=1.d0/(1.d0+dexp(7.5d2*(1.d0/7.5d1-1.d0/T_gr_K)))
      endif
      xk(23)=6.d-17*dsqrt(T_K/3.d2)*Z_metal*f_a
     &     /(1.d0+4.d-2*dsqrt(T_K+T_gr_K)
     &     +2.d-3*T_K+8.d-6*T_K**2)
C  24)   He+   +   H2    ->   H+    +   H     +   He
c   642
c   GA08 (24)
      xk(24)=3.70d-14*dexp(-35.d0/T_K)
C  25)   H2+   +   He    ->   HeH+  +   H
c   641
      xk(25)=1.30d-10
C  26)   H2+   +   H2    ->   H3+   +   H
c   640
      xk(26)=2.08d-9
C  27)   H3+   +   H-    -> 2 H2
c   4601
      xk(27)=2.30d-7*(T_K/300.d0)**(-0.50d0)
C  28)   He+   +   H     ->   H+    +   He
c   2938, 2939
c   GA08 (26)
      xk(28)=1.2d-15*dexp(T_K/300.d0)**0.25d0
C  29)   He+   +   H-    ->   H     +   He
c   3491
c   GA08 (28)
       xk(29)=2.32d-7*(T_K/300.d0)**(-0.52d0)
     &     *dexp(T_K/22400.d0)
C  30)   He+   +   H2    ->   H2+   +   He
c   3060
c   GA08 (25)
      xk(30)=7.20d-15
C  31)   HeH+  +   H     ->   H2+   +   He
c   550
      xk(31)=9.10d-10
C  32)   HeH+  +   H2    ->   H3+   +   He
c   643
      xk(32)=1.50d-9
C  33)
c   none
      xk(33)=0.d0
C  34)   H3+   +   e     ->   H2    +   H
c   3522
      xk(34)=2.34d-8*(T_K/300.d0)**(-0.52d0)
C  35)   H3+   +   e     -> 3 H
c   3521
      xk(35)=4.36d-8*(T_K/300.d0)**(-0.52d0)
C  36)   HeH+  +   e     ->   H     +   He
c   3523
      xk(36)=1.00d-8*(T_K/300.d0)**(-0.60d0)
c  37)   H2    +   He    -> 2 H     +   He
c  GA08 (11)
c     v=0
      xk37v0=10.d0**(-27.029d0+3.801d0*xlgT-29487.d0/T_K)
c     LTE
      xk37LTE=10.d0**(-2.729d0-1.75d0*xlgT-23474.d0/T_K)
c
      xlgk37=(xcr/(1.d0+xcr))*dlog10(xk37LTE)
     &     +(1.d0/(1.d0+xcr))*dlog10(xk37v0)
c      
      xk(37)=10.d0**xlgk37
c
c  38)   H-    +   H     -> 2 H     +   e
c  GA08 (15)
      if(T_eV < 0.1d0) then
         xk(38)=2.5634d-9*T_eV**1.78186d0
      else
         xk(38)=dexp(-20.372609d0+(1.13944933d0
     &      +(-1.4210135d-1+(8.4644554d-3+(-1.4327641d-3
     &      +(2.0122503d-4+(8.6639632d-5+(-2.5850097d-5
     &      +(2.4555012d-6-8.0683825d-8*xlnT_eV)
     &        *xlnT_eV)*xlnT_eV)*xlnT_eV)
     &      *xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)*xlnT_eV)
      endif
c  39)   H-    +   H2+   -> 3 H
c  GA08 (22)
      xk(39)=1.4d-7*T300**(-0.5d0)
c  40)   H2    +   e     ->   H-    +   H
c  GA08 (23)
      xk(40)=2.7d-8*(T_K**(-1.27d0))*dexp(-43000.d0/T_K)
c  41)   He    +   H+    ->   He+   +   H
c  GA08 (27)
      if(T_K < 1.d4) then         
         xk(41)=1.26d-9*(T_K**(-0.75d0))*dexp(-127500.d0/T_K)
      else
         xk(41)=4.0d-37*T_K**4.74d0
      endif
c  42)   He    +   H-    ->   He    +   H   +   e
c  GA08 (29)
      xk(42)=4.1d-17*(T_K**2)*dexp(-19870.d0/T_K)
c  43) 2 H     +   He    ->   H2    +   He
c  GA08 (32)
      xk(43)=6.9d-32*T_K**(-0.4d0)
c
      xk(44:46)=0.d0
c      
c     47-50: used for additional D reactions 
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     D    reactions                           C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  51)   D+  +   e    ->   D    +   ph. 
c  GA08 (33)
      xk(51)=xk(2)
C  52)   D   +   H+   ->   D+   +   H   
c  GA08 (34)
      if(T_K < 2.d5) then
         xk(52)=2.0d-10*(T_K**0.402d0)*dexp(-37.1d0/T_K)
     &        -3.31d-17*(T_K**1.48d0)
      else
         xk(52)=3.44d-10*(T_K**0.35d0)
      endif
C  53)   D+  +   H    ->   D    +   H+ 
c  GA08 (35)
      xk(53)=2.06d-10*(T_K**0.396d0)*dexp(-33.d0/T_K)
     &     +2.03d-9*(T_K**(-0.332d0))
C  54)   D   +   H    ->   HD   +   ph. 
c  GA08 (36)
      if(T_K < 200.d0) then 
         xk(54)=1.d-25*(2.80202d0-6.63697d0*xlnT
     &        +4.75619d0*xlnT**2-1.39325d0*xlnT**3
     &        +0.178259d0*xlnT**4-0.00817097d0*xlnT**5)
      else
         xk(54)=1.d-25*dexp(507.207d0-370.889d0*xlnT
     &        +104.854d0*xlnT**2-14.4192d0*xlnT**3
     &        +0.971469d0*xlnT**4-0.0258076d0*xlnT**5)
      endif
C  55)   D   +   H2   ->   H    +   HD 
c  GA08 (37)
      if(T_K < 2000.d0) then 
         xk(55)=10.d0**(-56.4737d0+5.88886d0*xlgT
     &        +7.19692d0*xlgT**2+2.25069d0*xlgT**3
     &        -2.16903d0*xlgT**4+0.317887d0*xlgT**5)
      else
         xk(55)=3.17d-10*dexp(-5207.d0/T_K)
      endif
C  56)   HD+ +   H    ->   H+   +   HD   
c  GA08 (38)
      xk(56)=xk(10)
C  57)   D+  +   H2   ->   H+   +   HD  
c  GA08 (39)
      xk(57)=(0.417d0+0.846d0*xlgT-0.137d0*(xlgT**2))*1.d-9
C  58)   HD  +   H    ->   H2   +   D
c  GA08 (40)
      if(T_K < 200.d0) then
         xk(58)=5.25d-11*dexp(-4430.d0/T_K)
      else
         xk(58)=5.25d-11*dexp(-4430.d0/T_K+173900.d0/T_K**2)
      endif
C  59)   HD  +   H+   ->   H2   +   D+ 
c  GA08 (41)
      xk(59)=1.1d-9*dexp(-488.d0/T_K)
C  60)   D   +   H+   ->   HD+  +   ph.
c  GA08 (42)
      xk(60)=3.9d-19*(T300**1.8d0)*dexp(20.d0/T_K)
C  61)   D+  +   H    ->   HD+  +   ph.   
c  GA08 (43)
      xk(61)=3.9d-19*(T300**1.8d0)*dexp(20.d0/T_K)
C  62)   HD+ +   e    ->   H    +   D 
c  GA08 (44)
      xk(62)=7.2d-8*T_K**(-0.5d0)
C  63)   D   +   e    ->   D-   +   ph.    
c  GA08 (51)
      xk(63)=xk(7)
C  64)   D+  +   D-   ->   2D 
c  GA08 (68)
      xk(64)=xk(15)
C  65)   H+  +   D-   ->   D    +   H 
c  GA08 (67)
      xk(65)=xk(15)
C  66)   H-  +   D    ->   H    +   D-   
c  GA08 (53)
      xk(66)=6.4d-9*(T300**0.41d0)
C  67)   D-  +   H    ->   D    +   H-  
c  GA08 (52)
      xk(67)=6.4d-9*(T300**0.41d0)
C  68)   D-  +   H    ->   HD   +   e
c  GA08 (55)
      xk(68)=0.5d0*xk(8)
C  69)   D   +   e    ->   D+   +  2 e
c  GA08 (45)
      xk(69)=xk(1)
C  70)   He+ +   D    ->   D+   +   He
c  GA08 (46)
      xk(70)=1.1d-15*T300**0.25d0
C  71)   He  +   D+   ->   D    +   He+
c  GA08 (47)
      if(T_K < 10000.d0) then
         xk(71)=1.85d-9*(T_K**(-0.75d0))*dexp(-127500.d0/T_K)
      else
         xk(71)=5.9d-37*T_K**4.74d0
      endif
C  72)   H2+ +   D    ->  HD+   +   H
c  GA08 (48)
      xk(72)=1.07d-9*(T300**0.062d0)*dexp(-T_K/41400.d0)
C  73)   HD+ +   D    ->  HD    +   D+
c  GA08 (49)
      xk(73)=xk(10)
C  74)   HD+ +   H    ->  H2+   +   D
c  GA08 (50)
      xk(74)=1.0d-9*dexp(-154.d0/T_K)
C  75)   HD  +   e    ->  H     +   D-
c  GA08 (57)
      xk(75)=1.35d-9*(T_K**(-1.27d0))*dexp(-43000.d0/T_K)
C  76)   HD  +   e    ->  D     +   H-
c  GA08 (58)
      xk(76)=1.35d-9*(T_K**(-1.27d0))*dexp(-43000.d0/T_K)
C  77)   H+  +   D-   ->  HD+   +   e
c  GA08 (60)
      xk(77)=1.1d-9*T300**(-0.4d0)
C  78)   D+  +   H-   ->  HD+   +   e
c  GA08 (61)
      xk(78)=1.1d-9*T300**(-0.4d0)
C  79)   D-  +   e    ->  D     + 2 e
c  GA08 (63)
      xk(79)=xk(14)
C  80)   D-  +   H    ->  D     +   H    +   e
c  GA08 (64)
      xk(80)=xk(38)
C  81)   D-  +   He   ->  D     +   He   +   e
c  GA08 (65)
      xk(81)=1.5d-17*(T_K**2)*dexp(-19870.d0/T_K)
C  82)   D+  +   H-   ->  D     +   H
c  GA08 (66)
      xk(82)=xk(15)
C  83)   H2+ +   D-   ->  H2    +   D
c  GA08 (69)
      xk(83)=1.7d-7*T300**(-0.5d0)
C  84)   H2+ +   D-   -> 2 H    +   D
c  GA08 (70)
      xk(84)=1.7d-7*T300**(-0.5d0)
C  85)   HD+ +   H-   ->   HD   +   H
c  GA08 (71)
      xk(85)=1.5d-7*T300**(-0.5d0)
C  86)   HD+ +   H-   ->   D    + 2 H
c  GA08 (72)
      xk(86)=1.5d-7*T300**(-0.5d0)
C  87)   HD+ +   D-   ->  HD    +   D
c  GA08 (73)
      xk(87)=1.9d-7*T300**(-0.5d0)
C  88)   HD+ +   D-   ->  2 D   +   H
c  GA08 (74)
      xk(88)=1.9d-7*T300**(-0.5d0)
C  89)   He+ +   D-   ->  He    +   D
c  GA08 (79)
      xk(89)=3.03d-7*(T300**(-0.52d0))*dexp(T_K/22400.d0)
C  90)   D   +   H2+  ->  H2    +   D+
c  GA08 (81)
      xk(90)=xk(10)
C  91)  H2+  +   D    ->  HD    +   H+
c  GA08 (82)
      xk(91)=1.0d-9
C  92)  HD+  +   H    ->  H2    +   D+
c  GA08 (83)
      xk(92)=1.0d-9
C  93)  H2   +   D+   ->  H2+   +   D
c  GA08 (90)
      xk(93)=xk(11)
C  94)  H2   +   D+   ->  HD+   +   H
c  GA08 (91)
      xk(94)=(1.04d-9+9.52d-9*(T_K/10000.d0) 
     &     -1.81d-9*(T_K/10000.d0)**2)*dexp(-21000.d0/T_K)
C  95)  HD   +   H+   ->  HD+   +   H
c  GA08 (92)
      xk(95)=xk(11)
C  96)  HD   +   H+   ->  H2+   +   D
c  GA08 (93)
      xk(96)=1.0d-9*dexp(-21600.d0/T_K)
C  97)  HD   +   D+   ->  HD+   +   D
c  GA08 (94)
      xk(97)=xk(11)
C  98)  HD   +   He+  ->  HD+   +   He
c  GA08 (101)
      xk(98)=xk(30)
C  99)  HD   +   He+  ->  He    +   H+   +   D
c  GA08 (102)
      xk(99)=1.85d-14*dexp(-35.d0/T_K)
C 100)  HD   +   He+  ->  He    +   H    +   D+
c  GA08 (103)
      xk(100)=1.85d-14*dexp(-35.d0/T_K)
C D50)  HD   +   H    -> 2 H    +   D       
c  GA08 (108)
      xk(50)=xk(13)
C D49)  HD   +   H2   ->  H     +   D    +   H2   
c  GA08 (109)
      xk(49)=xk(21)
C D48)  HD   +   He   ->  H     +   D    +   He
c  GA08 (110)
      xk(48)=xk(37)
C D47)  HD   +   e    ->  H     +   D    +   e
c  GA08 (111)
c     v=0
      xk47v0=5.09d-9*(T_K**0.128d0)*dexp(-103258.d0/T_K)
c     LTE
      xk47LTE=1.04d-9*(T_K**0.218d0)*dexp(-53070.7d0/T_K)
c
      xlgk47=(xcr/(1.d0+xcr))*dlog10(xk47LTE)
     &     +(1.d0/(1.d0+xcr))*dlog10(xk47v0)
c      
      xk(47)=10.d0**xlgk47
c
c
c     metal reactions
c
C  101)   H     +   CH    ->   C     +   H2
c   1
      xk(101)= 1.31d-10*dexp(-80.d0/T_K)
C  102)   H     +   CH    ->   C     + 2 H
c   4585
      xk(102)=6.00d-09*dexp(-40200.d0/T_K)
      
C  103)   H     +   CH2   ->   CH    +   H2
c   2
      xk(103)=6.64d-11
C  104)   H     +   CH3   ->   CH2   +   H2
c   4
      xk(104)=1.00d-10*dexp(-7600.d0/T_K)
C  105)   H     +   CH4   ->   H2    +   CH3
c   7
      xk(105)=5.94d-13*(T_K/300.d0)**3.d0*dexp(-4045.d0/T_K)
C  106)   H     +   OH    ->   H2    +   O
c   8
      xk(106)=6.99d-14*(T_K/300.d0)**2.80d0*dexp(-1950.d0/T_K)
C  107)   H     +   OH    ->   O     + 2 H
c   4586
      xk(107)=6.00d-9*dexp(-50900.d0/T_K)
C  108)   H     +   H2O   ->   OH    +   H2
c   10
      xk(108)=1.59d-11*(T_K/300.d0)**1.2d0*dexp(-9610.d0/T_K)
C  109)   H     +   H2O   ->   OH    + 2 H
c   4587
      xk(109)=5.80d-9*dexp(-52900.d0/T_K)
C  110)   H     +   C2    ->   CH    +   C
c   11
      xk(110)=4.67d-10*(T_K/300.d0)**0.50d0*dexp(-30450.d0/T_K)
C  111)   H     +   CO    ->   C     +   OH
c   15
      xk(111)=1.10d-10*(T_K/300.d0)**0.50d0*dexp(-77700.d0/T_K)
C  112)   H     +   H2CO  ->   HCO   +   H2
c   21
      xk(112)=4.85d-12*(T_K/300.d0)**1.9d0*dexp(-1379.d0/T_K)
C  113)   H     +   O2    ->   OH    +   O
c   25
      xk(113)=2.61d-10*dexp(-8156.d0/T_K)
C  114)   H     +   O2    -> 2 O     +   H
c   4591
      xk(114)=6.00d-9*dexp(-52300.d0/T_K)
C  115)   H     +   O2H   ->   H2O   +   O
c   27
      xk(115)=5.00d-11*dexp(-866.d0/T_K)
C  116)   H     +   O2H   ->   H2    +   O2
c   28
      xk(116)=2.06d-11*(T_K/300.d0)**0.84d0*dexp(-277.d0/T_K)
C  117)   H     +   O2H   -> 2 OH
c   26
      xk(117)=1.66d-10*dexp(-413.d0/T_K)
C  118)   H     +   H2O2  ->   O2H   +   H2
c   31
      xk(118)=8.00d-11*dexp(-4000.d0/T_K)
C  119)   H     +   CO2   ->   OH    +   CO
c   36
      xk(119)=3.38d-10*dexp(-13163.d0/T_K)
C  120)   C     +   H2    ->   CH    +   H
c   47
      xk(120)=6.64d-10*dexp(-11700.d0/T_K)
C  121)   C     +   OH    ->   CH    +   O
c   72
      xk(121)=2.25d-11*(T_K/300.d0)**0.50d0*dexp(-14800.d0/T_K)
C  122)   C     +   CO    ->   C2    +   O
c   79
      xk(122)=2.94d-11*(T_K/300.d0)**0.50d0*dexp(-58025.d0/T_K)
C  123)   O     +   H2    ->   OH    +   H
c   53
      xk(123)=3.14d-13*(T_K/300.d0)**2.7d0*dexp(-3150.d0/T_K)
C  124)   O     +   CH    ->   OH    +   C
c  138
      xk(124)=2.52d-11*dexp(-2381.d0/T_K)
C  125)   O     +   CH2   ->   CH    +   OH
c  240 
      xk(125)=4.98d-10*dexp(-6000.d0/T_K)
C  126)   O     +   CH2   ->   HCO   +   H
c  243
      xk(126)=5.01d-11
C  127)   O     +   CH4   ->   CH3   +   OH
c  311
      xk(127)=2.29d-12*(T_K/300.d0)**2.2d0*dexp(-3820.d0/T_K)
C  128)   O     +   H2O   -> 2 OH   
c  314
      xk(128)=1.85d-11*(T_K/300.d0)**0.95d0*dexp(-8571.d0/T_K)
C  129)   O     +   H2CO  ->   HCO   +   OH
c  339
      xk(129)=1.07d-11*(T_K/300.d0)**1.17d0*dexp(-1242.d0/T_K)
C  130)   O     +   H2O2  ->   O2H   +   OH
c  354
      xk(130)=8.54d-14*(T_K/300.d0)**3.25d0*dexp(-1200.d0/T_K)
C  131)   O     +   CO2   ->   CO    +   O2
c  370
      xk(131)=2.46d-11*dexp(-26567.d0/T_K)
C  132)   H+    +   O     ->   O+    +   H
c  2945, 2946  
      if(T_K <1.d4) then
         xk(132)=7.31d-10*(T_K/300.d0)**0.23d0*dexp(-225.9d0/T_K)
      else
         xk(132)=3.04d-10*(T_K/300.d0)**0.47d0*dexp(11.5d0/T_K)
      endif
C  133)   H2    +   CH    ->   CH2   +   H
c   48
      xk(133)=5.46d-10*dexp(-1943.d0/T_K)
C  134)   H2    +   CH    ->   C     +   H   +   H2
c   4594
      xk(134)=6.00d-9*dexp(-40200.d0/T_K)
C  135)   H2    +   CH2   ->   CH3   +   H
c   50
      xk(135)=5.18d-11*(T_K/300.d0)**0.17d0*dexp(-6400.d0/T_K)
C  136)   H2    +   CH3   ->   CH4   +   H
c   52
      xk(136)=6.86d-14*(T_K/300.d0)**2.74d0*dexp(-4740.d0/T_K)
C  137)   H2    +   OH    ->   H2O   +   H
c   55
      xk(137)=2.05d-12*(T_K/300.d0)**1.52d0*dexp(-1736.d0/T_K)
C  138)   H2    +   OH    ->   O     +   H     +   H2
c   4595
      xk(138)=6.00d-9*dexp(-50900.d0/T_K)
C  139)   H2    +   H2O   ->   OH    +   H     +   H2
c   4596
      xk(139)=5.80d-9*dexp(-52900.d0/T_K)
C  140)   H2    +   O2    ->   O2H   +   H
c   59
      xk(140)=2.40d-10*dexp(-28500.d0/T_K)
C  141)   H2    +   O2    -> 2 OH
c   58
      xk(141)=3.16d-10*dexp(-21890.d0/T_K)
C  142)   H2    +   O2    -> 2 O     +   H2
c   4598
      xk(142)=6.00d-9*dexp(-52300.d0/T_K)
C  143)   H2    +   O2H   ->   H2O2  +   H
c   61
      xk(143)=4.38d-12*dexp(-10751.d0/T_K)
C  144) 
c   none
      xk(144)=0.d0
C  145)   H3+   +   O2    ->   O2H+  +   H2
c   788
      xk(145)=9.30d-10*dexp(-100.d0/T_K)
C  146)   C+    +   H2    ->   CH+   +   H
c   645
      xk(146)=1.00d-10*dexp(-4640.d0/T_K)
C  147)   CH    +   CH4   ->   CH2   +   CH3
c   141
      xk(147)=2.28d-11*(T_K/300.d0)**0.70d0*dexp(-3000.d0/T_K)
C  148)   CH    +   OH    ->   HCO   +   H
c   143
      xk(148)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-5000.d0/T_K)
C  149)   CH    +   HCO   ->   CH2   +   CO
c   147
      xk(149)=2.87d-12*(T_K/300.d0)**0.70d0*dexp(-500.d0/T_K)
C  150)   CH    +   H2CO  ->   CH2   +   HCO
c   152
      xk(150)=9.21d-12*(T_K/300.d0)**0.70d0*dexp(-2000.d0/T_K)
C  151)   CH    +   O2    ->   HCO   +   O
c   155
      xk(151)=1.44d-11*(T_K/300.d0)**0.70d0*dexp(-3000.d0/T_K)
C  152)   CH    +   O2H   ->   CH2   +   O2
c   159
      xk(152)=2.94d-13*(T_K/300.d0)**0.50d0*dexp(-7550.d0/T_K)
C  153)   CH    +   O2H   ->   HCO   +   OH
c   160
      xk(153)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
C  154)   CH    +   CO2   ->   HCO   +   CO
c   161
      xk(154)=2.94d-13*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
C  155) 2 CH2             ->   CH    +   CH3
c   236
      xk(155)=4.00d-10*dexp(-5000.d0/T_K)
C  156)   CH2   +   CH4   -> 2 CH3
c   244
      xk(156)=7.13d-12*dexp(-5050.d0/T_K)
C  157)   CH2   +   HCO   ->   CH3   +   CO
c   251
      xk(157)=3.00d-11
C  158)   CH2   +   H2CO  ->   CH3   +   HCO
c   254
      xk(158)=3.30d-13*dexp(-3270.d0/T_K)
C  159)   CH2+  +   H     ->   CH+   +   H2
c   553
      xk(159)=1.00d-9*dexp(-7080.d0/T_K)
C  160)   CH3   +   H2CO  ->   CH4   +   HCO
c   299
      xk(160)=1.34d-15*(T_K/300.d0)**5.05d0*dexp(-1636.d0/T_K)
C  161)   CH3+  +   H     ->   CH2+  +   H2
c   554
      xk(161)=7.00d-10*dexp(-10560.d0/T_K)
C  162)   OH    +   CH2   ->   CH3   +   O
c   245
      xk(162)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
C  163)   OH    +   CH2   ->   H2O   +   CH
c   246
      xk(163)=1.44d-11*(T_K/300.d0)**0.50d0*dexp(-3000.d0/T_K)
C 164)   OH    +   CH2   ->   H2CO  +   H
c  247
      xk(164)=3.00d-11
C 165)   OH    +   CH3   ->   H2O   +   CH2
c  290
      xk(165)=1.20d-10*dexp(-1400.d0/T_K)
C 166)   OH    +   CH4   ->   H2O   +   CH3
c  422
      xk(166)=3.77d-13*(T_K/300.d0)**2.42d0*dexp(-1162.d0/T_K)
C 167) 2 OH              ->   H2O   +   O
c  426
      xk(167)=1.65d-12*(T_K/300.d0)**1.14d0*dexp(-50.d0/T_K)
C 168)   OH    +   CO    ->   CO2   +   H
c  437
      xk(168)=2.81d-13*dexp(-176.d0/T_K)
C 169)   OH    +   H2O2  ->   H2O   +   O2H
c  449
      xk(169)=5.26d-12*dexp(-307.d0/T_K)
C 170)   H2O   +   CH3   ->   CH4   +   OH
c  294
      xk(170)=2.30d-15*(T_K/300.d0)**3.47d0*dexp(-6681.d0/T_K)
C 171)   CO    +   O2    ->   CO2   +   O
c  504
      xk(171)=5.99d-12*dexp(-24075.d0/T_K)
C 172)   CO    +   O2H   ->   CO2   +   OH
c  505
      xk(172)=5.60d-10*dexp(-12160.d0/T_K)
C 173)   O2    +   CH2   ->   HCO   +   OH
c  258
      xk(173)=4.10d-11*dexp(-750.d0/T_K)
C 174)   O2    +   CH2   ->   H2CO  +   O
c  259
      xk(174)=3.65d-11*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
C 175)   O2    +   CH3   ->   O2H   +   CH2
c  303
      xk(175)=5.30d-12*dexp(-34975.d0/T_K)
C 176)   O2    +   CH3   ->   H2CO  +   OH
c  302
      xk(176)=5.64d-13*dexp(-4500.d0/T_K)
C 177)   O2    +   CH4   ->   CH3   +   O2H
c  424
      xk(177)=6.70d-11*dexp(-28640.d0/T_K)
C 178)   O2    +   HCO   ->   O2H   +   CO
c  520
      if(T_K > 200.d0) then
         xk(178)=1.58d-12*(T_K/300.d0)**1.24d0*dexp(353.d0/T_K)
      else
         xk(178)=5.58d-12
      endif
C 179)   O2H   +   CH3   ->   CH4   +   O2
c  306
      xk(179)=6.00d-12
C 180)   O2H   +   H2O   ->   H2O2  +   OH
c  460
      xk(180)=4.65d-11*dexp(-16500.d0/T_K)
C 181)   O2H   +   HCO   ->   H2CO  +   O2
c  521
      xk(181)=5.00d-11
C 182)   O2H   +   H2CO  ->   H2O2  +   HCO
c  531
      xk(182)=3.30d-12*dexp(-5870.d0/T_K)
C 183) 2 O2H             ->   H2O2  +   O2   
c  540  
      if(T_K > 200.d0) then
         xk(183)=3.81d-15*(T_K/300.d0)**3.71d0*dexp(1761.d0/T_K)
      else
         xk(183)=5.64d-12
      endif
C 184)   H     +   HCO   ->   CO    +   H2
c  18
      xk(184)=2.00d-10
C 185)   C     +   H     ->   CH    +   ph.
c  4081
      xk(185)=1.00d-17
C 186) 2 C               ->   C2    +   ph.
c  4106
      xk(186)=4.36d-18*(T_K/300.d0)**0.35d0*dexp(-161.3d0/T_K)
C 187)   C     +   O     ->   CO    +   ph.
c  4112, 4113
      if(T_K > 300.d0) then
         xk(187)=3.09d-17*(T_K/300.d0)**0.33d0*dexp(-1629.d0/T_K)
      else
         xk(187)=2.10d-19
      endif
C 188)   C     +   CH    ->   C2    +   H
c  63
      xk(188)=6.59d-11
C 189)   C     +   OH    ->   CO    +   H
c  73
      xk(189)=1.00d-10
C 190)   C     +   HCO   ->   CH    +   CO
c  83
      xk(190)=1.00d-10
C 191)   C     +   O2    ->   CO    +   O    
c  92, 93
      if(T_K < 295.d0) then
         xk(191)=4.70d-11*(T_K/300.d0)**(-0.34d0)
      else
         xk(191)=2.48d-12*(T_K/300.d0)**1.54d0*dexp(613.d0/T_K)
      endif
C 192)   CH3   +   HCO   ->   CH4   +   CO
c  297
      xk(192)=2.00d-10
C 193)   O     +   H     ->   OH    +   ph.
c  4083
      xk(193)=9.90d-19*(T_K/300.d0)**(-0.38d0)
C 194) 2 O               ->   O2    +   ph.
c  4136
      xk(194)=4.90d-20*(T_K/300.d0)**1.58d0
C 195)   O     +   CH    ->   CO    +   H
c  139, 140
      if(T_K < 2000.d0) then
         xk(195)=6.60d-11
      else
         xk(195)=1.02d-10*dexp(-914.d0/T_K)
      endif
C 196)   O     +   CH    ->   HCO+  +   e
c  4600
      xk(196)=2.00d-11*(T_K/300.d0)**0.44d0
C 197)   O     +   CH2   ->   CO    + 2 H
c  241
      xk(197)=1.33d-10
C 198)   O     +   CH3   ->   H2CO  +   H
c  286
      xk(198)=1.30d-10
C 199)   O     +   OH    ->   O2    +   H
c  312
      if(T_K > 158.d0) then 
         xk(199)=1.77d-11*dexp(178.d0/T_K)
      else
         xk(199)=5.46d-11
      endif
C 200)   O     +   C2    ->   CO    +   C
c  315, 316
      if(T_K < 300.d0) then
         xk(200)=1.d-10
      else
         xk(200)=6.d-10
      endif
C 201)   O     +   HCO   ->   CO2   +   H
c  332
      xk(201)=5.00d-11
C 202)   O     +   HCO   ->   OH    +   CO
c  333
      xk(202)=5.00d-11
C 203)   O     +   O2H   ->   OH    +   O2
c  351
      if(T_K > 200.d0) then
         xk(203)=3.17d-11*dexp(174.d0/T_K)
      else
         xk(203)=7.57d-11
      endif
C 204)   OH    +   HCO   ->   H2O   +   CO
c  439
      xk(204)=1.70d-10
C 205)   OH    +   H2CO  ->   H2O   +   HCO
c  443
      if(T_K > 200.d0) then
         xk(205)=2.22d-12*(T_K/300.d0)**1.42d0*dexp(416.d0/T_K)
      else
         xk(205)=9.99d-12
      endif
C 206)   OH    +   O2H   ->   H2O   +   O2
c  448
      if(T_K > 200.d0) then
         xk(206)=3.66d-11*(T_K/300.d0)**(-0.13d0)*dexp(244.d0/T_K)
      else
         xk(206)=1.31d-10
      endif
C 207) 2 HCO             ->   H2CO  +   CO
c  516
      xk(207)=3.00d-11
C 208)   H+    +   CH    ->   CH+   +   H
c  2940
      xk(208)=1.90d-9*(T_K/300.d0)**(-0.5d0)
C 209)   H+    +   CH2   ->   CH+   +   H2
c  552
      xk(209)=1.40d-9
C 210)   H+    +   CH2   ->   CH2+  +   H
c  2941
      xk(210)=1.40d-9
C 211)   H+    +   CH3   ->   CH3+  +   H
c  2943
      xk(211)=3.40d-9
C 212)   H+    +   CH4   ->   CH3+  +   H2
c  555
      xk(212)=2.30d-9
C 213)   H+    +   CH4   ->   CH4+  +   H
c  2948
      xk(213)=1.50d-9
C 214)   H+    +   OH    ->   OH+   +   H
c  2949
      xk(214)=2.10d-9*(T_K/300.d0)**(-0.5d0)
C 215)   H+    +   H2O   ->   H2O+  +   H
c  2951
      xk(215)=6.90d-9*(T_K/300.d0)**(-0.5d0)
C 216)   H+    +   C2    ->   C2+   +   H
c  2952
      xk(216)=3.10d-9
C 217)   H+    +   HCO   ->   CO+   +   H2
c  566
      xk(217)=9.40d-10*(T_K/300.d0)**(-0.5d0)
C 218)   H+    +   HCO   ->   H2+   +   CO
c  567
      xk(218)=9.40d-10*(T_K/300.d0)**(-0.5d0)
C 219)   H+    +   HCO   ->   HCO+  +   H
c  2963
      xk(219)=9.40d-10*(T_K/300.d0)**(-0.5d0)
C 220)   H+    +   H2CO  ->   H2CO+ +   H
c  2967
      xk(220)=2.96d-9*(T_K/300.d0)**(-0.5d0)
C 221)   H+    +   H2CO  ->   HCO+  +   H2
c  575
      xk(221)=3.57d-9*(T_K/300.d0)**(-0.5d0)
C 222)   H+    +   O2    ->   O2+   +   H
c  2972
      xk(222)=2.00d-9
C 223)   H+    +   CO2   ->   HCO+  +   O
c  606
      xk(223)=3.50d-9
C 224)   H-    +   C     ->   CH    +   e
c  4032
      xk(224)=1.00d-9
C 225)   H-    +   O     ->   OH    +   e
c  4039
      xk(225)=1.00d-9
C 226)   H-    +   CH    ->   CH2   +   e
c  4034
      xk(226)=1.00d-10
C 227)   H-    +   CH2   ->   CH3   +   e
c  4036
      xk(227)=1.00d-9
C 228)   H-    +   CH3   ->   CH4   +   e
c  4038
      xk(228)=1.00d-9
C 229)   H-    +   OH    ->   H2O   +   e
c  4042
      xk(229)=1.00d-10
C 230)   H-    +   CO    ->   HCO   +   e
c  4048
      xk(230)=5.00d-11
C 231)   H-    +   HCO   ->   H2CO  +   e
c  4049
      xk(231)=1.00d-9
C 232)   H2+   +   C     ->   CH+   +   H
c  644
      xk(232)=2.40d-9
C 233)   H2+   +   O     ->   OH+   +   H
c  655
      xk(233)=1.50d-9
C 234)   H2+   +   CH    ->   CH+   +   H2
c  3061
      xk(234)=7.10d-10*(T_K/300.d0)**(-0.5d0)
C 235)   H2+   +   CH    ->   CH2+  +   H
c  646
      xk(235)=7.10d-10*(T_K/300.d0)**(-0.5d0)
C 236)   H2+   +   CH2   ->   CH3+  +   H
c  650
      xk(236)=1.00d-9
C 237)   H2+   +   CH2   ->   CH2+  +   H2
c  3062
      xk(237)=1.00d-9
C 238)   H2+   +   CH4   ->   CH4+  +   H2
c  3065
      xk(238)=1.40d-9
C 239)   H2+   +   CH4   ->   CH5+  +   H
c  659
      xk(239)=1.14d-10
C 240)   H2+   +   CH4   ->   CH3+  +   H     +  H2   
c  662
      xk(240)=2.30d-9
C 241)   H2+   +   OH    ->   OH+   +   H2
c  3066
      xk(241)=7.60d-10*(T_K/300.d0)**(-0.5d0)
C 242)   H2+   +   OH    ->   H2O+  +   H
c  663
      xk(242)=7.60d-10*(T_K/300.d0)**(-0.5d0)
C 243)   H2+   +   H2O   ->   H2O+  +   H2
c  3068
      xk(243)=3.90d-9*(T_K/300.d0)**(-0.5d0)
C 244)   H2+   +   H2O   ->   H3O+  +   H
c  668
      xk(244)=3.40d-9*(T_K/300.d0)**(-0.5d0)
C 245)   H2+   +   C2    ->   C2+   +   H2
c  3069
      xk(245)=1.10d-9
C 246)   H2+   +   CO    ->   CO+   +   H2
c  3074
      xk(246)=6.44d-10
C 247)   H2+   +   CO    ->   HCO+  +   H
c  682
      xk(247)=2.16d-9
C 248)   H2+   +   HCO   ->   H3+   +   CO
c  689
      xk(248)=1.00d-9*(T_K/300.d0)**(-0.5d0)
C 249)   H2+   +   HCO   ->   HCO+  +   H2
c  3076
      xk(249)=1.00d-9*(T_K/300.d0)**(-0.5d0)
C 250)   H2+   +   H2CO  ->   H2CO+ +   H2
c  3078
      xk(250)=1.40d-9*(T_K/300.d0)**(-0.5d0)
C 251)   H2+   +   H2CO  ->   HCO+  +   H     +   H2
c  691
      xk(251)=1.40d-9*(T_K/300.d0)**(-0.5d0)
C 252)   H2+   +   O2    ->   O2H+  +   H
c  694
      xk(252)=1.90d-9
C 253)   H2+   +   O2    ->   O2+   +   H2
c  3080
      xk(253)=8.00d-10
C 254)   H2+   +   CO2   ->   HCO2+ +   H
c  713
      xk(254)=2.35d-9
C 255)   H3+   +   C     ->   CH+   +   H2
c  749
      xk(255)=2.00d-9
C 256)   H3+   +   O     ->   OH+   +   H2
c  754
      xk(256)=8.40d-10
C 257)   H3+   +   CH    ->   CH2+  +   H2
c  750
      xk(257)=1.20d-9*(T_K/300.d0)**(-0.5d0)
C 258)   H3+   +   CH2   ->   CH3+  +   H2
c  751
      xk(258)=1.70d-9
C 259)   H3+   +   CH3   ->   CH4+  +   H2
c  753
      xk(259)=2.10d-9
C 260)   H3+   +   CH4   ->   CH5+  +   H2
c  757
      xk(260)=2.40d-9
C 261)   H3+   +   OH    ->   H2O+  +   H2
c  758
      xk(261)=1.30d-9*(T_K/300.d0)**(-0.5d0)
C 262)   H3+   +   H2O   ->   H3O+  +   H2
c  760
      xk(262)=5.90d-9*(T_K/300.d0)**(-0.5d0)
C 263)   H3+   +   CO    ->   HCO+  +   H2
c  771
      xk(263)=1.70d-9
C 264)   H3+   +   HCO   ->   H2CO+ +   H2
c  777
      xk(264)=1.70d-9*(T_K/300.d0)**(-0.5d0)
C 265)   H3+   +   H2CO  ->   H3CO+ +   H2
c  781
      xk(265)=6.30d-9*(T_K/300.d0)**(-0.5d0)
C 266)   H3+   +   CO2   ->   HCO2+ +   H2
c  818
      xk(266)=2.00d-9
C 267)   He+   +   CH    ->   C+    +   H     +   He
c  913
      xk(267)=1.10d-9*(T_K/300.d0)**(-0.5d0)
C 268)   He+   +   CH    ->   CH+   +   He
c  3083
      xk(268)=5.00d-10*(T_K/300.d0)**(-0.5d0)
C 269)   He+   +   CH2   ->   CH+   +   H     +   Hec
c  915
      xk(269)=7.50d-10
C 270)   He+   +   CH2   ->   C+    +   H2    +   He
c  914
      xk(270)=7.50d-10
C 271)   He+   +   CH3   ->   CH+   +   H2    +   He
c  917
      xk(271)=1.80d-9
C 272)   He+   +   CH4   ->   CH3   +   H+    +   He
c  922
      xk(272)=4.80d-10
C 273)   He+   +   CH4   ->   CH+   +   H2    +   He    +   H
c  920
      xk(273)=2.40d-10
C 274)   He+   +   CH4   ->   CH2+  +   H2    +   He
c  921
      xk(274)=9.50d-10
C 275)   He+   +   CH4   ->   CH3+  +   He    +   H
c  923
      xk(275)=8.50d-11
C 276)   He+   +   CH4   ->   CH4+  +   He
c  3084
      xk(276)=5.10d-11
C 277)   He+   +   OH    ->   O+    +   H     +   He
c  924
      xk(277)=1.10d-9*(T_K/300.d0)**(-0.5d0)
C 278)   He+   +   H2O   ->   H+    +   OH    +   He
c  927
      xk(278)=2.04d-10*(T_K/300.d0)**(-0.5d0)
C 279)   He+   +   H2O   ->   OH+   +   H     +   He
c  928
      xk(279)=2.86d-10*(T_K/300.d0)**(-0.5d0)
C 280)   He+   +   H2O   ->   H2O+  +   He 
c  3086
      xk(280)=6.05d-11*(T_K/300.d0)**(-0.5d0)
C 281)   He+   +   C2    ->   C+    +   C     +   He
c  930
      xk(281)=1.60d-9
C 282)   He+   +   C2    ->   C2+   +   He
c  3087
      xk(282)=5.00d-10
C 283)   He+   +   CO    ->   C+    +   O     +   He
c  948, 949      
      xk(283)=1.60d-9
c      xk(283)=1.40d-9*(T_K/300.d0)**(-0.5d0)
C 284)   He+   +   HCO   ->   CO+   +   H     +   He
c  957
      xk(284)=4.90d-10*(T_K/300.d0)**(-0.5d0)
C 285)   He+   +   HCO   ->   CH+   +   O     +   He
c  955
      xk(285)=4.90d-10*(T_K/300.d0)**(-0.5d0)
C 286)   He+   +   HCO   ->   HeH+  +   CO
c  956
      xk(286)=3.00d-10*(T_K/300.d0)**(-0.5d0)
C 287)   He+   +   H2CO  ->   CO+   +   H2    +   He
c  965
      xk(287)=1.88d-9*(T_K/300.d0)**(-0.5d0)
C 288)   He+   +   H2CO  ->   HCO+  +   H     +   He
c  966
      xk(288)=1.14d-9*(T_K/300.d0)**(-0.5d0)
C 289)   He+   +   O2    ->   O+    +   O     +   He
c  977
      xk(289)=1.10d-9
C 290)   He+   +   O2    ->   O2+   +   He
c  3095
      xk(290)=3.30d-11
C 291)   He+   +   CO2   ->   O2+   +   C     +   He
c  1025
      xk(291)=1.10d-11
C 292)   He+   +   CO2   ->   O+    +   CO    +   He
c  1023
      xk(292)=1.00d-10
C 293)   He+   +   CO2   ->   CO+   +   O     +   He
c  1024
      xk(293)=8.70d-10
C 294)   He+   +   CO2   ->   C+    +   O2    +   He
c  1026
      xk(294)=4.00d-11
C 295)   C+    +   H     ->   CH+   +   ph.
c  4082
      xk(295)=1.70d-17
C 296)   C+    +   O     ->   CO+   +   ph.
c  4114, 4115
      if(T_K < 300.d0) then
         xk(296)=2.50d-18
      else
         xk(296)=3.14d-18*(T_K/300.d0)**(-0.15d0)*dexp(-68.d0/T_K)
      endif
C 297)   C+    +   H-    ->   H     +   C
c  3493  
      xk(297)=2.30d-7*(T_K/300.d0)**(-0.50d0)
C 298)   C+    +   H2    ->   CH2+  +   ph.
c  4089
      xk(298)=4.00d-16*(T_K/300.d0)**(-0.20d0)
C 299)   C+    +   CH    ->   C2+   +   H
c  1173
      xk(299)=3.80d-10*(T_K/300.d0)**(-0.50d0)
C 300)   C+    +   CH    ->   CH+   +   C
c  3100
      xk(300)=3.80d-10*(T_K/300.d0)**(-0.50d0)
C 301)   C+    +   CH2   ->   CH2+  +   C
c  3101
      xk(301)=5.20d-10
C 302)   C+    +   OH    ->   CO+   +   H
c  1186
      xk(302)=7.70d-10*(T_K/300.d0)**(-0.50d0)
C 303)   C+    +   H2O   ->   HCO+  +   H
c  1191
      xk(303)=9.00d-10*(T_K/300.d0)**(-0.50d0)
C 304)   C+    +   HCO   ->   HCO+  +   C
c  3111
      xk(304)=4.80d-10*(T_K/300.d0)**(-0.50d0)
C 305)   C+    +   HCO   ->   CH+   +   CO
c  1212
      xk(305)=4.80d-10*(T_K/300.d0)**(-0.50d0)
C 306)   C+    +   H2CO  ->   CH2+  +   CO
c  1225
      xk(306)=2.34d-9*(T_K/300.d0)**(-0.50d0)
C 307)   C+    +   H2CO  ->   HCO+  +   CH
c  1226
      xk(307)=7.80d-10*(T_K/300.d0)**(-0.50d0)
C 308)   C+    +   H2CO  ->   H2CO+ +   C
c  3114
      xk(308)=7.80d-10*(T_K/300.d0)**(-0.50d0)
C 309)   C+    +   O2    ->   CO+   +   O
c  1234
      xk(309)=3.80d-10
C 310)   C+    +   O2    ->   O+    +   CO
c  1235
      xk(310)=6.20d-10
C 311)   C+    +   CO2   ->   CO+   +   CO
c  1287
      xk(311)=1.10d-9
C 312)   CH+   +   H     ->   C+    +   H2
c  551
      xk(312)=7.50d-10
C 313)   CH+   +   C     ->   C2+   +   H
c  1174
      xk(313)=1.20d-9
C 314)   CH+   +   O     ->   CO+   +   H
c  1404
      xk(314)=3.50d-10
C 315)   CH+   +   H2    ->   CH2+  +   H
c  647
      xk(315)=1.20d-9
C 316)   CH+   +   CH    ->   C2+   +   H2
c  1397
      xk(316)=7.40d-10*(T_K/300.d0)**(-0.50d0)
C 317)   CH+   +   OH    ->   CO+   +   H2
c  1412
      xk(317)=7.50d-10*(T_K/300.d0)**(-0.50d0)
C 318)   CH+   +   H2O   ->   HCO+  +   H2
c  1420
      xk(318)=2.90d-9*(T_K/300.d0)**(-0.50d0)
C 319)   CH+   +   H2O   ->   H3O+  +   C
c  1417
      xk(319)=5.80d-10*(T_K/300.d0)**(-0.50d0)
C 320)   CH+   +   H2O   ->   H2CO+ +   H
c  1419
      xk(320)=5.80d-10*(T_K/300.d0)**(-0.50d0)
C 321)
c  none
      xk(321)=0.d0
C 322)   CH+   +   HCO   ->   CH2+  +   CO
c  1440
      xk(322)=4.60d-10*(T_K/300.d0)**(-0.50d0)
C 323)   CH+   +   HCO   ->   HCO+  +   CH
c  3174
      xk(323)=4.60d-10*(T_K/300.d0)**(-0.50d0)
C 324)   CH+   +   H2CO  ->   CH3+  +   CO
c  1444
      xk(324)=9.60d-10*(T_K/300.d0)**(-0.50d0)
C 325)   CH+   +   H2CO  ->   HCO+  +   CH2
c  1446
      xk(325)=9.60d-10*(T_K/300.d0)**(-0.50d0)
C 326)   CH+   +   H2CO  ->   H3CO+ +   C
c  1445
      xk(326)=9.60d-10*(T_K/300.d0)**(-0.50d0)
C 327)   CH+   +   O2    ->   HCO+  +   O
c  1452
      xk(327)=9.70d-10
C 328)   CH+   +   O2    ->   HCO   +   O+
c  1453
      xk(328)=1.00d-11            
C 329)   CH+   +   O2    ->   CO+   +   OH
c  1451
      xk(329)=1.00d-11      
C 330)   CH+   +   CO2   ->   HCO+  +   CO
c  1465
      xk(330)=1.60d-9
C 331)   CH2+  +   O     ->   HCO+  +   H
c  1554
      xk(331)=7.50d-10
C 332)   CH2+  +   H2    ->   CH3+  +   H
c  651
      xk(332)=1.60d-9
C 333)   CH2+  +   H2O   ->   H3CO+ +   H
c  1563
      xk(333)=1.20d-9*(T_K/300.d0)**(-0.50d0)
C 334)   CH2+  +   HCO   ->   CH3+  +   CO
c  1577
      xk(334)=4.50d-10*(T_K/300.d0)**(-0.50d0)
C 335)   CH2+  +   H2CO  ->   HCO+  +   CH3
c  1581
      xk(335)=2.81d-9*(T_K/300.d0)**(-0.50d0)
C 336)   CH2+  +   O2    ->   HCO+  +   OH
c  1586
      xk(336)=9.10d-10
C 337)   CH2+  +   CO2   ->   H2CO+ +   CO
c  1593
      xk(337)=1.60d-9
C 338)   CH3+  +   O     ->   HCO+  +   H2
c  1649
      xk(338)=4.00d-10
C 339)   CH3+  +   O     ->   H2CO+ +   H
c  1648
      xk(339)=4.00d-11
C 340)   CH3+  +   H2    ->   CH5+  +   ph.
c  4091
      xk(340)=1.30d-14*(T_K/300.d0)**(-1.00d0)
C 341)   CH3+  +   OH    ->   H2CO+ +   H2
c  1652
      xk(341)=7.20d-10*(T_K/300.d0)**(-0.50d0)
C 342)   CH3+  +   HCO   ->   CH4+  +   CO
c  1669
      xk(342)=4.40d-10*(T_K/300.d0)**(-0.50d0)
C 343)   CH3+  +   HCO   ->   HCO+  +   CH3
c  3227
      xk(343)=4.40d-10*(T_K/300.d0)**(-0.50d0)
C 344)   CH3+  +   H2CO  ->   HCO+  +   CH4
c  1674
      xk(344)=1.60d-9*(T_K/300.d0)**(-0.50d0)
C 345)   CH3+  +   O2    ->   H3CO+ +   O
c  1676
      xk(345)=5.00d-12
C 346)   O+    +   H     ->   H+    +   O
c  2944
      xk(346)=5.66d-10*(T_K/300.d0)**0.36d0*dexp(8.6d0/T_K)
C 347)   O+    +   H-    ->   H     +   O
c  3495
      xk(347)=2.30d-7*(T_K/300.d0)**(-0.50d0)
C 348)   O+    +   H2    ->   OH+   +   H
c  656
      xk(348)=1.70d-9
C 349)   O+    +   CH    ->   CH+   +   O
c  3162
      xk(349)=3.50d-10*(T_K/300.d0)**(-0.50d0)
C 350)   O+    +   CH    ->   CO+   +   H
c  1405
      xk(350)=3.50d-10*(T_K/300.d0)**(-0.50d0)
C 351)   O+    +   CH2   ->   CH2+  +   O
c  3203
      xk(351)=9.70d-10
C 352)   O+    +   CH4   ->   CH3+  +   OH
c  1711
      xk(352)=1.10d-10
C 353)   O+    +   CH4   ->   CH4+  +   O
c  3232
      xk(353)=8.90d-10
C 354)   O+    +   OH    ->   OH+   +   O
c  3233
      xk(354)=3.60d-10*(T_K/300.d0)**(-0.50d0)
C 355)   O+    +   OH    ->   O2+   +   H
c  1714
      xk(355)=3.60d-10*(T_K/300.d0)**(-0.50d0)
C 356)   O+    +   H2O   ->   H2O+  +   O
c  3235
      xk(356)=3.20d-9*(T_K/300.d0)**(-0.50d0)
C 357)   O+    +   C2    ->   C2+   +   O
c  3236
      xk(357)=4.80d-10
C 358)   O+    +   C2    ->   CO+   +   C
c  1720
      xk(358)=4.80d-10
C 359)   O+    +   HCO   ->   OH+   +   CO
c  1738
      xk(359)=4.30d-10*(T_K/300.d0)**(-0.50d0)
C 360)   O+    +   HCO   ->   HCO+  +   O
c  3245
      xk(360)=4.30d-10*(T_K/300.d0)**(-0.50d0)
C 361)   O+    +   H2CO  ->   HCO+  +   OH
c  1741
      xk(361)=1.40d-9*(T_K/300.d0)**(-0.50d0)
C 362)   O+    +   H2CO  ->   H2CO+ +   O
c  3246
      xk(362)=2.10d-9*(T_K/300.d0)**(-0.50d0)
C 363)   O+    +   O2    ->   O2+   +   O
c  3247
      xk(363)=1.90d-11
C 364)   O+    +   CO2   ->   O2+   +   CO
c  1768
      xk(364)=9.40d-10
C 365)   CH4+  +   H     ->   CH3+  +   H2
c  556
      xk(365)=1.00d-11
C 366)   CH4+  +   O     ->   CH3+  +   OH
c  1713
      xk(366)=1.00d-9
C 367)   CH4+  +   H2    ->   CH5+  +   H
c  660, 661
      if(T_K > 300.d0) then
         xk(367)=3.30d-11
      else
         xk(367)=3.40d-11*(T_K/300.d0)**(-1.35d0)*dexp(-23.d0/T_K)
      endif
C 368)   CH4+  +   CH4   ->   CH5+  +   CH3
c  1840
      xk(368)=1.50d-9
C 369)   CH4+  +   H2O   ->   H3O+  +   CH3
c  1845
      xk(369)=2.60d-9*(T_K/300.d0)**(-0.50d0)
C 370)   CH4+  +   CO    ->   HCO+  +   CH3
c  1863
      xk(370)=1.40d-9
C 371)   CH4+  +   H2CO  ->   H2CO+ +   CH4
c  3272
      xk(371)=1.62d-9*(T_K/300.d0)**(-0.50d0)
C 372)   CH4+  +   H2CO  ->   H3CO+ +   CH3
c  1869
      xk(372)=1.98d-9*(T_K/300.d0)**(-0.50d0)
C 373)   CH4+  +   O2    ->   O2+   +   CH4
c  3273
      xk(373)=3.90d-10
C 374)   CH4+  +   CO2   ->   HCO2+ +   CH3
c  1889
      xk(374)=1.20d-9
C 375)   OH+   +   C     ->   CH+   +   O
c  1185
      xk(375)=1.20d-9
C 376)   OH+   +   O     ->   O2+   +   H
c  1715
      xk(376)=7.10d-10
C 377)   OH+   +   H2    ->   H2O+  +   H
c  664
      xk(377)=1.01d-9
C 378)   OH+   +   CH    ->   CH+   +   OH
c  3164
      xk(378)=3.50d-10*(T_K/300.d0)**(-0.50d0)
C 379)   OH+   +   CH    ->   CH2+  +   O
c  1411
      xk(379)=3.50d-10*(T_K/300.d0)**(-0.50d0)
C 380)   OH+   +   CH2   ->   CH2+  +   OH
c  3205
      xk(380)=4.80d-10
C 381)   OH+   +   CH2   ->   CH3+  +   O
c  1558
      xk(381)=4.80d-10
C 382)   OH+   +   CH4   ->   H3O+  +   CH2
c  1842
      xk(382)=1.31d-9
C 383)   OH+   +   CH4   ->   CH5+  +   O
c  1841
      xk(383)=1.95d-10
C 384)   OH+   +   OH    ->   H2O+  +   O
c  1923
      xk(384)=7.00d-10*(T_K/300.d0)**(-0.50d0)
C 385)   OH+   +   H2O   ->   H2O+  +   OH
c  3279
      xk(385)=1.59d-9*(T_K/300.d0)**(-0.50d0)
C 386)   OH+   +   H2O   ->   H3O+  +   O
c  1927
      xk(386)=1.30d-9*(T_K/300.d0)**(-0.50d0)
C 387)   OH+   +   C2    ->   C2+   +   OH
c  3280
      xk(387)=4.80d-10
C 388)   OH+   +   CO    ->   HCO+  +   O
c  1937
      xk(388)=1.05d-9
C 389)   OH+   +   HCO   ->   H2O+  +   CO
c  1943
      xk(389)=2.80d-10*(T_K/300.d0)**(-0.50d0)
C 390)   OH+   +   HCO   ->   HCO+  +   OH
c  3287
      xk(390)=2.80d-10*(T_K/300.d0)**(-0.50d0)
C 391)   OH+   +   HCO   ->   H2CO+ +   O
c  1944
      xk(391)=2.80d-10*(T_K/300.d0)**(-0.50d0)
C 392)   OH+   +   H2CO  ->   H3CO+ +   O
c  1949
      xk(392)=1.12d-9*(T_K/300.d0)**(-0.50d0)
C 393)   OH+   +   H2CO  ->   H2CO+ +   OH
c  3289
      xk(393)=7.44d-10*(T_K/300.d0)**(-0.50d0)
C 394)   OH+   +   O2    ->   O2+   +   OH
c  3291
      xk(394)=5.90d-10
C 395)   OH+   +   CO2   ->   HCO2+ +   O
c  1962
      xk(395)=1.44d-9
C 396)   CH5+  +   H     ->   CH4+  +   H2
c  557
      xk(396)=1.50d-10
C 397)   CH5+  +   C     ->   CH+   +   CH4
c  1189
      xk(397)=1.20d-9
C 398)   CH5+  +   O     ->   H3O+  +   CH2
c  1717
      xk(398)=2.20d-10
C 399)   CH5+  +   O     ->   H3CO+ +   H2
c  1718
      xk(399)=4.40d-12
C 400)   CH5+  +   CH    ->   CH2+  +   CH4
c  1416
      xk(400)=6.90d-10*(T_K/300.d0)**(-0.50d0)
C 401)   CH5+  +   CH2   ->   CH3+  +   CH4
c  1562
      xk(401)=9.60d-10
C 402)   CH5+  +   OH    ->   H2O+  +   CH4
c  1926  
      xk(402)=7.00d-10*(T_K/300.d0)**(-0.50d0)
C 403)   CH5+  +   H2O   ->   H3O+  +   CH4
c  2021
      xk(403)=3.70d-9*(T_K/300.d0)**(-0.50d0)
C 404)   CH5+  +   CO    ->   HCO+  +   CH4
c  2028
      xk(404)=1.00d-9
C 405)   CH5+  +   HCO   ->   H2CO+ +   CH4
c  2032
      xk(405)=8.50d-10*(T_K/300.d0)**(-0.50d0)
C 406)   CH5+  +   H2CO  ->   H3CO+ +   CH4
c  2033
      xk(406)=4.50d-9*(T_K/300.d0)**(-0.50d0)
C 407)   CH5+  +   CO2   ->   HCO2+ +   CH4
c  2039
      xk(407)=3.20d-11
C 408)   H2O+  +   C     ->   CH+   +   OH
c  1190
      xk(408)=1.10d-9
C 409)   H2O+  +   O     ->   O2+   +   H2
c  1719
      xk(409)=4.00d-11
C 410)   H2O+  +   H2    ->   H3O+  +   H
c  669
      xk(410)=6.40d-10
C 411)   H2O+  +   CH    ->   CH+   +   H2O
c  3166
      xk(411)=3.40d-10*(T_K/300.d0)**(-0.50d0)
C 412)   H2O+  +   CH    ->   CH2+  +   OH
c  1418
      xk(412)=3.40d-10*(T_K/300.d0)**(-0.50d0)
C 413)   H2O+  +   CH2   ->   CH3+  +   OH
c  1564
      xk(413)=4.70d-10
C 414)   H2O+  +   CH2   ->   CH2+  +   H2O
c  3206
      xk(414)=4.70d-10
C 415)   H2O+  +   CH4   ->   H3O+  +   CH3
c  1846
      xk(415)=1.40d-9
C 416)   H2O+  +   OH    ->   H3O+  +   O
c  1928
      xk(416)=6.90d-10*(T_K/300.d0)**(-0.50d0)
C 417)   H2O+  +   H2O   ->   H3O+  +   OH
c  2043
      xk(417)=2.10d-9*(T_K/300.d0)**(-0.50d0)
C 418)   H2O+  +   C2    ->   C2+   +   H2O
c  3319
      xk(418)=4.70d-10
C 419)   H2O+  +   CO    ->   HCO+  +   OH
c  2058
      xk(419)=5.00d-10
C 420)   H2O+  +   HCO   ->   H3O+  +   CO
c  2062
      xk(420)=2.80d-10*(T_K/300.d0)**(-0.50d0)
C 421)   H2O+  +   HCO   ->   HCO+  +   H2O
c  3328
      xk(421)=2.80d-10*(T_K/300.d0)**(-0.50d0)
C 422)   H2O+  +   HCO   ->   H2CO+ +   OH
c  2063
      xk(422)=2.80d-10*(T_K/300.d0)**(-0.50d0)
C 423)   H2O+  +   H2CO  ->   H2CO+ +   H2O
c  3330
      xk(423)=1.41d-9*(T_K/300.d0)**(-0.50d0)
C 424)   H2O+  +   H2CO  ->   H3CO+ +   OH
c  2068
      xk(424)=6.62d-10*(T_K/300.d0)**(-0.50d0)
C 425)   H2O+  +   O2    ->   O2+   +   H2O
c  3332
      xk(425)=4.60d-10
C 426)   H3O+  +   C     ->   HCO+  +   H2
c  1194
      xk(426)=1.00d-11
C 427)   H3O+  +   H-    ->   OH    +   H2    +   H
c  4603
      xk(427)=2.30d-7*(T_K/300.d0)**(-0.50d0)
C 428)   H3O+  +   H-    ->   H2O   +   H2
c  4604
      xk(428)=2.30d-7*(T_K/300.d0)**(-0.50d0)
C 429)   H3O+  +   CH    ->   CH2+  +   H2O
c  1421
      xk(429)=6.80d-10*(T_K/300.d0)**(-0.50d0)
C 430)   H3O+  +   CH2   ->   CH3+  +   H2O
c  1565
      xk(430)=9.40d-10
C 431)   H3O+  +   H2CO  ->   H3CO+ +   H2O
c  2122
      xk(431)=3.40d-9*(T_K/300.d0)**(-0.50d0)
C 432)   C2+   +   C     ->   C+    +   C2
c  3104
      xk(432)=1.10d-10
C 433)   C2+   +   O     ->   CO+   +   C
c  1721
      xk(433)=3.10d-10
C 434)   C2+   +   CH    ->   CH+   +   C2
c  3168
      xk(434)=3.20d-10*(T_K/300.d0)**(-0.50d0)
C 435)   C2+   +   CH2   ->   CH2+  +   C2
c  3207
      xk(435)=4.50d-10
C 436)   C2+   +   OH    ->   OH+   +   C2
c  3281  
      xk(436)=6.50d-10*(T_K/300.d0)**(-0.50d0)
C 437)   C2+   +   HCO   ->   HCO+  +   C2
c  3356
      xk(437)=3.80d-10*(T_K/300.d0)**(-0.50d0)
C 438)   C2+   +   O2    ->   CO+   +   CO
c  2187
      xk(438)=8.00d-10
C 439)   CO+   +   H     ->   H+    +   CO
c  2960
      xk(439)=7.50d-10
C 440)   CO+   +   C     ->   C+    +   CO
c  3107
      xk(440)=1.10d-10
C 441)   CO+   +   O     ->   O+    +   CO
c  3241
      xk(441)=1.40d-10
C 442)   CO+   +   H2    ->   HCO+  +   H
c  683
      xk(442)=7.50d-10
C 443)   CO+   +   CH    ->   CH+   +   CO
c  3171
      xk(443)=3.20d-10*(T_K/300.d0)**(-0.50d0)
C 444)   CO+   +   CH    ->   HCO+  +   C
c  1436
      xk(444)=3.20d-10*(T_K/300.d0)**(-0.50d0)
C 445)   CO+   +   CH2   ->   CH2+  +   CO
c  3209
      xk(445)=4.30d-10
C 446)   CO+   +   CH2   ->   HCO+  +   CH
c  1573
      xk(446)=4.30d-10
C 447)   CO+   +   CH4   ->   CH4+  +   CO
c  3270
      xk(447)=7.93d-10
C 448)   CO+   +   CH4   ->   HCO+  +   CH3
c  1864
      xk(448)=4.55d-10
C 449)   CO+   +   OH    ->   OH+   +   CO
c  3285
      xk(449)=3.10d-10*(T_K/300.d0)**(-0.50d0)
C 450)   CO+   +   OH    ->   HCO+  +   O
c  1938
      xk(450)=3.10d-10*(T_K/300.d0)**(-0.50d0)
C 451)   CO+   +   H2O   ->   H2O+  +   CO
c  3324
      xk(451)=1.72d-9*(T_K/300.d0)**(-0.50d0)
C 452)   CO+   +   H2O   ->   HCO+  +   OH
c  2059
      xk(452)=8.84d-10*(T_K/300.d0)**(-0.50d0)
C 453)   CO+   +   C2    ->   C2+   +   CO
c  3354
      xk(453)=8.40d-10
C 454)   CO+   +   HCO   ->   HCO+  +   CO
c  3410
      xk(454)=7.40d-10*(T_K/300.d0)**(-0.50d0)
C 455)   CO+   +   H2CO  ->   HCO+  +   HCO
c  2450
      xk(455)=1.65d-9*(T_K/300.d0)**(-0.50d0)
C 456)   CO+   +   H2CO  ->   H2CO+ +   CO
c  3412
      xk(456)=1.35d-9*(T_K/300.d0)**(-0.50d0)
C 457)   CO+   +   O2    ->   O2+   +   CO
c  3413
      xk(457)=1.20d-10
C 458)   HCO+  +   C     ->   CH+   +   CO
c  1213
      xk(458)=1.10d-9
C 459)   HCO+  +   H-    ->   CO    +   H2
c  4605
      xk(459)=2.30d-7*(T_K/300.d0)**(-0.50d0)
C 460)   HCO+  +   CH    ->   CH2+  +   CO
c  1441
      xk(460)=6.30d-10*(T_K/300.d0)**(-0.50d0)
C 461)   HCO+  +   CH2   ->   CH3+  +   CO
c  1578
      xk(461)=8.60d-10
C 462)   HCO+  +   OH    ->   H2O+  +   CO
c  1945
      xk(462)=6.20d-10*(T_K/300.d0)**(-0.50d0)
C 463)   HCO+  +   OH    ->   HCO2+ +   H
c  1942
      xk(463)=1.00d-9*(T_K/300.d0)**(-0.50d0)
C 464)   HCO+  +   H2O   ->   H3O+  +   CO
c  2064
      xk(464)=2.50d-9*(T_K/300.d0)**(-0.50d0)
C 465)   HCO+  +   HCO   ->   H2CO+ +   CO
c  2587
      xk(465)=7.30d-10*(T_K/300.d0)**(-0.50d0)
C 466)   HCO+  +   H2CO  ->   H3CO+ +   CO  
c  2592
      xk(466)=3.30d-9*(T_K/300.d0)**(-0.50d0)
C 467)   H2CO+ +   CH    ->   CH+   +   H2CO
c  3176
      xk(467)=3.10d-10*(T_K/300.d0)**(-0.50d0)
C 468)   H2CO+ +   CH    ->   CH2+  +   HCO
c  1448
      xk(468)=3.10d-10*(T_K/300.d0)**(-0.50d0)
C 469)   H2CO+ +   CH2   ->   CH3+  +   HCO
c  1582
      xk(469)=4.30d-10
C 470)   H2CO+ +   CH2   ->   CH2+  +   H2CO
c  3212
      xk(470)=4.30d-10
C 471)   H2CO+ +   CH4   ->   H3CO+ +   CH3
c  1870
      xk(471)=9.35d-11
C 472)   H2CO+ +   H2O   ->   H3O+  +   HCO
c  2069
      xk(472)=2.60d-9*(T_K/300.d0)**(-0.50d0)
C 473)   H2CO+ +   HCO   ->   HCO+  +   H2CO
c  3442
      xk(473)=3.60d-10*(T_K/300.d0)**(-0.50d0)
C 474)   H2CO+ +   HCO   ->   H3CO+ +   CO
c  2593
      xk(474)=3.60d-10*(T_K/300.d0)**(-0.50d0)
C 475)   H2CO+ +   H2CO  ->   H3CO+ +   HCO
c  2712
      xk(475)=3.20d-9*(T_K/300.d0)**(-0.50d0)
C 476)   H2CO+ +   O2    ->   HCO+  +   O2H
c  2715
      xk(476)=7.70d-11
C 477)   H3CO+ +   CH    ->   CH2+  +   H2CO
c  1450
      xk(477)=6.20d-10*(T_K/300.d0)**(-0.50d0)
C 478)   H3CO+ +   H2O   ->   H3O+  +   H2CO
c  2076  
      xk(478)=2.30d-10*(T_K/300.d0)**(-0.50d0)
C 479)   O2+   +   C     ->   CO+   +   O
c  1237
      xk(479)=5.20d-11
C 480)   O2+   +   C     ->   C+    +   O2
c  3119
      xk(480)=5.20d-11
C 481)   O2+   +   CH    ->   CH+   +   O2
c  3177
      xk(481)=3.10d-10*(T_K/300.d0)**(-0.50d0)
C 482)   O2+   +   CH    ->   HCO+  +   O
c  1454
      xk(482)=3.10d-10*(T_K/300.d0)**(-0.50d0)
C 483)   O2+   +   CH2   ->   H2CO+ +   O
c  1585
      xk(483)=4.30d-10
C 484)   O2+   +   CH2   ->   CH2+  +   O2
c  3213
      xk(484)=4.30d-10
C 485)   O2+   +   C2    ->   CO+   +   CO
c  2188
      xk(485)=4.10d-10
C 486)   O2+   +   C2    ->   C2+   +   O2
c  3358  
      xk(486)=4.10d-10
C 487)   O2+   +   HCO   ->   O2H+  +   CO
c  2599
      xk(487)=3.60d-10*(T_K/300.d0)**(-0.50d0)
C 488)   O2+   +   HCO   ->   HCO+  +   O2
c  3443
      xk(488)=3.60d-10*(T_K/300.d0)**(-0.50d0)
C 489)   O2+   +   H2CO  ->   HCO+  +   O2    +   H
c  2714
      xk(489)=2.30d-10*(T_K/300.d0)**(-0.50d0)
C 490)   O2+   +   H2CO  ->   H2CO+ +   O2
c  3461
      xk(490)=2.07d-9*(T_K/300.d0)**(-0.50d0)
C 491)   O2H+  +   C     ->   CH+   +   O2
c  1243
      xk(491)=1.00d-9
C 492)   O2H+  +   O     ->   OH+   +   O2
c  1749
      xk(492)=6.20d-10
C 493)   O2H+  +   H2    ->   H3+   +   O2
c  697
      xk(493)=6.40d-10
C 494)   O2H+  +   CH    ->   CH2+  +   O2
c  1461
      xk(494)=6.20d-10*(T_K/300.d0)**(-0.50d0)
C 495)   O2H+  +   CH2   ->   CH3+  +   O2
c  1589
      xk(495)=8.50d-10
C 496)   O2H+  +   OH    ->   H2O+  +   O2
c  1959
      xk(496)=6.10d-10*(T_K/300.d0)**(-0.50d0)
C 497)   O2H+  +   H2O   ->   H3O+  +   O2
c  2085
      xk(497)=8.20d-10*(T_K/300.d0)**(-0.50d0)
C 498)   O2H+  +   CO    ->   HCO+  +   O2
c  2453
      xk(498)=8.40d-10
C 499)   O2H+  +   HCO   ->   H2CO+ +   O2
c  2605
      xk(499)=7.10d-10*(T_K/300.d0)**(-0.50d0)
C 500)   O2H+  +   H2CO  ->   H3CO+ +   O2
c  2720
      xk(500)=9.80d-10*(T_K/300.d0)**(-0.50d0)
C 501)   O2H+  +   CO2   ->   HCO2+ +   O2
c  2788
      xk(501)=1.10d-9
C 502)   HCO2+ +   C     ->   CH+   +   CO2
c  1296
      xk(502)=1.00d-9
C 503)   HCO2+ +   O     ->   HCO+  +   O2
c  1773
      xk(503)=1.00d-9
C 504)   HCO2+ +   CH4   ->   CH5+  +   CO2
c  1892
      xk(504)=7.80d-10
C 505)   HCO2+ +   H2O   ->   H3O+  +   CO2
c  2102
      xk(505)=2.30d-9*(T_K/300.d0)**(-0.50d0)
C 506)   HCO2+ +   CO    ->   HCO+  +   CO2
c  2459
      xk(506)=7.80d-10
C 507)   C+    +   e     ->   C     +   ph.
c  4012, 4013, 4014  
      if(T_K < 7950.d0) then
         xk(507)=4.67d-12*(T_K/300.d0)**(-0.6d0)
      elseif(T_K < 21140.d0) then 
         xk(507)=1.23d-17*(T_K/300.d0)**2.49d0*dexp(21845.6/T_K)
      else
         xk(507)=9.62d-8*(T_K/300.d0)**(-1.37d0)*dexp(-115786.2/T_K)
      endif
C 508)   CH+   +   e     ->   C     +   H
c  3524
      xk(508)=1.50d-7*(T_K/300.d0)**(-0.42d0)
C 509)   CH2+  +   e     ->   C     +   H2
c  3527
      xk(509)=7.68d-8*(T_K/300.d0)**(-0.60d0)
C 510)   CH2+  +   e     ->   CH    +   H
c  3525
      xk(510)=1.60d-7*(T_K/300.d0)**(-0.60d0)
C 511)   CH3+  +   e     ->   CH2   +   H
c  3529
      xk(511)=7.75d-8*(T_K/300.d0)**(-0.50d0)
C 512)   CH3+  +   e     ->   CH    +   H2
c  3530
      xk(512)=1.95d-7*(T_K/300.d0)**(-0.50d0)
C 513)   CH3+  +   e     ->   CH    + 2 H
c  3531
      xk(513)=2.00d-7*(T_K/300.d0)**(-0.40d0)
C 514)   CH3+  +   e     ->   CH3   +   ph.
c  4017
      xk(514)=1.10d-10*(T_K/300.d0)**(-0.50d0)
C 515)   O+    +   e     ->   O     +   ph.
c  4019
      xk(515)=3.24d-12*(T_K/300.d0)**(-0.66d0)
C 516)   CH4+  +   e     ->   CH3   +   H
c  3538
      xk(516)=1.75d-7*(T_K/300.d0)**(-0.50d0)
C 517)   CH4+  +   e     ->   CH2   + 2 H
c  3539
      xk(517)=1.75d-7*(T_K/300.d0)**(-0.50d0)
C 518)   OH+   +   e     ->   O     +   H
c  3540
      xk(518)=3.75d-8*(T_K/300.d0)**(-0.50d0)
C 519)   CH5+  +   e     ->   CH4   +   H
c  3544
      xk(519)=1.40d-8*(T_K/300.d0)**(-0.52d0)
C 520)   CH5+  +   e     ->   CH3   +   H2
c  3543
      xk(520)=1.40d-8*(T_K/300.d0)**(-0.52d0)
C 521)   H2O+  +   e     ->   OH    +   H
c  3550
      xk(521)=8.60d-8*(T_K/300.d0)**(-0.50d0)
C 522)   H2O+  +   e     ->   O     +   H2
c  3549
      xk(522)=3.90d-8*(T_K/300.d0)**(-0.50d0)
C 523)   H3O+  +   e     ->   H2O   +   H
c  3563
      xk(523)=1.08d-7*(T_K/300.d0)**(-0.50d0)
C 524)   H3O+  +   e     ->   OH    + 2 H
c  3561
      xk(524)=2.58d-7*(T_K/300.d0)**(-0.50d0)
C 525)   C2+   +   e     -> 2 C
c  3567
      xk(525)=3.00d-7*(T_K/300.d0)**(-0.50d0)
C 526)   CO+   +   e     ->   O     +   C
c  3587
      xk(526)=2.00d-7*(T_K/300.d0)**(-0.48d0)
C 527)   HCO+  +   e     ->   CO    +   H
c  3601
      xk(527)=2.40d-7*(T_K/300.d0)**(-0.69d0)
C 528)   H2CO+ +   e     ->   HCO   +   H
c  3613
      xk(528)=1.00d-7*(T_K/300.d0)**(-0.50d0)
C 529)   H2CO+ +   e     ->   CO    + 2 H
c  3612
      xk(529)=5.00d-7*(T_K/300.d0)**(-0.50d0)
C 530)   H2CO+ +   e     ->   H2CO  +   ph.
c  4023
      xk(530)=1.10d-10*(T_K/300.d0)**(-0.70d0)
C 531)   H3CO+ +   e     ->   CO    +   H     +   H2
c  3624
      xk(531)=2.00d-7*(T_K/300.d0)**(-0.50d0)
C 532)   H3CO+ +   e     ->   HCO   + 2 H
c  3625
      xk(532)=2.00d-7*(T_K/300.d0)**(-0.50d0)
C 533)   H3CO+ +   e     ->   H2CO  +   H
c  3626
      xk(533)=2.00d-7*(T_K/300.d0)**(-0.50d0)
C 534)   O2+   +   e     -> 2 O
c  3633
      xk(534)=1.95d-7*(T_K/300.d0)**(-0.70d0)
C 535)   O2H+  +   e     ->   O2    +   H
c  3641
      xk(535)=3.00d-7*(T_K/300.d0)**(-0.50d0)
C 536)   HCO2+ +   e     ->   CO2   +   H
c  3731
      xk(536)=6.00d-8*(T_K/300.d0)**(-0.64d0)
C 537)   HCO2+ +   e     ->   CO    +   OH
c  3730
      xk(537)=3.20d-7*(T_K/300.d0)**(-0.64d0)
C 538)   H2    +   C     ->   CH2   +   ph.
c  4088
      xk(538)=1.00d-17
C 539)   OH    +   CH3   ->   CH4   +   O
c  289
      xk(539)=3.27d-14*(T_K/300.d0)**2.20d0*dexp(-2240.d0/T_K)
C 540)   H+    +   H2CO  ->   CO+   +   H2   +   H
c  574
      xk(540)=1.06d-9*(T_K/300.d0)**(-0.5d0)
C 541)   He+   +   C     ->   C+    +   He
c  3082
      xk(541)=6.30d-15*(T_K/300.d0)**0.75d0
C 542)   He+   +   H2CO  ->   H2CO+ +   He
c  3093
      xk(542)=9.69d-10*(T_K/300.d0)**(-0.5d0)
C 543)   He+   +   H2CO  ->   CH2+  +   O    +   He
c  964
      xk(543)=1.71d-9*(T_K/300.d0)**(-0.5d0)
C 544)   H     +   CR    ->   H+    +   e
c  4385
c      xk(544)=4.60d-1*zeta
      xk(544)=5.98d-18*(zeta/1.36d-17)
C 545)   He    +   CR    ->   He+   +   e
c  4390
c      xk(545)=5.00d-1*zeta
      xk(545)=6.50d-18*(zeta/1.36d-17)
C 546)   C     +   CR    ->   C+    +   e
c  4391
c      xk(546)=1.77d0*zeta
      xk(546)=2.30d-17*(zeta/1.36d-17)
C 547)   O     +   CR    ->   O+    +   e
c  4393
c      xk(547)=2.62d0*zeta
      xk(547)=3.40d-17*(zeta/1.36d-17)
C 548)   H2    +   CR    ->   H+    +   H     +   e
c  4386
c      xk(548)=1.69d-2*zeta
      xk(548)=2.20d-19*(zeta/1.36d-17)
C 549)   H2    +   CR    ->   H2+   +   e
c  4389
c      xk(549)=9.23d-1*zeta
      xk(549)=1.20d-17*(zeta/1.36d-17)
C 550)   H2    +   CR    -> 2 H
c  4387
c      xk(550)=1.00d-1*zeta
       xk(550)=1.30d-18*(zeta/1.36d-17)
C 551)   H2    +   CR    ->   H+    +   H-
c  4388
c      xk(551)=3.00d-4*zeta
       xk(551)=3.90d-21*(zeta/1.36d-17)
C 552)   CO    +   CR    ->   CO+   +   e
c  4394
c      xk(552)=3.00d0*zeta
       xk(552)=3.90d-17*(zeta/1.36d-17)
C 553)   C     +   ph.   ->   C+    +   e
c  4173
      xk(553)=3.00d-10*dexp(-3.d0*A_v)*(G_0/1.71d0)
C 554)   H-    +   ph.   ->   H     +   e
c  4169
      xk(554)=2.40d-7*dexp(-0.5d0*A_v)*(G_0/1.71d0)
C 555)   H2+   +   ph.   ->   H+    +   H
c  4170
      xk(555)=5.70d-10*dexp(-1.9d0*A_v)*(G_0/1.71d0)
C 556)   H3+   +   ph.   ->   H2+   +   H
c  4171
      xk(556)=5.00d-15*dexp(-2.3d0*A_v)*(G_0/1.71d0)
C 557)   H3+   +   ph.   ->   H+    +   H2
c  4172
      xk(557)=5.00d-15*dexp(-1.8d0*A_v)*(G_0/1.71d0)
C 558)   CH    +   ph.   ->   CH+   +   e
c  4175
      xk(558)=7.60d-10*dexp(-2.8d0*A_v)*(G_0/1.71d0)
C 559)   CH    +   ph.   ->   C     +   H
c  4176
      xk(559)=8.60d-10*dexp(-1.2d0*A_v)*(G_0/1.71d0)
C 560)   CH+   +   ph.   ->   C+    +   H
c  4177
      xk(560)=2.50d-10*dexp(-2.5d0*A_v)*(G_0/1.71d0)
C 561)   CH2   +   ph.   ->   CH2+  +   e
c  4178
      xk(561)=1.00d-9*dexp(-2.3d0*A_v)*(G_0/1.71d0)
C 562)   CH2   +   ph.   ->   CH    +   H
c  4179
      xk(562)=7.20d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 563)   CH2+  +   ph.   ->   CH+   +   H
c  4180
      xk(563)=4.60d-11*dexp(-1.8d0*A_v)*(G_0/1.71d0)
C 564)   CH3   +   ph.   ->   CH    +   H2
c  4188
      xk(564)=2.50d-10*dexp(-1.9d0*A_v)*(G_0/1.71d0)
C 565)   CH3   +   ph.   ->   CH2   +   H
c  4186
      xk(565)=2.50d-10*dexp(-1.9d0*A_v)*(G_0/1.71d0)
C 566)   CH3   +   ph.   ->   CH3+  +   e
c  4185
      xk(566)=1.00d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
C 567)   CH3+  +   ph.   ->   CH2+  +   H
c  4187
      xk(567)=1.00d-9*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 568)   CH3+  +   ph.   ->   CH+   +   H2
c  4189
      xk(568)=1.00d-9*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 569)   CH4   +   ph.   ->   CH3   +   H
c  4193
      xk(569)=2.20d-10*dexp(-2.2d0*A_v)*(G_0/1.71d0)
C 570)   CH4   +   ph.   ->   CH2   +   H2
c  4194
      xk(570)=9.80d-10*dexp(-2.2d0*A_v)*(G_0/1.71d0)
C 571)   CH4   +   ph.   ->   CH    +   H     +   H2
c  4195
      xk(571)=2.20d-10*dexp(-2.2d0*A_v)*(G_0/1.71d0)
C 572)   OH    +   ph.   ->   O     +   H
c  4200
      xk(572)=3.50d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 573)   OH    +   ph.   ->   OH+   +   e
c  4198
      xk(573)=1.60d-12*dexp(-3.1d0*A_v)*(G_0/1.71d0)
C 574)   OH+   +   ph.   ->   H+    +   O
c  none
      xk(574)=0.d0
C 575)   H2O   +   ph.   ->   OH    +   H
c  4205
      xk(575)=5.90d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 576)   H2O   +   ph.   ->   H2O+  +   e
c  4206
      xk(576)=3.30d-11*dexp(-3.9d0*A_v)*(G_0/1.71d0)
C 577)   C2    +   ph.   -> 2 C
c  4210
      xk(577)=1.50d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
C 578)   C2    +   ph.   ->   C2+   +   e
c  4211
      xk(578)=4.10d-10*dexp(-3.5d0*A_v)*(G_0/1.71d0)
C 579)   C2+   +   ph.   ->   C+    +   C
c  4212
      xk(579)=1.00d-11*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 580)   CO    +   ph.   ->   C     +   O
c  4225
      shield=3.5d0*A_v+xNc_H2/1.6d21
      xk(580)=2.00d-10*dexp(-shield) *(G_0/1.71d0)
C 581)   HCO   +   ph.   ->   H     +   CO
c  4231
      xk(581)=1.10d-9*dexp(-0.8d0*A_v)*(G_0/1.71d0)
C 582)   HCO   +   ph.   ->   HCO+  +   e
c  4232
      xk(582)=5.60d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
C 583)   H2CO  +   ph.   ->   CO    +   H2
c  4243
      xk(583)=7.00d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 584)   H2CO  +   ph.   ->   CO    + 2 H
c  4242
      xk(584)=7.00d-10*dexp(-1.7d0*A_v)*(G_0/1.71d0)
C 585)   H2CO  +   ph.   ->   H2CO+ +   e
c  4241
      xk(585)=4.70d-10*dexp(-2.8d0*A_v)*(G_0/1.71d0)
C 586)   H2CO  +   ph.   ->   HCO+  +   e     +   H
c  4244
      xk(586)=1.40d-11*dexp(-3.1d0*A_v)*(G_0/1.71d0)
C 587)   O2    +   ph.   -> 2 O
c  4247
      xk(587)=6.90d-10*dexp(-1.8d0*A_v)*(G_0/1.71d0)
C 588)   O2    +   ph.   ->   O2+   +   e
c  4248
      xk(588)=5.60d-11*dexp(-3.7d0*A_v)*(G_0/1.71d0)
C 589)   CO2   +   ph.   ->   CO    +   O 
c  4284
      xk(589)=1.40d-9*dexp(-2.5d0*A_v)*(G_0/1.71d0)
C 590)   H2   +   ph.   ->   2H
c  none
      xm_p=1.67d-24
      xk_B=1.38d-16
      D5=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))*1.d-5
      x=xNc_H2/8.465d13
      fsh=0.9379d0/(1.d0+x/D5)**1.879d0
     &     +(0.03465d0/(1.d0+x)**0.473d0)
     &     *dexp(-2.293d-4*(1.d0+x)**0.5d0)
c      if(xNc_H2.le.1.d14) then
c         fsh=1.d0
c      else
c         fsh=(1.d14/xNc_H2)**0.75d0
c      endif
      D_0=G_0/1.71d0
      xk(590)=7.70d-11*fsh*dexp(-2.5d0*A_v)*D_0

C 591)   HD   +   ph.   ->    H     +   D
c  none
      D5=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p))*1.d-5
      x=xNc_HD/8.465d13
      fsh=0.9379d0/(1.d0+x/D5)**1.879d0
     &     +(0.03465d0/(1.d0+x)**0.473d0)
     &     *dexp(-2.293d-4*(1.d0+x)**0.5d0)
      x=xNc_H/2.848d23      
      fH=dexp(-0.149d0*x)/(1.d0+x)**1.620d0
      x=xNc_H/2.848d23  
      fH2=dexp(-5.2d-3*x)/(1.d0+x)**0.238d0
      fsh=fsh*fH*fH2
      xk(591)=7.70d-11*fsh*dexp(-2.5d0*A_v)*D_0
c
C supplemented reactions 
C 601)   H    +   HCO   ->   O   +   CH2
c  17
      xk(601)=6.61d-11*dexp(-51598.d0/T_K)
C 602)   H    +   H2O2  ->  H2O   +   OH
c  30
      xk(602)=1.70d-11*dexp(-1800.d0/T_K)
C 603)   C    +   CH2   ->   2 CH
c  64
      xk(603)=2.69d-12*dexp(-23550.d0/T_K)
C 604)   CH   +    O2   ->   CO   +   OH
c  154
      xk(604)=2.60d-11
C 605)   H    +   HCO   ->   O   +   CH2
c  17
      xk(605)=8.00d-11
C 606)   CH2  +   O2   ->   CO2    +   2 H
c  256
      xk(606)=3.65d-11*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
C 607)   CH2  +   O2   ->   CO2    +   H2
c  257
      xk(607)=2.92d-11*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
C 608)   CH2  +   O2   ->   CO     +   H2O
c  260
      xk(608)=2.48d-10*(T_K/300.d0)**(-3.3d0)*dexp(-1443.d0/T_K)
C 609)   2 CH3         ->   CH4    +   CH2
c  283
      xk(609)=7.13d-12*dexp(-5052.d0/T_K)
C 610)   CH3  +   O    ->   CO     +   H2    +   H
c  287
      xk(610)=3.60d-11*dexp(-202.d0/T_K)
C 611)   CH3  +   OH   ->   H2CO   +   H2
c  291
      xk(611)=1.70d-12
C 612)   CH3  +   O2   ->   HCO    +   H2O
c  301
      xk(612)=1.66d-12
C 613)   O    +  H2CO  ->   CO     +   OH  +   H
c  338
      xk(613)=1.00d-10
C 614)   C2   +   O2   ->   2 CO
c  464
      xk(614)=1.50d-11*dexp(-4300.d0/T_K)
C 615)   2 HCO         ->   2 CO   +   H2
c  515
      xk(615)=3.63d-11
C 616)   HCO   +   O2   ->   CO2   +   OH
c  519
      xk(616)=7.60d-13
C 617)    O    +   OH   ->   O2    +   H
c  544
      xk(617)=3.50d-11
C
      xk(618)=0.d0
      xk(619)=0.d0
      xk(620)=0.d0
C 621)    H3+  +    O    ->  H2O+  +   H
c  755
      xk(621)=3.60d-10
C
      xk(622)=0.d0
      xk(623)=0.d0
      xk(624)=0.d0
      xk(625)=0.d0
C 626)   O+    +  CO     ->  CO+   +   O
c  3242
      xk(626)=4.90d-12*(T_K/300.d0)**0.5d0*dexp(-4580.d0/T_K)
C
      xk(627)=0.d0
      xk(628)=0.d0
      xk(629)=0.d0
      xk(630)=0.d0
      xk(631)=0.d0
C 632)   CH2+  +  e      ->  C     +  2 H
c  3526
      xk(632)=4.03d-7*(T_K/300.d0)**(-0.6d0)
C 633)   CH5+  +  e      ->  CH3   +  2 H
c  3545
      xk(633)=1.96d-7*(T_K/300.d0)**(-0.52d0)
C 634)   CH5+  +  e      ->  CH2   +    H2    +   H
c  3546
      xk(634)=4.76d-8*(T_K/300.d0)**(-0.52d0)
C 635)   CH5+  +  e      ->  CH    +  2 H2
c  3547
      xk(635)=8.40d-9*(T_K/300.d0)**(-0.52d0)
C 636)   H2O+  +  e      ->  O     +  2 H
c  3548
      xk(636)=3.05d-7*(T_K/300.d0)**(-0.5d0)
C 637)   H3O+  +  e      ->  O     +    H2   +    H
c  3560
      xk(637)=5.60d-9*(T_K/300.d0)**(-0.5d0)
C 638)   H3O+  +  e      ->  OH    +    H2
c  3562
      xk(638)=6.02d-8*(T_K/300.d0)**(-0.5d0)
C 
      xk(639)=0.d0
C 640)   HCO2+ +  e      ->  CO    +    O   +    H
c  3732
      xk(640)=8.10d-7*(T_K/300.d0)**(-0.64d0)
C 641)   H+    +  He     ->  HeH+  +    ph.
c  4080
      xk(641)=5.26d-20*(T_K/300.d0)**(-0.51d0)
C 642)   H     +  OH     ->  H2O   +    ph.
c  4084
      xk(642)=5.26d-18*(T_K/300.d0)**(-5.22d0)*dexp(-90.d0/T_K)
C 643)   H2    +  CH     ->  CH3   +    ph.
c  4090
      xk(643)=5.09d-18*(T_K/300.d0)**(-0.71d0)*dexp(-11.6d0/T_K)
C 644)   C+    +  C      ->  C2+   +    ph.
c  4107
      xk(644)=4.01d-18*(T_K/300.d0)**0.17d0*dexp(-101.5d0/T_K)
C 645)   C     +  O+     ->  CO+   +    ph.
c  4111
      if(T_K > 2000.d0) then
         xk(645)=4.69d-11*(T_K/300.d0)**(-3.08d0)
     &        *dexp(2114.d0/T_K)
      else
         xk(645)=3.91d-11
      endif
C 646)   CH2+  +  ph.    ->  CH    +    H+
c  4181
      xk(646)=4.60d-11*dexp(-1.8d0*A_v)*(G_0/1.71d0)
C 647)   CH2+  +  ph.    ->  C+    +    H2
c  4182
      xk(647)=4.60d-11*dexp(-1.8d0*A_v)*(G_0/1.71d0)
C 648)   CH4+  +  ph.    ->  CH2+  +    H2
c  4196
      xk(648)=3.40d-11*dexp(-2.d0*A_v)*(G_0/1.71d0)
C 649)   CH4+  +  ph.    ->  CH3+  +    H
c  4197
      xk(649)=8.00d-12*dexp(-2.d0*A_v)*(G_0/1.71d0)
C 650)   OH+   +  ph.    ->  O+    +    H
c  4201
      xk(650)=1.00d-12*dexp(-1.8d0*A_v)*(G_0/1.71d0)
C 651)   H2O+   +  ph.   ->  OH+   +    H
c  4207
      xk(651)=1.00d-12*dexp(-2.d0*A_v)*(G_0/1.71d0)
C 652)   CO+   +  ph.   ->   C+    +    O
c  4226
      xk(652)=1.30d-10*dexp(-2.1d0*A_v)*(G_0/1.71d0)
C 653)   HCO+  +  ph.   ->   CO+   +    H
c  4233
      xk(653)=5.40d-12*dexp(-2.d0*A_v)*(G_0/1.71d0)
C 654)   O2+   +  ph.   ->   O+    +    O
c  4249
      xk(654)=2.40d-11*dexp(-2.d0*A_v)*(G_0/1.71d0)
C 655)   H2O2  +  ph.   -> 2 OH 
c  4258
      xk(655)=8.30d-10*dexp(-1.8d0*A_v)*(G_0/1.71d0)
c albedo
      omega=0.6d0
C 656)   C     + CR ph. ->   C+    +    e
c  4396
      xk(656)=zeta*255.d0/(1.d0-omega)
C 657)   CH    + CR ph. ->   C     +    H
c  4397
      xk(657)=zeta*365.d0/(1.d0-omega)
C 658)   CH+   + CR ph. ->   C+    +    H
c  4398
      xk(658)=zeta*88.d0/(1.d0-omega)
C 659)   CH2   + CR ph. ->   CH2+  +    e
c  4399
      xk(659)=zeta*250.d0/(1.d0-omega)
C 660)   CH2   + CR ph. ->   CH    +    H
c  4400
      xk(660)=zeta*250.d0/(1.d0-omega)
C 661)   CH3   + CR ph. ->   CH3+  +    e
c  4403
      xk(661)=zeta*250.d0/(1.d0-omega)
C 662)   CH3   + CR ph. ->   CH2   +    H
c  4404
      xk(662)=zeta*250.d0/(1.d0-omega)
C 663)   CH3   + CR ph. ->   CH    +    H2
c  4405
      xk(663)=zeta*250.d0/(1.d0-omega)
C 664)   CH4   + CR ph. ->   CH2   +    H2
c  4408
      xk(664)=zeta*1169.5d0/(1.d0-omega)
C 665)   OH    + CR ph. ->   O     +    H
c  4409
      xk(665)=zeta*254.5d0/(1.d0-omega)
C 666)   H2O   + CR ph. ->   OH    +    H
c  4413
      xk(666)=zeta*485.5d0/(1.d0-omega)
C 667)   C2    + CR ph. -> 2 C
c  4416
      xk(667)=zeta*119.5d0/(1.d0-omega)
C 668)   CO    + CR ph. ->   O     +    C
c  4427
      xk(668)=zeta*(T_K/300.d0)**1.17d0*105.d0/(1.d0-omega)
C 669)  HCO    + CR ph. ->   CO    +    H
c  4432
      xk(669)=zeta*210.5d0/(1.d0-omega)
C 670)  HCO    + CR ph. ->   HCO+  +    e
c  4433
      xk(670)=zeta*584.5d0/(1.d0-omega)
C 671)  H2CO   + CR ph. ->   CO    +    H2
c  4440
      xk(671)=zeta*1329.5d0/(1.d0-omega)
C 672)  O2     + CR ph. -> 2 O 
c  4445
      xk(672)=zeta*375.5d0/(1.d0-omega)
C 673)  O2     + CR ph. ->   O2+   +    e
c  4446
      xk(673)=zeta*58.5d0/(1.d0-omega)
C 674)  H2O2   + CR ph. -> 2 OH
c  4452
      xk(674)=zeta*750.d0/(1.d0-omega)
C 675)  CO2    + CR ph. ->   CO    +    O 
c  4474
      xk(675)=zeta*854.d0/(1.d0-omega)

      return
      END

      SUBROUTINE react_rat(xk,xnH,y,r_f_tot)
      IMPLICIT REAL*8(a-h,o-z)
****************************************************************
*     dy(i)/dt=r_f_tot(i)                                      * 
*     This subroutine returns reaction rate for each spieces,  *
*     r_f_tot(i).                                              *
****************************************************************
C     N_sp = number of spiecies
C     N_react = number of reactions
      PARAMETER(N_sp=50,N_react=675)   
      DIMENSION  y(N_sp),r_f(N_react,N_sp),r_f_tot(N_sp),xk(N_react)
****************************************************************
*     SPECIES                                                  *
*     1 : H      2 : H2     3 : e      4 : H+     5 : H2+      *
*     6 : H3+    7 : H-                                        *
*     8 : He     9 : He+    10: He++   11: HeH+                *
*     12: D      13: HD     14: D+     15: HD+    16: D-      *
*     17: C      18: C2     19: CH     20: CH2    21: CH3      *
*     22: CH4    23: C+     24: C2+    25: CH+    26: CH2+     *
*     27: CH3+   28: CH4+   29: CH5+                           *
*     30: O      31: O2     32: OH     33: CO     34: H2O      *
*     35: HCO    36: O2H    37: CO2    38: H2CO   39: H2O2     *
*     40: O+     41: O2+    42: OH+    43: CO+    44: H2O+     *   
*     45: HCO+   46: O2H+   47: H3O+   48: H2CO+  49: HCO2+    *
*     50: H3CO+                                                *
****************************************************************

      do isp=1,N_sp
         do ire=1,N_react
            r_f(ire,isp)=0.d0
         enddo
      enddo

********* primordial gas reactions *********

C   1)   H     +   e     ->   H+    + 2 e 
C         (1         3          4       2*3) 
      rate=xk(1)*y(1)*y(3)*xnH
      r_f(1,1)=-rate
      r_f(1,3)=rate
      r_f(1,4)=rate
C   2)   H+    +   e     ->   H     +   ph.
C         (4         3          1)
      rate=xk(2)*y(4)*y(3)*xnH
      r_f(2,4)=-rate
      r_f(2,3)=-rate
      r_f(2,1)=rate
C   3)   He    +   e     ->   He+   + 2 e   
C         (8         3          9       2*3)
      rate=xk(3)*y(8)*y(3)*xnH
      r_f(3,8)=-rate
      r_f(3,9)=rate
      r_f(3,3)=rate
C   4)   He+   +   e     ->   He    +   ph.     
C         (9         3          8)
      rate=xk(4)*y(9)*y(3)*xnH
      r_f(4,9)=-rate
      r_f(4,3)=-rate
      r_f(4,8)=rate
C   5)   He+   +   e     ->   He++  + 2 e   
C         (9         3          10      2*3)
      rate=xk(5)*y(9)*y(3)*xnH
      r_f(5,9)=-rate
      r_f(5,3)=rate
      r_f(5,10)=rate
C   6)   He++  +   e     ->   He+   +   ph. 
C         (10        3          9)
      rate=xk(6)*y(10)*y(3)*xnH
      r_f(6,10)=-rate
      r_f(6,3)=-rate
      r_f(6,9)=rate
C   7)   H     +   e     ->   H-    +   ph. 
C         (1         3          7)
      rate=xk(7)*y(1)*y(3)*xnH
      r_f(7,1)=-rate
      r_f(7,3)=-rate
      r_f(7,7)=rate
C   8)   H-    +   H     ->   H2    +   e       
C         (7         1          2         3)
      rate=xk(8)*y(7)*y(1)*xnH
      r_f(8,7)=-rate
      r_f(8,1)=-rate
      r_f(8,2)=rate
      r_f(8,3)=rate
C   9)   H     +   H+    ->   H2+   +   ph.    
C         (1         4          5)
      rate=xk(9)*y(1)*y(4)*xnH
      r_f(9,1)=-rate
      r_f(9,4)=-rate
      r_f(9,5)=rate
C  10)   H2+   +   H     ->   H2    +   H+  
C         (5         1          2         4)
      rate=xk(10)*y(5)*y(1)*xnH
      r_f(10,5)=-rate
      r_f(10,1)=-rate
      r_f(10,2)=rate
      r_f(10,4)=rate
C  11)   H2    +   H+    ->   H2+   +   H  
C         (2         4          5         1)
      rate=xk(11)*y(2)*y(4)*xnH
      r_f(11,2)=-rate
      r_f(11,4)=-rate
      r_f(11,5)=rate
      r_f(11,1)=rate
C  12)   H2    +   e     -> 2 H     +   e        
C         (2         3        2*1         3)
      rate=xk(12)*y(2)*y(3)*xnH
      r_f(12,2)=-rate
      r_f(12,1)=2.d0*rate
C  13)   H2    +   H     -> 3 H    
C         (2         1        3*1)
      rate=xk(13)*y(2)*y(1)*xnH
      r_f(13,2)=-rate
      r_f(13,1)=2.d0*rate
C  14)   H-    +   e     ->   H     + 2 e   
C         (7         3          1       2*3)
      rate=xk(14)*y(7)*y(3)*xnH
      r_f(14,7)=-rate
      r_f(14,3)=rate
      r_f(14,1)=rate
C  15)   H-    +   H+    -> 2 H      
C         (7         4        2*1)
      rate=xk(15)*y(7)*y(4)*xnH
      r_f(15,7)=-rate
      r_f(15,4)=-rate
      r_f(15,1)=2.d0*rate
C  16)   H-    +   H+    ->   H2+   +   e          
C         (7         4          5         3)
      rate=xk(16)*y(7)*y(4)*xnH
      r_f(16,7)=-rate
      r_f(16,4)=-rate
      r_f(16,5)=rate
      r_f(16,3)=rate
C  17)   H2+   +   e     -> 2 H   
C         (5         3        2*1)
      rate=xk(17)*y(5)*y(3)*xnH
      r_f(17,5)=-rate
      r_f(17,3)=-rate
      r_f(17,1)=2.d0*rate
C  18)   H2+   +   H-    ->   H2    +   H    
C         (5         7          2         1)
      rate=xk(18)*y(5)*y(7)*xnH
      r_f(18,5)=-rate
      r_f(18,7)=-rate
      r_f(18,2)=rate
      r_f(18,1)=rate
C  19) 3 H               ->   H2    +   H  
C       (3*1                    2         1)
      rate=xk(19)*(y(1)**3)*(xnH**2)
      r_f(19,1)=-2.d0*rate
      r_f(19,2)=rate
C  20) 2 H     +   H2    -> 2 H2  
C       (2*1         2        2*2)
      rate=xk(20)*(y(1)**2)*y(2)*(xnH**2) 
      r_f(20,1)=-2.d0*rate
      r_f(20,2)=rate
C  21) 2 H2              -> 2 H     +   H2        
C       (2*2                  2*1         2)
      rate=xk(21)*(y(2)**2)*xnH
      r_f(21,2)=-rate
      r_f(21,1)=2.d0*rate
C  22) 2 H               ->   H+    +   e     +   H     
C       (2*1                    4         3         1)
      rate=xk(22)*(y(1)**2)*xnH
      r_f(22,1)=-rate
      r_f(22,4)=rate
      r_f(22,3)=rate
C  23) 2 H     +   grain ->   H2  
C       (2*1                    2) special
      rate=xk(23)*y(1)*xnH
      r_f(23,1)=-2.d0*rate
      r_f(23,2)=rate
C  24)   He+   +   H2    ->   H+    +   H     +   He                 
C         (9         2          4         1         8)
      rate=xk(24)*y(9)*y(2)*xnH
      r_f(24,9)=-rate
      r_f(24,2)=-rate
      r_f(24,4)=rate
      r_f(24,1)=rate
      r_f(24,8)=rate
C  25)   H2+   +   He    ->   HeH+  +   H                  
C         (5         8          11        1)
      rate=xk(25)*y(5)*y(8)*xnH
      r_f(25,5)=-rate
      r_f(25,8)=-rate
      r_f(25,11)=rate
      r_f(25,1)=rate
C  26)   H2+   +   H2    ->   H3+   +   H                 
C         (5         2          6         1)
      rate=xk(26)*y(5)*y(2)*xnH
      r_f(26,5)=-rate
      r_f(26,2)=-rate
      r_f(26,6)=rate
      r_f(26,1)=rate
C  27)   H3+   +   H-    -> 2 H2                    
C         (6         7        2*2)
      rate=xk(27)*y(6)*y(7)*xnH
      r_f(27,6)=-rate
      r_f(27,7)=-rate
      r_f(27,2)=2.d0*rate
C  28)   He+   +   H     ->   H+    +   He                 
C         (9         1          4         8)
      rate=xk(28)*y(9)*y(1)*xnH
      r_f(28,9)=-rate
      r_f(28,1)=-rate
      r_f(28,4)=rate
      r_f(28,8)=rate
C  29)   He+   +   H-    ->   H     +   He                    
C         (9         7          1         8)
      rate=xk(29)*y(9)*y(7)*xnH
      r_f(29,9)=-rate
      r_f(29,7)=-rate
      r_f(29,1)=rate
      r_f(29,8)=rate
C  30)   He+   +   H2    ->   H2+   +   He               
C         (9         2          5         8)
      rate=xk(30)*y(9)*y(2)*xnH
      r_f(30,9)=-rate
      r_f(30,2)=-rate
      r_f(30,5)=rate
      r_f(30,8)=rate
C  31)   HeH+  +   H     ->   H2+   +   He                 
C         (11        1          5         8)
      rate=xk(31)*y(11)*y(1)*xnH
      r_f(31,11)=-rate
      r_f(31,1)=-rate
      r_f(31,5)=rate
      r_f(31,8)=rate
C  32)   HeH+  +   H2    ->   H3+   +   He                   
C         (11        2          6         8)
      rate=xk(32)*y(11)*y(2)*xnH
      r_f(32,11)=-rate
      r_f(32,2)=-rate
      r_f(32,6)=rate
      r_f(32,8)=rate
C  33)   H2+   +   e     ->   H2    +   ph.                
C         (5         3          2)
      rate=xk(33)*y(5)*y(3)*xnH
      r_f(33,5)=-rate
      r_f(33,3)=-rate
      r_f(33,2)=rate
C  34)   H3+   +   e     ->   H2    +   H                
C         (6         3          2         1)
      rate=xk(34)*y(6)*y(3)*xnH
      r_f(34,6)=-rate
      r_f(34,3)=-rate
      r_f(34,2)=rate
      r_f(34,1)=rate
C  35)   H3+   +   e     -> 3 H                        
C         (6         3        3*1)
      rate=xk(35)*y(6)*y(3)*xnH
      r_f(35,6)=-rate
      r_f(35,3)=-rate
      r_f(35,1)=3.d0*rate
C  36)   HeH+  +   e     ->   H     +   He              
C         (11        3          1         8)
      rate=xk(36)*y(11)*y(3)*xnH
      r_f(36,11)=-rate
      r_f(36,3)=-rate
      r_f(36,1)=rate
      r_f(36,8)=rate
c  37)   H2    +   He    -> 2 H     +   He
c       (2         8        2x1         8)
      rate=xk(37)*y(2)*y(8)*xnH
      r_f(37,2)=-rate
      r_f(37,1)=2.d0*rate
c  38)   H-    +   H     -> 2 H     +   e
c       (7         1        2x1         3)
      rate=xk(38)*y(7)*y(1)*xnH
      r_f(38,7)=-rate
      r_f(38,1)=rate
      r_f(38,3)=rate
c  39)   H-    +   H2+   -> 3 H
c       (7         5        3x1     )
      rate=xk(39)*y(7)*y(5)*xnH
      r_f(39,7)=-rate
      r_f(39,5)=-rate
      r_f(39,1)=3.d0*rate
c  40)   H2    +   e     ->   H-    +   H
c       (2         3          7         1)
      rate=xk(40)*y(2)*y(3)*xnH
      r_f(40,2)=-rate
      r_f(40,3)=-rate
      r_f(40,7)=rate
      r_f(40,1)=rate
c  41)   He    +   H+    ->   He+   +   H
c       (8         4          9         1)
      rate=xk(41)*y(8)*y(4)*xnH
      r_f(41,8)=-rate
      r_f(41,4)=-rate
      r_f(41,9)=rate
      r_f(41,1)=rate
c  42)   He    +   H-    ->   He    +   H   +   e
c       (8         7          8         1       3)
      rate=xk(42)*y(8)*y(7)*xnH
      r_f(42,7)=-rate
      r_f(42,1)=rate
      r_f(42,3)=rate
c  43) 2 H     +   He    ->   H2    +   He
c     (2x1         8          2         8)
      rate=xk(43)*(y(1)**2)*y(8)*(xnH**2)
      r_f(43,1)=-2.d0*rate
      r_f(43,2)=rate
c
c     D reactions
c
C  51)   D+  +   e    ->   D    +   ph. 
c     (  14      3         12   )
      rate=xk(51)*y(14)*y(3)*xnH
      r_f(51,14)=-rate
      r_f(51,3)=-rate
      r_f(51,12)=rate
C  52)   D   +   H+   ->   D+   +   H   
c     (  12      4         14       1 )
      rate=xk(52)*y(12)*y(4)*xnH
      r_f(52,12)=-rate
      r_f(52,4)=-rate
      r_f(52,14)=rate
      r_f(52,1)=rate
C  53)   D+  +   H    ->   D    +   H+  
c     (  14      1         12       4 )
      rate=xk(53)*y(14)*y(1)*xnH
      r_f(53,14)=-rate
      r_f(53,1)=-rate
      r_f(53,12)=rate
      r_f(53,4)=rate
C  54)   D   +   H    ->   HD   +   ph. 
c     (  12      1         13  )
      rate=xk(54)*y(12)*y(1)*xnH
      r_f(54,12)=-rate
      r_f(54,1)=-rate
      r_f(54,13)=rate
C  55)   D   +   H2   ->   H    +   HD 
c     (  12      2         1        13  )
      rate=xk(55)*y(12)*y(2)*xnH
      r_f(55,12)=-rate
      r_f(55,2)=-rate
      r_f(55,1)=rate
      r_f(55,13)=rate
C  56)   HD+ +   H    ->   H+   +   HD   
c     (  15      1         4        13  )
      rate=xk(56)*y(15)*y(1)*xnH
      r_f(56,15)=-rate
      r_f(56,1)=-rate
      r_f(56,4)=rate
      r_f(56,13)=rate
C  57)   D+  +   H2   ->   H+   +   HD  
c     (  14      2         4        13  )
      rate=xk(57)*y(14)*y(2)*xnH
      r_f(57,14)=-rate
      r_f(57,2)=-rate
      r_f(57,4)=rate
      r_f(57,13)=rate
C  58)   HD  +   H    ->   H2   +   D
c     (  13      1         2        12  )
      rate=xk(58)*y(13)*y(1)*xnH
      r_f(58,13)=-rate
      r_f(58,1)=-rate
      r_f(58,2)=rate
      r_f(58,12)=rate
C  59)   HD  +   H+   ->   H2   +   D+ 
c     (  13      4         2        14  )
      rate=xk(59)*y(13)*y(4)*xnH
      r_f(59,13)=-rate
      r_f(59,4)=-rate
      r_f(59,2)=rate
      r_f(59,14)=rate
C  60)   D   +   H+   ->   HD+  +   ph.
c     (  12      4         15   )
      rate=xk(60)*y(12)*y(4)*xnH
      r_f(60,12)=-rate
      r_f(60,4)=-rate
      r_f(60,15)=rate
C  61)   D+  +   H    ->   HD+  +   ph.   
c     (  14      1         15   )
      rate=xk(61)*y(14)*y(1)*xnH
      r_f(61,14)=-rate
      r_f(61,1)=-rate
      r_f(61,15)=rate
C  62)   HD+ +   e    ->   H    +   D  
c     (  15      3         1        12  )
      rate=xk(62)*y(15)*y(3)*xnH
      r_f(62,15)=-rate
      r_f(62,3)=-rate
      r_f(62,1)=rate
      r_f(62,12)=rate
C  63)   D   +   e    ->   D-   +   ph.
c     (  12      3         16  )
      rate=xk(63)*y(12)*y(3)*xnH
      r_f(63,12)=-rate
      r_f(63,3)=-rate
      r_f(63,16)=rate
C  64)   D+  +   D-   ->   2D 
c     (  14      16        2*12)
      rate=xk(64)*y(14)*y(16)*xnH
      r_f(64,14)=-rate
      r_f(64,16)=-rate
      r_f(64,12)=2.d0*rate
C  65)   H+  +   D-   ->   D    +   H 
c     (  4       16        12       1  )
      rate=xk(65)*y(4)*y(16)*xnH
      r_f(65,4)=-rate
      r_f(65,16)=-rate
      r_f(65,12)=rate
      r_f(65,1)=rate
C  66)   H-  +   D    ->   H    +   D-   
c     (  7       12        1        16 )
      rate=xk(66)*y(7)*y(12)*xnH
      r_f(66,7)=-rate
      r_f(66,12)=-rate
      r_f(66,1)=rate
      r_f(66,16)=rate
C  67)   D-  +   H    ->   D    +   H-  
c     (  16      1         12       7  )
      rate=xk(67)*y(16)*y(1)*xnH
      r_f(67,16)=-rate
      r_f(67,1)=-rate
      r_f(67,12)=rate
      r_f(67,7)=rate
C  68)   D-  +   H    ->   HD   +   e
c     (  16      1         13       3  )
      rate=xk(68)*y(16)*y(1)*xnH
      r_f(68,16)=-rate
      r_f(68,1)=-rate
      r_f(68,13)=rate
      r_f(68,3)=rate
C  69)   D   +   e    ->   D+   +  2 e
c     ( 12       3         14      2x3    )
      rate=xk(69)*y(12)*y(3)*xnH
      r_f(69,12)=-rate
      r_f(69,14)=rate
      r_f(69,3)=rate
C  70)   He+ +   D    ->   D+   +   He
c     (  9      12         14       8   )
      rate=xk(70)*y(9)*y(12)*xnH
      r_f(70,9)=-rate
      r_f(70,12)=-rate
      r_f(70,14)=rate
      r_f(70,8)=rate
C  71)   He  +   D+   ->   D    +   He+
c     (  8       14        12       9   )
      rate=xk(71)*y(8)*y(14)*xnH
      r_f(71,8)=-rate
      r_f(71,14)=-rate
      r_f(71,12)=rate
      r_f(71,9)=rate
C  72)   H2+ +   D    ->  HD+   +   H
c     (  5      12        15        1   )
      rate=xk(72)*y(5)*y(12)*xnH
      r_f(72,5)=-rate
      r_f(72,12)=-rate
      r_f(72,15)=rate
      r_f(72,1)=rate
C  73)   HD+ +   D    ->  HD    +   D+
c     (  15     12        13        14   )
      rate=xk(73)*y(15)*y(12)*xnH
      r_f(73,15)=-rate
      r_f(73,12)=-rate
      r_f(73,13)=rate
      r_f(73,14)=rate
C  74)   HD+ +   H    ->  H2+   +   D
c     (  15      1        5        12   )
      rate=xk(74)*y(15)*y(1)*xnH
      r_f(74,15)=-rate
      r_f(74,1)=-rate
      r_f(74,5)=rate
      r_f(74,12)=rate
C  75)   HD  +   e    ->  H     +   D-
c     (  13      3        1         16  )
      rate=xk(75)*y(13)*y(3)*xnH
      r_f(75,13)=-rate
      r_f(75,3)=-rate
      r_f(75,1)=rate
      r_f(75,16)=rate
C  76)   HD  +   e    ->  D     +   H-
c     (  13      3       12         7   )
      rate=xk(76)*y(13)*y(3)*xnH
      r_f(76,13)=-rate
      r_f(76,3)=-rate
      r_f(76,12)=rate
      r_f(76,7)=rate
C  77)   H+  +   D-   ->  HD+   +   e
c     (  4       16       15        3   )
      rate=xk(77)*y(4)*y(16)*xnH
      r_f(77,4)=-rate
      r_f(77,16)=-rate
      r_f(77,15)=rate
      r_f(77,3)=rate
C  78)   D+  +   H-   ->  HD+   +   e
c     (  14      7        15        3  )
      rate=xk(78)*y(14)*y(7)*xnH
      r_f(78,14)=-rate
      r_f(78,7)=-rate
      r_f(78,15)=rate
      r_f(78,3)=rate
C  79)   D-  +   e    ->  D     + 2 e
c     (  16      3       12       2x3     )
      rate=xk(79)*y(16)*y(3)*xnH
      r_f(79,16)=-rate
      r_f(79,12)=rate
      r_f(79,3)=rate
C  80)   D-  +   H    ->  D     +   H    +   e
c     (  16      1       12         1        3 )
      rate=xk(80)*y(16)*y(1)*xnH
      r_f(80,16)=-rate
      r_f(80,12)=rate
      r_f(80,3)=rate
C  81)   D-  +   He   ->  D     +   He   +   e
c     (  16      8       12         8        3 )
      rate=xk(81)*y(16)*y(8)*xnH
      r_f(81,16)=-rate
      r_f(81,12)=rate
      r_f(81,3)=rate
C  82)   D+  +   H-   ->  D     +   H
c     (  14      7       12         1   )
      rate=xk(82)*y(14)*y(7)*xnH
      r_f(82,14)=-rate
      r_f(82,7)=-rate
      r_f(82,12)=rate
      r_f(82,1)=rate
C  83)   H2+ +   D-   ->  H2    +   D
c     (  5       16       2        12    )
      rate=xk(83)*y(5)*y(16)*xnH
      r_f(83,5)=-rate
      r_f(83,16)=-rate
      r_f(83,2)=rate
      r_f(83,12)=rate
C  84)   H2+ +   D-   -> 2 H    +   D
c     (  5       16      2x1       12   )
      rate=xk(84)*y(5)*y(16)*xnH
      r_f(84,5)=-rate
      r_f(84,16)=-rate
      r_f(84,1)=2.d0*rate
      r_f(84,12)=rate
C  85)   HD+ +   H-   ->   HD   +   H
c     (  15      7         13       1   )
      rate=xk(85)*y(15)*y(7)*xnH
      r_f(85,15)=-rate
      r_f(85,7)=-rate
      r_f(85,13)=rate
      r_f(85,1)=rate
C  86)   HD+ +   H-   ->   D    + 2 H
c     (  15      7        12      2x1     )
      rate=xk(86)*y(15)*y(7)*xnH
      r_f(86,15)=-rate
      r_f(86,7)=-rate
      r_f(86,12)=rate
      r_f(86,1)=2.d0*rate
C  87)   HD+ +   D-   ->  HD    +   D
c     (  15      16       13       12   )
      rate=xk(87)*y(15)*y(16)*xnH
      r_f(87,15)=-rate
      r_f(87,16)=-rate
      r_f(87,13)=rate
      r_f(87,12)=rate
C  88)   HD+ +   D-   ->  2 D   +   H
c     (  15      16       2x12      1    )
      rate=xk(88)*y(15)*y(16)*xnH
      r_f(88,15)=-rate
      r_f(88,16)=-rate
      r_f(88,12)=2.d0*rate
      r_f(88,1)=rate
C  89)   He+ +   D-   ->  He    +   D
c     (  9       16       8        12   )
      rate=xk(89)*y(9)*y(16)*xnH
      r_f(89,9)=-rate
      r_f(89,16)=-rate
      r_f(89,8)=rate
      r_f(89,12)=rate
C  90)   D   +   H2+  ->  H2    +   D+
c     ( 12       5        2         14  )
      rate=xk(90)*y(12)*y(5)*xnH
      r_f(90,12)=-rate
      r_f(90,5)=-rate
      r_f(90,2)=rate
      r_f(90,14)=rate
C  91)  H2+  +   D    ->  HD    +   H+
c     ( 5       12        13        4    )
      rate=xk(91)*y(5)*y(12)*xnH
      r_f(91,5)=-rate
      r_f(91,12)=-rate
      r_f(91,13)=rate
      r_f(91,4)=rate
C  92)  HD+  +   H    ->  H2    +   D+
c     ( 15       1        2         14  )
      rate=xk(92)*y(15)*y(1)*xnH
      r_f(92,15)=-rate
      r_f(92,1)=-rate
      r_f(92,2)=rate
      r_f(92,14)=rate
C  93)  H2   +   D+   ->  H2+   +   D
c     ( 2        14       5        12   )
      rate=xk(93)*y(2)*y(14)*xnH
      r_f(93,2)=-rate
      r_f(93,14)=-rate
      r_f(93,5)=rate
      r_f(93,12)=rate
C  94)  H2   +   D+   ->  HD+   +   H
c     ( 2        14       15        1   )
      rate=xk(94)*y(2)*y(14)*xnH
      r_f(94,2)=-rate
      r_f(94,14)=-rate
      r_f(94,15)=rate
      r_f(94,1)=rate
C  95)  HD   +   H+   ->  HD+   +   H
c     ( 13       4        15        1   )
      rate=xk(95)*y(13)*y(4)*xnH
      r_f(95,13)=-rate
      r_f(95,4)=-rate
      r_f(95,15)=rate
      r_f(95,1)=rate
C  96)  HD   +   H+   ->  H2+   +   D
c     ( 13       4        5        12   )
      rate=xk(96)*y(13)*y(4)*xnH
      r_f(96,13)=-rate
      r_f(96,4)=-rate
      r_f(96,5)=rate
      r_f(96,12)=rate
C  97)  HD   +   D+   ->  HD+   +   D
c     ( 13       14       15       12   )
      rate=xk(97)*y(13)*y(14)*xnH
      r_f(97,13)=-rate
      r_f(97,14)=-rate
      r_f(97,15)=rate
      r_f(97,12)=rate
C  98)  HD   +   He+  ->  HD+   +   He
c     ( 13       9        15        8   )
      rate=xk(98)*y(13)*y(9)*xnH
      r_f(98,13)=-rate
      r_f(98,9)=-rate
      r_f(98,15)=rate
      r_f(98,8)=rate
C  99)  HD   +   He+  ->  He    +   H+   +   D
c     ( 13       9        8         4       12 )
      rate=xk(99)*y(13)*y(9)*xnH
      r_f(99,13)=-rate
      r_f(99,9)=-rate
      r_f(99,8)=rate
      r_f(99,4)=rate
      r_f(99,12)=rate
C 100)  HD   +   He+  ->  He    +   H    +   D+
c     ( 13       9        8         1        14)
      rate=xk(100)*y(13)*y(9)*xnH
      r_f(100,13)=-rate
      r_f(100,9)=-rate
      r_f(100,8)=rate
      r_f(100,1)=rate
      r_f(100,14)=rate
C D50)  HD   +   H    -> 2 H    +   D       
c     ( 13       1       2x1       12 )
      rate=xk(50)*y(13)*y(1)*xnH
      r_f(50,13)=-rate
      r_f(50,1)=rate
      r_f(50,12)=rate
C D49)  HD   +   H2   ->  H     +   D    +   H2   
c     ( 13       2        1        12        2)
      rate=xk(49)*y(13)*y(2)*xnH
      r_f(49,13)=-rate
      r_f(49,1)=rate
      r_f(49,12)=rate
C D48)  HD   +   He   ->  H     +   D    +   He
c     ( 13       8        1        12        8 )
      rate=xk(48)*y(13)*y(8)*xnH
      r_f(48,13)=-rate
      r_f(48,1)=rate
      r_f(48,12)=rate
C D47)  HD   +   e    ->  H     +   D    +   e
c     ( 13       3        1        12        3 )
      rate=xk(47)*y(13)*y(3)*xnH
      r_f(47,13)=-rate
      r_f(47,1)=rate
      r_f(47,12)=rate
c
c
c     C,O reactions
C  101)   H     +   CH    ->   C     +   H2                    
C         (1         19         17        2)
      rate=xk(101)*y(1)*y(19)*xnH
      r_f(101,1)=-rate
      r_f(101,19)=-rate
      r_f(101,17)=rate
      r_f(101,2)=rate
C  102)   H     +   CH    ->   C     + 2 H                      
C         (1         19         17      2*1)
      rate=xk(102)*y(1)*y(19)*xnH
      r_f(102,19)=-rate
      r_f(102,17)=rate
      r_f(102,1)=rate
C  103)   H     +   CH2   ->   CH    +   H2                     
C         (1         20         19        2)
      rate=xk(103)*y(1)*y(20)*xnH
      r_f(103,1)=-rate
      r_f(103,20)=-rate
      r_f(103,19)=rate
      r_f(103,2)=rate
C  104)   H     +   CH3   ->   CH2   +   H2                     
C         (1         21         20        2)
      rate=xk(104)*y(1)*y(21)*xnH
      r_f(104,1)=-rate
      r_f(104,21)=-rate
      r_f(104,20)=rate
      r_f(104,2)=rate
C  105)   H     +   CH4   ->   H2    +   CH3                     
C         (1         22         2         21)
      rate=xk(105)*y(1)*y(22)*xnH
      r_f(105,1)=-rate
      r_f(105,22)=-rate
      r_f(105,2)=rate
      r_f(105,21)=rate
C  106)   H     +   OH    ->   H2    +   O                        
C         (1         32         2         30) 
      rate=xk(106)*y(1)*y(32)*xnH
      r_f(106,1)=-rate
      r_f(106,32)=-rate
      r_f(106,2)=rate
      r_f(106,30)=rate
C  107)   H     +   OH    ->   O     + 2 H                      
C         (1         32         30      2*1)
      rate=xk(107)*y(1)*y(32)*xnH
      r_f(107,32)=-rate
      r_f(107,30)=rate
      r_f(107,1)=rate
C  108)   H     +   H2O   ->   OH    +   H2                        
C         (1         34         32        2)
      rate=xk(108)*y(1)*y(34)*xnH
      r_f(108,1)=-rate
      r_f(108,34)=-rate
      r_f(108,32)=rate
      r_f(108,2)=rate
C  109)   H     +   H2O   ->   OH    + 2 H                          
C         (1         34         32      2*1)
      rate=xk(109)*y(1)*y(34)*xnH
      r_f(109,34)=-rate
      r_f(109,32)=rate
      r_f(109,1)=rate
C  110)   H     +   C2    ->   CH    +   C                              
C         (1         18         19        17)               
      rate=xk(110)*y(1)*y(18)*xnH
      r_f(110,1)=-rate
      r_f(110,18)=-rate
      r_f(110,19)=rate
      r_f(110,17)=rate
C  111)   H     +   CO    ->   C     +   OH                            
C         (1         33         17        32)   
      rate=xk(111)*y(1)*y(33)*xnH
      r_f(111,1)=-rate
      r_f(111,33)=-rate
      r_f(111,17)=rate
      r_f(111,32)=rate
C  112)   H     +   H2CO  ->   HCO   +   H2                  
C         (1         38         35        2)
      rate=xk(112)*y(1)*y(38)*xnH
      r_f(112,1)=-rate
      r_f(112,38)=-rate
      r_f(112,35)=rate
      r_f(112,2)=rate
C  113)   H     +   O2    ->   OH    +   O                         
C         (1         31         32        30)      
      rate=xk(113)*y(1)*y(31)*xnH
      r_f(113,1)=-rate
      r_f(113,31)=-rate
      r_f(113,32)=rate
      r_f(113,30)=rate
C  114)   H     +   O2    -> 2 O     +   H              
C         (1         31       2*30        1)
      rate=xk(114)*y(1)*y(31)*xnH
      r_f(114,31)=-rate
      r_f(114,30)=2.d0*rate
C  115)   H     +   O2H   ->   H2O   +   O              
C         (1         36         34        30)
      rate=xk(115)*y(1)*y(36)*xnH
      r_f(115,1)=-rate
      r_f(115,36)=-rate
      r_f(115,34)=rate
      r_f(115,30)=rate
C  116)   H     +   O2H   ->   H2    +   O2                    
C         (1         36         2         31)
      rate=xk(116)*y(1)*y(36)*xnH
      r_f(116,1)=-rate
      r_f(116,36)=-rate
      r_f(116,2)=rate
      r_f(116,31)=rate
C  117)   H     +   O2H   -> 2 OH                               
C         (1         36       2*32)
      rate=xk(117)*y(1)*y(36)*xnH
      r_f(117,1)=-rate
      r_f(117,36)=-rate
      r_f(117,32)=2.d0*rate
C  118)   H     +   H2O2  ->   O2H   +   H2                      
C         (1         39         36        2)
      rate=xk(118)*y(1)*y(39)*xnH
      r_f(118,1)=-rate
      r_f(118,39)=-rate
      r_f(118,36)=rate
      r_f(118,2)=rate
C  119)   H     +   CO2   ->   OH    +   CO                        
C         (1         37         32        33)
      rate=xk(119)*y(1)*y(37)*xnH
      r_f(119,1)=-rate
      r_f(119,37)=-rate
      r_f(119,32)=rate
      r_f(119,33)=rate
C  120)   C     +   H2    ->   CH    +   H                    
C         (17        2          19        1)
      rate=xk(120)*y(17)*y(2)*xnH
      r_f(120,17)=-rate
      r_f(120,2)=-rate
      r_f(120,19)=rate
      r_f(120,1)=rate
C  121)   C     +   OH    ->   CH    +   O                     
C         (17        32         19        30)
      rate=xk(121)*y(17)*y(32)*xnH
      r_f(121,17)=-rate
      r_f(121,32)=-rate
      r_f(121,19)=rate
      r_f(121,30)=rate
C  122)   C     +   CO    ->   C2    +   O                             
C         (17        33         18        30)
      rate=xk(122)*y(17)*y(33)*xnH
      r_f(122,17)=-rate
      r_f(122,33)=-rate
      r_f(122,18)=rate
      r_f(122,30)=rate
C  123)   O     +   H2    ->   OH    +   H                        
C         (30        2          32        1)
      rate=xk(123)*y(30)*y(2)*xnH
      r_f(123,30)=-rate
      r_f(123,2)=-rate
      r_f(123,32)=rate
      r_f(123,1)=rate
C  124)   O     +   CH    ->   OH    +   C                         
C         (30        19         32        17) 
      rate=xk(124)*y(30)*y(19)*xnH
      r_f(124,30)=-rate
      r_f(124,19)=-rate
      r_f(124,32)=rate
      r_f(124,17)=rate
C  125)   O     +   CH2   ->   CH    +   OH                       
C         (30        20         19        32)
      rate=xk(125)*y(30)*y(20)*xnH
      r_f(125,30)=-rate
      r_f(125,20)=-rate
      r_f(125,19)=rate
      r_f(125,32)=rate
C  126)   O     +   CH2   ->   HCO   +   H                   
C         (30        20         35        1)
      rate=xk(126)*y(30)*y(20)*xnH
      r_f(126,30)=-rate
      r_f(126,20)=-rate
      r_f(126,35)=rate
      r_f(126,1)=rate
C  127)   O     +   CH4   ->   CH3   +   OH                   
C         (30        22         21        32)
      rate=xk(127)*y(30)*y(22)*xnH
      r_f(127,30)=-rate
      r_f(127,22)=-rate
      r_f(127,21)=rate
      r_f(127,32)=rate
C   128)   O     +   H2O   -> 2 OH   
C         (30        34       2*32)
      rate=xk(128)*y(30)*y(34)*xnH
      r_f(128,30)=-rate
      r_f(128,34)=-rate
      r_f(128,32)=2.d0*rate
C   129)   O     +   H2CO  ->   HCO   +   OH                   
C         (30        38         35        32)
      rate=xk(129)*y(30)*y(38)*xnH
      r_f(129,30)=-rate
      r_f(129,38)=-rate
      r_f(129,35)=rate
      r_f(129,32)=rate
C   130)   O     +   H2O2  ->   O2H   +   OH                 
C         (30        39         36        32)
      rate=xk(130)*y(30)*y(39)*xnH
      r_f(130,30)=-rate
      r_f(130,39)=-rate
      r_f(130,36)=rate
      r_f(130,32)=rate
C   131)   O     +   CO2   ->   CO    +   O2                       
C         (30        37         33        31)
      rate=xk(131)*y(30)*y(37)*xnH
      r_f(131,30)=-rate
      r_f(131,37)=-rate
      r_f(131,33)=rate
      r_f(131,31)=rate
C   132)   H+    +   O     ->   O+    +   H                            
C         (4         30         40        1)
      rate=xk(132)*y(4)*y(30)*xnH
      r_f(132,4)=-rate
      r_f(132,30)=-rate
      r_f(132,40)=rate
      r_f(132,1)=rate
C   133)   H2    +   CH    ->   CH2   +   H                      
C         (2         19         20        1)
      rate=xk(133)*y(2)*y(19)*xnH
      r_f(133,2)=-rate
      r_f(133,19)=-rate
      r_f(133,20)=rate
      r_f(133,1)=rate
C   134)   H2    +   CH    ->   C     +   H     +   H2           
C         (2         19         17        1         2)
      rate=xk(134)*y(2)*y(19)*xnH
      r_f(134,19)=-rate
      r_f(134,17)=rate
      r_f(134,1)=rate
C   135)   H2    +   CH2   ->   CH3   +   H                   
C         (2         20         21        1)
      rate=xk(135)*y(2)*y(20)*xnH
      r_f(135,2)=-rate
      r_f(135,20)=-rate
      r_f(135,21)=rate
      r_f(135,1)=rate
C   136)   H2    +   CH3   ->   CH4   +   H                       
C         (2         21         22        1)
      rate=xk(136)*y(2)*y(21)*xnH
      r_f(136,2)=-rate
      r_f(136,21)=-rate
      r_f(136,22)=rate
      r_f(136,1)=rate
C   137)   H2    +   OH    ->   H2O   +   H               
C         (2         32         34        1)
      rate=xk(137)*y(2)*y(32)*xnH
      r_f(137,2)=-rate
      r_f(137,32)=-rate
      r_f(137,34)=rate
      r_f(137,1)=rate
C   138)   H2    +   OH    ->   O     +   H     +   H2                 
C         (2         32         30        1         2)
      rate=xk(138)*y(2)*y(32)*xnH
      r_f(138,32)=-rate
      r_f(138,30)=rate
      r_f(138,1)=rate
C   139)   H2    +   H2O   ->   OH    +   H     +   H2             
C         (2         34         32        1         2)
      rate=xk(139)*y(2)*y(34)*xnH
      r_f(139,34)=-rate
      r_f(139,32)=rate
      r_f(139,1)=rate
C   140)   H2    +   O2    ->   O2H   +   H                     
C         (2         31         36        1)
      rate=xk(140)*y(2)*y(31)*xnH
      r_f(140,2)=-rate
      r_f(140,31)=-rate
      r_f(140,36)=rate
      r_f(140,1)=rate
C   141)   H2    +   O2    -> 2 OH                         
C         (2         31       2*32)  
      rate=xk(141)*y(2)*y(31)*xnH
      r_f(141,2)=-rate
      r_f(141,31)=-rate
      r_f(141,32)=2.d0*rate
C   142)   H2    +   O2    -> 2 O     +   H2                       
C         (2         31       2*30        2)
      rate=xk(142)*y(2)*y(31)*xnH
      r_f(142,31)=-rate
      r_f(142,30)=2.d0*rate
C   143)   H2    +   O2H   ->   H2O2  +   H                       
C         (2         36         39        1)
      rate=xk(143)*y(2)*y(36)*xnH
      r_f(143,2)=-rate
      r_f(143,36)=-rate
      r_f(143,39)=rate
      r_f(143,1)=rate
C   144)   H2    +   CO2   ->   H2O   +   CO                          
C         (2         37         34        33)   
      rate=xk(144)*y(2)*y(37)*xnH
      r_f(144,2)=-rate
      r_f(144,37)=-rate
      r_f(144,34)=rate
      r_f(144,33)=rate  
C   145)   H3+   +   O2    ->   O2H+  +   H2                   
C         (6         31         46        2)
      rate=xk(145)*y(6)*y(31)*xnH
      r_f(145,6)=-rate
      r_f(145,31)=-rate
      r_f(145,46)=rate
      r_f(145,2)=rate
C   146)   C+    +   H2    ->   CH+   +   H                            
C         (23        2          25        1)
      rate=xk(146)*y(23)*y(2)*xnH
      r_f(146,23)=-rate
      r_f(146,2)=-rate
      r_f(146,25)=rate
      r_f(146,1)=rate
C   147)   CH    +   CH4   ->   CH2   +   CH3                  
C         (19        22         20        21)
      rate=xk(147)*y(19)*y(22)*xnH
      r_f(147,19)=-rate
      r_f(147,22)=-rate
      r_f(147,20)=rate
      r_f(147,21)=rate
C   148)   CH    +   OH    ->   HCO   +   H                    
C         (19        32         35        1)
      rate=xk(148)*y(19)*y(32)*xnH
      r_f(148,19)=-rate
      r_f(148,32)=-rate
      r_f(148,35)=rate
      r_f(148,1)=rate
C   149)   CH    +   HCO   ->   CH2   +   CO                           
C         (19        35         20        33)
      rate=xk(149)*y(19)*y(35)*xnH
      r_f(149,19)=-rate
      r_f(149,35)=-rate
      r_f(149,20)=rate
      r_f(149,33)=rate
C   150)   CH    +   H2CO  ->   CH2   +   HCO                        
C         (19        38         20        35)
      rate=xk(150)*y(19)*y(38)*xnH
      r_f(150,19)=-rate
      r_f(150,38)=-rate
      r_f(150,20)=rate
      r_f(150,35)=rate
C   151)   CH    +   O2    ->   HCO   +   O                         
C         (19        31         35        30)
      rate=xk(151)*y(19)*y(31)*xnH
      r_f(151,19)=-rate
      r_f(151,31)=-rate
      r_f(151,35)=rate
      r_f(151,30)=rate
C   152)   CH    +   O2H   ->   CH2   +   O2                      
C         (19        36         20        31)
      rate=xk(152)*y(19)*y(36)*xnH
      r_f(152,19)=-rate
      r_f(152,36)=-rate
      r_f(152,20)=rate
      r_f(152,31)=rate
C   153)   CH    +   O2H   ->   HCO   +   OH                     
C         (19        36         35        32)
      rate=xk(153)*y(19)*y(36)*xnH
      r_f(153,19)=-rate
      r_f(153,36)=-rate
      r_f(153,35)=rate
      r_f(153,32)=rate
C   154)   CH    +   CO2   ->   HCO   +   CO                      
C         (19        37         35        33)
      rate=xk(154)*y(19)*y(37)*xnH
      r_f(154,19)=-rate
      r_f(154,37)=-rate
      r_f(154,35)=rate
      r_f(154,33)=rate
C   155) 2 CH2             ->   CH    +   CH3                  
C       (2*20                   19        21)
      rate=xk(155)*(y(20)**2)*xnH
      r_f(155,20)=-2.d0*rate
      r_f(155,19)=rate
      r_f(155,21)=rate
C   156)   CH2   +   CH4   -> 2 CH3                            
C         (20        22       2*21)
      rate=xk(156)*y(20)*y(22)*xnH
      r_f(156,20)=-rate
      r_f(156,22)=-rate
      r_f(156,21)=2.d0*rate
C   157)   CH2   +   HCO   ->   CH3   +   CO                   
C         (20        35         21        33)
      rate=xk(157)*y(20)*y(35)*xnH
      r_f(157,20)=-rate
      r_f(157,35)=-rate
      r_f(157,21)=rate
      r_f(157,33)=rate
C   158)   CH2   +   H2CO  ->   CH3   +   HCO                   
C         (20        38         21        35)
      rate=xk(158)*y(20)*y(38)*xnH
      r_f(158,20)=-rate
      r_f(158,38)=-rate
      r_f(158,21)=rate
      r_f(158,35)=rate
C   159)   CH2+  +   H     ->   CH+   +   H2                       
C         (26        1          25        2)
      rate=xk(159)*y(26)*y(1)*xnH
      r_f(159,26)=-rate
      r_f(159,1)=-rate
      r_f(159,25)=rate
      r_f(159,2)=rate
C   160)   CH3   +   H2CO  ->   CH4   +   HCO                  
C         (21        38         22        35)
      rate=xk(160)*y(21)*y(38)*xnH
      r_f(160,21)=-rate
      r_f(160,38)=-rate
      r_f(160,22)=rate
      r_f(160,35)=rate
C   161)   CH3+  +   H     ->   CH2+  +   H2                   
C         (27        1          26        2)
      rate=xk(161)*y(27)*y(1)*xnH
      r_f(161,27)=-rate
      r_f(161,1)=-rate
      r_f(161,26)=rate
      r_f(161,2)=rate
C   162)   OH    +   CH2   ->   CH3   +   O                       
C         (32        20         21        30) 
      rate=xk(162)*y(32)*y(20)*xnH
      r_f(162,32)=-rate
      r_f(162,20)=-rate
      r_f(162,21)=rate
      r_f(162,30)=rate
C   163)   OH    +   CH2   ->   H2O   +   CH                    
C         (32        20         34        19)   
      rate=xk(163)*y(32)*y(20)*xnH
      r_f(163,32)=-rate
      r_f(163,20)=-rate
      r_f(163,34)=rate
      r_f(163,19)=rate
C  164)   OH    +   CH2   ->   H2CO  +   H                     
C         (32        20         38        1)
      rate=xk(164)*y(32)*y(20)*xnH
      r_f(164,32)=-rate
      r_f(164,20)=-rate
      r_f(164,38)=rate
      r_f(164,1)=rate
C  165)   OH    +   CH3   ->   H2O   +   CH2                         
C         (32        21         34        20)
      rate=xk(165)*y(32)*y(21)*xnH
      r_f(165,32)=-rate
      r_f(165,21)=-rate
      r_f(165,34)=rate
      r_f(165,20)=rate
C  166)   OH    +   CH4   ->   H2O   +   CH3                      
C         (32        22         34        21) 
      rate=xk(166)*y(32)*y(22)*xnH
      r_f(166,32)=-rate
      r_f(166,22)=-rate
      r_f(166,34)=rate
      r_f(166,21)=rate
C  167) 2 OH              ->   H2O   +   O                             
C       (2*32                   34        30)
      rate=xk(167)*(y(32)**2)*xnH
      r_f(167,32)=-2.d0*rate
      r_f(167,34)=rate
      r_f(167,30)=rate
C  168)   OH    +   CO    ->   CO2   +   H                       
C         (32        33         37        1)
      rate=xk(168)*y(32)*y(33)*xnH
      r_f(168,32)=-rate
      r_f(168,33)=-rate
      r_f(168,37)=rate
      r_f(168,1)=rate
C  169)   OH    +   H2O2  ->   H2O   +   O2H                        
C         (32        39         34        36)
      rate=xk(169)*y(32)*y(39)*xnH
      r_f(169,32)=-rate
      r_f(169,39)=-rate
      r_f(169,34)=rate
      r_f(169,36)=rate
C  170)   H2O   +   CH3   ->   CH4   +   OH                         
C         (34        21         22        32)
      rate=xk(170)*y(34)*y(21)*xnH
      r_f(170,34)=-rate
      r_f(170,21)=-rate
      r_f(170,22)=rate
      r_f(170,32)=rate
C  171)   CO    +   O2    ->   CO2   +   O                              
C         (33        31         37        30)
      rate=xk(171)*y(33)*y(31)*xnH
      r_f(171,33)=-rate
      r_f(171,31)=-rate
      r_f(171,37)=rate
      r_f(171,30)=rate
C  172)   CO    +   O2H   ->   CO2   +   OH                         
C         (33        36         37        32)
      rate=xk(172)*y(33)*y(36)*xnH
      r_f(172,33)=-rate
      r_f(172,36)=-rate
      r_f(172,37)=rate
      r_f(172,32)=rate
C  173)   O2    +   CH2   ->   HCO   +   OH                        
C         (31        20         35        32)
      rate=xk(173)*y(31)*y(20)*xnH
      r_f(173,31)=-rate
      r_f(173,20)=-rate
      r_f(173,35)=rate
      r_f(173,32)=rate
C  174)   O2    +   CH2   ->   H2CO  +   O                          
C         (31        20         38        30)
      rate=xk(174)*y(31)*y(20)*xnH
      r_f(174,31)=-rate
      r_f(174,20)=-rate
      r_f(174,38)=rate
      r_f(174,30)=rate
C  175)   O2    +   CH3   ->   O2H   +   CH2                         
C         (31        21         36        20)
      rate=xk(175)*y(31)*y(21)*xnH
      r_f(175,31)=-rate
      r_f(175,21)=-rate
      r_f(175,36)=rate
      r_f(175,20)=rate
C  176)   O2    +   CH3   ->   H2CO  +   OH                           
C         (31        21         38        32)
      rate=xk(176)*y(31)*y(21)*xnH
      r_f(176,31)=-rate
      r_f(176,21)=-rate
      r_f(176,38)=rate
      r_f(176,32)=rate
C  177)   O2    +   CH4   ->   CH3   +   O2H                         
C         (31        22         21        36)
      rate=xk(177)*y(31)*y(22)*xnH
      r_f(177,31)=-rate
      r_f(177,22)=-rate
      r_f(177,21)=rate
      r_f(177,36)=rate
C  178)   O2    +   HCO   ->   O2H   +   CO                         
C         (31        35         36        33)
      rate=xk(178)*y(31)*y(35)*xnH
      r_f(178,31)=-rate
      r_f(178,35)=-rate
      r_f(178,36)=rate
      r_f(178,33)=rate
C  179)   O2H   +   CH3   ->   CH4   +   O2                         
C         (36        21         22        31)
      rate=xk(179)*y(36)*y(21)*xnH
      r_f(179,36)=-rate
      r_f(179,21)=-rate
      r_f(179,22)=rate
      r_f(179,31)=rate
C  180)   O2H   +   H2O   ->   H2O2  +   OH                        
C         (36        34         39        32)
      rate=xk(180)*y(36)*y(34)*xnH
      r_f(180,36)=-rate
      r_f(180,34)=-rate
      r_f(180,39)=rate
      r_f(180,32)=rate
C  181)   O2H   +   HCO   ->   H2CO  +   O2                           
C         (36        35         38        31)
      rate=xk(181)*y(36)*y(35)*xnH
      r_f(181,36)=-rate
      r_f(181,35)=-rate
      r_f(181,38)=rate
      r_f(181,31)=rate
C  182)   O2H   +   H2CO  ->   H2O2  +   HCO                            
C         (36        38         39        35)
      rate=xk(182)*y(36)*y(38)*xnH
      r_f(182,36)=-rate
      r_f(182,38)=-rate
      r_f(182,39)=rate
      r_f(182,35)=rate
C  183) 2 O2H             ->   H2O2  +   O2                           
C       (2*36                   39        31)
      rate=xk(183)*(y(36)**2)*xnH
      r_f(183,36)=-2.d0*rate
      r_f(183,39)=rate
      r_f(183,31)=rate
C  184)   H     +   HCO   ->   CO    +   H2                      
C         (1         35         33        2)
      rate=xk(184)*y(1)*y(35)*xnH
      r_f(184,1)=-rate
      r_f(184,35)=-rate
      r_f(184,33)=rate
      r_f(184,2)=rate
C  185)   C     +   H     ->   CH    +   ph.                  
C         (17        1          19)
      rate=xk(185)*y(17)*y(1)*xnH
      r_f(185,17)=-rate
      r_f(185,1)=-rate
      r_f(185,19)=rate
C  186) 2 C               ->   C2    +   ph.         
C       (2*17                   18)
      rate=xk(186)*(y(17)**2)*xnH
      r_f(186,17)=-2.d0*rate
      r_f(186,18)=rate
C  187)   C     +   O     ->   CO    +   ph.             
C         (17        30         33)
      rate=xk(187)*y(17)*y(30)*xnH
      r_f(187,17)=-rate
      r_f(187,30)=-rate
      r_f(187,33)=rate
C  188)   C     +   CH    ->   C2    +   H                  
C         (17        19         18        1)
      rate=xk(188)*y(17)*y(19)*xnH
      r_f(188,17)=-rate
      r_f(188,19)=-rate
      r_f(188,18)=rate
      r_f(188,1)=rate
C  189)   C     +   OH    ->   CO    +   H              
C         (17        32         33        1)
      rate=xk(189)*y(17)*y(32)*xnH
      r_f(189,17)=-rate
      r_f(189,32)=-rate
      r_f(189,33)=rate
      r_f(189,1)=rate
C  190)   C     +   HCO   ->   CH    +   CO                       
C         (17        35         19        33)
      rate=xk(190)*y(17)*y(35)*xnH
      r_f(190,17)=-rate
      r_f(190,35)=-rate
      r_f(190,19)=rate
      r_f(190,33)=rate
C  191)   C     +   O2    ->   CO    +   O                         
C         (17        31         33        30)
      rate=xk(191)*y(17)*y(31)*xnH
      r_f(191,17)=-rate
      r_f(191,31)=-rate
      r_f(191,33)=rate
      r_f(191,30)=rate
C  192)   CH3   +   HCO   ->   CH4   +   CO                        
C         (21        35         22        33)
      rate=xk(192)*y(21)*y(35)*xnH
      r_f(192,21)=-rate
      r_f(192,35)=-rate
      r_f(192,22)=rate
      r_f(192,33)=rate
C  193)   O     +   H     ->   OH    +   ph.        
C         (30        1          32)
      rate=xk(193)*y(30)*y(1)*xnH
      r_f(193,30)=-rate
      r_f(193,1)=-rate
      r_f(193,32)=rate
C  194) 2 O               ->   O2    +   ph.             
C       (2*30                   31)
      rate=xk(194)*(y(30)**2)*xnH
      r_f(194,30)=-2.d0*rate
      r_f(194,31)=rate
C  195)   O     +   CH    ->   CO    +   H                
C         (30        19         33        1)
      rate=xk(195)*y(30)*y(19)*xnH
      r_f(195,30)=-rate
      r_f(195,19)=-rate
      r_f(195,33)=rate
      r_f(195,1)=rate
C  196)   O     +   CH    ->   HCO+  +   e                
C         (30        19         45        3)
      rate=xk(196)*y(30)*y(19)*xnH
      r_f(196,30)=-rate
      r_f(196,19)=-rate
      r_f(196,45)=rate
      r_f(196,3)=rate
C  197)   O     +   CH2   ->   CO    + 2 H                   
C         (30        20         33      2*1)
      rate=xk(197)*y(30)*y(20)*xnH
      r_f(197,30)=-rate
      r_f(197,20)=-rate
      r_f(197,33)=rate
      r_f(197,1)=2.d0*rate
C  198)   O     +   CH3   ->   H2CO  +   H                 
C         (30        21         38        1)
      rate=xk(198)*y(30)*y(21)*xnH
      r_f(198,30)=-rate
      r_f(198,21)=-rate
      r_f(198,38)=rate
      r_f(198,1)=rate
C  199)   O     +   OH    ->   O2    +   H             
C         (30        32         31        1)
      rate=xk(199)*y(30)*y(32)*xnH
      r_f(199,30)=-rate
      r_f(199,32)=-rate
      r_f(199,31)=rate
      r_f(199,1)=rate
C  200)   O     +   C2    ->   CO    +   C                                
C         (30        18         33        17)
      rate=xk(200)*y(30)*y(18)*xnH
      r_f(200,30)=-rate
      r_f(200,18)=-rate
      r_f(200,33)=rate
      r_f(200,17)=rate
C  201)   O     +   HCO   ->   CO2   +   H                    
C         (30        35         37        1)
      rate=xk(201)*y(30)*y(35)*xnH
      r_f(201,30)=-rate
      r_f(201,35)=-rate
      r_f(201,37)=rate
      r_f(201,1)=rate
C  202)   O     +   HCO   ->   OH    +   CO                  
C         (30        35         32        33)
      rate=xk(202)*y(30)*y(35)*xnH
      r_f(202,30)=-rate
      r_f(202,35)=-rate
      r_f(202,32)=rate
      r_f(202,33)=rate
C  203)   O     +   O2H   ->   OH    +   O2                      
C         (30        36         32        31)
      rate=xk(203)*y(30)*y(36)*xnH
      r_f(203,30)=-rate
      r_f(203,36)=-rate
      r_f(203,32)=rate
      r_f(203,31)=rate
C  204)   OH    +   HCO   ->   H2O   +   CO                  
C         (32        35         34        33)
      rate=xk(204)*y(32)*y(35)*xnH
      r_f(204,32)=-rate
      r_f(204,35)=-rate
      r_f(204,34)=rate
      r_f(204,33)=rate
C  205)   OH    +   H2CO  ->   H2O   +   HCO                    
C         (32        38         34        35)
      rate=xk(205)*y(32)*y(38)*xnH
      r_f(205,32)=-rate
      r_f(205,38)=-rate
      r_f(205,34)=rate
      r_f(205,35)=rate
C  206)   OH    +   O2H   ->   H2O   +   O2                    
C         (32        36         34        31)
      rate=xk(206)*y(32)*y(36)*xnH
      r_f(206,32)=-rate
      r_f(206,36)=-rate
      r_f(206,34)=rate
      r_f(206,31)=rate
C  207) 2 HCO             ->   H2CO  +   CO               
C       (2*35                   38        33)
      rate=xk(207)*(y(35)**2)*xnH
      r_f(207,35)=-2.d0*rate
      r_f(207,38)=rate
      r_f(207,33)=rate
C 208)   H+    +   CH    ->   CH+   +   H                       
C         (4         19         25        1)
      rate=xk(208)*y(4)*y(19)*xnH
      r_f(208,4)=-rate
      r_f(208,19)=-rate
      r_f(208,25)=rate
      r_f(208,1)=rate
C 209)   H+    +   CH2   ->   CH+   +   H2                        
C         (4         20         25        2)
      rate=xk(209)*y(4)*y(20)*xnH
      r_f(209,4)=-rate
      r_f(209,20)=-rate
      r_f(209,25)=rate
      r_f(209,2)=rate
C 210)   H+    +   CH2   ->   CH2+  +   H                        
C         (4         20         26        1)
      rate=xk(210)*y(4)*y(20)*xnH
      r_f(210,4)=-rate
      r_f(210,20)=-rate
      r_f(210,26)=rate
      r_f(210,1)=rate
C 211)   H+    +   CH3   ->   CH3+  +   H                        
C         (4         21         27        1)
      rate=xk(211)*y(4)*y(21)*xnH
      r_f(211,4)=-rate
      r_f(211,21)=-rate
      r_f(211,27)=rate
      r_f(211,1)=rate
C 212)   H+    +   CH4   ->   CH3+  +   H2                    
C         (4         22         27        2)
      rate=xk(212)*y(4)*y(22)*xnH
      r_f(212,4)=-rate
      r_f(212,22)=-rate
      r_f(212,27)=rate
      r_f(212,2)=rate
C 213)   H+    +   CH4   ->   CH4+  +   H                    
C         (4         22         28        1)
      rate=xk(213)*y(4)*y(22)*xnH
      r_f(213,4)=-rate
      r_f(213,22)=-rate
      r_f(213,28)=rate
      r_f(213,1)=rate
C 214)   H+    +   OH    ->   OH+   +   H                  
C         (4         32         42        1)
      rate=xk(214)*y(4)*y(32)*xnH
      r_f(214,4)=-rate
      r_f(214,32)=-rate
      r_f(214,42)=rate
      r_f(214,1)=rate
C 215)   H+    +   H2O   ->   H2O+  +   H                   
C         (4         34         44        1)
      rate=xk(215)*y(4)*y(34)*xnH
      r_f(215,4)=-rate
      r_f(215,34)=-rate
      r_f(215,44)=rate
      r_f(215,1)=rate
C 216)   H+    +   C2    ->   C2+   +   H                 
C         (4         18         24        1)
      rate=xk(216)*y(4)*y(18)*xnH
      r_f(216,4)=-rate
      r_f(216,18)=-rate
      r_f(216,24)=rate
      r_f(216,1)=rate
C 217)   H+    +   HCO   ->   CO+   +   H2                      
C         (4         35         43        2)
      rate=xk(217)*y(4)*y(35)*xnH
      r_f(217,4)=-rate
      r_f(217,35)=-rate
      r_f(217,43)=rate
      r_f(217,2)=rate
C 218)   H+    +   HCO   ->   H2+   +   CO                     
C         (4         35         5         33)
      rate=xk(218)*y(4)*y(35)*xnH
      r_f(218,4)=-rate
      r_f(218,35)=-rate
      r_f(218,5)=rate
      r_f(218,33)=rate
C 219)   H+    +   HCO   ->   HCO+  +   H                  
C         (4         35         45        1)
      rate=xk(219)*y(4)*y(35)*xnH
      r_f(219,4)=-rate
      r_f(219,35)=-rate
      r_f(219,45)=rate
      r_f(219,1)=rate
C 220)   H+    +   H2CO  ->   H2CO+ +   H              
C         (4         38         48        1)
      rate=xk(220)*y(4)*y(38)*xnH
      r_f(220,4)=-rate
      r_f(220,38)=-rate
      r_f(220,48)=rate
      r_f(220,1)=rate
C 221)   H+    +   H2CO  ->   HCO+  +   H2                
C         (4         38         45        2)
      rate=xk(221)*y(4)*y(38)*xnH
      r_f(221,4)=-rate
      r_f(221,38)=-rate
      r_f(221,45)=rate
      r_f(221,2)=rate
C 222)   H+    +   O2    ->   O2+   +   H                 
C         (4         31         41        1)
      rate=xk(222)*y(4)*y(31)*xnH
      r_f(222,4)=-rate
      r_f(222,31)=-rate
      r_f(222,41)=rate
      r_f(222,1)=rate
C 223)   H+    +   CO2   ->   HCO+  +   O                    
C         (4         37         45        30)
      rate=xk(223)*y(4)*y(37)*xnH
      r_f(223,4)=-rate
      r_f(223,37)=-rate
      r_f(223,45)=rate
      r_f(223,30)=rate
C 224)   H-    +   C     ->   CH    +   e               
C         (7         17         19        3)
      rate=xk(224)*y(7)*y(17)*xnH
      r_f(224,7)=-rate
      r_f(224,17)=-rate
      r_f(224,19)=rate
      r_f(224,3)=rate
C 225)   H-    +   O     ->   OH    +   e                
C         (7         30         32        3)
      rate=xk(225)*y(7)*y(30)*xnH
      r_f(225,7)=-rate
      r_f(225,30)=-rate
      r_f(225,32)=rate
      r_f(225,3)=rate
C 226)   H-    +   CH    ->   CH2   +   e                  
C         (7         19         20        3)
      rate=xk(226)*y(7)*y(19)*xnH
      r_f(226,7)=-rate
      r_f(226,19)=-rate
      r_f(226,20)=rate
      r_f(226,3)=rate
C 227)   H-    +   CH2   ->   CH3   +   e                
C         (7         20         21        3)
      rate=xk(227)*y(7)*y(20)*xnH
      r_f(227,7)=-rate
      r_f(227,20)=-rate
      r_f(227,21)=rate
      r_f(227,3)=rate
C 228)   H-    +   CH3   ->   CH4   +   e              
C         (7         21         22        3)
      rate=xk(228)*y(7)*y(21)*xnH
      r_f(228,7)=-rate
      r_f(228,21)=-rate
      r_f(228,22)=rate
      r_f(228,3)=rate
C 229)   H-    +   OH    ->   H2O   +   e               
C         (7         32         34        3)
      rate=xk(229)*y(7)*y(32)*xnH
      r_f(229,7)=-rate
      r_f(229,32)=-rate
      r_f(229,34)=rate
      r_f(229,3)=rate
C 230)   H-    +   CO    ->   HCO   +   e                 
C         (7         33         35        3)
      rate=xk(230)*y(7)*y(33)*xnH
      r_f(230,7)=-rate
      r_f(230,33)=-rate
      r_f(230,35)=rate
      r_f(230,3)=rate
C 231)   H-    +   HCO   ->   H2CO  +   e               
C         (7         35         38        3)
      rate=xk(231)*y(7)*y(35)*xnH
      r_f(231,7)=-rate
      r_f(231,35)=-rate
      r_f(231,38)=rate
      r_f(231,3)=rate
C 232)   H2+   +   C     ->   CH+   +   H                
C         (5         17         25        1)
      rate=xk(232)*y(5)*y(17)*xnH
      r_f(232,5)=-rate
      r_f(232,17)=-rate
      r_f(232,25)=rate
      r_f(232,1)=rate
C 233)   H2+   +   O     ->   OH+   +   H              
C         (5         30         42        1)
      rate=xk(233)*y(5)*y(30)*xnH
      r_f(233,5)=-rate
      r_f(233,30)=-rate
      r_f(233,42)=rate
      r_f(233,1)=rate
C 234)   H2+   +   CH    ->   CH+   +   H2               
C         (5         19         25        2)
      rate=xk(234)*y(5)*y(19)*xnH
      r_f(234,5)=-rate
      r_f(234,19)=-rate
      r_f(234,25)=rate
      r_f(234,2)=rate
C 235)   H2+   +   CH    ->   CH2+  +   H              
C         (5         19         26        1)
      rate=xk(235)*y(5)*y(19)*xnH
      r_f(235,5)=-rate
      r_f(235,19)=-rate
      r_f(235,26)=rate
      r_f(235,1)=rate
C 236)   H2+   +   CH2   ->   CH3+  +   H            
C         (5         20         27        1)
      rate=xk(236)*y(5)*y(20)*xnH
      r_f(236,5)=-rate
      r_f(236,20)=-rate
      r_f(236,27)=rate
      r_f(236,1)=rate
C 237)   H2+   +   CH2   ->   CH2+  +   H2              
C         (5         20         26        2)
      rate=xk(237)*y(5)*y(20)*xnH
      r_f(237,5)=-rate
      r_f(237,20)=-rate
      r_f(237,26)=rate
      r_f(237,2)=rate
C 238)   H2+   +   CH4   ->   CH4+  +   H2             
C         (5         22         28        2)
      rate=xk(238)*y(5)*y(22)*xnH
      r_f(238,5)=-rate
      r_f(238,22)=-rate
      r_f(238,28)=rate
      r_f(238,2)=rate
C 239)   H2+   +   CH4   ->   CH5+  +   H                
C         (5         22         29        1)
      rate=xk(239)*y(5)*y(22)*xnH
      r_f(239,5)=-rate
      r_f(239,22)=-rate
      r_f(239,29)=rate
      r_f(239,1)=rate
C 240)   H2+   +   CH4   ->   CH3+  +   H     +  H2          
C         (5         22         27        1        2)
      rate=xk(240)*y(5)*y(22)*xnH
      r_f(240,5)=-rate
      r_f(240,22)=-rate
      r_f(240,27)=rate
      r_f(240,1)=rate
      r_f(240,2)=rate
C 241)   H2+   +   OH    ->   OH+   +   H2            
C         (5         32         42        2)
      rate=xk(241)*y(5)*y(32)*xnH
      r_f(241,5)=-rate
      r_f(241,32)=-rate
      r_f(241,42)=rate
      r_f(241,2)=rate
C 242)   H2+   +   OH    ->   H2O+  +   H               
C         (5         32         44        1)
      rate=xk(242)*y(5)*y(32)*xnH
      r_f(242,5)=-rate
      r_f(242,32)=-rate
      r_f(242,44)=rate
      r_f(242,1)=rate
C 243)   H2+   +   H2O   ->   H2O+  +   H2              
C         (5         34         44        2)
      rate=xk(243)*y(5)*y(34)*xnH
      r_f(243,5)=-rate
      r_f(243,34)=-rate
      r_f(243,44)=rate
      r_f(243,2)=rate
C 244)   H2+   +   H2O   ->   H3O+  +   H               
C         (5         34         47        1)
      rate=xk(244)*y(5)*y(34)*xnH
      r_f(244,5)=-rate
      r_f(244,34)=-rate
      r_f(244,47)=rate
      r_f(244,1)=rate
C 245)   H2+   +   C2    ->   C2+   +   H2              
C         (5         18         24        2)
      rate=xk(245)*y(5)*y(18)*xnH
      r_f(245,5)=-rate
      r_f(245,18)=-rate
      r_f(245,24)=rate
      r_f(245,2)=rate
C 246)   H2+   +   CO    ->   CO+   +   H2                 
C         (5         33         43        2)
      rate=xk(246)*y(5)*y(33)*xnH
      r_f(246,5)=-rate
      r_f(246,33)=-rate
      r_f(246,43)=rate
      r_f(246,2)=rate
C 247)   H2+   +   CO    ->   HCO+  +   H                 
C         (5         33         45        1)
      rate=xk(247)*y(5)*y(33)*xnH
      r_f(247,5)=-rate
      r_f(247,33)=-rate
      r_f(247,45)=rate
      r_f(247,1)=rate
C 248)   H2+   +   HCO   ->   H3+   +   CO                  
C         (5         35         6         33)
      rate=xk(248)*y(5)*y(35)*xnH
      r_f(248,5)=-rate
      r_f(248,35)=-rate
      r_f(248,6)=rate
      r_f(248,33)=rate
C 249)   H2+   +   HCO   ->   HCO+  +   H2               
C         (5         35         45        2)
      rate=xk(249)*y(5)*y(35)*xnH
      r_f(249,5)=-rate
      r_f(249,35)=-rate
      r_f(249,45)=rate
      r_f(249,2)=rate
C 250)   H2+   +   H2CO  ->   H2CO+ +   H2               
C         (5         38         48        2)
      rate=xk(250)*y(5)*y(38)*xnH
      r_f(250,5)=-rate
      r_f(250,38)=-rate
      r_f(250,48)=rate
      r_f(250,2)=rate
C 251)   H2+   +   H2CO  ->   HCO+  +   H     +   H2           
C         (5         38         45        1         2)
      rate=xk(251)*y(5)*y(38)*xnH
      r_f(251,5)=-rate
      r_f(251,38)=-rate
      r_f(251,45)=rate
      r_f(251,1)=rate
      r_f(251,2)=rate
C 252)   H2+   +   O2    ->   O2H+  +   H                   
C         (5         31         46        1)
      rate=xk(252)*y(5)*y(31)*xnH
      r_f(252,5)=-rate
      r_f(252,31)=-rate
      r_f(252,46)=rate
      r_f(252,1)=rate
C 253)   H2+   +   O2    ->   O2+   +   H2                   
C         (5         31         41        2)
      rate=xk(253)*y(5)*y(31)*xnH
      r_f(253,5)=-rate
      r_f(253,31)=-rate
      r_f(253,41)=rate
      r_f(253,2)=rate
C 254)   H2+   +   CO2   ->   HCO2+ +   H                     
C         (5         37         49        1)
      rate=xk(254)*y(5)*y(37)*xnH
      r_f(254,5)=-rate
      r_f(254,37)=-rate
      r_f(254,49)=rate
      r_f(254,1)=rate
C 255)   H3+   +   C     ->   CH+   +   H2                   
C         (6         17         25        2)
      rate=xk(255)*y(6)*y(17)*xnH
      r_f(255,6)=-rate
      r_f(255,17)=-rate
      r_f(255,25)=rate
      r_f(255,2)=rate
C  256)   H3+   +   O     ->   OH+   +   H2                
C         (6         30         42        2) 
      rate=xk(256)*y(6)*y(30)*xnH
      r_f(256,6)=-rate
      r_f(256,30)=-rate
      r_f(256,42)=rate
      r_f(256,2)=rate
C  257)   H3+   +   CH    ->   CH2+  +   H2                
C         (6         19         26        2)
      rate=xk(257)*y(6)*y(19)*xnH
      r_f(257,6)=-rate
      r_f(257,19)=-rate
      r_f(257,26)=rate
      r_f(257,2)=rate
C  258)   H3+   +   CH2   ->   CH3+  +   H2               
C         (6         20         27        2)
      rate=xk(258)*y(6)*y(20)*xnH
      r_f(258,6)=-rate
      r_f(258,20)=-rate
      r_f(258,27)=rate
      r_f(258,2)=rate
C  259)   H3+   +   CH3   ->   CH4+  +   H2              
C         (6         21         28        2)
      rate=xk(259)*y(6)*y(21)*xnH
      r_f(259,6)=-rate
      r_f(259,21)=-rate
      r_f(259,28)=rate
      r_f(259,2)=rate
C  260)   H3+   +   CH4   ->   CH5+  +   H2               
C         (6         22         29        2)
      rate=xk(260)*y(6)*y(22)*xnH
      r_f(260,6)=-rate
      r_f(260,22)=-rate
      r_f(260,29)=rate
      r_f(260,2)=rate
C  261)   H3+   +   OH    ->   H2O+  +   H2               
C         (6         32         44        2)
      rate=xk(261)*y(6)*y(32)*xnH
      r_f(261,6)=-rate
      r_f(261,32)=-rate
      r_f(261,44)=rate
      r_f(261,2)=rate
C  262)   H3+   +   H2O   ->   H3O+  +   H2             
C         (6         34         47        2)
      rate=xk(262)*y(6)*y(34)*xnH
      r_f(262,6)=-rate
      r_f(262,34)=-rate
      r_f(262,47)=rate
      r_f(262,2)=rate
C  263)   H3+   +   CO    ->   HCO+  +   H2            
C         (6         33         45        2)
      rate=xk(263)*y(6)*y(33)*xnH
      r_f(263,6)=-rate
      r_f(263,33)=-rate
      r_f(263,45)=rate
      r_f(263,2)=rate
C  264)   H3+   +   HCO   ->   H2CO+ +   H2              
C         (6         35         48        2)
      rate=xk(264)*y(6)*y(35)*xnH
      r_f(264,6)=-rate
      r_f(264,35)=-rate
      r_f(264,48)=rate
      r_f(264,2)=rate
C  265)   H3+   +   H2CO  ->   H3CO+ +   H2              
C         (6         38         50        2)
      rate=xk(265)*y(6)*y(38)*xnH
      r_f(265,6)=-rate
      r_f(265,38)=-rate
      r_f(265,50)=rate
      r_f(265,2)=rate
C  266)   H3+   +   CO2   ->   HCO2+ +   H2               
C         (6         37         49        2)
      rate=xk(266)*y(6)*y(37)*xnH
      r_f(266,6)=-rate
      r_f(266,37)=-rate
      r_f(266,49)=rate
      r_f(266,2)=rate
C  267)   He+   +   CH    ->   C+    +   H     +   He             
C         (9         19         23        1         8)
      rate=xk(267)*y(9)*y(19)*xnH
      r_f(267,9)=-rate
      r_f(267,19)=-rate
      r_f(267,23)=rate
      r_f(267,1)=rate
      r_f(267,8)=rate
C  268)   He+   +   CH    ->   CH+   +   He                
C         (9         19         25        8)
      rate=xk(268)*y(9)*y(19)*xnH
      r_f(268,9)=-rate
      r_f(268,19)=-rate
      r_f(268,25)=rate
      r_f(268,8)=rate
C  269)   He+   +   CH2   ->   CH+   +   H     +   He            
C         (9         20         25        1         8)
      rate=xk(269)*y(9)*y(20)*xnH
      r_f(269,9)=-rate
      r_f(269,20)=-rate
      r_f(269,25)=rate
      r_f(269,1)=rate
      r_f(269,8)=rate
C  270)   He+   +   CH2   ->   C+    +   H2    +   He            
C         (9         20         23        2         8)
      rate=xk(270)*y(9)*y(20)*xnH
      r_f(270,9)=-rate
      r_f(270,20)=-rate
      r_f(270,23)=rate
      r_f(270,2)=rate
      r_f(270,8)=rate
C  271)   He+   +   CH3   ->   CH+   +   H2    +   He           
C         (9         21         25        2         8)
      rate=xk(271)*y(9)*y(21)*xnH
      r_f(271,9)=-rate
      r_f(271,21)=-rate
      r_f(271,25)=rate
      r_f(271,2)=rate
      r_f(271,8)=rate
C  272)   He+   +   CH4   ->   CH3   +   H+    +   He               
C         (9         22         21        4         8)
      rate=xk(272)*y(9)*y(22)*xnH
      r_f(272,9)=-rate
      r_f(272,22)=-rate
      r_f(272,21)=rate
      r_f(272,4)=rate
      r_f(272,8)=rate
C  273)   He+   +   CH4   ->   CH+   +   H2    +   He    +   H         
C         (9         22         25        2         8         1)
      rate=xk(273)*y(9)*y(22)*xnH
      r_f(273,9)=-rate
      r_f(273,22)=-rate
      r_f(273,25)=rate
      r_f(273,2)=rate
      r_f(273,8)=rate
      r_f(273,1)=rate
C  274)   He+   +   CH4   ->   CH2+  +   H2    +   He            
C         (9         22         26        2         8)
      rate=xk(274)*y(9)*y(22)*xnH
      r_f(274,9)=-rate
      r_f(274,22)=-rate
      r_f(274,26)=rate
      r_f(274,2)=rate
      r_f(274,8)=rate
C  275)   He+   +   CH4   ->   CH3+  +   He    +   H           
C         (9         22         27        8         1)
      rate=xk(275)*y(9)*y(22)*xnH
      r_f(275,9)=-rate
      r_f(275,22)=-rate
      r_f(275,27)=rate
      r_f(275,8)=rate
      r_f(275,1)=rate
C  276)   He+   +   CH4   ->   CH4+  +   He                
C         (9         22         28        8)
      rate=xk(276)*y(9)*y(22)*xnH
      r_f(276,9)=-rate
      r_f(276,22)=-rate
      r_f(276,28)=rate
      r_f(276,8)=rate
C  277)   He+   +   OH    ->   O+    +   H     +   He          
C         (9         32         40        1         8)
      rate=xk(277)*y(9)*y(32)*xnH
      r_f(277,9)=-rate
      r_f(277,32)=-rate
      r_f(277,40)=rate
      r_f(277,1)=rate
      r_f(277,8)=rate
C  278)   He+   +   H2O   ->   H+    +   OH    +   He           
C         (9         34         4         32        8)
      rate=xk(278)*y(9)*y(34)*xnH
      r_f(278,9)=-rate
      r_f(278,34)=-rate
      r_f(278,4)=rate
      r_f(278,32)=rate
      r_f(278,8)=rate
C  279)   He+   +   H2O   ->   OH+   +   H     +   He            
C         (9         34         42        1         8) 
      rate=xk(279)*y(9)*y(34)*xnH
      r_f(279,9)=-rate
      r_f(279,34)=-rate
      r_f(279,42)=rate
      r_f(279,1)=rate
      r_f(279,8)=rate
C  280)   He+   +   H2O   ->   H2O+  +   He                 
C         (9         34         44        8)
      rate=xk(280)*y(9)*y(34)*xnH
      r_f(280,9)=-rate
      r_f(280,34)=-rate
      r_f(280,44)=rate
      r_f(280,8)=rate
C  281)   He+   +   C2    ->   C+    +   C     +   He       
C         (9         18         23        17        8)
      rate=xk(281)*y(9)*y(18)*xnH
      r_f(281,9)=-rate
      r_f(281,18)=-rate
      r_f(281,23)=rate
      r_f(281,17)=rate
      r_f(281,8)=rate
C  282)   He+   +   C2    ->   C2+   +   He                  
C         (9         18         24        8)
      rate=xk(282)*y(9)*y(18)*xnH
      r_f(282,9)=-rate
      r_f(282,18)=-rate
      r_f(282,24)=rate
      r_f(282,8)=rate
C  283)   He+   +   CO    ->   C+    +   O     +   He            
C         (9         33         23        30        8)
      rate=xk(283)*y(9)*y(33)*xnH
      r_f(283,9)=-rate
      r_f(283,33)=-rate
      r_f(283,23)=rate
      r_f(283,30)=rate
      r_f(283,8)=rate
C  284)   He+   +   HCO   ->   CO+   +   H     +   He           
C         (9         35         43        1         8)
      rate=xk(284)*y(9)*y(35)*xnH
      r_f(284,9)=-rate
      r_f(284,35)=-rate
      r_f(284,43)=rate
      r_f(284,1)=rate
      r_f(284,8)=rate
C  285)   He+   +   HCO   ->   CH+   +   O     +   He          
C         (9         35         25        30        8)
      rate=xk(285)*y(9)*y(35)*xnH
      r_f(285,9)=-rate
      r_f(285,35)=-rate
      r_f(285,25)=rate
      r_f(285,30)=rate
      r_f(285,8)=rate
C  286)   He+   +   HCO   ->   HeH+  +   CO                        
C         (9         35         11        33)
      rate=xk(286)*y(9)*y(35)*xnH
      r_f(286,9)=-rate
      r_f(286,35)=-rate
      r_f(286,11)=rate
      r_f(286,33)=rate
C  287)   He+   +   H2CO  ->   CO+   +   H2    +   He             
C         (9         38         43        2         8)
      rate=xk(287)*y(9)*y(38)*xnH
      r_f(287,9)=-rate
      r_f(287,38)=-rate
      r_f(287,43)=rate
      r_f(287,2)=rate
      r_f(287,8)=rate
C  288)   He+   +   H2CO  ->   HCO+  +   H     +   He           
C         (9         38         45        1         8)
      rate=xk(288)*y(9)*y(38)*xnH
      r_f(288,9)=-rate
      r_f(288,38)=-rate
      r_f(288,45)=rate
      r_f(288,1)=rate
      r_f(288,8)=rate
C  289)   He+   +   O2    ->   O+    +   O     +   He           
C         (9         31         40        30        8)
      rate=xk(289)*y(9)*y(31)*xnH
      r_f(289,9)=-rate
      r_f(289,31)=-rate
      r_f(289,40)=rate
      r_f(289,30)=rate
      r_f(289,8)=rate
C  290)   He+   +   O2    ->   O2+   +   He                  
C         (9         31         41        8)
      rate=xk(290)*y(9)*y(31)*xnH
      r_f(290,9)=-rate
      r_f(290,31)=-rate
      r_f(290,41)=rate
      r_f(290,8)=rate
C  291)   He+   +   CO2   ->   O2+   +   C     +   He           
C         (9         37         41        17        8)
      rate=xk(291)*y(9)*y(37)*xnH
      r_f(291,9)=-rate
      r_f(291,37)=-rate
      r_f(291,41)=rate
      r_f(291,17)=rate
      r_f(291,8)=rate
C  292)   He+   +   CO2   ->   O+    +   CO    +   He            
C         (9         37         40        33        8)
      rate=xk(292)*y(9)*y(37)*xnH
      r_f(292,9)=-rate
      r_f(292,37)=-rate
      r_f(292,40)=rate
      r_f(292,33)=rate
      r_f(292,8)=rate
C  293)   He+   +   CO2   ->   CO+   +   O     +   He             
C         (9         37         43        30        8)  
      rate=xk(293)*y(9)*y(37)*xnH
      r_f(293,9)=-rate
      r_f(293,37)=-rate
      r_f(293,43)=rate
      r_f(293,30)=rate
      r_f(293,8)=rate
C  294)   He+   +   CO2   ->   C+    +   O2    +   He        
C         (9         37         23        31        8)
      rate=xk(294)*y(9)*y(37)*xnH
      r_f(294,9)=-rate
      r_f(294,37)=-rate
      r_f(294,23)=rate
      r_f(294,31)=rate
      r_f(294,8)=rate
C  295)   C+    +   H     ->   CH+   +   ph.              
C         (23        1          25)
      rate=xk(295)*y(23)*y(1)*xnH
      r_f(295,23)=-rate
      r_f(295,1)=-rate
      r_f(295,25)=rate
C  296)   C+    +   O     ->   CO+   +   ph.                
C         (23        30         43)
      rate=xk(296)*y(23)*y(30)*xnH
      r_f(296,23)=-rate
      r_f(296,30)=-rate
      r_f(296,43)=rate
C  297)   C+    +   H-    ->   H     +   C                
C         (23        7          1         17)
      rate=xk(297)*y(23)*y(7)*xnH
      r_f(297,23)=-rate
      r_f(297,7)=-rate
      r_f(297,1)=rate
      r_f(297,17)=rate
C  298)   C+    +   H2    ->   CH2+  +   ph.              
C         (23        2          26)
      rate=xk(298)*y(23)*y(2)*xnH
      r_f(298,23)=-rate
      r_f(298,2)=-rate
      r_f(298,26)=rate
C  299)   C+    +   CH    ->   C2+   +   H                 
C         (23        19         24        1)
      rate=xk(299)*y(23)*y(19)*xnH
      r_f(299,23)=-rate
      r_f(299,19)=-rate
      r_f(299,24)=rate
      r_f(299,1)=rate
C  300)   C+    +   CH    ->   CH+   +   C                
C         (23        19         25        17)
      rate=xk(300)*y(23)*y(19)*xnH
      r_f(300,23)=-rate
      r_f(300,19)=-rate
      r_f(300,25)=rate
      r_f(300,17)=rate
C  301)   C+    +   CH2   ->   CH2+  +   C                 
C         (23        20         26        17)
      rate=xk(301)*y(23)*y(20)*xnH
      r_f(301,23)=-rate
      r_f(301,20)=-rate
      r_f(301,26)=rate
      r_f(301,17)=rate
C  302)   C+    +   OH    ->   CO+   +   H                           
C         (23        32         43        1)
      rate=xk(302)*y(23)*y(32)*xnH
      r_f(302,23)=-rate
      r_f(302,32)=-rate
      r_f(302,43)=rate
      r_f(302,1)=rate
C  303)   C+    +   H2O   ->   HCO+  +   H                   
C         (23        34         45        1)
      rate=xk(303)*y(23)*y(34)*xnH
      r_f(303,23)=-rate
      r_f(303,34)=-rate
      r_f(303,45)=rate
      r_f(303,1)=rate
C  304)   C+    +   HCO   ->   HCO+  +   C                     
C         (23        35         45        17)
      rate=xk(304)*y(23)*y(35)*xnH
      r_f(304,23)=-rate
      r_f(304,35)=-rate
      r_f(304,45)=rate
      r_f(304,17)=rate
C  305)   C+    +   HCO   ->   CH+   +   CO                    
C         (23        35         25        33)
      rate=xk(305)*y(23)*y(35)*xnH
      r_f(305,23)=-rate
      r_f(305,35)=-rate
      r_f(305,25)=rate
      r_f(305,33)=rate
C  306)   C+    +   H2CO  ->   CH2+  +   CO                      
C         (23        38         26        33)
      rate=xk(306)*y(23)*y(38)*xnH
      r_f(306,23)=-rate
      r_f(306,38)=-rate
      r_f(306,26)=rate
      r_f(306,33)=rate
C  307)   C+    +   H2CO  ->   HCO+  +   CH                       
C         (23        38         45        19)
      rate=xk(307)*y(23)*y(38)*xnH
      r_f(307,23)=-rate
      r_f(307,38)=-rate
      r_f(307,45)=rate
      r_f(307,19)=rate
C  308)   C+    +   H2CO  ->   H2CO+ +   C                        
C         (23        38         48        17)   
      rate=xk(308)*y(23)*y(38)*xnH
      r_f(308,23)=-rate
      r_f(308,38)=-rate
      r_f(308,48)=rate
      r_f(308,17)=rate                       
C  309)   C+    +   O2    ->   CO+   +   O                           
C         (23        31         43        30)
      rate=xk(309)*y(23)*y(31)*xnH
      r_f(309,23)=-rate
      r_f(309,31)=-rate
      r_f(309,43)=rate
      r_f(309,30)=rate
C  310)   C+    +   O2    ->   O+    +   CO                   
C         (23        31         40        33)
      rate=xk(310)*y(23)*y(31)*xnH
      r_f(310,23)=-rate
      r_f(310,31)=-rate
      r_f(310,40)=rate
      r_f(310,33)=rate
C  311)   C+    +   CO2   ->   CO+   +   CO                        
C         (23        37         43        33)
      rate=xk(311)*y(23)*y(37)*xnH
      r_f(311,23)=-rate
      r_f(311,37)=-rate
      r_f(311,43)=rate
      r_f(311,33)=rate
C  312)   CH+   +   H     ->   C+    +   H2                 
C         (25        1          23        2)
      rate=xk(312)*y(25)*y(1)*xnH
      r_f(312,25)=-rate
      r_f(312,1)=-rate
      r_f(312,23)=rate
      r_f(312,2)=rate
C  313)   CH+   +   C     ->   C2+   +   H                       
C         (25        17         24        1)
      rate=xk(313)*y(25)*y(17)*xnH
      r_f(313,25)=-rate
      r_f(313,17)=-rate
      r_f(313,24)=rate
      r_f(313,1)=rate
C  314)   CH+   +   O     ->   CO+   +   H                     
C         (25        30         43        1)
      rate=xk(314)*y(25)*y(30)*xnH
      r_f(314,25)=-rate
      r_f(314,30)=-rate
      r_f(314,43)=rate
      r_f(314,1)=rate
C  315)   CH+   +   H2    ->   CH2+  +   H                       
C         (25        2          26        1)
      rate=xk(315)*y(25)*y(2)*xnH
      r_f(315,25)=-rate
      r_f(315,2)=-rate
      r_f(315,26)=rate
      r_f(315,1)=rate
C  316)   CH+   +   CH    ->   C2+   +   H2                     
C         (25        19         24        2) 
      rate=xk(316)*y(25)*y(19)*xnH
      r_f(316,25)=-rate
      r_f(316,19)=-rate
      r_f(316,24)=rate
      r_f(316,2)=rate
C  317)   CH+   +   OH    ->   CO+   +   H2                        
C         (25        32         43        2)
      rate=xk(317)*y(25)*y(32)*xnH
      r_f(317,25)=-rate
      r_f(317,32)=-rate
      r_f(317,43)=rate
      r_f(317,2)=rate
C  318)   CH+   +   H2O   ->   HCO+  +   H2                          
C         (25        34         45        2)
      rate=xk(318)*y(25)*y(34)*xnH
      r_f(318,25)=-rate
      r_f(318,34)=-rate
      r_f(318,45)=rate
      r_f(318,2)=rate
C  319)   CH+   +   H2O   ->   H3O+  +   C                     
C         (25        34         47        17)
      rate=xk(319)*y(25)*y(34)*xnH
      r_f(319,25)=-rate
      r_f(319,34)=-rate
      r_f(319,47)=rate
      r_f(319,17)=rate
C  320)   CH+   +   H2O   ->   H2CO+ +   H                  
C         (25        34         48        1)
      rate=xk(320)*y(25)*y(34)*xnH
      r_f(320,25)=-rate
      r_f(320,34)=-rate
      r_f(320,48)=rate
      r_f(320,1)=rate
C  321)   CH+   +   CO    ->   HCO+  +   C                      
C         (25        33         45        17)
      rate=xk(321)*y(25)*y(33)*xnH
      r_f(321,25)=-rate
      r_f(321,33)=-rate
      r_f(321,45)=rate
      r_f(321,17)=rate
C  322)   CH+   +   HCO   ->   CH2+  +   CO                   
C         (25        35         26        33)
      rate=xk(322)*y(25)*y(35)*xnH
      r_f(322,25)=-rate
      r_f(322,35)=-rate
      r_f(322,26)=rate
      r_f(322,33)=rate
C  323)   CH+   +   HCO   ->   HCO+  +   CH                 
C         (25        35         45        19)
      rate=xk(323)*y(25)*y(35)*xnH
      r_f(323,25)=-rate
      r_f(323,35)=-rate
      r_f(323,45)=rate
      r_f(323,19)=rate
C  324)   CH+   +   H2CO  ->   CH3+  +   CO               
C         (25        38         27        33)
      rate=xk(324)*y(25)*y(38)*xnH
      r_f(324,25)=-rate
      r_f(324,38)=-rate
      r_f(324,27)=rate
      r_f(324,33)=rate
C  325)   CH+   +   H2CO  ->   HCO+  +   CH2                   
C         (25        38         45        20)
      rate=xk(325)*y(25)*y(38)*xnH
      r_f(325,25)=-rate
      r_f(325,38)=-rate
      r_f(325,45)=rate
      r_f(325,20)=rate
C  326)   CH+   +   H2CO  ->   H3CO+ +   C                          
C         (25        38         50        17)    
      rate=xk(326)*y(25)*y(38)*xnH
      r_f(326,25)=-rate
      r_f(326,38)=-rate
      r_f(326,50)=rate
      r_f(326,17)=rate                              
C  327)   CH+   +   O2    ->   HCO+  +   O                         
C         (25        31         45        30)
      rate=xk(327)*y(25)*y(31)*xnH
      r_f(327,25)=-rate
      r_f(327,31)=-rate
      r_f(327,45)=rate
      r_f(327,30)=rate
C  328)   CH+   +   O2    ->   HCO   +   O+                     
C         (25        31         35        40)
      rate=xk(328)*y(25)*y(31)*xnH
      r_f(328,25)=-rate
      r_f(328,31)=-rate
      r_f(328,35)=rate
      r_f(328,40)=rate
C  329)   CH+   +   O2    ->   CO+   +   OH                         
C         (25        31         43        32)
      rate=xk(329)*y(25)*y(31)*xnH
      r_f(329,25)=-rate
      r_f(329,31)=-rate
      r_f(329,43)=rate
      r_f(329,32)=rate
C  330)   CH+   +   CO2   ->   HCO+  +   CO                    
C         (25        37         45        33)
      rate=xk(330)*y(25)*y(37)*xnH
      r_f(330,25)=-rate
      r_f(330,37)=-rate
      r_f(330,45)=rate
      r_f(330,33)=rate
C  331)   CH2+  +   O     ->   HCO+  +   H            
C         (26        30         45        1)
      rate=xk(331)*y(26)*y(30)*xnH
      r_f(331,26)=-rate
      r_f(331,30)=-rate
      r_f(331,45)=rate
      r_f(331,1)=rate
C  332)   CH2+  +   H2    ->   CH3+  +   H               
C         (26        2          27        1)
      rate=xk(332)*y(26)*y(2)*xnH
      r_f(332,26)=-rate
      r_f(332,2)=-rate
      r_f(332,27)=rate
      r_f(332,1)=rate
C  333)   CH2+  +   H2O   ->   H3CO+ +   H                   
C         (26        34         50        1)
      rate=xk(333)*y(26)*y(34)*xnH
      r_f(333,26)=-rate
      r_f(333,34)=-rate
      r_f(333,50)=rate
      r_f(333,1)=rate
C  334)   CH2+  +   HCO   ->   CH3+  +   CO                 
C         (26        35         27        33)
      rate=xk(334)*y(26)*y(35)*xnH
      r_f(334,26)=-rate
      r_f(334,35)=-rate
      r_f(334,27)=rate
      r_f(334,33)=rate
C  335)   CH2+  +   H2CO  ->   HCO+  +   CH3                       
C         (26        38         45        21)
      rate=xk(335)*y(26)*y(38)*xnH
      r_f(335,26)=-rate
      r_f(335,38)=-rate
      r_f(335,45)=rate
      r_f(335,21)=rate
C  336)   CH2+  +   O2    ->   HCO+  +   OH                    
C         (26        31         45        32)
      rate=xk(336)*y(26)*y(31)*xnH
      r_f(336,26)=-rate
      r_f(336,31)=-rate
      r_f(336,45)=rate
      r_f(336,32)=rate
C  337)   CH2+  +   CO2   ->   H2CO+ +   CO             
C         (26        37         48        33)
      rate=xk(337)*y(26)*y(37)*xnH
      r_f(337,26)=-rate
      r_f(337,37)=-rate
      r_f(337,48)=rate
      r_f(337,33)=rate
C  338)   CH3+  +   O     ->   HCO+  +   H2              
C         (27        30         45        2)
      rate=xk(338)*y(27)*y(30)*xnH
      r_f(338,27)=-rate
      r_f(338,30)=-rate
      r_f(338,45)=rate
      r_f(338,2)=rate
C  339)   CH3+  +   O     ->   H2CO+ +   H                  
C         (27        30         48        1)
      rate=xk(339)*y(27)*y(30)*xnH
      r_f(339,27)=-rate
      r_f(339,30)=-rate
      r_f(339,48)=rate
      r_f(339,1)=rate
C  340)   CH3+  +   H2    ->   CH5+  +   ph.                 
C         (27        2          29)
      rate=xk(340)*y(27)*y(2)*xnH
      r_f(340,27)=-rate
      r_f(340,2)=-rate
      r_f(340,29)=rate
C  341)   CH3+  +   OH    ->   H2CO+ +   H2                 
C         (27        32         48        2)
      rate=xk(341)*y(27)*y(32)*xnH
      r_f(341,27)=-rate
      r_f(341,32)=-rate
      r_f(341,48)=rate
      r_f(341,2)=rate
C  342)   CH3+  +   HCO   ->   CH4+  +   CO                         
C         (27        35         28        33)
      rate=xk(342)*y(27)*y(35)*xnH
      r_f(342,27)=-rate
      r_f(342,35)=-rate
      r_f(342,28)=rate
      r_f(342,33)=rate
C  343)   CH3+  +   HCO   ->   HCO+  +   CH3           
C         (27        35         45        21)
      rate=xk(343)*y(27)*y(35)*xnH
      r_f(343,27)=-rate
      r_f(343,35)=-rate
      r_f(343,45)=rate
      r_f(343,21)=rate
C  344)   CH3+  +   H2CO  ->   HCO+  +   CH4         
C         (27        38         45        22) 
      rate=xk(344)*y(27)*y(38)*xnH
      r_f(344,27)=-rate
      r_f(344,38)=-rate
      r_f(344,45)=rate
      r_f(344,22)=rate
C  345)   CH3+  +   O2    ->   H3CO+ +   O                 
C         (27        31         50        30)
      rate=xk(345)*y(27)*y(31)*xnH
      r_f(345,27)=-rate
      r_f(345,31)=-rate
      r_f(345,50)=rate
      r_f(345,30)=rate
C  346)   O+    +   H     ->   H+    +   O            
C         (40        1          4         30)
      rate=xk(346)*y(40)*y(1)*xnH
      r_f(346,40)=-rate
      r_f(346,1)=-rate
      r_f(346,4)=rate
      r_f(346,30)=rate
C  347)   O+    +   H-    ->   H     +   O            
C         (40        7          1         30)
      rate=xk(347)*y(40)*y(7)*xnH
      r_f(347,40)=-rate
      r_f(347,7)=-rate
      r_f(347,1)=rate
      r_f(347,30)=rate
C  348)   O+    +   H2    ->   OH+   +   H                
C         (40        2          42        1)
      rate=xk(348)*y(40)*y(2)*xnH
      r_f(348,40)=-rate
      r_f(348,2)=-rate
      r_f(348,42)=rate
      r_f(348,1)=rate
C  349)   O+    +   CH    ->   CH+   +   O             
C         (40        19         25        30)
      rate=xk(349)*y(40)*y(19)*xnH
      r_f(349,40)=-rate
      r_f(349,19)=-rate
      r_f(349,25)=rate
      r_f(349,30)=rate
C  350)   O+    +   CH    ->   CO+   +   H            
C         (40        19         43        1)
      rate=xk(350)*y(40)*y(19)*xnH
      r_f(350,40)=-rate
      r_f(350,19)=-rate
      r_f(350,43)=rate
      r_f(350,1)=rate
C  351)   O+    +   CH2   ->   CH2+  +   O            
C         (40        20         26        30)
      rate=xk(351)*y(40)*y(20)*xnH
      r_f(351,40)=-rate
      r_f(351,20)=-rate
      r_f(351,26)=rate
      r_f(351,30)=rate
C  352)   O+    +   CH4   ->   CH3+  +   OH            
C         (40        22         27        32)
      rate=xk(352)*y(40)*y(22)*xnH
      r_f(352,40)=-rate
      r_f(352,22)=-rate
      r_f(352,27)=rate
      r_f(352,32)=rate
C  353)   O+    +   CH4   ->   CH4+  +   O             
C         (40        22         28        30)
      rate=xk(353)*y(40)*y(22)*xnH
      r_f(353,40)=-rate
      r_f(353,22)=-rate
      r_f(353,28)=rate
      r_f(353,30)=rate
C  354)   O+    +   OH    ->   OH+   +   O               
C         (40        32         42        30)
      rate=xk(354)*y(40)*y(32)*xnH
      r_f(354,40)=-rate
      r_f(354,32)=-rate
      r_f(354,42)=rate
      r_f(354,30)=rate
C  355)   O+    +   OH    ->   O2+   +   H             
C         (40        32         41        1)
      rate=xk(355)*y(40)*y(32)*xnH
      r_f(355,40)=-rate
      r_f(355,32)=-rate
      r_f(355,41)=rate
      r_f(355,1)=rate
C  356)   O+    +   H2O   ->   H2O+  +   O             
C         (40        34         44        30)
      rate=xk(356)*y(40)*y(34)*xnH
      r_f(356,40)=-rate
      r_f(356,34)=-rate
      r_f(356,44)=rate
      r_f(356,30)=rate
C  357)   O+    +   C2    ->   C2+   +   O            
C         (40        18         24        30)
      rate=xk(357)*y(40)*y(18)*xnH
      r_f(357,40)=-rate
      r_f(357,18)=-rate
      r_f(357,24)=rate
      r_f(357,30)=rate
C  358)   O+    +   C2    ->   CO+   +   C               
C         (40        18         43        17)
      rate=xk(358)*y(40)*y(18)*xnH
      r_f(358,40)=-rate
      r_f(358,18)=-rate
      r_f(358,43)=rate
      r_f(358,17)=rate
C  359)   O+    +   HCO   ->   OH+   +   CO            
C         (40        35         42        33)
      rate=xk(359)*y(40)*y(35)*xnH
      r_f(359,40)=-rate
      r_f(359,35)=-rate
      r_f(359,42)=rate
      r_f(359,33)=rate
C  360)   O+    +   HCO   ->   HCO+  +   O              
C         (40        35         45        30)
      rate=xk(360)*y(40)*y(35)*xnH
      r_f(360,40)=-rate
      r_f(360,35)=-rate
      r_f(360,45)=rate
      r_f(360,30)=rate
C  361)   O+    +   H2CO  ->   HCO+  +   OH             
C         (40        38         45        32)
      rate=xk(361)*y(40)*y(38)*xnH
      r_f(361,40)=-rate
      r_f(361,38)=-rate
      r_f(361,45)=rate
      r_f(361,32)=rate
C  362)   O+    +   H2CO  ->   H2CO+ +   O                 
C         (40        38         48        30)
      rate=xk(362)*y(40)*y(38)*xnH
      r_f(362,40)=-rate
      r_f(362,38)=-rate
      r_f(362,48)=rate
      r_f(362,30)=rate
C  363)   O+    +   O2    ->   O2+   +   O             
C         (40        31         41        30) 
      rate=xk(363)*y(40)*y(31)*xnH
      r_f(363,40)=-rate
      r_f(363,31)=-rate
      r_f(363,41)=rate
      r_f(363,30)=rate
C  364)   O+    +   CO2   ->   O2+   +   CO                
C         (40        37         41        33)
      rate=xk(364)*y(40)*y(37)*xnH
      r_f(364,40)=-rate
      r_f(364,37)=-rate
      r_f(364,41)=rate
      r_f(364,33)=rate
C  365)   CH4+  +   H     ->   CH3+  +   H2              
C         (28        1          27        2) 
      rate=xk(365)*y(28)*y(1)*xnH
      r_f(365,28)=-rate
      r_f(365,1)=-rate
      r_f(365,27)=rate
      r_f(365,2)=rate
C  366)   CH4+  +   O     ->   CH3+  +   OH              
C         (28        30         27        32)
      rate=xk(366)*y(28)*y(30)*xnH
      r_f(366,28)=-rate
      r_f(366,30)=-rate
      r_f(366,27)=rate
      r_f(366,32)=rate
C  367)   CH4+  +   H2    ->   CH5+  +   H                 
C         (28        2          29        1)
      rate=xk(367)*y(28)*y(2)*xnH
      r_f(367,28)=-rate
      r_f(367,2)=-rate
      r_f(367,29)=rate
      r_f(367,1)=rate
C  368)   CH4+  +   CH4   ->   CH5+  +   CH3              
C         (28        22         29        21)
      rate=xk(368)*y(28)*y(22)*xnH
      r_f(368,28)=-rate
      r_f(368,22)=-rate
      r_f(368,29)=rate
      r_f(368,21)=rate
C  369)   CH4+  +   H2O   ->   H3O+  +   CH3             
C         (28        34         47        21)
      rate=xk(369)*y(28)*y(34)*xnH
      r_f(369,28)=-rate
      r_f(369,34)=-rate
      r_f(369,47)=rate
      r_f(369,21)=rate
C  370)   CH4+  +   CO    ->   HCO+  +   CH3              
C         (28        33         45        21)
      rate=xk(370)*y(28)*y(33)*xnH
      r_f(370,28)=-rate
      r_f(370,33)=-rate
      r_f(370,45)=rate
      r_f(370,21)=rate
C  371)   CH4+  +   H2CO  ->   H2CO+ +   CH4                
C         (28        38         48        22)
      rate=xk(371)*y(28)*y(38)*xnH
      r_f(371,28)=-rate
      r_f(371,38)=-rate
      r_f(371,48)=rate
      r_f(371,22)=rate
C  372)   CH4+  +   H2CO  ->   H3CO+ +   CH3               
C         (28        38         50        21)
      rate=xk(372)*y(28)*y(38)*xnH
      r_f(372,28)=-rate
      r_f(372,38)=-rate
      r_f(372,50)=rate
      r_f(372,21)=rate
C  373)   CH4+  +   O2    ->   O2+   +   CH4              
C         (28        31         41        22)
      rate=xk(373)*y(28)*y(31)*xnH
      r_f(373,28)=-rate
      r_f(373,31)=-rate
      r_f(373,41)=rate
      r_f(373,22)=rate
C  374)   CH4+  +   CO2   ->   HCO2+ +   CH3            
C         (28        37         49        21)
      rate=xk(374)*y(28)*y(37)*xnH
      r_f(374,28)=-rate
      r_f(374,37)=-rate
      r_f(374,49)=rate
      r_f(374,21)=rate
C  375)   OH+   +   C     ->   CH+   +   O                  
C         (42        17         25        30)
      rate=xk(375)*y(42)*y(17)*xnH
      r_f(375,42)=-rate
      r_f(375,17)=-rate
      r_f(375,25)=rate
      r_f(375,30)=rate
C  376)   OH+   +   O     ->   O2+   +   H               
C         (42        30         41        1)
      rate=xk(376)*y(42)*y(30)*xnH
      r_f(376,42)=-rate
      r_f(376,30)=-rate
      r_f(376,41)=rate
      r_f(376,1)=rate
C  377)   OH+   +   H2    ->   H2O+  +   H                
C         (42        2          44        1)
      rate=xk(377)*y(42)*y(2)*xnH
      r_f(377,42)=-rate
      r_f(377,2)=-rate
      r_f(377,44)=rate
      r_f(377,1)=rate
C  378)   OH+   +   CH    ->   CH+   +   OH                
C         (42        19         25        32)
      rate=xk(378)*y(42)*y(19)*xnH
      r_f(378,42)=-rate
      r_f(378,19)=-rate
      r_f(378,25)=rate
      r_f(378,32)=rate
C  379)   OH+   +   CH    ->   CH2+  +   O                
C         (42        19         26        30)
      rate=xk(379)*y(42)*y(19)*xnH
      r_f(379,42)=-rate
      r_f(379,19)=-rate
      r_f(379,26)=rate
      r_f(379,30)=rate
C  380)   OH+   +   CH2   ->   CH2+  +   OH                  
C         (42        20         26        32)
      rate=xk(380)*y(42)*y(20)*xnH
      r_f(380,42)=-rate
      r_f(380,20)=-rate
      r_f(380,26)=rate
      r_f(380,32)=rate
C  381)   OH+   +   CH2   ->   CH3+  +   O                  
C         (42        20         27        30)
      rate=xk(381)*y(42)*y(20)*xnH
      r_f(381,42)=-rate
      r_f(381,20)=-rate
      r_f(381,27)=rate
      r_f(381,30)=rate
C  382)   OH+   +   CH4   ->   H3O+  +   CH2                  
C         (42        22         47        20)
      rate=xk(382)*y(42)*y(22)*xnH
      r_f(382,42)=-rate
      r_f(382,22)=-rate
      r_f(382,47)=rate
      r_f(382,20)=rate
C  383)   OH+   +   CH4   ->   CH5+  +   O                
C         (42        22         29        30)
      rate=xk(383)*y(42)*y(22)*xnH
      r_f(383,42)=-rate
      r_f(383,22)=-rate
      r_f(383,29)=rate
      r_f(383,30)=rate
C  384)   OH+   +   OH    ->   H2O+  +   O               
C         (42        32         44        30)
      rate=xk(384)*y(42)*y(32)*xnH
      r_f(384,42)=-rate
      r_f(384,32)=-rate
      r_f(384,44)=rate
      r_f(384,30)=rate
C  385)   OH+   +   H2O   ->   H2O+  +   OH                  
C         (42        34         44        32)
      rate=xk(385)*y(42)*y(34)*xnH
      r_f(385,42)=-rate
      r_f(385,34)=-rate
      r_f(385,44)=rate
      r_f(385,32)=rate
C  386)   OH+   +   H2O   ->   H3O+  +   O            
C         (42        34         47        30)
      rate=xk(386)*y(42)*y(34)*xnH
      r_f(386,42)=-rate
      r_f(386,34)=-rate
      r_f(386,47)=rate
      r_f(386,30)=rate
C  387)   OH+   +   C2    ->   C2+   +   OH             
C         (42        18         24        32)
      rate=xk(387)*y(42)*y(18)*xnH
      r_f(387,42)=-rate
      r_f(387,18)=-rate
      r_f(387,24)=rate
      r_f(387,32)=rate
C  388)   OH+   +   CO    ->   HCO+  +   O               
C         (42        33         45        30)
      rate=xk(388)*y(42)*y(33)*xnH
      r_f(388,42)=-rate
      r_f(388,33)=-rate
      r_f(388,45)=rate
      r_f(388,30)=rate
C  389)   OH+   +   HCO   ->   H2O+  +   CO                
C         (42        35         44        33)
      rate=xk(389)*y(42)*y(35)*xnH
      r_f(389,42)=-rate
      r_f(389,35)=-rate
      r_f(389,44)=rate
      r_f(389,33)=rate
C  390)   OH+   +   HCO   ->   HCO+  +   OH                  
C         (42        35         45        32)    
      rate=xk(390)*y(42)*y(35)*xnH
      r_f(390,42)=-rate
      r_f(390,35)=-rate
      r_f(390,45)=rate
      r_f(390,32)=rate
C  391)   OH+   +   HCO   ->   H2CO+ +   O                       
C         (42        35         48        30)
      rate=xk(391)*y(42)*y(35)*xnH
      r_f(391,42)=-rate
      r_f(391,35)=-rate
      r_f(391,48)=rate
      r_f(391,30)=rate
C  392)   OH+   +   H2CO  ->   H3CO+ +   O                
C         (42        38         50        30)
      rate=xk(392)*y(42)*y(38)*xnH
      r_f(392,42)=-rate
      r_f(392,38)=-rate
      r_f(392,50)=rate
      r_f(392,30)=rate
C  393)   OH+   +   H2CO  ->   H2CO+ +   OH                
C         (42        38         48        32) 
      rate=xk(393)*y(42)*y(38)*xnH
      r_f(393,42)=-rate
      r_f(393,38)=-rate
      r_f(393,48)=rate
      r_f(393,32)=rate
C  394)   OH+   +   O2    ->   O2+   +   OH               
C         (42        31         41        32)
      rate=xk(394)*y(42)*y(31)*xnH
      r_f(394,42)=-rate
      r_f(394,31)=-rate
      r_f(394,41)=rate
      r_f(394,32)=rate
C  395)   OH+   +   CO2   ->   HCO2+ +   O                            
C         (42        37         49        30)
      rate=xk(395)*y(42)*y(37)*xnH
      r_f(395,42)=-rate
      r_f(395,37)=-rate
      r_f(395,49)=rate
      r_f(395,30)=rate
C  396)   CH5+  +   H     ->   CH4+  +   H2                         
C         (29        1          28        2)
      rate=xk(396)*y(29)*y(1)*xnH
      r_f(396,29)=-rate
      r_f(396,1)=-rate
      r_f(396,28)=rate
      r_f(396,2)=rate
C  397)   CH5+  +   C     ->   CH+   +   CH4                             
C         (29        17         25        22)
      rate=xk(397)*y(29)*y(17)*xnH
      r_f(397,29)=-rate
      r_f(397,17)=-rate
      r_f(397,25)=rate
      r_f(397,22)=rate
C  398)   CH5+  +   O     ->   H3O+  +   CH2                           
C         (29        30         47        20)
      rate=xk(398)*y(29)*y(30)*xnH
      r_f(398,29)=-rate
      r_f(398,30)=-rate
      r_f(398,47)=rate
      r_f(398,20)=rate
C  399)   CH5+  +   O     ->   H3CO+ +   H2                         
C         (29        30         50        2)
      rate=xk(399)*y(29)*y(30)*xnH
      r_f(399,29)=-rate
      r_f(399,30)=-rate
      r_f(399,50)=rate
      r_f(399,2)=rate
C  400)   CH5+  +   CH    ->   CH2+  +   CH4                          
C         (29        19         26        22)
      rate=xk(400)*y(29)*y(19)*xnH
      r_f(400,29)=-rate
      r_f(400,19)=-rate
      r_f(400,26)=rate
      r_f(400,22)=rate
C  401)   CH5+  +   CH2   ->   CH3+  +   CH4             
C         (29        20         27        22)
      rate=xk(401)*y(29)*y(20)*xnH
      r_f(401,29)=-rate
      r_f(401,20)=-rate
      r_f(401,27)=rate
      r_f(401,22)=rate
C  402)   CH5+  +   OH    ->   H2O+  +   CH4              
C         (29        32         44        22)
      rate=xk(402)*y(29)*y(32)*xnH
      r_f(402,29)=-rate
      r_f(402,32)=-rate
      r_f(402,44)=rate
      r_f(402,22)=rate
C  403)   CH5+  +   H2O   ->   H3O+  +   CH4               
C         (29        34         47        22)
      rate=xk(403)*y(29)*y(34)*xnH
      r_f(403,29)=-rate
      r_f(403,34)=-rate
      r_f(403,47)=rate
      r_f(403,22)=rate
C  404)   CH5+  +   CO    ->   HCO+  +   CH4                 
C         (29        33         45        22)
      rate=xk(404)*y(29)*y(33)*xnH
      r_f(404,29)=-rate
      r_f(404,33)=-rate
      r_f(404,45)=rate
      r_f(404,22)=rate
C  405)   CH5+  +   HCO   ->   H2CO+ +   CH4             
C         (29        35         48        22)
      rate=xk(405)*y(29)*y(35)*xnH
      r_f(405,29)=-rate
      r_f(405,35)=-rate
      r_f(405,48)=rate
      r_f(405,22)=rate
C  406)   CH5+  +   H2CO  ->   H3CO+ +   CH4                
C         (29        38         50        22)
      rate=xk(406)*y(29)*y(38)*xnH
      r_f(406,29)=-rate
      r_f(406,38)=-rate
      r_f(406,50)=rate
      r_f(406,22)=rate
C  407)   CH5+  +   CO2   ->   HCO2+ +   CH4               
C         (29        37         49        22)
      rate=xk(407)*y(29)*y(37)*xnH
      r_f(407,29)=-rate
      r_f(407,37)=-rate
      r_f(407,49)=rate
      r_f(407,22)=rate
C  408)   H2O+  +   C     ->   CH+   +   OH               
C         (44        17         25        32)
      rate=xk(408)*y(44)*y(17)*xnH
      r_f(408,44)=-rate
      r_f(408,17)=-rate
      r_f(408,25)=rate
      r_f(408,32)=rate
C  409)   H2O+  +   O     ->   O2+   +   H2                 
C         (44        30         41        2)
      rate=xk(409)*y(44)*y(30)*xnH
      r_f(409,44)=-rate
      r_f(409,30)=-rate
      r_f(409,41)=rate
      r_f(409,2)=rate
C  410)   H2O+  +   H2    ->   H3O+  +   H               
C         (44        2          47        1)
      rate=xk(410)*y(44)*y(2)*xnH
      r_f(410,44)=-rate
      r_f(410,2)=-rate
      r_f(410,47)=rate
      r_f(410,1)=rate
C  411)   H2O+  +   CH    ->   CH+   +   H2O               
C         (44        19         25        34)
      rate=xk(411)*y(44)*y(19)*xnH
      r_f(411,44)=-rate
      r_f(411,19)=-rate
      r_f(411,25)=rate
      r_f(411,34)=rate
C  412)   H2O+  +   CH    ->   CH2+  +   OH                  
C         (44        19         26        32)
      rate=xk(412)*y(44)*y(19)*xnH
      r_f(412,44)=-rate
      r_f(412,19)=-rate
      r_f(412,26)=rate
      r_f(412,32)=rate
C  413)   H2O+  +   CH2   ->   CH3+  +   OH                
C         (44        20         27        32)
      rate=xk(413)*y(44)*y(20)*xnH
      r_f(413,44)=-rate
      r_f(413,20)=-rate
      r_f(413,27)=rate
      r_f(413,32)=rate
C  414)   H2O+  +   CH2   ->   CH2+  +   H2O                 
C         (44        20         26        34)
      rate=xk(414)*y(44)*y(20)*xnH
      r_f(414,44)=-rate
      r_f(414,20)=-rate
      r_f(414,26)=rate
      r_f(414,34)=rate
C  415)   H2O+  +   CH4   ->   H3O+  +   CH3              
C         (44        22         47        21) 
      rate=xk(415)*y(44)*y(22)*xnH
      r_f(415,44)=-rate
      r_f(415,22)=-rate
      r_f(415,47)=rate
      r_f(415,21)=rate
C  416)   H2O+  +   OH    ->   H3O+  +   O               
C         (44        32         47        30)
      rate=xk(416)*y(44)*y(32)*xnH
      r_f(416,44)=-rate
      r_f(416,32)=-rate
      r_f(416,47)=rate
      r_f(416,30)=rate
C  417)   H2O+  +   H2O   ->   H3O+  +   OH                
C         (44        34         47        32)
      rate=xk(417)*y(44)*y(34)*xnH
      r_f(417,44)=-rate
      r_f(417,34)=-rate
      r_f(417,47)=rate
      r_f(417,32)=rate
C  418)   H2O+  +   C2    ->   C2+   +   H2O                 
C         (44        18         24        34)
      rate=xk(418)*y(44)*y(18)*xnH
      r_f(418,44)=-rate
      r_f(418,18)=-rate
      r_f(418,24)=rate
      r_f(418,34)=rate
C  419)   H2O+  +   CO    ->   HCO+  +   OH                
C         (44        33         45        32)
      rate=xk(419)*y(44)*y(33)*xnH
      r_f(419,44)=-rate
      r_f(419,33)=-rate
      r_f(419,45)=rate
      r_f(419,32)=rate
C  420)   H2O+  +   HCO   ->   H3O+  +   CO               
C         (44        35         47        33)
      rate=xk(420)*y(44)*y(35)*xnH
      r_f(420,44)=-rate
      r_f(420,35)=-rate
      r_f(420,47)=rate
      r_f(420,33)=rate
C  421)   H2O+  +   HCO   ->   HCO+  +   H2O              
C         (44        35         45        34)
      rate=xk(421)*y(44)*y(35)*xnH
      r_f(421,44)=-rate
      r_f(421,35)=-rate
      r_f(421,45)=rate
      r_f(421,34)=rate
C  422)   H2O+  +   HCO   ->   H2CO+ +   OH               
C         (44        35         48        32)
      rate=xk(422)*y(44)*y(35)*xnH
      r_f(422,44)=-rate
      r_f(422,35)=-rate
      r_f(422,48)=rate
      r_f(422,32)=rate
C  423)   H2O+  +   H2CO  ->   H2CO+ +   H2O                    
C         (44        38         48        34)
      rate=xk(423)*y(44)*y(38)*xnH
      r_f(423,44)=-rate
      r_f(423,38)=-rate
      r_f(423,48)=rate
      r_f(423,34)=rate
C  424)   H2O+  +   H2CO  ->   H3CO+ +   OH                
C         (44        38         50        32)
      rate=xk(424)*y(44)*y(38)*xnH
      r_f(424,44)=-rate
      r_f(424,38)=-rate
      r_f(424,50)=rate
      r_f(424,32)=rate
C  425)   H2O+  +   O2    ->   O2+   +   H2O                         
C         (44        31         41        34)
      rate=xk(425)*y(44)*y(31)*xnH
      r_f(425,44)=-rate
      r_f(425,31)=-rate
      r_f(425,41)=rate
      r_f(425,34)=rate
C  426)   H3O+  +   C     ->   HCO+  +   H2                
C         (47        17         45        2)
      rate=xk(426)*y(47)*y(17)*xnH
      r_f(426,47)=-rate
      r_f(426,17)=-rate
      r_f(426,45)=rate
      r_f(426,2)=rate
C  427)   H3O+  +   H-    ->   OH    +   H2    +   H             
C         (47        7          32        2         1)
      rate=xk(427)*y(47)*y(7)*xnH
      r_f(427,47)=-rate
      r_f(427,7)=-rate
      r_f(427,32)=rate
      r_f(427,2)=rate
      r_f(427,1)=rate
C  428)   H3O+  +   H-    ->   H2O   +   H2                      
C         (47        7          34        2)
      rate=xk(428)*y(47)*y(7)*xnH
      r_f(428,47)=-rate
      r_f(428,7)=-rate
      r_f(428,34)=rate
      r_f(428,2)=rate
C  429)   H3O+  +   CH    ->   CH2+  +   H2O               
C         (47        19         26        34)
      rate=xk(429)*y(47)*y(19)*xnH
      r_f(429,47)=-rate
      r_f(429,19)=-rate
      r_f(429,26)=rate
      r_f(429,34)=rate
C  430)   H3O+  +   CH2   ->   CH3+  +   H2O                 
C         (47        20         27        34)
      rate=xk(430)*y(47)*y(20)*xnH
      r_f(430,47)=-rate
      r_f(430,20)=-rate
      r_f(430,27)=rate
      r_f(430,34)=rate
C  431)   H3O+  +   H2CO  ->   H3CO+ +   H2O            
C         (47        38         50        34)
      rate=xk(431)*y(47)*y(38)*xnH
      r_f(431,47)=-rate
      r_f(431,38)=-rate
      r_f(431,50)=rate
      r_f(431,34)=rate
C  432)   C2+   +   C     ->   C+    +   C2                                 
C         (24        17         23        18) 
      rate=xk(432)*y(24)*y(17)*xnH
      r_f(432,24)=-rate
      r_f(432,17)=-rate
      r_f(432,23)=rate
      r_f(432,18)=rate
C  433)   C2+   +   O     ->   CO+   +   C                 
C         (24        30         43        17)
      rate=xk(433)*y(24)*y(30)*xnH
      r_f(433,24)=-rate
      r_f(433,30)=-rate
      r_f(433,43)=rate
      r_f(433,17)=rate
C  434)   C2+   +   CH    ->   CH+   +   C2              
C         (24        19         25        18)
      rate=xk(434)*y(24)*y(19)*xnH
      r_f(434,24)=-rate
      r_f(434,19)=-rate
      r_f(434,25)=rate
      r_f(434,18)=rate
C  435)   C2+   +   CH2   ->   CH2+  +   C2            
C         (24        20         26        18)
      rate=xk(435)*y(24)*y(20)*xnH
      r_f(435,24)=-rate
      r_f(435,20)=-rate
      r_f(435,26)=rate
      r_f(435,18)=rate
C  436)   C2+   +   OH    ->   OH+   +   C2                
C         (24        32         42        18)
      rate=xk(436)*y(24)*y(32)*xnH
      r_f(436,24)=-rate
      r_f(436,32)=-rate
      r_f(436,42)=rate
      r_f(436,18)=rate
C  437)   C2+   +   HCO   ->   HCO+  +   C2            
C         (24        35         45        18)
      rate=xk(437)*y(24)*y(35)*xnH
      r_f(437,24)=-rate
      r_f(437,35)=-rate
      r_f(437,45)=rate
      r_f(437,18)=rate
C  438)   C2+   +   O2    ->   CO+   +   CO                
C         (24        31         43        33)
      rate=xk(438)*y(24)*y(31)*xnH
      r_f(438,24)=-rate
      r_f(438,31)=-rate
      r_f(438,43)=rate
      r_f(438,33)=rate
C  439)   CO+   +   H     ->   H+    +   CO             
C         (43        1          4         33)
      rate=xk(439)*y(43)*y(1)*xnH
      r_f(439,43)=-rate
      r_f(439,1)=-rate
      r_f(439,4)=rate
      r_f(439,33)=rate
C  440)   CO+   +   C     ->   C+    +   CO            
C         (43        17         23        33)
      rate=xk(440)*y(43)*y(17)*xnH
      r_f(440,43)=-rate
      r_f(440,17)=-rate
      r_f(440,23)=rate
      r_f(440,33)=rate
C  441)   CO+   +   O     ->   O+    +   CO                 
C         (43        30         40        33)
      rate=xk(441)*y(43)*y(30)*xnH
      r_f(441,43)=-rate
      r_f(441,30)=-rate
      r_f(441,40)=rate
      r_f(441,33)=rate
C  442)   CO+   +   H2    ->   HCO+  +   H                   
C         (43        2          45        1)
      rate=xk(442)*y(43)*y(2)*xnH
      r_f(442,43)=-rate
      r_f(442,2)=-rate
      r_f(442,45)=rate
      r_f(442,1)=rate
C  443)   CO+   +   CH    ->   CH+   +   CO                  
C         (43        19         25        33)
      rate=xk(443)*y(43)*y(19)*xnH
      r_f(443,43)=-rate
      r_f(443,19)=-rate
      r_f(443,25)=rate
      r_f(443,33)=rate
C  444)   CO+   +   CH    ->   HCO+  +   C              
C         (43        19         45        17)
      rate=xk(444)*y(43)*y(19)*xnH
      r_f(444,43)=-rate
      r_f(444,19)=-rate
      r_f(444,45)=rate
      r_f(444,17)=rate
C  445)   CO+   +   CH2   ->   CH2+  +   CO               
C         (43        20         26        33)
      rate=xk(445)*y(43)*y(20)*xnH
      r_f(445,43)=-rate
      r_f(445,20)=-rate
      r_f(445,26)=rate
      r_f(445,33)=rate
C  446)   CO+   +   CH2   ->   HCO+  +   CH               
C         (43        20         45        19)
      rate=xk(446)*y(43)*y(20)*xnH
      r_f(446,43)=-rate
      r_f(446,20)=-rate
      r_f(446,45)=rate
      r_f(446,19)=rate
C  447)   CO+   +   CH4   ->   CH4+  +   CO                  
C         (43        22         28        33)
      rate=xk(447)*y(43)*y(22)*xnH
      r_f(447,43)=-rate
      r_f(447,22)=-rate
      r_f(447,28)=rate
      r_f(447,33)=rate
C  448)   CO+   +   CH4   ->   HCO+  +   CH3                  
C         (43        22         45        21)
      rate=xk(448)*y(43)*y(22)*xnH
      r_f(448,43)=-rate
      r_f(448,22)=-rate
      r_f(448,45)=rate
      r_f(448,21)=rate
C  449)   CO+   +   OH    ->   OH+   +   CO                  
C         (43        32         42        33)
      rate=xk(449)*y(43)*y(32)*xnH
      r_f(449,43)=-rate
      r_f(449,32)=-rate
      r_f(449,42)=rate
      r_f(449,33)=rate
C  450)   CO+   +   OH    ->   HCO+  +   O                
C         (43        32         45        30)
      rate=xk(450)*y(43)*y(32)*xnH
      r_f(450,43)=-rate
      r_f(450,32)=-rate
      r_f(450,45)=rate
      r_f(450,30)=rate
C  451)   CO+   +   H2O   ->   H2O+  +   CO              
C         (43        34         44        33) 
      rate=xk(451)*y(43)*y(34)*xnH
      r_f(451,43)=-rate
      r_f(451,34)=-rate
      r_f(451,44)=rate
      r_f(451,33)=rate
C  452)   CO+   +   H2O   ->   HCO+  +   OH                 
C         (43        34         45        32) 
      rate=xk(452)*y(43)*y(34)*xnH
      r_f(452,43)=-rate
      r_f(452,34)=-rate
      r_f(452,45)=rate
      r_f(452,32)=rate
C  453)   CO+   +   C2    ->   C2+   +   CO                 
C         (43        18         24        33)
      rate=xk(453)*y(43)*y(18)*xnH
      r_f(453,43)=-rate
      r_f(453,18)=-rate
      r_f(453,24)=rate
      r_f(453,33)=rate
C  454)   CO+   +   HCO   ->   HCO+  +   CO             
C         (43        35         45        33)
      rate=xk(454)*y(43)*y(35)*xnH
      r_f(454,43)=-rate
      r_f(454,35)=-rate
      r_f(454,45)=rate
      r_f(454,33)=rate
C  455)   CO+   +   H2CO  ->   HCO+  +   HCO                             
C         (43        38         45        35)
      rate=xk(455)*y(43)*y(38)*xnH
      r_f(455,43)=-rate
      r_f(455,38)=-rate
      r_f(455,45)=rate
      r_f(455,35)=rate
C  456)   CO+   +   H2CO  ->   H2CO+ +   CO              
C         (43        38         48        33)
      rate=xk(456)*y(43)*y(38)*xnH
      r_f(456,43)=-rate
      r_f(456,38)=-rate
      r_f(456,48)=rate
      r_f(456,33)=rate
C  457)   CO+   +   O2    ->   O2+   +   CO            
C         (43        31         41        33)
      rate=xk(457)*y(43)*y(31)*xnH
      r_f(457,43)=-rate
      r_f(457,31)=-rate
      r_f(457,41)=rate
      r_f(457,33)=rate
C  458)   HCO+  +   C     ->   CH+   +   CO             
C         (45        17         25        33)    
      rate=xk(458)*y(45)*y(17)*xnH
      r_f(458,45)=-rate
      r_f(458,17)=-rate
      r_f(458,25)=rate
      r_f(458,33)=rate
C  459)   HCO+  +   H-    ->   CO    +   H2          
C         (45        7          33        2)
      rate=xk(459)*y(45)*y(7)*xnH
      r_f(459,45)=-rate
      r_f(459,7)=-rate
      r_f(459,33)=rate
      r_f(459,2)=rate
C  460)   HCO+  +   CH    ->   CH2+  +   CO              
C         (45        19         26        33)
      rate=xk(460)*y(45)*y(19)*xnH
      r_f(460,45)=-rate
      r_f(460,19)=-rate
      r_f(460,26)=rate
      r_f(460,33)=rate
C  461)   HCO+  +   CH2   ->   CH3+  +   CO               
C         (45        20         27        33)
      rate=xk(461)*y(45)*y(20)*xnH
      r_f(461,45)=-rate
      r_f(461,20)=-rate
      r_f(461,27)=rate
      r_f(461,33)=rate
C  462)   HCO+  +   OH    ->   H2O+  +   CO                   
C         (45        32         44        33)
      rate=xk(462)*y(45)*y(32)*xnH
      r_f(462,45)=-rate
      r_f(462,32)=-rate
      r_f(462,44)=rate
      r_f(462,33)=rate
C  463)   HCO+  +   OH    ->   HCO2+ +   H                
C         (45        32         49        1)
      rate=xk(463)*y(45)*y(32)*xnH
      r_f(463,45)=-rate
      r_f(463,32)=-rate
      r_f(463,49)=rate
      r_f(463,1)=rate
C  464)   HCO+  +   H2O   ->   H3O+  +   CO               
C         (45        34         47        33)
      rate=xk(464)*y(45)*y(34)*xnH
      r_f(464,45)=-rate
      r_f(464,34)=-rate
      r_f(464,47)=rate
      r_f(464,33)=rate
C  465)   HCO+  +   HCO   ->   H2CO+ +   CO          
C         (45        35         48        33)
      rate=xk(465)*y(45)*y(35)*xnH
      r_f(465,45)=-rate
      r_f(465,35)=-rate
      r_f(465,48)=rate
      r_f(465,33)=rate
C  466)   HCO+  +   H2CO  ->   H3CO+ +   CO                    
C         (45        38         50        33)
      rate=xk(466)*y(45)*y(38)*xnH
      r_f(466,45)=-rate
      r_f(466,38)=-rate
      r_f(466,50)=rate
      r_f(466,33)=rate
C  467)   H2CO+ +   CH    ->   CH+   +   H2CO                            
C         (48        19         25        38)
      rate=xk(467)*y(48)*y(19)*xnH
      r_f(467,48)=-rate
      r_f(467,19)=-rate
      r_f(467,25)=rate
      r_f(467,38)=rate
C  468)   H2CO+ +   CH    ->   CH2+  +   HCO                 
C         (48        19         26        35)
      rate=xk(468)*y(48)*y(19)*xnH
      r_f(468,48)=-rate
      r_f(468,19)=-rate
      r_f(468,26)=rate
      r_f(468,35)=rate
C  469)   H2CO+ +   CH2   ->   CH3+  +   HCO                 
C         (48        20         27        35)
      rate=xk(469)*y(48)*y(20)*xnH
      r_f(469,48)=-rate
      r_f(469,20)=-rate
      r_f(469,27)=rate
      r_f(469,35)=rate
C  470)   H2CO+ +   CH2   ->   CH2+  +   H2CO             
C         (48        20         26        38)
      rate=xk(470)*y(48)*y(20)*xnH
      r_f(470,48)=-rate
      r_f(470,20)=-rate
      r_f(470,26)=rate
      r_f(470,38)=rate
C  471)   H2CO+ +   CH4   ->   H3CO+ +   CH3                     
C         (48        22         50        21)
      rate=xk(471)*y(48)*y(22)*xnH
      r_f(471,48)=-rate
      r_f(471,22)=-rate
      r_f(471,50)=rate
      r_f(471,21)=rate
C  472)   H2CO+ +   H2O   ->   H3O+  +   HCO              
C         (48        34         47        35)
      rate=xk(472)*y(48)*y(34)*xnH
      r_f(472,48)=-rate
      r_f(472,34)=-rate
      r_f(472,47)=rate
      r_f(472,35)=rate
C  473)   H2CO+ +   HCO   ->   HCO+  +   H2CO              
C         (48        35         45        38)
      rate=xk(473)*y(48)*y(35)*xnH
      r_f(473,48)=-rate
      r_f(473,35)=-rate
      r_f(473,45)=rate
      r_f(473,38)=rate
C  474)   H2CO+ +   HCO   ->   H3CO+ +   CO       
C         (48        35         50        33)
      rate=xk(474)*y(48)*y(35)*xnH
      r_f(474,48)=-rate
      r_f(474,35)=-rate
      r_f(474,50)=rate
      r_f(474,33)=rate
C  475)   H2CO+ +   H2CO  ->   H3CO+ +   HCO                   
C         (48        38         50        35)
      rate=xk(475)*y(48)*y(38)*xnH
      r_f(475,48)=-rate
      r_f(475,38)=-rate
      r_f(475,50)=rate
      r_f(475,35)=rate
C  476)   H2CO+ +   O2    ->   HCO+  +   O2H                  
C         (48        31         45        36)
      rate=xk(476)*y(48)*y(31)*xnH
      r_f(476,48)=-rate
      r_f(476,31)=-rate
      r_f(476,45)=rate
      r_f(476,36)=rate
C  477)   H3CO+ +   CH    ->   CH2+  +   H2CO           
C         (50        19         26        38)
      rate=xk(477)*y(50)*y(19)*xnH
      r_f(477,50)=-rate
      r_f(477,19)=-rate
      r_f(477,26)=rate
      r_f(477,38)=rate
C  478)   H3CO+ +   H2O   ->   H3O+  +   H2CO                 
C         (50        34         47        38)
      rate=xk(478)*y(50)*y(34)*xnH
      r_f(478,50)=-rate
      r_f(478,34)=-rate
      r_f(478,47)=rate
      r_f(478,38)=rate
C  479)   O2+   +   C     ->   CO+   +   O                       
C         (41        17         43        30)
      rate=xk(479)*y(41)*y(17)*xnH
      r_f(479,41)=-rate
      r_f(479,17)=-rate
      r_f(479,43)=rate
      r_f(479,30)=rate
C  480)   O2+   +   C     ->   C+    +   O2                    
C         (41        17         23        31)
      rate=xk(480)*y(41)*y(17)*xnH
      r_f(480,41)=-rate
      r_f(480,17)=-rate
      r_f(480,23)=rate
      r_f(480,31)=rate
C  481)   O2+   +   CH    ->   CH+   +   O2                 
C         (41        19         25        31)
      rate=xk(481)*y(41)*y(19)*xnH
      r_f(481,41)=-rate
      r_f(481,19)=-rate
      r_f(481,25)=rate
      r_f(481,31)=rate
C  482)   O2+   +   CH    ->   HCO+  +   O                  
C         (41        19         45        30)
      rate=xk(482)*y(41)*y(19)*xnH
      r_f(482,41)=-rate
      r_f(482,19)=-rate
      r_f(482,45)=rate
      r_f(482,30)=rate
C  483)   O2+   +   CH2   ->   H2CO+ +   O                   
C         (41        20         48        30)
      rate=xk(483)*y(41)*y(20)*xnH
      r_f(483,41)=-rate
      r_f(483,20)=-rate
      r_f(483,48)=rate
      r_f(483,30)=rate
C  484)   O2+   +   CH2   ->   CH2+  +   O2                  
C         (41        20         26        31)
      rate=xk(484)*y(41)*y(20)*xnH
      r_f(484,41)=-rate
      r_f(484,20)=-rate
      r_f(484,26)=rate
      r_f(484,31)=rate
C  485)   O2+   +   C2    ->   CO+   +   CO                      
C         (41        18         43        33)
      rate=xk(485)*y(41)*y(18)*xnH
      r_f(485,41)=-rate
      r_f(485,18)=-rate
      r_f(485,43)=rate
      r_f(485,33)=rate
C  486)   O2+   +   C2    ->   C2+   +   O2                  
C         (41        18         24        31)
      rate=xk(486)*y(41)*y(18)*xnH
      r_f(486,41)=-rate
      r_f(486,18)=-rate
      r_f(486,24)=rate
      r_f(486,31)=rate
C  487)   O2+   +   HCO   ->   O2H+  +   CO               
C         (41        35         46        33)
      rate=xk(487)*y(41)*y(35)*xnH
      r_f(487,41)=-rate
      r_f(487,35)=-rate
      r_f(487,46)=rate
      r_f(487,33)=rate
C  488)   O2+   +   HCO   ->   HCO+  +   O2                 
C         (41        35         45        31)
      rate=xk(488)*y(41)*y(35)*xnH
      r_f(488,41)=-rate
      r_f(488,35)=-rate
      r_f(488,45)=rate
      r_f(488,31)=rate
C  489)   O2+   +   H2CO  ->   HCO+  +   O2    +   H                
C         (41        38         45        31        1)
      rate=xk(489)*y(41)*y(38)*xnH
      r_f(489,41)=-rate
      r_f(489,38)=-rate
      r_f(489,45)=rate
      r_f(489,31)=rate
      r_f(489,1)=rate
C  490)   O2+   +   H2CO  ->   H2CO+ +   O2                   
C         (41        38         48        31)
      rate=xk(490)*y(41)*y(38)*xnH
      r_f(490,41)=-rate
      r_f(490,38)=-rate
      r_f(490,48)=rate
      r_f(490,31)=rate
C  491)   O2H+  +   C     ->   CH+   +   O2                
C         (46        17         25        31)
      rate=xk(491)*y(46)*y(17)*xnH
      r_f(491,46)=-rate
      r_f(491,17)=-rate
      r_f(491,25)=rate
      r_f(491,31)=rate
C  492)   O2H+  +   O     ->   OH+   +   O2             
C         (46        30         42        31)
      rate=xk(492)*y(46)*y(30)*xnH
      r_f(492,46)=-rate
      r_f(492,30)=-rate
      r_f(492,42)=rate
      r_f(492,31)=rate
C  493)   O2H+  +   H2    ->   H3+   +   O2                   
C         (46        2          6         31)   
      rate=xk(493)*y(46)*y(2)*xnH
      r_f(493,46)=-rate
      r_f(493,2)=-rate
      r_f(493,6)=rate
      r_f(493,31)=rate          
C  494)   O2H+  +   CH    ->   CH2+  +   O2              
C         (46        19         26        31)    
      rate=xk(494)*y(46)*y(19)*xnH
      r_f(494,46)=-rate
      r_f(494,19)=-rate
      r_f(494,26)=rate
      r_f(494,31)=rate         
C  495)   O2H+  +   CH2   ->   CH3+  +   O2               
C         (46        20         27        31)     
      rate=xk(495)*y(46)*y(20)*xnH
      r_f(495,46)=-rate
      r_f(495,20)=-rate
      r_f(495,27)=rate
      r_f(495,31)=rate            
C  496)   O2H+  +   OH    ->   H2O+  +   O2                  
C         (46        32         44        31)     
      rate=xk(496)*y(46)*y(32)*xnH
      r_f(496,46)=-rate
      r_f(496,32)=-rate
      r_f(496,44)=rate
      r_f(496,31)=rate           
C  497)   O2H+  +   H2O   ->   H3O+  +   O2                
C         (46        34         47        31)     
      rate=xk(497)*y(46)*y(34)*xnH
      r_f(497,46)=-rate
      r_f(497,34)=-rate
      r_f(497,47)=rate
      r_f(497,31)=rate          
C  498)   O2H+  +   CO    ->   HCO+  +   O2               
C         (46        33         45        31)      
      rate=xk(498)*y(46)*y(33)*xnH
      r_f(498,46)=-rate
      r_f(498,33)=-rate
      r_f(498,45)=rate
      r_f(498,31)=rate             
C  499)   O2H+  +   HCO   ->   H2CO+ +   O2                
C         (46        35         48        31)      
      rate=xk(499)*y(46)*y(35)*xnH
      r_f(499,46)=-rate
      r_f(499,35)=-rate
      r_f(499,48)=rate
      r_f(499,31)=rate             
C  500)   O2H+  +   H2CO  ->   H3CO+ +   O2                   
C         (46        38         50        31)   
      rate=xk(500)*y(46)*y(38)*xnH
      r_f(500,46)=-rate
      r_f(500,38)=-rate
      r_f(500,50)=rate
      r_f(500,31)=rate
C  501)   O2H+  +   CO2   ->   HCO2+ +   O2                   
C         (46        37         49        31)
      rate=xk(501)*y(46)*y(37)*xnH
      r_f(501,46)=-rate
      r_f(501,37)=-rate
      r_f(501,49)=rate
      r_f(501,31)=rate
C  502)   HCO2+ +   C     ->   CH+   +   CO2                         
C         (49        17         25        37)
      rate=xk(502)*y(49)*y(17)*xnH
      r_f(502,49)=-rate
      r_f(502,17)=-rate
      r_f(502,25)=rate
      r_f(502,37)=rate
C  503)   HCO2+ +   O     ->   HCO+  +   O2                     
C         (49        30         45        31)
      rate=xk(503)*y(49)*y(30)*xnH
      r_f(503,49)=-rate
      r_f(503,30)=-rate
      r_f(503,45)=rate
      r_f(503,31)=rate
C  504)   HCO2+ +   CH4   ->   CH5+  +   CO2                     
C         (49        22         29        37)
      rate=xk(504)*y(49)*y(22)*xnH
      r_f(504,49)=-rate
      r_f(504,22)=-rate
      r_f(504,29)=rate
      r_f(504,37)=rate
C  505)   HCO2+ +   H2O   ->   H3O+  +   CO2                 
C         (49        34         47        37)
      rate=xk(505)*y(49)*y(34)*xnH
      r_f(505,49)=-rate
      r_f(505,34)=-rate
      r_f(505,47)=rate
      r_f(505,37)=rate
C  506)   HCO2+ +   CO    ->   HCO+  +   CO2                  
C         (49        33         45        37)
      rate=xk(506)*y(49)*y(33)*xnH
      r_f(506,49)=-rate
      r_f(506,33)=-rate
      r_f(506,45)=rate
      r_f(506,37)=rate
C  507)   C+    +   e     ->   C     +   ph.              
C         (23        3          17)
      rate=xk(507)*y(23)*y(3)*xnH
      r_f(507,23)=-rate
      r_f(507,3)=-rate
      r_f(507,17)=rate
C  508)   CH+   +   e     ->   C     +   H                
C         (25        3          17        1)
      rate=xk(508)*y(25)*y(3)*xnH
      r_f(508,25)=-rate
      r_f(508,3)=-rate
      r_f(508,17)=rate
      r_f(508,1)=rate
C  509)   CH2+  +   e     ->   C     +   H2              
C         (26        3          17        2)
      rate=xk(509)*y(26)*y(3)*xnH
      r_f(509,26)=-rate
      r_f(509,3)=-rate
      r_f(509,17)=rate
      r_f(509,2)=rate
C  510)   CH2+  +   e     ->   CH    +   H                       
C         (26        3          19        1)
      rate=xk(510)*y(26)*y(3)*xnH
      r_f(510,26)=-rate
      r_f(510,3)=-rate
      r_f(510,19)=rate
      r_f(510,1)=rate
C  511)   CH3+  +   e     ->   CH2   +   H                         
C         (27        3          20        1)
      rate=xk(511)*y(27)*y(3)*xnH
      r_f(511,27)=-rate
      r_f(511,3)=-rate
      r_f(511,20)=rate
      r_f(511,1)=rate
C  512)   CH3+  +   e     ->   CH    +   H2                  
C         (27        3          19        2)
      rate=xk(512)*y(27)*y(3)*xnH
      r_f(512,27)=-rate
      r_f(512,3)=-rate
      r_f(512,19)=rate
      r_f(512,2)=rate
C  513)   CH3+  +   e     ->   CH    + 2 H               
C         (27        3          19      2*1)
      rate=xk(513)*y(27)*y(3)*xnH
      r_f(513,27)=-rate
      r_f(513,3)=-rate
      r_f(513,19)=rate
      r_f(513,1)=2.d0*rate
C  514)   CH3+  +   e     ->   CH3   +   ph.              
C         (27        3          21)
      rate=xk(514)*y(27)*y(3)*xnH
      r_f(514,27)=-rate
      r_f(514,3)=-rate
      r_f(514,21)=rate
C  515)   O+    +   e     ->   O     +   ph.               
C         (40        3          30)
      rate=xk(515)*y(40)*y(3)*xnH
      r_f(515,40)=-rate
      r_f(515,3)=-rate
      r_f(515,30)=rate
C  516)   CH4+  +   e     ->   CH3   +   H                
C         (28        3          21        1)
      rate=xk(516)*y(28)*y(3)*xnH
      r_f(516,28)=-rate
      r_f(516,3)=-rate
      r_f(516,21)=rate
      r_f(516,1)=rate
C  517)   CH4+  +   e     ->   CH2   + 2 H                     
C         (28        3          20      2*1)
      rate=xk(517)*y(28)*y(3)*xnH
      r_f(517,28)=-rate
      r_f(517,3)=-rate
      r_f(517,20)=rate
      r_f(517,1)=2.d0*rate
C  518)   OH+   +   e     ->   O     +   H                  
C         (42        3          30        1)
      rate=xk(518)*y(42)*y(3)*xnH
      r_f(518,42)=-rate
      r_f(518,3)=-rate
      r_f(518,30)=rate
      r_f(518,1)=rate
C  519)   CH5+  +   e     ->   CH4   +   H                
C         (29        3          22        1)
      rate=xk(519)*y(29)*y(3)*xnH
      r_f(519,29)=-rate
      r_f(519,3)=-rate
      r_f(519,22)=rate
      r_f(519,1)=rate
C  520)   CH5+  +   e     ->   CH3   +   H2                     
C         (29        3          21        2)
      rate=xk(520)*y(29)*y(3)*xnH
      r_f(520,29)=-rate
      r_f(520,3)=-rate
      r_f(520,21)=rate
      r_f(520,2)=rate
C  521)   H2O+  +   e     ->   OH    +   H                    
C         (44        3          32        1)
      rate=xk(521)*y(44)*y(3)*xnH
      r_f(521,44)=-rate
      r_f(521,3)=-rate
      r_f(521,32)=rate
      r_f(521,1)=rate
C  522)   H2O+  +   e     ->   O     +   H2                    
C         (44        3          30        2)
      rate=xk(522)*y(44)*y(3)*xnH
      r_f(522,44)=-rate
      r_f(522,3)=-rate
      r_f(522,30)=rate
      r_f(522,2)=rate
C  523)   H3O+  +   e     ->   H2O   +   H                      
C         (47        3          34        1)
      rate=xk(523)*y(47)*y(3)*xnH
      r_f(523,47)=-rate
      r_f(523,3)=-rate
      r_f(523,34)=rate
      r_f(523,1)=rate
C  524)   H3O+  +   e     ->   OH    + 2 H                  
C         (47        3          32      2*1)
      rate=xk(524)*y(47)*y(3)*xnH
      r_f(524,47)=-rate
      r_f(524,3)=-rate
      r_f(524,32)=rate
      r_f(524,1)=2.d0*rate
C  525)   C2+   +   e     -> 2 C                            
C         (24        3        2*17)  
      rate=xk(525)*y(24)*y(3)*xnH
      r_f(525,24)=-rate
      r_f(525,3)=-rate
      r_f(525,17)=2.d0*rate
C  526)   CO+   +   e     ->   O     +   C                     
C         (43        3          30        17)
      rate=xk(526)*y(43)*y(3)*xnH
      r_f(526,43)=-rate
      r_f(526,3)=-rate
      r_f(526,30)=rate
      r_f(526,17)=rate
C  527)   HCO+  +   e     ->   CO    +   H                   
C         (45        3          33        1) 
      rate=xk(527)*y(45)*y(3)*xnH
      r_f(527,45)=-rate
      r_f(527,3)=-rate
      r_f(527,33)=rate
      r_f(527,1)=rate
C  528)   H2CO+ +   e     ->   HCO   +   H                             
C         (48        3          35        1)
      rate=xk(528)*y(48)*y(3)*xnH
      r_f(528,48)=-rate
      r_f(528,3)=-rate
      r_f(528,35)=rate
      r_f(528,1)=rate
C  529)   H2CO+ +   e     ->   CO    + 2 H                             
C         (48        3          33      2*1)
      rate=xk(529)*y(48)*y(3)*xnH
      r_f(529,48)=-rate
      r_f(529,3)=-rate
      r_f(529,33)=rate
      r_f(529,1)=2.d0*rate
C  530)   H2CO+ +   e     ->   H2CO  +   ph.                       
C         (48        3          38)
      rate=xk(530)*y(48)*y(3)*xnH
      r_f(530,48)=-rate
      r_f(530,3)=-rate
      r_f(530,38)=rate
C  531)   H3CO+ +   e     ->   CO    +   H     +   H2             
C         (50        3          33        1         2)
      rate=xk(531)*y(50)*y(3)*xnH
      r_f(531,50)=-rate
      r_f(531,3)=-rate
      r_f(531,33)=rate
      r_f(531,1)=rate
      r_f(531,2)=rate
C  532)   H3CO+ +   e     ->   HCO   + 2 H                   
C         (50        3          35      2*1)
      rate=xk(532)*y(50)*y(3)*xnH
      r_f(532,50)=-rate
      r_f(532,3)=-rate
      r_f(532,35)=rate
      r_f(532,1)=2.d0*rate
C  533)   H3CO+ +   e     ->   H2CO  +   H                       
C         (50        3          38        1) 
      rate=xk(533)*y(50)*y(3)*xnH
      r_f(533,50)=-rate
      r_f(533,3)=-rate
      r_f(533,38)=rate
      r_f(533,1)=rate
C  534)   O2+   +   e     -> 2 O                               
C         (41        3        2*30)
      rate=xk(534)*y(41)*y(3)*xnH
      r_f(534,41)=-rate
      r_f(534,3)=-rate
      r_f(534,30)=2.d0*rate
C  535)   O2H+  +   e     ->   O2    +   H                       
C         (46        3          31        1)
      rate=xk(535)*y(46)*y(3)*xnH
      r_f(535,46)=-rate
      r_f(535,3)=-rate
      r_f(535,31)=rate
      r_f(535,1)=rate
C  536)   HCO2+ +   e     ->   CO2   +   H                      
C         (49        3          37        1)
      rate=xk(536)*y(49)*y(3)*xnH
      r_f(536,49)=-rate
      r_f(536,3)=-rate
      r_f(536,37)=rate
      r_f(536,1)=rate
C  537)   HCO2+ +   e     ->   CO    +   OH                        
C         (49        3          33        32)
      rate=xk(537)*y(49)*y(3)*xnH
      r_f(537,49)=-rate
      r_f(537,3)=-rate
      r_f(537,33)=rate
      r_f(537,32)=rate
C  538)    H2    +   C     ->   CH2   +   ph.
C         (2        17         20)      
      rate=xk(538)*y(2)*y(17)*xnH
      r_f(538,2)=-rate
      r_f(538,17)=-rate
      r_f(538,20)=rate
C  539)    OH    +   CH3   ->   CH4   +   O
C         (32       21         22        30)
      rate=xk(539)*y(32)*y(21)*xnH
      r_f(539,32)=-rate
      r_f(539,21)=-rate
      r_f(539,22)=rate
      r_f(539,30)=rate
C  540)    H+    +   H2CO  ->   CO+   +   H2   +   H
C        (4         38         43        2        1)
      rate=xk(540)*y(4)*y(38)*xnH
      r_f(540,4)=-rate
      r_f(540,38)=-rate
      r_f(540,43)=rate
      r_f(540,2)=rate
      r_f(540,1)=rate
C  541)    He+   +   C     ->   C+    +   He
C        (9         17         23        8)
      rate=xk(541)*y(9)*y(17)*xnH
      r_f(541,9)=-rate
      r_f(541,17)=-rate
      r_f(541,23)=rate
      r_f(541,8)=rate
C  542)    He+   +   H2CO  ->   H2CO+ +   He
C        (9         38         48        8)
      rate=xk(542)*y(9)*y(38)*xnH
      r_f(542,9)=-rate
      r_f(542,38)=-rate
      r_f(542,48)=rate
      r_f(542,8)=rate
C  543)    He+   +   H2CO  ->   CH2+  +   O    +   He
C        (9         38         26        30       8)
      rate=xk(543)*y(9)*y(38)*xnH
      r_f(543,9)=-rate
      r_f(543,38)=-rate
      r_f(543,26)=rate
      r_f(543,30)=rate
      r_f(543,8)=rate
C  544)   H     +   CR    ->   H+    +   e                    
C         (1                    4         3)
      rate=xk(544)*y(1)
      r_f(544,1)=-rate
      r_f(544,4)=rate
      r_f(544,3)=rate
C  545)   He    +   CR    ->   He+   +   e                  
C         (8                    9         3)
      rate=xk(545)*y(8)
      r_f(545,8)=-rate
      r_f(545,9)=rate
      r_f(545,3)=rate
C  546)   C     +   CR    ->   C+    +   e                         
C         (17                   23        3)
      rate=xk(546)*y(17)
      r_f(546,17)=-rate
      r_f(546,23)=rate
      r_f(546,3)=rate
C  547)   O     +   CR    ->   O+    +   e                     
C         (30                   40        3)
      rate=xk(547)*y(30)
      r_f(547,30)=-rate
      r_f(547,40)=rate
      r_f(547,3)=rate
C  548)   H2    +   CR    ->   H+    +   H     +   e               
C         (2                    4         1         3)
      rate=xk(548)*y(2)
      r_f(548,2)=-rate
      r_f(548,4)=rate
      r_f(548,1)=rate
      r_f(548,3)=rate      
C  549)   H2    +   CR    ->   H2+   +   e                   
C         (2                    5         3)
      rate=xk(549)*y(2)
      r_f(549,2)=-rate
      r_f(549,5)=rate
      r_f(549,3)=rate
C  550)   H2    +   CR    -> 2 H                      
C         (2                  2*1)
      rate=xk(550)*y(2)
      r_f(550,2)=-rate
      r_f(550,1)=2.d0*rate
C  551)   H2    +   CR    ->   H+    +   H-               
C         (2                    4         7)
      rate=xk(551)*y(2)
      r_f(551,2)=-rate
      r_f(551,4)=rate
      r_f(551,7)=rate
C  552)   CO    +   CR    ->   CO+   +   e                      
C         (33                   43        3)
      rate=xk(552)*y(33)
      r_f(552,33)=-rate
      r_f(552,43)=rate
      r_f(552,3)=rate
C  553)   C     +   ph.   ->   C+    +   e                   
C         (17                   23        3)
      rate=xk(553)*y(17)
      r_f(553,17)=-rate
      r_f(553,23)=rate
      r_f(553,3)=rate
C  554)   H-    +   ph.   ->   H     +   e             
C         (7                    1         3)
      rate=xk(554)*y(7)
      r_f(554,7)=-rate
      r_f(554,1)=rate
      r_f(554,3)=rate
C  555)   H2+   +   ph.   ->   H+    +   H             
C         (5                    4         1)
      rate=xk(555)*y(5)
      r_f(555,5)=-rate
      r_f(555,4)=rate
      r_f(555,1)=rate
C  556)   H3+   +   ph.   ->   H2+   +   H               
C         (6                    5         1)
      rate=xk(556)*y(6)
      r_f(556,6)=-rate
      r_f(556,5)=rate
      r_f(556,1)=rate
C  557)   H3+   +   ph.   ->   H+    +   H2               
C         (6                    4         2)
      rate=xk(557)*y(6)
      r_f(557,6)=-rate
      r_f(557,4)=rate
      r_f(557,2)=rate
C  558)   CH    +   ph.   ->   CH+   +   e                 
C         (19                   25        3)
      rate=xk(558)*y(19)
      r_f(558,19)=-rate
      r_f(558,25)=rate
      r_f(558,3)=rate
C  559)   CH    +   ph.   ->   C     +   H                         
C         (19                   17        1)
      rate=xk(559)*y(19)
      r_f(559,19)=-rate
      r_f(559,17)=rate
      r_f(559,1)=rate
C  560)   CH+   +   ph.   ->   C+    +   H                      
C         (25                   23        1)
      rate=xk(560)*y(25)
      r_f(560,25)=-rate
      r_f(560,23)=rate
      r_f(560,1)=rate
C  561)   CH2   +   ph.   ->   CH2+  +   e               
C         (20                   26        3)
      rate=xk(561)*y(20)
      r_f(561,20)=-rate
      r_f(561,26)=rate
      r_f(561,3)=rate
C  562)   CH2   +   ph.   ->   CH    +   H              
C         (20                   19        1)
      rate=xk(562)*y(20)
      r_f(562,20)=-rate
      r_f(562,19)=rate
      r_f(562,1)=rate
C  563)   CH2+  +   ph.   ->   CH+   +   H         
C         (26                   25        1)
      rate=xk(563)*y(26)
      r_f(563,26)=-rate
      r_f(563,25)=rate
      r_f(563,1)=rate
C  564)   CH3   +   ph.   ->   CH    +   H2                  
C         (21                   19        2)
      rate=xk(564)*y(21)
      r_f(564,21)=-rate
      r_f(564,19)=rate
      r_f(564,2)=rate
C  565)   CH3   +   ph.   ->   CH2   +   H                     
C         (21                   20        1)
      rate=xk(565)*y(21)
      r_f(565,21)=-rate
      r_f(565,20)=rate
      r_f(565,1)=rate
C  566)   CH3   +   ph.   ->   CH3+  +   e                         
C         (21                   27        3)
      rate=xk(566)*y(21)
      r_f(566,21)=-rate
      r_f(566,27)=rate
      r_f(566,3)=rate
C  567)   CH3+  +   ph.   ->   CH2+  +   H                        
C         (27                   26        1)
      rate=xk(567)*y(27)
      r_f(567,27)=-rate
      r_f(567,26)=rate
      r_f(567,1)=rate
C  568)   CH3+  +   ph.   ->   CH+   +   H2                       
C         (27                   25        2)
      rate=xk(568)*y(27)
      r_f(568,27)=-rate
      r_f(568,25)=rate
      r_f(568,2)=rate
C  569)   CH4   +   ph.   ->   CH3   +   H                     
C         (22                   21        1)
      rate=xk(569)*y(22)
      r_f(569,22)=-rate
      r_f(569,21)=rate
      r_f(569,1)=rate
C  570)   CH4   +   ph.   ->   CH2   +   H2                      
C         (22                   20        2)
      rate=xk(570)*y(22)
      r_f(570,22)=-rate
      r_f(570,20)=rate
      r_f(570,2)=rate
C  571)   CH4   +   ph.   ->   CH    +   H     +   H2                   
C         (22                   19        1         2)
      rate=xk(571)*y(22)
      r_f(571,22)=-rate
      r_f(571,19)=rate
      r_f(571,1)=rate
      r_f(571,2)=rate
C  572)   OH    +   ph.   ->   O     +   H                         
C         (32                   30        1)
      rate=xk(572)*y(32)
      r_f(572,32)=-rate
      r_f(572,30)=rate
      r_f(572,1)=rate
C  573)   OH    +   ph.   ->   OH+   +   e                          
C         (32                   42        3) 
      rate=xk(573)*y(32)
      r_f(573,32)=-rate
      r_f(573,42)=rate
      r_f(573,3)=rate
C  574)   OH+   +   ph.   ->   H+    +   O                      
C         (42                   4         30)
      rate=xk(574)*y(42)
      r_f(574,42)=-rate
      r_f(574,4)=rate
      r_f(574,30)=rate
C  575)   H2O   +   ph.   ->   OH    +   H                            
C         (34                   32        1)
      rate=xk(575)*y(34)
      r_f(575,34)=-rate
      r_f(575,32)=rate
      r_f(575,1)=rate
C  576)   H2O   +   ph.   ->   H2O+  +   e                           
C         (34                   44        3)
      rate=xk(576)*y(34)
      r_f(576,34)=-rate
      r_f(576,44)=rate
      r_f(576,3)=rate
C  577)   C2    +   ph.   -> 2 C                                    
C         (18                 2*17)
      rate=xk(577)*y(18)
      r_f(577,18)=-rate
      r_f(577,17)=2.d0*rate
C  578)   C2    +   ph.   ->   C2+   +   e                   
C         (18                   24        3)
      rate=xk(578)*y(18)
      r_f(578,18)=-rate
      r_f(578,24)=rate
      r_f(578,3)=rate
C  579)   C2+   +   ph.   ->   C+    +   C                 
C         (24                   23        17)
      rate=xk(579)*y(24)
      r_f(579,24)=-rate
      r_f(579,23)=rate
      r_f(579,17)=rate
C  580)   CO    +   ph.   ->   C     +   O                          
C         (33                   17        30)
      rate=xk(580)*y(33)
      r_f(580,33)=-rate
      r_f(580,17)=rate
      r_f(580,30)=rate
C  581)   HCO   +   ph.   ->   H     +   CO                     
C         (35                   1         33)
      rate=xk(581)*y(35)
      r_f(581,35)=-rate
      r_f(581,1)=rate
      r_f(581,33)=rate
C  582)   HCO   +   ph.   ->   HCO+  +   e                    
C         (35                   45        3)
      rate=xk(582)*y(35)
      r_f(582,35)=-rate
      r_f(582,45)=rate
      r_f(582,3)=rate
C  583)   H2CO  +   ph.   ->   CO    +   H2                  
C         (38                   33        2)
      rate=xk(583)*y(38)
      r_f(583,38)=-rate
      r_f(583,33)=rate
      r_f(583,2)=rate
C  584)   H2CO  +   ph.   ->   CO    + 2 H                
C         (38                   33      2*1)
      rate=xk(584)*y(38)
      r_f(584,38)=-rate
      r_f(584,33)=rate
      r_f(584,1)=2.d0*rate
C  585)   H2CO  +   ph.   ->   H2CO+ +   e                
C         (38                   48        3)
      rate=xk(585)*y(38)
      r_f(585,38)=-rate
      r_f(585,48)=rate
      r_f(585,3)=rate
C  586)   H2CO  +   ph.   ->   HCO+  +   e     +   H             
C         (38                   45        3         1)
      rate=xk(586)*y(38)
      r_f(586,38)=-rate
      r_f(586,45)=rate
      r_f(586,3)=rate
      r_f(586,1)=rate
C  587)   O2    +   ph.   -> 2 O                                
C         (31                 2*30)
      rate=xk(587)*y(31)
      r_f(587,31)=-rate
      r_f(587,30)=2.d0*rate
C  588)   O2    +   ph.   ->   O2+   +   e                    
C         (31                   41        3)
      rate=xk(588)*y(31)
      r_f(588,31)=-rate
      r_f(588,41)=rate
      r_f(588,3)=rate
C  589)   CO2   +   ph.   ->   CO    +   O                    
C         (37                   33        30)
      rate=xk(589)*y(37)
      r_f(589,37)=-rate
      r_f(589,33)=rate
      r_f(589,30)=rate
C  590)   H2   +   ph.   -> 2 H
C        (2                   2*1)
      rate=xk(590)*y(2)
      r_f(590,2)=-rate
      r_f(590,1)=2.d0*rate
C 591)   HD   +   ph.   ->    H     +   D
C        (13                  1         12)
      rate=xk(591)*y(13)
      r_f(591,13)=-rate
      r_f(591,1)=rate
      r_f(591,12)=rate
c
C supplemented reactions 
C 601)   H    +   HCO   ->   O   +   CH2
C       (1        35         30      20  )
      rate=xk(601)*y(1)*y(35)*xnH
      r_f(601,1)=-rate
      r_f(601,35)=-rate
      r_f(601,30)=rate
      r_f(601,20)=rate
C 602)   H    +   H2O2  ->  H2O   +   OH
C       (1        39        34        32 )
      rate=xk(602)*y(1)*y(39)*xnH
      r_f(602,1)=-rate
      r_f(602,39)=-rate
      r_f(602,34)=rate
      r_f(602,32)=rate
C 603)   C    +   CH2   ->   2 CH
C       (17       20         2x19         )
      rate=xk(603)*y(17)*y(20)*xnH
      r_f(603,17)=-rate
      r_f(603,20)=-rate
      r_f(603,19)=2.d0*rate
C 604)   CH   +    O2   ->   CO   +   OH
C       (19        31        33       32  )
      rate=xk(604)*y(19)*y(31)*xnH
      r_f(604,19)=-rate
      r_f(604,31)=-rate
      r_f(604,33)=rate
      r_f(604,32)=rate
C 605)   H    +   HCO   ->   O   +   CH2
C       (1        35         30      20  )
      rate=xk(605)*y(1)*y(35)*xnH
      r_f(605,1)=-rate
      r_f(605,35)=-rate
      r_f(605,30)=rate
      r_f(605,20)=rate
C 606)   CH2  +   O2   ->   CO2    +   2 H
C       (20       31        37         2x1 )
      rate=xk(606)*y(20)*y(31)*xnH
      r_f(606,20)=-rate
      r_f(606,31)=-rate
      r_f(606,37)=rate
      r_f(606,1)=2.d0*rate
C 607)   CH2  +   O2   ->   CO2    +   H2
C       (20       31        37         2  )
      rate=xk(607)*y(20)*y(31)*xnH
      r_f(607,20)=-rate
      r_f(607,31)=-rate
      r_f(607,37)=rate
      r_f(607,2)=rate
C 608)   CH2  +   O2   ->   CO     +   H2O
C       (20       31        33         34 )
      rate=xk(608)*y(20)*y(31)*xnH
      r_f(608,20)=-rate
      r_f(608,31)=-rate
      r_f(608,33)=rate
      r_f(608,34)=rate
C 609)   2 CH3         ->   CH4    +   CH2
C       (2x21               22         20 )
      rate=xk(609)*(y(21)**2)*xnH
      r_f(609,21)=-2.d0*rate
      r_f(609,22)=rate
      r_f(609,20)=rate
C 610)   CH3  +   O    ->   CO     +   H2    +   H
C       (21       30        33         2         1 )
      rate=xk(610)*y(21)*y(30)*xnH
      r_f(610,21)=-rate
      r_f(610,30)=-rate
      r_f(610,33)=rate
      r_f(610,2)=rate
      r_f(610,1)=rate
C 611)   CH3  +   OH   ->   H2CO   +   H2
C       (21       32        38         2 )
      rate=xk(611)*y(21)*y(32)*xnH
      r_f(611,21)=-rate
      r_f(611,32)=-rate
      r_f(611,38)=rate
      r_f(611,2)=rate
C 612)   CH3  +   O2   ->   HCO    +   H2O
C       (21       31        35         34 )
      rate=xk(612)*y(21)*y(31)*xnH
      r_f(612,21)=-rate
      r_f(612,31)=-rate
      r_f(612,35)=rate
      r_f(612,34)=rate
C 613)   O    +  H2CO  ->   CO     +   OH  +   H
C       (30      38         33         32      1 )
      rate=xk(613)*y(30)*y(38)*xnH
      r_f(613,30)=-rate
      r_f(613,38)=-rate
      r_f(613,33)=rate
      r_f(613,32)=rate
      r_f(613,1)=rate
C 614)   C2   +   O2   ->   2 CO
C       (18       31        2x33          )
      rate=xk(614)*y(18)*y(31)*xnH
      r_f(614,18)=-rate
      r_f(614,31)=-rate
      r_f(614,33)=2.d0*rate
C 615)   2 HCO         ->   2 CO   +   H2
C       (2x35               2x33       2   )
      rate=xk(615)*(y(35)**2)*xnH
      r_f(615,35)=-2.d0*rate
      r_f(615,33)=2.d0*rate
      r_f(615,2)=rate
C 616)   HCO   +   O2   ->   CO2   +   OH
C       (35        31        37        32 )
      rate=xk(616)*y(35)*y(31)*xnH
      r_f(616,35)=-rate
      r_f(616,31)=-rate
      r_f(616,37)=rate
      r_f(616,32)=rate
C 617)    O    +   OH   ->   O2    +   H
C        (30       32        31        1  )
      rate=xk(617)*y(30)*y(32)*xnH
      r_f(617,30)=-rate
      r_f(617,32)=-rate
      r_f(617,31)=rate
      r_f(617,1)=rate
C 621)    H3+  +    O    ->  H2O+  +   H
C        (6         30       44        1 )
      rate=xk(621)*y(6)*y(30)*xnH
      r_f(621,6)=-rate
      r_f(621,30)=-rate
      r_f(621,44)=rate
      r_f(621,1)=rate
C 626)   O+    +  CO     ->  CO+   +   O
C       (40       33         43        30  )
      rate=xk(626)*y(40)*y(33)*xnH
      r_f(626,40)=-rate
      r_f(626,33)=-rate
      r_f(626,43)=rate
      r_f(626,30)=rate
C 632)   CH2+  +  e      ->  C     +  2 H
C       (26       3          17       2x1 )
      rate=xk(632)*y(26)*y(3)*xnH
      r_f(632,26)=-rate
      r_f(632,3)=-rate
      r_f(632,17)=rate
      r_f(632,1)=2.d0*rate
C 633)   CH5+  +  e      ->  CH3   +  2 H
C       (29       3          21       2x1 )
      rate=xk(633)*y(29)*y(3)*xnH
      r_f(633,29)=-rate
      r_f(633,3)=-rate
      r_f(633,21)=rate
      r_f(633,1)=2.d0*rate
C 634)   CH5+  +  e      ->  CH2   +    H2    +   H
C       (29       3          20         2         1)
      rate=xk(634)*y(29)*y(3)*xnH
      r_f(634,29)=-rate
      r_f(634,3)=-rate
      r_f(634,20)=rate
      r_f(634,2)=rate
      r_f(634,1)=rate
C 635)   CH5+  +  e      ->  CH    +  2 H2
C       (29       3          19       2x2 )
      rate=xk(635)*y(29)*y(3)*xnH
      r_f(635,29)=-rate
      r_f(635,3)=-rate
      r_f(635,19)=rate
      r_f(635,2)=2.d0*rate
C 636)   H2O+  +  e      ->  O     +  2 H
C       (44       3          30       2x1 )
      rate=xk(636)*y(44)*y(3)*xnH
      r_f(636,44)=-rate
      r_f(636,3)=-rate
      r_f(636,30)=rate
      r_f(636,1)=2.d0*rate
C 637)   H3O+  +  e      ->  O     +    H2   +    H
C       (47       3          30         2         1)
      rate=xk(637)*y(47)*y(3)*xnH
      r_f(637,47)=-rate
      r_f(637,3)=-rate
      r_f(637,30)=rate
      r_f(637,2)=rate
      r_f(637,1)=rate
C 638)   H3O+  +  e      ->  OH    +    H2
C       (47       3          32         2)
      rate=xk(638)*y(47)*y(3)*xnH
      r_f(638,47)=-rate
      r_f(638,3)=-rate
      r_f(638,32)=rate
      r_f(638,2)=rate
C 640)   HCO2+ +  e      ->  CO    +    O   +    H
C       (49       3          33         30       1)
      rate=xk(640)*y(49)*y(3)*xnH
      r_f(640,49)=-rate
      r_f(640,3)=-rate
      r_f(640,33)=rate
      r_f(640,30)=rate
      r_f(640,1)=rate
C 641)   H+    +  He     ->  HeH+  +    ph.
C       (4        8          11       )
      rate=xk(641)*y(4)*y(8)*xnH
      r_f(641,4)=-rate
      r_f(641,8)=-rate
      r_f(641,11)=rate
C 642)   H     +  OH     ->  H2O   +    ph.
C       (1        32         34        )
      rate=xk(642)*y(1)*y(32)*xnH
      r_f(642,1)=-rate
      r_f(642,32)=-rate
      r_f(642,34)=rate
C 643)   H2    +  CH     ->  CH3   +    ph.
C       (2        19         21        )
      rate=xk(643)*y(2)*y(19)*xnH
      r_f(643,2)=-rate
      r_f(643,19)=-rate
      r_f(643,21)=rate
C 644)   C+    +  C      ->  C2+   +    ph.
C       (23       17         24         )
      rate=xk(644)*y(23)*y(17)*xnH
      r_f(644,23)=-rate
      r_f(644,17)=-rate
      r_f(644,24)=rate
C 645)   C     +  O+     ->  CO+   +    ph.
C       (17       40         43         )
      rate=xk(645)*y(17)*y(40)*xnH
      r_f(645,17)=-rate
      r_f(645,40)=-rate
      r_f(645,43)=rate
C 646)   CH2+  +  ph.    ->  CH    +    H+
C       (26                  19         4)
      rate=xk(646)*y(26)
      r_f(646,26)=-rate
      r_f(646,19)=rate
      r_f(646,4)=rate
C 647)   CH2+  +  ph.    ->  C+    +    H2
C       (26                  23         2)
      rate=xk(647)*y(26)
      r_f(647,26)=-rate
      r_f(647,23)=rate
      r_f(647,2)=rate
C 648)   CH4+  +  ph.    ->  CH2+  +    H2
C       (28                  26         2)
      rate=xk(648)*y(28)
      r_f(648,28)=-rate
      r_f(648,26)=rate
      r_f(648,2)=rate
C 649)   CH4+  +  ph.    ->  CH3+  +    H
C       (28                  27         1)
      rate=xk(649)*y(28)
      r_f(649,28)=-rate
      r_f(649,27)=rate
      r_f(649,1)=rate
C 650)   OH+   +  ph.    ->  O+    +    H
C       (42                  40         1)
      rate=xk(650)*y(42)
      r_f(650,42)=-rate
      r_f(650,40)=rate
      r_f(650,1)=rate
C 651)   H2O+   +  ph.   ->  OH+   +    H
C       (44                  42         1)
      rate=xk(651)*y(44)
      r_f(651,44)=-rate
      r_f(651,42)=rate
      r_f(651,1)=rate
C 652)   CO+   +  ph.   ->   C+    +    O
C       (43                  23         30)
      rate=xk(652)*y(43)
      r_f(652,43)=-rate
      r_f(652,23)=rate
      r_f(652,30)=rate
C 653)   HCO+  +  ph.   ->   CO+   +    H
C       (45                  43         1)
      rate=xk(653)*y(45)
      r_f(653,45)=-rate
      r_f(653,43)=rate
      r_f(653,1)=rate
C 654)   O2+   +  ph.   ->   O+    +    O
C       (41                  40         30)
      rate=xk(654)*y(41)
      r_f(654,41)=-rate
      r_f(654,40)=rate
      r_f(654,30)=rate
C 655)   H2O2  +  ph.   -> 2 OH 
C       (39                2x32         )
      rate=xk(655)*y(39)
      r_f(655,39)=-rate
      r_f(655,32)=2.d0*rate
C 656)   C     + CR ph. ->   C+    +    e
C       (17                  23         3)
      rate=xk(656)*y(17)
      r_f(656,17)=-rate
      r_f(656,23)=rate
      r_f(656,3)=rate
C 657)   CH    + CR ph. ->   C     +    H
C       (19                  17         1)
      rate=xk(657)*y(19)
      r_f(657,19)=-rate
      r_f(657,17)=rate
      r_f(657,1)=rate
C 658)   CH+   + CR ph. ->   C+    +    H
C       (25                  23         1)
      rate=xk(658)*y(25)
      r_f(658,25)=-rate
      r_f(658,23)=rate
      r_f(658,1)=rate
C 659)   CH2   + CR ph. ->   CH2+  +    e
C       (20                  26         3)
      rate=xk(659)*y(20)
      r_f(659,20)=-rate
      r_f(659,26)=rate
      r_f(659,3)=rate
C 660)   CH2   + CR ph. ->   CH    +    H
C       (20                  19         1)
      rate=xk(660)*y(20)
      r_f(660,20)=-rate
      r_f(660,19)=rate
      r_f(660,1)=rate
C 661)   CH3   + CR ph. ->   CH3+  +    e
C       (21                  27         3)
      rate=xk(661)*y(21)
      r_f(661,21)=-rate
      r_f(661,27)=rate
      r_f(661,3)=rate
C 662)   CH3   + CR ph. ->   CH2   +    H
C       (21                  20         1)
      rate=xk(662)*y(21)
      r_f(662,21)=-rate
      r_f(662,20)=rate
      r_f(662,1)=rate
C 663)   CH3   + CR ph. ->   CH    +    H2
C       (21                  19         2)
      rate=xk(663)*y(21)
      r_f(663,21)=-rate
      r_f(663,19)=rate
      r_f(663,2)=rate
C 664)   CH4   + CR ph. ->   CH2   +    H2
C       (22                  20         2)
      rate=xk(664)*y(22)
      r_f(664,22)=-rate
      r_f(664,20)=rate
      r_f(664,2)=rate
C 665)   OH    + CR ph. ->   O     +    H
C       (32                  30         1)
      rate=xk(665)*y(32)
      r_f(665,32)=-rate
      r_f(665,30)=rate
      r_f(665,1)=rate
C 666)   H2O   + CR ph. ->   OH    +    H
C       (34                  32         1)
      rate=xk(666)*y(34)
      r_f(666,34)=-rate
      r_f(666,32)=rate
      r_f(666,1)=rate
C 667)   C2    + CR ph. -> 2 C
C       (18                2x17         )
      rate=xk(667)*y(18)
      r_f(667,18)=-rate
      r_f(667,17)=2.d0*rate
C 668)   CO    + CR ph. ->   O     +    C
C       (33                  30         17)
      rate=xk(668)*y(33)
      r_f(668,33)=-rate
      r_f(668,30)=rate
      r_f(668,17)=rate
C 669)  HCO    + CR ph. ->   CO    +    H
C      (35                   33         1)
      rate=xk(669)*y(35)
      r_f(669,35)=-rate
      r_f(669,33)=rate
      r_f(669,1)=rate
C 670)  HCO    + CR ph. ->   HCO+  +    e
C      (35                   45         3)
      rate=xk(670)*y(35)
      r_f(670,35)=-rate
      r_f(670,45)=rate
      r_f(670,3)=rate
C 671)  H2CO   + CR ph. ->   CO    +    H2
C      (38                   33         2)
      rate=xk(671)*y(38)
      r_f(671,38)=-rate
      r_f(671,33)=rate
      r_f(671,2)=rate
C 672)  O2     + CR ph. -> 2 O 
C      (31                 2x30        )
      rate=xk(672)*y(31)
      r_f(672,31)=-rate
      r_f(672,30)=2.d0*rate
C 673)  O2     + CR ph. ->   O2+   +    e
C      (31                   41         3)
      rate=xk(673)*y(31)
      r_f(673,31)=-rate
      r_f(673,41)=rate
      r_f(673,3)=rate
C 674)  H2O2   + CR ph. -> 2 OH
C      (39                 2x32        )
      rate=xk(674)*y(39)
      r_f(674,39)=-rate
      r_f(674,32)=2.d0*rate
C 675)  CO2    + CR ph. ->   CO    +    O 
C      (37                   33         30)
      rate=xk(675)*y(37)
      r_f(675,37)=-rate
      r_f(675,33)=rate
      r_f(675,30)=rate

************************************************************

      do isp=1,N_sp      
         r_f_tot(isp)=0.d0
         do ire=1,N_react
            r_f_tot(isp)=r_f_tot(isp)+r_f(ire,isp)
         enddo
      enddo

      return 
      END


