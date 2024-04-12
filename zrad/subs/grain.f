      SUBROUTINE grtemp(xnH,T,esc_cnt,T_rad,xJ_ex,T_gr)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (yHe=8.333d-2,xk_vis=40.d0)
      DATA xm_p/1.67d-24/     
      funcL(T_d)=4.9d-2*vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
     &     *dsqrt(T*1.d-3)*(0.354+0.5d0*yHe)*(T-T_d)
c      funcL(T_d)=4.9d-2*vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
c     &     *dsqrt(T*1.d-3)*(1.d0-0.8d0*dexp(-75.d0/T))*(T-T_d)
c     for stellar radiation
c      funcD(T_d)=vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)/
c     &     vol_gr((1.d0+4.d0*yHe)*xnH*xm_p,0.d0)
      func(T_d)=xkp_gr((1.d0+4.d0*yHe)*xnH*xm_p,T_d)
     &     *(T_d**4)*esc_cnt
c     collision with gas particles
     &     -xnH*funcL(T_d)
c     stellar radiation
c     &     -4.41d3*xk_vis*funcD(T_d)*xJ_ex
c     CMB
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
      do 11 j=1,100
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
11    continue
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
      END



      FUNCTION xkp_gr(rho,T)
      IMPLICIT REAL*8(a-h,o-z)
C     Planck mean mass absorption coefficient ("opacity") 
C     owing to Dust Grain (solar metalicity)
c     USES vaptemp,linear
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
      END



      FUNCTION vol_gr(rho,T_K)
C     gives grain volume per unit mass of gas (solar metalicity)
      IMPLICIT REAL*8(a-h,o-z)
c     USES vaptemp
      x(t,t1)=dmin1(1.d0,ddim(1.025d0,t/t1)/0.05d0)

      call vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)

      vol_gr=(0.78d-4-0.46d-4*x(T_K,T_tr))*x(T_K,T_ir)
     &     +7.19d-4*x(T_K,T_ol)+2.16d-4*x(T_K,T_pyr)
     &     +1.18d-4*x(T_K,T_tr)+23.53d-4*x(T_K,T_ro)
     &     +6.02d-4*x(T_K,T_vo)+12.93d-4*x(T_K,T_ice)
      
      return
      END



      SUBROUTINE vaptemp(rho,T_ice,T_vo,T_ro,T_tr,T_ir,T_pyr,T_ol)
C     Vaporization Temperatures (K) of Grain Species
C     T_ice   for   Water ice
C     T_vo    for   Volatile organics
C     T_ro    for   Refractory organics
C     T_tr    for   Troilite 
C     T_ir    for   Metallic iron
C     T_pyr   for   Orthopyroxene
C     T_ol    for   Olivine 
      IMPLICIT REAL*8(a-h,o-z)
c     USES linear
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
      
      END


      SUBROUTINE phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)
      IMPLICIT REAL*8(a-h,o-z)
      DATA xm_p/1.67d-24/, T_ro/575.d0/
      if(T_gr_K > T_ro) then
         Gmm_pe=0.d0
         return
      endif
c
      yHe=8.333d-2
      rho=(1.d0+4.d0*yHe)*xm_p*xnH 
      x=G_0*dsqrt(T_K)/(y_e*xnH)
      eps=4.9d-2/(1.d0+4.0d-3*(x**0.73d0))
     &     +3.7d-2*(T_K*1.d-4)**0.7d0/(1.d0+2.0d-4*x)
c
      Gam_pe=1.0d-24*eps*G_0*dexp(-1.8d0*A_v)
c
      beta=0.74d0/(T_K**0.068d0)
      xLd_pe=4.65d-30*(T_K**0.94d0)*(x**beta)*y_e
c
      Gmm_pe=Z_metal*(Gam_pe*xnH-xLd_pe*xnH**2)/rho
      
      return
      END













