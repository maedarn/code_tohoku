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
      Z_gas=1.0d0
      Z_dust=1.0d0
      G_0=1.d0
      iit=15
      itr=10000
c      zeta=1.d-17
c       G_0=100.d0
c      zeta=0.d0
       T_rad=10.d0
c
      yHe=8.333d-2
      yD=3.d-5
c
      open(8,file='data/nT.dat',status='unknown')
      open(11,file='data/1.dat',status='unknown')
      open(12,file='data/2.dat',status='unknown')
c      open(13,file='3.dat',status='unknown')
      open(14,file='data/4.dat',status='unknown')
      open(15,file='data/5.dat',status='unknown')
      open(16,file='data/data_y.dat',status='unknown')
c      open(15,file='data/data_5',status='unknown')
c      open(16,file='data/ar.dat',status='unknown')
c      
************************************************************

      do j = 1 , iit
c      sisu = -2.0d0
c      sisu =sisu + 0.05d0*dble(j)
c      xnH=10.0d0**sisu
      read(16,*) xNH,y(3),y(1),y(2),y(4)
     &       ,y(13),y(33),y(32),y(34),y(23),y(17),y(30)
      write(*,*)y(3),y(1),y(2),y(4)
      rho = (1.0d0+4.d0*yHe)*xm_p*xnH
      i_ev=1
C     initial condition 
      !xnH=1.d-1
c      xnH=1.d5
      !rho=(1.d0+4.d0*yHe)*xm_p*xnH*6.d0
c      T_K=3000.d0
      !T_K=30000.d0
      T_K_middle=3000d0
      T_K_High=1.d6
      T_K_Low =1.d0
      !T1=1.d6
      !T2=1.d0
c      T_K=8000.d0
      Z_metal_gas=Z_gas
      do i=1,N_sp
         y(i)=0.d0
      enddo
c     set abundance of non-zero species
c      y_Hp=3.d-4
c      y_H2=2.d-5
c      y_Hp=0.5d0
c      y_H2=0.d0
      y_Hp=1.d-4
      y_H2=1.d-6
      y_Dp=0.d0
      y_HD=0.d0
c      y_HD=1.d-9
c
      y_H=1.d0-y_Hp-2.d0*y_H2-y_HD
c      y(1)=y_H
c      y(2)=y_H2
c      y(4)=y_Hp
      y(8)=yHe
c
      y_D=yD-y_Dp-y_HD
      y(12)=y_D
      y(13)=y_HD
      y(14)=y_Dp
c     
c     for fiducial dust model
      yC=0.927d-4*Z_metal_gas
      yO=3.568d-4*Z_metal_gas 
c
c     no dust model (Anders & Grevesse 1989)
cc      yC=3.58d-4*Z_metal_gas
c      yC=3.97d-4*Z_metal_gas
c      yO=8.49d-4*Z_metal_gas
c
      y_C=0.d0*yC
      y_Cp=1.0d0*yC
c      y(17)=y_C
c      y(23)=y_Cp
c      y(30)=yO
c      
      y(3)=y_Hp+y_Dp+y_Cp
c     dust Z
      Z_metal=Z_dust
c    
      xmu=(1.d0+4.d0*yHe)
     &     /(y(1)+y(2)+y(3)+y(4)+y(8)+y(9)+y(10))
c
      gamma=1.d0+(1.d0+4.d0*yHe)
     &     /(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))
     &     +c_H2(T_K)*y(2)))
c
      e=xk_B*T_K/((gamma-1.d0)*xmu*xm_p)
************************************************************
      T_gr_K=10.d0
      xLmbd_ch=0.d0
c
      t=0.d0
      dt=1.d-1
      itchem=0
      t_chem=1.d-1
      esc_cnt=1.d0
C     temporal evolution
      do i=1,itr
c     collapse timescale, radius, bulk velocity
c     1) self-gravitating
c         xlmbd_J=3.d18!dsqrt(pi*xk_B*T_K/(G*xmu*xm_p*rho))
c         radius=xlmbd_J
c         xNcolJ=xnH*xlmbd_J
c     2) constant column density
         xNcol=1.d19
         radius=xNcol/xnH
c
         t_col=0.d0!dsqrt(3.d0*pi/(32.d0*G*rho))
         v_bulk=0.d0!radius/3.d0/t_col
         v_turb=0.d0
         v_bulk=dsqrt(v_bulk**2+v_turb**2)
c     continuum cooling
         A_v=Z_metal*xnH*radius*5.3d-22

         T_K_middle=0.5d0*(T_K_High + T_K_Low)
c         write(*,*) T_K_High,T_K_middle,T_K_Low,xLmbd_net_High,
c     &   xLmbd_net,xLmbd_net_Low
         T_K=T_K_middle
         call  rad_cool(Z_metal,T_K,T_gr_K,radius,A_v,esc_cnt,
     &        dt,xmu,gamma,y,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,
     &        i,i_ev)
c         tau_cont=tau_cnt
c
c     PE heating/cooling
         call  phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)
c
c     CR heating
         Gmm_CR=Gamma_CR(zeta,y_a,y_m,yHe)
c
c         v_D=dsqrt(2.d0*xk_B*T_K/(2.d0*xm_p))
c         xNc_H2=dmin1(1.d0,v_D/v_bulk)*y_H2*xnH*radius
c         v_D=dsqrt(2.d0*xk_B*T_K/(1.d0*xm_p))
c         xNc_H=dmin1(1.d0,v_D/v_bulk)*y_H*xnH*radius
c         v_D=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p))
c         xNc_HD=dmin1(1.d0,v_D/v_bulk)*y_HD*xnH*radius
c
c         call chemcool(xnH,T_K,T_gr_K,Z_metal,
c     &        y,dt,t_chem,xmu,gamma,xLmbd_ch)
c
c         y_H=y(1)
c         y_H2=y(2)
c         y_e=y(3)
c         y_HD=y(13)
c
         xLmbd_net=xLmbd_cnt+xLmbd_line+xLmbd_gr+xLmbd_ch
     &        -Gmm_pe-Gmm_CR

         T_K=T_K_High
         call  rad_cool(Z_metal,T_K,T_gr_K,radius,A_v,esc_cnt,
     &        dt,xmu,gamma,y,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,
     &        i,i_ev)
         call  phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)
         Gmm_CR=Gamma_CR(zeta,y_a,y_m,yHe)
c         v_D=dsqrt(2.d0*xk_B*T_K_High/(2.d0*xm_p))
c         xNc_H2=dmin1(1.d0,v_D/v_bulk)*y_H2*xnH*radius
c         v_D=dsqrt(2.d0*xk_B*T_K/(1.d0*xm_p))
c         xNc_H=dmin1(1.d0,v_D/v_bulk)*y_H*xnH*radius
c         v_D=dsqrt(2.d0*xk_B*T_K/(3.d0*xm_p))
c         xNc_HD=dmin1(1.d0,v_D/v_bulk)*y_HD*xnH*radius
         xLmbd_net_High=xLmbd_cnt+xLmbd_line+xLmbd_gr+xLmbd_ch
     &        -Gmm_pe-Gmm_CR

         T_K=T_K_Low
         call  rad_cool(Z_metal,T_K,T_gr_K,radius,A_v,esc_cnt,
     &        dt,xmu,gamma,y,xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,
     &        i,i_ev)
         call  phelectr(xnH,T_K,T_gr_K,y_e,Z_metal,G_0,A_v,Gmm_pe)
         Gmm_CR=Gamma_CR(zeta,y_a,y_m,yHe)
         xLmbd_net_Low=xLmbd_cnt+xLmbd_line+xLmbd_gr+xLmbd_ch
     &        -Gmm_pe-Gmm_CR

c         write(*,*) xLmbd_net_Low,xLmbd_net,xLmbd_net_High,T_K_middle
         if(((xLmbd_net*xLmbd_net_High).le.0.d0)) then
            T_K_Low=T_K_middle
         endif
         if(((xLmbd_net*xLmbd_net_Low).lt.0.d0)) then
            T_K_High=T_K_middle
         endif

c         if(((xLmbd_net*xLmbd_net_High).gt.0.d0) .and.
c     &      ((xLmbd_net*xLmbd_net_Low).le.0.d0)) then
c             T_K_High=T_K;T_K=0.5d0*(T_K+T_K_Low);T_K_Low=T_K_Low
c             write(*,*)T_K_High,T_K,T_K_Low
c         endif
c         if(((xLmbd_net*xLmbd_net_High).lt.0.d0) .and.
c     &      ((xLmbd_net*xLmbd_net_Low).ge.0.d0)) then
c             T_K_Low=T_K;T_K=0.5d0*(T_K+T_K_High);T_K_High=T_K_High
c             write(*,*) T_K_High,T_K,T_K_Low
c         endif
c         xMJ=rho*radius**3/2.d33

c     update values of rho & e
c         call update(i_ev)
c         
c     data output
         T_K=T_K_middle
         if((mod(i,10).eq.0).and.(i_ev==0))  then 
c            write(*,101) t/t_col, T_K, y_H2, y_e
 101        format(4e10.3)
         endif
         if((mod(i,10).eq.0).and.(i_ev==1))  then             
            write(*,102) xnH, T_K, T_gr_K, y_H2, y_e
 102        format(5e10.3)
            write(8,108) xnH, T_K
            write(11,201) xnH,T_K,T_gr_K!,p
c            write(12,202) xnH,y
             write(12,202) xnH,y(2),y(23),y(17),y(33)             
c            write(13,203) xnH,Gmm_cmp,xLmbd_cnt,xLmbd_line,xLmbd_gr,
c     &           xLmbd_ch,Gmm_pe,Gmm_CR
c             write(*,*)'xLmbd_net',xLmbd_net
             !xLmbd_CI=1.d-20!dsign(dmax1(dabs(xLmbd_CI),1.d-10),xLmbd_CI)
             !write(*,*)'xLmbd_CI_2',xLmbd_CI
             write(14,204) xnH, xLmbd_CII+xLmbd_CI,xLmbd_OI,
     &                     xLmbd_CII,xLmbd_CI
c            write(14,204) xnH, xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,
c     &           xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
c            write(15,108) xMJ,T_K 
c            write(16,*) xnH,T_K,T_gr_K,A_v,xNc_H2,y
         endif
 108     format(2e10.3)
 201     format(4e10.3)
 202     format(51e10.3)
 203     format(8e11.3)
 204     format(10e11.3)
c     
c     end of calculation
c         if(xnH.gt.1.d23) then
         if(xnH.gt.1.d14) then
            close(11)
            close(12)
            close(13)
            close(14)
c            close(15)
c            close(16)
            stop
         endif
c     
c     dt for next step
c         if(xLmbd_net.ne.0.d0) then
c            t_cool=e/dabs(xLmbd_net)
c         else
c            t_cool=1.d50
c         endif
c     
c         if(itchem.eq.0) then
c            dt=dmax1(1.d0*t_chem,1.2d0*dt)
c            if(dmin1(2.d-2*t_col,2.d-2*t_cool).le.dt) itchem=1
c         elseif(xnH < 1.d0) then
c            dt=dmin1(2.d-2*t_col,2.d-2*t_cool)
c         elseif(xnH < 1.d14) then
c            dt=dmin1(2.d-2*t_col,2.d-2*t_cool)
c         else
c            dt=dmin1(5.d-2*t_col,5.d-2*t_cool)
c         endif
c         print *, itchem, 'dt=', t_chem, 
c     &        dmin1(2.d-2*t_col,2.d-2*t_cool)
c         pause
c     
      enddo
      write(15,*) xnH,T_K_middle
      enddo
      close(15)
      close(16)
c
      CONTAINS
C
      SUBROUTINE update(i_ev)
      implicit real*8(a-h,o-z)    
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
         drho=0.d0!(rhoo/t_col)*dt
         Gmm_cmp=0.d0!(gamma-1.d0)*eo/t_col !Heeting by compression
         de=(Gmm_cmp-xLmbd_net)*dt         
         if(i_ev == 1) then
            rho=rhoo+drho
            xnH=rho/((1.d0+4.d0*yHe)*xm_p)
         endif
         e=eo+de           
c
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
c
         if(it > 1) ftn=ft-err*(ft-fto)/(err-erro)            
         erro=err
         fto=ft
         ft=ftn
      enddo
      t=t+dt
      if((i_ev==0).and.(t/t_col > 1.d0)) then
         i_ev=1
         print *, '***********************'
c         pause
      endif         
      return
      END SUBROUTINE update
C
C
      END PROGRAM 
C
C
C
      SUBROUTINE coltime(gamma_eff,ft)
      implicit real*8(a-h,o-z)      
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
c      
      if(f > 0.95d0) then
         ft=1.d0/dsqrt(0.05d0)
      else
         ft=1.d0/dsqrt(1.d0-f)
      endif
c
      return
      END
C
C
C
      FUNCTION Q_bg(T_nu)
      IMPLICIT REAL*8(a-h,o-z)
      COMMON /radgr/ T_gr_K,tau_cont
      PARAMETER (B_ex=5.d-4)
      common /radtemp/ T_rad
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
c      Q_bg=Q_bg_CMB
c      Q_bg=0.d0
      return
      END
c
c
      SUBROUTINE rad_cool(Z_metal,Tp,Tp_gr,
     &     radius,A_v,esc_cnt,dt,xmu,gamma,y,
     &     xLmbd_ch,xLmbd_line,xLmbd_cnt,xLmbd_gr,i,i_ev)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER(maxit=100,N_sp=50)
      common /radtemp/ T_rad 
      DIMENSION y(N_sp)
      COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc  
      DATA xm_p/1.67d-24/,xk_B/1.38d-16/
c      func(t,xLmbd)=Tp-t-(gamma-1.d0)*xmu*xm_p*xLmbd*dt/xk_B

      Tp0=Tp
      Tp1=Tp
      if(tau_cnt < 10.d0) then
         T_K=Tp1
         call line_cool(y,radius,xLmbd_line,i)
      else
         xLmbd_line=0.d0
      endif
      call cnt(Z_metal,xnH,Tp1,Tp_gr,
     &     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
c      if(dabs(Tp1/T_rad-1.d0) > 1.d0) return
c      xLmbd=xLmbd_ch+xLmbd_line+xLmbd_cnt+xLmbd_gr
c      fL=func(Tp1,xLmbd)

c      Tp2=1.01d0*Tp1
c      call cnt(Z_metal,xnH,Tp2,Tp_gr,
c     &     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
c      xLmbd=xLmbd_ch+xLmbd_line+xLmbd_cnt+xLmbd_gr
c      f=func(Tp2,xLmbd)

c      if(dabs(fL).lt.dabs(f)) then
c         Tpsec=Tp1
c         TpL=Tp2
c         swap=fL
c         fL=f
c         f=swap
c      else
c         TpL=Tp1
c         Tpsec=Tp2
c      endif
      
c      do 11 it=1,maxit
c         dTp=(TpL-Tpsec)*f/(f-fL)
c         TpL=Tpsec
c         fL=f
c         Tpsec=Tpsec+dTp
c         call cnt(Z_metal,xnH,Tpsec,Tp_gr,
c     &     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt1,xLmbd_gr1)
c         xLmbd1=xLmbd_ch+xLmbd_line+xLmbd_cnt1+xLmbd_gr1
c         f=func(Tpsec,xLmbd1)
c         err=dabs(f/Tpsec)
c         if(err < 1.d-6) go to 12
c 11   continue
c 12   continue
c      xLmbd_cnt=xLmbd_cnt1
c      xLmbd_gr=xLmbd_gr1
      T_K=Tp0
      return
      END
c
c
      SUBROUTINE cnt(Z_metal,xnH,T_K,T_gr_K,
     &     radius,A_v,esc_cnt,tau_cnt,xLmbd_cnt,xLmbd_gr)
      IMPLICIT REAL*8(a-h,o-z)
      DATA sigma_B/5.67d-5/,xm_p/1.67d-24/ 
c      PARAMETER (xk_vis=40.d0,B_ex=5.d-4)
c      PARAMETER (xk_vis=40.d0,B_ex=0.d0)
      common /radtemp/ T_rad
      yHe=8.333d-2
c      xJ_ex=B_ex*dexp(-A_v/2.5d0)
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
c     &     -xk_vis*Z_metal*xJ_ex      
      return
      END
c      
c     
      SUBROUTINE line_cool(y,radius,xLmbd_line,i)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER (N_sp=50,N_p=100,yHe=8.333d-2)
      DIMENSION y(N_sp),func(N_p),
     &     esc_CI(N_p),esc_CII(N_p),esc_OI(N_p),
     &     esc_CImeta(N_p),esc_CIImeta(N_p),esc_OImeta(N_p)
      common /radtemp/ T_rad
      COMMON /lines/ xnH,T_K,y_e,y_a,y_m,v_bulk,tau_cnt,xNc   
      COMMON /mol_lines/ A0,DT0,sigma,eta_T,xmu_mol,J_max
      common /indcool/ xLmbd_H2, xLmbd_HD, xLmbd_CO, xLmbd_OH,
     &        xLmbd_H2O, xLmbd_CII, xLmbd_CI, xLmbd_OI, xLmbd_Lya
      EXTERNAL pop_CI,pop_OI,pop_CII,
     &     pop_CImeta,pop_OImeta,pop_CIImeta
      DATA xm_p/1.67d-24/,xk_B/1.38d-16/,pi/3.14159265358979d0/
      SAVE esc_CI,esc_CII,esc_OI,
     &     esc_CImeta,esc_CIImeta,esc_OImeta
c      
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
c
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
c
      xNc_H=xnH*radius
      xNc_H2=y_m*xNc_H
      xNc_CO=y_CO*xNc_H
      xNc_OH=y_OH*xNc_H
      xNc_H2O=y_H2O*xNc_H
      xNc_CII=y_CII*xNc_H
      xNc_CI=y_CI*xNc_H
      xNc_OI=y_OI*xNc_H
      xNc_HD=y_HD*xNc_H
c
      v_th=dsqrt(2.d0*xk_B*T_K/xm_p)
c     
c     H2 cooling
      v_D=v_th/dsqrt(2.d0)
      xNc_H2=dmin1(1.d0,v_D/v_bulk)*xNc_H2
      xLmbd_H2=y_m
     &     *xLd_H2(xnH,T_K,y_a,y_m,y_e,y_Hp,xNc_H2,tau_cnt)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c
c     HD cooling
      v_D=v_th/dsqrt(3.d0)
      xNc_HD=dmin1(1.d0,v_D/v_bulk)*xNc_HD      
      xLmbd_HD=y_HD*xLd_HD(xnH,T_K,y_a,y_m,xNc_HD,tau_cnt)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c
c     CO cooling
      xmu_mol=28.d0
      v_D=v_th/dsqrt(xmu_mol)
      xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CO
      call COcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_CO)
      call COcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_CO)
      xLmbd_CO=y_CO
     &     *(xLd_CO-xLdr_CO)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c
c     OH cooling
      xmu_mol=17.d0
      v_D=v_th/dsqrt(xmu_mol)
      xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OH
      call OHcool(xnH,T_K,y_m,xNc,tau_cnt,xLd_OH)
      call OHcool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_OH)
      xLmbd_OH=y_OH
     &     *(xLd_OH-xLdr_OH)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c
c     H2O cooling
      xmu_mol=18.d0
      v_D=v_th/dsqrt(xmu_mol)
      xNc=dmin1(1.d0,v_D/v_bulk)*xNc_H2O
      call H2Ocool(xnH,T_K,y_m,xNc,tau_cnt,xLd_H2O)
      call H2Ocool(xnH,T_rad,y_m,xNc,tau_cnt,xLdr_H2O)
      xLmbd_H2O=y_H2O
     &     *(xLd_H2O-xLdr_H2O)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c
c     CII cooling
      xmu_C=12.d0
      xmu_O=16.d0
      v_D=v_th/dsqrt(xmu_C)
      xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CII
c     fine structure line
      xLmbd_CII=y_CII
     &     *xLd(pop_CII,esc_CII,1)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c     metastable line
      if(T_K > 2.d2) then
         xLmbd_CIImeta=y_CII
     &        *xLd(pop_CIImeta,esc_CIImeta,1)
     &        *xnH/((1.d0+4.d0*yHe)*xm_p)
         xLmbd_CII=xLmbd_CII+xLmbd_CIImeta
      endif
c
c     CI cooling
      xNc=dmin1(1.d0,v_D/v_bulk)*xNc_CI
c     fine structure line
      xLmbd_CI=y_CI
     &     *xLd(pop_CI,esc_CI,3)
     &     *xnH/((1.d0+4.d0*yHe)*xm_p)
c     
      if(T_K > 3.d3) then
         xLmbd_CImeta=y_CI
     &        *xLd(pop_CImeta,esc_CImeta,3)
     &        *xnH/((1.d0+4.d0*yHe)*xm_p)
         xLmbd_CI=xLmbd_CI+xLmbd_CImeta
      endif
c     
c     OI cooling
      v_D=v_th/dsqrt(xmu_O)
      xNc=dmin1(1.d0,v_D/v_bulk)*xNc_OI
c     fine structure line
      if(T_K > 15.d0) then
         xLmbd_OI=y_OI
     &        *xLd(pop_OI,esc_OI,3)
     &        *xnH/((1.d0+4.d0*yHe)*xm_p)
      else
         xLmbd_OI=0.d0
      endif
c     metastable line
      if(T_K > 3.5d3) then
         xLmbd_OImeta=y_OI
     &        *xLd(pop_OImeta,esc_OImeta,3)
     &        *xnH/((1.d0+4.d0*yHe)*xm_p)
         xLmbd_OI=xLmbd_OI+xLmbd_OImeta
      endif

c     HI Lya cooling 
      if(xnH < 1.d5) then
         T_5=1.d-5*T_K
         xLmbd_Lya=y_e*y_a*7.50d-19/(1.d0+dsqrt(T_5))
     &        *dexp(-1.18348d5/T_K)
     &        *xnH/((1.d0+4.d0*yHe)*xm_p)
      else
         xLmbd_Lya=0.d0
      endif
c     total cooling rate by lines
      xLmbd_line=xLmbd_H2+xLmbd_HD+xLmbd_CO+xLmbd_OH+xLmbd_H2O
     &     +xLmbd_CII+xLmbd_CI+xLmbd_OI+xLmbd_Lya         
      return
      END
c
c    
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





















