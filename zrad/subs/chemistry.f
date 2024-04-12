
      SUBROUTINE chemcool(xnH,Tp,Tp_gr,Z_metal,
     &     y,dt,tchem,xmu,gamma,xLmbdch)
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER(maxit=100,N_sp=50)      
      DIMENSION y(N_sp),ytmp(N_sp)
      DATA xm_p/1.67d-24/,xk_B/1.38d-16/

C     input    Tp,xnH,dt,y,xmu,gamma  
C     output   y,xmu,gamma,xLmbdch

      func(t,xLmbdch)=Tp-t-(gamma-1.d0)*xmu*xm_p*xLmbdch*dt/xk_B

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
      
      do 11 i=1,maxit
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
 11   continue 
 12   continue 
      do 13 j=1,N_sp
         y(j)=ytmp(j)
 13   continue
      xmu=xmu1
      gamma=gamma1
      xLmbdch=xLmbdch1
      tchem=tchem1
      return
      END


      SUBROUTINE chemreact(xnH,T_K,T_gr_K,Z_metal,
     &     y,dt,tchem,xmu,gamma,xLmbdch)
******************************************
*      Version 4+      (13 Sep 1998)     *
*      no radiation                      *
******************************************
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER(N_sp=50,N_react=675,eps=1.d-4,eps_y=1.d-10)
      PARAMETER (xnH_eq=3.d16)
      DIMENSION y(N_sp),y_init(N_sp),y_tmp(N_sp),dy(N_sp),ddy(N_sp)
      DIMENSION xk(N_react)
      DIMENSION r_f(N_sp),r_f_fw(N_sp),r_f_bw(N_sp)
      DIMENSION dr_fdy(N_sp,N_sp),A(N_sp,N_sp)  
      common /xcrit/xH,xH2,xHe

      do isp=1,N_sp
         y_init(isp)=y(isp)
         dy(isp)=0.d0
      enddo

      if(xnH.gt.xnH_eq) go to 50         

**************************************
*    Non Equilibrium (Low Density)   *
**************************************     
      xH=y(1)
      xH2=2.d0*y(2)
      xHe=y(8)
c
      call react_coef(xnH,T_K,T_gr_K,Z_metal,xk)
c
      do 20 itr=1,20
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

C ------ set the matrix A
         do isp=1,N_sp
            do jsp=1,N_sp
               if(isp.ne.jsp) then
                  A(isp,jsp)=-dt*dr_fdy(isp,jsp)
               else
                  A(isp,jsp)=1.d0-dt*dr_fdy(isp,jsp)
               endif
            enddo
         enddo   
         
C ------- set the vector ddy
         do isp=1,N_sp
            ddy(isp)=r_f(isp)*dt-dy(isp)
         enddo

C ------- solve linear equations 
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
 20   continue
      go to 100

 50   continue
***********************************************
*        Equillibrium ( High Density )        *
***********************************************
      call equichem
     &     (xnH,T_K,y_e,y_HI,y_HII,y_H2,y_HeI,y_HeII,y_HeIII)
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
c       
      yHe=y(8)+y(9)+y(10)

      xmu=(1.d0+4.d0*yHe)
     &     /(y(1)+y(2)+y(3)+y(4)+y(8)+y(9)+y(10))

      gamma=1.d0+(1.d0+4.d0*yHe)
     &     /(xmu*(1.5d0*(y(1)+y(3)+y(4)+y(8)+y(9)+y(10))
     &     +c_H2(T_K)*y(2)))
      
      xm_p=1.67d-24
C     H_2 formation/dissociation cooling
      if(xnH.lt.1.d13) then
         xn_cr=1.d6/dsqrt(T_K)
     &        /(1.6d0*y(1)*dexp(-(4.d2/T_K)**2)
     &        +1.4d0*y(2)*dexp(-1.2d4/(T_K+1.2d3)))
         crit=1.d0/(1.d0+xn_cr/xnH)
C     2 H + grain -> H2
         rtgrain=xk(23)*y(1)        
C     H- + H -> H2 +  e
         rtHm=xk(8)*y(7)*y(1)
C     H2+ + H  -> H2 + H+
         rtH2p=xk(10)*y(5)*y(1)
C     3 H  -> H2 +  H 
C     2H  + H2  -> 2H2
c     2 H + He  -> H2 + He :43
         rt3b=(xk(19)*(y(1)**3)+xk(20)*(y(1)**2)*y(2)
     &        +xk(43)*(y(1)**2)*y(8))*xnH
C     H2 + H  -> 3 H
C     2 H2  -> 2 H + H2
c     H2 + He -> 2 H + He :37
         rtdis=xk(13)*y(2)*y(1)+xk(21)*(y(2)**2)
     &        +xk(37)*y(2)*y(8)
C     H2 + ph. -> 2H
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
      
      return
      END


      SUBROUTINE equichem
     &     (xnH,T_K,y_e,y_HI,y_HII,y_H2,y_HeI,y_HeII,y_HeIII)
***********************************************
*        Equillibrium  Chemistry    (Saha)    *
***********************************************
      IMPLICIT REAL*8(a-h,o-z)
      PARAMETER(yHe=8.333d-2)
      DATA pi/3.14159d0/,rmasse/9.109534d-28/,h/6.626176d-27/,
     &     rmassp/1.67d-24/

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
      fH2=(zH0**2/zH2)*(1.60219d-12*pi*rmassp/h**2)**1.5d0
     &     *(T_eV**1.5d0)*dexp(-chiH2/T_eV)/xnH

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
c     Newton-Raphson Loop
         do 11 itr=1,1000
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
 11      continue
 12      continue
      
         y_HI=5.d-1*d
         y_HII=fH0*y_HI/y_e
         y_H2=y_HI**2/fH2
         y_HeI=yHe/(1.d0+fHe0/y_e*(1.d0+fHep/y_e))
         y_HeII=fHe0*y_HeI/y_e
         y_HeIII=fHep*y_HeII/y_e
         return
      endif
      END


















