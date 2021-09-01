program main
  implicit none
  double precision B0,B1,rho0,rho1,T0,T1,p0,p1,v0,v1,rhopre,vpre,ppre,bpre,time
  double precision :: gamma,nu,n0,n1,M,cs,TtoK=1.d3,PtoPkb,dm,rhoinv,vnew,pnew,rhonew,bnew
  double precision AbsL,Laml,Lamc,CooL,N,T,metal,ekin,emag,dx,dt,Lcf,Tnew,nnew
  double precision Lamm,Lamh,Tsisu,CooLonly,nst,nls,CooLst,CooLls,CooLmd,CooLonlyKI
  double precision, parameter :: kb=8.63359d0
  integer i,tstep,j,nloop
  character(21) :: fname='/Users/maeda/Desktop/'
  !input

  tstep=40
  nloop=50
  metal=1.d0

  open(100,file=fname//'cooling.dat')
  !write(100,*) n0,p0*PtoPkb,rho0,v0,p0,b0,n0,T0*TtoK,time
  !write(100,*) n1,p1*PtoPkb,rho1,v1,p1,b1,n1,T1*TtoK,time
  do i=1,tstep
     Tsisu = (2.0d0-3.d0)! -10^3
     Tsisu =Tsisu + 0.1d0*dble(i)
     T1 = 10.0d0**Tsisu
     !-------------------------( Cooling & Heating )
     
     !nst=0.001d0
     !nls=1.d6
     !do j=1,nloop
     
     if(T1*TtoK<1.4577d4) then
     Laml = 1.0d7 * dexp( -1.184d5/(T1*TtoK+1.0d3) )
     !Laml = dmin1(Laml,5.d3)
     Lamc = metal * 1.4d-2 * dsqrt(T1*TtoK) * dexp( -9.2d1/(T1*TtoK) )
     !CooL = ( n1**2.d0 * (Laml + Lamc) - n1 ) * AbsL
     CooLonly = Laml + Lamc
     elseif((1.4577d4<T1*TtoK).and.(T1*TtoK<1.9449d4)) then
     Lamm = 5.d3 + metal * 1.4d-2 * dsqrt(T1*TtoK) * dexp( -9.2d1/(T1*TtoK) )
     CooL = ( n1**2.d0 * (Lamm) - n1 ) * AbsL
     CooLonly = Lamm
     elseif(1.9449d4<T1*TtoK) then
     Lamh = 3.75d4*(1.d0-dtanh((T1*TtoK-2.d5)/(2.d5)))*dexp(-5.d4/(T1*TtoK))+1.d3*dexp(-5.d4/(T1*TtoK))
     CooL = ( n1**2.d0 * (Lamh) - n1 ) * AbsL
     CooLonly = Lamh
     endif
     Laml = 1.0d7 * dexp( -1.184d5/(T1*TtoK+1.0d3) )
     Lamc = metal * 1.4d-2 * dsqrt(T1*TtoK) * dexp( -9.2d1/(T1*TtoK) )
     CooLonlyKI=Laml + Lamc

     
     write(100,*) T1*TtoK,',',CooLonly,',',CooLonlyKI
  end do
  close(100)
end program main
