program main
  implicit none
  double precision B0,B1,rho0,rho1,T0,T1,p0,p1,v0,v1,rhopre,vpre,ppre,bpre,time
  double precision :: gamma,nu,n0,n1,M,cs,TtoK=1.d3,PtoPkb,dm,rhoinv,vnew,pnew,rhonew,bnew
  double precision AbsL,Laml,Lamc,CooL,N,T,metal,ekin,emag,dx,dt,Lcf,Tnew,nnew
  double precision, parameter :: kb=8.63359d0
  integer i,tstep
  !input
  rho0=1.0d0
  v0=48.89d0
  B0=1.57734d0
  p0=8.810807d3*kb*1.d-3
  gamma=5.0d0/3.d0
  nu=1.27d0
  PtoPkb= 1.52d0/1.38d0 * 1.d2
  metal=1.d0
  dx=1.d0
  dt=0.001d0
  tstep=1000
  M=10.d0
  !input

  n0=rho0/nu !(rho[mp/cc]/{mu*mp})=1/cc
  T0=P0/(kb*n0)
  cs=dsqrt(gamma*P0/rho0)
  M=v0/cs

  !shock jump (M >> 1)
  n1=4.0d0*n0
  rho1=4.0d0*rho0
  v1=0.250d0*v0
  B1=4.0d0*B0
  P1=p0 * 2.d0*gamma/(1.d0+gamma)*M**2
  T1=P1/(kb*n1)
  ekin=rho0*v0*v0
  emag=B1*B1*0.5d0
  time=0.d0

  AbsL = 3.94656d1
  TtoK   = 1.0d3

  rhopre=rho1
  vpre=v1
  ppre=p1
  bpre=B1

  open(100,file='cooling.dat')
  write(100,*) n0,p0*PtoPkb,rho0,v0,p0,b0,n0,T0*TtoK,time
  write(100,*) n1,p1*PtoPkb,rho1,v1,p1,b1,n1,T1*TtoK,time
  do i=1,tstep
     call Stblty(dt,n1,T1,gamma,dx,metal,AbsL,TtoK,p1)
     dt=dt*0.1d0

     time=time+dt
     !cooling function
     !-------------------------( Cooling & Heating )
     Laml = 1.0d7 * dexp( -1.184d5/(T1*TtoK+1.0d3) )
     Laml = dmin1(Laml,5.d3)
     Lamc = metal * 1.4d-2 * dsqrt(T1*TtoK) * dexp( -9.2d1/(T1*TtoK) )
     CooL = ( n1**2.d0 * (Laml + Lamc) - n1 ) * AbsL

     write(*,*)p1,rho1,cool

!     if( dt .le. 0.2d0*p1/(gammi1*dabs(CooL)) ) then
!        p1 = p1 - gammi1*CooL*dt !explicit
!     else
!        Call IMC( P1,n1,dt,N(i) ) !implicit
!        U(i,5) = Pn(i)
!     end if

     dm=dx*(rho1-rhopre+1.d-7)

     rhoinv=1.d0/rho1+(v1-vpre)/dm * dt
     rhonew=1.d0/rhoinv
     vnew = v1-(p1-ppre+b1**2/2.d0-bpre**2/2.d0)/dm * dt
     pnew = p1-(gamma-1.d0)*p1*(v1-vpre)/dm - (gamma-1.d0)*CooL/rho1/T1
     bnew = rhonew*B1/rho1

     rhopre=rho1
     vpre=v1
     ppre=p1
     bpre=B1
     rho1=rhonew
     v1=vnew
     p1=pnew
     b1=bnew
     n1=rhonew/nu
     T1=Pnew/(kb*nnew)

     write(100,*) n1,p1*PtoPkb,rho1,v1,p1,b1,n1,T1*TtoK,time
  end do
  close(100)

end program main

SUBROUTINE Stblty(tLMT,N,Tn,gamma,dx,metal,AbsL,TtoK,p1)
!USE comvar
implicit none
double precision  tLMT,alpha,tauC,N,Tn,dx,TtoK,Laml,metal
double precision  CooL,gamma,gammi1,AbsL,Lamc,p1


!do i = 1, Ncell
   gammi1 = gamma-1.d0
   !*** Stability for Conduction ***!
   !alpha = gammi1*Kcond*dsqrt(Tn)
   !alpha = Nn*kb*dx**2/alpha
   !*** Avoid over cooling ***!
   !Call Fcool(CooL,N,Tn)
   Laml = 1.0d7 * dexp( -1.184d5/(Tn*TtoK+1.0d3) )
   Laml = dmin1(Laml,5.d3)
   Lamc = metal * 1.4d-2 * dsqrt(Tn*TtoK) * dexp( -9.2d1/(Tn*TtoK) )
   CooL = ( n**2.d0 * (Laml + Lamc) - n ) * AbsL

   tauC = p1/gammi1/dabs(CooL)

   !tLMT =  dmin1( tLMT, alpha )
   tLMT =  dmin1( tLMT, tauC  )
!end do
if(tLMT.lt.0.d0) write(*,*) tLMT,'err at Stblty'

END SUBROUTINE Stblty

!SUBROUTINE IMC( P,n,dt,NN,gamma )
!USE comvar
!double precision P,n,T,dt,NN,gamma,gammi1
!double precision Pu,Pd,Pm,fev,iud
!double precision CooL,Pmold,nkbi
!integer kkk
!gammi1=gamma-1.d0
!Pu    = 1.d10
!Pd    = 1.d-10
!Pm    = 1.d1
!Pmold = -1.d0
!nkbi  = 1.d0/(kb*n)
!do kkk = 1,30
!  T = Pm*nkbi
!  call Fcool(CooL,NN,T)
!  fev = Pm + CooL*dt*gammi1 - P
!  iud = dsign(1.d0,fev)
!  Pd  = dmax1(-iud*Pm,Pd)
!  Pu  = dmin1(2.d0*Pm/(iud+1.d0+1.d-10),Pu)
!  Pm  = 1.d1**(0.5d0*(dlog10(Pu)+dlog10(Pd)))
!  if(dabs(Pm-Pmold).lt.1.d-5) goto 835
!  Pmold = Pm
!end do
!835 continue
!P = Pm
!END SUBROUTINE IMC
