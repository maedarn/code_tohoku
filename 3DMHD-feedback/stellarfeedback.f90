SUBROUTINE feedback(dt,mode)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION  :: dt,dxi,intt=0.d0,raninp,ranout,Rst,rloop
DOUBLE PRECISION  :: Q,lQ,aq=-39.3178d0,bq=221.997d0,cq=-227.456d0,dq=117.410d0,eq=-30.1511d0,fq=3.06810d0
DOUBLE PRECISION  :: LSFE=0.02d0,mstarmn=1.d0,pstar,Mstar1,Mstarlp
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt,jtime=0,idum,nid=0
INTEGER, parameter :: nstar=1000000
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU
character(3) Nfinal,itime
integer :: i,j,k,nlpmass=100000,itrrst=10
double precision :: div(-1:ndx,-1:ndy,-1:ndz,1:4)
double precision :: rhoth=5.d0*1.d5*1.27d0,SFE,SFratio=0.5d0
double precision :: Mmassive=10.d0,Msun=2.4d-2 !3*3*3*1.7*10^(18+18+18-24)/2*10^33=3*3*3*1.7/2 *10^(-3)
DOUBLE PRECISION ran
double precision :: Ustar(1:10,1:nstar)


if(mode==1) then

intt=intt+dt

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
div(i,j,k,1)=(U(i-1,j,k,2)-U(i+1,j,k,2))*0.5d0/dx1+(U(i,j-1,k,3)-U(i,j+1,k,3))*0.5d0/dy1+(U(i,j,k-1,4)-U(i,j,k+1,4))*0.5d0/dz1
div(i,j,k,2)=U(i-1,j,k,2)*U(i+1,j,k,2)
div(i,j,k,3)=U(i,j-1,k,3)*U(i,j+1,k,3)
div(i,j,k,4)=U(i,j,k-1,4)*U(i,j,k+1,4)
end do; end do; end do

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
if((div(i,j,k,1)<0.d0).and.(div(i,j,k,2)<0.d0).and.(div(i,j,k,3)<0.d0).and.(div(i,j,k,4)<0.d0).and.(U(i,j,k,1)>rhoth)) then
SFE=LSFE*dsqrt(U(i,j,k,1)/(4.04d0*1.d3))*dt
!SFE=LSFE*U(i,j,k,1)/dsqrt(G*U(i,j,k,1))*dx1*dy1*dz1*dt
!pstar=(1-dexp(-SFE/mstarmn))
!call ran0(ran,1)
!if(ran<pstar) then
Mstar1=U(i,j,k)*SFE*Msun
Mstarlp=0.d0
idum=1
do i=1,nlpmass
if(Mstarlp>Mstar1) then
goto 2021
endif

call ran0(raninp,idum)
call IMF(raninp,ranout)
Mstarlp=Mstarlp+ranout

!if(ran<pstar) then
if(Mmassive<ranout) then
nid=nid+1
Ustar(1,nid)=dble(IST*Ncellx + i)
Ustar(2,nid)=dble(JST*Ncelly + j)
Ustar(3,nid)=dble(KST*Ncellz + k)
Ustar(4,nid)=U(i,j,k,2)
Ustar(5,nid)=U(i,j,k,3)
Ustar(6,nid)=U(i,j,k,4)
Ustar(7,nid)=ranout !U(i,j,k,1)*dx1*dy1*dz1*SFratio
Ustar(8,nid)=intt
Ustar(9,nid)=dble(nid)
lQ=aq+bq*dlog10(ranout)+cq*dlog10(ranout)**2+dq*dlog10(ranout)**3+eq*dlog10(ranout)**4+fq*dlog10(ranout)**5
!Q=10.d0**lQ
Ustar(10,nid)=10.d0**lQ
endif
enddo
2021 continue
endif

end do; end do; end do
endif

if(mode==2) then
i=mod(int(Ustar(1,nid)),Ncellx)
j=mod(int(Ustar(1,nid)),Ncelly)
k=mod(int(Ustar(1,nid)),Ncellz)

Rst=10.d0*(Ustar(10,nid)/10.d49)**(1.d0/3.d0)*(U(i,j,k)/10.d0)**(-2.d0/3.d0) !loop


endif

if(mode==3) then
Ustar(1,nid)=Ustar(1,nid)+Ustar(4,nid)*dt/dx
Ustar(2,nid)=Ustar(2,nid)+Ustar(5,nid)*dt/dy
Ustar(3,nid)=Ustar(3,nid)+Ustar(6,nid)*dt/dz
endif



end SUBROUTINE feedback

program main
  implicit none
  double precision :: tstart,xstart,tend,vstart
  double precision :: t,x,v,dt
  double precision :: t0,t1,t2,t3
  double precision :: kx0,kx1,kx2,kx3
  double precision :: kv0,kv1,kv2,kv3
  double precision :: f0,f1,f2,f3
  double precision :: v0,v1,v2,v3
  double precision :: funcf,funcv
  double precision,allocatable,dimension(:) :: tint,xint,vint,xexa,vexa
  integer :: i,iend
  character(21) :: fname='/Users/maeda/Desktop/'

  tstart=0.d0
  tend  =100.d0
  xstart=0.d0
  vstart=1.d0
  tstart=0.d0
  dt=0.1d0
  iend=100000
  
  allocate(tint(1:iend),xint(1:iend),vint(1:iend),xexa(1:iend),vexa(1:iend))
  tint(:)=0.d0
  xint(:)=0.d0
  vint(:)=0.d0
 
 
  t=tstart
  x=xstart
  v=vstart
  do i=1,iend
    !write(10,'(2f13.5)') t,x
 
    t0=t
    kv0=v
    kx0=x
    f0=funcf(t0,kx0,kv0)
    v0=funcv(t0,kx0,kv0)
 
    t1=t+dt/2.0
    kv1=v+f0*dt/2.0
    kx1=x+v0*dt/2.0
    f1=funcf(t1,kx1,kv1)
    v1=funcv(t1,kx1,kv1)
 
    t2=t+dt/2.0
    kv2=v+f1*dt/2.0
    kx2=x+v1*dt/2.0
    f2=funcf(t2,kx2,kv2)
    v2=funcv(t2,kx2,kv2)
 
    t3=t+dt
    kv3=v+f2*dt
    kx3=x+v2*dt
    f3=funcf(t3,kx3,kv3)
    v3=funcv(t3,kx3,kv3)
 
    v=v+(f0+f1*2.0+f2*2.0+f3)*dt/6.0
    x=x+(v0+v1*2.0+v2*2.0+v3)*dt/6.0

    t=t+dt

    tint(i)=t
    xint(i)=x
    vint(i)=v
    xexa(i)=vstart*dsin(t)
    vexa(i)=-vstart*dcos(t)
 
  end do

  open(10,file=fname//'runge4.dat')
  do i=1,iend
  write(10,*) tint(i),',',xint(i),',',vint(i),',',xexa(i),',',vexa(i)
  end do
  close(10)

  deallocate(tint,xint,vint,xexa,vexa)
 
  !close(10)
 
end program main
 
function funcf(t,x,v)
  double precision :: t,x,v
  double precision :: funcf
 
  !funcf=(1.0-x**2)
  funcf=-x!9.8d0!(1.0-x**2)
 
  return
end function funcf



function funcv(t,x,v)
  double precision :: t,x,v
  double precision :: funcv
 
  !funcv=(1.0-x**2)
  funcv=v!(1.0-x**2)
 
  return
end function funcv



SUBROUTINE IMF(raninp,ranout)
integer :: nmesh,cnc1,cnc2,nran,idum
integer :: i
double precision :: raninp,ranout
double precision :: mmax,mmin,nsisu,scale,slide,dmy,Ncumcnc1,Ncumcnc2
double precision :: M1,M2,coefflow,coeffmid,coeffhig,pwllow,pwlmid,pwlhig,coefflow2,coeffmid2,coeffhig2
double precision, allocatable, dimension(:) :: M,lgM,dN,N,lgN,Ncum,Mran,ranini,ts1,ts2,ts3,ts4
double precision :: dm,lgmmax,lgmmin,ddN,Ntot
character(21) :: dir='/Users/maeda/Desktop/'

!-------paramerter------
!nmesh = 1000
!nran=10000000
mmin = 0.01d0
mmax = 150.d0
Ntot = 100.d0
!slide=5.d0
M1=0.08d0
M2=0.5d0
pwllow=-0.3d0
pwlmid=-1.3d0
pwlhig=-2.3d0
coeffhig=1.d0
coeffmid=(coeffhig*(M2)**pwlhig)/((M2)**pwlmid)
coefflow=(coeffmid*(M1)**pwlmid)/((M1)**pwllow)
coefflow2=1.d0
coeffmid2=(coefflow2*(M1)**(pwllow))/((M1)**(pwlmid))
coeffhig2=(coeffmid2*(M2)**(pwlmid))/((M2)**(pwlhig))
!-------paramerter------

!allocate(M(0:nmesh),lgM(0:nmesh),dN(0:nmesh),N(0:nmesh),lgN(0:nmesh),Ncum(0:nmesh),Mran(1:nran),ranini(1:nran))
!allocate(ts1(1:nran),ts2(1:nran),ts3(1:nran),ts4(1:nran))

!lgmmax = dlog10(mmax)
!lgmmin = dlog10(mmin)
!dm = (lgmmax-lgmmin)/dble(nmesh)
!cnc1=int((dlog10(M1)-dlog10(mmin))/dm)
!cnc2=int((dlog10(M2)-dlog10(mmin))/dm)



Ntot=coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))+&
     coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
     coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncumcnc2=coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2     )**(pwlhig+1.d0)/(pwlhig+1.d0))
Ncumcnc1=coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1     )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
         coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2     )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncumcnc2=Ncumcnc2/Ntot
Ncumcnc1=Ncumcnc1/Ntot



if (raninp>Ncumcnc1) then
ranout=((coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
          coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))-coeffmid*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-coeffhig*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*raninp)/coefflow*(pwllow+1.d0))**(1.d0/(pwllow+1.d0))

else if ((Ncumcnc1>raninp).and.(raninp>Ncumcnc2)) then
ranout=((coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*raninp)*(pwlmid+1.d0)/coeffmid)**(1.d0/(pwlmid+1.d0))

else if (Ncumcnc2>raninp) then
ranout=((coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-Ntot*raninp)*(pwlhig+1.d0)/coeffhig)**(1.d0/(pwlhig+1.d0))

endif


!deallocate(M,lgM,dN,N,lgN,Ncum,Mran,ranini)
!deallocate(ts1,ts2,ts3,ts4)
end SUBROUTINE IMF
