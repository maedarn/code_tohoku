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
  double precision :: rhoc,rc,cs0,fratio,r0,Mint
  double precision :: Msun=2.4d-2,AU=2.d0*1.d5
  double precision,allocatable,dimension(:) :: tint,xint,vint,xexa,vexa
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  integer :: i,iend,irc
  character(21) :: fname='/Users/maeda/Desktop/'

  tstart=0.d0
  tend  =10.d0
  xstart=0.d0
  vstart=0.d0
  tstart=0.d0
  dt=0.01d0
  iend=5000
  cs0=0.19d0/0.953d0
  rhoc=1.d0*1.d4*2.d0
  rc=6.45d0*cs0/dsqrt(G4pi*rhoc)
  fratio=20.d0
  r0=cs0/dsqrt(G4pi*rhoc)

  write(*,*) rc*AU,rc**3*rhoc*Msun,rc,rhoc,r0
  
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
  irc=0
  Mint=0.d0
  open(10,file=fname//'runge4BE.dat')
  do i=1,iend
  if(tint(i)*r0<rc)then
  Mint=Mint+4.d0*3.141592d0*rhoc*dexp(-xint(i))*tint(i)*r0*tint(i)*r0*dt*r0*Msun*fratio
  write(*,*) Mint,rhoc*dexp(-xint(i)),tint(i)*r0,4.d0*3.141592d0*rhoc*dexp(-xint(i))*tint(i)*r0*tint(i)*r0*dt*r0*Msun
  irc=i
  write(10,*) tint(i)*r0,',',xint(i),',',vint(i)/r0,',',fratio*rhoc*dexp(-xint(i))
  else
  write(10,*) tint(i)*r0,',',xint(irc),',',vint(irc)/r0,',',fratio*rhoc*dexp(-xint(irc))
  endif
  end do
  close(10)

  open(10,file=fname//'BE.DAT')
  do i=1,iend
  if(tint(i)*r0<rc)then
  irc=i
  write(10,*) tint(i)*r0,fratio*rhoc*dexp(-xint(i))
  else
  write(10,*) tint(i)*r0,fratio*rhoc*dexp(-xint(irc))
  endif
  end do
  close(10)


  deallocate(tint,xint,vint,xexa,vexa)
 
  !close(10)
 
end program main
 
function funcf(t,x,v)
  double precision :: t,x,v
  double precision :: funcf
 
  !funcf=(1.0-x**2)
  funcf=-2.d0*t*v+t*t*dexp(-x)
 
  return
end function funcf



function funcv(t,x,v)
  double precision :: t,x,v
  double precision :: funcv
 
  !funcv=(1.0-x**2)
  funcv=v!(1.0-x**2)
 
  return
end function funcv

