MODULE comvar
implicit none
integer, parameter :: ndx=128,ndy=128
double precision :: Phi(0:ndx+1,0:ndy+1)
double precision :: dx,dy
END MODULE comvar

program main
  use comvar
  implicit none
  double precision :: tstart,xstart,ystart,tend,vxstart,vystart
  double precision :: t,x,y,vx,vy,dt
  double precision :: t0,t1,t2,t3!,dx,dy
  double precision :: kx0,kx1,kx2,kx3
  double precision :: kv0,kv1,kv2,kv3
  double precision :: f0,f1,f2,f3
  double precision :: v0,v1,v2,v3
  double precision :: funcf,funcv,rrsph3,dinit1
  double precision,allocatable,dimension(:) :: tint,xint,yint,vxint,vyint,xexa,yexa,vxexa,vyexa
  !double precision :: Msun=2.4d-2,G=1.11142d-4,Lbox,Msunex,Msun=2.47077d-2
  double precision :: Msun=2.4d-2,G=1.11142d-4,Lbox,Msunex
  !double precision,allocatable,dimension(:,:) :: Phi
  
  integer :: i,iend,mode!,ndx,ndy
  character(21) :: fname='/Users/maeda/Desktop/'

  Msunex=1.672622d0*1.d-24 * (3.085677857083d0*1.d18)**3 / (1.98892d0 * 1.d33)
  !write(*,*) Msunex,vxstart,'Msunex'

  tstart=0.d0
  tend  =100.d0
  xstart=0.5d0
  vxstart=1.2d0/0.953d0
  !vxstart=dsqrt(G*100.d0/Msun/0.3d0)/0.953d0
  ystart=0.2d0
  vystart=0.d0
  tstart=0.d0
  dt=0.001d0
  iend=10000
  !ndx=128
  !ndy=128
  Lbox=1.d0
  dx=Lbox/(dble(ndx))
  dy=Lbox/(dble(ndy))
  rrsph3=0.1d0
  dinit1=100.d0/Msun/(4.d0*3.141592d0/3.d0 *  rrsph3**3)

  write(*,*) Msunex,vxstart,'Msunex'

  
  allocate(tint(1:iend),xint(1:iend),vxint(1:iend),xexa(1:iend),vxexa(1:iend),yint(1:iend),vyint(1:iend),yexa(1:iend),vyexa(1:iend))
  !allocate(Phi(0:ndx+1,0:ndy+1))
  tint(:)=0.d0
  xint(:)=0.d0
  vxint(:)=0.d0
  yint(:)=0.d0
  vyint(:)=0.d0

  call ini(rrsph3,dinit1)
 
 
  t=tstart
  x=xstart
  vx=vxstart
  y=ystart
  vy=vystart

  write(*,*)xstart,ystart,vxstart,vystart,dx,dy
  write(*,*)x,y,vx,vy
  do i=1,iend
    !write(10,'(2f13.5)') t,x
 
    !call rng(vx,x,y,dt*0.5d0,1)
    call rng(vx,x,y,dt,1)
    call rng(vy,y,x,dt,2)
    !call rng(vx,x,y,dt*0.5d0,1)
    t=t+dt

    !write(*,*)x,y,vx,vy
    !write(*,*)y

    tint(i)=t
    xint(i)=x
    vxint(i)=vx
    yint(i)=y
    vyint(i)=vy

    xexa(i)=vxstart*dsin(t)
    vxexa(i)=-vxstart*dcos(t)
    yexa(i)=vystart*dsin(t)
    vyexa(i)=-vystart*dcos(t)
 
  end do

  open(10,file=fname//'runge4.dat')
  do i=1,iend
  write(10,*) tint(i),',',xint(i),',',vxint(i),',',xexa(i),',',vxexa(i),',',yint(i),',',vyint(i),',',yexa(i),',',vyexa(i)&
              ,',',vxint(i)*vxint(i)+vyint(i)*vyint(i)!-vxstart*vxstart
  end do
  close(10)

  deallocate(tint,xint,vxint,xexa,vxexa,yint,vyint,yexa,vyexa)
 
  !close(10)
 
end program main
 
function funcf(x,v,xdm,mode)
  use comvar
  double precision :: t,x,v,xdm,xmod,ymod,delx1,delx2,delx
  double precision :: funcf
  integer :: mode,i,j

  !write(*,*)'mode',mode
  !funcf=(1.0-x**2)
  if(mode==1) then
  !funcf=-x!9.8d0!(1.0-x**2)
  i=int(x/dx)
  j=int(xdm/dy)
  xmod=x-i*dx
  ymod=xdm-j*dy
  xmod=xmod/dx
  ymod=ymod/dy
  !write(*,*)i,j,x,xdm,ymod,dy,'x'
  !funcf=-(Phi(i+1,j)-Phi(i,j))/dx
  !funcf=-(Phi(i+1,j)-Phi(i-1,j))/dx*0.5d0
   !funcf=(-(Phi(i+1,j)-Phi(i,j))/dx -(Phi(i+1,j+1)-Phi(i,j+1))/dx)*0.5d0
  delx1=-(Phi(i+1,j)-Phi(i-1,j))*0.5d0/dx * (1.d0-ymod) -(Phi(i+1,j+1)-Phi(i-1,j+1))*0.5d0/dx * (ymod)
  !delx1=-(Phi(i+1,j)-Phi(i-1,j))/dx*0.5d0
  !delx2=-(Phi(i+2,j)-Phi(i  ,j))/dx*0.5d0
  delx2=-(Phi(i+2,j)-Phi(i  ,j))*0.5d0/dx * (1.d0-ymod) -(Phi(i+2,j+1)-Phi(i  ,j+1))*0.5d0/dx * (ymod)
  delx=delx1* (1.d0-xmod)+delx2* xmod
  funcf=delx
  !funcf=-(Phi(i+1,j)-Phi(i,j))/dx * ymod/dy -(Phi(i+1,j+1)-Phi(i,j+1))/dx * (1.d0-ymod)/dy
  !funcf=-(Phi(i+1,j)-Phi(i,j))/dx * (1.d0-ymod) -(Phi(i+1,j+1)-Phi(i,j+1))/dx * (ymod)
  !write(*,*)funcf,'x'
  endif
  if(mode==2) then
  !funcf=-y!9.8d0!(1.0-x**2)
  i=int(xdm/dx)
  j=int(x/dy)
  ymod=x-j*dy
  xmod=xdm-i*dx
  xmod=xmod/dx
  ymod=ymod/dy
  !write(*,*)i,j,x,xmod/dx,'y'
  !funcf=-(Phi(i,j+1)-Phi(i,j))/dx
  !funcf=-(Phi(i,j+1)-Phi(i,j-1))/dy*0.5d0
   !funcf=(-(Phi(i,j+1)-Phi(i,j))/dx-(Phi(i+1,j+1)-Phi(i+1,j))/dx)*0.5d0
  delx1=-(Phi(i,j+1)-Phi(i,j-1))*0.5d0/dy * (1.d0-xmod) -(Phi(i+1,j+1)-Phi(i+1,j-1))*0.5d0/dy * (xmod)
  delx2=-(Phi(i,j+2)-Phi(i,j  ))*0.5d0/dy * (1.d0-xmod) -(Phi(i+1,j+2)-Phi(i+1,j  ))*0.5d0/dy * (xmod)
  delx=delx1* (1.d0-ymod)+delx2* ymod
  funcf=delx
  !funcf=-(Phi(i,j+1)-Phi(i,j))/dy* xmod/dx-(Phi(i+1,j+1)-Phi(i+1,j))/dy * (1.d0-xmod)/dx
  !funcf=-(Phi(i,j+1)-Phi(i,j))/dy* (1.d0-xmod)-(Phi(i+1,j+1)-Phi(i+1,j))/dy * (xmod)
  !write(*,*)funcf,'y'
  endif
 
  return
end function funcf



function funcv(x,v,xdm,mode)
  double precision :: t,x,v,xdm
  double precision :: funcv
  integer :: mode
 
  !funcv=(1.0-x**2)
  funcv=v!(1.0-x**2)
 
  return
end function funcv


subroutine rng(v,x,xdm,dt,mode)
double precision :: t,x,v,dt
double precision :: t0,t1,t2,t3
double precision :: kx0,kx1,kx2,kx3
double precision :: kv0,kv1,kv2,kv3
double precision :: f0,f1,f2,f3
double precision :: v0,v1,v2,v3
double precision :: funcf,funcv,xdm
integer :: mode

   !t0=t
   kv0=v
   kx0=x
   f0=funcf(kx0,kv0,xdm,mode)
   v0=funcv(kx0,kv0,xdm,mode)

   !t1=t+dt/2.0
   kv1=v+f0*dt/2.0
   kx1=x+v0*dt/2.0
   f1=funcf(kx1,kv1,xdm,mode)
   v1=funcv(kx1,kv1,xdm,mode)

   !t2=t+dt/2.0
   kv2=v+f1*dt/2.0
   kx2=x+v1*dt/2.0
   f2=funcf(kx2,kv2,xdm,mode)
   v2=funcv(kx2,kv2,xdm,mode)

   !t3=t+dt
   kv3=v+f2*dt
   kx3=x+v2*dt
   f3=funcf(kx3,kv3,xdm,mode)
   v3=funcv(kx3,kv3,xdm,mode)

   v=v+(f0+f1*2.0+f2*2.0+f3)*dt/6.0
   x=x+(v0+v1*2.0+v2*2.0+v3)*dt/6.0
end subroutine rng

subroutine ini(rrsph3,dinit1)
use comvar
!implicit none
integer :: i,j
!double precision :: Phi(0:ndx+1,0:ndy+1)
double precision :: rsph3,dinit1,rrsph3,G=1.11142d-4, G4pi!=12.56637d0*G

G4pi=12.56637d0*G
do j = 0, ndy+1; do i = 0, ndx+1
   
   rsph3 =dsqrt( (dble(ndx/2-i)*dx)**2 + (dble(ndy/2-j)*dy)**2)

   if(rsph3 .le. rrsph3 ) then
      Phi(i,j)=G4pi/6.d0*dinit1*(rsph3)**2
   else
      Phi(i,j)=-G4pi/rsph3/3.d0*dinit1*(rrsph3)**3.d0+G4pi/2.d0*dinit1*(rrsph3)**2.d0
   end if
end do
!write(*,*)Phi(ndy/2,j),'Phi'
end do

end subroutine ini
