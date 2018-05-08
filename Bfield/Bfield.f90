program main
  implicit none
  integer i,n,j,l
  double precision B0,B1,x,y,dx,dy,Bx,By,dBx,dBy,a,k,omega,t,pi,xpiku,ypiku
  CHARACTER(33) :: dir='/Users/maeda/Desktop/code/Bfield/'
  character(4) :: filenm='BF01'
  character filename*128
100 format(E19.10e3,E19.10e3)
  open(10,file=dir//filenm//'.dat')

  x=0.0d0
  y=0.0d0
  k=2.0d0
  a=0.0d0
  omega=1.0d0
  t=0.0d0
  B0=10.0d0
  B1=3.0d0
  dx=1.0d0/100.0d0
  dy=1.0d0/100.0d0
  xpiku=0.0d0
  ypiku=0.0d0
  pi=3.14d0

  do l=1,5
     xpiku=dsqrt(1.0d0*pi/k)+xpiku
     ypiku=dsqrt(1.0d0*pi/k)+ypiku
     write(10,100) xpiku,ypiku
     !write(10,*)
  end do

  do j= 1,5
     t=0.2*int(j)
     write (filename, '("BF",i1.1, ".dat")') j
     open (17, file=filename, status='replace')
     do n=1,10
        a=int(n)
        do i=1,100

           Bx=-B1*dcos(k*(x+y)-omega*t)
           By=B0+B1*dsin(k*(x+y)-omega*t)

           !dBx=-B1*k*dsin(k*(x+y))+1.d-10
           !dBy=B1*k*dcos(k*(x+y))+1.d-10

           x=Bx*dx+x
           y=By*dy+y

           write(17,100) x , y
        end do
        x=0.0d0 + a
        y=0.0d0
        write(17,*) 
     end do
     a=0.0d0
     x=0.0d0
     y=0.0d0
  end do
end program main

