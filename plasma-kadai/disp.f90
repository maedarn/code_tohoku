program main
implicit none
double precision :: pi, kappa, theta, dtheta, hosei, coeff ,phi2
double precision, allocatable, dimension(:) :: phi
integer i,n,step

pi=3.141592653589d0
!*******parameter*****
step=1000
kappa=1.8d0
coeff= -1.0d0
hosei= 0.d0*pi
!*********************

allocate(phi(step))
dtheta=2.d0*pi/dble(step)
theta=0.d0

do i = 1,step
theta=theta+dtheta
phi(i)=(kappa*(kappa-2.d0)*dsin(theta)-0.5d0*kappa*(kappa-1.d0)*dsin(2.d0*theta)) / &
(0.5d0*(kappa-1.d0)*(kappa-2.d0)-kappa*(kappa-2.d0)*dcos(theta)+0.5d0*kappa*(kappa-1.d0)*dcos(2.d0*theta))
end do

theta=0.d0
open (10, file='phi20exa.dat')
do i=1,step
   theta=theta+dtheta
   if((phi(i-1)*phi(i)<0.d0).and.(dabs(phi(i))>1.d0))then
   hosei=hosei+coeff*pi
   endif
   !theta=theta+hosei
   !phi2=(kappa*(kappa-2.d0)*dsin(theta)-0.5d0*kappa*(kappa-1.d0)*dsin(2.d0*theta)) / &
   !(0.5d0*(kappa-1.d0)*(kappa-2.d0)-kappa*(kappa-2.d0)*dcos(theta)+0.5d0*kappa*(kappa-1.d0)*dcos(2.d0*theta))
   write(10,*) phi(i),theta,datan(phi(i)),-(datan(phi(i))+hosei)/kappa/theta,&
   dtan(datan(phi(i))+hosei),phi(i)+hosei,datan(phi(i)+hosei)!,phi2
end do
close(10)

deallocate(phi)
end program
