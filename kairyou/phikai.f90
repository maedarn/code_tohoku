program main
  implicit none
  integer i,j,len
  real(8) ,allocatable , dimension(:) :: x , y , Phi
  real(8) G,R,pi,Rcen,Lbox,dx,rho,G4pi
  !***********************
  len=256
  R=40.d0
  Rcen=50.d0
  Lbox=1.d2
  !***********************
  G = 1.11142d-4
  G4pi=12.56637d0*G
  pi = 3.1415926535d0
  rho = 3.0d0/G4pi/4.d1/4.d1

  ALLOCATE(x(1:len))
  ALLOCATE(y(1:len))
  ALLOCATE(Phi(1:len))
  dx = Lbox/dble(len)
  x(1) = dx/2.0d0
  do i = 2 , len
     x(i) = x(i-1) + dx
  end do
  y(:) = x(:)

  do i=1,len
     if(dabs(x(i) - Rcen) .gt. R ) then
        Phi(i) = - G * rho * 4.0d0/3.0d0 * pi * R * R * R /dabs(x(i) - Rcen)
     else
        Phi(i) = 2.0d0 * pi * G * rho / 3.0d0 * (x(i) - Rcen)**2 &
             - 2.0d0 * pi * G * rho * R**2
     end if
  end do

  open(unit=66,file='phikai1D.dat')
  do i=1,len
     write(66,*) sngl(x(i)) , sngl(Phi(i))
  end do
  close(66)
  !open(unit=67,file='phikai1D.dat')
  !do j=1,len
  !   do i=1,len
  !      write(67,*) sngl(x(i)) , sngl(x(i)) , sngl(Phi(i))
  !   end do
  !   write(67,*)
  !end do
  !close(67)
end program main
