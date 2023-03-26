PROGRAM MAIN
  implicit none
  doubleprecision :: pi,rhoc,dr
  DOUBLE PRECISION, parameter :: G=1.11142d-4
  DOUBLE PRECISION, parameter :: cs=0.2d0
  INTEGER :: i,j,loop
  double precision, allocatable, dimension(:) :: rho,M,lnrho,r
  character(21) :: dir='/Users/maeda/Desktop/'

  pi=3.1415926535d0
  loop=500
  rhoc=1.d4
  dr=0.05d0

  allocate(rho(0:loop),M(0:loop),lnrho(0:loop),r(0:loop))
  rho(0)=rhoc
  lnrho(0)=dlog(rhoc)
  M(0)=0.d0
  r(0)=0.d0
  do i=1,loop
  r(i)=r(i-1)+dr
  M(i)=4.d0*pi*r(i)*r(i)*dr*rho(i-1)+M(i-1)
  lnrho(i)=-G*M(i-1)/cs/cs/r(i)/r(i)*dr+lnrho(i-1)
  rho(i)=dexp(lnrho(i))
  enddo

  open(20,file = dir//'BE-2.dat')
  do i=0,loop
  write(20,*) r(i),',',M(i),',',rho(i),',',lnrho(i)
  end do
  close(20)
  
  deallocate(rho,M,lnrho,r)
END PROGRAM

