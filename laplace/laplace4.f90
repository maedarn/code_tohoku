program main
  implicit none
  integer :: itrmax = 150 , j , i , itr
  integer , parameter :: n1=100 ,n2=100
  real(8) ::  c , d , omg , dx1=1.0d0 , dx2=1.0d0 , x(n1) , y(n2)
  real(8) ::  phi(1:n1,1:n2)
  integer k
  real(8) :: pi , h=0 , rhs=0 ,er=0 ,er0=1.0d-6
  pi = acos(-1.0d0)
  do k = 1,n1
     x(k) = k*dx1
  end do
  do k = 1,n2
       y(k) = k*dx2
  end do
  do k = 1, n2
     phi(1,k) = 1 * sin(pi * k / n2)
  end do
  do k = 1, n2
     phi(n1,k) = 0.0d0
  end do
  do k = 1, n1
     phi(k,1) = 0.0d0
  end do
  do k = 1, n1
     phi(k,n2) = 0.0d0
  end do
  phi(2:n1-1,2:n2-1) = 0.0d0
  c = -dx2**2 / (dx1**2 + dx2**2) * 0.5d0
  d = (dx1**2 / dx2**2)**2 * c
  omg = 2.0d0/(1.0d0 + sin(pi/(n2-1)))
  do itr = 1 , itrmax
     er = 0
     do j = 2 , n2 - 1
        do i = 2 , n1 - 1
           rhs = - c * (phi(i-1,j) + phi(i+1,j)) - d * (phi(i,j-1) + phi(i,j+1))
           h = (rhs - phi(i,j))**2
           phi(i,j) = phi(i,j) + omg * (rhs - phi(i,j))
           er = er + h
        end do
     end do
     write(*,*) 'ite , er = ' , itr , er
     if(er<er0) exit
  end do
  open(50,file='laplace5.dat')
  do j=1,n2
     do i=1,n1
        write(50,'(3e12.4)') x(i),y(j),phi(i,j)
     end do
     write(50,'(3e12.4)')
  end do
end program main
