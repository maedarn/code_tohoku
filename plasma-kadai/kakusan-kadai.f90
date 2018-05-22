program main
  implicit none
  character filename*128,filename2*128
  integer :: j , i , istep , nstep = 160
  integer , parameter :: n1=100 ,n2=100
  real(8) ::  d1 , d2 , dx1=1.0d-2 , dx2=1.0d-2 , x(0:n1) , y(n2) , dt = 1.d-4 , a = 0.5d0+0.02d0
  real(8) ::  phi(0:n1) , phi2(0:n1)=0.0d0
  integer k
  real(8) :: pi , rhs=0.0d0 ,er=0.0d0 ,er0=1.0d-6
  pi = acos(-1.0d0)
  d1 = a*dt/dx1/dx1

  x(:)=0.0d0
  do k = 1,n1
     x(k) = k*dx1
  end do

  !初期条件
  do i= 1 , 50
     phi(i) = 10.0d0*dsin(pi*dble(i)*0.01d0)+5*dsin(7.0d0*pi*dble(i)*0.01d0)+0.001d0*dsin(pi*dble(i)*0.5d0)
  end do
  do i= 51 , n1-1
     phi(i) = 10.0d0*dsin(pi*dble(i)*0.01d0)+5*dsin(15.0d0*pi*dble(i)*0.01d0)+0.001d0*dsin(pi*dble(i)*0.5d0)
  end do

  write (filename2, '("kakusan-3", i2.2, ".dat")') 0
  open (16, file=filename2, status='replace')
  !do j=1,n2
  do i=0,n1
     write(16,'(3e12.4)') x(i),phi(i)
     !write(16,*) ''
  end do
  !end do
  close (16)

  !ディリクレ境界
  phi(0) = 0.0d0
  phi(n1) = 0.0d0
  !do k = 1, n2
   !  phi(k) = 1 * sin(pi * k / n2)
  !end do
  !ノイマン境界
  !do k = 1, n2
   !  phi(n1,k) = phi(n1-1,k) !0.0d0
  !end do

  do istep = 1 , nstep
     er = 0.0d0
     rhs = 0.0d0
     !do j = 2 , n2-1
        do i = 1 , n1-1
           phi2(i) = phi(i) &
                + d1 * (phi(i-1) - 2.0d0*phi(i) + phi(i+1))
        end do
     !end do
     !ノイマン境界条件
     !do k = 1, n2
      !  phi(n1,k) = phi(n1-1,k)
     !end do
     do i = 1 , n1-1
        !do j = 2 , n2-1
           rhs = (phi2(i) - phi(i))**2
           er = er + rhs
        !enddo
     enddo
     phi(1:n1-1) = phi2(1:n1-1)
     write(*,*) 'ite , er = ' , istep , er

     if(mod(istep,20)==0) then
        write (filename, '("kakusan-3", i2.2, ".dat")') istep/20
        open (17, file=filename, status='replace')
        !do j=1,n2
           do i=0,n1
              write(17,'(3e12.4)') x(i),phi(i)
              !write(17,*) ''
           end do
        !end do
        close (17)
     endif
     !if(er<er0) exit
  end do
end program main
