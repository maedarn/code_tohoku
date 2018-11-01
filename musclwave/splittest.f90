program main
  implicit none
  integer i
  DOUBLE PRECISION , dimension(1:128) :: x,Phi, Phi1step ,cal,Phi2step,Phi3step
  DOUBLE PRECISION , dimension(1:128) :: Phidt,Phigrd,Phiexa,Phi2dt,Phi3dt
  double precision :: cg = 1.0d0,dx,dt,a
  !real(4) a

  dx = 100.0/128.d0
  dt=7.812d-2
  open(18,file='/Users/maeda/Desktop/kaiseki/testcode2/phi00800.dat')
  open(19,file='/Users/maeda/Desktop/kaiseki/testcode2/phi00799.dat')
  open(21,file='/Users/maeda/Desktop/kaiseki/testcode2/phi00798.dat')
  open(22,file='/Users/maeda/Desktop/kaiseki/testcode2/phi00797.dat')
  !write(*,*) 'aa'
  do i=1,128
     read(18,*) a , Phi(i) , Phi1step(i)
  end do
  do i=1,128
     read(19,*) a , Phidt(i),Phi2step(i)
  end do
  do i=1,128
     read(21,*) a , Phi2dt(i),Phi3step(i)
  end do
  do i=1,128
     read(22,*) a , Phi3dt(i)
  end do
  close(18)
  close(19)
  close(21)
  close(22)
  do i=3,127
     cal(i)= ((Phi(i) - Phidt(i)) / cg / dt - (Phi(i+1) - Phi(i-1)) * 0.5d0 /dx)
     !cal(i)= ((3.0d0*Phidt(i) -4.0d0* Phi2dt(i) + Phi3dt(i) )*0.5d0 / cg / dt +&
     !     (3.0d0*Phidt(i) - 4.0d0*Phidt(i+1)+Phidt(i+2)) * 0.5d0 /dx)
     !cal(i)= ((3.0d0*Phidt(i) -4.0d0* Phi2dt(i) + Phi3dt(i) )*0.5d0 / cg / dt -&
     !(3.0d0*Phidt(i) - 4.0d0*Phidt(i-1)+Phidt(i-2)) * 0.5d0 /dx)
     !write(201,*) sngl(x(i)) , Phi1step(i)
     !write(201,*) i , Phi1step(i),Phiv(i) , Phidt(i),Phiv(i+1) , Phiv(i-1),Phiv(i+1) - Phiv(i-1),Phiv(i) - Phidt(i)
  end do
  open(20,file='/Users/maeda/Desktop/kaiseki/testcode2/cal.dat')
  do i=2,127
     write(20,*) i,cal(i),Phi1step(i),Phi(i),Phidt(i),Phi2dt(i),Phi(i) - Phidt(i),&
          Phi(i+1) - Phi(i-1),Phi2step(i),Phi3step(i),((Phi(i) - Phidt(i)) / cg / dt)
  end do
  close(20)
  write(*,*) Phi(1)
end program main
