program main
  implicit none
  double precision x, y, A ,theta,dx
  integer i
  CHARACTER(29) :: dir='/Users/maeda/Desktop/code/cy/'
  character(4) :: filenm='cl02'

  A=2.0d0
  dx=3.1414592d0 * 2.0d0 * 1.0d0 / 100.0d0

  theta=0.0d0

100 format(E19.10e3,E19.10e3)
  open(10,file=dir//filenm//'.dat')
  do i=1,100
     theta = theta + dx

     x=A*(theta - dsin(theta))
     y=A*(1- dcos(theta))


     if(mod(i,10)==0)  then
        write(*,*) x,y
     end if


     write(10,100) x , y

  end do

end program main
