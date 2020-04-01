program main
  implicit none
  integer i

  do i=1,20
     open(10,file='number.dat'number.dat)
     write(10,*) i
     close(10)
  end do
end program main
