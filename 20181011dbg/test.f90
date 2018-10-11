program test
  implicit none
  integer i,j
  do i=1,10
     j=2/(i-4)
     write(*,*) j
     !print *,j
  end do
  stop
end program test
