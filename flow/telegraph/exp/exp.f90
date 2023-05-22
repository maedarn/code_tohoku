program lowexp2
implicit none
real(8) :: x, y, ans
integer :: i, j, n
write(*,*) 'Please input number n:'
read(*,*) n
!write(*,*) 'Please input number x'
!read(*,*) x
x=-2.d0*0.5d0/0.2d0 * 0.5d0 * 0.5d0
ans = 1.d0
y = 1.d0
j = 1.d0
!x=-1.d-2
do i=1,n
 y = y * x/dble(i)
 ans = ans + y
end do
write(*,*) dexp(x), ans,x
stop
end program lowexp2
