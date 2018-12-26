program main
  implicit none
  double precision , allocatable , dimension(:,:) :: u
  integer n,time,step,i,j,val
  character(3) name
  character(5) name2

  n=128
  val=3

  allocate(u(1:n,1:val))

  do i=0,24999
     if(mod(i,100)==0) then
        write(name,'(i3.3)') i/100
        write(name2,'(i5.5)') i
        open(120,file='/Users/maeda/Desktop/kaiseki/testcode4/phi'//name2//'.dat')
        do j=1,n
           read(120,*) u(j,1),u(j,2),u(j,3)
        end do
        close(120)
        open(130,file='/Users/maeda/Desktop/kaiseki/testcode4/PHI'//name//'.dat')
        do j=1,n
           write(130,*) u(j,1),u(j,2),u(j,3)
        end do
        close(130)
     end if
  end do
end program main
