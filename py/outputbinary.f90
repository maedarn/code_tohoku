program main
  implicit none
  real(4) , allocatable , dimension(:,:) :: V
  real(4) , allocatable , dimension(:,:,:) :: Un,Un2
  real(8) , allocatable , dimension(:,:,:,:) :: U
  integer i,j,k,n

  !allocate(V(3,3))
  allocate(Un(4,4,4),Un2(4,4,4))

  do k=1,4
     do j=1,4
        do i=1,4
           !V(i,j)=10*i+j
           Un(i,j,k)=100*i+10*j+k
           Un2(i,j,k)=100*i+10*j+k+1000
        end do
     end do
  end do

!多分 open(11, 云々 ,access='stream') でデータべた書きにしないと駄目なのでは？
!unformatted だけだとバイナリで書き出すがレコード区切りが入る。 stream は Fortran2003 以降の属性だが、昔風には direct access 形式でもまぁ普通は行けなくもない。
  open(unit=8,file='txt.dat',FORM='UNFORMATTED' &
       !,access='stream',convert='little_endian')!convert=”big_endian”
       ,convert='little_endian')!convert=”big_endian”
  do k=1,4
     do j=1,4
        !do i=1,3
           !   write(8) V(i,j)
           !   write(8) (V(i,j),i=1,3)
           !   write(*,*) V(i,j)
           write(8) (Un(i,j,k),Un2(i,j,k),i=1,4)
           !write(8) Un(i,j,k),Un2(i,j,k)
        !end do
     end do
  end do

  !deallocate(V)
  deallocate(Un,Un2)
end program main
