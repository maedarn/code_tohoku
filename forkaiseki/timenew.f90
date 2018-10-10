program main
  implicit none
  integer i,j,in,out,ratio
  double precision , allocatable :: t(:) , tn(:)
  double precision dt,tt

  !*******************
  in = 20
  out = 200
  !*******************

  ALLOCATE(t(0:20),tn(0:200))

  open(unit=150,file='time.DAT')

  do i = 0 , in
     read(150,*) t(i)
     !write(*,*) t(i)
  end do

  close(150)

  !tt=0
  !dt = (t(1) - t(0)) * dble(in) / dble(out)
  !write(*,*) dt
  !do k = 0 , in
  !do i = 0 , out
  !tn(i) = t(k) + dt
  !tt = tt + dt
  !k = int(tt)
  dt = 2.5d-2
  ratio = out/in
  tn(0) = 0.0d0
  do  j = 1 , out

     !do i = 1 , ratio-1
     !   tn( j * 10 + i ) = t(j)  + dt
     !end do
     tn(j) = tn(j-1) + dt



  end do

  open(unit=151,file='timein.DAT')

  do i = 0 , out
     write(151,*) tn(i)
  end do

  close(151)

end program main
