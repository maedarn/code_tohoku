program main
  implicit none
  integer i,j,in,out,ratio,mode
  double precision , allocatable :: t(:) , tn(:)
  double precision dt,tt

  !*******************
  in = 91
  out = 910
  mode=2
  !*******************

  ALLOCATE(t(0:in),tn(0:out/mode))

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
  dt = 2.5d-2*mode
  ratio = out/in/mode
  tn(0) = 0.0d0
  do  j = 1 , out/mode

     !do i = 1 , ratio-1
     !   tn( j * 10 + i ) = t(j)  + dt
     !end do
     tn(j) = tn(j-1) + dt



  end do

  open(unit=151,file='timein.DAT')

  do i = 0 , out/mode
     write(151,*) tn(i)
  end do

  close(151)

end program main
