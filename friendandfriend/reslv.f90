module common
  implicit none
  integer,parameter :: nx=512,ny=512,nz=512
  integer,parameter :: rank=nx*ny*nz/4
  integer,dimension(1:rank) :: rl_table,next_table,tail_table,T
end module common

program main
  use common
  integer lab,i,ik
  rl_table(1)=1
  rl_table(2)=1
  rl_table(3)=3

  next_table(1)=2
  next_table(2)=-1
  next_table(3)=-1

  tail_table(1)=2
  tail_table(2)=2
  tail_table(3)=3
  !lab=1
  !do i=1,2
  !   rl_table(lab)=lab
  !   next_table(lab)= -1
  !   tail_table(lab)= lab
  !   lab=lab+1
  !end do
  call resolve(2,3)
  open(190,file='/Users/maeda/Desktop/kaiseki/image/kakunin2.DAT',FORM='FORMATTED')
  do ik=1,10
     write(190,*) ik,rl_table(ik),next_table(ik),tail_table(ik)
  end do
  close(190)
end program main

subroutine resolve(a,b)
  !implicit none
  use common
  implicit none
  integer u,v,i,ik,a,b
  !write(*,*) 'ok1',rl_table(v)
  v=rl_table(a)
  u=rl_table(b)
  write(*,*) 'ok2'
  !write(*,*)'aa',u,v,next_table(tail_table(u)),next_table(tail_table(v)),tail_table(u)
  next_table(tail_table(u))=v
  tail_table(u)=tail_table(v)
  !write(*,*)tail_table(u),tail_table(v),next_table(v),next_table(next_table(v))
  i=v
  write(*,*)'aa',u,v,next_table(tail_table(u)),tail_table(u)
  do while(i.ne.-1)
     !write(*,*) u,v,i,next_table(i),next_table(v)!,next_table(31)
     rl_table(i)=u
     i=next_table(i)
     write(*,*) i
  end do
  !write(*,*)'bb'
end subroutine resolve
