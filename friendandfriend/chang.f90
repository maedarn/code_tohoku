
  t(:)=rl(:)
  m=1
  do i=1,max
     if(t(i)==i) then
        a(i)=m
        m=m+1
     end if
  end do
  write(*,*)'separate region number = ' , m


