!!!!!!!!!!!!
!個数密度1/1.27の加熱冷却が釣り合う密度を求める。koyama inutsuka 2002
!!!!!!!!!!!!
program main
  implicit none
  double precision :: n , p1 , p2 , kai1 , kai2 , kai3 , root1 , root2 , root3 , sisu1 , sisu2 , sisu3 , sisu4 , sisu5 , sisu6 , &
       kaikai1 , kaikai2 , pmid , sa  , pos , inv,nmid,n1,n2,p
  integer :: i

  inv = 1.0d0/1.27d0
  !n = inv
  p=2.5d5
  n1 = 1.0d0
  n2 = 1.d5
  do i = 1 , 50
     !pmid = (p1 + p2) * 0.5d0
     nmid = (n1 + n2) * 0.5d0
     sisu1 = dexp(-118400.0d0/(p/(n1+pos) + 1000.0d0))
     sisu2 = dexp(-92.0d0/(p/(n1+pos)))
     sisu3 = dexp(-118400.0d0/(p/(n2+pos) + 1000.0d0))
     sisu4 = dexp(-92.0d0/(p/(n2+pos)))
     sisu5 = dexp(-118400.0d0/(p/(nmid+pos) + 1000.0d0))
     sisu6 = dexp(-92.0d0/(p/(nmid+pos)))
     root1 = dsqrt(p/(n1+pos))
     root2 = dsqrt(p/(n2+pos))
     root3 = dsqrt(p/(nmid+pos))
     kai1 = 1 - (1.0d7 * sisu1 + 1.4d-2 * root1 *sisu2) * n1
     kai2 = 1 - (1.0d7 * sisu3 + 1.4d-2 * root2 *sisu4) * n2
     kai3 = 1 - (1.0d7 * sisu5 + 1.4d-2 * root3 *sisu6) * nmid

     sa = dabs(n1 - n2)
     if(sa < 1.0d-6) goto 1000



     kaikai1 = kai1 * kai3
     kaikai2 = kai2 * kai3

     !write(*,*) kaikai1 , kaikai2

     if(kaikai1 < 0) then
        n2 = nmid
     else if(kaikai2 < 0) then
        n1 = nmid
     else
        write(*,*) 'please stop'
     end if
     if(i == 50) then
        write(*,*) 'itrmax'
     end if

  end do

1000 continue
  write(* , *) nmid , p


end program main
