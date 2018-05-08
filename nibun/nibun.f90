program main
  implicit none
  double precision :: n , p1 , p2 , kai1 , kai2 , kai3 , root1 , root2 , root3 , sisu1 , sisu2 , sisu3 , sisu4 , sisu5 , sisu6 , &
       kaikai1 , kaikai2 , pmid , sa , nsisu , pos , inv
  integer :: i , j

  !character filename*6
  pos = 1d-10
  inv = 1.0d0 / 1.27d0
  n = 0.0d0
  do j = 1 , 100
     nsisu = -2.0d0
     nsisu =nsisu + 0.05d0*dble(j)
     n = 10.0d0**nsisu
     p1 = 1.0d0
     p2 = 1.d6
     write(*,*) n
     do i = 1 , 50
        pmid = (p1 + p2) * 0.5d0
        sisu1 = dexp(-118400.0d0/(p1/(n+pos) + 1000.0d0))
        sisu2 = dexp(-92.0d0/(p1/(n+pos)))
        sisu3 = dexp(-118400.0d0/(p2/(n+pos) + 1000.0d0))
        sisu4 = dexp(-92.0d0/(p2/(n+pos)))
        sisu5 = dexp(-118400.0d0/(pmid/(n+pos) + 1000.0d0))
        sisu6 = dexp(-92.0d0/(pmid/(n+pos)))
        root1 = dsqrt(p1/(n+pos))
        root2 = dsqrt(p2/(n+pos))
        root3 = dsqrt(pmid/(n+pos))
        kai1 = 1 - (1.0d7 * sisu1 + 1.4d-2 * root1 *sisu2) * n
        kai2 = 1 - (1.0d7 * sisu3 + 1.4d-2 * root2 *sisu4) * n
        kai3 = 1 - (1.0d7 * sisu5 + 1.4d-2 * root3 *sisu6) * n

        sa = dabs(p1 - p2)
        if(sa < 1.0d-6) goto 1000



        kaikai1 = kai1 * kai3
        kaikai2 = kai2 * kai3

        !write(*,*) kaikai1 , kaikai2

        if(kaikai1 < 0) then
           p2 = pmid
        else if(kaikai2 < 0) then
           p1 = pmid
        else
           write(*,*) 'please stop'
        end if
        if(i == 50) then
           write(*,*) 'itrmax'
        end if

     end do

1000 open(2 , file = 'KIcool.dat' , STATUS = 'UNKNOWN')
     write(2 , *) n , pmid
  end do
  close(2)
end program main
