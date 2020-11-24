program main
    implicit none

    integer date_time(8)
    character(len=10) sys_time(3)
    
    integer i,j,k
    double precision a,b,c,dat

    a=34.d0
    b=21.d0
    c=0.26d0
    
    call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    !write(*, '(a, 1x, i6)') "milli second     : ", date_time(7)
    call print_time8()
    do k=1,1064
    do j=1,1064
    do i=1,1064
    dat=b+c
    enddo
    enddo
    enddo
    !write(*, '(a, 1x, i6)') "milli second     : ", date_time(7)
    call print_time8()
    !call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    !write(*, '(a, 1x, i6)') "milli second     : ", date_time(7)
    !call print_time8()
    do k=1,1064
    do j=1,1064
    do i=1,1064
    dat=b*c
    enddo
    enddo
    enddo
    !write(*, '(a, 1x, i6)') "milli second     : ", date_time(7)
    call print_time8()
    do k=1,1064
    do j=1,1064
    do i=1,1064
    dat=b*c+b*c+b*c+b*c+b*c+b*c+b*c+b*c+b*c+b*c+b*c+b*c+b*c
    enddo
    enddo
    enddo
    call print_time8()

    
    write(*,*) "date_and_time"

    !write(*, '(a)') "<< argument 1-3 >>"
    !write(*, '(a, 1x, a)') "year, month, date : ", sys_time(1)
    !write(*, '(a, 1x, a)') "time              : ", sys_time(2)
    !write(*, '(a, 1x, a)') "time zone         : ", sys_time(3)
    !write(*,*)

    !write(*, '(a)') "<< argument 4   array>>"
    !write(*, '(a, 1x, i6)') "year             : ", date_time(1)
    !write(*, '(a, 1x, i6)') "month            : ", date_time(2)
    !write(*, '(a, 1x, i6)') "date             : ", date_time(3)
    !write(*, '(a, 1x, i6)') "time zone[min.]  : ", date_time(4)
    !write(*, '(a, 1x, i6)') "hour             : ", date_time(5)
    !write(*, '(a, 1x, i6)') "minute           : ", date_time(6)
    !write(*, '(a, 1x, i6)') "second           : ", date_time(7)
    !write(*, '(a, 1x, i6)') "milli second     : ", date_time(8)
end program main

subroutine print_time8()
    integer(kind=8) :: icount, icount_rate, icount_max
integer :: ipre=0

    call system_clock(COUNT=icount, COUNT_RATE=icount_rate, &
        & COUNT_MAX=icount_max)

    write(*, '(a, 1x, i20)') "icount       ", icount
    !write(*, '(a, 1x, i20)') "icount_rate  ", icount_rate
    !write(*, '(a, 1x, i20)') "icount_max   ", icount_max
    write(*,*) icount-ipre
    ipre=icount
end subroutine print_time8
