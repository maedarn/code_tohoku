program main
  implicit none
  integer :: i,j,k,isw,jsw,times1=0,times2=0
  integer ifl,li,aa

  li=0
  do jsw=2,1,-1
     isw=jsw
     do k=1,4
        do j=1,4
           do i=isw,4,2
              ifl=li*int(1/i)
              write(*,*) i,j,k,isw,jsw,ifl
           end do
           isw=3-isw
           times1=times1+1
        end do
        isw=3-isw
        times2=times2+1
     end do
     write(*,*)'----', times1,times2,aa*li,'----'
  end do
  !write(*,*) '----',times1,times2,'----'
end program main

