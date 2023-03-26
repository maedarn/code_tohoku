PROGRAM MAIN
  implicit none
  doubleprecision :: phase1,pi
  INTEGER :: i,j,idum
  double precision, allocatable, dimension(:) :: Mran
  pi=3.1415926535d0
  idum=1
  allocate(Mran(1:4))
  call ran0(phase1,idum); phase1=2.d0*pi*phase1
  write(*,*)'fst',idum
  do j=1,2
     do i=1,2
        !idum=1+idum
        write(*,*)'pre',idum
        call ran0(phase1,idum); !phase1=2.d0*pi*phase1
        write(*,*)'mid',idum
        call ran0(Mran(i),idum); !phase1=2.d0*pi*phase1
        write(*,*)'pst',idum
        write(*,*) Mran(i), phase1
     enddo
  enddo

  deallocate(Mran)
END PROGRAM

SUBROUTINE ran0(ran,idum)
INTEGER idum,IA,IM,IQ,IR,MASK
doubleprecision :: ran,AM
PARAMETER (IA=16807,IM=2147483647,AM=1./IM, IQ=127773,IR=2836,MASK=123459876)
INTEGER k
!idum=ieor(idum,MASK)
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if(idum.lt.0) idum=idum+IM
ran=AM*idum
!idum=ieor(idum,MASK)
END SUBROUTINE ran0
