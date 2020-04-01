program main
  implicit none
  integer nsisu,smesh,sloop,nb,val,i,j
  double precision,allocatable,dimension(:) :: lgM
  double precision,allocatable,dimension(:,:) :: U
  integer,allocatable,dimension(:) :: numM

  nb=13
  val=22
  sloop=6
  smesh=2

  allocate(lgM(0:sloop),numM(sloop))
  allocate(U(nb,val))

  nsisu=2
  lgM(:)=0.d0
  numM(:)=0
  U(:,:)=0.d0
  lgM(0)=10.0d0**nsisu
  do j=1,sloop
     nsisu=2
     nsisu =nsisu + 1.d0/dble(smesh)*dble(j)
     lgM(j) = 10.0d0**nsisu
     write(*,*) lgM(j)
  end do

  open(110,file='mdivhighden098001.DAT',FORM='FORMATTED')
  do j=1,nb
     read(110,*) (U(j,i),i=1,val)
  end do
  close(110)

  do i=1,nb
     do j=1,sloop
        if((lgM(j-1) < U(i,2)).and.(U(i,2) < lgM(j))) then
           numM(j)=numM(j)+1
        end if
     end do
  end do

  open(120,file='musfn-new',FORM='FORMATTED')
  do j=1,sloop
     write(120,*) 10.d0**(0.5d0*(dlog10(lgM(j-1))+dlog10(lgM(j)))),numM(j)
  end do
  close(120)

end program main
