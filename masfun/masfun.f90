program main
  implicit none
  integer smesh,sloop,nb,val,i,j
  double precision nsisu,a
  double precision,allocatable,dimension(:) :: lgM
  double precision,allocatable,dimension(:,:) :: U
  integer,allocatable,dimension(:) :: numM
  integer,allocatable,dimension(:) :: LL
  

  nb=30
  val=2
  sloop=8
  smesh=2
  a=0.d0

  allocate(lgM(0:sloop),numM(sloop))
  allocate(U(nb,val))
  allocate(LL(nb))

  nsisu=2.d0
  lgM(:)=0.d0
  numM(:)=0
  U(:,:)=0.d0
  lgM(0)=10.0d0**nsisu
  do j=1,sloop
     nsisu=2.d0
     nsisu =nsisu + 1.d0/dble(smesh)*dble(j)
     lgM(j) = 10.0d0**nsisu
     write(*,*) lgM(j)
  end do

  open(110,file='/Users/maeda/Desktop/Dropbox/analysis/samplecnvK/N-Masfuncrrnew082001.DAT',FORM='FORMATTED')
  do j=1,nb
     read(110,*) U(j,1),LL(j)
     !write(*,*) U(j,i)
     U(j,2)=dble(LL(j))
  end do
  close(110)

  do i=1,nb
     do j=1,sloop
        if((lgM(j-1) < U(i,2)).and.(U(i,2) < lgM(j))) then
           numM(j)=numM(j)+1
        end if
     end do
  end do

  open(120,file='musfn-new4.dat',FORM='FORMATTED')
  do j=1,sloop
     write(120,*) 10.d0**(0.5d0*(dlog10(lgM(j-1))+dlog10(lgM(j)))),numM(j)
  end do
  close(120)

  open(130,file='musfn-new5-25.dat',FORM='FORMATTED')
  do j=1,sloop
     !write(130,*) 10.d0**(0.5d0*(dlog10(lgM(j-1))+dlog10(lgM(j)))),numM(j)
     write(130,*) lgM(j-1),numM(j)
     write(130,*) lgM(j),numM(j)
     write(130,*) lgM(j),a
    ! write(130,*) lgM(j),a
  end do
  close(130)

  open(140,file='musfn-new5-cum.dat',FORM='FORMATTED')
  do i=sloop-1,1,-1
     numM(i)=numM(i)+numM(i+1)
  end do
  do j=1,sloop
     !write(130,*) 10.d0**(0.5d0*(dlog10(lgM(j-1))+dlog10(lgM(j)))),numM(j)
     write(140,*) lgM(j-1),numM(j)
     write(140,*) lgM(j),numM(j)
     write(140,*) lgM(j),a
  end do
  close(140)

end program main
