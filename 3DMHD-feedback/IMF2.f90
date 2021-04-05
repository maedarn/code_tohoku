program IMF
implicit none
integer :: nmesh,cnc
integer :: i
double precision :: mmax,mmin,nsisu,scale,slide
double precision, allocatable, dimension(:) :: M,lgM,dN,N,lgN,Ncum
double precision :: dm,lgmmax,lgmmin,ddN,Ntot
character(21) :: dir='/Users/maeda/Desktop/'

!-------paramerter------
nmesh = 1000
mmin = 0.1d0
mmax = 100.d0
Ntot = 100.d0
slide=5.d0
dm=0.1d0
!-------paramerter------


allocate(M(0:nmesh),lgM(0:nmesh),dN(0:nmesh),N(0:nmesh),lgN(0:nmesh),Ncum(0:nmesh))


lgmmax = dlog10(mmax)
lgmmin = dlog10(mmin)
!dm = (lgmmax-lgmmin)/dble(nmesh)
!cnc=int((dlog10(1.d0)-dlog10(mmin))/dm)
cnc=10

write(*,*)lgmmax,lgmmin,dm,cnc

!nsisu=lgmmin
!lgM(:)=0.d0
!lgM(0)=10.0d0**nsisu
!do i=1,nmesh
!   nsisu=lgmmin
!   nsisu =nsisu+ 1.d0/dble(lgmmax-lgmmin)*dble(i)
!   lgM(i) = 10.0d0**nsisu
!   write(*,*) lgM(i),nsisu
!end do

!lgM(0)=lgmmin
!do i= 1,nmesh
!lgM(i)=lgM(i-1)+dm
!write(*,*)lgM(i)
!enddo

M(0)=0.1d0
do i= 1,nmesh
!M(i)=10.d0**(lgM(i))
M(i)=M(i-1)+dm
!write(*,*)M(i)
write(*,*) M(i),i
enddo



ddN=dexp(-(dlog10(1.d0/0.079d0))**2.d0 / 2.d0/0.69d0/0.69d0)/dlog(10.d0)-(1.d0)**(-1.3d0)/dlog(10.d0)

write(*,*) cnc
do i= 0,cnc
dN(i)=dexp(-(dlog10(M(i)/0.079d0))**2.d0 / 2.d0/0.69d0/0.69d0)/M(i)/dlog(10.d0)
enddo

do i= cnc,nmesh
!dN(i)=(M(i))**(-1.3d0)+ddN
dN(i)=(M(i))**(-1.3d0)/M(i)/dlog(10.d0)!+ddN
enddo

!N(0)=dN(0)
do i= 1,nmesh
N(i)=dN(i)*dm!+N(i-1)
!N(i)=dN(i)*dm!+N(i-1)
!write(*,*)M(i)
enddo

!Ncum(:)=0.d0
!Ncum(nmesh)=10.d0**(N(nmesh))
!do i= nmesh-1,0,-1
!Ncum(i)=Ncum(i+1)+10.d0**(N(i))
!enddo
!write(*,*)scale,ddN
!N(:)=N(:)*Ntot/scale
!N(:)=N(:)+slide
!Ncum(:)=Ncum(:)/Ncum(0)
!Ntot=Ncum(0)
!do i= nmesh,0,-1
!Ncum(i)=Ncum(i)/Ntot
!enddo
write(*,*) 'here'

open(20,file = dir//'IMF.dat')
do i=0,nmesh
write(*,*) 'here', M(i),',',dN(i),',',N(i)
write(20,*) M(i),',',dN(i),',',N(i)
end do
close(20)

deallocate(M,lgM,dN,N,lgN,Ncum)
end program IMF
