program IMF
implicit none
integer :: nmesh,cnc1,cnc2,nran,idum
integer :: i
double precision :: mmax,mmin,nsisu,scale,slide,dmy
double precision :: M1,M2,coefflow,coeffmid,coeffhig,pwllow,pwlmid,pwlhig,coefflow2,coeffmid2,coeffhig2
double precision, allocatable, dimension(:) :: M,lgM,dN,N,lgN,Ncum,Mran,ranini,ts1,ts2,ts3,ts4
double precision :: dm,lgmmax,lgmmin,ddN,Ntot
character(21) :: dir='/Users/maeda/Desktop/'

!-------paramerter------
nmesh = 1000
nran=10000000
mmin = 0.01d0
mmax = 150.d0
Ntot = 100.d0
slide=5.d0
M1=0.08d0
M2=0.5d0
pwllow=-0.3d0
pwlmid=-1.3d0
pwlhig=-2.3d0
coeffhig=1.d0
coeffmid=(coeffhig*(M2)**pwlhig)/((M2)**pwlmid)
coefflow=(coeffmid*(M1)**pwlmid)/((M1)**pwllow)
coefflow2=1.d0
coeffmid2=(coefflow2*(M1)**(pwllow))/((M1)**(pwlmid))
coeffhig2=(coeffmid2*(M2)**(pwlmid))/((M2)**(pwlhig))
!-------paramerter------

allocate(M(0:nmesh),lgM(0:nmesh),dN(0:nmesh),N(0:nmesh),lgN(0:nmesh),Ncum(0:nmesh),Mran(1:nran),ranini(1:nran))
allocate(ts1(1:nran),ts2(1:nran),ts3(1:nran),ts4(1:nran))

lgmmax = dlog10(mmax)
lgmmin = dlog10(mmin)
dm = (lgmmax-lgmmin)/dble(nmesh)
cnc1=int((dlog10(M1)-dlog10(mmin))/dm)
cnc2=int((dlog10(M2)-dlog10(mmin))/dm)

write(*,*)lgmmax,lgmmin,dm,cnc1,cnc2

lgM(0)=lgmmin
do i= 1,nmesh
lgM(i)=lgM(i-1)+dm
!write(*,*)lgM(i)
enddo

do i= 0,nmesh
M(i)=10.d0**(lgM(i))
!write(*,*)M(i)
enddo


write(*,*) cnc1
do i= 0,cnc1
dN(i)=coefflow*M(i)**(pwllow+1.d0)
enddo

write(*,*) cnc2
do i= cnc1,cnc2
dN(i)=coeffmid*M(i)**(pwlmid+1.d0)
enddo

do i= cnc2,nmesh
dN(i)=coeffhig*M(i)**(pwlhig+1.d0)
enddo


Ntot=coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))+&
     coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
     coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))



Ncum(:)=0.d0
Mran(:)=0.d0

do i= nmesh,cnc2,-1
Ncum(i)=coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M(i))**(pwlhig+1.d0)/(pwlhig+1.d0))

!Mran(i)=((coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-Ncum(i))*(pwlhig+1.d0)/coeffhig)**(1.d0/(pwlhig+1.d0))
enddo

do i= cnc2,cnc1,-1
Ncum(i)=coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M(i))**(pwlmid+1.d0)/(pwlmid+1.d0))+&
        coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))

!Mran(i)=((coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
!          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
!         -Ncum(i))*(pwlmid+1.d0)/coeffmid)**(1.d0/(pwlmid+1.d0))
enddo

do i= cnc1,0,-1
Ncum(i)=coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(M(i))**(pwllow+1.d0)/(pwllow+1.d0))+&
        coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
        coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))

!Mran(i)=((coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
!          coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))-coeffmid*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
!          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-coeffhig*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
!         -Ncum(i))/coefflow*(pwllow+1.d0))**(1.d0/(pwllow+1.d0))
enddo

do i= nmesh,0,-1
Ncum(i)=Ncum(i)/Ntot
enddo


!goto 2021
!--------------
Mran(:)=0.d0
idum=1
do i=1,nran
call ran0(Mran(i),idum)
enddo

ranini(:)=Mran(:)!*Ntot


write(*,*)Ncum(cnc1),Ncum(cnc2),M(cnc1),M(cnc2)
do i=1,nran
if (ranini(i)>Ncum(cnc1)) then
!Mran(i)=((coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
!          coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
!          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
!          -ranini(i))*(pwllow+1.d0)/coefflow)**(1.d0/(pwllow+1.d0))

Mran(i)=((coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
          coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))-coeffmid*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-coeffhig*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*ranini(i))/coefflow*(pwllow+1.d0))**(1.d0/(pwllow+1.d0))

!write(*,*)Mran(i),ranini(i),Ncum(cnc1),Ncum(cnc2),'1'
else if ((Ncum(cnc1)>ranini(i)).and.(ranini(i)>Ncum(cnc2))) then
Mran(i)=((coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*ranini(i))*(pwlmid+1.d0)/coeffmid)**(1.d0/(pwlmid+1.d0))
!write(*,*)Mran(i),ranini(i),Ncum(cnc1),Ncum(cnc2),'2'
else if (Ncum(cnc2)>ranini(i)) then
Mran(i)=((coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-Ntot*ranini(i))*(pwlhig+1.d0)/coeffhig)**(1.d0/(pwlhig+1.d0))
!write(*,*)Mran(i),ranini(i),Ncum(cnc1),Ncum(cnc2),'3'
endif
enddo
!--------------

goto 2021
!--------------
Ntot=coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))+&
     coeffmid2*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
     coeffhig2*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncum(:)=0.d0

do i= 0,cnc1
Ncum(i)=coefflow2*((M(i))**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))
enddo
do i= cnc1,cnc2
Ncum(i)=coeffmid2*((M(i))**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
        coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))
enddo

do i= cnc2,nmesh
Ncum(i)=coeffhig2*((M(i))**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))+&
        coeffmid2*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
        coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))
enddo


do i= 0,nmesh
Ncum(i)=Ncum(i)/Ntot
enddo


!do i=1,nran
!if (ranini(i)<Ncum(cnc1)) then
!Mran(i)=(((ranini(i)+coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coefflow2)*(pwllow+1.d0))**(1.d0/(pwllow+1.d0))
!else if ((Ncum(cnc2)>ranini(i)).and.(ranini(i)>Ncum(cnc1))) then
!Mran(i)=(((ranini(i)+coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))-&
!          coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
!          coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coeffmid2)/(pwlmid+1.d0))**(1.d0/(pwlmid+1.d0))
!else if (Ncum(cnc2)<ranini(i)) then
!Mran(i)=(((ranini(i)+coeffhig2*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))-&
!           coeffmid2*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))-&
!           coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coeffhig2)&
!           *(pwlhig+1.d0))**(1.d0/(pwlhig+1.d0))
!endif
!enddo

ts1(:)=0.d0
ts2(:)=0.d0
ts3(:)=0.d0
ts4(:)=0.d0

do i=1,nran
if (Ncum(i)<Ncum(cnc1)) then
Mran(i)=(((Ncum(i)+coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coefflow2)*(pwllow+1.d0))!**(1.d0/(pwllow+1.d0))
else if ((Ncum(cnc2)>Ncum(i)).and.(Ncum(i)>Ncum(cnc1))) then
!Mran(i)=(((Ncum(i)+coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))-&
!          coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
!          coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coeffmid2)/(pwlmid+1.d0))**(1.d0/(pwlmid+1.d0))
Mran(i)=(((Ncum(i)+coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))-&
coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coeffmid2)/(pwlmid+1.d0))**(1.d0/(pwlmid+1.d0))

ts1(i)=coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))
ts2(i)=-coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))!+coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0))
ts3(i)=coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0))
ts4(i)=+coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))-&
coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0))

else if (Ncum(cnc2)<Ncum(i)) then
!Mran(i)=(((Ncum(i)+coeffhig2*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))-&
!           coeffmid2*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+coeffmid2*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))-&
!           coefflow2*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+coefflow2*((mmin)**(pwllow+1.d0)/(pwllow+1.d0)))/coeffhig2)&
!           *(pwlhig+1.d0))**(1.d0/(pwlhig+1.d0))
Mran(i)=0.d0
endif
enddo
!-------------
2021 continue

open(20,file = dir//'IMF.dat')
do i=0,nmesh
write(20,*) 10.d0**lgM(i),',',Ncum(i),',',dN(i)
end do
close(20)

open(21,file = dir//'Mran.dat')
do i=1,nran
write(21,*) Mran(i),',',ranini(i)!,Ncum(i)!,ts1(i),ts2(i),ts3(i),ts4(i)!ranini(i)
end do
close(21)

deallocate(M,lgM,dN,N,lgN,Ncum,Mran,ranini)
deallocate(ts1,ts2,ts3,ts4)
end program IMF

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
