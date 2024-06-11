module comvar
real(4), dimension(:,:,:), allocatable :: Mtl
integer, parameter :: ifile = 30, imtrx=189,num_lin=30*(3+32)/2
integer, dimension(1:ifile) :: in_mtl
double precision, dimension(1:ifile) :: Zmetals
double precision :: Zi=1.d0,alphaB,ni,ne,kb=1.38d-16,Zsolar=1.d0,xHe=0.1d0
double precision :: Sigmo,Tth1=2.d4,Tth2=3.5d4,asigmo=1.d1
end module comvar

program main
use comvar
implicit none
integer, parameter :: ix = 100
integer :: i,j,k,jm,j_lin
double precision, parameter :: Tmin=1.d4,Tmax=1.d8
double precision :: dlogT,Tmp1,dlogT1
double precision, dimension(:), allocatable :: Tmp,Rho
double precision, dimension(:), allocatable :: Lambda,Lambdaff,Lambdarr,Lambdameta,LambdaHe
double precision, dimension(:,:), allocatable :: Lambda_CIE


allocate(Tmp(0:ix),Rho(0:ix))
allocate(Lambda(0:ix),Lambdaff(0:ix),Lambdarr(0:ix),Lambdameta(0:ix),LambdaHe(0:ix))
allocate(Lambda_CIE(0:ix,1:num_lin))
allocate(Mtl(1:ifile,1:ifile+3,1:imtrx))
Mtl(:,:,:)=0.e0
Lambda_CIE(:,:)=0.d0
Lambdaff(:)=0.d0
Lambdarr(:)=0.d0
Lambdameta(:)=0.d0
LambdaHe(:)=0.d0

Zmetals(1)=1.d0
Zmetals(2)=8.33d-2
!Zmetals(2)=1.d-1
Zmetals(3)=2.04d-9
Zmetals(4)=2.63d-11
Zmetals(5)=6.17d-10
Zmetals(6)=2.45d-4
Zmetals(7)=6.03d-5
Zmetals(8)=4.57d-4
Zmetals(9)=3.02d-8
Zmetals(10)=1.95d-4
Zmetals(11)=2.14d-6
Zmetals(12)=3.39d-5
Zmetals(13)=2.95d-6
Zmetals(14)=3.24d-5
Zmetals(15)=3.20d-7
Zmetals(16)=1.38d-5
Zmetals(17)=1.91d-7
Zmetals(18)=2.51d-6
Zmetals(19)=1.32d-7
Zmetals(20)=2.29d-6
Zmetals(21)=1.48d-9
Zmetals(22)=1.05d-7
Zmetals(23)=1.08d-8
Zmetals(24)=4.68d-7
Zmetals(25)=2.88d-7
Zmetals(26)=2.82d-5
Zmetals(27)=8.32d-8
Zmetals(28)=1.78d-6
Zmetals(29)=1.62d-8
Zmetals(30)=3.98d-8

do i=1,ifile
in_mtl(i)=i+3
enddo


dlogT = (-dlog10(Tmin) + dlog10(Tmax))/dble(ix)
Tmp(0) = Tmin
dlogT1=dlog10(Tmin)
do i=1,ix
dlogT1=dlogT1+dlogT
!Tmp(i)=Tmp(i-1)+10.d0**dlogT
Tmp(i)=10.d0**dlogT1
enddo

!***free-free***
!Draine p.95
ni=1.d0
ne=1.d0
do i=0,ix
Lambdaff(i)=1.422d-25*(1.d0+(0.44d0)/(1.d0+0.058d0*(dlog(Tmp(i)/10.d0**(5.4d0)/Zi**2))))*Zi**2.d0 * (Tmp(i)/1.d4)**0.5d0 * ni * ne
enddo
!***free-free***

!***recomb.-rad.(case B)***
!Draine p.139,319,321
ni=1.d0
ne=1.d0
do i=0,ix
alphaB=2.54d0*1.d-13*Zi*(Tmp(i)/1.d4/Zi/Zi)**(-0.8163d0-0.0208d0*dlog(Tmp(i)/1.d4/Zi/Zi))
!Lambdarr(i)=alphaB * ne * np * (0.684d0-0.0416d0*dlog(Tmp(ix)/1.d4/Zi/Zi))*kb*Tmp(ix)
Lambdarr(i)=alphaB * ne * ni * (0.684d0-0.0416d0*dlog(Tmp(i)/1.d4/Zi/Zi))*kb*Tmp(i)
enddo
!***recomb.-rad.(case B)***
!***metals***
call metal_cool_corona()
!***metals***


j_lin=0
do k=1,ifile
jm=in_mtl(k)
do j=2,jm
j_lin=j_lin+1
do i=0,ix
Tmp1=Tmp(i)
Sigmo=1.d0/(1.d0+dexp(-asigmo*(Tmp(i)-0.5d0*(Tth1+Tth2))/(Tth2-Tth1)))
!Mtl(j,:,k)
call linear(Mtl(k,1,:),Mtl(k,j,:),imtrx,sngl(Tmp1),Lambda_CIE(i,j_lin))
if((k==2).and.(j==jm))then
LambdaHe(i)=Lambda_CIE(i,j_lin)*Zmetals(k)*Sigmo+LambdaHe(i)
endif
if((k>2).and.(j==jm)) then
Lambdameta(i)=Lambda_CIE(i,j_lin)*Zmetals(k)*Sigmo*Zsolar+Lambdameta(i)
endif
enddo
enddo
enddo


open(10,file='cool_rate.dat',access='stream',FORM='UNFORMATTED')
open(11,file='cool_rate_nobi.dat',FORM='FORMATTED')
do i=0,ix
write(10) sngl(Tmp(i)), (sngl(Lambda_CIE(i,j)),j=1,j_lin)
write(11,*) sngl(Tmp(i)), (sngl(Lambda_CIE(i,j)),j=1,j_lin)
enddo
close(10)
close(11)

open(13,file='test_cool_rate.dat',access='stream',FORM='UNFORMATTED')
do i=1,imtrx
write(13) (Mtl(1,j,i),j=1,4)
enddo
close(13)

open(14,file='cool_rate_net.dat',access='stream',FORM='UNFORMATTED')
do i=0,ix
write(14) sngl(Tmp(i)),sngl(Lambdameta(i)),sngl(LambdaHe(i)),sngl(Lambdaff(i)),sngl(Lambdarr(i))
enddo
close(14)

deallocate(Tmp,Rho)
deallocate(Lambda,Lambdaff,Lambdarr,Lambdameta,LambdaHe)
deallocate(Mtl)
end program main


subroutine metal_cool_corona()
!http://wise-obs.tau.ac.il/~orlyg/ion_by_ion/
use comvar
!integer, parameter :: ifile = 30, imtrx=189
integer :: i,j,k,imaxlp
!integer, dimension(1:ifile) :: in_mtl
!H,He,Li,Be,B,C,N,O,F,Ne
!Na,Mg,Al,Si,P,S,Cl,Ar,K,Ca
!Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn
integer :: icrn1=4,icrn2=5,icrn3=6,icrn4=7,icrn5=8,icrn6=9,icrn7=10,icrn8=11,icrn9=12,icrn10=13, &
icrn11=14,icrn12=15,icrn13=16,icrn14=17,icrn15=18,icrn16=19,icrn17=20,icrn18=21,icrn19=22,icrn20=23, &
icrn21=24,icrn22=25,icrn23=26,icrn24=27,icrn25=28,icrn26=29,icrn27=30,icrn28=31,icrn29=32,icrn30=33
!real(4), dimension(:,:), allocatable :: Mtl1,Mtl2,Mtl3,Mtl4,Mtl5,Mtl6,Mtl7,Mtl8,Mtl9
character(2) :: nm



do k=1,ifile
write(nm,'(I2.2)') k
open(14,file=nm//'.dat')
imaxlp=k+3
do j=1,imtrx
!read(14,'E8.2E2') (Mtl(i,j,k),i=1,imaxlp)
read(14,*) (Mtl(k,i,j),i=1,imaxlp)
enddo
close(14)
enddo

end subroutine metal_cool_corona



SUBROUTINE linear(xa,ya,m,x,y)
integer m,i,ms
real(4) :: xa(1:m),ya(1:m)
real(4) :: x
double precision :: y1,y2,t,y

!ms=m
do i=1,m
   if(x-xa(i).le.0.d0) then
      ms=i
      !write(*,*) 'lin'
      go to 12
   endif
enddo
ms=m
12   continue

if(ms.eq.1) ms=2
y1=dble(ya(ms-1))
y2=dble(ya(ms))
!x=xa(ms)
t=(dble(x)-dble(xa(ms-1)))/(dble(xa(ms))-dble(xa(ms-1)))
y=(1.d0-t)*y1+t*y2
END SUBROUTINE linear
