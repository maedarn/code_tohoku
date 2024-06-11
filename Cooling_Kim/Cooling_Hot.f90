program main
implicit none
integer, parameter :: ix = 100
integer :: i
double precision, parameter :: Tmin=3.d4,Tmax=3.d6
double precision :: dlogT
double precision, dimension(:), allocatable :: Tmp,Rho
double precision, dimension(:), allocatable :: Lambda,Lambdaff,Lambdarr
double precision :: Zi=1.d0,alphaB


allocate(Tmp(0:ix),Rho(0:ix))
allocate(Lambda(0:ix),Lambdaff(0:ix),Lambdarr(0:ix))

dlogT = (-dlog10(Tmin) + dlog10(Tmax))/dble(ix)

Tmp(0) = Tmin
do i=1,ix
Tmp(i)=Tmp(i-1)+10.d0**dlogT
enddo

!***free-free***
!Draine p.95
do i=0,ix
Lambdaff(i)=1.422d-25*(1.d0+(0.44d0)/(1.d0+0.058d0*(dlog(Tmp(ix)/10.d0**(5.4d0)/Zi**2))))*Zi**2.d0 * (Tmp(ix)/1.d4)**0.5d0 * ni * ne
enddo
!***free-free***

!***recomb.-rad.(case B)***
!Draine p.139,319,321
do i=0,ix
alphaB=2.54d0*1.d-13*Zi*(Tmp(ix)/1.d4/Zi/Zi)**(-0.8163d0-0.0208d0*dlog(Tmp(ix)/1.d4/Zi/Zi))
!Lambdarr(i)=alphaB * ne * np * (0.684d0-0.0416d0*dlog(Tmp(ix)/1.d4/Zi/Zi))*kb*Tmp(ix)
Lambdarr(i)=alphaB * ne * ni * (0.684d0-0.0416d0*dlog(Tmp(ix)/1.d4/Zi/Zi))*kb*Tmp(ix)
enddo
!***recomb.-rad.(case B)***

!***metals***
call metal_cool_corona()
!***metals***



deallocate(Tmp,Rho)
deallocate(Lambda,Lambdaff,Lambdarr)
deallocate(Mlt)
end program main



subroutine metal_cool_corona()
!http://wise-obs.tau.ac.il/~orlyg/ion_by_ion/
integer, parameter :: ifile = 30, imtrx=189
integer :: i,j,k,imaxlp
integer, dimension(1:ifile) :: in_mtl
!H,He,Li,Be,B,C,Ni,O,F,Ne
!Na,Mg,Al,Si,P,Sl,Cl,Ar,K,Ca
!Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn
integer :: icrn1=4,icrn2=5,icrn3=6,icrn4=7,icrn5=8,icrn6=9,icrn7=10,icrn8=11,icrn9=12,icrn10=13, &
icrn11=14,icrn12=15,icrn13=16,icrn14=17,icrn15=18,icrn16=19,icrn17=20,icrn18=21,icrn19=22,icrn20=23, &
icrn21=24,icrn22=25,icrn23=26,icrn24=27,icrn25=28,icrn26=29,icrn27=30,icrn28=31,icrn29=32,icrn30=33
!real(4), dimension(:,:), allocatable :: Mtl1,Mtl2,Mtl3,Mtl4,Mtl5,Mtl6,Mtl7,Mtl8,Mtl9
real(4), dimension(:,:,:), allocatable :: Mtl
character(2) :: nm

do i=1,ifile
in_mtl(i)=i+4
enddo

allocate(Mtl(1:ifile,1:imtrx,1:ifile))
Mtl(:,:,:)=0.e0


do k=1,ifile
write(nm,'(I2.2)') k
open(14,file=nm//'.dat')
imaxlp=k+3
do j=1,imtrx
read(14,*) (Mtl(i,j,k),i=1,imaxlp)
enddo
close(14)
enddo

end subroutine metal_cool_corona



SUBROUTINE linear(xa,ya,m,x,y)
integer m,i,ms
double precision :: xa(m),ya(m)
double precision :: y1,y2,t,y

ms=m
do i=1,m
   if(x-xa(i).le.0.d0) then
      ms=i
   endif
enddo
!ms=m

if(ms.eq.1) ms=2
y1=ya(ms-1)
y2=ya(ms)
t=(x-xa(ms-1))/(xa(ms)-xa(ms-1))
y=(1.d0-t)*y1+t*y2
END SUBROUTINE linear



SUBROUTINE bilinear(x1a,x2a,ya,m,n,x1,x2,y)
double precision :: x1a(m),x2a(n),ya(m,n)
do i=1,m
   if(x1-x1a(i).le.0.d0) then
      ms=i
   endif
enddo
ms=m
do i=1,n
   if(x2-x2a(i).le.0.d0) then
      ns=i
   endif
enddo
ns=n
if(ms.eq.1) ms=2
if(ns.eq.1) ns=2
y1=ya(ms-1,ns-1)
y2=ya(ms,ns-1)
y3=ya(ms,ns)
y4=ya(ms-1,ns)
t=(x1-x1a(ms-1))/(x1a(ms)-x1a(ms-1))
u=(x2-x2a(ns-1))/(x2a(ns)-x2a(ns-1))
y=(1.d0-t)*(1.d0-u)*y1+t*(1.d0-u)*y2+t*u*y3+(1.d0-t)*u*y4
END SUBROUTINE bilinear
