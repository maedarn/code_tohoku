program main
implicit none
double precision :: e, a1, a3, R, z, rho, Gpi, dz, dR, G=1.11142d-4, Phimean, vesc, wvesc, wPhimean, rtot, theta
double precision, allocatable, dimension(:,:) :: Phi
integer :: i, j, nr, nz
!----parameter----
a1=2.0d0
!a1=9.9999d0
a3=15.0d0
e=sqrt(1.0d0-(a1/a3)**2.0d0)
nr=10
nz=100
rho=1.0d4
!rho=3.d0*1.0d4/4.d0/4.d0/30.d0
Gpi=12.56637d0*G/4.d0
!----parameter----

allocate(Phi(1:nr,1:nz))
dz=a3/dble(nz)
Phi(:,:)=0.d0

do j=1,nz
 z=dz*dble(j)
 !R=sqrt(dabs(a1**2.d0 - z**2.d0/(1.d0-e**2.d0)))
 R=sqrt(dabs(a1*a1 - z*z*a1*a1/a3/a3))
 !do i=1,nr
 i=1
 !Phi(i,j) = -Gpi*rho*(dlog((1.d0+e)/(1.d0-e))*a1*a1/e - &
 !            (1.d0-e*e)*((1.d0/(1.d0-e*e))-dlog((1.d0+e)/(1.d0-e))/(2.d0*e))*R*R/e/e - &
 !            2.d0*(1.d0-e*e)*(dlog((1.d0+e)/(1.d0-e))/(2.d0*e)-1.d0)*z*z/e/e)
 Phi(i,j) = -Gpi*rho*((dlog(1.d0+e)-dlog(1.d0-e))*a1*a1/e &
            -(1.d0-e*e)/e/e * ((1.d0/(1.d0-e*e))-(dlog(1.d0+e)-dlog(1.d0-e))/2.d0/e)*R*R &
            -2.d0*(1.d0-e*e)/e/e * ((dlog(1.d0+e)-dlog(1.d0-e))/2.d0/e - 1.d0 )*z*z)
 !enddo
enddo

Phimean=0.d0
wPhimean=0.d0
i=1
do j=1,nz
z=dz*dble(j)
R=sqrt(dabs(a1*a1 - z*z*a1*a1/a3/a3))
Phimean=Phimean+Phi(i,j)
wPhimean=wPhimean+Phi(i,j)*R
rtot=R+rtot
enddo
Phimean=Phimean/dble(nz)
wPhimean=wPhimean/rtot


open(3,file='phifil.DAT')

do j=1,nz
i=1
z=dz*dble(j)
!R=sqrt(dabs(a1**2.d0 - z**2.d0/(1.d0-e**2.d0)))
R=sqrt(dabs(a1**2.d0 - z**2.d0/a3**2.d0 * a1**2.d0))
vesc=sqrt(2.d0*dabs(Phimean))
wvesc=sqrt(2.d0*dabs(wPhimean))
theta=datan(R/z)/3.14159265d0*2.d0*90.d0
write(3,*) R,z,Phi(i,j),-Gpi*rho*a1**2.d0*4.d0/3.d0,-Gpi*rho*a3**2.d0*4.d0/3.d0,e,Phimean,vesc,wPhimean,wvesc,theta !,sqrt(R*R+z*z)
enddo

close(3)

deallocate(Phi)
end program main
