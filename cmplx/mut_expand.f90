program main
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'
INTEGER :: lsph,msph,lsphmax=1
double precision :: xg_mpi(0:NPE-1),yg_mpi(0:NPE-1),zg_mpi(0:NPE-1)
double precision :: mass1,pi=3.14159265358979d0
double precision :: xg,yg,zg,rg,xr,yr,zr,thetag,phig,plgndr,Ylm_coeff,Ylm1,Ylm0_coeff
COMPLEX(16) :: Ylm,Ylmcnj,Ylmm,Ylmreal,Ylmreal0,Ylmrealb,Ylmreals,i_cmplx,sgn_pm

xg=0.5d0
yg=0.5d0
zg=0.5d0


!multipole moment
i_cmplx=(0.d0,1.d0)
do lsph = 0, lsphmax
  do msph = -lsph, lsph
    Ylm0_coeff=0.d0
    if(msph==0) then; Ylm0_coeff=1.d0; endif
    sgn_pm=0.5d0+dsign(0.5d0,dble(msph)+0.1d0)
    do k=1,ndz-2; do j=1,ndy-2; do i=1,ndx-2
    !spherical polar coordinate
    rg=dsqrt((x(i)-xg)**2.d0+(y(j)-yg)**2.d0+(z(k)-zg)**2.d0)
    thetag=dacos((z(k)-zg)/rg)
    phig=dsign(1.d0,(y(j)-yg))*dacos((x(i)-xg)/dsqrt((x(i)-xg)**2+(y(j)-yg)**2))

    !Real spherical harmonics
    Ylm1=plgndr(lsph,iabs(msph),thetag)*dsqrt((2.d0*dble(lsph)+1.d0)/4.d0/pi*dble(factorial(lsph-iabs(msph)))/dble(factorial(lsph+iabs(msph)))) !Ylm(m=|m|)
    Ylm=cmplx(Ylm1*dcos(dble(iabs(msph))*phig),Ylm1*sin(dble(iabs(msph))*phig)) !Ylm(m=|m|)
    Ylmcnj=conjg(Ylm)
    Ylmrealb=(-1.d0)**msph/dsqrt(2.d0)*(Ylm+Ylmcnj)         !m>0
    Ylmreals=(-1.d0)**msph/dsqrt(2.d0)*(Ylm-Ylmcnj)/i_cmplx !m<0
    Ylmreal0=Ylm                                            !m=0
    Ylmreal=(Ylmrealb*sgn_pm+Ylmrealb*(1.d0-sgn_pm))*(1.d0-Ylm0_coeff)+Ylmreal0*Ylm0_coeff

    Qlm(lsph,msph)=U(i,j,k,1)*dx1*dy1*dz1*rg**dble(lsph)*dble(Ylmreal)*dsqrt(4.d0*pi/(2.d0*dble(lsph)+1.d0))
    enddo; enddo; enddo
  end do
end do

do lsph = 0, lsphmax
  do msph = -lsph, lsph
  Phi(i,j,k)=G*Qlm(lsph,msph)*dsqrt(4.d0*pi/(2.d0*dble(lsph)+1.d0))*Ylm(i,j,k,l,m)/(rg**dble(lsph+1))
  end do
end do

end program main

! numerical in F77
FUNCTION plgndr(l,m,phi_cos)
INTEGER l,m
double precision plgndr,phi_cos
INTEGER i,ll
double precision fact,pll,pmm,pmmp1,somx2

if((m.lt.0).or.(m.gt.l).or.(abs(phi_cos).gt.1.d0)) goto 295 !'bad arguments in plgndr'

pmm=1.d0
 if(m.gt.0) then
 somx2=sqrt((1.d0-phi_cos)*(1.d0+phi_cos))
 fact=1.d0
  do  i=1,m
   pmm=-pmm*fact*somx2 !numerical receipt (6.8.8)
   fact=fact+2.d0
  enddo
 endif
 if(l.eq.m) then
 plgndr=pmm
 else
 pmmp1=phi_cos*dble(2*m+1)*pmm !numerical receipt (6.8.9)
 m+1
  if(l.eq.m+1) then
  plgndr=pmmp1
  else
   do ll=m+2,l
    pll=(phi_cos*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m) !numerical receipt (6.8.7)
    pmm=pmmp1
    pmmp1=pll
   enddo
  plgndr=pll
  endif
 endif

295 continue
return
END FUNCTION plgndr

! calcurate n!
recursive function factorial(n) result(factorial_n)
    implicit none
    integer(int32), intent(in) :: n
    integer(int64) :: factorial_n

    if (n > 0) then
        factorial_n = n*factorial(n - 1)
        return
   end if

   factorial_n = 1
end function

