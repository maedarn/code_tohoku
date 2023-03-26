program main
implicit none
double precision a,b,plgndr
integer factorial,c
COMPLEX(16) :: Ylm,Ylmcnj,Ylm_m,i_cmplx

a=2.d0
b=3.d0

i_cmplx=(0.d0,1.d0)
Ylm=cmplx(a,b)
Ylmcnj=conjg(Ylm)
Ylm_m=Ylm/i_cmplx

c=factorial(3)

write(*,*)Ylm,Ylmcnj,Ylm_m,c,plgndr(2,1,1.d0)
end program main


recursive function factorial(nin) result(factorial_n)
    integer, intent(in) :: nin
    integer :: factorial_n

    if (nin > 0) then
        factorial_n = nin*factorial(nin - 1)
        return
   end if

   factorial_n = 1
end function


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
