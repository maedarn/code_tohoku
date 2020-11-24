program telegraphexa
implicit none
integer, parameter :: ndx=130
integer :: i,j,loop
double precision, dimension(-1:ndx) :: Phiexa,x,Phiadv
double precision :: f,omega,t=0.d0,cg,Tdiff,xpc,xmc,tmx,bessi0,bessi1,kappa,alxp,alxm
double precision :: Lbox,dx,dalph,int0=0.d0,int1=0.d0,alph,elm,amp,inf0,inf1,k
!character(58) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-1D-1/'
character(21) :: dir='/Users/maeda/Desktop/'

!****************
omega=20.d0
cg=1.d0
Tdiff=2.0d0
Lbox=1.d0
t=0.5d0
loop=100
elm=1.d-8
amp=1.d0
k=3.1415926535d0*2.d0/Lbox * 1.d0
kappa=Tdiff*0.5d0
!****************

dx=Lbox/dble(ndx-2)
x(0)=-dx/2.d0
x(1)=dx/2.d0
x(-1)=-dx/2.d0-dx
do i=2,ndx
x(i)=x(i-1)+dx
enddo

!f=dsin(omega*t)
Phiexa(:)=0.d0
Phiadv(:)=0.d0

do i=1,ndx-2
xpc=x(i)+cg*t
xmc=x(i)-cg*t

dalph=(xpc-xmc)/dble(loop)
int0=0.d0
alph=xmc
do j = 0,loop-1
alph=alph+dalph
alxp=cg*t+x(i)-alph
alxm=cg*t-x(i)+alph
int0=int0+1.d0/cg * (inf1(alph,k,amp)+kappa*inf0(alph,k,amp)) * bessi0(2.d0*(dsqrt(cg*alxp*alxm+elm)))*dalph   !alph=x
enddo

int1=0.d0
alph=xmc
do j = 0,loop-1
alph=alph+dalph
alxp=cg*t+x(i)-alph
alxm=cg*t-x(i)+alph
int1=int1+dsqrt(kappa*kappa/cg/cg) * cg * t * inf0(alph,k,amp) * bessi1(2.d0*(dsqrt(cg*alxp*alxm+elm)))/(dsqrt(alxp*alxm+elm))*dalph   !alph=x
enddo


Phiexa(i) = dexp(-t*kappa)*(0.5d0*(inf0(xpc,k,amp)+inf0(xmc,k,amp)) + int0 + int1)

!Phiadv(i) = amp*dsin(omega*(t-x(i)/cg))
Phiadv(i) = amp*dsin(-cg*k*x(i)+t)
enddo

open(22,file = dir//'phiexa.txt')
do i=1,ndx-2
write(22,*) x(i),',',Phiexa(i),',',Phiadv(i)
end do
close(22)
end program telegraphexa

FUNCTION bessi0(x)
DOUBLE PRECISION bessi0,x
DOUBLE PRECISION ax
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y

p1=1.0d0
p2=3.5156229d0
p3=3.0899424d0
p4=1.2067492d0
p5=0.2659732d0
p6=0.360768d-1
p7=0.45813d-2

q1=0.39894228d0
q2=0.1328592d-1
q3=0.225319d-2
q4=-0.157565d-2
q5=0.916281d-2
q6=-0.2057706d-1
q7=0.2635537d-1
q8=-0.1647633d-1
q9=0.392377d-2

if (abs(x).lt.3.75) then
y=(x/3.75d0)**2d0
bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
else
ax=dabs(x)
y=3.75/ax
bessi0=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
endif

END FUNCTION bessi0


FUNCTION bessi1(x)
double precision bessi1,x
double precision ax
DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
!SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
!DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
!* 0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
!DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
!* -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
!* -0.2895312d-1,0.1787654d-1,-0.420059d-2/
p1=0.5d0
p2=0.87890594d0
p3=0.51498869d0
p4=0.15084934d0
p5=0.2658733d-1
p6=0.301532d-2
p7=0.32411d-3

q1=0.39894228d0
q2=-0.3988024d-1
q3=-0.362018d-2
q4=0.163801d-2
q5=-0.1031555d-1
q6=0.2282967d-1
q7=-0.2895312d-1
q8=0.1787654d-1
q9=-0.420059d-2


if (dabs(x) .lt. 3.75d0) then
y=(x/3.75d0)**2d0
bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
else
ax=dabs(x)
y=3.75d0/ax
bessi1=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))

if(x.lt.0.d0) then
bessi1=-bessi1
endif

endif

END FUNCTION bessi1


FUNCTION inf0(x,k,amp)
double precision inf0,x
double precision k,amp

inf0=amp*dsin(k*x)
END FUNCTION inf0


FUNCTION inf1(x,k,amp)
double precision inf1,x
double precision k,amp

inf1=amp*dcos(k*x)
END FUNCTION inf1
