program telegraphexa
implicit none
integer, parameter :: ndx=130
integer :: i,j,loop
double precision, dimension(-1:ndx) :: Phiexa,x,Phiadv
double precision :: f,omega,t=0.d0,cg,Tdiff,tdash,xdash,tmx
double precision :: Lbox,dx,dalph,int=0.d0,alph,elm,bessi1,amp
!character(58) :: dir='/Users/maeda/Desktop/Dropbox/analysis/telegraph-test-1D-1/'
character(21) :: dir='/Users/maeda/Desktop/'

!****************
omega=20.d0
cg=1.d0
Tdiff=2.0d0
Lbox=1.d0
t=1.d0
loop=100
elm=1.d-8
amp=1.d0
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
tdash=t/Tdiff
xdash=x(i)/cg/Tdiff
tmx=tdash-xdash
if(tmx>0.d0) then

dalph=(tdash-xdash)/dble(loop)
int=0.d0
alph=0.d0
do j = 0,loop-1
alph=alph+dalph
int=int+amp*dsin(omega*alph)*dexp(alph*0.5d0)/(dsqrt((tdash-alph)**2.d0-xdash**2.d0+elm))&
*bessi1(0.5d0*dsqrt((tdash-alph)**2.d0-xdash**2.d0+elm))*dalph
enddo

Phiexa(i) = dexp(-tdash*0.5d0)*amp*dsin(omega*tmx) &
+ 0.5d0*xdash*dexp(-tdash*0.5d0)*int

!Phiadv(i) = amp*dsin(omega*(t-x(i)/cg))
Phiadv(i) = amp*dsin(omega*tmx)
elseif(tmx<0.d0) then
Phiexa(i) = 0.d0
Phiadv(i) = 0.d0
endif

enddo

open(22,file = dir//'phiexa.txt')
do i=1,ndx-2
write(22,*) x(i),',',Phiexa(i),',',Phiadv(i)
end do
close(22)
end program telegraphexa


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
