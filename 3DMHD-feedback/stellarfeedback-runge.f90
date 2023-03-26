SUBROUTINE feedback(dt,mode)
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
DOUBLE PRECISION  :: dt,dxi,intt=0.d0,raninp,ranout,Rst,rloop
DOUBLE PRECISION  :: Q,lQ
DOUBLE PRECISION  :: mstarmn=1.d0,pstar,Mstar1,Mstarlp,T0   = 1.0d3
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt,jtime=0!,nid=0
!INTEGER, parameter :: nstar=1000000
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU
character(3) Nfinal,itime
integer :: i,j,k,itrrst=10,rixmn,rjymn,rkzmn,rixmx,rjymx,rkzmx,nrixmn,nrjymn,nrkzmn,nrixmx,nrjymx,nrkzmx
integer :: rsix,rsjy,rskz,nstloop,NRANKdm,nstdm
DOUBLE PRECISION  :: rsixc,rsjyc,rskzc
double precision :: div(-1:ndx,-1:ndy,-1:ndz,1:4)
double precision :: SFE,SFratio=0.5d0,Mmassivetot=0.d0,rdumy
DOUBLE PRECISION ran
!double precision :: Ustar(1:10,1:nstar)

!write(*,*) NRANK,U(1,1,1,1),U(1,1,1,5),'pint-feed1'

if(mode==1) then

Ustar(9,0)=0.d0
nidnw=0
intt=intt+dt
N_MPI(20)=4; N_MPI(1)=1; N_MPI(2)=2; N_MPI(3)=3; N_MPI(4)=4; iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
div(i,j,k,1)=(U(i+1,j,k,2)-U(i-1,j,k,2))*0.5d0/dx1+(U(i,j+1,k,3)-U(i,j-1,k,3))*0.5d0/dy1+(U(i,j,k+1,4)-U(i,j,k-1,4))*0.5d0/dz1
div(i,j,k,2)=(U(i+1,j,k,2)-U(i-1,j,k,2))*0.5d0/dx1
div(i,j,k,3)=(U(i,j+1,k,3)-U(i,j-1,k,3))*0.5d0/dy1
div(i,j,k,4)=(U(i,j,k+1,4)-U(i,j,k-1,4))*0.5d0/dz1
end do; end do; end do

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
if((div(i,j,k,1)<0.d0).and.(div(i,j,k,2)<0.d0).and.(div(i,j,k,3)<0.d0).and.(div(i,j,k,4)<0.d0).and.(U(i,j,k,1)>rhoth)) then
!if(U(i,j,k,1)>rhoth) then
SFE=LSFE*dsqrt(U(i,j,k,1)/(4.04d0*1.d3))*dt
!SFE=LSFE*U(i,j,k,1)/dsqrt(G*U(i,j,k,1))*dx1*dy1*dz1*dt
!pstar=(1-dexp(-SFE/mstarmn))
!call ran0(ran,1)
!if(ran<pstar) then
Mstar1=U(i,j,k,1)*SFE*Msun
Mstarlp=0.d0
Mmassivetot=0.d0
do nstdm=1,nlpmass
if(Mstarlp>Mstar1) then
goto 2021
endif
idum=1+NRANK
call ran0(raninp,idum)
call IMF(raninp,ranout)
Mstarlp=Mstarlp+ranout

!if(ran<pstar) then
if(Mmassive<ranout) then
Mmassivetot=Mmassivetot+ranout
!nid=nid+1
nidnw=nidnw+1
Ustar(1,nid+nidnw)=dble(IST*Ncellx + i)*dx1+dx1*0.5d0
Ustar(2,nid+nidnw)=dble(JST*Ncelly + j)*dy1+dy1*0.5d0
Ustar(3,nid+nidnw)=dble(KST*Ncellz + k)*dz1+dz1*0.5d0
Ustar(4,nid+nidnw)=U(i,j,k,2)
Ustar(5,nid+nidnw)=U(i,j,k,3)
Ustar(6,nid+nidnw)=U(i,j,k,4)
Ustar(7,nid+nidnw)=ranout !U(i,j,k,1)*dx1*dy1*dz1*SFratio
Ustar(8,nid+nidnw)=intt
Ustar(9,nid+nidnw)=dble(nid)
lQ=aq+bq*dlog10(ranout)+cq*dlog10(ranout)**2+dq*dlog10(ranout)**3+eq*dlog10(ranout)**4+fq*dlog10(ranout)**5
!Q=10.d0**lQ
!Ustar(10,nid+nidnw)=10.d0**lQ
Ustar(10,nid+nidnw)=lQ

Ustar(11,nid+nidnw)=dble(IST*Ncellx + i)
Ustar(12,nid+nidnw)=dble(JST*Ncelly + j)
Ustar(13,nid+nidnw)=dble(KST*Ncellz + k)

Rst=10.d0*(10**(Ustar(10,nid+nidnw)-49.d0))**(1.d0/3.d0)*(U(i,j,k,1)/1.27d0/10.d0)**(-2.d0/3.d0) !loop
Ustar(14,nid+nidnw)=Rst

Ustar(14,0)=1.d0

Ustar(9,0)=dble(nidnw)

write(*,*) NRANK,ranout,Ustar(10,nid+nidnw),Ustar(7,nid+nidnw),Ustar(1,nid+nidnw),Ustar(2,nid+nidnw),Ustar(3,nid+nidnw),nid,nidnw,'pint-feed2'

endif
enddo
2021 continue
U(i,j,k,1)=U(i,j,k,1)-Mmassivetot/Msun
endif
end do; end do; end do
call collectST()
endif

if(mode==2) then

do nstloop=1,nid
!i=mod(int(Ustar(1,nid)),dx1)
!j=mod(int(Ustar(2,nid)),dy1)
!k=mod(int(Ustar(3,nid)),dz1)

!i=int(dmod(dmod(Ustar(1,nstloop),dble(Ncellx)*dx1),dx1))
!j=int(dmod(dmod(Ustar(2,nstloop),dble(Ncelly)*dy1),dy1))
!k=int(dmod(dmod(Ustar(3,nstloop),dble(Ncellz)*dz1),dz1))
 !k=mod(int(Ustar(3,nid)),Ncellz)
!Rst=10.d0*(10**(Ustar(10,nstloop)-49.d0))**(1.d0/3.d0)*(U(i,j,k,1)/1.27d0/10.d0)**(-2.d0/3.d0) !loop

!rixmn=Ustar(1,nid)-Rst
!rjymn=Ustar(2,nid)-Rst
!rkzmn=Ustar(3,nid)-Rst
!rixmx=Ustar(1,nid)+Rst
!rjymx=Ustar(2,nid)+Rst
!rkzmx=Ustar(3,nid)+Rst
!write(*,*) NRANK,Ustar(10,nstloop),Rst,dx1,nid,'pint-feed3'

!nrixmn=int(dmod(Ustar(1,nstloop)-Rst,dble(Ncellx)*dx1))
!nrjymn=int(dmod(Ustar(2,nstloop)-Rst,dble(Ncelly)*dy1))
!nrkzmn=int(dmod(Ustar(3,nstloop)-Rst,dble(Ncellz)*dz1))
!nrixmx=int(dmod(Ustar(1,nstloop)+Rst,dble(Ncellx)*dx1))
!nrjymx=int(dmod(Ustar(2,nstloop)+Rst,dble(Ncelly)*dy1))
!nrkzmx=int(dmod(Ustar(3,nstloop)+Rst,dble(Ncellz)*dz1))

!rixmn=int((Ustar(1,nstloop)-Rst)/(dble(Ncellx)*dx1))
!rjymn=int((Ustar(2,nstloop)-Rst)/(dble(Ncelly)*dy1))
!rkzmn=int((Ustar(3,nstloop)-Rst)/(dble(Ncellz)*dz1))
!rixmx=int((Ustar(1,nstloop)+Rst)/(dble(Ncellx)*dx1))
!rjymx=int((Ustar(2,nstloop)+Rst)/(dble(Ncelly)*dy1))
!rkzmx=int((Ustar(3,nstloop)+Rst)/(dble(Ncellz)*dz1))

!nrixmn=int((Ustar(1,nstloop)-Rst)/(dble(Ncellx)*dx1)))
!nrjymn=int((Ustar(2,nstloop)-Rst)/(dble(Ncelly)*dy1)))
!nrkzmn=int((Ustar(3,nstloop)-Rst)/(dble(Ncellz)*dz1)))
!nrixmx=int((Ustar(1,nstloop)+Rst)/(dble(Ncellx)*dx1)))
!nrjymx=int((Ustar(2,nstloop)+Rst)/(dble(Ncelly)*dy1)))
!nrkzmx=int((Ustar(3,nstloop)+Rst)/(dble(Ncellz)*dz1)))

!rixmn=int((Ustar(1,nstloop)-Rst)/(dble(Ncellx)*dx1))
!rjymn=int((Ustar(2,nstloop)-Rst)/(dble(Ncelly)*dy1))
!rkzmn=int((Ustar(3,nstloop)-Rst)/(dble(Ncellz)*dz1))
!rixmx=int((Ustar(1,nstloop)+Rst)/(dble(Ncellx)*dx1))
!rjymx=int((Ustar(2,nstloop)+Rst)/(dble(Ncelly)*dy1))
!rkzmx=int((Ustar(3,nstloop)+Rst)/(dble(Ncellz)*dz1))

rixcn=int((Ustar(1,nstloop))/(dble(Ncellx)*dx1))
rjycn=int((Ustar(2,nstloop))/(dble(Ncelly)*dy1))
rkzcn=int((Ustar(3,nstloop))/(dble(Ncellz)*dz1))

write(*,*) rixcn,rjycn,rkxcn,'CN1'

if(rixcn.eq.-1 ) rixcn = NSPLTx-1
if(rjycn.eq.-1 ) rjycn = NSPLTy-1
if(rkzcn.eq.-1 ) rkzcn = NSPLTz-1
if(rixcn.eq.NSPLTx) rixcn = 0
if(rjycn.eq.NSPLTy) rjycn = 0
if(rkzcn.eq.NSPLTz) rkzcn = 0

!rixcnm1=rixcn-1
!rjycnm1=rixcn-1
!rkzcnm1=rixcn-1

!rixcnp1=rixcn+1
!rjycnp1=rixcn+1
!rkzcnp1=rixcn+1

!rixmn=int((Ustar(1,nstloop)-Rst)/(dble(Ncellx)*dx1))
!rjymn=int((Ustar(2,nstloop)-Rst)/(dble(Ncelly)*dy1))
!rkzmn=int((Ustar(3,nstloop)-Rst)/(dble(Ncellz)*dz1))
!rixmx=int((Ustar(1,nstloop)+Rst)/(dble(Ncellx)*dx1))
!rjymx=int((Ustar(2,nstloop)+Rst)/(dble(Ncelly)*dy1))
!rkzmx=int((Ustar(3,nstloop)+Rst)/(dble(Ncellz)*dz1))

!if(rixmn.eq.-1 ) rixmn = NSPLTx-1
!if(rjymn.eq.-1 ) rjymn = NSPLTy-1
!if(rkzmn.eq.-1 ) rkzmn = NSPLTz-1
!if(rixmn.eq.NSPLTx) rixmn = 0
!if(rjymn.eq.NSPLTy) rjymn = 0
!if(rkzmn.eq.NSPLTz) rkzmn = 0

!if(rixmx.eq.-1 ) rixmx = NSPLTx-1
!if(rjymx.eq.-1 ) rjymx = NSPLTy-1
!if(rkzmx.eq.-1 ) rkzmx = NSPLTz-1
!if(rixmx.eq.NSPLTx) rixmx = 0
!if(rjymx.eq.NSPLTy) rjymx = 0
!if(rkzmx.eq.NSPLTz) rkzmx = 0


!if(rixmn.eq.-1 ) rixmn = NSPLTx-1
!if(rjymn.eq.-1 ) rjymn = NSPLTy-1
!if(rkzmn.eq.-1 ) rkzmn = NSPLTz-1
!if(rixmx.eq.NSPLTx) rixmx = 0
!if(rjymx.eq.NSPLTy) rjymx = 0
!if(rkzmx.eq.NSPLTz) rkzmx = 0


!if((IST==rixmn).or.(IST==rixmx).or.(JST==rjymn).or.(JST==rjymx).or.(KST==rkzmn).or.(KST==rkzmx)) then

!if((IST==rixcn).or.(IST==rixcnm1).or.(IST==rixcnp1).or.(JST==rjycn).or.(JST==rjycnm1)..or.(JST==rjycnp1)&
!.or.(KST==rkzcn).or.(KST==rkzcnm1).or.(KST==rkzcnp1)) then

if((IST==rixcn).and.(JST==rjycn).and.(KST==rkzcn)) then

i=int(dmod(dmod(Ustar(1,nstloop),dble(Ncellx)*dx1),dx1))
j=int(dmod(dmod(Ustar(2,nstloop),dble(Ncelly)*dy1),dy1))
k=int(dmod(dmod(Ustar(3,nstloop),dble(Ncellz)*dz1),dz1))

write(*,*) NRANK,IST,JST,KST,i,j,k,'CN2'

Rst=10.d0*(10**(Ustar(10,nstloop)-49.d0))**(1.d0/3.d0)*(U(i,j,k,1)/1.27d0/10.d0)**(-2.d0/3.d0) !loop
Ustar(14,nstloop)=Rst
write(*,*) NRANK,Ustar(10,nstloop),Rst,dx1,nid,'pint-feed3'
!call BC_ST_RS()
end if

rixmn=int((Ustar(1,nstloop)-Ustar(14,nstloop))/(dble(Ncellx)*dx1))
rjymn=int((Ustar(2,nstloop)-Ustar(14,nstloop))/(dble(Ncelly)*dy1))
rkzmn=int((Ustar(3,nstloop)-Ustar(14,nstloop))/(dble(Ncellz)*dz1))
rixmx=int((Ustar(1,nstloop)+Ustar(14,nstloop))/(dble(Ncellx)*dx1))
rjymx=int((Ustar(2,nstloop)+Ustar(14,nstloop))/(dble(Ncelly)*dy1))
rkzmx=int((Ustar(3,nstloop)+Ustar(14,nstloop))/(dble(Ncellz)*dz1))

if(rixmn.eq.-1 ) rixmn = NSPLTx-1
if(rjymn.eq.-1 ) rjymn = NSPLTy-1
if(rkzmn.eq.-1 ) rkzmn = NSPLTz-1
if(rixmn.eq.NSPLTx) rixmn = 0
if(rjymn.eq.NSPLTy) rjymn = 0
if(rkzmn.eq.NSPLTz) rkzmn = 0

if(rixmx.eq.-1 ) rixmx = NSPLTx-1
if(rjymx.eq.-1 ) rjymx = NSPLTy-1
if(rkzmx.eq.-1 ) rkzmx = NSPLTz-1
if(rixmx.eq.NSPLTx) rixmx = 0
if(rjymx.eq.NSPLTy) rjymx = 0
if(rkzmx.eq.NSPLTz) rkzmx = 0

!rixcnm1=rixcn-1
!rjycnm1=rixcn-1
!rkzcnm1=rixcn-1

!rixcnp1=rixcn+1
!rjycnp1=rixcn+1
!rkzcnp1=rixcn+1

!if(rixcnm1.eq.-1 ) rixcnm1 = NSPLTx-1
!if(rjycnm1.eq.-1 ) rjycnm1 = NSPLTy-1
!if(rkzcnm1.eq.-1 ) rkzcnm1 = NSPLTz-1
!if(rixcnm1.eq.NSPLTx) rixcnm1 = 0
!if(rjycnm1.eq.NSPLTy) rjycnm1 = 0
!if(rkzcnm1.eq.NSPLTz) rkzcnm1 = 0

!if(rixcnp1.eq.-1 ) rixcnp1 = NSPLTx-1
!if(rjycnp1.eq.-1 ) rjycnp1 = NSPLTy-1
!if(rkzcnp1.eq.-1 ) rkzcnp1 = NSPLTz-1
!if(rixcnp1.eq.NSPLTx) rixcnp1 = 0
!if(rjycnp1.eq.NSPLTy) rjycnp1 = 0
!if(rkzcnp1.eq.NSPLTz) rkzcnp1 = 0

!if(i<1) rixcnm1 = NSPLTx-1
!if(j.eq.-1 ) rjycnm1 = NSPLTy-1
!if(k.eq.-1 ) rkzcnm1 = NSPLTz-1
!if(i.eq.NSPLTx) rixcnm1 = 0
!if(j.eq.NSPLTy) rjycnm1 = 0
!if(k.eq.NSPLTz) rkzcnm1 = 0

!k=mod(int(Ustar(3,nid)),Ncellz)
!if((IST==rixcn).or.(IST==rixmn).or.(IST==rixmx).or.(JST==rjycn).or.(JST==rjymn).or.(JST==rjymx)&
!.or.(KST==rkzcn).or.(KST==rkzmn).or.(KST==rkzmx)) then

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
rdumy2=y(j)-Ustar(2,nstloop)
rdumy3=z(k)-Ustar(3,nstloop)
!if(rdumy2>Lbox*0.5d0) rdumy2=rdumy2-Lbox
!if(rdumy3>Lbox*0.5d0) rdumy3=rdumy3-Lbox
rdumy2=rdumy2-Lbox*dble(int(rdumy2/(Lbox*0.5d0)))
rdumy3=rdumy3-Lbox*dble(int(rdumy3/(Lbox*0.5d0)))

rdumy=dsqrt((x(i)-Ustar(1,nstloop))**2+(rdumy2)**2+(rdumy3)**2)

if(Ustar(14,nstloop)>rdumy) then
write(*,*)x(i),y(j),z(k),NRANK,'ion'
!ni=int(dmod(dabs(x(i)-Ustar(1,nstloop)),dble(Ncellx)*dx1))!*dsign(x(i)-Ustar(1,nstloop))
!nj=int(dmod(dabs(y(j)-Ustar(2,nstloop)),dble(Ncelly)*dy1))!*dsign(y(j)-Ustar(2,nstloop))
!nk=int(dmod(dabs(z(k)-Ustar(3,nstloop)),dble(Ncellz)*dz1))!*dsign(z(k)-Ustar(3,nstloop))
!if(ni.eq.-1 ) ni = NSPLTx-1
!if(nj.eq.-1 ) nj = NSPLTy-1
!if(nk.eq.-1 ) nk = NSPLTz-1
!if(ni.eq.NSPLTx) ni = 0
!if(nj.eq.NSPLTy) nj = 0
!if(nk.eq.NSPLTz) nk = 0
!if((IST==ni).and.(JST==nj).and.(KST==nk)) then
!T(j,j,k,)=1.d4/T0
U(j,j,k,5)=(1.d4/T0)*kb*U(j,j,k,1)/1.27d0
endif
!endif
end do; end do; end do
!endif

enddo

endif


if(mode==3) then
do nstloop=1,nid
rsix=int((Ustar(1,nstloop))/(dble(Ncellx)*dx1))
rsjy=int((Ustar(2,nstloop))/(dble(Ncelly)*dy1))
rskz=int((Ustar(3,nstloop))/(dble(Ncellz)*dz1))

rsixc=dmod(dmod(Ustar(1,nstloop),dble(Ncellx)*dx1),dx1)
rsjyc=dmod(dmod(Ustar(2,nstloop),dble(Ncelly)*dy1),dy1)
rskzc=dmod(dmod(Ustar(3,nstloop),dble(Ncellz)*dz1),dz1)

if((IST==rsix).and.(JST==rsjy).and.(KST==rskz)) then
!NRANKdm = rsix + rsjy*NSPLTx + rskz*NSPLTx*NSPLTy
call runge(Ustar(1,nstloop),Ustar(4,nstloop),dt,nstloop,rsixc,rsjyc,rskzc,1)
call runge(Ustar(2,nstloop),Ustar(5,nstloop),dt,nstloop,rsixc,rsjyc,rskzc,2)
call runge(Ustar(3,nstloop),Ustar(6,nstloop),dt,nstloop,rsixc,rsjyc,rskzc,3)
!Ustar(1,nstloop)=Ustar(1,nstloop)+Ustar(4,nstloop)*dt/dx
!Ustar(2,nstloop)=Ustar(2,nstloop)+Ustar(5,nstloop)*dt/dy
!Ustar(3,nstloop)=Ustar(3,nstloop)+Ustar(6,nstloop)*dt/dz
endif
enddo

call BC_ST()
endif

end SUBROUTINE feedback

!SUBROUTINE runge(xint,vint,tint,dt,nid)
SUBROUTINE runge(xint,vint,dt,nid,rsixc,rsjyc,rskzc,mode)
  !implicit none
  double precision :: tstart,xstart,tend,vstart
  double precision :: tnx,xnx,vnx,dt,rsixc,rsjyc,rskzc
  double precision :: tstp0,tstp1,tstp2,tstp3
  double precision :: kx0,kx1,kx2,kx3
  double precision :: kv0,kv1,kv2,kv3
  double precision :: f0,f1,f2,f3
  double precision :: v0,v1,v2,v3
  double precision :: funcf,funcv
  !double precision,allocatable,dimension(:) :: tint,xint,vint,xexa,vexa
  integer :: i,iend,nid,mode
  character(21) :: fname='/Users/maeda/Desktop/'

  !tstart=0.d0
  !tend  =100.d0
  !xstart=0.d0
  !vstart=1.d0
  !tstart=0.d0
  !dt=0.1d0
  !nid=100000
  
  !allocate(tint(1:nid),xint(1:nid),vint(1:nid),xexa(1:nid),vexa(1:nid))
  !allocate(xint(1:nid),vint(1:nid))
  !tint(:)=0.d0
  !xint(:)=0.d0
  !vint(:)=0.d0
 
 
  !t=tstart
  !x=xstart
  !v=vstart
  !do i=1,nid
    !write(10,'(2f13.5)') t,x
 
    !tstp0=tint(i)
    kv0=vint
    kx0=xint
    f0=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx0,kv0)
    v0=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx0,kv0)
 
    !tstp1=tint(i)+dt/2.0
    kv1=vint+f0*dt/2.0
    kx1=xint+v0*dt/2.0
    f1=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx1,kv1)
    v1=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx1,kv1)
 
    !tstp2=tint(i)+dt/2.0
    kv2=vint+f1*dt/2.0
    kx2=xint+v1*dt/2.0
    f2=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx2,kv2)
    v2=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx2,kv2)
 
    !tstp3=tint(i)+dt
    kv3=vint+f2*dt
    kx3=xint+v2*dt
    f3=funcf(int(rsixc),int(rsjyc),int(rskzc),mode,kx3,kv3)
    v3=funcv(int(rsixc),int(rsjyc),int(rskzc),mode,kx3,kv3)
 
    vnx=vint+(f0+f1*2.0+f2*2.0+f3)*dt/6.0
    xnx=xint+(v0+v1*2.0+v2*2.0+v3)*dt/6.0

    !tnx=tint(i)+dt

    !tint(i)=tnx
    xint=xnx
    vint=vnx
    !xexa(i)=vstart*dsin(t)
    !vexa(i)=-vstart*dcos(t)
 
  !end do

  !open(10,file=fname//'runge4.dat')
  !do i=1,nid
  !write(10,*) tint(i),',',xint(i),',',vint(i),',',xexa(i),',',vexa(i)
  !end do
  !close(10)

  !deallocate(xint,vint)
  !deallocate(tint,xint,vint,xexa,vexa)
 
  !close(10)
 
end SUBROUTINE runge



 
function funcf(i,j,k,mode,kx,kv)
  USE comvar
  USE slfgrv
  !double precision :: t,x,v
  double precision :: kx,kv
  double precision :: funcf
  integer :: i,j,k,mode
 
  !funcf=(1.0-x**2)
  !funcf=-x!9.8d0!(1.0-x**2)

  !U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
  !U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
  !U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0

  !U(i,j,k,2) = -  ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
  !U(i,j,k,3) = -  ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dyi *0.5d0
  !U(i,j,k,4) = -  ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dzi *0.5d0
  if (mode==1)then
  funcf = -  ( Phi(i+1,j,k)-2.d0*Phi(i,j,k)+Phi(i-1,j,k) ) / dx1 / dx1
  endif
  if (mode==2)then
  funcf = -  ( Phi(i,j+1,k)-2.d0*Phi(i,j,k)+Phi(i,j-1,k) ) / dy1 / dy1
  endif
  if (mode==3)then
  funcf = -  ( Phi(i,j,k+1)-2.d0*Phi(i,j,k)+Phi(i,j,k-1) ) / dz1 / dz1
  endif
 
  return
end function funcf





function funcv(i,j,k,mode,kx,kv)
  !double precision :: t,x,v
  double precision :: funcv
  integer :: i,j,k,mode
  double precision :: kx,kv
 

  funcv = kv
  return
end function funcv



SUBROUTINE IMF(raninp,ranout)
integer :: nmesh,cnc1,cnc2,nran!,idum
integer :: i
double precision :: raninp,ranout
double precision :: mmax,mmin,nsisu,scale,slide,dmy,Ncumcnc1,Ncumcnc2
double precision :: M1,M2,coefflow,coeffmid,coeffhig,pwllow,pwlmid,pwlhig,coefflow2,coeffmid2,coeffhig2
double precision, allocatable, dimension(:) :: M,lgM,dN,N,lgN,Ncum,Mran,ranini,ts1,ts2,ts3,ts4
double precision :: dm,lgmmax,lgmmin,ddN,Ntot
character(21) :: dir='/Users/maeda/Desktop/'

!-------paramerter------
!nmesh = 1000
!nran=10000000
mmin = 0.01d0
mmax = 150.d0
Ntot = 100.d0
!slide=5.d0
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

!allocate(M(0:nmesh),lgM(0:nmesh),dN(0:nmesh),N(0:nmesh),lgN(0:nmesh),Ncum(0:nmesh),Mran(1:nran),ranini(1:nran))
!allocate(ts1(1:nran),ts2(1:nran),ts3(1:nran),ts4(1:nran))

!lgmmax = dlog10(mmax)
!lgmmin = dlog10(mmin)
!dm = (lgmmax-lgmmin)/dble(nmesh)
!cnc1=int((dlog10(M1)-dlog10(mmin))/dm)
!cnc2=int((dlog10(M2)-dlog10(mmin))/dm)



Ntot=coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0)-(mmin)**(pwllow+1.d0)/(pwllow+1.d0))+&
     coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
     coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncumcnc2=coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2     )**(pwlhig+1.d0)/(pwlhig+1.d0))
Ncumcnc1=coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0)-(M1     )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
         coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2     )**(pwlhig+1.d0)/(pwlhig+1.d0))


Ncumcnc2=Ncumcnc2/Ntot
Ncumcnc1=Ncumcnc1/Ntot



if (raninp>Ncumcnc1) then
ranout=((coefflow*((M1  )**(pwllow+1.d0)/(pwllow+1.d0))+&
          coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))-coeffmid*((M1  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-coeffhig*((M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*raninp)/coefflow*(pwllow+1.d0))**(1.d0/(pwllow+1.d0))

else if ((Ncumcnc1>raninp).and.(raninp>Ncumcnc2)) then
ranout=((coeffmid*((M2  )**(pwlmid+1.d0)/(pwlmid+1.d0))+&
          coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0)-(M2  )**(pwlhig+1.d0)/(pwlhig+1.d0)) &
          -Ntot*raninp)*(pwlmid+1.d0)/coeffmid)**(1.d0/(pwlmid+1.d0))

else if (Ncumcnc2>raninp) then
ranout=((coeffhig*((mmax)**(pwlhig+1.d0)/(pwlhig+1.d0))-Ntot*raninp)*(pwlhig+1.d0)/coeffhig)**(1.d0/(pwlhig+1.d0))

endif


!deallocate(M,lgM,dN,N,lgN,Ncum,Mran,ranini)
!deallocate(ts1,ts2,ts3,ts4)
end SUBROUTINE IMF

!do i = 0, Ncell+1
!   N(i) = U(i,1)/1.27d0
!  Tn(i) = U(i,5)/( kb*(ndp(i)+ndH(i)+ndH2(i)+ndHe(i)+ndHep(i)) )
!  Pn(i) = U(i,5)
!end do
!do i = 0, Ncell
!  rtTx  = dsqrt(  ( dx(i)*Tn(i+1)+dx(i+1)*Tn(i) ) /( dx(i) + dx(i+1) )  )
!  Qx(i) = Kcond * rtTx * ( Tn(i+1)-Tn(i) )*2.d0 /( dx(i) + dx(i+1) )
!end do
!do i = 1, Ncell
!!----- Cooling ---------------------------------------------------------
!  Call Fcool( CooL,N(i),Tn(i))
!  gammi1 =   3.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+5.d0*ndH2(i)
!  gammi1 = ( 2.d0*(ndH(i)+ndp(i)+ndHe(i)+ndHep(i))+2.d0*ndH2(i) )/gammi1
!  if( dt .le. 0.2d0*Pn(i)/(gammi1*dabs(CooL)) ) then
!    U(i,5) = U(i,5) - gammi1*CooL*dt*0.5d0 !explicit
!  else
!    Call IMC( U(i,5),ndH(i)+ndp(i)+ndHe(i)+ndHep(i)+ndH2(i),dt*0.5d0,N(i) ) !implicit
!  end if
!!----- Conduction ------------------------------------------------------
!  U(i,5) = U(i,5) + gammi1*dt*0.5d0*(Qx(i)-Qx(i-1))/dx(i)
!  U(i,5) = dmax1(U(i,5),pmin)
!  U(i,5) = dmin1(U(i,5),pmax)
!end do


SUBROUTINE collectST()
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: stMPI(1:valstar,0:numstar,0:NPE-1),u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

u1(:,:)=Ustar(:,:)!0.d0


do i=0,numstar!nid+nidnw
do j=1,valstar
  stMPI(j,i,NRANK)=Ustar(j,i)
  !write(*,*) NRANK,Ustar(j,i),'totst-1'
end do;end do

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

do Nroot=0,NPE-1
  CALL MPI_BCAST(stMPI(1,0,Nroot),(valstar)*(numstar+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

totst=0
do Nroot=0,NPE-1
totst=totst+int(stMPI(9,0,Nroot))
enddo
!write(*,*) NRANK,totst,'totst'
!nidnw=0

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

stcnt=nid
do Nroot=0,NPE-1
if (int(stMPI(9,0,Nroot)).ne.0) then
do i=1,int(stMPI(9,0,Nroot))
stcnt=stcnt+1
do iii=1,valstar
u1(iii,stcnt)=stMPI(iii,nid+i,Nroot)
!write(*,*) NRANK,u1(iii,stcnt),int(stMPI(9,0,Nroot)),iii,stcnt,nid+i,'totst12'
enddo
enddo
endif
enddo

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

nid=nid+totst
if(totst.ne.0) then
do i=nid-totst+1,nid; do j=1,valstar!numstar
  Ustar(j,i)=u1(j,i)
  Ustar(9,i)=dble(i)
!write(*,*) NRANK,Ustar(j,i),u1(j,i),nid,j,i,'totst2'
end do;end do
endif

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!nid=nid+totst
!stMPI(9,0,:)=0.d0

END SUBROUTINE collectST


SUBROUTINE BC_ST()
USE comvar
USE mpivar
USE slfgrv
USE fedvar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: stMPI(1:valstar,0:numstar,0:NPE-1),u1(1:valstar,0:numstar)
integer :: num1,totst,stcnt,nwid,NRANKdm
integer :: rsix,rsjy,rskz

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

do j=1,valstar; do i=0,numstar
  stMPI(j,i,NRANK)=Ustar(j,i)
  u1(j,i)=Ustar(j,i)
end do;end do

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
do Nroot=0,NPE-1
  CALL MPI_BCAST(stMPI(1,0,Nroot),(valstar)*(numstar+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

do i=1,numstar
rsix=int((Ustar(1,i))/(dble(Ncellx)*dx1))
rsjy=int((Ustar(2,i))/(dble(Ncelly)*dy1))
rskz=int((Ustar(3,i))/(dble(Ncellz)*dz1))
NRANKdm = rsix + rsjy*NSPLTx + rskz*NSPLTx*NSPLTy
do j=1,valstar
u1(j,i)=stMPI(j,i,NRANKdm)
end do
end do

END SUBROUTINE BC_ST


!SUBROUTINE BC_ST_RS()
!USE comvar
!USE mpivar
!USE slfgrv
!USE fedvar
!INCLUDE 'mpif.h'
!INTEGER :: MSTATUS(MPI_STATUS_SIZE)
!double precision :: stMPI(1:valstar,0:numstar,0:NPE-1),u1(1:valstar,0:numstar)
!integer :: num1,totst,stcnt,nwid,NRANKdm
!integer :: rsix,rsjy,rskz

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!j=14
!do i=0,numstar
!  stMPI(j,i,NRANK)=Ustar(j,i)
!  u1(j,i)=Ustar(j,i)
!end do

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!do Nroot=0,NPE-1
!  CALL MPI_BCAST(stMPI(14,0,Nroot),(numstar+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
!end do

!do i=1,numstar
!do Nroot=1,NPE-1
!  stMPI(14,i,Nroot)=stMPI(14,i,Nroot)+tMPI(14,i,Nroot-1)
!end do
!end do


!do i=1,numstar
!  Ustar(14,i)=stMPI(14,i,NPE-1)
!end do

!do i=1,numstar
!rsix=int((Ustar(1,i))/(dble(Ncellx)*dx1))
!rsjy=int((Ustar(2,i))/(dble(Ncelly)*dy1))
!rskz=int((Ustar(3,i))/(dble(Ncellz)*dz1))
!NRANKdm = rsix + rsjy*NSPLTx + rskz*NSPLTx*NSPLTy
!do j=1,valstar
!u1(j,i)=stMPI(j,i,NRANKdm)
!end do
!end do

!END SUBROUTINE BC_ST_RS
