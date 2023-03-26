
PROGRAM MAIN
INCLUDE 'mpif.h'
real*8, dimension(:), allocatable :: x_i,y_i,z_i,dx_i,dy_i,dz_i,x,y,z,dx,dy,dz
real*8, dimension(:,:,:,:), allocatable :: U,Ak
real*8 :: theta,pi,amp,xpi,ypi,zpi,phase1,phase2,phase3,kx,ky,kzz,kw,radius,mass,vir,sigmav,kw1,kmin,theta1,theta2,theta3
real*8 :: ampn(2048,3),ampn0(2048,3),md(2048,3),md0(2048,3)
real*8 :: ql1x,ql1y,ql1z,ql2x,ql2y,ql2z,dinit1,amp1,Lbox
character*3 :: NPENUM
DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
double precision :: Msun=2.4d-2
INTEGER :: MSTATUS(MPI_STATUS_SIZE)


CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPE  ,IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,NRANK,IERR)
if(NPE.eq.4)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 1; end if
if(NPE.eq.8)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 2; end if
if(NPE.eq.16)   then; NSPLTx = 4; NSPLTy = 2; NSPLTz = 2; end if
if(NPE.eq.32)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 2; end if
if(NPE.eq.64)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.128)  then; NSPLTx = 8; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.256)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 4; end if
if(NPE.eq.512)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 8; end if
if(NPE.eq.1024) then; NSPLTx = 8; NSPLTy = 8; NSPLTz =16; end if

IST = mod(NRANK,NSPLTx); KST = NRANK/(NSPLTx*NSPLTy); JST = NRANK/NSPLTx-NSPLTy*KST
LEFT = NRANK - 1            ; if(IST.eq.0       ) LEFT = NRANK + (NSPLTx-1)
RIGT = NRANK + 1            ; if(IST.eq.NSPLTx-1) RIGT = NRANK - (NSPLTx-1)
BOTM = NRANK - NSPLTx       ; if(JST.eq.0       ) BOTM = NRANK + NSPLTx*(NSPLTy-1)
TOP  = NRANK + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = NRANK - NSPLTx*(NSPLTy-1)
DOWN = NRANK - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = NRANK + NSPLTx*NSPLTy*(NSPLTz-1)
UP   = NRANK + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = NRANK - NSPLTx*NSPLTy*(NSPLTz-1)
!***** x-direction shock tube test *****!

Ncellx = 128; Ncelly = 128; Ncellz = 128
ql1x = 40.d0; ql2x = 40.d0; ql1y = 40.d0; ql2y = 40.d0; ql1z = 40.d0; ql2z = 40.d0


!WNM 0.3311311e+00  0.2201655e+04  0.6648891e+04
!dinit1 = 0.3311311d0*1.27d0
!TUE 0.2968247e+01  0.2202584e+04  0.7420489e+03
pi = 3.14159265358979323846d0
dinit1 = 86.d0!*2.d0
radius=ql1x*0.5d0
mass=1.d5*Msun
vir=1.d0
!amp1=vir*3.d0/5.d0*G*mass/radius/dsqrt(3.d0)
amp1=dsqrt(vir*3.d0/5.d0*G*mass/radius)*5.d-1
!sigmav=1.d0/8.d0/3.141592d0/2.d0/3.141592d0
sigmav=1.d0/64.d0/3.141592d0/2.d0/3.141592d0
Lbox=ql1x+ql2x
kmin=1.d0!2.d0*pi/Lbox
!kmin=2.d0*pi*Lbox/Lbox


Ncellx = Ncellx/NSPLTx; Ncelly = Ncelly/NSPLTy; Ncellz = Ncellz/NSPLTz
ALLOCATE(dx_i(-1:Ncellx*NSPLTx+2)); ALLOCATE(dy_i(-1:Ncelly*NSPLTy+2)); ALLOCATE(dz_i(-1:Ncellz*NSPLTz+2))
ALLOCATE( x_i(-1:Ncellx*NSPLTx+2)); ALLOCATE( y_i(-1:Ncelly*NSPLTy+2)); ALLOCATE( z_i(-1:Ncellz*NSPLTz+2))
ALLOCATE(dx(-1:Ncellx+2)); ALLOCATE(dy(-1:Ncelly+2)); ALLOCATE(dz(-1:Ncellz+2))
ALLOCATE( x(-1:Ncellx+2)); ALLOCATE( y(-1:Ncelly+2)); ALLOCATE( z(-1:Ncellz+2))
ALLOCATE(U(Ncellx+1,Ncelly+1,Ncellz+1,3))
ALLOCATE(Ak(Ncellx+1,Ncelly+1,Ncellz+1,3))

do i = -1, Ncellx*NSPLTx+2
  dx_i(i) = (ql1x+ql2x)/(Ncellx*NSPLTx)
end do
do j = -1, Ncelly*NSPLTy+2
  dy_i(j) = (ql1y+ql2y)/(Ncelly*NSPLTy)
end do
do k = -1, Ncellz*NSPLTz+2
  dz_i(k) = (ql1z+ql2z)/(Ncellz*NSPLTz)
end do

x_i(-1) = -dx_i(0)
do i = 0, Ncellx*NSPLTx+2
   x_i(i) = x_i(i-1) + dx_i(i)
end do
y_i(-1) = -dy_i(0)
do j = 0, Ncelly*NSPLTy+2
   y_i(j) = y_i(j-1) + dy_i(j)
end do
z_i(-1) = -dz_i(0)
do k = 0, Ncellz*NSPLTz+2
   z_i(k) = z_i(k-1) + dz_i(k)
end do

do i = -1, Ncellx+2
  ix    =  IST*Ncellx + i
  x(i)  =  x_i(ix)
  dx(i) =  dx_i(ix)
end do
do j = -1, Ncelly+2
  jy    =  JST*Ncelly + j
  y(j)  =  y_i(jy)
  dy(j) =  dy_i(jy)
end do
do k = -1, Ncellz+2
  kz    =  KST*Ncellz + k
  z(k)  =  z_i(kz)
  dz(k) =  dz_i(kz)
end do

IF(NRANK.EQ.0) THEN
  400 format(D25.17)
  open(4,file='/work/maedarn/3DMHD/samplecnv/cdnt.DAT')
    write(4,400) ( 0.5d0 * ( x_i(i-1)+x_i(i) ), i=1, Ncellx*NSPLTx )
    write(4,400) ( 0.5d0 * ( y_i(j-1)+y_i(j) ), j=1, Ncelly*NSPLTy )
    write(4,400) ( 0.5d0 * ( z_i(k-1)+z_i(k) ), k=1, Ncellz*NSPLTz )
  close(4)
END IF





  U(:,:,:,:) = 0.d0
  idum = 1!;jdum = 20;kdum = 300;
  kmax = 64
  do l =-kmax,kmax; if(NRANK.eq.0) write(*,*) l; do ll= -kmax,kmax; do lll=-kmax,kmax
    kw  = dsqrt( dble(l)**2.d0+dble(ll)**2.d0+dble(lll)**2.d0)
    if((l.ne.0).or.(ll.ne.0).or.(lll.ne.0)) then; if(kw.le.dble(kmax)) then; if(kw.ge.dble(2)) then
      !amp = dsqrt( kw**(-11.d0/3.d0) )
      !amp = dsqrt( kw**(-4.d0 -2.d0) )
      !kx=2.d0*pi*dble(l); ky=2.d0*pi*dble(ll); kzz=2.d0*pi*dble(lll)
      kx=2.d0*pi*dble(l); ky=2.d0*pi*dble(ll); kzz=2.d0*pi*dble(lll)
      kw1=(kx*kx+ky*ky+kz*kz)/Lbox/Lbox
      !kw1=
      !kw1=dsqrt(kw1)
      !amp=kw1*sigmav*dexp(-kw1**2.d0/2.d0/sigmav) ** dsqrt( kw1**(-4.d0 -2.d0) )
      !amp=kw1*sigmav*dexp(-kw1**2.d0/2.d0/sigmav) ** dsqrt( kw1**(-4.d0 -2.d0) )
      !amp= dsqrt( (kw1+2.d0*pi)**(-4.d0 -2.d0) )
      !amp= dsqrt( (kw1+kmin**2.d0)**((-4.d0 -2.d0)/2.d0) )
      !amp= dsqrt( (kw +1.d0**2.d0)**((-4.d0 -2.d0)/2.d0) )
      amp= dsqrt( (kw)**((-4.d0 -2.d0)/2.d0) )
      !amp=dqsrt(amp)
      !kx=2.d0*pi*dble(l); ky=2.d0*pi*dble(ll); kzz=2.d0*pi*dble(lll)
      call ran0(phase1,idum); phase1=2.d0*pi*phase1
      call ran0(phase2,idum); phase2=2.d0*pi*phase2
      call ran0(phase3,idum); phase3=2.d0*pi*phase3
      !if(NRANK.eq.0) write(*,*) kw,amp,kw1,dsqrt( (kw1+kmin)**(-4.d0 -2.d0) ),kmin,dsqrt( (kmin*2.d0)**(-2.d0 -1.d0) ),'amp'
      do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
        xpi=0.5d0*(x(i-1)+x(i))/(ql1x+ql2x);ypi=0.5d0*(y(j-1)+y(j))/(ql1y+ql2y);zpi=0.5d0*(z(k-1)+z(k))/(ql1z+ql2z)
        !kw1=kx*kx*xpi*xpi+ky*ky*ypi*ypi+kz*kz*zpi*zpi
        !amp=kw1*sigmav*dsqrt(-kw1**2.d0/2.d0/sigmav) ** dsqrt( kw1**(-4.d0 -2.d0) )
        Ak(i,j,k,1) = amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase1)!sin(x+pi/2) ik * Ak
        Ak(i,j,k,2) = amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase2)
        Ak(i,j,k,3) = amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase3)

        !U(i,j,k,1) = U(i,j,k,1)+(ky/Lbox*Ak(i,j,k,3)-kz/Lbox*Ak(i,j,k,2))*amp1
        !U(i,j,k,2) = U(i,j,k,2)+(kz/Lbox*Ak(i,j,k,1)-kx/Lbox*Ak(i,j,k,3))*amp1
        !U(i,j,k,3) = U(i,j,k,3)+(kx/Lbox*Ak(i,j,k,2)-ky/Lbox*Ak(i,j,k,1))*amp1
        U(i,j,k,1) = U(i,j,k,1)+(ky*Ak(i,j,k,3)-kz*Ak(i,j,k,2))*amp1/(2.d0*pi)
        U(i,j,k,2) = U(i,j,k,2)+(kz*Ak(i,j,k,1)-kx*Ak(i,j,k,3))*amp1/(2.d0*pi)
        U(i,j,k,3) = U(i,j,k,3)+(kx*Ak(i,j,k,2)-ky*Ak(i,j,k,1))*amp1/(2.d0*pi)
      end do;end do;end do
    end if; end if; end if
  end do;end do;end do

!****ABE****
!Turbulence
!idum = 1;pi22 = 6.28318530717959d0;vdisx=0.d0;vdisy=0.d0;vdisz=0.d0
!amp = (0.078d0 + 1.003077577d0/8.d0 + 0.9d0)/1.1d0
!do kx = -10,10; do ky = -10,10; do kz = -10,10
!  if((kx.eq.0).and.(ky.eq.0).and.(kz.eq.0)) goto 5328
!    kw = sqrt(real(kx)**2+real(ky)**2+real(kz)**2)
!    kr = (kx*x(i)+ky*y(j)+kz*z(k))/(2.d0*rboundary)
!   ! kw**(-11/3) : Kolmogorov,  kw**(-12/3) : lawson's law
!    call ran0(rndm,idum)
!    vdisx = vdisx +amp*sqrt( kw**(-12.d0/3.d0) )*sin(pi22*(kr+rndm))
!    call ran0(rndm,idum)
!    vdisy = vdisy +amp*sqrt( kw**(-12.d0/3.d0) )*sin(pi22*(kr+rndm))
!    call ran0(rndm,idum)
!    vdisz = vdisz +amp*sqrt( kw**(-12.d0/3.d0) )*sin(pi22*(kr+rndm))
!  5328 continue
!enddo;enddo;enddo
!vdis2 = vdis2 + vdisx**2+vdisy**2+vdisz**2
!clnmb = clnmb + 1
!vx(i,j,k) = vx(i,j,k) + vdisx!*0.5d0*( 1.d0-tanh((sqrt(r2_sc)-r_sc)/0.1d0) ) *exp(-r2_sc/(r_sc**2))
!vy(i,j,k) = vy(i,j,k) + vdisy!*0.5d0*( 1.d0-tanh((sqrt(r2_sc)-r_sc)/0.1d0) ) *exp(-r2_sc/(r_sc**2))
!vz(i,j,k) = vz(i,j,k) + vdisz!*0.5d0*( 1.d0-tanh((sqrt(r2_sc)-r_sc)/0.1d0) ) *exp(-r2_sc/(r_sc**2))
!****ABE****

do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
do l =-kmax,kmax; if(NRANK.eq.0) write(*,*) l; do ll= -kmax,kmax; do lll=-kmax,kmax
  kw  = dsqrt( dble(l)**2.d0+dble(ll)**2.d0+dble(lll)**2.d0)
  amp= dsqrt( (kw)**((-4.d0 -2.d0)/2.d0) )
  if((l.ne.0).or.(ll.ne.0).or.(lll.ne.0)) then; if(kw.le.dble(kmax)) then; if(kw.ge.dble(2)) then
   ! kr = (kx*x(i)+ky*y(j)+kz*z(k))/(2.d0*rboundary)
   ! kw**(-11/3) : Kolmogorov,  kw**(-12/3) : lawson's law
    xpi=0.5d0*(x(i-1)+x(i))/(ql1x+ql2x);ypi=0.5d0*(y(j-1)+y(j))/(ql1y+ql2y);zpi=0.5d0*(z(k-1)+z(k))/(ql1z+ql2z)
    call ran0(phase1,idum); phase1=2.d0*pi*phase1
    call ran0(phase2,idum); phase2=2.d0*pi*phase2
    call ran0(phase3,idum); phase3=2.d0*pi*phase3
    Ak(i,j,k,1) =Ak(i,j,k,1)+ amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase1)!sin(x+pi/2) ik * Ak
    Ak(i,j,k,2) =Ak(i,j,k,2)+ amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase2)
    Ak(i,j,k,3) =Ak(i,j,k,3)+ amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase3)
    !call ran0(rndm,idum)
   ! vdisx = vdisx +amp*sqrt( kw**(-12.d0/3.d0) )*sin(pi22*(kr+rndm))
    !call ran0(rndm,idum)
   ! vdisy = vdisy +amp*sqrt( kw**(-12.d0/3.d0) )*sin(pi22*(kr+rndm))
    !call ran0(rndm,idum)
   ! vdisz = vdisz +amp*sqrt( kw**(-12.d0/3.d0) )*sin(pi22*(kr+rndm))
   end if; end if; end if
enddo;enddo;enddo
end do;end do;end do

do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
do l =-kmax,kmax; if(NRANK.eq.0) write(*,*) l; do ll= -kmax,kmax; do lll=-kmax,kmax
  kw  = dsqrt( dble(l)**2.d0+dble(ll)**2.d0+dble(lll)**2.d0)
  amp= dsqrt( (kw)**((-4.d0 -2.d0)/2.d0) )
  if((l.ne.0).or.(ll.ne.0).or.(lll.ne.0)) then; if(kw.le.dble(kmax)) then; if(kw.ge.dble(2)) then

    U(i,j,k,1) = U(i,j,k,1)+(ky*Ak(i,j,k,3)-kz*Ak(i,j,k,2))*amp1/(2.d0*pi)
    U(i,j,k,2) = U(i,j,k,2)+(kz*Ak(i,j,k,1)-kx*Ak(i,j,k,3))*amp1/(2.d0*pi)
    U(i,j,k,3) = U(i,j,k,3)+(kx*Ak(i,j,k,2)-ky*Ak(i,j,k,1))*amp1/(2.d0*pi)
    
   end if; end if; end if
enddo;enddo;enddo
end do;end do;end do

!vdis2 = vdis2 + vdisx**2+vdisy**2+vdisz**2
!clnmb = clnmb + 1
!vx(i,j,k) = vx(i,j,k) + vdisx!*0.5d0*( 1.d0-tanh((sqrt(r2_sc)-r_sc)/0.1d0) ) *exp(-r2_sc/(r_sc**2))
!vy(i,j,k) = vy(i,j,k) + vdisy!*0.5d0*( 1.d0-tanh((sqrt(r2_sc)-r_sc)/0.1d0) ) *exp(-r2_sc/(r_sc**2))
!vz(i,j,k) = vz(i,j,k) + vdisz!*0.5d0*( 1.d0-tanh((sqrt(r2_sc)-r_sc)/0.1d0) ) *exp(-r2_sc/(r_sc**2))

U(:,:,:,:) = 0.d0
idum = 1!;jdum = 20;kdum = 300;
kmax = 64
do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1

do l =-kmax,kmax; if(NRANK.eq.0) write(*,*) l; do ll= -kmax,kmax; do lll=-kmax,kmax
  kw  = dsqrt( dble(l)**2.d0+dble(ll)**2.d0+dble(lll)**2.d0)
  if((l.ne.0).or.(ll.ne.0).or.(lll.ne.0)) then; if(kw.le.dble(kmax)) then; if(kw.ge.dble(2)) then
    kx=2.d0*pi*dble(l); ky=2.d0*pi*dble(ll); kzz=2.d0*pi*dble(lll)
    kw1=(kx*kx+ky*ky+kz*kz)/Lbox/Lbox
    amp= dsqrt( (kw)**((-4.d0 -2.d0)/2.d0) )
    call ran0(phase1,idum); phase1=2.d0*pi*phase1
    call ran0(phase2,idum); phase2=2.d0*pi*phase2
    call ran0(phase3,idum); phase3=2.d0*pi*phase3
    !do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
      xpi=0.5d0*(x(i-1)+x(i))/(ql1x+ql2x);ypi=0.5d0*(y(j-1)+y(j))/(ql1y+ql2y);zpi=0.5d0*(z(k-1)+z(k))/(ql1z+ql2z)
      Ak(i,j,k,1) = amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase1)!sin(x+pi/2) ik * Ak
      Ak(i,j,k,2) = amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase2)
      Ak(i,j,k,3) = amp*dsin(kx*xpi+ky*ypi+kzz*zpi+phase3)
      U(i,j,k,1) = U(i,j,k,1)+(ky*Ak(i,j,k,3)-kz*Ak(i,j,k,2))*amp1/(2.d0*pi)
      U(i,j,k,2) = U(i,j,k,2)+(kz*Ak(i,j,k,1)-kx*Ak(i,j,k,3))*amp1/(2.d0*pi)
      U(i,j,k,3) = U(i,j,k,3)+(kx*Ak(i,j,k,2)-ky*Ak(i,j,k,1))*amp1/(2.d0*pi)
  end if; end if; end if
end do;end do;end do

!U(i,j,k,1) = U(i,j,k,1)+(ky*Ak(i,j,k,3)-kz*Ak(i,j,k,2))*amp1/(2.d0*pi)
!U(i,j,k,2) = U(i,j,k,2)+(kz*Ak(i,j,k,1)-kx*Ak(i,j,k,3))*amp1/(2.d0*pi)
!U(i,j,k,3) = U(i,j,k,3)+(kx*Ak(i,j,k,2)-ky*Ak(i,j,k,1))*amp1/(2.d0*pi)
end do;end do;end do


!ampn0(NRANK)=0.d0
!do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
!  ampn0(NRANK) = ampn0(NRANK) + U(i,j,k,1)**2
!end do;end do;end do
!CALL MPI_GATHER(ampn0(NRANK),1,MPI_REAL8,ampn,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
!amp = 0.d0
!if(NRANK.eq.0) then
!  do NN=0,NPE-1
!    amp = amp + ampn(NN) !集める
!  end do; amp = amp/NPE/Ncellx/Ncelly/Ncellz !平均
!  amp = 0.1d0/dsqrt(amp)
!end if
!CALL MPI_BCAST(amp,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
!do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
!  U(i,j,k)=dinit1*dexp(amp*U(i,j,k))
!end do;end do;end do

!ampn0(NRANK,:)=0.d0; md0(NRANK,:)=0.d0
ampn0(:,:)=0.d0; md0(:,:)=0.d0
ampn(:,:)=0.d0; md(:,:)=0.d0
do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
  ampn0(NRANK,1) = ampn0(NRANK,1) + U(i,j,k,1)**2
    md0(NRANK,1) =   md0(NRANK,1) + U(i,j,k,1)
  ampn0(NRANK,2) = ampn0(NRANK,2) + U(i,j,k,2)**2
    md0(NRANK,2) =   md0(NRANK,2) + U(i,j,k,2)
  ampn0(NRANK,3) = ampn0(NRANK,3) + U(i,j,k,3)**2
    md0(NRANK,3) =   md0(NRANK,3) + U(i,j,k,3)
end do;end do;end do
CALL MPI_GATHER(ampn0(NRANK,1),1,MPI_REAL8,ampn,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_GATHER(  md0(NRANK,1),1,MPI_REAL8,  md,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_GATHER(ampn0(NRANK,2),1,MPI_REAL8,ampn,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_GATHER(  md0(NRANK,2),1,MPI_REAL8,  md,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_GATHER(ampn0(NRANK,3),1,MPI_REAL8,ampn,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_GATHER(  md0(NRANK,3),1,MPI_REAL8,  md,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
amp = 0.d0; theta=0.d0
theta1=0.d0
theta2=0.d0
theta3=0.d0
if(NRANK.eq.0) then
  do NN=0,NPE-1
    amp   = amp   + ampn(NN,1)+ ampn(NN,2)+ ampn(NN,3)
    !theta = theta +   md(NN,1)+   md(NN,2)+   md(NN,3)
    theta1 = theta1 +   md(NN,1)
    theta2 = theta2 +   md(NN,2)
    theta3 = theta3 +   md(NN,3)
  end do; amp = amp/NPE/Ncellx/Ncelly/Ncellz; theta = theta/NPE/Ncellx/Ncelly/Ncellz
  theta1 = theta1/NSPLTx/Ncellx
  theta2 = theta2/NSPLTy/Ncelly
  theta3 = theta3/NSPLTz/Ncellz
  write(*,*) 'mean=',theta1,theta2,theta3
  write(*,*) 'disp=',dsqrt(amp-theta1**2-theta2**2-theta3**2)
  write(*,*) dsqrt(amp-theta**2)/theta
end if


!--------------------------------

!***** WRITE Initial Conditions *****!
  WRITE(NPENUM,'(I3.3)') NRANK
  open(unit=8,file='/work/maedarn/3DMHD/VTF/V'//NPENUM//'.dat',FORM='UNFORMATTED')
  do k = 1, Ncellz
  do j = 1, Ncelly
    write(8) (sngl(U(i,j,k,1)),sngl(U(i,j,k,2)),sngl(U(i,j,k,3)),i=1,Ncellx)
  end do
  end do
  close(8)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END PROGRAM


SUBROUTINE ran0(ran,idum)
INTEGER idum,IA,IM,IQ,IR,MASK
REAL*8 ran,AM
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,MASK=123459876)
INTEGER k
!idum=ieor(idum,MASK)
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if(idum.lt.0) idum=idum+IM
ran=AM*idum
!idum=ieor(idum,MASK)
END SUBROUTINE ran0
