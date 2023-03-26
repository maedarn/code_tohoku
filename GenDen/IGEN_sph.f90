
PROGRAM MAIN
INCLUDE 'mpif.h'
real*8, dimension(:), allocatable :: x_i,y_i,z_i,dx_i,dy_i,dz_i,x,y,z,dx,dy,dz
real*8, dimension(:,:,:,:), allocatable :: U
real*8 :: theta,pi,amp,xpi,ypi,zpi,phase1,phase2,phase3,kx,ky,kzz,kw,radius,mass,vir
real*8 :: ampn(2048),ampn0(2048),md(2048),md0(2048)
real*8 :: ql1x,ql1y,ql1z,ql2x,ql2y,ql2z,dinit1,amp1
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

Ncellx = 256; Ncelly = 256; Ncellz = 256
ql1x = 5.d0; ql2x = 5.d0; ql1y = 5.d0; ql2y = 5.d0; ql1z = 5.d0; ql2z = 5.d0


!WNM 0.3311311e+00  0.2201655e+04  0.6648891e+04
!dinit1 = 0.3311311d0*1.27d0
!TUE 0.2968247e+01  0.2202584e+04  0.7420489e+03
dinit1 = 86.d0*1.27d0
radius=20.d0
mass=1.d0*Msun
vir=2.d0
amp1=vir*3.d0/5.d0*G*mass/radius/dsqrt(3.d0)

Ncellx = Ncellx/NSPLTx; Ncelly = Ncelly/NSPLTy; Ncellz = Ncellz/NSPLTz
ALLOCATE(dx_i(-1:Ncellx*NSPLTx+2)); ALLOCATE(dy_i(-1:Ncelly*NSPLTy+2)); ALLOCATE(dz_i(-1:Ncellz*NSPLTz+2))
ALLOCATE( x_i(-1:Ncellx*NSPLTx+2)); ALLOCATE( y_i(-1:Ncelly*NSPLTy+2)); ALLOCATE( z_i(-1:Ncellz*NSPLTz+2))
ALLOCATE(dx(-1:Ncellx+2)); ALLOCATE(dy(-1:Ncelly+2)); ALLOCATE(dz(-1:Ncellz+2))
ALLOCATE( x(-1:Ncellx+2)); ALLOCATE( y(-1:Ncelly+2)); ALLOCATE( z(-1:Ncellz+2))
ALLOCATE(U(Ncellx+1,Ncelly+1,Ncellz+1,3))

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


  pi = 3.14159265358979323846d0
  U(:,:,:,:) = 0.d0
  idum = 1;jdum = 2;kdum = 3; kmax = 64
  do l =-kmax,kmax; if(NRANK.eq.0) write(*,*) l; do ll= -kmax,kmax; do lll=-kmax,kmax
    kw  = dsqrt( dble(l**2+ll**2+lll**2) )
    if((l.ne.0).or.(ll.ne.0).or.(lll.ne.0)) then; if(kw.le.dble(kmax)) then; if(kw>dble(2)) then
      !amp = dsqrt( kw**(-11.d0/3.d0) )
      amp = dsqrt( kw**(-4.d0) )
      amp=dqsrt(amp)
      kx=2.d0*pi*dble(l); ky=2.d0*pi*dble(ll); kzz=2.d0*pi*dble(lll)
      call ran0(phase1,idum); phase1=2.d0*pi*phase1
      call ran0(phase2,jdum); phase2=2.d0*pi*phase2
      call ran0(phase3,kdum); phase3=2.d0*pi*phase3
      do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
        xpi=0.5d0*(x(i-1)+x(i))/(ql1x+ql2x);ypi=0.5d0*(y(j-1)+y(j))/(ql1y+ql2y);zpi=0.5d0*(z(k-1)+z(k))/(ql1z+ql2z)
        U(i,j,k,1) = U(i,j,k,1)+amp*dcos(kx*xpi+ky*ypi+kzz*zpi+phase1)
        U(i,j,k,2) = U(i,j,k,2)+amp*dcos(kx*xpi+ky*ypi+kzz*zpi+phase2)
        U(i,j,k,3) = U(i,j,k,3)+amp*dcos(kx*xpi+ky*ypi+kzz*zpi+phase3)
      end do;end do;end do
    end if; end if; end if
  end do;end do;end do

  
  ampn0(NRANK)=0.d0
  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
    ampn0(NRANK) = ampn0(NRANK) + U(i,j,k,1)**2+ U(i,j,k,2)**2+ U(i,j,k,3)**2
  end do;end do;end do
  CALL MPI_GATHER(ampn0(NRANK),1,MPI_REAL8,ampn,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
  amp = 0.d0
  if(NRANK.eq.0) then
    do NN=0,NPE-1
      amp = amp + ampn(NN) !集める
    end do; amp = amp/NPE/Ncellx/Ncelly/Ncellz !平均
  !amp = 2.d0*3.d0/5.d0*G/radius*mass!0.1d0/dsqrt(amp)
  end if
  CALL MPI_BCAST(amp,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
    U(i,j,k)=dinit1*dexp(amp*U(i,j,k))
  end do;end do;end do


  ampn0(NRANK)=0.d0; md0(NRANK)=0.d0
  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
    ampn0(NRANK) = ampn0(NRANK) + U(i,j,k)**2
      md0(NRANK) =   md0(NRANK) + U(i,j,k)
  end do;end do;end do
  CALL MPI_GATHER(ampn0(NRANK),1,MPI_REAL8,ampn,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
  CALL MPI_GATHER(  md0(NRANK),1,MPI_REAL8,  md,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
  amp = 0.d0; theta=0.d0
  if(NRANK.eq.0) then
    do NN=0,NPE-1
      amp   = amp   + ampn(NN)
      theta = theta +   md(NN)
    end do; amp = amp/NPE/Ncellx/Ncelly/Ncellz; theta = theta/NPE/Ncellx/Ncelly/Ncellz
    write(*,*) 'mean=',theta
    write(*,*) 'disp=',dsqrt(amp-theta**2)
    write(*,*) dsqrt(amp-theta**2)/theta
  end if


!--------------------------------

!***** WRITE Initial Conditions *****!
  WRITE(NPENUM,'(I3.3)') NRANK
  open(unit=8,file='/work/maedarn/3DMHD/samplecnv/DTF/D'//NPENUM//'.dat',FORM='UNFORMATTED')
  do k = 1, Ncellz
  do j = 1, Ncelly
    write(8) (sngl(U(i,j,k)),i=1,Ncellx)
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
