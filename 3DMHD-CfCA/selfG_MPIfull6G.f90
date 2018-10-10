MODULE comvar
!INTEGER, parameter :: ndx=130, ndy=130, ndz=130, ndmax=130, Dim=3 !1024^3
INTEGER, parameter :: ndx=66, ndy=66, ndz=66, ndmax=66, Dim=3 !512^3
!INTEGER, parameter :: ndx=34, ndy=34, ndz=34, ndmax=34, Dim=3
DOUBLE PRECISION, dimension(-1:ndx) :: x,dx
DOUBLE PRECISION, dimension(-1:ndy) :: y,dy
DOUBLE PRECISION, dimension(-1:ndz) :: z,dz
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U, Bcc, Blg, Vfc, EMF
DOUBLE PRECISION, dimension(:,:,:),   allocatable :: dnc, xlag, dxlagM

DOUBLE PRECISION, parameter :: kb=8.63359d0, Kcond=1.6384d-2
DOUBLE PRECISION  :: gamma,gammi1,gammi2,gammi3,gampl1,gampl2,gampl3
DOUBLE PRECISION  :: CFL,facdep,tfinal,time,phr(-1:400)
DOUBLE PRECISION  :: pmin,pmax,rmin,rmax
INTEGER :: Ncellx,Ncelly,Ncellz,iwx,iwy,iwz,maxstp,nitera
INTEGER :: ifchem,ifthrm,ifrad,ifgrv
END MODULE comvar

MODULE mpivar
INTEGER :: NPE,NRANK, NSPLTx,NSPLTy,NSPLTz, IST,JST,KST, LEFT,RIGT,BOTM,TOP,UP,DOWN
INTEGER :: BCx1,BCx2,BCy1,BCy2,BCz1,BCz2, N_MPI(20)
DOUBLE PRECISION  :: BBRV(10,2,2),BBRV_cm(8)
REAL*4, dimension(:,:,:), allocatable :: DTF
END MODULE mpivar

MODULE chmvar
DOUBLE PRECISION, parameter :: mH=1.d0, mHe=4.d0, mH2=2.d0, mC=12.d0, mCO=28.d0
DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-4, xo=3.2d-4, dv=2.d0, Tgr=5.d-3
DOUBLE PRECISION, dimension(:,:,:)  , allocatable :: ndp,ndH,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: Ntot,NH2,NnC,NCO,tCII
DOUBLE PRECISION  :: ndpmin,ndHmin,ndH2min,ndHemin,ndHepmin,ndCmin,ndCpmin,ndCOmin
END MODULE chmvar

MODULE slfgrv
DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
INTEGER :: point1(0:15),point2(0:15),NGL,NGcr,Nmem1,Nmem2
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Phi
DOUBLE PRECISION :: Lbox

INTEGER :: pointb1(0:15),pointb2(0:15)
DOUBLE PRECISION, dimension(:,:), allocatable :: bphi1,bphi2
END MODULE slfgrv

!======================================================================*
!                                 MAIN                                 *
!======================================================================*

PROGRAM MAIN_3DMHD
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'

CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPE  ,IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,NRANK,IERR)
!NPE=NPE/40
!NRANK=NRANK/40
!write(*,*) NPE,NRANK
!write(*,*) 'OK'
!----- Prepare MPI SPLIT -----------------------------------------------!

if(NPE.eq.4)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 1; end if
if(NPE.eq.8)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 2; end if
if(NPE.eq.16)   then; NSPLTx = 4; NSPLTy = 2; NSPLTz = 2; end if
if(NPE.eq.32)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 2; end if
if(NPE.eq.64)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.128)  then; NSPLTx = 8; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.256)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 4; end if
if(NPE.eq.512)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 8; end if
if(NPE.eq.1024) then; NSPLTx = 8; NSPLTy = 8; NSPLTz =16; end if

!write(*,*) 'OK1'

IST = mod(NRANK,NSPLTx); KST = NRANK/(NSPLTx*NSPLTy); JST = NRANK/NSPLTx-NSPLTy*KST
!write(*,*) 'OK11'
LEFT = NRANK - 1            ; if(IST.eq.0       ) LEFT = NRANK + (NSPLTx-1)
RIGT = NRANK + 1            ; if(IST.eq.NSPLTx-1) RIGT = NRANK - (NSPLTx-1)
BOTM = NRANK - NSPLTx       ; if(JST.eq.0       ) BOTM = NRANK + NSPLTx*(NSPLTy-1)
TOP  = NRANK + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = NRANK - NSPLTx*(NSPLTy-1)
DOWN = NRANK - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = NRANK + NSPLTx*NSPLTy*(NSPLTz-1)
UP   = NRANK + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = NRANK - NSPLTx*NSPLTy*(NSPLTz-1)
!----------------------------------------------------------------------!
!write(*,*) 'OK2'
ALLOCATE( U(-1:ndx,-1:ndy,-1:ndz,8) )
ALLOCATE(ndH(-1:ndx,-1:ndy,-1:ndz),ndp(-1:ndx,-1:ndy,-1:ndz),ndH2(-1:ndx,-1:ndy,-1:ndz),ndHe(-1:ndx,-1:ndy,-1:ndz), &
       ndHep(-1:ndx,-1:ndy,-1:ndz),ndC(-1:ndx,-1:ndy,-1:ndz),ndCp(-1:ndx,-1:ndy,-1:ndz),ndCO(-1:ndx,-1:ndy,-1:ndz), &
         nde(-1:ndx,-1:ndy,-1:ndz),ndtot(-1:ndx,-1:ndy,-1:ndz),Ntot(-1:ndx,-1:ndy,-1:ndz,2),                        &
         NH2(-1:ndx,-1:ndy,-1:ndz,2),NnC(-1:ndx,-1:ndy,-1:ndz,2),NCO(-1:ndx,-1:ndy,-1:ndz,2),tCII(-1:ndx,-1:ndy,-1:ndz,2) )
ALLOCATE(DTF(-1:(ndx-2)*NSPLTx+2,-1:ndy,-1:ndz))
ALLOCATE(Phi(-1:ndx,-1:ndy,-1:ndz))

!write(*,*) 'OK3'

call INITIA
!write(*,*) 'OK'
call EVOLVE

write(*,*) 'OK4'

DEALLOCATE(U)
DEALLOCATE(ndH,ndp,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)
DEALLOCATE(DTF)
DEALLOCATE(Phi)

CALL MPI_FINALIZE(IERR)

END PROGRAM MAIN_3DMHD

!======================================================================*
!                     Prepare an Initial State                         *
!======================================================================*
SUBROUTINE INITIA
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'

integer :: Np1x, Np2x, Np1y, Np2y, Np1z, Np2z, nunit, ix, jy, kz,b,c
double precision ::  ql1x,ql2x,ql1y,ql2y,ql1z,ql2z,dinit1,dinit2,pinit1,pinit2, &
           vinitx1,vinitx2,vinity1,vinity2,vinitz1,vinitz2,           &
           binitx1,binitx2,binity1,binity2,binitz1,binitz2
double precision, dimension(:), allocatable :: x_i,y_i,z_i,dx_i,dy_i,dz_i
double precision :: theta,pi,amp,xpi,ypi,zpi,phase1,phase2,phase3,kx,ky,kzz,kw
double precision :: Hini,pini,H2ini,Heini,Hepini,Cini,COini,Cpini,dBC
double precision :: ampn(2048),ampn0(2048)
character*3 :: NPENUM
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision, dimension(:,:), allocatable :: plane,rand
integer i3,i4,i2y,i2z,rsph2
double precision cenx,ceny,cenz,rsph

open(8,file='/work/maedarn/3DMHD/test/INPUT3D.DAT')
  read(8,*)  Np1x,Np2x
  read(8,*)  Np1y,Np2y
  read(8,*)  Np1z,Np2z
  read(8,*)  ql1x,ql2x
  read(8,*)  ql1y,ql2y
  read(8,*)  ql1z,ql2z
  read(8,*)  vinitx1,vinitx2
  read(8,*)  vinity1,vinity2
  read(8,*)  vinitz1,vinitz2
  read(8,*)  binitx1,binitx2
  read(8,*)  binity1,binity2
  read(8,*)  binitz1,binitz2
  read(8,*)  CFL,facdep
  read(8,*)  maxstp,nitera,tfinal
  read(8,*)  BCx1,BCx2,BCy1,BCy2,BCz1,BCz2
  read(8,*)  ifchem,ifthrm,ifrad,ifgrv
close(8)

!WNM ntot = 1.024
!goto 10000
 pinit1=8.810807d3*kb*1.d-3; pinit2=pinit1
 Hini=0.9219098d0; pini=0.9503446d-2; H2ini=0.9465513d-8; Heini=0.9155226d-1; Hepini=0.5655353d-3
 Cini=0.1565848d-8; COini=0.2202631d-20; Cpini=0.1433520d-3
 dinit1=mH*Hini+mH*pini+mH2*H2ini+mHe*Heini+mHe*Hepini; dinit2=dinit1
 BBRV_cm(1)=Hini; BBRV_cm(2)=pini; BBRV_cm(3)=H2ini; BBRV_cm(4)=Heini
 BBRV_cm(5)=Hepini; BBRV_cm(6)=Cini; BBRV_cm(7)=COini; BBRV_cm(8)=Cpini
!10000 continue

IF(BCx1.eq.4) THEN; IF(IST.EQ.0)        LEFT = MPI_PROC_NULL; END IF
IF(BCx2.eq.4) THEN; IF(IST.EQ.NSPLTx-1) RIGT = MPI_PROC_NULL; END IF

Ncellx = Np1x + Np2x; Ncelly = Np1y + Np2y; Ncellz = Np1z + Np2z
gamma  = ( 5.d0*(Hini+pini+Heini+Hepini)+7.d0*H2ini )/( 3.d0*(Hini+pini+Heini+Hepini)+5.d0*H2ini )
gammi1 = gamma - 1.0d0; gammi2 = gamma - 2.0d0; gammi3 = gamma - 3.0d0
gampl1 = gamma + 1.0d0; gampl2 = gamma + 2.0d0; gampl3 = gamma + 3.0d0
pi     = 3.14159265358979323846d0


pmin = 1.829797d0 * 8.6336d0   !p/kb=1.d3
pmax = 1.d10  !604.5288d0 !p/kb =7.d4
rmin = 0.1949628d0
rmax = 1.d10  !4168.669d0*1.27d0

ndHmin  = rmin*0.91d0; ndpmin  = 1.d-20; ndH2min = 1.d-20; ndHemin = rmin*0.09d0
ndHepmin= 1.d-20; ndCpmin = 1.d-20; ndCmin = 1.d-20; ndCOmin = 1.d-20

!***** for constrained boundary *****!
BBRV_cm(1)=0.91d0/1.27d0 !H
BBRV_cm(2)=0.91d0*pini/Hini/1.27d0 !p
BBRV_cm(3)=0.91d0*H2ini/Hini/1.27d0 !H2
BBRV_cm(4)=0.09d0/1.27d0 !He
BBRV_cm(5)=0.09d0*Hepini/Heini/1.27d0 !Hep
BBRV_cm(8)=xc/1.27d0 !Cp
BBRV_cm(6)=xc*Cini/Cpini/1.27d0 !C
BBRV_cm(7)=xc*COini/Cpini/1.27d0 !CO

dBC = mH*BBRV_cm(1) + mH*BBRV_cm(2) + mH2*BBRV_cm(3) + mHe*BBRV_cm(4) + mHe*BBRV_cm(5)
BBRV(1,1,1) = dBC;     BBRV(1,2,1) = dBC;         BBRV(1,1,2) =  dBC;     BBRV(1,2,2) =  dBC
BBRV(2,1,1) = vinitx1; BBRV(2,2,1) = dBC*vinitx1; BBRV(2,1,2) = -vinitx1; BBRV(2,2,2) = -dBC*vinitx1
BBRV(3,1,1) = vinity1; BBRV(3,2,1) = dBC*vinity1; BBRV(3,1,2) =  vinity1; BBRV(3,2,2) =  dBC*vinity1
BBRV(4,1,1) = vinitz1; BBRV(4,2,1) = dBC*vinitz1; BBRV(4,1,2) =  vinitz1; BBRV(4,2,2) =  dBC*vinitz1
BBRV(5,1,1) = pinit1;  BBRV(5,2,1) = pinit1/gammi1 + 0.5d0*(dBC*vinitx1**2+binitx1**2+binity1**2+binitz1**2)
BBRV(5,1,2) = pinit1;  BBRV(5,2,2) = pinit1/gammi1 + 0.5d0*(dBC*vinitx1**2+binitx1**2+binity1**2+binitz1**2)
BBRV(6,1,1) = binitx1; BBRV(6,2,1) = binitx1; BBRV(6,1,2) = binitx2; BBRV(6,2,2) = binitx2
BBRV(7,1,1) = binity1; BBRV(7,2,1) = binity1; BBRV(7,1,2) = binity2; BBRV(7,2,2) = binity2
BBRV(8,1,1) = binitz1; BBRV(8,2,1) = binitz1; BBRV(8,1,2) = binitz2; BBRV(8,2,2) = binitz2

!***** x-direction shock tube test *****!

Ncellx = Ncellx/NSPLTx; Ncelly = Ncelly/NSPLTy; Ncellz = Ncellz/NSPLTz
dinit1 = mH*Hini + mH*pini + mH2*H2ini + mHe*Heini + mHe*Hepini

do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   i2 = IST*Ncellx+i
   !i2y = JST*Ncelly+j
   !i2z = KST*Ncellz+k
  if(i2.le.Np1x) then
    U(i,j,k,1) = dinit1
    U(i,j,k,2) = vinitx1
    U(i,j,k,3) = vinity1
    U(i,j,k,4) = vinitz1
    U(i,j,k,5) = pinit1
    U(i,j,k,6) = binitx1
    U(i,j,k,7) = binity1
    U(i,j,k,8) = binitz1
    ndH(i,j,k)   = Hini
    ndp(i,j,k)   = pini
    ndH2(i,j,k)  = H2ini
    ndHe(i,j,k)  = Heini
    ndHep(i,j,k) = Hepini
    ndC(i,j,k)   = Cini
    ndCO(i,j,k)  = COini
    ndCp(i,j,k)  = Cpini
    nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
    ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
    Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
    Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
  end if
  if(i2.gt.Np1x) then
    U(i,j,k,1) = dinit1
    U(i,j,k,2) = vinitx2
    U(i,j,k,3) = vinity1
    U(i,j,k,4) = vinitz1
    U(i,j,k,5) = pinit1
    U(i,j,k,6) = binitx1
    U(i,j,k,7) = binity2
    U(i,j,k,8) = binitz1
    ndH(i,j,k)   = Hini
    ndp(i,j,k)   = pini
    ndH2(i,j,k)  = H2ini
    ndHe(i,j,k)  = Heini
    ndHep(i,j,k) = Hepini
    ndC(i,j,k)   = Cini
    ndCO(i,j,k)  = COini
    ndCp(i,j,k)  = Cpini
    nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
    ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
    Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
    Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
  end if
end do; end do; end do
write(*,*) NRANK,'INIT'
ALLOCATE(dx_i(-1:Ncellx*NSPLTx+2)); ALLOCATE(dy_i(-1:Ncelly*NSPLTy+2)); ALLOCATE(dz_i(-1:Ncellz*NSPLTz+2))
ALLOCATE( x_i(-1:Ncellx*NSPLTx+2)); ALLOCATE( y_i(-1:Ncelly*NSPLTy+2)); ALLOCATE( z_i(-1:Ncellz*NSPLTz+2))

do i = -1, Np1x
  dx_i(i) = ql1x/dble(Np1x)
end do
do i = Np1x+1, Ncellx*NSPLTx+2
  dx_i(i) = ql2x/dble(Np2x)
end do
do j = -1, Ncelly*NSPLTy+2
  dy_i(j) = ql1y/dble(Np1y)
end do
do k = -1, Ncellz*NSPLTz+2
  dz_i(k) = ql1z/dble(Np1z)
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
  open(4,file='/work/maedarn/3DMHD/test/cdnt.DAT')
    write(4,400) ( 0.5d0 * ( x_i(i-1)+x_i(i) ), i=1, Ncellx*NSPLTx )
    write(4,400) ( 0.5d0 * ( y_i(j-1)+y_i(j) ), j=1, Ncelly*NSPLTy )
    write(4,400) ( 0.5d0 * ( z_i(k-1)+z_i(k) ), k=1, Ncellz*NSPLTz )
  close(4)
END IF

open(2,file='/work/maedarn/3DMHD/test/tsave.DAT')
  read(2,'(1p1d25.17)') amp
  read(2,'(i8)') nunit
  close(2)

  !********************sphere***********************
  !goto 6001
  DTF(:,:,:) = 0.0d0
  dinit1=1.0d0
  cenx=dble(Np1x)+0.5d0
  ceny=dble(Np1y)+0.5d0
  cenz=dble(Np1z)+0.5d0
  !rsph = ql1x-ql1x/5.0d0
  rsph2=int(dble(Np1x)*0.8d0)
  do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   i2 = IST*Ncellx+i
   i2y = JST*Ncelly+j
   i2z = KST*Ncellz+k
   cenx=dble(Np1x)+0.5d0
   ceny=dble(Np1y)+0.5d0
   cenz=dble(Np1z)+0.5d0
   rsph=dsqrt( (cenx-dble(i2))**2 + (ceny-dble(i2y))**2 + (cenz-dble(i2z))**2 )
   if(rsph.le.dble(rsph2)) then
      U(i,j,k,1) = dinit1
      U(i,j,k,2) = 0.0d0
      U(i,j,k,3) = 0.0d0
      U(i,j,k,4) = 0.0d0
      U(i,j,k,5) = pinit1
      U(i,j,k,6) = 0.0d0
      U(i,j,k,7) = 0.0d0
      U(i,j,k,8) = 0.0d0
      ndH(i,j,k)   = Hini
      ndp(i,j,k)   = pini
      ndH2(i,j,k)  = H2ini
      ndHe(i,j,k)  = Heini
      ndHep(i,j,k) = Hepini
      ndC(i,j,k)   = Cini
      ndCO(i,j,k)  = COini
      ndCp(i,j,k)  = Cpini
      nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
      ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
      Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
      Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
   else
      U(i,j,k,1) = 0.0d0
      U(i,j,k,2) = 0.0d0
      U(i,j,k,3) = 0.0d0
      U(i,j,k,4) = 0.0d0
      U(i,j,k,5) = 0.0d0
      U(i,j,k,6) = 0.0d0
      U(i,j,k,7) = 0.0d0
      U(i,j,k,8) = 0.0d0
      ndH(i,j,k)   = 0.0d0
      ndp(i,j,k)   = 0.0d0
      ndH2(i,j,k)  = 0.0d0
      ndHe(i,j,k)  = 0.0d0
      ndHep(i,j,k) = 0.0d0
      ndC(i,j,k)   = 0.0d0
      ndCO(i,j,k)  = 0.0d0
      ndCp(i,j,k)  = 0.0d0
      nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
      ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
      Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
      Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
   end if
end do
end do
end do
 !6001 continue
 !********************sphere***********************

  !********purtube yz plane***********!
  goto 1333
  ALLOCATE (plane(-1:Ncelly*NSPLTy+2,-1:Ncellz*NSPLTz+2))
  open(unit=28,file='/work/maedarn/3DMHD/test/delta2.dat',FORM='UNFORMATTED')
  do c=-1,Ncellz*NSPLTz+2
     do b=-1,Ncelly*NSPLTy+2
        read(28) plane(b,c)
     end do
  end do
  close(28)


do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   i2 = IST*Ncellx+i
   i3 = JST*Ncelly+j
   i4 = KST*Ncellz+k
  if(x_i(i2).le.plane(i3,i4)) then
    U(i,j,k,1) = dinit1
    U(i,j,k,2) = vinitx1
    U(i,j,k,3) = vinity1
    U(i,j,k,4) = vinitz1
    U(i,j,k,5) = pinit1
    U(i,j,k,6) = binitx1
    U(i,j,k,7) = binity1
    U(i,j,k,8) = binitz1
    ndH(i,j,k)   = Hini
    ndp(i,j,k)   = pini
    ndH2(i,j,k)  = H2ini
    ndHe(i,j,k)  = Heini
    ndHep(i,j,k) = Hepini
    ndC(i,j,k)   = Cini
    ndCO(i,j,k)  = COini
    ndCp(i,j,k)  = Cpini
    nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
    ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
    Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
    Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
  end if
  if(x_i(i2).gt.plane(i3,i4)) then
    U(i,j,k,1) = dinit1
    U(i,j,k,2) = vinitx2
    U(i,j,k,3) = vinity1
    U(i,j,k,4) = vinitz1
    U(i,j,k,5) = pinit1
    U(i,j,k,6) = binitx1
    U(i,j,k,7) = binity2
    U(i,j,k,8) = binitz1
    ndH(i,j,k)   = Hini
    ndp(i,j,k)   = pini
    ndH2(i,j,k)  = H2ini
    ndHe(i,j,k)  = Heini
    ndHep(i,j,k) = Hepini
    ndC(i,j,k)   = Cini
    ndCO(i,j,k)  = COini
    ndCp(i,j,k)  = Cpini
    nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
    ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
    Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
    Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
  end if
end do; end do; end do


1333 continue
  !********purtube yz plane***********!


!/work/maedarn/3DMHD/test/
!***** Alfven wave propagation *****!
goto 111
do k = 1, Ncellz+1; do j = 1, Ncelly+1; do i = 1, Ncellx+1
  xpi = 0.5d0*( x(i)+x(i-1) ); amp = 1.d-3
  U(i,j,k,3) =  amp*dcos(2.d0*pi*xpi)
  U(i,j,k,4) =  amp*dcos(2.d0*pi*xpi)
  U(i,j,k,7) = -amp*dsqrt(dinit1)*dcos(2.d0*pi*xpi)
  U(i,j,k,8) =  amp*dsqrt(dinit1)*dcos(2.d0*pi*xpi)
end do; end do; end do
111 continue
!***** Blast wave *****!
goto 112
do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
  xpi = 0.5d0*( x(i)+x(i-1) ); ypi = 0.5d0*( y(j)+y(j-1) )
  amp = dsqrt( (xpi-0.5d0)**2 + (ypi-0.5d0)**2 )
  if(amp.lt.0.125d0) U(i,j,k,5) =  1.d2
end do; end do; end do
112 continue

!**** read inhomogeneous density field ****!
  !DTF(:,:,:) = dinit1
  goto 119
  do MRANK = 0, NPE-1
    IS = mod(MRANK,NSPLTx); KS = MRANK/(NSPLTx*NSPLTy); JS = MRANK/NSPLTx-NSPLTy*KS
    if((JS.eq.JST).and.(KS.eq.KST)) then
      WRITE(NPENUM,'(I3.3)') MRANK
      open(unit=8,file='/work/maedarn/3DMHD/test/DTF/D'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
      do k = 1, Ncellz
      do j = 1, Ncelly
        read(8) (DTF(i,j,k),i=Ncellx*IS+1,Ncellx*IS+Ncellx)
      end do
      end do
      close(8)
    end if
  end do

  CALL MPI_TYPE_VECTOR(Ncellz+4,2*(Ncellx*NSPLTx+4),(Ncellx*NSPLTx+4)*(Ncelly+4),MPI_REAL4,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  CALL MPI_SENDRECV(DTF(-1,Ncelly-1,-1),1,VECU,TOP ,1, &
                    DTF(-1,      -1,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(DTF(-1,1       ,-1),1,VECU,BOTM,1, &
                    DTF(-1,Ncelly+1,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_TYPE_FREE(VECU,IERR)
  CALL MPI_TYPE_VECTOR(1,2*(Ncellx*NSPLTx+4)*(Ncelly+4),2*(Ncellx*NSPLTx+4)*(Ncelly+4),MPI_REAL4,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  CALL MPI_SENDRECV(DTF(-1,-1,Ncellz-1),1,VECU,UP  ,1, &
                    DTF(-1,-1,      -1),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(DTF(-1,-1,1       ),1,VECU,DOWN,1, &
                    DTF(-1,-1,Ncellz+1),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_TYPE_FREE(VECU,IERR)

  if(nunit.ne.1) goto 119
  write(*,*) '119',NRANK
  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
    ix = Ncellx*IST+i
    U(i,j,k,1)   = dble(DTF(ix,j,k))
    U(i,j,k,2)   = -vinitx1*dtanh(0.5d0*(x(i)-ql1x))
    ndH(i,j,k)   = U(i,j,k,1)*BBRV_cm(1)
    ndp(i,j,k)   = U(i,j,k,1)*BBRV_cm(2)
    ndH2(i,j,k)  = U(i,j,k,1)*BBRV_cm(3)
    ndHe(i,j,k)  = U(i,j,k,1)*BBRV_cm(4)
    ndHep(i,j,k) = U(i,j,k,1)*BBRV_cm(5)
    ndC(i,j,k)   = U(i,j,k,1)*BBRV_cm(6)
    ndCO(i,j,k)  = U(i,j,k,1)*BBRV_cm(7)
    ndCp(i,j,k)  = U(i,j,k,1)*BBRV_cm(8)
  end do;end do;end do
  119 continue
!--------------------------------

DEALLOCATE(dx_i); DEALLOCATE(dy_i); DEALLOCATE(dz_i); DEALLOCATE(x_i); DEALLOCATE(y_i); DEALLOCATE(z_i)

!***** Read Initial Conditions *****!
if(nunit.eq.1) goto 120
  WRITE(NPENUM,'(I3.3)') NRANK
  open(unit=8,file='/work/maedarn/3DMHD/test/000'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
  do k = 1, Ncellz+1
  do j = 1, Ncelly+1
    read(8) (U(i,j,k,1),U(i,j,k,2),U(i,j,k,3),U(i,j,k,4),U(i,j,k,5),U(i,j,k,6),U(i,j,k,7),U(i,j,k,8), &
             ndH(i,j,k),ndp(i,j,k),ndH2(i,j,k),ndHe(i,j,k),          &
             ndHep(i,j,k),ndC(i,j,k),ndCO(i,j,k),ndCp(i,j,k),Phi(i,j,k),i=1,Ncellx+1)
  end do
  end do
  close(8)
120  continue

IF(BCx1.eq.4) THEN; IF(IST.EQ.0)        LEFT = MPI_PROC_NULL; END IF
IF(BCx2.eq.4) THEN; IF(IST.EQ.NSPLTx-1) RIGT = MPI_PROC_NULL; END IF

!call CC(4,0.d0)

do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
  nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
  ndtot(i,j,k) = ndp(i,j,k)+ndH(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
  Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; NCO(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
  Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; NCO(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
end do; end do; end do

!if(ifrad.eq.2) then; do l=1,20; call SHIELD(); end do; end if
write(*,*) '==ok1=='
if(ifgrv.eq.2) then
  N_MPI(20)=1; N_MPI(1)=1; iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)
  Lbox=ql1x+ql2x; call GRAVTY(0.d0,1); call GRAVTY(0.d0,2)
end if
write(*,*) '==ok2=='
END SUBROUTINE INITIA


SUBROUTINE ran0(ran,idum)
INTEGER idum,IA,IM,IQ,IR,MASK
DOUBLE PRECISION ran,AM
PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,MASK=123459876)
INTEGER k
!idum=ieor(idum,MASK)
k=idum/IQ
idum=IA*(idum-k*IQ)-IR*k
if(idum.lt.0) idum=idum+IM
ran=AM*idum
!idum=ieor(idum,MASK)
END SUBROUTINE ran0


!=====================================================================*
!                  Integration of The Evolution                       *
!=====================================================================*

SUBROUTINE EVOLVE
USE comvar
USE mpivar
USE chmvar
INCLUDE 'mpif.h'

double precision  :: t(1000),dt, stt, tLMT, dt_mpi(0:1024), dt_gat(0:1024), time_CPU(3)
double precision  :: tsave,dtsave,tsave2D,dtsave2D
integer :: nunit, st, st_mpi(0:1024), st_gat(0:2047), Time_signal
character*7 stb(3)
character*3 fnunit,fnpe


!===========time===========
goto 2200
open(2,file='/work/maedarn/3DMHD/test/tsave.DAT')
  read(2,*) time
  read(2,*) nunit
close(2)
open(2,file='/work/maedarn/3DMHD/test/tsave2D.DAT')
  read(2,*) nunit2D
close(2)
open(3,file='/work/maedarn/3DMHD/test/time.DAT')
do i = 1, nunit
  read(3,'(1p1d25.17)') t(i)
end do
close(3)
!IF(NRANK.EQ.0) THEN
!  open(2,file='/work/maedarn/3DMHD/test/test.DAT')
!END IF

write(fnunit,'(I3.3)') nunit;  write(fnpe,'(I3.3)') NRANK
open(5,file='/work/maedarn/3DMHD/test/info'//fnunit//'.DAT')
!open(5,file='/work/maedarn/3DMHD/test/info'//fnunit//fnpe//'.DAT')

st    = 1
ifEVO = 1
dt    = 0.d0
dtsave = 0.25d0
dtsave2D = 0.025d0
itime  = 1 + int( (time + 1.d-8)/dtsave )

stb(1)='FastSpd'
stb(2)='Conduct'
time_CPU(1) = 0.d0
time_CPU(2) = 0.d0
time_CPU(3) = 0.d0
Time_signal = 0
2200 continue
!===========time===========

!do in10 = 1, maxstp

!  time_CPU(1) = MPI_WTIME()
!  tsave = dtsave * dble(itime)
!  if(time.ge.tfinal) goto 9000
!  if(time.ge.tsave ) goto 7777



!call SAVEU(nunit,dt,stb,st,t,0)




do in20 = 1, nitera

!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),Bcc(1,1,1,2),U(1,1,1,7),'point'
!    tsave2D = dtsave2D * nunit2D
!    if(time.ge.tsave2D) call SAVEU2D(nunit2D)
!    if(time.ge.tfinal) goto 9000
!    if(time.ge.tsave ) goto 7777
!***** Determine time-step dt *****
!    dt_mpi(NRANK) = tfinal
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point1'
!    call Couran(tLMT)
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),tLMT,'point1'
!    dt_mpi(NRANK) = dmin1( dt_mpi(NRANK), CFL * tLMT )
!    st_mpi(NRANK) = 1
!    stt= dt_mpi(NRANK)

!    call Stblty(tLMT)
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),tLMT,'point2'
!    dt_mpi(NRANK) = dmin1( dt_mpi(NRANK), tLMT    )
!    if(dt_mpi(NRANK).lt.stt) st_mpi(NRANK) = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for MPI
!    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!    CALL MPI_GATHER(dt_mpi(NRANK),1,MPI_REAL8,   &
!                    dt_gat       ,1,MPI_REAL8,   &
!                    0            ,MPI_COMM_WORLD,IERR)
!    CALL MPI_GATHER(st_mpi(NRANK),1,MPI_INTEGER, &
!                    st_gat       ,1,MPI_INTEGER, &
!                    0            ,MPI_COMM_WORLD,IERR)
!    IF(NRANK.EQ.0)  THEN
!      dt  = tfinal
!      dtt = tfinal
!      do i_t = 0, NPE-1
!        dt  = dmin1( dt, dt_gat(i_t) )
!        if(dt.lt.dtt) st = st_gat(i_t)
!        dtt = dt
!      end do
!END IF
!    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!    CALL MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if((mod(in20,10).eq.1).and.(NRANK.eq.0)) write(*,*) in20,time,dt
!    if(NRANK.eq.0) write(*,*) in20,time,dt
!    if(time+dt.gt.tfinal) dt = tfinal - time
!    if(time+dt.gt.tsave ) dt = tsave  - time

!if(NRANK==40) write(*,*) NRANK,in20,dt,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point3'
!***** Source parts 1*****
    !if(ifgrv.eq.2) then; call GRAVTY(dt,3); end if
!   call SOURCE(0.5d0*dt)
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point4'
!***** Godunov parts *****
       goto 2201
    if(ifEVO.eq.1) then
      iwx=1; iwy=0; iwz=0; call MHD(x,dx,dt); iwx=0; iwy=1; iwz=0; call MHD(y,dy,dt); iwx=0; iwy=0; iwz=1; call MHD(z,dz,dt)
      ifEVO = 2; goto 1000
    end if
    if(ifEVO.eq.2) then
      iwx=0; iwy=1; iwz=0; call MHD(y,dy,dt); iwx=0; iwy=0; iwz=1; call MHD(z,dz,dt); iwx=1; iwy=0; iwz=0; call MHD(x,dx,dt)
      ifEVO = 3; goto 1000
    end if
    if(ifEVO.eq.3) then
      iwx=0; iwy=0; iwz=1; call MHD(z,dz,dt); iwx=1; iwy=0; iwz=0; call MHD(x,dx,dt); iwx=0; iwy=1; iwz=0; call MHD(y,dy,dt)
      ifEVO = 4; goto 1000
    end if
    if(ifEVO.eq.4) then
      iwx=1; iwy=0; iwz=0; call MHD(x,dx,dt); iwx=0; iwy=0; iwz=1; call MHD(z,dz,dt); iwx=0; iwy=1; iwz=0; call MHD(y,dy,dt)
      ifEVO = 5; goto 1000
    end if
    if(ifEVO.eq.5) then
      iwx=0; iwy=1; iwz=0; call MHD(y,dy,dt); iwx=1; iwy=0; iwz=0; call MHD(x,dx,dt); iwx=0; iwy=0; iwz=1; call MHD(z,dz,dt)
      ifEVO = 6; goto 1000
    end if
    if(ifEVO.eq.6) then
      iwx=0; iwy=0; iwz=1; call MHD(z,dz,dt); iwx=0; iwy=1; iwz=0; call MHD(y,dy,dt); iwx=1; iwy=0; iwz=0; call MHD(x,dx,dt)
      ifEVO = 1; goto 1000
    end if
1000 continue
    DEALLOCATE(Bcc)
if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point5'
!***** CT part *****
    ALLOCATE(Vfc(-1:ndx,-1:ndy,-1:ndz,3))
    call CC(1,dt)
    ALLOCATE(dnc(-1:ndx,-1:ndy,-1:ndz)); ALLOCATE(EMF(-1:ndx,-1:ndy,-1:ndz,3))
    iwx=1; iwy=0; iwz=0; call CC(2,dt); call CCT(dx,dy,dt) !calculate Ez
    iwx=0; iwy=1; iwz=0; call CC(2,dt); call CCT(dy,dz,dt) !calculate Ex
    iwx=0; iwy=0; iwz=1; call CC(2,dt); call CCT(dz,dx,dt) !calculate Ey
    DEALLOCATE(Vfc); DEALLOCATE(dnc)

    call CC(3,dt)
    DEALLOCATE(EMF)
    call CC(4,dt)
    if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point6'
2201 continue


!***** Source parts 2*****
!    call SOURCE(0.5d0*dt)
    if(ifgrv.eq.2) then; call GRAVTY(0.0d0,2)!; call GRAVTY(dt,3);
    end if
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point7'
!    call DISSIP()
!    time = time + dt
!  end do
!  itime = itime - 1
!  7777   continue
!  itime = itime + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for MPI
!  IF(NRANK.EQ.0) THEN
!    time_CPU(2) = MPI_WTIME()
!    time_CPU(2) = ( time_CPU(2)-time_CPU(1) )/3.6d3
!    time_CPU(3) = time_CPU(3)+time_CPU(2)
!    IF(time_CPU(3)+time_CPU(2).GT.11.7d0) Time_signal=1
!  END IF
!  CALL MPI_BCAST(Time_signal,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
!  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!  IF(Time_signal.EQ.1) GOTO 9000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do

!9000 continue
!IF(NRANK.EQ.0) write(*,*) 'MPI time1 = ',MPI_WTIME()
!call SAVEU(nunit,dt,stb,st,t,1)

END SUBROUTINE EVOLVE

SUBROUTINE GRAVTY(dt,mode)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION  :: dt,dxi
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

if(mode.eq.1) then
  call pinter(Nmem1,Nmem2,Ncellx,Ncelly,Ncellz)
end if

if(mode.eq.2) then
  N_MPI(20)=1; N_MPI(1)=1
  iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)
  call PB()
  call mglin(Nmem1,Nmem2,2,5,5)
  DEALLOCATE(bphi1,bphi2)

  N_ol = 2
                      !count, blocklength, stride
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
  !*****  BC for the leftsides of domains  *****
  CALL MPI_SENDRECV(Phi(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
                    Phi(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !*****  BC for the rightsides of domains *****
  CALL MPI_SENDRECV(Phi(1            ,-1,-1),1,VECU,LEFT,1, &
                    Phi(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_TYPE_FREE(VECU,IERR)
  LEFT = LEFTt; RIGT = RIGTt

  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
  !*****  BC for the downsides of domains  ****
  CALL MPI_SENDRECV(Phi(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
                    Phi(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !*****  BC for the upsides of domains  ****
  CALL MPI_SENDRECV(Phi(-1,1            ,-1),1,VECU,BOTM,1, &
                    Phi(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt

  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
  !*****  BC for the downsides of domains  ****
  CALL MPI_SENDRECV(Phi(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                    Phi(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !*****  BC for the upsides of domains  ****
  CALL MPI_SENDRECV(Phi(-1,-1,1            ),1,VECU,DOWN,1, &
                    Phi(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt

  if(IST.eq.0       ) then; do k=1,Ncellz; do j=1,Ncelly
    Phi(0       ,j,k) = Phi(1     ,j,k); Phi(-1       ,j,k) = Phi(1     ,j,k) !grad=0
  end do; end do; end if
  if(IST.eq.NSPLTx-1) then; do k=1,Ncellz; do j=1,Ncelly
    Phi(Ncellx+1,j,k) = Phi(Ncellx,j,k); Phi(Ncellx+2,j,k) = Phi(Ncellx,j,k)
  end do; end do; end if
end if

!if(mode.eq.3) then !acceraration because of gravity
!  dxi = 1.d0/(12.d0*dx(0))
!  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
!    U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
!    U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
!    U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0
!  end do;end do;end do
!end if

END SUBROUTINE GRAVTY

SUBROUTINE pinter(Need1,Need2,Ncellx,Ncelly,Ncellz)
USE slfgrv
USE mpivar

!***  finest grid pointer : point(NGL) ***
!*** coasest grid pointer : point(0)   ***

NGL = min0(Ncellx*NSPLTx,Ncelly*NSPLTy,Ncellz*NSPLTz)
NGL = int(dlog(dble(NGL))/dlog(2.d0)+1.d-3)


!*** Unsplit Pointer ***!
NGcr = max0(NSPLTx,NSPLTy,NSPLTz)
NGcr = int(dlog(dble(NGcr))/dlog(2.d0)+1.d-3) + 1

point1(1) = 1
nl = 1
nx=3; ny=3; nz=3
2 continue
point1(nl+1)=point1(nl)+(nx)*(ny)*(nz)
write(*,*) NRANK,point1(nl+1),nx,nl ,'point1'
nx=nx*2-1; ny=ny*2-1; nz=nz*2-1
nl = nl+1
if(nl.ne.NGcr+1) goto 2

Need1 = point1(NGcr+1)

!*** MPI Split Pointer ***!
!nl = NGcr-1
nl=NGcr
point2(nl) = 1
nx=(2**NGcr)/NSPLTx+2; ny=(2**NGcr)/NSPLTy+2; nz=(2**NGcr)/NSPLTz+2
3 continue
point2(nl+1)=point2(nl)+(nx)*(ny)*(nz)
write(*,*) NRANK,point2(nl+1),nx,nl,'point2'
nx=nx*2-2; ny=ny*2-2; nz=nz*2-2
nl = nl+1
if(nl.ne.NGL+1) goto 3

Need2 = point2(NGL+1)

if(NRANK.eq.0) write(*,*) 'NGL=',NGL
if(NRANK.eq.0) write(*,*) 'NGcr=',NGcr
if(NRANK.eq.0) write(*,*) 'need1=',Need1
if(NRANK.eq.0) write(*,*) 'need2=',Need2

END SUBROUTINE pinter

!***********************************
! set BC in interp (except in addint) and the first slvsml when fixed boundary except phi=0
!
SUBROUTINE mglin(Need1,Need2,ncycle,NPRE,NPOST)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION :: cphi1(Need1),crho1(Need1),cres1(Need1),crhs1(Need1) !Unsplit
DOUBLE PRECISION :: cphi2(Need2),crho2(Need2),cres2(Need2),crhs2(Need2) !MPI Split
DOUBLE PRECISION :: tMPI(Need1,0:NPE-1)
character(3) nnn
character(3) name1
double precision :: shd1=1.0d4
!***********initialaze**********
!cphi1(:)=0.0d0
!crho1(:)=0.0d0
!cres1(:)=0.0d0
!crhs1(:)=0.0d0
!cphi2(:)=0.0d0
!crho2(:)=0.0d0
!cres2(:)=0.0d0
!crhs2(:)=0.0d0
!***********initialaze**********

write(name1,'(I3.3)') NRANK
open(250+NRANK,file='GINI'//name1//'.dat')
do k=0,Ncellz+1; kk=(Ncellx+2)*(Ncelly+2)*k+point2(NGL)
do j=0,Ncelly+1; jj=(Ncellx+2)*j+kk
do i=0,Ncellx+1; nc = i+jj
   write(250+NRANK,*) nc,jj,kk,i,j,k,  U(i,j,k,1)*G4pi
end do
end do
end do
close(250+NRANK)
!Pre-BC for rho is necessary
!write(*,*) 'mglin' , NRANK

write(*,*) point2(NGL)+(Ncellz+2)*(Ncelly+2)*(Ncellx+2),Need2,'neeeed'

do k=0,Ncellz+1; kk=(Ncellx+2)*(Ncelly+2)*k+point2(NGL)
do j=0,Ncelly+1; jj=(Ncellx+2)*j+kk
do i=0,Ncellx+1; nc = i+jj
   crhs2(nc) = U(i,j,k,1)*G4pi
end do;end do;end do

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'mglin2' , NRANK

ncx=(Ncellx+1)/2+1; ncy=(Ncelly+1)/2+1; ncz=(Ncellz+1)/2+1
ngrid=NGL-1
write(*,*) NRANK,ngrid,ncx,point2(ngrid),point2(ngrid+1),point2(ngrid+1)-point2(ngrid),'=======ppp-======'
call rstrctMPI(crho2(point2(ngrid)),crhs2(point2(ngrid+1)),ncx,ncy,ncz,0)

do while(ngrid.ne.NGcr) !真で繰り返す
  ncx=ncx/2+1; ncy=ncy/2+1; ncz=ncz/2+1
  ngrid=ngrid-1
  call rstrctMPI(crho2(point2(ngrid)),crho2(point2(ngrid+1)),ncx,ncy,ncz,0)
  write(*,*) ncx, ngrid,NRANK,crho2(point2(ngrid)),crho2(point2(ngrid)+1),crho2(point2(ngrid+1)-1),point2(ngrid),point2(ngrid+1),&
  point2(ngrid+1)-point2(ngrid),'gggrid'
end do
ncxcr=ncx; ncycr=ncy; nczcr=ncz

ncx=2**NGcr+1;ncy=ncx;ncz=ncx
write(*,*) ncxcr,ncx, ngrid,NGcr,NRANK,point2(NGcr+1)-point2(NGcr),point2(NGcr+1),point2(NGcr),crho2(point2(NGcr+1)), &
     point1(NGcr),point1(NGcr+1)-point1(NGcr),crho2(point2(NGcr)), 'precll'

!call saveu1(crho1(point1(NGcr)),crho2(point2(NGcr)),ncx,ncy,ncz,ncxcr,ncycr,nczcr)
call collect(crho1(point1(NGcr)),crho2(point2(NGcr)),ncx,ncy,ncz,ncxcr,ncycr,nczcr,1)

do while(ngrid.ne.1)
  ncx=ncx/2+1; ncy=ncy/2+1; ncz=ncz/2+1
  ngrid=ngrid-1
  call rstrct(crho1(point1(ngrid)),crho1(point1(ngrid+1)),ncx,ncy,ncz,0,NRANK)
end do

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

call slvsmlb(cphi1(point1(1)),crho1(point1(1))) !BC set is necessary
!write(*,*) 'mglin3' , NRANK






!Here nc=3
ngrid = NGL
do j=2,ngrid
   write(*,*) 'mglin0' , NRANK , j
   IF(j.le.NGcr) THEN !*** generate candidate sol. from j-1 to j (upward) ***

      if(cphi1(point1(j-1)+1) > shd1) then
         write(*,*) j,cphi1(point1(j-1)+1),NRANK,'==============1==============='
      end if

     ncx=ncx*2-1; ncy=ncy*2-1; ncz=ncz*2-1
     !********:fordebug:************
     CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
     !********:fordebug:************
    call interp(cphi1(point1(j)),cphi1(point1(j-1)),ncx,ncy,ncz,pointb1(j),1,NRANK)  !BC set is necessary
    call copy(crhs1(point1(j)),crho1(point1(j)),ncx,ncy,ncz)

    if(cphi1(point1(j)+1) > shd1) then
       write(*,*) j,cphi1(point1(j)+1),NRANK,'==============2==============='
    end if
!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

if(j.eq.NGcr) then

   !********:fordebug:************
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   !********:fordebug:************
    if(cphi1(point1(NGcr)+1) > shd1) then
       write(*,*) j,cphi1(point1(j)+1),NRANK,'==============3==============='
    end if

      ncx=ncxcr; ncy=ncycr; ncz=nczcr
      call divide( cphi2(point2(NGcr)),cphi1(point1(NGcr)),ncx,ncy,ncz,2**NGcr+1,2**NGcr+1,2**NGcr+1 )
      call copyMPI(crhs2(point2(NGcr)),crho2(point2(NGcr)),ncx,ncy,ncz)

      if(cphi2(point1(NGcr)+1) > shd1) then
         write(*,*) j,cphi1(point1(j)+1),NRANK,'==============4==============='
      end if
    end if
 ELSE

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

if(cphi2(point2(j-1)+1) > shd1) then
   write(*,*) j,cphi2(point2(j-1)+1),NRANK,'==============5==============='
end if

    ncx=ncx*2-1; ncy=ncy*2-1; ncz=ncz*2-1
    call interpMPI(cphi2(point2(j)),cphi2(point2(j-1)),ncx,ncy,ncz,pointb2(j),1)  !BC set is necessary
    if(j.ne.ngrid) call copyMPI(crhs2(point2(j)),crho2(point2(j)),ncx,ncy,ncz)

    if(cphi2(point2(j)+1) > shd1) then
       write(*,*) j,cphi2(point2(j)+1),NRANK,'==============6==============='
    end if
  END IF
  write(*,*) 'mglin1' , NRANK , j

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

  do jcycle=1,ncycle !V-cycle
    write(*,*) 'mglin2' , NRANK , j , jcycle
    nfx=ncx; nfy=ncy; nfz=ncz
    do jj=j,2,-1          !*** DOWNWARD *****************************
                          !    phi + rhs --> res -->      (level N  )
                          !                          rhs  (level N-1)

       write(*,*) 'mglinpost1' , NRANK , j , jcycle ,jj
       IF(jj.lt.NGcr) THEN !*** generate residual from jj to jj-1 ****

         ! if(cphi2(point2(j-1)+1) > shd1) then
         !    write(*,*) j,cphi2(point2(j-1)+1),NRANK,'==============5==============='
         ! end if

        do jpre=1,NPRE
           mode=2; if((jj.ne.j).and.(jpre.eq.1)) mode=1
           !********:fordebug:************
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
           !********:fordebug:************

           if(cphi1(point1(jj)+1) > shd1) then
              write(*,*) j,jpre,jcycle,cphi1(point1(jj)+1),crhs1(point1(jj)+1),NRANK,'==============7==============='
           end if

           call relax(cphi1(point1(jj)),crhs1(point1(jj)),nfx,nfy,nfz,mode,NRANK)

           if(cphi1(point1(jj)+1) > shd1) then
              write(*,*) j,jpre,jcycle,cphi1(point1(jj)+1),crhs1(point1(jj)+1),NRANK,'==============8==============='
           end if

       end do
       !********:fordebug:************
       CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
       !********:fordebug:************
        call resid(cres1(point1(jj)),cphi1(point1(jj)),crhs1(point1(jj)),nfx,nfy,nfz)
        nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1

        if(cres1(point1(jj)+1) > shd1) then
           write(*,*) j,jcycle,cphi1(point1(jj)+1),cres1(point1(jj)+1) ,NRANK,'==============9==============='
        end if

        call rstrct(crhs1(point1(jj-1)),cres1(point1(jj)),nfx,nfy,nfz,1,NRANK)  !fill0 at BC below this subroutine is necessary
        call  fill0(cphi1(point1(jj-1)),nfx,nfy,nfz)
      ELSE
        NPRE1 = NPRE; if(j.ge.NGL) NPRE1 = 2
        do jpre=1,NPRE1
           mode=2; if((jj.ne.j).and.(jpre.eq.1)) mode=1
           !********:fordebug:************
           CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
           !********:fordebug:************

           if(cphi2(point2(jj)+1) > shd1) then
              write(*,*) j,jpre,jcycle,cphi2(point2(jj)+1),crhs2(point2(jj)+1),NRANK,'==============10==============='
           end if

           call relaxMPI(cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz,mode)

           if(cphi2(point2(jj)+1) > shd1) then
              write(*,*) j,jpre,jcycle,cphi2(point2(jj)+1),crhs2(point2(jj)+1),NRANK,'==============11==============='
           end if
       end do
       !write(*,*) 'mglinpost2' , NRANK , j , jcycle ,jj
       call residMPI(cres2(point2(jj)),cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz)
       !write(*,*) 'mglinpost3' , NRANK , j , jcycle ,jj

       if(cres2(point2(jj)+1) > shd1) then
           write(*,*) j,jcycle,cphi2(point1(jj)+1),cres2(point1(jj)+1) ,NRANK,'==============12==============='
        end if
       if(jj.eq.NGcr) then

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

       !   write(*,*) 'mglinpost4' , NRANK , j , jcycle ,jj
          nfx=2**NGcr+1;nfy=nfx;nfz=nfx
          call collect( cphi1(point1(NGcr)),cphi2(point2(NGcr)),nfx,nfy,nfz,ncxcr,ncycr,nczcr,1 ) !necessary at upward loop below
          call collect( cres1(point1(NGcr)),cres2(point2(NGcr)),nfx,nfy,nfz,ncxcr,ncycr,nczcr,1 )

          if(cres1(point1(NGcr)+1) > shd1 .or. cphi1(point1(NGcr)+1) > shd1) then
             write(*,*) j,jcycle,cphi1(point1(NGcr)+1),cres1(point1(NGcr)+1),cphi2(point2(NGcr)+1) ,NRANK,'==============13==============='
          end if
          nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
       !   write(*,*) 'mglinpost5' , NRANK , j , jcycle ,jj
          call rstrct(crhs1(point1(jj-1)),cres1(point1(jj)),nfx,nfy,nfz,1,NRANK)  !fill0 at BC below this subroutine is necessary
          call  fill0(cphi1(point1(jj-1)),nfx,nfy,nfz)
       !   write(*,*) 'mglinpost6' , NRANK , j , jcycle ,jj
        else
          nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
          call rstrctMPI(crhs2(point2(jj-1)),cres2(point2(jj)),nfx,nfy,nfz,1)  !fill0 at BC below this subroutine is necessary
          call  fill0MPI(cphi2(point2(jj-1)),nfx,nfy,nfz)
        end if
      END IF
    end do

    call slvsml(cphi1(point1(1)),crhs1(point1(1)))  !BC set is unnecessary

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

    nfx=3; nfy=3; nfz=3
    do jj=2,j             !*** UPWARD **********************************
                          !    phi --> phi +rhs --> phi      (level N  )
                          !  + phi                           (level N-1)
      IF(jj.le.NGcr) THEN !*** generate new solution from jj-1 to jj ***
         nfx=2*nfx-1; nfy=2*nfy-1; nfz=2*nfz-1


         if(cphi1(point1(jj)+1) > shd1) then
              write(*,*) j,jj,jcycle,cphi1(point1(jj-1)+1),cres1(point1(jj)+1),NRANK,'==============14==============='
           end if

         call addint(cphi1(point1(jj)),cphi1(point1(jj-1)),cres1(point1(jj)),nfx,nfy,nfz,NRANK)  !BC set is unnecessary
        !********:fordebug:************
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        !********:fordebug:************

        if(cphi1(point1(jj)+1) > shd1) then
              write(*,*) j,jj,jcycle,cphi1(point1(jj)+1),cres1(point1(jj)+1),NRANK,'==============15==============='
           end if

        if(jj.eq.NGcr) then
          nfx=ncxcr; nfy=ncycr; nfz=nczcr
          call divide( cphi2(point2(NGcr)),cphi1(point1(NGcr)),nfx,nfy,nfz,2**NGcr+1,2**NGcr+1,2**NGcr+1 ) !境界の値は0?!!!!!!!!!!!!!!!!!!!!!!!!!
          !********:fordebug:************
          CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          !********:fordebug:************

          if(cphi2(point2(jj)+1) > shd1) then
              write(*,*) j,jj,jcycle,cphi2(point2(jj)+1),cres2(point2(jj)+1),NRANK,'==============16==============='
           end if

           do jpost=1,NPOST

              if(cphi2(point2(jj)+1) > shd1) then
                 write(*,*) j,jj,jpost,jcycle,cphi2(point2(jj)+1),crhs2(point2(jj)+1),NRANK,'==============17==============='
              end if

              call relaxMPI(cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz,2)

              if(cphi2(point2(jj)+1) > shd1) then
                 write(*,*) j,jj,jpost,jcycle,cphi2(point2(jj)+1),crhs2(point2(jj)+1),NRANK,'==============18==============='
              end if
          end do
       else
          !********:fordebug:************
          CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
          !********:fordebug:************
          do jpost=1,NPOST
             if(cphi1(point1(jj)+1) > shd1) then
                 write(*,*) j,jj,jpost,jcycle,cphi1(point1(jj)+1),crhs1(point1(jj)+1),NRANK,'==============19==============='
              end if
              call relax(cphi1(point1(jj)),crhs1(point1(jj)),nfx,nfy,nfz,2,NRANK)
              if(cphi1(point1(jj)+1) > shd1) then
                 write(*,*) j,jj,jpost,jcycle,cphi1(point1(jj)+1),crhs1(point1(jj)+1),NRANK,'==============20==============='
              end if
          end do
        end if
     ELSE
        !********:fordebug:************
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        !********:fordebug:************
        if(cphi2(point2(jj-1)+1) > shd1) then
           write(*,*) j,jj,jcycle,cphi2(point2(jj-1)+1),cres2(point2(jj)+1),NRANK,'==============17==============='
        end if
        nfx=2*nfx-1; nfy=2*nfy-1; nfz=2*nfz-1
        call addintMPI(cphi2(point2(jj)),cphi2(point2(jj-1)),cres2(point2(jj)),nfx,nfy,nfz)  !BC set is unnecessary
        NPOST1 = NPOST; if(j.ge.NGL) NPOST1 = 2
        !********:fordebug:************
        CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        !********:fordebug:************
        do jpost=1,NPOST1
           if(cphi2(point2(jj)+1) > shd1) then
              write(*,*) j,jj,jpost,jcycle,cphi2(point2(jj)+1),crhs2(point2(jj)+1),NRANK,'==============22==============='
           end if
           call relaxMPI(cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz,2)
           if(cphi2(point2(jj)+1) > shd1) then
              write(*,*) j,jj,jpost,jcycle,cphi2(point2(jj)+1),crhs2(point2(jj)+1),NRANK,'==============23==============='
           end if
        end do
      END IF
   end do
   !********:fordebug:************
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   !********:fordebug:************
     write(*,*) 'mglin4' , NRANK , j , jcycle
  end do
end do

!ncx = Ncellx+1
WRITE(nnn,'(I3.3)') NRANK
open(521,file='final'//nnn//'.dat')
do k=1,ncz; kk=(ncx+1)*(ncy+1)*k+point2(NGL)
do j=1,ncy; jj=(ncx+1)*j+kk
do i=1,ncx; ii = i+jj
   Phi(i,j,k) = cphi2(ii)
   write(521,*) Phi(i,j,k)
end do; end do; end do
close(521)
   !call saveu(Phi(1,1,1),nxc,nyc,nzc,2)

END SUBROUTINE mglin


SUBROUTINE BCsgr_MPI(u,nx,ny,nz,lx1,lx2,ly1,ly2,lz1,lz2)
USE mpivar
INCLUDE 'mpif.h'
DOUBLE PRECISION :: u(0:nx,0:ny,0:nz)
DOUBLE PRECISION  :: VECU
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
INTEGER :: LEFTt,RIGTt
character(4) name
character(3) nR
integer :: countin=0
!write(*,*) 'BC1',NRANK

write (name,'(i4.4)') countin
write (nR,'(i3.3)') NRANK
!open(510+NRANK,file='BCsgr'//name//nR//'.dat')
!do kkkk=0,nz
!   do jjjj=0,ny
!      do iiii=0,nx
!         write(510+NRANK,*) iiii,jjjj,kkkk,nx,u(iiii,jjjj,kkkk)
!      end do
!   end do
!end do
!close(510+NRANK)
!write(*,*) 'rMPI' , countin , NRANK
countin =countin + 1


CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

!*** for X ***!
LEFTt = LEFT; IF(IST.eq.0)        LEFTt = MPI_PROC_NULL
RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGTt = MPI_PROC_NULL
write(*,*) 'BC2',NRANK,LEFTt,RIGTt
CALL MPI_TYPE_VECTOR((ny+1)*(nz+1),1,nx+1,MPI_REAL8,VECU,IERR); CALL MPI_TYPE_COMMIT(VECU,IERR)
if(lx1.eq.1) CALL MPI_SENDRECV(u(nx-1,0,0),1,VECU,RIGTt,1, &
                               u(   0,0,0),1,VECU,LEFTt,1, MPI_COMM_WORLD,MSTATUS,IERR)
if(lx2.eq.1) CALL MPI_SENDRECV(u(   1,0,0),1,VECU,LEFTt,1, &
                               u(nx  ,0,0),1,VECU,RIGTt,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
!write(*,*) 'BC2',NRANK
!*** for Y ***!
!IF(JST.eq.0)        BOTM = MPI_PROC_NULL
!IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
CALL MPI_TYPE_VECTOR(nz+1,nx+1,(nx+1)*(ny+1),MPI_REAL8,VECU,IERR); CALL MPI_TYPE_COMMIT(VECU,IERR)
if(ly1.eq.1) CALL MPI_SENDRECV(u(0,ny-1,0),1,VECU,TOP ,1, &
                               u(0,   0,0),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
if(ly2.eq.1) CALL MPI_SENDRECV(u(0,   1,0),1,VECU,BOTM,1, &
                               u(0,ny  ,0),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
!write(*,*) 'BC3',NRANK
!*** for z ***!
!IF(KST.eq.0)        DOWN = MPI_PROC_NULL
!IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
CALL MPI_TYPE_VECTOR(1,(nx+1)*(ny+1),(nx+1)*(ny+1),MPI_REAL8,VECU,IERR); CALL MPI_TYPE_COMMIT(VECU,IERR)
if(lz1.eq.1) CALL MPI_SENDRECV(u(0,0,nz-1),1,VECU,UP  ,1, &
                               u(0,0,   0),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
if(lz2.eq.1) CALL MPI_SENDRECV(u(0,0,   1),1,VECU,DOWN,1, &
                               u(0,0,nz  ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'BC4',NRANK
END SUBROUTINE BCsgr_MPI


SUBROUTINE collect(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2,pos)
USE mpivar
INCLUDE 'mpif.h'
integer :: countin=0, nx1,ny1,nz1,nx2,ny2,nz2,pos!,MSTATUS(MPI_STATUS_SIZE)
double precision :: u2(0:nx2,0:ny2,0:nz2) , u1(nx1,ny1,nz1)
!double precision ,intent(in) :: u2(0:nx2,0:ny2,0:nz2)!,u2(0:nx2,0:ny2,0:nz2)
!double precision :: u3(1:(nx2+1)*(ny2+1)*(nz2+1))
double precision :: tMPI(0:nx2,0:ny2,0:nz2,0:NPE-1)
!integer :: countin=0!,nx1,ny1,nz1,nx2,ny2,nz2
character*4 name
character*3 nR
write(*,*) 'cl=======1'
goto 3001
write (name,'(I4.4)') countin
write (nR,'(i3.3)') NRANK
open(10+NRANK,file='/work/maedarn/3DMHD/test/collectpre'//name//nR//'.dat')
!goto 3000
do k=0,nz2
   do j=0,ny2
      do i=0,nx2
         write(10+NRANK,*) pos,nx1,nx2,i,j,k,u2(i,j,k)
      end do
   end do
end do
!3000 continue

!do i=1,(nx2+1)*(ny2+1)*(nz2+1)
!   write(131+NRANK,*) i , u3(i)
!end do

!do k=0,nz2
!   do j=0,ny2
!      do i=0,nx2
!         u2(i,j,k)=u3( (i+1) + nx2*(j) +  nx2*ny2*(k) )
!      end do
!   end do
!end do
close(10+NRANK)
3001 continue

do k=0,nz2; do j=0,ny2; do i=0,nx2
  tMPI(i,j,k,NRANK)=u2(i,j,k)
end do;end do;end do
countin =countin + 1
!**********NEW FORM***********
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!**********NEW FORM***********

do Nroot=0,NPE-1
  CALL MPI_BCAST(tMPI(0,0,0,Nroot),(nx2+1)*(ny2+1)*(nz2+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

!**********NEW FORM***********
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!**********NEW FORM***********

do Nroot=0,NPE-1
 ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
 nxed=nx2-1; IF(ISTt.eq.NSPLTx-1) nxed=nx2
 nyed=ny2-1; IF(JSTt.eq.NSPLTy-1) nyed=ny2
 nzed=nz2-1; IF(KSTt.eq.NSPLTz-1) nzed=nz2
 do kk=1,nzed;k=KSTt*(nz2-1)+kk
  do jj=1,nyed;j=JSTt*(ny2-1)+jj
   do ii=1,nxed;i=ISTt*(nx2-1)+ii
    u1(i,j,k) = tMPI(ii,jj,kk,Nroot)
 end do;end do;end do;end do

 !write (name,'(i4.4)') countin
 !write (nR,'(i3.3)') NRANK
 !goto 1202
! open(251+NRANK,file='collect'//name//nR//'.dat')

 !do k=1,nz1; do j=1,ny1; do i=1,nx1
 !        write(251+NRANK,*) nx1, i,j,k, u1(i,j,k)
 !     end do
 !  end do
!end do
 !do kk=1,nzed;k=KSTt*(nz2-1)+kk
 ! do jj=1,nyed;j=JSTt*(ny2-1)+jj
 !  do ii=1,nxed;i=ISTt*(nx2-1)+ii
 !         write(251+NRANK,*) ii,jj,kk, u1(i,j,k)!,tMPI(ii,jj,kk,Nroot)
 !      end do
 !   end do
 !end do
 !close(251+NRANK)
 !write(*,*) 'collect' , countin , NRANK
! countin =countin + 1
 !1202 continue
 write(*,*) 'cl=======2'
END SUBROUTINE collect


SUBROUTINE divide(u2,u1,nx2,ny2,nz2,nx1,ny1,nz1)
USE mpivar
double precision :: u2(0:nx2,0:ny2,0:nz2),u1(nx1,ny1,nz1)
integer :: countin=0
character(4) name
character(3) nR
!write(*,*) 'rsmp=======1'
write (name,'(i4.4)') countin
write (nR,'(i3.3)') NRANK
open(171+NRANK,file='dvd'//name//nR//'.dat')

write(*,*) 'di======1'
!**********NEW FORM***********
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!**********NEW FORM***********

iist=0; IF(IST.eq.0) iist=1
jjst=0; IF(JST.eq.0) jjst=1
kkst=0; IF(KST.eq.0) kkst=1
do kk=kkst,nz2;k=KST*(nz2-1)+kk
  do jj=jjst,ny2;j=JST*(ny2-1)+jj
    do ii=iist,nx2;i=IST*(nx2-1)+ii
      u2(ii,jj,kk) = u1(i,j,k)
   end do;end do;end do

   do k=0,nz2; do j=0,ny2; do i=0,nx2
         write(171+NRANK,*) i,j,k, u2(i,j,k)
      end do
   end do
end do
close(171+NRANK)
countin = countin+1
write(*,*) 'di=======2'
END SUBROUTINE divide


SUBROUTINE rstrctMPI(uc,uf,nx,ny,nz,mode)
USE mpivar
double precision :: uc(0:nx,0:ny,0:nz),uf(0:2*nx-1,0:2*ny-1,0:2*nz-1)
double precision, parameter :: w = 1.d0/12.d0
integer :: countin=0
character(4) name
character(3) nR
write(*,*) 'rsmp=======1'

goto 3002
write (name,'(i4.4)') countin
write (nR,'(i3.3)') NRANK
open(171+NRANK,file='rstrctpre'//name//nR//'.dat')
do k=0,2*nz-1; do j=0,2*ny-1; do i=0,2*nx-1
         write(171+NRANK,*) i,j,k, uf(i,j,k)
      end do
   end do
end do
close(171+NRANK)
countin = countin+1
3002 continue

ixst=1; ixed=nx-1; iyst=1; iyed=ny-1; izst=1; ized=nz-1
IF(IST.eq.0) ixst = 2; IF(JST.eq.0) iyst = 2; IF(KST.eq.0) izst = 2 !境界は導出済

do kc=izst,ized; kf=2*kc-1
  do jc=iyst,iyed; jf=2*jc-1
    do ic=ixst,ixed; if=2*ic-1
      uc(ic,jc,kc)=0.5d0*uf(if,jf,kf)+ &
      w*(uf(if+1,jf,kf)+uf(if-1,jf,kf)+uf(if,jf+1,kf)+uf(if,jf-1,kf)+uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do


IF(IST.eq.0) THEN
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(1 ,jc,kc)=uf(1 ,jf,kf)
  end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
  nf=2*nx-1
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(nx,jc,kc)=uf(nf,jf,kf)
  end do; end do
END IF
IF(JST.eq.0) THEN
  do kc=1,nz; kf=2*kc-1; do ic=1,nx; if=2*ic-1
    uc(ic,1 ,kc)=uf(if,1 ,kf)
  end do; end do
END IF
IF(JST.eq.NSPLTy-1) THEN
  nf=2*ny-1
  do kc=1,nz; kf=2*kc-1; do ic=1,nx; if=2*ic-1
    uc(ic,ny,kc)=uf(if,nf,kf)
  end do; end do
END IF
IF(KST.eq.0) THEN
  do jc=1,ny; jf=2*jc-1; do ic=1,nx; if=2*ic-1
    uc(ic,jc,1 )=uf(if,jf,1 )
  end do; end do
END IF
IF(KST.eq.NSPLTz-1) THEN
  nf=2*nz-1
  do jc=1,ny; jf=2*jc-1; do ic=1,nx; if=2*ic-1
    uc(ic,jc,nz)=uf(if,jf,nf)
  end do; end do
END IF

IF(mode.eq.1) THEN

IF(IST.eq.0) THEN
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(1 ,jc,kc)=0.d0
  end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
  nf=2*nx-1
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(nx,jc,kc)=0.d0
  end do; end do
END IF

END IF

CALL BCsgr_MPI(uc,nx,ny,nz,1,1,1,1,1,1)

!open(271+NRANK,file='rstrct'//name//nR//'.dat')
!do k=0,nz; do j=0,ny; do i=0,nx
!         write(271+NRANK,*) i,j,k, uc(i,j,k)
!      end do
!   end do
!end do
!close(271+NRANK)
write(*,*) 'rsmp=======2'
END SUBROUTINE rstrctMPI


SUBROUTINE rstrct(uc,uf,nx,ny,nz,mode,NRANK)
double precision uc(nx,ny,nz),uf(2*nx-1,2*ny-1,2*nz-1)
double precision, parameter :: w = 1.d0/12.d0
integer :: countin=0,NRANK
character(4) name
character(3) nR
write(*,*) 'rs=======1'

goto 3003
write (name,'(i4.4)') countin
write (nR,'(i3.3)') NRANK
!open(171+NRANK,file='rstrctnompi'//name//nR//'.dat')

open(371+NRANK,file='rstrctnompipre'//name//nR//'.dat')
do k=1,2*nz-1; do j=1,2*ny-1; do i=1,2*nx-1
         write(371+NRANK,*)nx, i,j,k, uf(i,j,k)
      end do
   end do
end do
close(371+NRANK)

3003 continue

do kc=2,nz-1; kf=2*kc-1
  do jc=2,ny-1; jf=2*jc-1
    do ic=2,nx-1; if=2*ic-1
      uc(ic,jc,kc)=0.5d0*uf(if,jf,kf)+ &
      w*(uf(if+1,jf,kf)+uf(if-1,jf,kf)+uf(if,jf+1,kf)+uf(if,jf-1,kf)+uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do

nf=2*nx-1
do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
  uc(1 ,jc,kc)=uf(1 ,jf,kf)
  uc(nx,jc,kc)=uf(nf,jf,kf)
end do; end do
nf=2*nz-1
do jc=1,ny; jf=2*jc-1; do ic=1,nx; if=2*ic-1
  uc(ic,jc,1 )=uf(if,jf,1 )
  uc(ic,jc,nz)=uf(if,jf,nf)
end do; end do
nf=2*ny-1
do kc=1,nz; kf=2*kc-1; do ic=1,nx; if=2*ic-1
  uc(ic,1 ,kc)=uf(if,1 ,kf)
  uc(ic,ny,kc)=uf(if,nf,kf)
end do; end do

IF(mode.eq.1) THEN

do kc=1,nz; do jc=1,ny
  uc(1 ,jc,kc) = 0.d0
  uc(nx,jc,kc) = 0.d0
end do; end do

goto 3004
open(171+NRANK,file='rstrctnompi'//name//nR//'.dat')
do k=1,nz; do j=1,ny; do i=1,nx
         write(171+NRANK,*) i,j,k, uc(i,j,k)
      end do
   end do
end do
close(171+NRANK)
countin = countin+1
3004 continue
END IF
write(*,*) 'rs=======2'
END SUBROUTINE rstrct


SUBROUTINE interpMPI(uf,uc,nx,ny,nz,np,mode)
USE mpivar
USE slfgrv
double precision uf(0:nx,0:ny,0:nz),uc(0:nx/2+1,0:ny/2+1,0:nz/2+1)
write(*,*) 'inmp=======1'

do kc=1,nz/2+1; kf=2*kc-1
  do jc=1,ny/2+1; jf=2*jc-1
    do ic=1,nx/2+1
      uf(2*ic-1,jf,kf)=uc(ic,jc,kc)
    end do
  end do
end do

do kf=1,nz,2
  do jf=1,ny,2
    do if=2,nx,2
      uf(if,jf,kf)=.5d0*(uf(if+1,jf,kf)+uf(if-1,jf,kf))
    end do
  end do
end do
do kf=1,nz,2
  do jf=2,ny,2
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf+1,kf)+uf(if,jf-1,kf))
    end do
  end do
end do
do kf=2,nz,2
  do jf=1,ny
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do

IF(mode.eq.1) THEN

IF(IST.eq.0) THEN
do k=0,nz; kk=(ny+1)*k
do j=0,ny; nn=j+kk
  uf(1 ,j,k)=bphi2(np+nn,1)
end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
do k=0,nz; kk=(ny+1)*k
do j=0,ny; nn=j+kk
  uf(nx,j,k)=bphi2(np+nn,2)
end do; end do
END IF

END IF

CALL BCsgr_MPI(uf,nx,ny,nz,1,0,1,0,1,0)
write(*,*) 'inmp=======2'
END SUBROUTINE interpMPI


SUBROUTINE interp(uf,uc,nx,ny,nz,np,mode,NRANK)
USE slfgrv
double precision uf(nx,ny,nz),uc(nx/2+1,ny/2+1,nz/2+1)
integer :: countin=0,NRANK
character(4) name
character(3) nR
write(*,*) 'in=======1'
write (name,'(i4.4)') countin
write (nR,'(i3.3)') NRANK


do kc=1,nz/2+1; kf=2*kc-1
  do jc=1,ny/2+1; jf=2*jc-1
    do ic=1,nx/2+1
      uf(2*ic-1,jf,kf)=uc(ic,jc,kc)
    end do
  end do
end do

do kf=1,nz,2
  do jf=1,ny,2
    do if=2,nx,2
      uf(if,jf,kf)=.5d0*(uf(if+1,jf,kf)+uf(if-1,jf,kf))
    end do
  end do
end do
do kf=1,nz,2
  do jf=2,ny,2
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf+1,kf)+uf(if,jf-1,kf))
    end do
  end do
end do
do kf=2,nz,2
  do jf=1,ny
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do

IF(mode.eq.1) THEN
  do kf=1,nz; kk=ny*(kf-1)
  do jf=1,ny; nn=jf-1+kk
    uf(1 ,jf,kf) = bphi1(np+nn,1)
    uf(nx,jf,kf) = bphi1(np+nn,2)
  end do; end do
END IF

open(171+NRANK,file='interp'//name//nR//'.dat')
do k=1,nz; do j=1,ny; do i=1,nx
         write(171+NRANK,*) mode,i,j,k, uf(i,j,k)
      end do
   end do
end do
close(171+NRANK)
countin = countin+1
write(*,*) 'in=======1'
END SUBROUTINE interp


SUBROUTINE addintMPI(uf,uc,res,nx,ny,nz)
  double precision res(0:nx,0:ny,0:nz),uc(0:nx/2+1,0:ny/2+1,0:nz/2+1),uf(0:nx,0:ny,0:nz)
  integer :: up=0
call interpMPI(res,uc,nx,ny,nz,0,0)
!integer :: up=0
!do k=0,nz; do j=0,ny; do i=0,nx
!  uf(i,j,k)=uf(i,j,k)+res(i,j,k)
!end do; end do; end do
write(*,*) 'addmp=======1'
!isw=2 !for speed up
!isw=2-mod(up,2)
!up=up+1
!write(*,*) isw ,'iswiswmpi'
do jsw=2,1,-1
isw=jsw
do k=0,nz
  do j=0,ny
    do i=isw-1,nx,2
      uf(i,j,k)=uf(i,j,k)+res(i,j,k)
    end do
    isw=3-isw
  end do
  isw=3-isw
end do
end do
write(*,*) 'addmp=======2'
END SUBROUTINE addintMPI


SUBROUTINE addint(uf,uc,res,nx,ny,nz,NRANK)
  double precision res(nx,ny,nz),uc(nx/2+1,ny/2+1,nz/2+1),uf(nx,ny,nz)
  integer NRANK
  integer :: up=0
  write(*,*) 'add=======1'
call interp(res,uc,nx,ny,nz,0,0,NRANK)
!do k=1,nz; do j=1,ny; do i=1,nx
!uf(i,j,k)=uf(i,j,k)+res(i,j,k)
!end do; end do; end do

!isw=1 ! for speed up
!isw=2-mod(up,2)
!up=up+1
!write(*,*) isw ,'iswisw'
do jsw=1,2
   isw=jsw
do k=1,nz
  do j=1,ny
    do i=isw,nx,2
      uf(i,j,k)=uf(i,j,k)+res(i,j,k)
    end do
    isw=3-isw
  end do
end do
end do
write(*,*) 'add=======2'
END SUBROUTINE addint


SUBROUTINE slvsml(u,rhs)
USE slfgrv
double precision  u(3,3,3),rhs(3,3,3)
double precision h
write(*,*) 'sl=======1'
h=0.5d0*Lbox

u(1,1,1)=0.d0;u(3,1,1)=0.d0
u(1,2,1)=0.d0;u(3,2,1)=0.d0
u(1,3,1)=0.d0;u(3,3,1)=0.d0

u(1,1,2)=0.d0;u(3,1,2)=0.d0
u(1,2,2)=0.d0;u(3,2,2)=0.d0
u(1,3,2)=0.d0;u(3,3,2)=0.d0

u(1,1,3)=0.d0;u(3,1,3)=0.d0
u(1,2,3)=0.d0;u(3,2,3)=0.d0
u(1,3,3)=0.d0;u(3,3,3)=0.d0

u(2,1,1) = 0.025d0*h*h*(-8.d0*rhs(2,1,1)-2.d0*rhs(2,1,2)-2.d0*rhs(2,1,3)-2.d0*rhs(2,2,1)-rhs(2,2,2)-rhs(2,2,3)-2.d0*rhs(2,3,1)-rhs(2,3,2)-rhs(2,3,3))
u(2,2,1) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-rhs(2,1,2)-rhs(2,1,3)-8.d0*rhs(2,2,1)-2.d0*rhs(2,2,2)-2.d0*rhs(2,2,3)-2.d0*rhs(2,3,1)-rhs(2,3,2)-rhs(2,3,3))
u(2,3,1) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-rhs(2,1,2)-rhs(2,1,3)-2.d0*rhs(2,2,1)-rhs(2,2,2)-rhs(2,2,3)-8.d0*rhs(2,3,1)-2.d0*rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,1,2) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-8.d0*rhs(2,1,2)-2.d0*rhs(2,1,3)-rhs(2,2,1)-2.d0*rhs(2,2,2)-rhs(2,2,3)-rhs(2,3,1)-2.d0*rhs(2,3,2)-rhs(2,3,3))
u(2,2,2) = 0.025d0*h*h*(-rhs(2,1,1)-2.d0*rhs(2,1,2)-rhs(2,1,3)-2.d0*rhs(2,2,1)-8.d0*rhs(2,2,2)-2.d0*rhs(2,2,3)-rhs(2,3,1)-2.d0*rhs(2,3,2)-rhs(2,3,3))
u(2,3,2) = 0.025d0*h*h*(-rhs(2,1,1)-2.d0*rhs(2,1,2)-rhs(2,1,3)-rhs(2,2,1)-2.d0*rhs(2,2,2)-rhs(2,2,3)-2.d0*rhs(2,3,1)-8.d0*rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,1,3) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-2.d0*rhs(2,1,2)-8.d0*rhs(2,1,3)-rhs(2,2,1)-rhs(2,2,2)-2.d0*rhs(2,2,3)-rhs(2,3,1)-rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,2,3) = 0.025d0*h*h*(-rhs(2,1,1)-rhs(2,1,2)-2.d0*rhs(2,1,3)-2.d0*rhs(2,2,1)-2.d0*rhs(2,2,2)-8.d0*rhs(2,2,3)-rhs(2,3,1)-rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,3,3) = 0.025d0*h*h*(-rhs(2,1,1)-rhs(2,1,2)-2.d0*rhs(2,1,3)-rhs(2,2,1)-rhs(2,2,2)-2.d0*rhs(2,2,3)-2.d0*rhs(2,3,1)-2.d0*rhs(2,3,2)-8.d0*rhs(2,3,3))
write(*,*) 'sl=======2'
END SUBROUTINE slvsml


SUBROUTINE slvsmlb(u,rhs)
USE slfgrv
double precision  u(3,3,3),rhs(3,3,3)
double precision h
write(*,*) 'rslvb=======1'
h=0.5d0*Lbox

u(1,1,1)=bphi1(1,1);  u(3,1,1)=bphi1(1,2)
u(1,2,1)=bphi1(2,1);  u(3,2,1)=bphi1(2,2)
u(1,3,1)=bphi1(3,1);  u(3,3,1)=bphi1(3,2)

u(1,1,2)=bphi1(4,1);  u(3,1,2)=bphi1(4,2)
u(1,2,2)=bphi1(5,1);  u(3,2,2)=bphi1(5,2)
u(1,3,2)=bphi1(6,1);  u(3,3,2)=bphi1(6,2)

u(1,1,3)=bphi1(7,1);  u(3,1,3)=bphi1(7,2)
u(1,2,3)=bphi1(8,1);  u(3,2,3)=bphi1(8,2)
u(1,3,3)=bphi1(9,1);  u(3,3,3)=bphi1(9,2)

u(2,1,1) = 0.025d0* (8.d0*u(1,1,1) + 2.d0*u(1,1,2) + 2.d0*u(1,1,3) + 2.d0*u(1,2,1) + u(1,2,2) + u(1,2,3) + 2.d0*u(1,3,1) + &
u(1,3,2) + u(1,3,3) + 8.d0*u(3,1,1) + 2.d0*u(3,1,2) + 2.d0*u(3,1,3) + 2.d0*u(3,2,1) + u(3,2,2) + u(3,2,3) + 2.d0*u(3,3,1) + u(3,3,2) + u(3,3,3) &
 - 8.d0*h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,2,1) = 0.025d0* (2.d0*u(1,1,1) + u(1,1,2) + u(1,1,3) + 8.d0*u(1,2,1) + 2.d0*u(1,2,2) + 2.d0*u(1,2,3) + 2.d0*u(1,3,1) + &
u(1,3,2) + u(1,3,3) + 2.d0*u(3,1,1) + u(3,1,2) + u(3,1,3) + 8.d0*u(3,2,1) + 2.d0*u(3,2,2) + 2.d0*u(3,2,3) + 2.d0*u(3,3,1) + u(3,3,2) + u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - 8.d0*h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,3,1) = 0.025d0* (2.d0*u(1,1,1) + u(1,1,2) + u(1,1,3) + 2.d0*u(1,2,1) + u(1,2,2) + u(1,2,3) + 8.d0*u(1,3,1) + 2.d0*u(1,3,2) + &
2.d0*u(1,3,3) + 2.d0*u(3,1,1) + u(3,1,2) + u(3,1,3) + 2.d0*u(3,2,1) + u(3,2,2) + u(3,2,3) + 8.d0*u(3,3,1) + 2.d0*u(3,3,2) + 2.d0*u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
8.d0*h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,1,2) = 0.025d0* (2.d0*u(1,1,1) + 8.d0*u(1,1,2) + 2.d0*u(1,1,3) + u(1,2,1) + 2.d0*u(1,2,2) + u(1,2,3) + u(1,3,1) + 2.d0*u(1,3,2) + &
u(1,3,3) + 2.d0*u(3,1,1) + 8.d0*u(3,1,2) + 2.d0*u(3,1,3) + u(3,2,1) + 2.d0*u(3,2,2) + u(3,2,3) + u(3,3,1) + 2.d0*u(3,3,2) + u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - 8.d0*h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,2,2) = 0.025d0* ( u(1,1,1) + 2.d0*u(1,1,2) + u(1,1,3) + 2.d0*u(1,2,1) + 8.d0*u(1,2,2) + 2.d0*u(1,2,3) + u(1,3,1) + 2.d0*u(1,3,2) + &
u(1,3,3) + u(3,1,1) + 2.d0*u(3,1,2) + u(3,1,3) + 2.d0*u(3,2,1) + 8.d0*u(3,2,2) + 2.d0*u(3,2,3) + u(3,3,1) + 2.d0*u(3,3,2) + u(3,3,3) &
- h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - 8.d0*h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,3,2) = 0.025d0* ( u(1,1,1) + 2.d0*u(1,1,2) + u(1,1,3) + u(1,2,1) + 2.d0*u(1,2,2) + u(1,2,3) + 2.d0*u(1,3,1) + 8.d0*u(1,3,2) + &
2.d0*u(1,3,3) + u(3,1,1) + 2.d0*u(3,1,2) + u(3,1,3) + u(3,2,1) + 2.d0*u(3,2,2) + u(3,2,3) + 2.d0*u(3,3,1) + 8.d0*u(3,3,2) + 2.d0*u(3,3,3) &
- h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - 8.d0*h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,1,3) = 0.025d0* (2.d0*u(1,1,1) + 2.d0*u(1,1,2) + 8.d0*u(1,1,3) + u(1,2,1) + u(1,2,2) + 2.d0*u(1,2,3) + u(1,3,1) + u(1,3,2) + &
2.d0*u(1,3,3) + 2.d0*u(3,1,1) + 2.d0*u(3,1,2) + 8.d0*u(3,1,3) + u(3,2,1) + u(3,2,2) + 2.d0*u(3,2,3) + u(3,3,1) + u(3,3,2) + 2.d0*u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - 8.d0*h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,2,3) = 0.025d0* ( u(1,1,1) + u(1,1,2) + 2.d0*u(1,1,3) + 2.d0*u(1,2,1) + 2.d0*u(1,2,2) + 8.d0*u(1,2,3) + u(1,3,1) + u(1,3,2) + &
2.d0*u(1,3,3) + u(3,1,1) + u(3,1,2) + 2.d0*u(3,1,3) + 2.d0*u(3,2,1) + 2.d0*u(3,2,2) + 8.d0*u(3,2,3) + u(3,3,1) + u(3,3,2) + 2.d0*u(3,3,3) &
- h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - 8.d0*h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,3,3) = 0.025d0* ( u(1,1,1) + u(1,1,2) + 2.d0*u(1,1,3) + u(1,2,1) + u(1,2,2) + 2.d0*u(1,2,3) + 2.d0*u(1,3,1) + 2.d0*u(1,3,2) + &
8.d0*u(1,3,3) + u(3,1,1) + u(3,1,2) + 2.d0*u(3,1,3) + u(3,2,1) + u(3,2,2) + 2.d0*u(3,2,3) + 2.d0*u(3,3,1) + 2.d0*u(3,3,2) + 8.d0*u(3,3,3) &
- h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - 8.d0*h*h*rhs(2,3,3))

do k=1,3
   do j=1,3
      do i=1,3
         write(*,*) i,j,k, u(i,j,k),rhs(i,j,k)
      end do
   end do
end do
write(*,*) 'slvb=======2'
END SUBROUTINE slvsmlb


SUBROUTINE relaxMPI(u,rhs,nx,ny,nz,mode)
USE mpivar
USE slfgrv

!=====fordbug
INCLUDE 'mpif.h'

double precision, parameter :: w=1.d0/6.d0
double precision u(0:nx,0:ny,0:nz),rhs(0:nx,0:ny,0:nz)
double precision h,h2
integer :: countin=0
character(4) name
character(3) nR
h=Lbox/dble((nx-1)*NSPLTx)
h2=h*h

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

write(*,*) 'rxmp=======1' ,NRANK,h2,u(0,0,0)

!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************
!call saveu(u(0:nx,0:ny,0:nz),nx,nz,ny,0,0,0,1)
!goto 3005
write (name,'(i4.4)') countin
write (nR,'(i3.3)') NRANK
!open(751+NRANK,file='rMPIpre'//name//nR//'.dat')
!do kkkk=0,nz
!   do jjjj=0,ny
!      do iiii=0,nx
!         write(751+NRANK,*) iiii,jjjj,kkkk,nx,u(iiii,jjjj,kkkk),rhs(iiii,jjjj,kkkk)
!      end do
!   end do
!end do
!close(751+NRANK)
write(*,*) 'rMPI' , countin , NRANK
countin =countin + 1
!3005 continue
!********:fordebug:************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********:fordebug:************

li=0; IF(IST.eq.0) li=1
lj=0; IF(JST.eq.0) lj=1
lk=0; IF(KST.eq.0) lk=1
do jsw=2,1,-1 !Red-Black Gauss-Seidel
  isw=jsw
  if(mode.eq.1) then
    do k=1,nz-1
      do j=1,ny-1
        do i=isw,nx-1,2
          ifl=li*int(1/i)
          u(i,j,k)=-w*h2*rhs(i,j,k) &
          *dble(1-ifl) + (0.5d0+dsign(0.5d0,dble(ifl)-0.5d0))*u(i,j,k)
!		  if(NRANK.eq.0) write(*,*) 'i=',i,'ifl=',ifl
!		  if(NRANK.eq.0) write(*,*) (1-ifl),(0.5d0+dsign(0.5d0,ifl-0.5d0))
        end do
        isw=3-isw
      end do
      isw=3-isw
   end do

   !********:fordebug:************
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   !********:fordebug:************
   write(*,*) 'rPRE1',NRANK , mode
   mode=2
  else
    do k=1,nz-1
      do j=1,ny-1
        do i=isw,nx-1,2
          ifl=li*int(1/i)
          u(i,j,k)=w*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-h2*rhs(i,j,k)) &
          *dble(1-ifl) + (0.5d0+dsign(0.5d0,dble(ifl)-0.5d0))*u(i,j,k)
        end do
        isw=3-isw
      end do
      isw=3-isw
    end do
 end if
 write(*,*)'rPRE2', NRANK , mode
 CALL BCsgr_MPI(u,nx,ny,nz,1,1,1,1,1,1)
 write(*,*)'rPRE3', NRANK , mode
end do
write(*,*) 'rxmp=======2'
END SUBROUTINE relaxMPI


SUBROUTINE relax(u,rhs,nx,ny,nz,mode,a)
USE slfgrv
double precision, parameter :: w=1.d0/6.d0
double precision u(nx,ny,nz),rhs(nx,ny,nz)
double precision h,h2
integer :: countin=0,a
character(3) nR
character(4) name
write(*,*) 'rx=======1'
goto 1010
write (name,'(i4.4)') countin
write (nR,'(i3.3)') a
open(511+a,file='rpre'//name//nR//'.dat')
do kkkk=1,nz
   do jjjj=1,ny
      do iiii=1,nx
         write(511+a,*) iiii,jjjj,kkkk,nx,u(iiii,jjjj,kkkk),rhs(iiii,jjjj,kkkk)
      end do
   end do
end do
close(511+a)
!write(*,*) 'rMPI' , countin , NRANK
countin =countin + 1
1010 continue
write(*,*) a,'relax=======1'
h=Lbox/dble(nx-1)
h2=h*h

do jsw=1,2 !Red-Black Gauss-Seidel
  isw=jsw
  if(mode.eq.1) then
    do k=1,nz
      do j=1,ny
        do i=isw+1,nx-1,2
          u(i,j,k)=-w*h2*rhs(i,j,k)
        end do
        isw=3-isw
      end do
    end do
    mode=2
  else
    do k=1,nz
      km = k-1; if(k.eq.1) km = nz
!      km = (k-1)*(1+isign(1,k-2)   )/2 + nz*(1-isign(1,k-2)   )/2 !original
      kp = k+1; if(k.eq.nz) kp = 1
!      kp = (k+1)*(1+isign(1,nz-1-k))/2 +  1*(1-isign(1,nz-1-k))/2 !original
      do j=1,ny
        jm = j-1; if(j.eq.1) jm = ny
!        jm = (j-1)*(1+isign(1,j-2)   )/2 + ny*(1-isign(1,j-2)   )/2 !original
        jp = j+1; if(j.eq.ny) jp = 1
!        jp = (j+1)*(1+isign(1,ny-1-j))/2 +  1*(1-isign(1,ny-1-j))/2 !original
        do i=isw+1,nx-1,2
           u(i,j,k)=w*(u(i+1,j,k)+u(i-1,j,k)+u(i,jp,k)+u(i,jm,k)+u(i,j,kp)+u(i,j,km)-h2*rhs(i,j,k))
          !qq write(*,*) a,i,i-1,'relaxinGG'
        end do
        isw=3-isw
      end do
    end do
  end if
end do
write(*,*) a,'relax=======2'
END SUBROUTINE relax


SUBROUTINE residMPI(res,u,rhs,nx,ny,nz)
USE mpivar
USE slfgrv
double precision res(0:nx,0:ny,0:nz),rhs(0:nx,0:ny,0:nz),u(0:nx,0:ny,0:nz)
double precision h,h2i
write(*,*) 'resmp=======1'
h=Lbox/dble((nx-1)*NSPLTx)
h2i=1.d0/(h*h)

!is=1; IF(IST.eq.0) is=2
!js=1; IF(JST.eq.0) js=2
!ks=1; IF(KST.eq.0) ks=2
!do k=ks,nz-1; do j=js,ny-1; do i=is,nx-1
!  res(i,j,k)=-h2i*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-6.d0*u(i,j,k))+rhs(i,j,k) 
!end do; end do; end do

do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
  res(i,j,k)=-h2i*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-6.d0*u(i,j,k))+rhs(i,j,k) 
end do; end do; end do

IF(IST.eq.0) THEN
!do k=1,nz,2; do j=1,ny,2 !for speed up
do k=1,nz; do j=1,ny
  res(1 ,j,k)=0.d0 
end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
!do k=1,nz,2; do j=1,ny,2 !for speed up
do k=1,nz; do j=1,ny
  res(nx,j,k)=0.d0
end do; end do
END IF

CALL BCsgr_MPI(res,nx,ny,nz,1,1,1,1,1,1)
write(*,*) 'resmp=======1'
END SUBROUTINE residMPI


SUBROUTINE resid(res,u,rhs,nx,ny,nz)
USE slfgrv
double precision res(nx,ny,nz),rhs(nx,ny,nz),u(nx,ny,nz)
double precision h,h2i
write(*,*) 'res=======1'
h=Lbox/dble(nx-1)
h2i=1.d0/(h*h)
do k=1,nz; do j=1,ny; do i=2,nx-1
  res(i,j,k)=-h2i*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-6.d0*u(i,j,k))+rhs(i,j,k) 
end do; end do; end do

!do k=1,nz,2; do j=1,ny,2 !for speed up
do k=1,nz; do j=1,ny
  res(1 ,j,k)=0.d0
  res(nx,j,k)=0.d0
end do; end do
write(*,*) 'res=======2'
END SUBROUTINE resid


SUBROUTINE copy(aout,ain,nx,ny,nz)
  double precision ain(nx,ny,nz),aout(nx,ny,nz)
  write(*,*) 'cp=======1'
do k=1,nz; do j=1,ny; do i=1,nx; aout(i,j,k)=ain(i,j,k); end do; end do; end do
END SUBROUTINE copy

SUBROUTINE copyMPI(aout,ain,nx,ny,nz)
  double precision ain(0:nx,0:ny,0:nz),aout(0:nx,0:ny,0:nz)
  write(*,*) 'cpmp=======1'
do k=0,nz; do j=0,ny; do i=0,nx; aout(i,j,k)=ain(i,j,k); end do; end do; end do
END SUBROUTINE copyMPI

SUBROUTINE fill0(u,nx,ny,nz)
  double precision u(nx,ny,nz)
  write(*,*) 'fi=======1'
  !===============
  u(:,:,:)=0.0d0
  !===============
do j=1,ny; do i=1,nx
  u(i,j,1 )=0.d0
  u(i,j,nz)=0.d0
end do; end do
do k=1,nz; do i=1,nx
  u(i,1 ,k)=0.d0
  u(i,ny,k)=0.d0
end do; end do
do k=1,nz; do j=1,ny
  u(1 ,j,k)=0.d0
  u(nx,j,k)=0.d0
end do; end do
END SUBROUTINE fill0

SUBROUTINE fill0MPI(u,nx,ny,nz)
  double precision u(0:nx,0:ny,0:nz)
  write(*,*) 'fimp=======1'
  !===============
  u(:,:,:)=0.0d0
  !===============
do j=0,ny; do i=0,nx
  u(i,j,0 )=0.d0
  u(i,j,1 )=0.d0
  u(i,j,nz)=0.d0
end do; end do
do k=0,nz; do i=0,nx
  u(i,0 ,k)=0.d0
  u(i,1 ,k)=0.d0
  u(i,ny,k)=0.d0
end do; end do
do k=0,nz; do j=0,ny
  u(0 ,j,k)=0.d0
  u(1 ,j,k)=0.d0
  u(nx,j,k)=0.d0
end do; end do
END SUBROUTINE fill0MPI


SUBROUTINE PB()
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

double precision,  parameter :: pi = 3.14159265359d0
DOUBLE PRECISION, dimension(:,:,:), allocatable :: temp1,temp2

DOUBLE PRECISION, dimension(:,:,:),  allocatable :: data
complex*16, dimension(:,:), allocatable :: speq
DOUBLE PRECISION :: rho

DOUBLE PRECISION, dimension(:,:),  allocatable :: dat1,dat2
complex*16, dimension(:), allocatable :: spe1,spe2
double precision :: kap,temp1r,temp1i,temp2r,temp2i,facG,fac,dxx,dyy,dzz,zp1,zp2
double precision, dimension(:,:,:), allocatable :: fint0,fint1

character*4 fnum
integer :: countin=0
character(4) name
character(3) nR



iwx=1;iwy=1;iwz=1;N_MPI(20)=1;N_MPI(1)=1;CALL BC_MPI(2,1)

!*** Pointer for boundary ***!
pointb1(1) = 1
nl = 1
nc=3
2 continue
pointb1(nl+1)=pointb1(nl)+nc**2
nc=nc*2-1
nl = nl+1
if(nl.ne.NGcr+1) goto 2
Needb = pointb1(NGcr+1)
ALLOCATE(bphi1(Needb,2))

!nl = NGcr-1
nl = NGcr
pointb2(nl) = 1
nx=(2**NGcr)/NSPLTx+2; ny=(2**NGcr)/NSPLTy+2; nz=(2**NGcr)/NSPLTz+2
3 continue
pointb2(nl+1)=pointb2(nl)+max0(nx*ny,ny*nz,nz*nx)
nx=nx*2-2; ny=ny*2-2; nz=nz*2-2
nl = nl+1
if(nl.ne.NGL+1) goto 3
Needb = pointb2(NGL+1)
ALLOCATE(bphi2(Needb,2))


!MIYAMA method ---------------------------------------------------

ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1),speq(Ncellz*NSPLTz,Ncellx))
ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
         dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))

!nccy = Ncelly/NSPLTy; nccz = Ncellz/NSPLTz
nccy = Ncelly; nccz = Ncellz


do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
do i=1,Ncellx
  data(jy,kz,i) = U(i,j,k,1)
end do;end do;end do

!****************new*****************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!****************new*****************

                    !count, blocklength, stride
CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)

do Nlp = 1,NSPLTy*NSPLTz-1

  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

  Nis = JSs + NSPLTy*KSs
  kls = Nis + 1
  Nir = JST + NSPLTy*KST
  klr = Nir + 1

  write(*,*) NRANK,Nlp,isend,KSs,JSs,irecv,KSr,JSr,Nis,kls,Nir,klr , 'MMMMMMMMMMMPPPPPPPPPPPIIIIIIIIIIII'

  if(kls.gt.Ncellx) then; isend = MPI_PROC_NULL; kls = Ncellx+1; end if
  if(klr.gt.Ncellx) then; irecv = MPI_PROC_NULL; klr = Ncellx+1; end if
  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls),1,VECU,isend,1, & !send
                    data(JSr*Ncelly+1,KSr*Ncellz+1,klr),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv
end do

CALL MPI_TYPE_FREE(VECU,IERR)


dxx = dy(1); dyy = dz(1); dzz = dx(1)
facG = -G4pi*dzz

dat1(:,:)=0.d0; dat2(:,:)=0.d0; spe1(:)=(0.d0,0.d0); spe2(:)=(0.d0,0.d0)

!Nir = JST + NSPLTy*KST
!klr = Nir + 1


nn1 = Ncelly*NSPLTy; nn2 = Ncellz*NSPLTz

if(klr.le.Ncellx) then

  call rlft3(data(1,1,klr),speq(1,klr),nn1,nn2,1) !dataをfourier

  kz = klr
  zp1 = x(kz)-0.5d0*dzz
  zp2 = Lbox - zp1
  temp1r = dat1(1,1) - data(1,1,klr) * 0.5d0*zp1 * facG
  temp1i = dat1(2,1) - data(2,1,klr) * 0.5d0*zp1 * facG
  temp2r = dat2(1,1) - data(1,1,klr) * 0.5d0*zp2 * facG
  temp2i = dat2(2,1) - data(2,1,klr) * 0.5d0*zp2 * facG

  write (name,'(i4.4)') countin
  write (nR,'(i3.3)') NRANK
  open(171+NRANK,file='PBa1'//name//nR//'.dat')
  open(271+NRANK,file='PBb1'//name//nR//'.dat')

  do m=1,nn2/2+1; do l=1,nn1/2
    kap = 4.d0*( dsin(pi*(dble(l-1)/dble(nn1)))**2/dxx**2 + dsin(pi*(dble(m-1)/dble(nn2)))**2/dyy**2 )
    kap = dsqrt(kap) - 1.d-4!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(271+NRANK,*) kap,nn1,nn2,dyy,dxx,dat1(2*l-1,m),dat1(2*l  ,m),dat2(2*l-1,m),dat2(2*l  ,m),data(2*l-1,m,klr) &
    ,-zp1*kap,zp1,kz,dzz,0.5d0*dexp(-zp1*kap)/kap *facG
    dat1(2*l-1,m) = dat1(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*dexp(-zp1*kap)/kap *facG
    dat1(2*l  ,m) = dat1(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*dexp(-zp1*kap)/kap *facG
    dat2(2*l-1,m) = dat2(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*dexp(-zp2*kap)/kap *facG
    dat2(2*l  ,m) = dat2(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*dexp(-zp2*kap)/kap *facG
    write(171+NRANK,*) kap,nn1,nn2,dyy,dxx,dat1(2*l-1,m),dat1(2*l  ,m),dat2(2*l-1,m),dat2(2*l  ,m),-zp1*kap,zp1,kz,dzz,0.5d0*dexp(-zp1*kap)/kap *facG
  end do;end do

  !do k=0,2*nz-1; do j=0,2*ny-1; do i=0,2*nx-1
  !       write(171+NRANK,*) i,j,k, uf(i,j,k)
  !    end do
  ! end do
!end do
  close(171+NRANK)
  close(271+NRANK)
  countin = countin+1

  l=nn1/2+1
  do m=1,nn2/2+1
     kap = 4.d0*( dsin(pi*(dble(l-1)/dble(nn1)))**2/dxx**2 + dsin(pi*(dble(m-1)/dble(nn2)))**2/dyy**2 )
     kap = dsqrt(kap)- 1.d-4
     !write(*,*) NRANK,-zp1*kap
    spe1(m) = spe1(m) + speq(m,klr)* 0.5d0*dexp(-zp1*kap)/kap *facG
    spe2(m) = spe2(m) + speq(m,klr)* 0.5d0*dexp(-zp2*kap)/kap *facG
  end do

  do m2=nn2/2+2,nn2; m=nn2+2-m2; do l=1,nn1/2
    kap = 4.d0*( dsin(pi*(dble(l-1)/dble(nn1)))**2/dxx**2 + dsin(pi*(dble(m-1)/dble(nn2)))**2/dyy**2 )
    kap = dsqrt(kap) - 1.d-4
    !write(*,*) NRANK,-zp1*kap,'BP1'
    dat1(2*l-1,m2) = dat1(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*dexp(-zp1*kap)/kap *facG
    dat1(2*l  ,m2) = dat1(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*dexp(-zp1*kap)/kap *facG
    dat2(2*l-1,m2) = dat2(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*dexp(-zp2*kap)/kap *facG
    dat2(2*l  ,m2) = dat2(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*dexp(-zp2*kap)/kap *facG
  end do;end do

  l=nn1/2+1
  do m2=nn2/2+2,nn2; m=nn2+2-m2
    kap = 4.d0*( dsin(pi*(dble(l-1)/dble(nn1)))**2/dxx**2 + dsin(pi*(dble(m-1)/dble(nn2)))**2/dyy**2 )
    kap = dsqrt(kap)- 1.d-4 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !write(*,*) NRANK,-zp1*kap
    spe1(m2) = spe1(m2) + speq(m2,klr)* 0.5d0*dexp(-zp1*kap)/kap *facG
    spe2(m2) = spe2(m2) + speq(m2,klr)* 0.5d0*dexp(-zp2*kap)/kap *facG
  end do

  dat1(1,1) = temp1r
  dat1(2,1) = temp1i
  dat2(1,1) = temp2r
  dat2(2,1) = temp2i

end if

CALL MPI_ALLREDUCE(dat1(1,1),data(1,1,1),Ncelly*NSPLTy*Ncellz*NSPLTz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
CALL MPI_ALLREDUCE(spe1(1)  ,speq(  1,1),Ncellz*NSPLTz,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,IERR)
CALL MPI_ALLREDUCE(dat2(1,1),data(1,1,2),Ncelly*NSPLTy*Ncellz*NSPLTz,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,IERR)
CALL MPI_ALLREDUCE(spe2(1)  ,speq(  1,2),Ncellz*NSPLTz,MPI_COMPLEX16,MPI_SUM,MPI_COMM_WORLD,IERR)


call rlft3(data(1,1,1),speq(1,1),nn1,nn2,-1)
call rlft3(data(1,1,2),speq(1,2),nn1,nn2,-1)

fac = 2*(nn1/2)**2
do j=1,nn2; do i=1,nn1
  data(i,j,1) = data(i,j,1)/fac
  data(i,j,2) = data(i,j,2)/fac
end do; end do

!If(NRANK.eq.0) then; open(3,file='Pot.dat')
!do j=1,nn2
!  do i=1,nn1
!    write(3,*) i,j,data(i,j,1),data(i,j,2)
!  end do; write(3,*) ' '
!end do
!close(3); END IF


ncx=Ncellx+1; ncy=Ncelly+1; ncz=Ncellz+1
do k=0,ncz; kk= (ncy+1)*k
do j=0,ncy; n = j+kk
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 1
  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 1
  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz

  bphi2(pointb2(NGL)+n,1) = dble(data(jb,kbb,1))
  bphi2(pointb2(NGL)+n,2) = dble(data(jb,kbb,2))
end do; end do


!write(fnum,'(I3.3)') NRANK
!open(3,file='/work/inouety/SGte/bphi'//fnum//'.dat')
!ncx=Ncellx+1; ncy=Ncelly+1; ncz=Ncellz+1
!do k=0,ncz; kk= (ncy+1)*k
!do j=0,ncy; n = j+kk
!  write(3,*) JST*Ncelly+j,KST*Ncellz+k,bphi2(pointb2(NGL)+n,2)
!end do;write(3,*) ' '; end do
!close(3)


DEALLOCATE(data,speq)
DEALLOCATE(dat1,spe1,dat2,spe2)
!-----------------------------------------------------------------

open(171+NRANK,file='Phi1'//name//nR//'.dat')
open(271+NRANK,file='Phi2'//name//nR//'.dat')
open(371+NRANK,file='Phi3'//name//nR//'.dat')
open(471+NRANK,file='Phi4'//name//nR//'.dat')

ncx = Ncellx+1; ncy = Ncelly+1; ncz = Ncellz+1
w  = 0.125d0
do lc=NGL-1,NGcr-1,-1
  nfx = ncx; nfy = ncy; nfz = ncz
  ncx = ncx/2+1; ncy = ncy/2+1; ncz = ncz/2+1
  do l = 1, 2
  if(l.eq.1) then; mfx=nfy; mfy=nfz; mcx=ncy; mcy=ncz; end if
  if(l.eq.2) then; mfx=nfy; mfy=nfz; mcx=ncy; mcy=ncz; end if

  do jc = 2,mcy-1; jf = 2*jc-1
  do ic = 2,mcx-1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 4.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1    ,l)+bphi2(pointb2(lc+1)+kf+1    ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    !write(171+NRANK,*) lc,kf,kc,jc,ic,jf,if,l,bphi2(pointb2(lc)+kc,l)
 end do
  end do

  do jc = 2,mcy-1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    ! write(271+NRANK,*) 0,lc,kf,kc,jc,ic,jf,if,l,bphi2(pointb2(lc)+kc,l)
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
   ! write(271+NRANK,*) 1,lc,kf,kc,jc,ic,jf,if,l,bphi2(pointb2(lc)+kc,l)
 end do
  do ic = 2,mcx-1; if = 2*ic-1
    jc = 1; jf = 2*jc-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    jc = mcy; jf = 2*jc-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )
     !write(371+NRANK,*) lc,kf,kc,jc,ic,jf,if,l,bphi2(pointb2(lc)+kc,l)
 end do

    jc = 1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    jc = 1 ; jf = 2*jc-1
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    jc = mcy; jf = 2*jc-1
    ic = 1 ; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )
    jc = mcy; jf = 2*jc-1
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1 ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )

     !write(471+NRANK,*) 1,lc,kf,kc,jc,ic,jf,if,l,bphi2(pointb2(lc)+kc,l)
 end do
end do


!close(171+NRANK)
!close(271+NRANK)
!close(371+NRANK)
!close(471+NRANK)

nz=(2**NGcr)/NSPLTz+1; ny=(2**NGcr)/NSPLTy+1; nx=(2**NGcr)/NSPLTx+1; n1=2**NGcr+1
ALLOCATE( temp1(n1,n1,n1), temp2(0:nx,0:ny,0:nz) )


do k=0,nz; kk = (ny+1)*k + pointb2(NGcr)
do j=0,ny; n = j + kk
   temp2(1 ,j,k) = bphi2(n,1)
    write(471+NRANK,*) j,k,n, temp2(1 ,j,k),bphi2(n,1)
end do;end do
close(471+NRANK)
call collect( temp1,temp2,n1,n1,n1,nx,ny,nz,0 )
do k=1,n1; kk = n1*(k-1) + pointb1(NGcr)
do j=1,n1; n = j-1 + kk
   bphi1(n,1) = temp1(1 ,j,k)
   write(371+NRANK,*) j,k,n, bphi1(n,1) , temp1(1 ,j,k)
end do;end do
close(371+NRANK)
do k=0,nz; kk = (ny+1)*k + pointb2(NGcr)
do j=0,ny; n = j + kk
   temp2(nx,j,k) = bphi2(n,2)
   write(171+NRANK,*) j,k,n,nx, temp2(nx ,j,k),bphi2(n,2)
end do;end do
call collect( temp1,temp2,n1,n1,n1,nx,ny,nz,0 )
do k=1,n1; kk = n1*(k-1) + pointb1(NGcr)
do j=1,n1; n = j-1 + kk
   bphi1(n,2) = temp1(n1,j,k)
   write(271+NRANK,*) j,k,n,n1, bphi1(n,2) , temp1(n1 ,j,k)
end do;end do
close(171+NRANK)
close(271+NRANK)
DEALLOCATE(temp1,temp2)

open(171+NRANK,file='Phi5'//name//nR//'.dat')
open(271+NRANK,file='Phi6'//name//nR//'.dat')
open(371+NRANK,file='Phi7'//name//nR//'.dat')
open(471+NRANK,file='Phi8'//name//nR//'.dat')

nc = (2**NGcr)+1
w  = 0.125d0
do lc=NGcr-1,1,-1
  nf = nc
  nc = nc/2+1
  do l = 1, 2

  do jc = 2,nc-1; jf = 2*jc-1
  do ic = 2,nc-1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 4.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
     write(171+NRANK,*) lc,kf,kc,jc,ic,kf,kc,l,bphi1(pointb1(lc)+kc,l)
 end do
  end do

  do jc = 2,nc-1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
     write(271+NRANK,*) lc,kf,kc,jc,ic,kf,kc,l,bphi1(pointb1(lc)+kc,l)
 end do
  do ic = 2,nc-1; if = 2*ic-1
    jc = 1; jf = 2*jc-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
    jc = nc; jf = 2*jc-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )
   write(371+NRANK,*) lc,kf,kc,jc,ic,kf,kc,l,bphi1(pointb1(lc)+kc,l)
 end do

    jc = 1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
    jc = 1 ; jf = 2*jc-1
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
    jc = nc; jf = 2*jc-1
    ic = 1 ; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )
    jc = nc; jf = 2*jc-1
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )

    write(471+NRANK,*) lc,kf,kc,jc,ic,kf,kc,l,bphi1(pointb1(lc)+kc,l)
    ! write(*,*) bphi1(pointb1(lc)+kc,l),'===========bphi1============='
  end do
end do

close(171+NRANK)
close(271+NRANK)
close(371+NRANK)
close(471+NRANK)
!do i=1,100
!   write(*,*) bphi1(i,1),bphi1(i,2),'===========bphi1============='
!end do
END SUBROUTINE PB


SUBROUTINE rlft3(data,speq,nn1,nn2,isign) ! numelical in F77
INTEGER isign,nn1,nn2
COMPLEX*16 data(nn1/2,nn2),speq(nn2)
INTEGER i1,i2,j1,j2,nn(2)
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
COMPLEX*16 c1,c2,h1,h2,w
c1=cmplx(0.5d0,0.0d0)
c2=cmplx(0.0d0,-0.5d0*isign)
theta=6.28318530717959d0/dble(isign*nn1) !(yの数) 2*pi/grid(y)
wpr=-2.0d0*dsin(0.5d0*theta)**2
wpi=dsin(theta)
nn(1)=nn1/2
nn(2)=nn2

if(isign.eq.1) then
  call fourn(data,nn,2,isign)
    do i2=1,nn2
      speq(i2)=data(1,i2)
    enddo
endif
wr=1.0d0
wi=0.0d0
do i1=1,nn1/4+1
  j1=nn1/2-i1+2
  do i2=1,nn2
    j2=1
    if(i2.ne.1) j2=nn2-i2+2
    if(i1.eq.1) then
      h1=c1*(data(1,i2)+conjg(speq(j2)))
      h2=c2*(data(1,i2)-conjg(speq(j2)))
      data(1,i2)=h1+h2
      speq(j2)=conjg(h1-h2)
    else
      h1=c1*(data(i1,i2)+conjg(data(j1,j2)))
      h2=c2*(data(i1,i2)-conjg(data(j1,j2)))
      data(i1,i2)=h1+w*h2
      data(j1,j2)=conjg(h1-w*h2)
    endif
  enddo
  wtemp=wr
  wr=wr*wpr-wi*wpi+wr
  wi=wi*wpr+wtemp*wpi+wi
  w=cmplx(wr,wi)
enddo

if(isign.eq.-1) call fourn(data,nn,2,isign)
END SUBROUTINE

SUBROUTINE fourn(data,nn,ndim,isign)
INTEGER isign,ndim,nn(ndim)
DOUBLE PRECISION data(*)
INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2,ip3,k1,k2,n,nprev,nrem,ntot
DOUBLE PRECISION tempi,tempr
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
ntot=1
do idim=1,ndim
  ntot=ntot*nn(idim) !z*y/2
enddo
nprev=1
do idim=1,ndim !2
  n=nn(idim)
  nrem=ntot/(n*nprev)
  ip1=2*nprev
  ip2=ip1*n
  ip3=ip2*nrem
  i2rev=1
  do i2=1,ip2,ip1
    if(i2.lt.i2rev)then
      do i1=i2,i2+ip1-2,2
        do i3=i1,ip3,ip2
          i3rev=i2rev+i3-i2
          tempr=data(i3)
          tempi=data(i3+1)
          data(i3)=data(i3rev)
          data(i3+1)=data(i3rev+1)
          data(i3rev)=tempr
          data(i3rev+1)=tempi
        enddo
      enddo
    endif
    ibit=ip2/2
1   if((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
      i2rev=i2rev-ibit
      ibit=ibit/2
      goto 1
    endif
    i2rev=i2rev+ibit
  enddo
  ifp1=ip1
2 if(ifp1.lt.ip2)then
    ifp2=2*ifp1
    theta=isign*6.28318530717959d0/(ifp2/ip1)
    wpr=-2.d0*sin(0.5d0*theta)**2
    wpi=sin(theta)
    wr=1.d0
    wi=0.d0
    do i3=1,ifp1,ip1
      do i1=i3,i3+ip1-2,2
        do i2=i1,ip3,ifp2
          k1=i2
          k2=k1+ifp1
          tempr=wr*data(k2)-wi*data(k2+1)
          tempi=wr*data(k2+1)+wi*data(k2)
          data(k2)=data(k1)-tempr
          data(k2+1)=data(k1+1)-tempi
          data(k1)=data(k1)+tempr
          data(k1+1)=data(k1+1)+tempi
        enddo
      enddo
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    enddo
    ifp1=ifp2
    goto 2
  endif
  nprev=n*nprev
enddo
END SUBROUTINE

subroutine saveu(Uin,nx,ny,nz,nix,niy,niz,mode)
  USE comvar
  USE mpivar
  USE chmvar
  USE slfgrv
  INCLUDE 'mpif.h'

integer :: nunit=0,st,msig,nx,ny,nz,mode,nix,niy,niz
double precision  :: dt,t(1000)
!double precision,dimension(:,:,:),allocatable :: Uin
character*7 stb(3)
character*3 filenm
CHARACTER*3 NPENUM
double precision,dimension(nix:nx,niy:ny,niz:nz) :: Uin
WRITE(NPENUM,'(I3.3)') NRANK
write(filenm,'(I3.3)') nunit

if (mode==1) then
   !nix=0
   !niy=0
   !niz=0
   !ALLOCATE(Uin(0:nx,0:ny,0:nz))
   !double precision,dimension(nix:nx,niy:ny,niz:nz) :: Uin
   open(10,FILE='/work/maedarn/3DMHD/testrelaxMPI/'//filenm//NPENUM//'.dat',FORM='UNFORMATTED')
end if

!if (mode==2) then
!   nix=1
!   niy=1
!   niz=1
!   ALLOCATE(Uin(nix:nx,niy:ny,niz:nz))
!   open(10,FILE='/work/maedarn/3DMHD/testPhif/'//filenm//NPENUM//'.dat',FORM='UNFORMATTED')
!end if

WRITE(NPENUM,'(I3.3)') NRANK
write(filenm,'(I3.3)') nunit
!open(10,file='/work/maedarn/3DMHD/test/'//filenm//NPENUM//'.dat')
!open(10,FILE='/work/maedarn/3DMHD/testg/'//filenm//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
!100 format(D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3)
!  k=1;j=1;i=1
  do k = niz,nz
  do j = niy,ny
!    write(10) (sngl(U(i,j,k,1)),sngl(U(i,j,k,2)),sngl(U(i,j,k,3)),sngl(U(i,j,k,4)),sngl(U(i,j,k,5)), &
!               sngl(Bcc(i,j,k,1)),sngl(Bcc(i,j,k,2)),sngl(Bcc(i,j,k,3)), &
!               sngl(ndH(i,j,k)),sngl(ndp(i,j,k)),sngl(ndH2(i,j,k)),sngl(ndHe(i,j,k)), &
!               sngl(ndHep(i,j,k)),sngl(ndC(i,j,k)),sngl(ndCO(i,j,k)),sngl(ndCp(i,j,k)), &
!               sngl(Phi(i,j,k)),i=1,Ncellx+1 )
   ! write(10,100) (sngl(U(i,j,k,1)),sngl(U(i,j,k,2)),sngl(U(i,j,k,3)),sngl(U(i,j,k,4)),sngl(U(i,j,k,5)), &
   !            sngl(U(i,j,k,6)),sngl(U(i,j,k,7)),sngl(U(i,j,k,8)), &
   !            sngl(ndH(i,j,k)),sngl(ndp(i,j,k)),sngl(ndH2(i,j,k)),sngl(ndHe(i,j,k)), &
   !            sngl(ndHep(i,j,k)),sngl(ndC(i,j,k)),sngl(ndCO(i,j,k)),sngl(ndCp(i,j,k)), &
     !            sngl(Phi(i,j,k)),i=1,Ncellx+1 )

     write(10,*) (Uin(i,j,k),i=nix,nx)
  end do
  end do

!  101 format(E19.10,E19.10,E19.10,E19.10,E19.10,E19.10,E19.10,E19.10,E19.10)
!  write(10,101) (0.5d0*(x(i)+x(i-1)),U(i,j,k,1),U(i,j,k,2),U(i,j,k,3),U(i,j,k,4),U(i,j,k,5), &
!                 Bcc(i,j,k,1),Bcc(i,j,k,2),Bcc(i,j,k,3), &
!                 ndH(i,j,k),ndp(i,j,k),ndH2(i,j,k),ndHe(i,j,k), &
!                 ndHep(i,j,k),ndC(i,j,k),ndCO(i,j,k),ndCp(i,j,k),i=1,Ncellx )
  close(10)
 ! DEALLOCATE(Uin)
  nunit = nunit + 1


goto 2202
IF(NRANK.EQ.0) THEN
  t(nunit) = time
  open(3,file='/work/maedarn/3DMHD/test/time.DAT')
  do i = 1, nunit
    write(3,'(1p1d25.17)') t(i)
  end do
  close(3)
  write(5,'(1p1d25.17,a,i8,a,1p1e11.3,a)') time,'  Step =',(nunit-1)*nitera,'  dt =', dt,stb(st)
END IF

nunit = nunit + 1

!if(msig.eq.1) then
  IF(NRANK.EQ.0) THEN
    write(5,'(a,1p1e11.3,1p1e11.3)') 'Done ! Time =', time, Tfinal
!    close(5)
    open(2,file='/work/maedarn/3DMHD/test/tsave.DAT')
    write(2,'(1p1d25.17)') time
    write(2,'(i8)') nunit-1
    close(2)
  END IF

  open(8,file='/work/maedarn/3DMHD/test/000'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
  do k = 1, Ncellz+1
  do j = 1, Ncelly+1
    write(8) (U(i,j,k,1),U(i,j,k,2),U(i,j,k,3),U(i,j,k,4),U(i,j,k,5),U(i,j,k,6),U(i,j,k,7),U(i,j,k,8), &
              ndH(i,j,k),ndp(i,j,k),ndH2(i,j,k),ndHe(i,j,k), &
              ndHep(i,j,k),ndC(i,j,k),ndCO(i,j,k),ndCp(i,j,k),Phi(i,j,k), i=1,Ncellx+1 )
  end do
  end do
  close(8)
  !end if
  2202 continue
end subroutine saveu

SUBROUTINE BC_MPI(N_ol,mode)
USE comvar
USE mpivar
INCLUDE 'mpif.h'

integer :: N_ol,mode
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!********************************  BC for the leftsides of domains  *****
  DO K = 1, N_MPI(20)
    CALL MPI_SENDRECV(U(Ncellx+1-N_ol,-1,-1,N_MPI(K)),1,VECU,RIGT,1, &
                      U(       1-N_ol,-1,-1,N_MPI(K)),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    IF((BCx1.eq.4).and.(IST.eq.0)) THEN
      if(N_MPI(K).eq.1) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
        xp = x(0)-0.5d0*dx(0)+IX*dx(0) - time*BBRV(2,1,1)
        nn = int( abs(xp)/(dx(0)*Ncellx*NSPLTx) ) + 1
        xp = xp + dx(0)*Ncellx*NSPLTx * nn
        II = int( xp/dx(0) ) + 1
        U(IX,JY,KZ,1) = DTF(II,JY,KZ)
      END DO;END DO;END DO; goto 1
      end if
      if(N_MPI(K).eq.6) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0 !!! divfree
        U(IX,JY,KZ,6) = U(1,JY,KZ,6)
      END DO;END DO;END DO; goto 1
      end if
      if((N_MPI(K).eq.5).and.(mode.eq.2)) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
        U(IX,JY,KZ,5) = BBRV(5,1,1)*1.5d0 + 0.5d0*( U(IX,JY,KZ,1)*BBRV(2,1,1)**2 &
                      + Bcc(IX,JY,KZ,1)**2+Blg(IX,JY,KZ,1)**2+Blg(IX,JY,KZ,2)**2 )
      END DO;END DO;END DO; goto 1
      end if
      if((N_MPI(K).eq.2).and.(mode.eq.2)) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
        U(IX,JY,KZ,2) = U(IX,JY,KZ,1)*BBRV(2,1,1)
      END DO;END DO;END DO; goto 1
      end if
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
        U(IX,JY,KZ,N_MPI(K)) = BBRV(N_MPI(K),mode,1)
      END DO;END DO;END DO
      1 continue
    END IF
  END DO
!********************************  BC for the rightsides of domains  ****
  DO K = 1, N_MPI(20)
    CALL MPI_SENDRECV(U(1            ,-1,-1,N_MPI(K)),1,VECU,LEFT,1, &
                      U(Ncellx+1     ,-1,-1,N_MPI(K)),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    IF((BCx2.eq.4).and.(IST.eq.NSPLTx-1)) THEN
      if(N_MPI(K).eq.1) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
        xp = 0.5d0*(x(IX-1)+x(IX)) - time*BBRV(2,1,2)
        nn = int( abs(xp)/(dx(0)*Ncellx*NSPLTx) )
        xp = xp - dx(0)*Ncellx*NSPLTx * nn
        II = int( xp/dx(0) ) + 1
        U(IX,JY,KZ,1) = DTF(II,JY,KZ)
      END DO;END DO;END DO; goto 2
      end if
      if(N_MPI(K).eq.6) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+2, Ncellx+N_ol !!! divfree
        U(IX,JY,KZ,6) = U(Ncellx+1,JY,KZ,6)
      END DO;END DO;END DO; goto 2
      end if
      if((N_MPI(K).eq.5).and.(mode.eq.2)) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
        U(IX,JY,KZ,5) = BBRV(5,1,2)*1.5d0 + 0.5d0*( U(IX,JY,KZ,1)*BBRV(2,1,2)**2 &
                      + Bcc(IX,JY,KZ,1)**2+Blg(IX,JY,KZ,1)**2+Blg(IX,JY,KZ,2)**2 )
      END DO;END DO;END DO; goto 2
      end if
      if((N_MPI(K).eq.2).and.(mode.eq.2)) then
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
        U(IX,JY,KZ,2) = U(IX,JY,KZ,1)*BBRV(2,1,2)
      END DO;END DO;END DO; goto 2
      end if
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
        U(IX,JY,KZ,N_MPI(K)) = BBRV(N_MPI(K),mode,2)
      END DO;END DO;END DO
      2 continue
    END IF
  END DO
!************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF


IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!*************************************  BC for the downsides of domains  ****
  DO K = 1, N_MPI(20)
    CALL MPI_SENDRECV(U(-1,Ncelly+1-N_ol,-1,N_MPI(K)),1,VECU,TOP ,1, &
                      U(-1,       1-N_ol,-1,N_MPI(K)),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END DO
!**************************************  BC for the upsides of domains  ****
  DO K = 1, N_MPI(20)
    CALL MPI_SENDRECV(U(-1,1            ,-1,N_MPI(K)),1,VECU,BOTM,1, &
                      U(-1,Ncelly+1     ,-1,N_MPI(K)),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!*************************************  BC for the downsides of domains  ****
  DO K = 1, N_MPI(20)
    CALL MPI_SENDRECV(U(-1,-1,Ncellz+1-N_ol,N_MPI(K)),1,VECU,UP  ,1, &
                      U(-1,-1,       1-N_ol,N_MPI(K)),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END DO
!**************************************  BC for the upsides of domains  ****
  DO K = 1, N_MPI(20)
    CALL MPI_SENDRECV(U(-1,-1,1            ,N_MPI(K)),1,VECU,DOWN,1, &
                      U(-1,-1,Ncellz+1     ,N_MPI(K)),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF

END SUBROUTINE

subroutine saveu1(Uin1,Uin,nx1,ny1,nz1,nx2,ny2,nz2)
 ! USE comvar
  USE mpivar
 ! USE chmvar
 ! USE slfgrv
  INCLUDE 'mpif.h'

integer :: nunit=0,st,msig,nx1,ny1,nz1,mode,nx2,ny2,nz2
!double precision  :: dt,t(1000)
!double precision,dimension(:,:,:),allocatable :: Uin
character*7 stb(3)
character*3 filenm
CHARACTER*3 NPENUM
double precision :: Uin(0:nx2,0:ny2,0:nz2),Uin1(nx1,ny1,nz1)
double precision :: tMPI(0:nx2,0:ny2,0:nz2,0:NPE-1)
WRITE(NPENUM,'(I3.3)') NRANK
write(filenm,'(I3.3)') nunit
write(*,*) '==================saveU============'

open(10+NRANK,FILE='/work/maedarn/3DMHD/test/saveu'//filenm//NPENUM//'.dat')
!100 format(D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3)
!  k=1;j=1;i=1
  do k = 0,nz2
     do j = 0,ny2
        do i=0,nx2
     write(10+NRANK,*) i,j,k,Uin(i,j,k)
  end do
end do
end do

close(10+NRANK)

do k=0,nz2; do j=0,ny2; do i=0,nx2
  tMPI(i,j,k,NRANK)=Uin(i,j,k)
end do;end do;end do
 ! DEALLOCATE(Uin)
  nunit = nunit + 1
!**********NEW FORM***********
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!**********NEW FORM***********

do Nroot=0,NPE-1
  CALL MPI_BCAST(tMPI(0,0,0,Nroot),(nx2+1)*(ny2+1)*(nz2+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

!**********NEW FORM***********
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!**********NEW FORM***********

do Nroot=0,NPE-1
 ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
 nxed=nx2-1; IF(ISTt.eq.NSPLTx-1) nxed=nx2
 nyed=ny2-1; IF(JSTt.eq.NSPLTy-1) nyed=ny2
 nzed=nz2-1; IF(KSTt.eq.NSPLTz-1) nzed=nz2
 do kk=1,nzed;k=KSTt*(nz2-1)+kk
  do jj=1,nyed;j=JSTt*(ny2-1)+jj
   do ii=1,nxed;i=ISTt*(nx2-1)+ii
    Uin1(i,j,k) = tMPI(ii,jj,kk,Nroot)
 end do;end do;end do;end do

 !write (name,'(i4.4)') countin
 !write (nR,'(i3.3)') NRANK
 !goto 1202
 open(251+NRANK,file='saveuaf'//filenm//NPENUM//'.dat')

 do k=1,nz1; do j=1,ny1; do i=1,nx1
         write(251+NRANK,*) nx1, i,j,k, Uin1(i,j,k)
      end do
   end do
end do
 !do kk=1,nzed;k=KSTt*(nz2-1)+kk
 ! do jj=1,nyed;j=JSTt*(ny2-1)+jj
 !  do ii=1,nxed;i=ISTt*(nx2-1)+ii
 !         write(251+NRANK,*) ii,jj,kk, u1(i,j,k)!,tMPI(ii,jj,kk,Nroot)
 !      end do
 !   end do
 !end do
 close(251+NRANK)
end subroutine saveu1
