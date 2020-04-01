MODULE comvar
!INTEGER, parameter :: ndx=130, ndy=130, ndz=130, ndmax=130, Dim=3 !1024^3
!INTEGER, parameter :: ndx=66, ndy=66, ndz=66, ndmax=66, Dim=3 !512^3
!INTEGER, parameter :: ndx=34, ndy=34, ndz=34, ndmax=34, Dim=3
INTEGER, parameter :: ndx=18, ndy=18, ndz=18, ndmax=18, Dim=3
DOUBLE PRECISION, dimension(-1:ndx) :: x,dx
DOUBLE PRECISION, dimension(-1:ndy) :: y,dy
DOUBLE PRECISION, dimension(-1:ndz) :: z,dz
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U, Bcc, Blg, Vfc, EMF
DOUBLE PRECISION, dimension(:,:,:),   allocatable :: dnc, xlag, dxlagM

DOUBLE PRECISION, parameter :: kb=8.63359d0, Kcond=1.6384d-2
DOUBLE PRECISION  :: gamma,gammi1,gammi2,gammi3,gampl1,gampl2,gampl3
DOUBLE PRECISION  :: CFL,facdep,tfinal,time,phr(-1:400)
DOUBLE PRECISION  :: pmin,pmax,rmin,rmax
double precision :: shusoku1=0.0d0 , phiratio=1.0d0/3.0d0,rhomean
INTEGER :: Ncellx,Ncelly,Ncellz,iwx,iwy,iwz,maxstp,nitera
INTEGER :: ifchem,ifthrm,ifrad,ifgrv

!DOUBLE PRECISION :: cg=1.0d0,sourratio=0.5d0
DOUBLE PRECISION, parameter :: cg=1.0d0,sourratio=0.5d0
double precision :: deltalength
!integer ifevogrv,ifevogrv2
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
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Phi ! , Phiexa
double precision, dimension(:,:,:), allocatable :: Phidt! , Phicgp , Phicgm
DOUBLE PRECISION :: Lbox
!double precision :: deltalength , cgcsratio= 1.0d0,cgratio1=0.2d0 !, shusoku1=0.0d0
double precision ::  cgcsratio= 1.0d0,cgratio1=0.2d0 !, shusoku1=0.0d0

DOUBLE PRECISION , dimension(:,:,:,:), allocatable ::  Phicgm , Phi1step , Phi2step , Phicgp
DOUBLE PRECISION , dimension(:,:,:), allocatable ::  Phigrd , Phiexa

INTEGER :: pointb1(0:15),pointb2(0:15)
DOUBLE PRECISION, dimension(:,:,:), allocatable :: bphi1l,bphi2l,bstep1l,bstep2l &
     !,bphi1r,bphi2r,bstep1r,bstep2r
     ,bphil,bphir
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: bphixl,bphixr,bphiyl,bphiyr,bphizl,bphizr,&
     bstepxl,bstepxr,bstepyl,bstepyr,bstepzl,bstepzr,source,sourcedt,sourcedt2!,bphil,bphir
integer , parameter :: bnd=3,loopbc=3
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

   !write(*,*) 'OK1',NRANK
   !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   !write(*,*) 'OK2',NRANK

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

!*********grvwave*********
ALLOCATE(Phiexa(-1:ndx,-1:ndy,-1:ndz))
ALLOCATE(Phigrd(-1:ndx,-1:ndy,-1:ndz))
!ALLOCATE(Phicgp(-1:ndx,-1:ndy,-1:ndz))
!ALLOCATE(Phicgm(-1:ndx,-1:ndy,-1:ndz))
!ALLOCATE(Phi1step(-1:ndx,-1:ndy,-1:ndz))
!ALLOCATE(Phi2step(-1:ndx,-1:ndy,-1:ndz))
ALLOCATE(Phicgp(-1:ndx,-1:ndy,-1:ndz,1:Dim))
ALLOCATE(Phicgm(-1:ndx,-1:ndy,-1:ndz,1:Dim))
ALLOCATE(Phi1step(-1:ndx,-1:ndy,-1:ndz,1:Dim))
ALLOCATE(Phi2step(-1:ndx,-1:ndy,-1:ndz,1:Dim))
ALLOCATE(sourcedt(-1:ndx,-1:ndy,-1:ndz,1:Dim))
ALLOCATE(sourcedt2(-1:ndx,-1:ndy,-1:ndz,1:Dim))
ALLOCATE(source(-1:ndx,-1:ndy,-1:ndz,1:Dim))
!ALLOCATE(bphi1l(1:ndy-2,1:ndz-2,-1:1))
!ALLOCATE(bphi2l(1:ndy-2,1:ndz-2,-1:1))
!ALLOCATE(bstep1l(1:ndy-2,1:ndz-2,-1:1))
!ALLOCATE(bstep2l(1:ndy-2,1:ndz-2,-1:1))
!ALLOCATE(bphi1r(1:ndy-2,1:ndz-2,ndx-2:ndx))
!ALLOCATE(bphi2r(1:ndy-2,1:ndz-2,ndx-2:ndx))
!ALLOCATE(bstep1r(1:ndy-2,1:ndz-2,ndx-2:ndx))
!ALLOCATE(bstep2r(1:ndy-2,1:ndz-2,ndx-2:ndx))

ALLOCATE(bphil(-1:ndy,-1:ndz,-1:1     ))
ALLOCATE(bphir(-1:ndy,-1:ndz,ndx-2:ndx))

!ALLOCATE(bphi1l(-1:ndy,-1:ndz,-1:1,1:Dim))
!ALLOCATE(bphi2l(-1:ndy,-1:ndz,-1:1,1:Dim))
!ALLOCATE(bstep1l(-1:ndy,-1:ndz,-1:1,1:Dim))
!ALLOCATE(bstep2l(-1:ndy,-1:ndz,-1:1,1:Dim))
!ALLOCATE(bphi1r(-1:ndy,-1:ndz,ndx-2:ndx,1:Dim))
!ALLOCATE(bphi2r(-1:ndy,-1:ndz,ndx-2:ndx,1:Dim))
!ALLOCATE(bstep1r(-1:ndy,-1:ndz,ndx-2:ndx,1:Dim))
!ALLOCATE(bstep2r(-1:ndy,-1:ndz,ndx-2:ndx,1:Dim))

ALLOCATE( bphixl(-1:ndy,-1:ndz,-1:1     ,1:Dim))
ALLOCATE(bstepxl(-1:ndy,-1:ndz,-1:1     ,1:Dim))
ALLOCATE( bphixr(-1:ndy,-1:ndz,ndx-2:ndx,1:Dim))
ALLOCATE(bstepxr(-1:ndy,-1:ndz,ndx-2:ndx,1:Dim))

ALLOCATE( bphiyl(-1:ndz,-1:ndx,-1:1     ,1:Dim))
ALLOCATE(bstepyl(-1:ndz,-1:ndx,-1:1     ,1:Dim))
ALLOCATE( bphiyr(-1:ndz,-1:ndx,ndy-2:ndy,1:Dim))
ALLOCATE(bstepyr(-1:ndz,-1:ndx,ndy-2:ndy,1:Dim))

ALLOCATE( bphizl(-1:ndx,-1:ndy,-1:1     ,1:Dim))
ALLOCATE(bstepzl(-1:ndx,-1:ndy,-1:1     ,1:Dim))
ALLOCATE( bphizr(-1:ndx,-1:ndy,ndz-2:ndz,1:Dim))
ALLOCATE(bstepzr(-1:ndx,-1:ndy,ndz-2:ndz,1:Dim))
!*********grvwave*********

!write(*,*) 'OK3'

call INITIA
!write(*,*) 'OK'
call EVOLVE

!write(*,*) 'OK'

DEALLOCATE(U)
DEALLOCATE(ndH,ndp,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)
DEALLOCATE(DTF)
DEALLOCATE(Phi)

!********gravwave**********
DEALLOCATE(Phiexa,Phigrd)
DEALLOCATE(Phicgp)
DEALLOCATE(Phicgm)
DEALLOCATE(Phi1step)
DEALLOCATE(Phi2step)
!DEALLOCATE(bphi1l,bphi2l,bstep1l,bstep2l)
!DEALLOCATE(bphi1r,bphi2r,bstep1r,bstep2r)
DEALLOCATE(bphil,bphir)
DEALLOCATE(bphixl,bphixr,bphiyl,bphiyr,bphizl,bphizr,bstepxl,bstepxr,bstepyl,bstepyr,bstepzl,bstepzr)
DEALLOCATE(source,sourcedt,sourcedt2)
!********gravwave**********

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
character*3 :: NPENUM,MPIname
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision, dimension(:,:), allocatable :: plane,rand
integer i3,i4,i2y,i2z,rsph2,pls
double precision cenx,ceny,cenz,rsph,rrsph,Hsheet,censh,minexa

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

!dx=dy=dz
deltalength = dx_i(0)
!dx=dy=dz

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

write(*,*) 'check-dx' ,dx(1),deltalength

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


  !write(*,*) 'sheet1'
   !********************sheet***********************
  !goto 6011
  DTF(:,:,:) = 0.0d0
  !dinit1=1.0d0/G4pi
  dinit1 = 2.0d0/G4pi/90.d0
  censh = ql1x + dx(1)/2.0d0 !x=serfase
  Hsheet = 1.0d1
  !rsph = ql1x-ql1x/5.0d0
  !rsph2=int(dble(Np1x)*0.8d0)
  !Hsheet = dble(Np1x) / 5.0d0
  rhomean = dinit1*(Hsheet*2.d0)/(ql1x+ql2x)
  do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   !i2 = IST*Ncellx+i
   !i2y = JST*Ncelly+j
   !i2z = KST*Ncellz+k

   !rsph=dsqrt( (cenx-dble(i2))**2 + (ceny-dble(i2y))**2 + (cenz-dble(i2z))**2 )
   if( dabs(x(i) - censh ) .le. Hsheet ) then
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
      !ndH(i,j,k)   = 0.0d0
      !ndp(i,j,k)   = 0.0d0
      !ndH2(i,j,k)  = 0.0d0
      !ndHe(i,j,k)  = 0.0d0
      !ndHep(i,j,k) = 0.0d0
      !ndC(i,j,k)   = 0.0d0
      !ndCO(i,j,k)  = 0.0d0
      !ndCp(i,j,k)  = 0.0d0
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
end do
end do
end do

!minexa = 1.0d2
write(MPIname,'(i3.3)') NRANK
open(142,file='/work/maedarn/3DMHD/test/phiexact'//MPIname//'.DAT')
saexact1 = G4pi * dinit1 * Hsheet * dabs( x_i(0) - censh)  - G4pi/2.0d0 * dinit1 * Hsheet**2
!saexact1 = G4pi * dinit1 * Hsheet * dabs( x_i(0) - ql1x)  - G4pi/2.0d0 * dinit1 * Hsheet**2


!write(*,*) x_i(0) - censh , x_i(0),saexact1 ,'bcx'
!ALLOCATE(Phiexa(-1:ndx,-1:ndy,-1:ndz))
do k=-1,Ncellz+2
   do j=-1,Ncelly+2
      do i= -1,Ncellx+2
         if( dabs(x(i) - censh ) .le. Hsheet ) then
            Phiexa(i,j,k) = G4pi/2.0d0 * dinit1 * (x(i) - censh )**2
            write(142,*) sngl(G4pi/2.0d0 * dinit1 * (x(i) - censh )**2)
         else
            Phiexa(i,j,k) = G4pi * dinit1 * Hsheet * dabs(x(i) - censh)  - G4pi/2.0d0 * dinit1 * Hsheet**2
            write(142,*) sngl(G4pi * dinit1 * Hsheet * dabs(x(i) - censh)  - G4pi/2.0d0 * dinit1 * Hsheet**2)
         end if
!         write(142,*) sngl(Phiexa(i,j,k))
         !minexa=dmin1(minexa,Phiexa(i,j,k))
      end do
   end do
end do

!do k=-1,Ncellz+2
!   do j=-1,Ncelly+2
!      do pls = 0,2,1
!         bphi1l(j,k,1-abs(pls)) = Phiexa(1-abs(pls),j,k)
!         bphi1r(j,k,Ncellx+abs(pls)) = Phiexa(Ncellx+abs(pls),j,k)
!         bphi2l(j,k,1-abs(pls)) = Phiexa(1-abs(pls),j,k)
!         bphi2r(j,k,Ncellx+abs(pls)) = Phiexa(Ncellx+abs(pls),j,k)
!      end do
!   end do
!end do


do k=-1,ndz
do j=-1,ndy
do i=0,ndx-1
   Phigrd(i,j,k)=(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx(1)
   !write(144,*) sngl(x(i)) , Phigrd(i) , Phiexa(i-1),Phiexa(i+1)
end do
Phigrd(-1,j,k)=(-Phiexa(0,j,k)+Phiexa(1,j,k))/dx(1)
Phigrd(ndx,j,k)=(Phiexa(ndx-1,j,k)-Phiexa(ndx-2,j,k))/dx(1)
end do
end do

write(*,*) NRANK,Phigrd(0,1,1),Phigrd(Ncellx-1,1,1),dinit1
!DEALLOCATE(Phiexa)
close(142)

dinit1=0.0d0
 !6011 continue
!********************sheet***********************
!write(*,*) 'sheet2'

  !********************sphere***********************
  goto 6001
  DTF(:,:,:) = 0.0d0
  !dinit1=1.0d0/G4pi
  dinit1=3.0d0/G4pi/4.d1/4.d1
  cenx=dble(Np1x)+0.5d0
  ceny=dble(Np1y)+0.5d0
  cenz=dble(Np1z)+0.5d0
  !rsph = ql1x-ql1x/5.0d0
  !rsph2=int(dble(Np1x)*0.8d0)
  rrsph = dble(Np1x)*0.2d0
  do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   i2 = IST*Ncellx+i
   i2y = JST*Ncelly+j
   i2z = KST*Ncellz+k
   cenx=dble(Np1x)+0.5d0
   ceny=dble(Np1y)+0.5d0
   cenz=dble(Np1z)+0.5d0
   rsph=dsqrt( (cenx-dble(i2))**2 + (ceny-dble(i2y))**2 + (cenz-dble(i2z))**2 )
   if(rsph .le. rrsph ) then
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

rhomean = dinit1*G4pi/3.d0 * (rrsph)**3 / ((ql1x+ql2x)**3)
dinit1=0.0d0
 6001 continue
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
  DTF(:,:,:) = dinit1
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

call CC(4,0.d0)

do k=1,Ncellz+1; do j=1,Ncelly+1; do i=1,Ncellx+1
  nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
  ndtot(i,j,k) = ndp(i,j,k)+ndH(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
  Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; NCO(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
  Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; NCO(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
end do; end do; end do

if(ifrad.eq.2) then; do l=1,20; call SHIELD(); end do; end if

   !write(*,*) 'pregrv1'
if(ifgrv.eq.2) then
  N_MPI(20)=1; N_MPI(1)=1; iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)
  Lbox=ql1x+ql2x!; call GRAVTY(0.d0,1); call GRAVTY(0.d0,2)
  !call SELFGRAVWAVE(0.0d0,1)
  call SELFGRAVWAVE(0.0d0,0) !密度場の生成の時
  !call SELFGRAVWAVE(0.0d0,6) !calculate cg
end if
 !write(*,*) 'posgrv1'
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

!---------------debug-------------------
!write(*,*) '-------------3-----------',NRANK
!---------------debug-------------------

do in10 = 1, maxstp

  time_CPU(1) = MPI_WTIME()
  tsave = dtsave * dble(itime)
  if(time.ge.tfinal) goto 9000
  if(time.ge.tsave ) goto 7777
  !call SAVEU(nunit,dt,stb,st,t,0)

  do in20 = 1, nitera
     !write(*,*) '---top---'
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),Bcc(1,1,1,2),U(1,1,1,7),'point'
    !tsave2D = dtsave2D * nunit2D
    !if(time.ge.tsave2D) call SAVEU2D(nunit2D)
    !if(time.ge.tfinal) goto 9000
    !if(time.ge.tsave ) goto 7777
!***** Determine time-step dt *****
    !dt_mpi(NRANK) = tfinal
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point1'
    !call Couran(tLMT)
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),tLMT,'point1'
    !dt_mpi(NRANK) = dmin1( dt_mpi(NRANK), CFL * tLMT )
    !st_mpi(NRANK) = 1
    !stt= dt_mpi(NRANK)

    !---------------debug-------------------
    write(*,*) '-------------4-----------',NRANK,in20,in10
    !---------------debug-------------------


    !if(ifgrv==2) then
    !call SELFGRAVWAVE(stt,5)
    !end if


    !---------------debug-------------------
    !write(*,*) '-------------5-----------',NRANK,in20,in10
    !---------------debug-------------------



    !--------for INIT---------------
    !call Stblty(tLMT)
    !--------for INIT---------------



    !---------------debug-------------------
    !write(*,*) '-------------6-----------',tLMT,NRANK,in20,in10
    !---------------debug-------------------
    dt=dx(1)/cg*0.1d0
    !if(ifgrv==2) then
    !   call SELFGRAVWAVE(tLMT,7)

       !---------------debug-------------------
       !write(*,*) '-------------99-----------',NRANK,in20,in10
       !---------------debug-------------------


     !  call SELFGRAVWAVE(tLMT,11)

       !---------------debug-------------------
       !write(*,*) '-------------88-----------',NRANK
       !---------------debug-------------------
    !endif

 !---------------debug-------------------
 !write(*,*) '-------------1-----------',NRANK
 !---------------debug-------------------
goto 342
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),tLMT,'point2'
    dt_mpi(NRANK) = dmin1( dt_mpi(NRANK), tLMT    )
    if(dt_mpi(NRANK).lt.stt) st_mpi(NRANK) = 2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for MPI
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    CALL MPI_GATHER(dt_mpi(NRANK),1,MPI_REAL8,   &
                    dt_gat       ,1,MPI_REAL8,   &
                    0            ,MPI_COMM_WORLD,IERR)
    CALL MPI_GATHER(st_mpi(NRANK),1,MPI_INTEGER, &
                    st_gat       ,1,MPI_INTEGER, &
                    0            ,MPI_COMM_WORLD,IERR)
    IF(NRANK.EQ.0)  THEN
      dt  = tfinal
      dtt = tfinal
      do i_t = 0, NPE-1
         dt  = dmin1( dt, dt_gat(i_t) )
         !write(*,*) '--------------dt--------------' , dt
        if(dt.lt.dtt) st = st_gat(i_t)
        dtt = dt
      end do
   END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
    CALL MPI_BCAST(dt,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    if((mod(in20,10).eq.1).and.(NRANK.eq.0)) write(*,*) in20,time,dt
    !if(NRANK.eq.0) write(*,*) in20,time,dt
    if(time+dt.gt.tfinal) dt = tfinal - time
    if(time+dt.gt.tsave ) dt = tsave  - time

    !write(*,*) '-------------dt--------------------cl-------------' , dt,NRANK
!if(NRANK==40) write(*,*) NRANK,in20,dt,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point3'
!***** Source parts 1*****
    !if(ifgrv.eq.2) then
       !call GRAVTY(dt,3)
       !call SELFGRAVWAVE(dt,3)
       !call SELFGRAVWAVE(dt*0.5d0,3)
    !end if
    !call SOURCE(0.5d0*dt)
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point4'
    !***** Godunov parts *****

342 continue


    !---------------------------skip-----------------------------
    !---------------------------skip-----------------------------
    goto 263
    !---------------------------skip-----------------------------
    !---------------------------skip-----------------------------

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
!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point5'
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

    !---------------------------skip-----------------------------
    !---------------------------skip-----------------------------
    263 continue
    !---------------------------skip-----------------------------
    !---------------------------skip-----------------------------

    !---------------debug-------------------
    !write(*,*) '-------------10-----------',NRANK
    !---------------debug-------------------

!if(NRANK==40) write(*,*) NRANK,in20,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point6'
!***** Source parts 2*****
    !call SOURCE(0.5d0*dt)
    if(ifgrv.eq.2) then
       !call GRAVTY(dt,2)
       !call GRAVTY(dt,3)
       !call GRAVTY(dt*0.5d0,3)
       !write(*,*) dt,NRANK,'--dt--dt--'
       call SELFGRAVWAVE(dt,2)
       !---debug---
       !call  SELFGRAVWAVE(0.0d0,4)
       !call SELFGRAVWAVE(dt*0.5d0,3)

       !---------------debug-------------------
       !write(*,*) '-------------13-----------',NRANK, nitera ,maxstp
       !---------------debug-------------------


    end if

    !*********************************!収束判定
    !call SELFGRAVWAVE(0.0d0,8) !収束判定
    !if(shusoku1 > 1.0d0) then
    !   !---------------debug-------------------
    !   write(*,*) '-------------goto-----------',NRANK
       !---------------debug-------------------
    !   goto 2419
    !end if
    !*********************************!収束判定

    !---------------debug-------------------
    !write(*,*) '-------------12-----------',NRANK
    !---------------debug-------------------
    !call DISSIP()
    !time = time + dt
 end do
  !itime = itime - 1
  7777   continue
  !itime = itime + 1

  if(ifgrv.eq.2) then
     call  SELFGRAVWAVE(0.0d0,4)
  end if


  !*********************************!収束判定
  !call SELFGRAVWAVE(0.0d0,8) !収束判定
  !if(shusoku1 > 1.0d0) then
     !---------------debug-------------------
  !   write(*,*) '-------------goto-----------',NRANK
     !---------------debug-------------------
  !  goto 2419
  !end if
  !*********************************!収束判定



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for MPI
  IF(NRANK.EQ.0) THEN
    time_CPU(2) = MPI_WTIME()
    time_CPU(2) = ( time_CPU(2)-time_CPU(1) )/3.6d3
    time_CPU(3) = time_CPU(3)+time_CPU(2)
    IF(time_CPU(3)+time_CPU(2).GT.11.7d0) Time_signal=1
 END IF
  CALL MPI_BCAST(Time_signal,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  IF(Time_signal.EQ.1) GOTO 9000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end do

!---------------debug-------------------
!write(*,*) '-------------22-----------',NRANK
!---------------debug-------------------
!**********************!収束判定
!if(ifgrv.eq.2) then !if文中には飛べない
2419 continue
if(ifgrv.eq.2) then
   call  SELFGRAVWAVE(0.0d0,4)
end if
!**********************!収束判定
write(*,*) 'save',NRANK

9000 continue
IF(NRANK.EQ.0) write(*,*) 'MPI time1 = ',MPI_WTIME()
!call SAVEU(nunit,dt,stb,st,t,1)



END SUBROUTINE EVOLVE


!----------------------------------------------------------- SAVE VARIABLES ---!

SUBROUTINE SAVEU(nunit,dt,stb,st,t,msig)
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'

integer :: nunit,st,msig
double precision  :: dt,t(1000)
character*7 stb(3)
character*3 filenm
CHARACTER*3 NPENUM


WRITE(NPENUM,'(I3.3)') NRANK
write(filenm,'(I3.3)') nunit
!open(10,file='/work/maedarn/3DMHD/test/'//filenm//NPENUM//'.dat')
open(10,FILE='/work/maedarn/3DMHD/test/'//filenm//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
100 format(D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3,D10.3)
  k=1;j=1;i=1
  do k = 1, Ncellz+1
  do j = 1, Ncelly+1
!    write(10) (sngl(U(i,j,k,1)),sngl(U(i,j,k,2)),sngl(U(i,j,k,3)),sngl(U(i,j,k,4)),sngl(U(i,j,k,5)), &
!               sngl(Bcc(i,j,k,1)),sngl(Bcc(i,j,k,2)),sngl(Bcc(i,j,k,3)), &
!               sngl(ndH(i,j,k)),sngl(ndp(i,j,k)),sngl(ndH2(i,j,k)),sngl(ndHe(i,j,k)), &
!               sngl(ndHep(i,j,k)),sngl(ndC(i,j,k)),sngl(ndCO(i,j,k)),sngl(ndCp(i,j,k)), &
!               sngl(Phi(i,j,k)),i=1,Ncellx+1 )
    write(10) (sngl(U(i,j,k,1)),sngl(U(i,j,k,2)),sngl(U(i,j,k,3)),sngl(U(i,j,k,4)),sngl(U(i,j,k,5)), &
               sngl(U(i,j,k,6)),sngl(U(i,j,k,7)),sngl(U(i,j,k,8)), &
               sngl(ndH(i,j,k)),sngl(ndp(i,j,k)),sngl(ndH2(i,j,k)),sngl(ndHe(i,j,k)), &
               sngl(ndHep(i,j,k)),sngl(ndC(i,j,k)),sngl(ndCO(i,j,k)),sngl(ndCp(i,j,k)), &
               sngl(Phi(i,j,k)),i=1,Ncellx+1 )
  end do
  end do

!  101 format(E19.10,E19.10,E19.10,E19.10,E19.10,E19.10,E19.10,E19.10,E19.10)
!  write(10,101) (0.5d0*(x(i)+x(i-1)),U(i,j,k,1),U(i,j,k,2),U(i,j,k,3),U(i,j,k,4),U(i,j,k,5), &
!                 Bcc(i,j,k,1),Bcc(i,j,k,2),Bcc(i,j,k,3), &
!                 ndH(i,j,k),ndp(i,j,k),ndH2(i,j,k),ndHe(i,j,k), &
!                 ndHep(i,j,k),ndC(i,j,k),ndCO(i,j,k),ndCp(i,j,k),i=1,Ncellx )
close(10)


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

END SUBROUTINE SAVEU


SUBROUTINE SAVEU2D(nunit2D)
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'

character*3 filenm
CHARACTER*3 NPENUM
!write(*,*) Ncellx,Ncelly
WRITE(NPENUM,'(I3.3)') NRANK
write(filenm,'(I3.3)') nunit2D
open(11,FILE='/work/maedarn/3DMHD/test/2D'//filenm//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
  k=1; do j=1,Ncelly
    write(11) (sngl(U(i,j,k,1)),sngl(U(i,j,k,2)),sngl(U(i,j,k,3)),sngl(U(i,j,k,4)),sngl(U(i,j,k,5)), &
               sngl(Bcc(i,j,k,1)),sngl(Bcc(i,j,k,2)),sngl(Bcc(i,j,k,3)), &
               sngl(ndH(i,j,k)),sngl(ndp(i,j,k)),sngl(ndH2(i,j,k)),sngl(ndHe(i,j,k)), &
               sngl(ndHep(i,j,k)),sngl(ndC(i,j,k)),sngl(ndCO(i,j,k)),sngl(ndCp(i,j,k)), &
               sngl(Phi(i,j,k)), i=1,Ncellx )
  end do
close(11)

nunit2D = nunit2D + 1

IF(NRANK.EQ.0) THEN
  open(3,file='/work/maedarn/3DMHD/test/tsave2D.DAT')
    write(3,*) nunit2D
    write(3,*) time
  close(3)
END IF

END SUBROUTINE SAVEU2D


!=====================================================================*
!           Boundary condition                                        *
!=====================================================================*
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
!---------------------------------------------------------------------------
!----------------------------------------- Boundary Cond. for CelCenMagFie -

SUBROUTINE BC_MPI_OT(N_ol,mode)
USE comvar
USE mpivar
USE chmvar
INCLUDE 'mpif.h'

integer :: N_ol,mode,NV(2)
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

klps = -1

IF(iwx.eq.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!***********************************  BC for Lagrange coordinate  **********
  IF(N_MPI(10).eq.1) THEN
    CALL MPI_SENDRECV(dxlagM(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
                      dxlagM(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    IF((BCx1.eq.4).and.(IST.eq.0)) THEN
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
        dxlagM(IX,JY,KZ) = dx(0)
      END DO;END DO;END DO
    END IF

    CALL MPI_SENDRECV(dxlagM(1            ,-1,-1),1,VECU,LEFT,1, &
                      dxlagM(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    IF((BCx2.eq.4).and.(IST.eq.NSPLTx-1)) THEN
      DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
        dxlagM(IX,JY,KZ) = dx(Ncellx+1)
      END DO;END DO;END DO
    END IF
  END IF
!********************************************************  BC for Bcc ******
  IF(N_MPI(11).eq.1) THEN
    DO K = 1, 3
      CALL MPI_SENDRECV(Bcc(Ncellx+1-N_ol,-1,-1,K),1,VECU,RIGT,1, &
                        Bcc(       1-N_ol,-1,-1,K),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
      IF((BCx1.eq.4).and.(IST.eq.0)) THEN
        IF(K.eq.1) THEN
          DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
            Bcc(IX,JY,KZ,K) = Bcc(1,JY,KZ,K)  !!! divfree
          END DO;END DO;END DO
        ELSE
          DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
            Bcc(IX,JY,KZ,K) = BBRV(5+K,1,1)
          END DO;END DO;END DO
        END IF
      END IF

      CALL MPI_SENDRECV(Bcc(1            ,-1,-1,K),1,VECU,LEFT,1, &
                        Bcc(Ncellx+1     ,-1,-1,K),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
      IF((BCx2.eq.4).and.(IST.eq.NSPLTx-1)) THEN
        IF(K.eq.1) THEN
          DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
            Bcc(IX,JY,KZ,K) = Bcc(Ncellx,JY,KZ,K)  !!! divfree
          END DO;END DO;END DO
        ELSE
          DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
            Bcc(IX,JY,KZ,K) = BBRV(5+K,1,2)
          END DO;END DO;END DO
        END IF
      END IF
    END DO
  END IF
!********************************************************  BC for Blg ******
  IF(N_MPI(12).eq.1) THEN
    DO K = 1, 2
      CALL MPI_SENDRECV(Blg(Ncellx+1-N_ol,-1,-1,K),1,VECU,RIGT,1, &
                        Blg(       1-N_ol,-1,-1,K),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
      IF((BCx1.eq.4).and.(IST.eq.0)) THEN
        DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
          Blg(IX,JY,KZ,K) = BBRV(6+K,1,1)
        END DO;END DO;END DO
      END IF

      CALL MPI_SENDRECV(Blg(1            ,-1,-1,K),1,VECU,LEFT,1, &
                        Blg(Ncellx+1     ,-1,-1,K),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
      IF((BCx2.eq.4).and.(IST.eq.NSPLTx-1)) THEN
        DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
          Blg(IX,JY,KZ,K) = BBRV(6+K,1,2)
        END DO;END DO;END DO
      END IF
    END DO
  END IF
!********************************************************  BC for Vfc ******
  IF(N_MPI(13).eq.1) THEN
    NV(1) = 2; NV(2) = 3
    DO K = 1, 2
      CALL MPI_SENDRECV(Vfc(Ncellx+1-N_ol,-1,-1,NV(K)),1,VECU,RIGT,1, &
                        Vfc(       1-N_ol,-1,-1,NV(K)),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
      IF((BCx1.eq.4).and.(IST.eq.0)) THEN
        DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
          Vfc(IX,JY,KZ,NV(K)) = Vfc(1,JY,KZ,NV(K))
        END DO;END DO;END DO
      END IF

      CALL MPI_SENDRECV(Vfc(1            ,-1,-1,NV(K)),1,VECU,LEFT,1, &
                        Vfc(Ncellx+1     ,-1,-1,NV(K)),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
      IF((BCx2.eq.4).and.(IST.eq.NSPLTx-1)) THEN
        DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
          Vfc(IX,JY,KZ,NV(K)) = Vfc(Ncellx,JY,KZ,NV(K))
        END DO;END DO;END DO
      END IF
    END DO
  END IF
!***********************************  BC for chemical species  **********
  IF(N_MPI(10).eq.10) THEN
    CALL MPI_SENDRECV(  ndH(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                        ndH(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndp(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                        ndp(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndH2(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                       ndH2(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndHe(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                       ndHe(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(ndHep(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                      ndHep(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndC(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                        ndC(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCO(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                       ndCO(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCp(Ncellx+1-N_ol,-1,klps),1,VECU,RIGT,1, &
                       ndCp(       1-N_ol,-1,klps),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    IF((BCx1.eq.4).and.(IST.eq.0)) THEN
      DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 0
!        ndH(IX,JY,KZ)   = BBRV_cm(1); ndp(IX,JY,KZ) = BBRV_cm(2); ndH2(IX,JY,KZ) = BBRV_cm(3); ndHe(IX,JY,KZ) = BBRV_cm(4)
!        ndHep(IX,JY,KZ) = BBRV_cm(5); ndC(IX,JY,KZ) = BBRV_cm(6); ndCO(IX,JY,KZ) = BBRV_cm(7); ndCp(IX,JY,KZ) = BBRV_cm(8) 
        ndH(IX,JY,KZ)   = BBRV_cm(1)*U(IX,JY,KZ,1); ndp(IX,JY,KZ)  = BBRV_cm(2)*U(IX,JY,KZ,1)
        ndH2(IX,JY,KZ)  = BBRV_cm(3)*U(IX,JY,KZ,1); ndHe(IX,JY,KZ) = BBRV_cm(4)*U(IX,JY,KZ,1)
        ndHep(IX,JY,KZ) = BBRV_cm(5)*U(IX,JY,KZ,1); ndC(IX,JY,KZ)  = BBRV_cm(6)*U(IX,JY,KZ,1)
        ndCO(IX,JY,KZ)  = BBRV_cm(7)*U(IX,JY,KZ,1); ndCp(IX,JY,KZ) = BBRV_cm(8)*U(IX,JY,KZ,1)
      END DO;END DO;END DO
    END IF

    CALL MPI_SENDRECV(  ndH(1            ,-1,klps),1,VECU,LEFT,1, &
                        ndH(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndp(1            ,-1,klps),1,VECU,LEFT,1, &
                        ndp(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndH2(1            ,-1,klps),1,VECU,LEFT,1, &
                       ndH2(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndHe(1            ,-1,klps),1,VECU,LEFT,1, &
                       ndHe(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(ndHep(1            ,-1,klps),1,VECU,LEFT,1, &
                      ndHep(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndC(1            ,-1,klps),1,VECU,LEFT,1, &
                        ndC(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCO(1            ,-1,klps),1,VECU,LEFT,1, &
                       ndCO(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCp(1            ,-1,klps),1,VECU,LEFT,1, &
                       ndCp(Ncellx+1     ,-1,klps),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
    IF((BCx2.eq.4).and.(IST.eq.NSPLTx-1)) THEN
      DO KZ=1,Ncellz; DO JY=1,Ncelly; DO IX=Ncellx+1,Ncellx+N_ol
!        ndH(IX,JY,KZ)   = BBRV_cm(1); ndp(IX,JY,KZ) = BBRV_cm(2); ndH2(IX,JY,KZ) = BBRV_cm(3); ndHe(IX,JY,KZ) = BBRV_cm(4)
!        ndHep(IX,JY,KZ) = BBRV_cm(5); ndC(IX,JY,KZ) = BBRV_cm(6); ndCO(IX,JY,KZ) = BBRV_cm(7); ndCp(IX,JY,KZ) = BBRV_cm(8)
        ndH(IX,JY,KZ)   = BBRV_cm(1)*U(IX,JY,KZ,1); ndp(IX,JY,KZ)  = BBRV_cm(2)*U(IX,JY,KZ,1)
        ndH2(IX,JY,KZ)  = BBRV_cm(3)*U(IX,JY,KZ,1); ndHe(IX,JY,KZ) = BBRV_cm(4)*U(IX,JY,KZ,1)
        ndHep(IX,JY,KZ) = BBRV_cm(5)*U(IX,JY,KZ,1); ndC(IX,JY,KZ)  = BBRV_cm(6)*U(IX,JY,KZ,1)
        ndCO(IX,JY,KZ)  = BBRV_cm(7)*U(IX,JY,KZ,1); ndCp(IX,JY,KZ) = BBRV_cm(8)*U(IX,JY,KZ,1)
      END DO; END DO; END DO
    END IF
  END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF

IF(iwy.eq.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!***********************************  BC for Lagrange coordinate  **********
  IF(N_MPI(10).eq.1) THEN
    CALL MPI_SENDRECV(dxlagM(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
                      dxlagM(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(dxlagM(-1,1            ,-1),1,VECU,BOTM,1, &
                      dxlagM(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END IF
!********************************************************  BC for Bcc ******
  IF(N_MPI(11).eq.1) THEN
    DO K = 1, 3
      CALL MPI_SENDRECV(Bcc(-1,Ncelly+1-N_ol,-1,K),1,VECU,TOP ,1, &
                        Bcc(-1,       1-N_ol,-1,K),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
      CALL MPI_SENDRECV(Bcc(-1,1            ,-1,K),1,VECU,BOTM,1, &
                        Bcc(-1,Ncelly+1     ,-1,K),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    END DO
  END IF
!********************************************************  BC for Blg ******
  IF(N_MPI(12).eq.1) THEN
    DO K = 1, 2
      CALL MPI_SENDRECV(Blg(-1,Ncelly+1-N_ol,-1,K),1,VECU,TOP ,1, &
                        Blg(-1,       1-N_ol,-1,K),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
      CALL MPI_SENDRECV(Blg(-1,1            ,-1,K),1,VECU,BOTM,1, &
                        Blg(-1,Ncelly+1     ,-1,K),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    END DO
  END IF
!********************************************************  BC for Vfc ******
  IF(N_MPI(13).eq.1) THEN
    NV(1) = 1; NV(2) = 3
    DO K = 1, 2
      CALL MPI_SENDRECV(Vfc(-1,Ncelly+1-N_ol,-1,NV(K)),1,VECU,TOP ,1, &
                        Vfc(-1,       1-N_ol,-1,NV(K)),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
      CALL MPI_SENDRECV(Vfc(-1,1            ,-1,NV(K)),1,VECU,BOTM,1, &
                        Vfc(-1,Ncelly+1     ,-1,NV(K)),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    END DO
  END IF
!***********************************  BC for chemical species  **********
  IF(N_MPI(10).eq.10) THEN
    CALL MPI_SENDRECV(  ndH(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                        ndH(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndp(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                        ndp(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndH2(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                       ndH2(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndHe(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                       ndHe(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(ndHep(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                      ndHep(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndC(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                        ndC(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCO(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                       ndCO(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCp(-1,Ncelly+1-N_ol,klps),1,VECU,TOP ,1, &
                       ndCp(-1,       1-N_ol,klps),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
    
    CALL MPI_SENDRECV(  ndH(-1,1            ,klps),1,VECU,BOTM,1, &
                        ndH(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndp(-1,1            ,klps),1,VECU,BOTM,1, &
                        ndp(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndH2(-1,1            ,klps),1,VECU,BOTM,1, &
                       ndH2(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndHe(-1,1            ,klps),1,VECU,BOTM,1, &
                       ndHe(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(ndHep(-1,1            ,klps),1,VECU,BOTM,1, &
                      ndHep(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndC(-1,1            ,klps),1,VECU,BOTM,1, &
                        ndC(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCO(-1,1            ,klps),1,VECU,BOTM,1, &
                       ndCO(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCp(-1,1            ,klps),1,VECU,BOTM,1, &
                       ndCp(-1,Ncelly+1     ,klps),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF

IF(iwz.eq.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!***********************************  BC for Lagrange coordinate  **********
  IF(N_MPI(10).eq.1) THEN
    CALL MPI_SENDRECV(dxlagM(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                      dxlagM(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(dxlagM(-1,-1,1            ),1,VECU,DOWN,1, &
                      dxlagM(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END IF
!********************************************************  BC for Bcc ******
  IF(N_MPI(11).eq.1) THEN
    DO K = 1, 3
      CALL MPI_SENDRECV(Bcc(-1,-1,Ncellz+1-N_ol,K),1,VECU,UP  ,1, &
                        Bcc(-1,-1,       1-N_ol,K),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
      CALL MPI_SENDRECV(Bcc(-1,-1,1            ,K),1,VECU,DOWN,1, &
                        Bcc(-1,-1,Ncellz+1     ,K),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    END DO
  END IF
!********************************************************  BC for Blg ******
  IF(N_MPI(12).eq.1) THEN
    DO K = 1, 2
      CALL MPI_SENDRECV(Blg(-1,-1,Ncellz+1-N_ol,K),1,VECU,UP  ,1, &
                        Blg(-1,-1,       1-N_ol,K),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
      CALL MPI_SENDRECV(Blg(-1,-1,1            ,K),1,VECU,DOWN,1, &
                        Blg(-1,-1,Ncellz+1     ,K),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    END DO
  END IF
!********************************************************  BC for Vfc ******
  IF(N_MPI(13).eq.1) THEN
    NV(1) = 1; NV(2) = 2
    DO K = 1, 2
      CALL MPI_SENDRECV(Vfc(-1,-1,Ncellz+1-N_ol,NV(K)),1,VECU,UP  ,1, &
                        Vfc(-1,-1,       1-N_ol,NV(K)),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
      CALL MPI_SENDRECV(Vfc(-1,-1,1            ,NV(K)),1,VECU,DOWN,1, &
                        Vfc(-1,-1,Ncellz+1     ,NV(K)),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    END DO
  END IF
!***********************************  BC for chemical species  **********
  IF(N_MPI(10).eq.10) THEN
    CALL MPI_SENDRECV(  ndH(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                        ndH(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndp(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                        ndp(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndH2(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                       ndH2(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndHe(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                       ndHe(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(ndHep(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                      ndHep(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndC(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                        ndC(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCO(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                       ndCO(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCp(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
                       ndCp(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)

    CALL MPI_SENDRECV(  ndH(-1,-1,1            ),1,VECU,DOWN,1, &
                        ndH(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndp(-1,-1,1            ),1,VECU,DOWN,1, &
                        ndp(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndH2(-1,-1,1            ),1,VECU,DOWN,1, &
                       ndH2(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndHe(-1,-1,1            ),1,VECU,DOWN,1, &
                       ndHe(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(ndHep(-1,-1,1            ),1,VECU,DOWN,1, &
                      ndHep(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV(  ndC(-1,-1,1            ),1,VECU,DOWN,1, &
                        ndC(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCO(-1,-1,1            ),1,VECU,DOWN,1, &
                       ndCO(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
    CALL MPI_SENDRECV( ndCp(-1,-1,1            ),1,VECU,DOWN,1, &
                       ndCp(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF

END SUBROUTINE


!*--------------------------------------------------------------*- f77 -*!
!*                 3D ideal MHD                                         *!
!*                 Lagrange and Eulerian remapping scheme               *!
!*                                                                      *!
!*                 E = P/(gamma-1) + rho*v^2/2 + B^2/2                  *!
!*                                                                      *!
!*======================================================================*!
SUBROUTINE MHD(xelr,dxelr,dt)
USE comvar
USE mpivar
USE chmvar

integer :: ix
double precision  :: dt,xelr(-1:ndmax),dxelr(-1:ndmax), ekin,emag, invd, vel
double precision  :: mmean

ALLOCATE(Vfc(-1:ndx,-1:ndy,-1:ndz,1)); ALLOCATE(Blg(-1:ndx,-1:ndy,-1:ndz,2))
ALLOCATE(xlag(-1:ndx,-1:ndy,-1:ndz)); ALLOCATE(dxlagM(-1:ndx,-1:ndy,-1:ndz))
do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
  ix = iwx*i    + iwy*j + iwz*k
  xlag(i,j,k) = xelr(ix)
end do; end do; end do

call RHS(dt,dxelr)
call MOC(dt,xelr)

!do k = 1, Ncellz; do j = 1, Ncelly;do i = 1, Ncellx
!  U(i,j,k,2) = U(i,j,k,2) * U(i,j,k,1); U(i,j,k,3) = U(i,j,k,3) * U(i,j,k,1)
!  U(i,j,k,4) = U(i,j,k,4) * U(i,j,k,1); U(i,j,k,5) = U(i,j,k,5) * U(i,j,k,1)
!end do; end do; end do

call REMAP(dt,xelr)
DEALLOCATE(xlag); DEALLOCATE(dxlagM); DEALLOCATE(Vfc)

if(iwx.eq.1) iBcc=1;if(iwy.eq.1) iBcc=2;if(iwz.eq.1) iBcc=3

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
  U(i,j,k,1) = dmax1( 1.27d0*rmin, U(i,j,k,1) )
  invd = 1.d0/U(i,j,k,1)
  U(i,j,k,2) = U(i,j,k,2)*invd
  U(i,j,k,3) = U(i,j,k,3)*invd
  U(i,j,k,4) = U(i,j,k,4)*invd
  gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1
  ekin = 0.5d0 * (   U(i,j,k,2)**2 + U(i,j,k,3)**2 + U(i,j,k,4)**2 ) * U(i,j,k,1)
  emag = 0.5d0 * ( Bcc(i,j,k,iBcc)**2+Blg(i,j,k,1)**2+Blg(i,j,k,2)**2 )
  U(i,j,k,5) = (  U(i,j,k,5)-ekin-emag  )*gammi1
  U(i,j,k,5) = dmax1( pmin, U(i,j,k,5) )
  U(i,j,k,5) = dmin1( pmax, U(i,j,k,5) )

    ndH(i,j,k) = dmax1( ndHmin  ,ndH(i,j,k)   )
    ndp(i,j,k) = dmax1( ndpmin  ,ndp(i,j,k)   )
   ndH2(i,j,k) = dmax1( ndH2min ,ndH2(i,j,k)  )
   ndHe(i,j,k) = dmax1( ndHemin ,ndHe(i,j,k)  )
  ndHep(i,j,k) = dmax1( ndHepmin,ndHep(i,j,k) )
    ndC(i,j,k) = dmax1( ndCmin  ,ndC(i,j,k)   )
   ndCp(i,j,k) = dmax1( ndCpmin ,ndCp(i,j,k)  )
   ndCO(i,j,k) = dmax1( ndCOmin ,ndCO(i,j,k)  )

  vel = dsqrt(2.d0*ekin*invd)
  U(i,j,k,2) = dsign(1.d0,U(i,j,k,2))*dmin1(dabs(U(i,j,k,2)),dabs(U(i,j,k,2))*30.d0/vel)
  U(i,j,k,3) = dsign(1.d0,U(i,j,k,3))*dmin1(dabs(U(i,j,k,3)),dabs(U(i,j,k,3))*30.d0/vel)
  U(i,j,k,4) = dsign(1.d0,U(i,j,k,4))*dmin1(dabs(U(i,j,k,4)),dabs(U(i,j,k,4))*30.d0/vel)

end do; end do; end do

DEALLOCATE(Blg)
END SUBROUTINE MHD


!======================================================================*
!           Right Hand Sides of Euler's Equations
!======================================================================*
SUBROUTINE RHS(dt,dxelr)
USE comvar
USE mpivar
USE chmvar

integer :: ix,jy,kz,Mnum,Lnum,Ncell,Ncm,Ncl
integer :: BT1,BT2,VN, itrn
double precision  :: dt,grdQ(-1:ndmax,16)
double precision  :: QL1,QL2,QL3,QL4,QR1,QR2,QR3,QR4
double precision  :: Va(0:ndmax), Pa(0:ndmax), dxelr(-1:ndmax)
double precision  :: depend1,depend2,cm
double precision  :: dm(-1:ndmax)
double precision  :: ndHm,ndpm,ndHem,ndHepm,ndH2m

itrn = 5

if(iwx.eq.1) then; Ncell = Ncellx; Ncm = Ncelly; Ncl = Ncellz; BT1 = 2; BT2 = 3; VN = 2; end if
if(iwy.eq.1) then; Ncell = Ncelly; Ncm = Ncellz; Ncl = Ncellx; BT1 = 3; BT2 = 1; VN = 3; end if
if(iwz.eq.1) then; Ncell = Ncellz; Ncm = Ncellx; Ncl = Ncelly; BT1 = 1; BT2 = 2; VN = 4; end if


N_MPI(20) = 3
N_MPI(1)  = 1
N_MPI(2)  = VN
N_MPI(3)  = 5
CALL BC_MPI(2,1)

N_MPI(10) = 10; CALL BC_MPI_OT(2,1); N_MPI(10) = 0 ! for chemical boundary

DO Lnum = 1, Ncl
DO Mnum = 1, Ncm

call VLIMIT(1    ,Mnum,Lnum,grdQ,dxelr,0,1)
call VLIMIT(VN   ,Mnum,Lnum,grdQ,dxelr,0,1)
call VLIMIT(5    ,Mnum,Lnum,grdQ,dxelr,0,1)
call VLIMIT(BT1+8,Mnum,Lnum,grdQ,dxelr,0,1)
call VLIMIT(BT2+8,Mnum,Lnum,grdQ,dxelr,0,1)

do i = 0, Ncell

  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)

  gamma =   3.d0*(ndH(ix ,jy ,kz )+ndp(ix ,jy ,kz )+ndHe(ix ,jy ,kz )+ndHep(ix ,jy ,kz ))+5.d0*ndH2(ix ,jy ,kz )
  gamma = ( 5.d0*(ndH(ix ,jy ,kz )+ndp(ix ,jy ,kz )+ndHe(ix ,jy ,kz )+ndHep(ix ,jy ,kz ))+7.d0*ndH2(ix ,jy ,kz ) )/gamma
  cm = dsqrt(  (gamma * U(ix ,jy ,kz ,5) + Bcc(ix ,jy ,kz ,BT1)**2 + Bcc(ix ,jy ,kz ,BT2)**2) / U(ix ,jy ,kz,1 )  )
  depend1 = 0.5d0*(dxelr(i  ) - cm * dt)*facdep
  depend1 = dmax1(0.d0, depend1)
  gamma =   3.d0*(ndH(ixp,jyp,kzp)+ndp(ixp,jyp,kzp)+ndHe(ixp,jyp,kzp)+ndHep(ixp,jyp,kzp))+5.d0*ndH2(ixp,jyp,kzp)
  gamma = ( 5.d0*(ndH(ixp,jyp,kzp)+ndp(ixp,jyp,kzp)+ndHe(ixp,jyp,kzp)+ndHep(ixp,jyp,kzp))+7.d0*ndH2(ixp,jyp,kzp) )/gamma
  cm = dsqrt(  (gamma * U(ixp,jyp,kzp,5) + Bcc(ixp,jyp,kzp,BT1)**2 + Bcc(ixp,jyp,kzp,BT2)**2) / U(ixp,jyp,kzp,1)  )
  depend2 = 0.5d0*(dxelr(i+1) - cm * dt)*facdep
  depend2 = dmax1(0.d0, depend2)

  QL1 = U(ix ,jy ,kz ,1 ) + depend1 * grdQ(i  ,1 )
  QL1 = (0.5d0-dsign(0.5d0,-QL1))*QL1 + (0.5d0+dsign(0.5d0,-QL1))*U(ix ,jy ,kz ,1)
  QR1 = U(ixp,jyp,kzp,1 ) - depend2 * grdQ(i+1,1 )
  QR1 = (0.5d0-dsign(0.5d0,-QR1))*QR1 + (0.5d0+dsign(0.5d0,-QR1))*U(ixp,jyp,kzp,1)
  QL2 = U(ix ,jy ,kz ,VN) + depend1 * grdQ(i  ,VN)
  QR2 = U(ixp,jyp,kzp,VN) - depend2 * grdQ(i+1,VN)
  QL3 = U(ix ,jy ,kz ,5 ) + depend1 * grdQ(i  ,5 )
  QL3 = (0.5d0-dsign(0.5d0,-QL3))*QL3 + (0.5d0+dsign(0.5d0,-QL3))*U(ix ,jy ,kz ,5)
  QR3 = U(ixp,jyp,kzp,5 ) - depend2 * grdQ(i+1,5 )
  QR3 = (0.5d0-dsign(0.5d0,-QR3))*QR3 + (0.5d0+dsign(0.5d0,-QR3))*U(ixp,jyp,kzp,5)
  QL4 = (Bcc(ix ,jy ,kz ,BT1)+depend1*grdQ(i  ,BT1+8))**2 + (Bcc(ix ,jy ,kz ,BT2)+depend1*grdQ(i  ,BT2+8))**2
  QR4 = (Bcc(ixp,jyp,kzp,BT1)-depend2*grdQ(i+1,BT1+8))**2 + (Bcc(ixp,jyp,kzp,BT2)-depend2*grdQ(i+1,BT2+8))**2

  ndHm =0.5d0*( ndH(ix,jy,kz)+ ndH(ixp,jyp,kzp)); ndpm =0.5d0*(  ndp(ix,jy,kz)+  ndp(ixp,jyp,kzp))
  ndHem=0.5d0*(ndHe(ix,jy,kz)+ndHe(ixp,jyp,kzp));ndHepm=0.5d0*(ndHep(ix,jy,kz)+ndHep(ixp,jyp,kzp))
  ndH2m=0.5d0*(ndH2(ix,jy,kz)+ndH2(ixp,jyp,kzp))
  gamma = ( 5.d0*(ndHm+ndpm+ndHem+ndHepm)+7.d0*ndH2m )/( 3.d0*(ndHm+ndpm+ndHem+ndHepm)+5.d0*ndH2m )
  gammi1 = gamma - 1.0d0
  gammi2 = gamma - 2.0d0
  gammi3 = gamma - 3.0d0
  gampl1 = gamma + 1.0d0
  gampl2 = gamma + 2.0d0
  gampl3 = gamma + 3.0d0
!---------------------------------------------------------*
  call RIEMAN( QL1,QL2,QL3,QL4,QR1,QR2,QR3,QR4,Pa(i),Va(i), itrn )
!---------------------------------------------------------*

end do

!***** store lagrangian mass *****
do i = 0, Ncell
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  dm(i)   = dxelr(i) * U(ix,jy,kz,1)
  xlag(ix,jy,kz) = xlag(ix,jy,kz) + Va(i) * dt
end do
!***** Conservation laws *****

do i = 1, Ncell
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

  gammi1 =   3.d0*(ndH(ix,jy,kz)+ndp(ix,jy,kz)+ndHe(ix,jy,kz)+ndHep(ix,jy,kz))+5.d0*ndH2(ix,jy,kz)
  gammi1 = ( 2.d0*(ndH(ix,jy,kz)+ndp(ix,jy,kz)+ndHe(ix,jy,kz)+ndHep(ix,jy,kz))+2.d0*ndH2(ix,jy,kz) )/gammi1
  U(ix,jy,kz,5 )  = U(ix,jy,kz,5)/( U(ix,jy,kz,1)*gammi1 ) &
                  + 0.5d0*(   U(ix,jy,kz,2)**2.d0+  U(ix,jy,kz,3)**2.d0+  U(ix,jy,kz,4)**2.d0 ) &
                  + 0.5d0*( Bcc(ix,jy,kz,1)**2.d0+Bcc(ix,jy,kz,2)**2.d0+Bcc(ix,jy,kz,3)**2.d0 )/U(ix,jy,kz,1)
  U(ix,jy,kz,1 )  = dm(i) / ( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
  U(ix,jy,kz,VN)  = U(ix,jy,kz,VN)- dt / dm(i) * ( Pa(i) - Pa(i-1) )
  U(ix,jy,kz,5 )  = U(ix,jy,kz,5 )- dt / dm(i) * ( Pa(i)*Va(i)-Pa(i-1)*Va(i-1) )
  Blg(ix,jy,kz,1) = Bcc(ix,jy,kz,BT1)*dxelr(i) / ( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
  Blg(ix,jy,kz,2) = Bcc(ix,jy,kz,BT2)*dxelr(i) / ( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )

    ndH(ix,jy,kz) = dxelr(i)*  ndH(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
    ndp(ix,jy,kz) = dxelr(i)*  ndp(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
   ndH2(ix,jy,kz) = dxelr(i)* ndH2(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
   ndHe(ix,jy,kz) = dxelr(i)* ndHe(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
  ndHep(ix,jy,kz) = dxelr(i)*ndHep(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
    ndC(ix,jy,kz) = dxelr(i)*  ndC(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
   ndCp(ix,jy,kz) = dxelr(i)* ndCp(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
   ndCO(ix,jy,kz) = dxelr(i)* ndCO(ix,jy,kz)/( xlag(ix,jy,kz) - xlag(ixm,jym,kzm) )
end do

do i = 1, Ncell+1
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i

  Vfc(ix,jy,kz,1) = Va(i-1)
end do

END DO
END DO

END SUBROUTINE RHS


!*======================================================================*
!*                 The Riemann Solver for Adiabatic Gases               *
!*                    with Tangential Magnetic Fields                   *
!*======================================================================*
SUBROUTINE RIEMAN( uLi1,uLi2,uLi3,uLi4,uRi1,uRi2,uRi3,uRi4,Pai,Vai, itrn )
USE comvar
integer :: itrn
double precision  :: Pai,Vai,uLi1,uLi2,uLi3,uLi4,uRi1,uRi2,uRi3,uRi4
double precision  :: D1i,D2i,V1i,V2i,P1i,P2i,B1i,B2i,GP1i,GP2i,GP1B1i,GP2B2i,pp1i,pp2i,qM1i,qM2i,qM1sqi,qM2sqi,PaOLDi

D1i = uLi1
V1i = uLi2
P1i = uLi3
B1i = 0.5d0 * uLi4
D2i = uRi1
V2i = uRi2
P2i = uRi3
B2i = 0.5d0 * uRi4

!*****    B1 = 0.5 * B1^2     B2 = 0.5 * B2^2

P1i    = P1i + B1i
P2i    = P2i + B2i
GP1i   = GAMmi1*P1i
GP2i   = GAMmi1*P2i
GP1B1i = GAMmi3*P1i - GAMmi2*B1i*2.d0
GP2B2i = GAMmi3*P2i - GAMmi2*B2i*2.d0

Pai  = ( P1i + P2i )*0.5d0

!do loop = 1, itrn

  pp1i = GAMpl3*Pai + GP1B1i 
  qM1sqi = pp1i+dsqrt( pp1i*pp1i-8.d0*(Pai-P1i)*(GAMpl1*Pai+GP1i) )
  qM1i   = dsqrt( qM1sqi*D1i )*0.5d0
  pp2i = GAMpl3*Pai + GP2B2i 
  qM2sqi = pp2i+dsqrt( pp2i*pp2i-8.d0*(Pai-P2i)*(GAMpl1*Pai+GP2i) )
  qM2i   = dsqrt( qM2sqi*D2i )*0.5d0
  PaOLDi = Pai
  Pai    = qM2i*P1i + qM1i*P2i + qM2i*qM1i*(V1i-V2i)
  Pai    = Pai/( qM2i + qM1i )
  if(Pai.le.0.0) Pai = 0.1d0*PaOLDi

  pp1i = GAMpl3*Pai + GP1B1i 
  qM1sqi = pp1i+dsqrt( pp1i*pp1i-8.d0*(Pai-P1i)*(GAMpl1*Pai+GP1i) )
  qM1i   = dsqrt( qM1sqi*D1i )*0.5d0
  pp2i = GAMpl3*Pai + GP2B2i 
  qM2sqi = pp2i+dsqrt( pp2i*pp2i-8.d0*(Pai-P2i)*(GAMpl1*Pai+GP2i) )
  qM2i   = dsqrt( qM2sqi*D2i )*0.5d0
  PaOLDi = Pai
  Pai    = qM2i*P1i + qM1i*P2i + qM2i*qM1i*(V1i-V2i)
  Pai    = Pai/( qM2i + qM1i )
  if(Pai.le.0.0) Pai = 0.1d0*PaOLDi

  pp1i = GAMpl3*Pai + GP1B1i 
  qM1sqi = pp1i+dsqrt( pp1i*pp1i-8.d0*(Pai-P1i)*(GAMpl1*Pai+GP1i) )
  qM1i   = dsqrt( qM1sqi*D1i )*0.5d0
  pp2i = GAMpl3*Pai + GP2B2i 
  qM2sqi = pp2i+dsqrt( pp2i*pp2i-8.d0*(Pai-P2i)*(GAMpl1*Pai+GP2i) )
  qM2i   = dsqrt( qM2sqi*D2i )*0.5d0
  PaOLDi = Pai
  Pai    = qM2i*P1i + qM1i*P2i + qM2i*qM1i*(V1i-V2i)
  Pai    = Pai/( qM2i + qM1i )
  if(Pai.le.0.0) Pai = 0.1d0*PaOLDi

  pp1i = GAMpl3*Pai + GP1B1i 
  qM1sqi = pp1i+dsqrt( pp1i*pp1i-8.d0*(Pai-P1i)*(GAMpl1*Pai+GP1i) )
  qM1i   = dsqrt( qM1sqi*D1i )*0.5d0
  pp2i = GAMpl3*Pai + GP2B2i 
  qM2sqi = pp2i+dsqrt( pp2i*pp2i-8.d0*(Pai-P2i)*(GAMpl1*Pai+GP2i) )
  qM2i   = dsqrt( qM2sqi*D2i )*0.5d0
  PaOLDi = Pai
  Pai    = qM2i*P1i + qM1i*P2i + qM2i*qM1i*(V1i-V2i)
  Pai    = Pai/( qM2i + qM1i )
  if(Pai.le.0.0) Pai = 0.1d0*PaOLDi

!end do

Vai = qM1i*V1i + qM2i*V2i + P1i - P2i
Vai = Vai / ( qM2i + qM1i )

END SUBROUTINE RIEMAN

!=====================================================================*
!                          Eulerian Remaping
!=====================================================================*
SUBROUTINE REMAP(dt,xelr)
USE comvar
USE mpivar
USE chmvar

integer :: Mnum,Lnum,Ncell,Ncm,Ncl,VN
double precision  :: xelr(-1:ndmax),dt
double precision  :: dxlag(-1:ndmax), F(0:ndmax,7), grdU(-1:ndmax,16)
double precision  :: depend, vadt
double precision  :: grdC(-1:ndmax,8),G(0:ndmax,8)

if(iwx.eq.1) then; Ncell = Ncellx; Ncm = Ncelly; Ncl = Ncellz; VN = 2; end if
if(iwy.eq.1) then; Ncell = Ncelly; Ncm = Ncellz; Ncl = Ncellx; VN = 3; end if
if(iwz.eq.1) then; Ncell = Ncellz; Ncm = Ncellx; Ncl = Ncelly; VN = 4; end if

N_MPI(12) = 1 !EB
CALL BC_MPI_OT(2,2)
N_MPI(12) = 0

N_MPI(20) = 5
N_MPI(1)  = 1
N_MPI(2)  = 2
N_MPI(3)  = 3
N_MPI(4)  = 4
N_MPI(5)  = 5
CALL BC_MPI(2,2)

!DO Lnum = 1, Ncl; DO Mnum = 1, Ncm; do i = 1, Ncell !Calculated in MOC
!  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
!  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
!  kz  = iwx*Lnum + iwy*Mnum + iwz*i
!  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
!  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
!  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
!  dxlagM(ix,jy,kz) = xlag(ix,jy,kz)-xlag(ixm,jym,kzm)
!end do; END DO; END DO
!N_MPI(10) = 1; CALL BC_MPI_OT(2,1); N_MPI(10) = 0

N_MPI(10) = 10; CALL BC_MPI_OT(2,1); N_MPI(10) = 0 ! for chemical boundary

DO Lnum = 1, Ncl
DO Mnum = 1, Ncm

do i = -1, Ncell+2
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  dxlag(i) = dxlagM(ix,jy,kz)
end do

call VLIMIT(1,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(2,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(3,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(4,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(5,Mnum,Lnum,grdU,dxlag,0,1)

call VLIMIT(12,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(13,Mnum,Lnum,grdU,dxlag,0,1)

call VLIMIT_OT(Mnum,Lnum,grdC,dxlag,0,1) !for chemical

do i = 0, Ncell
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  vadt= Vfc(ixp,jyp,kzp,1)*dt
  if(vadt.le.0.d0) then
    depend = 0.5d0 * (dxlag(i+1) + vadt) * facdep
    F(i,1) = -vadt * ( U(ixp,jyp,kzp,1) - depend*grdU(i+1,1) )
    F(i,2) = -vadt * ( U(ixp,jyp,kzp,2) - depend*grdU(i+1,2) )
    F(i,3) = -vadt * ( U(ixp,jyp,kzp,3) - depend*grdU(i+1,3) )
    F(i,4) = -vadt * ( U(ixp,jyp,kzp,4) - depend*grdU(i+1,4) )
    F(i,5) = -vadt * ( U(ixp,jyp,kzp,5) - depend*grdU(i+1,5) )
    F(i,6) = -vadt * ( Blg(ixp,jyp,kzp,1) - depend*grdU(i+1,12) )
    F(i,7) = -vadt * ( Blg(ixp,jyp,kzp,2) - depend*grdU(i+1,13) )
    G(i,1) = -vadt * (   ndH(ixp,jyp,kzp) - depend*grdC(i+1,1) )
    G(i,2) = -vadt * (   ndp(ixp,jyp,kzp) - depend*grdC(i+1,2) )
    G(i,3) = -vadt * (  ndH2(ixp,jyp,kzp) - depend*grdC(i+1,3) )
    G(i,4) = -vadt * (  ndHe(ixp,jyp,kzp) - depend*grdC(i+1,4) )
    G(i,5) = -vadt * ( ndHep(ixp,jyp,kzp) - depend*grdC(i+1,5) )
    G(i,6) = -vadt * (   ndC(ixp,jyp,kzp) - depend*grdC(i+1,6) )
    G(i,7) = -vadt * (  ndCp(ixp,jyp,kzp) - depend*grdC(i+1,7) )
    G(i,8) = -vadt * (  ndCO(ixp,jyp,kzp) - depend*grdC(i+1,8) )
  else  !*** va > 0
    depend = 0.5d0 * (dxlag(i  ) - vadt) * facdep
    F(i,1) = -vadt * ( U(ix ,jy ,kz ,1) + depend*grdU(i  ,1) )
    F(i,2) = -vadt * ( U(ix ,jy ,kz ,2) + depend*grdU(i  ,2) )
    F(i,3) = -vadt * ( U(ix ,jy ,kz ,3) + depend*grdU(i  ,3) )
    F(i,4) = -vadt * ( U(ix ,jy ,kz ,4) + depend*grdU(i  ,4) )
    F(i,5) = -vadt * ( U(ix ,jy ,kz ,5) + depend*grdU(i  ,5) )
    F(i,6) = -vadt * ( Blg(ix ,jy ,kz ,1) + depend*grdU(i  ,12) )
    F(i,7) = -vadt * ( Blg(ix ,jy ,kz ,2) + depend*grdU(i  ,13) )
    G(i,1) = -vadt * (   ndH(ix ,jy ,kz ) + depend*grdC(i  ,1) )
    G(i,2) = -vadt * (   ndp(ix ,jy ,kz ) + depend*grdC(i  ,2) )
    G(i,3) = -vadt * (  ndH2(ix ,jy ,kz ) + depend*grdC(i  ,3) )
    G(i,4) = -vadt * (  ndHe(ix ,jy ,kz ) + depend*grdC(i  ,4) )
    G(i,5) = -vadt * ( ndHep(ix ,jy ,kz ) + depend*grdC(i  ,5) )
    G(i,6) = -vadt * (   ndC(ix ,jy ,kz ) + depend*grdC(i  ,6) )
    G(i,7) = -vadt * (  ndCp(ix ,jy ,kz ) + depend*grdC(i  ,7) )
    G(i,8) = -vadt * (  ndCO(ix ,jy ,kz ) + depend*grdC(i  ,8) )
  end if
end do

do i = 1, Ncell
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  U(ix,jy,kz,1) = ( U(ix,jy,kz,1)*dxlag(i) + F(i,1) - F(i-1,1) )/( xelr(i) - xelr(i-1) )
  U(ix,jy,kz,2) = ( U(ix,jy,kz,2)*dxlag(i) + F(i,2) - F(i-1,2) )/( xelr(i) - xelr(i-1) )
  U(ix,jy,kz,3) = ( U(ix,jy,kz,3)*dxlag(i) + F(i,3) - F(i-1,3) )/( xelr(i) - xelr(i-1) )
  U(ix,jy,kz,4) = ( U(ix,jy,kz,4)*dxlag(i) + F(i,4) - F(i-1,4) )/( xelr(i) - xelr(i-1) )
  U(ix,jy,kz,5) = ( U(ix,jy,kz,5)*dxlag(i) + F(i,5) - F(i-1,5) )/( xelr(i) - xelr(i-1) )
  Blg(ix,jy,kz,1) = ( Blg(ix,jy,kz,1)*dxlag(i) + F(i,6) - F(i-1,6) )/( xelr(i) - xelr(i-1) )
  Blg(ix,jy,kz,2) = ( Blg(ix,jy,kz,2)*dxlag(i) + F(i,7) - F(i-1,7) )/( xelr(i) - xelr(i-1) )
    ndH(ix,jy,kz) = (   ndH(ix,jy,kz)*dxlag(i) + G(i,1) - G(i-1,1) )/( xelr(i) - xelr(i-1) )
    ndp(ix,jy,kz) = (   ndp(ix,jy,kz)*dxlag(i) + G(i,2) - G(i-1,2) )/( xelr(i) - xelr(i-1) )
   ndH2(ix,jy,kz) = (  ndH2(ix,jy,kz)*dxlag(i) + G(i,3) - G(i-1,3) )/( xelr(i) - xelr(i-1) )
   ndHe(ix,jy,kz) = (  ndHe(ix,jy,kz)*dxlag(i) + G(i,4) - G(i-1,4) )/( xelr(i) - xelr(i-1) )
  ndHep(ix,jy,kz) = ( ndHep(ix,jy,kz)*dxlag(i) + G(i,5) - G(i-1,5) )/( xelr(i) - xelr(i-1) )
    ndC(ix,jy,kz) = (   ndC(ix,jy,kz)*dxlag(i) + G(i,6) - G(i-1,6) )/( xelr(i) - xelr(i-1) )
   ndCp(ix,jy,kz) = (  ndCp(ix,jy,kz)*dxlag(i) + G(i,7) - G(i-1,7) )/( xelr(i) - xelr(i-1) )
   ndCO(ix,jy,kz) = (  ndCO(ix,jy,kz)*dxlag(i) + G(i,8) - G(i-1,8) )/( xelr(i) - xelr(i-1) )
end do

END DO
END DO

END SUBROUTINE REMAP


!======================================================================*
!                      Limiting (Van Albada)
!======================================================================*
SUBROUTINE VLIMIT(k,Mnum,Lnum,grdU,dxx,i_sta,i_end)
USE comvar

integer :: Mnum,Lnum,Ncell
double precision, parameter :: eps = 1.d-10
double precision  :: grdU(-1:ndmax,16), dxx(-1:ndmax)
double precision  :: delp,delm,flmt,T,gmm

if(iwx.eq.1) Ncell = Ncellx
if(iwy.eq.1) Ncell = Ncelly
if(iwz.eq.1) Ncell = Ncellz

if(k.le.8) then
do i = i_sta, Ncell+i_end
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

  delp = U(ixp,jyp,kzp,k)-U(ix ,jy ,kz ,k)
  delm = U(ix ,jy ,kz ,k)-U(ixm,jym,kzm,k)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
  grdU(i,k) = flmt*( U(ixp,jyp,kzp,k)-U(ixm,jym,kzm,k) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  T = 1.27d0*U(ix,jy,kz,5)/(kb*U(ix,jy,kz,1))
  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,k) = grdU(i,k)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))
 
end do
end if
if((k.ge.9).and.(k.le.11)) then
do i = i_sta, Ncell+i_end
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
  kbc = k-8

  delp = Bcc(ixp,jyp,kzp,kbc)-Bcc(ix ,jy ,kz ,kbc)
  delm = Bcc(ix ,jy ,kz ,kbc)-Bcc(ixm,jym,kzm,kbc)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
  grdU(i,k) = flmt*( Bcc(ixp,jyp,kzp,kbc)-Bcc(ixm,jym,kzm,kbc) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  T = 1.27d0*U(ix,jy,kz,5)/(kb*U(ix,jy,kz,1))
  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,k) = grdU(i,k)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))

end do
end if
if((k.ge.12).and.(k.le.13)) then
do i = i_sta, Ncell+i_end
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
  kbc = k-11

  delp = Blg(ixp,jyp,kzp,kbc)-Blg(ix ,jy ,kz ,kbc)
  delm = Blg(ix ,jy ,kz ,kbc)-Blg(ixm,jym,kzm,kbc)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
  grdU(i,k) = flmt*( Blg(ixp,jyp,kzp,kbc)-Blg(ixm,jym,kzm,kbc) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  T = 1.27d0*U(ix,jy,kz,5)/(kb*U(ix,jy,kz,1))
  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,k) = grdU(i,k)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))

end do
end if
if((k.ge.14).and.(k.le.16)) then
do i = i_sta, Ncell+i_end
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
  kbc = k-13

  delp = Vfc(ixp,jyp,kzp,kbc)-Vfc(ix ,jy ,kz ,kbc)
  delm = Vfc(ix ,jy ,kz ,kbc)-Vfc(ixm,jym,kzm,kbc)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
  grdU(i,k) = flmt*( Vfc(ixp,jyp,kzp,kbc)-Vfc(ixm,jym,kzm,kbc) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  T = 1.27d0*U(ix,jy,kz,5)/(kb*U(ix,jy,kz,1))
  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,k) = grdU(i,k)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))

end do
end if

END SUBROUTINE


SUBROUTINE VLIMIT_OT(Mnum,Lnum,grdU,dxx,i_sta,i_end)
USE comvar
USE chmvar
integer :: Mnum,Lnum,Ncell
double precision, parameter :: eps=1.d-10
double precision  :: grdU(-1:ndmax,8), dxx(-1:ndmax)
double precision  :: delp,delm,flmt,T,gmm

if(iwx.eq.1) Ncell = Ncellx
if(iwy.eq.1) Ncell = Ncelly
if(iwz.eq.1) Ncell = Ncellz

do i = i_sta, Ncell+i_end
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

  T = 1.27d0*U(ix,jy,kz,5)/(kb*U(ix,jy,kz,1))

  delp = ndH(ixp,jyp,kzp)-ndH(ix ,jy ,kz ); delm = ndH(ix ,jy ,kz )-ndH(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,1) = flmt*( ndH(ixp,jyp,kzp)-ndH(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,1) = grdU(i,1)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndp(ixp,jyp,kzp)-ndp(ix ,jy ,kz ); delm = ndp(ix ,jy ,kz )-ndp(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,2) = flmt*( ndp(ixp,jyp,kzp)-ndp(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,2) = grdU(i,2)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndH2(ixp,jyp,kzp)-ndH2(ix ,jy ,kz ); delm = ndH2(ix ,jy ,kz )-ndH2(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,3) = flmt*( ndH2(ixp,jyp,kzp)-ndH2(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,3) = grdU(i,3)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndHe(ixp,jyp,kzp)-ndHe(ix ,jy ,kz ); delm = ndHe(ix ,jy ,kz )-ndHe(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,4) = flmt*( ndHe(ixp,jyp,kzp)-ndHe(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,4) = grdU(i,4)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndHep(ixp,jyp,kzp)-ndHep(ix ,jy ,kz ); delm = ndHep(ix ,jy ,kz )-ndHep(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,5) = flmt*( ndHep(ixp,jyp,kzp)-ndHep(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,5) = grdU(i,5)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndC(ixp,jyp,kzp)-ndC(ix ,jy ,kz ); delm = ndC(ix ,jy ,kz )-ndC(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,6) = flmt*( ndC(ixp,jyp,kzp)-ndC(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,6) = grdU(i,6)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndCp(ixp,jyp,kzp)-ndCp(ix ,jy ,kz ); delm = ndCp(ix ,jy ,kz )-ndCp(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,7) = flmt*( ndCp(ixp,jyp,kzp)-ndCp(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,7) = grdU(i,7)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


  delp = ndCO(ixp,jyp,kzp)-ndCO(ix ,jy ,kz ); delm = ndCO(ix ,jy ,kz )-ndCO(ixm,jym,kzm)
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps))
  grdU(i,8) = flmt*( ndCO(ixp,jyp,kzp)-ndCO(ixm,jym,kzm) )/( dxx(i)+0.5d0*dxx(i-1)+0.5d0*dxx(i+1) )

  delp = 2.d0*delp/(dxx(i)+dxx(i+1)); delm = 2.d0*delm/(dxx(i)+dxx(i-1))
  gmm = (0.5d0+dsign(0.5d0,delp*delm))*dsign(1.d0,delp)*dmin1(dabs(delp),dabs(delm)) !minmod
  grdU(i,8) = grdU(i,8)*(0.5d0-dsign(0.5d0,T-3.d0)) + gmm*(0.5d0+dsign(0.5d0,T-3.d0))


end do

END SUBROUTINE VLIMIT_OT


!=====================================================================*
!     Courant Condition for Euler code
!=====================================================================*
SUBROUTINE Couran(tCFL)
USE comvar
USE mpivar
USE chmvar

double precision :: tCFL,c2

tCFL = tfinal
do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
  gamma =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  gamma = ( 5.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+7.d0*ndH2(i,j,k) )/gamma
  c2 = ( gamma * U(i,j,k,5) + Bcc(i,j,k,1)**2.d0+Bcc(i,j,k,2)**2.d0+Bcc(i,j,k,3)**2.d0 ) / U(i,j,k,1)
  if(U(i,j,k,1) .ne. 0.0d0) then
  tCFL = dmin1(tCFL, dx(i)/(dsqrt(c2) + dabs(U(i,j,k,2))) )
  tCFL = dmin1(tCFL, dy(j)/(dsqrt(c2) + dabs(U(i,j,k,3))) )
  tCFL = dmin1(tCFL, dz(k)/(dsqrt(c2) + dabs(U(i,j,k,4))) )
end if
end do; end do; end do
if(tCFL.lt.0.d0) write(5,*) time,NRANK,'err at Couran'

END SUBROUTINE Couran


!*-------------------------------------------------------------*- f90 -*!
!*          Solve Lorentz Force With Method of Characteristic          *!
!*                                   programed by T.Inoue              *!
!*                                        2006.11.                     *!
!*=====================================================================*!

SUBROUTINE MOC(dt,xelr)
USE comvar
USE mpivar

integer :: Mnum,Lnum,Ncell,Ncm,Ncl, BN,VT1,VT2,VN
double precision  :: dt, xelr(-1:ndmax), grdU(-1:ndmax,16), dxlag(-1:ndmax)
double precision  :: ca, depend1, depend2, deni
double precision  :: Bnap, Bnam, Vnap, Vnam
double precision  :: QL11,QL12,QL13,QR11,QR12,QR13
double precision  :: QL21,QL22,QL23,QR21,QR22,QR23
double precision  :: Vay(-1:ndmax), Bay(-1:ndmax), Vaz(-1:ndmax), Baz(-1:ndmax)

if(iwx.eq.1) then; Ncell=Ncellx; Ncm=Ncelly; Ncl=Ncellz; VN=2; VT1=3; VT2=4; BN=6; end if
if(iwy.eq.1) then; Ncell=Ncelly; Ncm=Ncellz; Ncl=Ncellx; VN=3; VT1=4; VT2=2; BN=7; end if
if(iwz.eq.1) then; Ncell=Ncellz; Ncm=Ncellx; Ncl=Ncelly; VN=4; VT1=2; VT2=3; BN=8; end if

N_MPI(20) = 4
N_MPI(1)  = 1
N_MPI(2)  = VT1
N_MPI(3)  = VT2
N_MPI(4)  = BN
CALL BC_MPI(2,1)

DO Lnum = 1, Ncl; DO Mnum = 1, Ncm; do i = 1, Ncell
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
  dxlagM(ix,jy,kz) = xlag(ix,jy,kz)-xlag(ixm,jym,kzm)
end do; END DO; END DO

N_MPI(10) = 1
N_MPI(12) = 1
CALL BC_MPI_OT(2,1)
N_MPI(10) = 0
N_MPI(12) = 0

DO Lnum = 1, Ncl
DO Mnum = 1, Ncm

do i = -1, Ncell+2
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  dxlag(i) = dxlagM(ix,jy,kz)
end do

!***** method of characteristic *****
call VLIMIT(1  ,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(VT1,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(VT2,Mnum,Lnum,grdU,dxlag,0,1)
call VLIMIT(12 ,Mnum,Lnum,grdU,dxlag,0,1) !BT1
call VLIMIT(13 ,Mnum,Lnum,grdU,dxlag,0,1) !BT2
!call VLIMIT(Bn ,Mnum,Lnum,grdU,dxlag,0,1)

do i = 0, Ncell

  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixpp= iwx*(i+2)+ iwy*Lnum + iwz*Mnum
  jypp= iwx*Mnum + iwy*(i+2)+ iwz*Lnum
  kzpp= iwx*Lnum + iwy*Mnum + iwz*(i+2)

  Bnap = U(ixp,jyp,kzp,Bn)+(xlag(ix,jy,kz)-xelr(i))*(U(ixpp,jypp,kzpp,Bn)-U(ix,jy,kz,Bn))/(xelr(i+1)-xelr(i-1))
  deni = dsqrt( (dxlag(i)+dxlag(i+1))/(dxlag(i)*U(ixp,jyp,kzp,1)+dxlag(i+1)*U(ix,jy,kz,1)) )
  if( Bnap.gt.0.d0) then
    ca = dabs( Bnap * deni )*0.5d0*dt
    ca = dmin1( dxlag(i  ), ca )
    depend1 = ( 0.5d0*dxlag(i  ) - ca )*facdep
    ca = dabs( Bnap * deni )*0.5d0*dt
    ca = dmin1( dxlag(i+1), ca )
    depend2 = ( 0.5d0*dxlag(i+1) - ca )*facdep
    QL11 =   U(ix ,jy ,kz ,VT1) + depend1 * grdU(i  ,VT1) ! Vt cc
    QR11 =   U(ixp,jyp,kzp,VT1) - depend2 * grdU(i+1,VT1)
    QL12 = Blg(ix ,jy ,kz ,1  ) + depend1 * grdU(i  ,12 ) ! Bt cc
    QR12 = Blg(ixp,jyp,kzp,1  ) - depend2 * grdU(i+1,12 )
    QL13 =   U(ix ,jy ,kz ,1  ) + depend1 * grdU(i  ,1  )
    QL13 =   (0.5d0-dsign(0.5d0,-QL13))*QL13 + (0.5d0+dsign(0.5d0,-QL13))*U(ix ,jy ,kz ,1)
    QR13 =   U(ixp,jyp,kzp,1  ) - depend2 * grdU(i+1,1  )
    QR13 =   (0.5d0-dsign(0.5d0,-QR13))*QR13 + (0.5d0+dsign(0.5d0,-QR13))*U(ixp,jyp,kzp,1)
    QL21 =   U(ix ,jy ,kz ,VT2) + depend1 * grdU(i  ,VT2) ! Vz cc
    QR21 =   U(ixp,jyp,kzp,VT2) - depend2 * grdU(i+1,VT2)
    QL22 = Blg(ix ,jy ,kz ,2  ) + depend1 * grdU(i  ,13 ) ! Bz cc
    QR22 = Blg(ixp,jyp,kzp,2  ) - depend2 * grdU(i+1,13 )
    QL23 = QL13
    QR23 = QR13
  else
    ca = dabs( Bnap * deni )*0.5d0*dt
    ca = dmin1( dxlag(i+1), ca )
    depend1 = ( 0.5d0*dxlag(i+1) - ca )*facdep
    ca = dabs( Bnap * deni )*0.5d0*dt
    ca = dmin1( dxlag(i  ), ca )
    depend2 = ( 0.5d0*dxlag(i  ) - ca )*facdep
    QL11 =   U(ixp,jyp,kzp,VT1) - depend1 * grdU(i+1,VT1) ! Vt cc
    QR11 =   U(ix ,jy ,kz ,VT1) + depend2 * grdU(i  ,VT1)
    QL12 = Blg(ixp,jyp,kzp,1  ) - depend1 * grdU(i+1,12 ) ! Bt cc
    QR12 = Blg(ix ,jy ,kz ,1  ) + depend2 * grdU(i  ,12 )
    QL13 =   U(ixp,jyp,kzp,1  ) - depend1 * grdU(i+1,1  )
    QL13 =   (0.5d0-dsign(0.5d0,-QL13))*QL13 + (0.5d0+dsign(0.5d0,-QL13))*U(ixp,jyp,kzp,1)
    QR13 =   U(ix ,jy ,kz ,1  ) + depend2 * grdU(i  ,1  )
    QR13 =   (0.5d0-dsign(0.5d0,-QR13))*QR13 + (0.5d0+dsign(0.5d0,-QR13))*U(ix ,jy ,kz ,1)
    QL21 =   U(ixp,jyp,kzp,VT2) - depend1 * grdU(i+1,VT2) ! Vz cc
    QR21 =   U(ix ,jy ,kz ,VT2) + depend2 * grdU(i  ,VT2)
    QL22 = Blg(ixp,jyp,kzp,2  ) - depend1 * grdU(i+1,13 ) ! Bz cc
    QR22 = Blg(ix ,jy ,kz ,2  ) + depend2 * grdU(i  ,13 )
    QL23 = QL13
    QR23 = QR13
  end if
!--------------------------------------------------------------*
  call CHAREQ( QL11,QL12,QL13,QR11,QR12,QR13,Bay(i),Vay(i) )
  call CHAREQ( QL21,QL22,QL23,QR21,QR22,QR23,Baz(i),Vaz(i) )
!--------------------------------------------------------------*
end do

do i = 1, Ncell
  ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
  jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
  kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
  ix  = iwx*i    + iwy*Lnum + iwz*Mnum
  jy  = iwx*Mnum + iwy*i    + iwz*Lnum
  kz  = iwx*Lnum + iwy*Mnum + iwz*i
  ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
  jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
  kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
  ixpp= iwx*(i+2)+ iwy*Lnum + iwz*Mnum
  jypp= iwx*Mnum + iwy*(i+2)+ iwz*Lnum
  kzpp= iwx*Lnum + iwy*Mnum + iwz*(i+2)

  Bnap = U(ixp,jyp,kzp,Bn)+(xlag(ix ,jy ,kz )-xelr(i  ))*(U(ixpp,jypp,kzpp,Bn)-U(ix ,jy ,kz ,Bn))/(xelr(i+1)-xelr(i-1))
  Bnam = U(ix ,jy ,kz ,Bn)+(xlag(ixm,jym,kzm)-xelr(i-1))*(U(ixp ,jyp ,kzp ,Bn)-U(ixm,jym,kzm,Bn))/(xelr(i  )-xelr(i-2))
  Vnap = Vfc(ixp,jyp,kzp,1 )
  Vnam = Vfc(ix ,jy ,kz ,1 )
!  U(ix,jy,kz,VN ) = U(ix,jy,kz,VN ) + dt / U(ix,jy,kz,1) / dxlag(i) * 0.5d0*( Bnap**2.d0 - Bnam**2.d0 )
!  U(ix,jy,kz,VT1) = U(ix,jy,kz,VT1) + dt / U(ix,jy,kz,1) / dxlag(i) * ( Bnap * Bay(i) - Bnam * Bay(i-1) )
!  U(ix,jy,kz,VT2) = U(ix,jy,kz,VT2) + dt / U(ix,jy,kz,1) / dxlag(i) * ( Bnap * Baz(i) - Bnam * Baz(i-1) )
!  U(ix,jy,kz,5  ) = U(ix,jy,kz,5  ) + dt / U(ix,jy,kz,1) / dxlag(i) * &
!     (  (0.5d0*Bnap*Vnap+Bay(i  )*Vay(i  )+Baz(i  )*Vaz(i  )) * Bnap  &
!       -(0.5d0*Bnam*Vnam+Bay(i-1)*Vay(i-1)+Baz(i-1)*Vaz(i-1)) * Bnam  )

  U(ix,jy,kz,VN ) = U(ix,jy,kz,VN )*U(ix,jy,kz,1) + dt / dxlag(i) * 0.5d0*( Bnap**2.d0 - Bnam**2.d0 )
  U(ix,jy,kz,VT1) = U(ix,jy,kz,VT1)*U(ix,jy,kz,1) + dt / dxlag(i) * ( Bnap * Bay(i) - Bnam * Bay(i-1) )
  U(ix,jy,kz,VT2) = U(ix,jy,kz,VT2)*U(ix,jy,kz,1) + dt / dxlag(i) * ( Bnap * Baz(i) - Bnam * Baz(i-1) )
  U(ix,jy,kz,5  ) = U(ix,jy,kz,5  )*U(ix,jy,kz,1) + dt / dxlag(i) * &
     (  (0.5d0*Bnap*Vnap+Bay(i  )*Vay(i  )+Baz(i  )*Vaz(i  ))*Bnap  &
       -(0.5d0*Bnam*Vnam+Bay(i-1)*Vay(i-1)+Baz(i-1)*Vaz(i-1))*Bnam  )

  Blg(ix,jy,kz,1) = Blg(ix,jy,kz,1) + dt / dxlag(i) * ( Bnap * Vay(i) - Bnam * Vay(i-1) )
  Blg(ix,jy,kz,2) = Blg(ix,jy,kz,2) + dt / dxlag(i) * ( Bnap * Vaz(i) - Bnam * Vaz(i-1) )
end do

END DO
END DO

END SUBROUTINE MOC


SUBROUTINE CHAREQ( QL1,QL2,QL3, QR1,QR2,QR3, Bay, Vay )
DOUBLE PRECISION :: QL1,QL2,QL3, QR1,QR2,QR3, Bay, Vay
double precision :: sqdyp, sqdym, vyp, vym, Byp, Bym

vyp   = QL1; vym   = QR1 ; Byp   = QL2; Bym   = QR2
sqdyp = dsqrt( QL3 ); sqdym = dsqrt( QR3 )

Vay = ( sqdyp*vyp + sqdym*vym + Bym - Byp )/(sqdyp+sqdym)
Bay = ( sqdyp*sqdym*(vym-vyp)+sqdyp*Bym+sqdym*Byp )/(sqdyp+sqdym)

END SUBROUTINE CHAREQ



SUBROUTINE CC(mode,dt)
USE comvar
USE mpivar
integer :: mode
double precision  :: dt,dS

if(mode.eq.1) then !FCVF & FCMFBC
  N_MPI(20) = 3
  N_MPI(1)  = 2; N_MPI(2) = 7; N_MPI(3) = 8; iwx=1; iwy=0; iwz=0; CALL BC_MPI(2,1)
  N_MPI(1)  = 3; N_MPI(2) = 6; N_MPI(3) = 8; iwx=0; iwy=1; iwz=0; CALL BC_MPI(2,1)
  N_MPI(1)  = 4; N_MPI(2) = 6; N_MPI(3) = 7; iwx=0; iwy=0; iwz=1; CALL BC_MPI(2,1)
  do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx+1
    Vfc(i,j,k,1) = ( dx(i-1)*U(i,j,k,2)+dx(i)*U(i-1,j,k,2) )/( dx(i-1)+dx(i) )
  end do; end do; end do
  do k = 1, Ncellz; do j = 1, Ncelly+1; do i = 1, Ncellx
    Vfc(i,j,k,2) = ( dy(j-1)*U(i,j,k,3)+dy(j)*U(i,j-1,k,3) )/( dy(j-1)+dy(j) )
  end do; end do; end do
  do k = 1, Ncellz+1; do j = 1, Ncelly; do i = 1, Ncellx
    Vfc(i,j,k,3) = ( dz(k-1)*U(i,j,k,4)+dz(k)*U(i,j,k-1,4) )/( dz(k-1)+dz(k) )
  end do; end do; end do
  N_MPI(13) = 1; iwx=1; iwy=1; iwz=1; CALL BC_MPI_OT(2,1); N_MPI(13) = 0

  N_MPI(20) = 1; N_MPI(1)  = 5; iwx=1; iwy=1; iwz=1; CALL BC_MPI(2,1) !for VLIMIT depends on T
end if

if(mode.eq.2) then !CSD
  if(iwx.eq.1) then
    N_MPI(20) = 1; N_MPI(1)  = 1; iwx=1; iwy=1; iwz=1; CALL BC_MPI(2,1); iwy=0; iwz=0
    do k = 1, Ncellz; do j = 1, Ncelly+1; do i = 1, Ncellx+1
      dnc(i,j,k) = 0.25d0*( U(i,j,k,1)+U(i-1,j,k,1)+U(i,j-1,k,1)+U(i-1,j-1,k,1) )
    end do; end do; end do
  end if
  if(iwy.eq.1) then
    do k = 1, Ncellz+1; do j = 1, Ncelly+1; do i = 1, Ncellx
      dnc(i,j,k) = 0.25d0*( U(i,j,k,1)+U(i,j-1,k,1)+U(i,j,k-1,1)+U(i,j-1,k-1,1) )
    end do; end do; end do
  end if
  if(iwz.eq.1) then
    do k = 1, Ncellz+1; do j = 1, Ncelly; do i = 1, Ncellx+1
      dnc(i,j,k) = 0.25d0*( U(i,j,k,1)+U(i-1,j,k,1)+U(i,j,k-1,1)+U(i-1,j,k-1,1) )
    end do; end do; end do
  end if
end if

if(mode.eq.3) then !update FCMF
  do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx+1
    U(i,j,k,6) = U(i,j,k,6) + dt/dy(j)*( EMF(i,j+1,k,3) - EMF(i,j,k,3) ) - dt/dz(k)*( EMF(i,j,k+1,2) - EMF(i,j,k,2) )
  end do; end do; end do
  do k = 1, Ncellz; do j = 1, Ncelly+1; do i = 1, Ncellx
    U(i,j,k,7) = U(i,j,k,7) + dt/dz(k)*( EMF(i,j,k+1,1) - EMF(i,j,k,1) ) - dt/dx(i)*( EMF(i+1,j,k,3) - EMF(i,j,k,3) )
  end do; end do; end do
  do k = 1, Ncellz+1; do j = 1, Ncelly; do i = 1, Ncellx
    U(i,j,k,8) = U(i,j,k,8) + dt/dx(i)*( EMF(i+1,j,k,2) - EMF(i,j,k,2) ) - dt/dy(j)*( EMF(i,j+1,k,1) - EMF(i,j,k,1) )
  end do; end do; end do
end if

if(mode.eq.4) then !CCMF
  ALLOCATE(Bcc(-1:ndx,-1:ndy,-1:ndz,3))
  do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
    Bcc(i,j,k,1) = 0.5d0 * ( U(i,j,k,6) + U(i+1,j,k,6) )
    Bcc(i,j,k,2) = 0.5d0 * ( U(i,j,k,7) + U(i,j+1,k,7) )
    Bcc(i,j,k,3) = 0.5d0 * ( U(i,j,k,8) + U(i,j,k+1,8) )
  end do; end do; end do
  N_MPI(10) = 0; N_MPI(11) = 1; N_MPI(12) = 0; N_MPI(13) = 0
  iwx=1; iwy=1; iwz=1; CALL BC_MPI_OT(2,1)
  N_MPI(11) = 0
end if

END SUBROUTINE CC


SUBROUTINE CCT(dxe,dye,dt)
USE comvar
USE mpivar

double precision  ::  dxe(-1:ndmax),dye(-1:ndmax),dt
double precision  ::  grdQ(-1:ndmax,16)

integer ::  Vx,Vy,Bx,By,Ez,ifW(4,2)
double precision  ::  C1p,C1m,C2p,C2m,W(4,2),v1a,v2a,B1a,B2a,facm
double precision  ::  C1pLL,C1pRL,C1pLR,C1pRR,C2mLL,C2mRL,C2mLR,C2mRR
double precision  ::  C1mLL,C1mRL,C1mLR,C1mRR,C2pLL,C2pRL,C2pLR,C2pRR
double precision  ::  a1p,a1m,a2p,a2m,da1p,da1m,da2p,da2m,S1,S2
double precision  ::  dQ(0:ndx,0:ndy,0:ndz,4),dsqinv
double precision  ::  Wf1p,Wf1m,Wf2p,Wf2m,safe

facm = facdep*0.5d0

if(iwx.eq.1) then
  Vx = 1; Vy = 2; Bx = 6; By = 7; Ez = 3
  Ncx = Ncellx; Ncy = Ncelly; Ncz = Ncellz
  do k = 1, Ncellz
  do j = 1, Ncelly+1
    call VLIMIT(15,j,k,grdQ,dxe,0,1) !vy
    call VLIMIT(7 ,j,k,grdQ,dxe,0,1) !By
    do i = 0, Ncellx+1
      dQ(i,j,k,2) = grdQ(i,15) !save as vy
      dQ(i,j,k,4) = grdQ(i,7)  !save as By
  end do; end do; end do
  iwx = 0; iwy = 1
  do k = 1, Ncellz
  do i = 1, Ncellx+1
    call VLIMIT(14,k,i,grdQ,dye,0,1) !vx
    call VLIMIT(6 ,k,i,grdQ,dye,0,1) !Bx
    do j = 0, Ncelly+1
      dQ(i,j,k,1) = grdQ(j,14) !save as vx
      dQ(i,j,k,3) = grdQ(j,6)  !save as Bx
  end do; end do; end do
  iwx = 1; iwy = 0
end if
if(iwy.eq.1) then
  Vx = 2; Vy = 3; Bx = 7; By = 8; Ez = 1
  Ncx = Ncelly; Ncy = Ncellz; Ncz = Ncellx
  do i = 1, Ncellx
  do k = 1, Ncellz+1
    call VLIMIT(16,k,i,grdQ,dxe,0,1) !vz
    call VLIMIT(8 ,k,i,grdQ,dxe,0,1) !Bz
    do j = 0, Ncelly+1
      dQ(i,j,k,2) = grdQ(j,16) !save as vy
      dQ(i,j,k,4) = grdQ(j,8)  !save as By
  end do; end do; end do
  iwy = 0; iwz = 1
  do i = 1, Ncellx
  do j = 1, Ncelly+1
    call VLIMIT(15,i,j,grdQ,dye,0,1) !vy
    call VLIMIT(7 ,i,j,grdQ,dye,0,1) !By
    do k = 0, Ncellz+1
      dQ(i,j,k,1) = grdQ(k,15) !save as vx
      dQ(i,j,k,3) = grdQ(k,7)  !save as Bx
  end do; end do; end do
  iwy = 1; iwz = 0
end if
if(iwz.eq.1) then
  Vx = 3; Vy = 1; Bx = 8; By = 6; Ez = 2
  Ncx = Ncellz; Ncy = Ncellx; Ncz = Ncelly
  do j = 1, Ncelly
  do i = 1, Ncellx+1
    call VLIMIT(14,i,j,grdQ,dxe,0,1) !vx
    call VLIMIT(6 ,i,j,grdQ,dxe,0,1) !Bx
    do k = 0, Ncellz+1
      dQ(i,j,k,2) = grdQ(k,14) !save as vy
      dQ(i,j,k,4) = grdQ(k,6)  !save as By
  end do; end do; end do
  iwz = 0; iwx = 1
  do j = 1, Ncelly
  do k = 1, Ncellz+1
    call VLIMIT(16,j,k,grdQ,dye,0,1) !vz
    call VLIMIT(8 ,j,k,grdQ,dye,0,1) !Bz
    do i = 0, Ncellx+1
      dQ(i,j,k,1) = grdQ(i,16) !save as vx
      dQ(i,j,k,3) = grdQ(i,8)  !save as Bx
  end do; end do; end do
  iwz = 1; iwx = 0
end if

!***** Calculate 4 candidates *****
!
!   z   y     !   x   z     !   y   x     
!   |  /      !   |  /      !   |  /      
!   | /       !   | /       !   | /       
!   |/        !   |/        !   |/        
!   ------>x  !   ------>y  !   ------>z  

do k = 1, Ncz
do j = 1, Ncy+1
do i = 1, Ncx+1
   ix0   = iwx*i     + iwy*k     + iwz*j
   jy0   = iwx*j     + iwy*i     + iwz*k
   kz0   = iwx*k     + iwy*j     + iwz*i
 dsqinv = 1.d0/dsqrt(dnc(ix0,jy0,kz0))

!!!!! for C1 > 0, C2 > 0
   S1   =  1.d0; S2   =  1.d0
   ix1   = iwx*i     + iwy*k     + iwz*(j-1)
   jy1   = iwx*(j-1) + iwy*i     + iwz*k
   kz1   = iwx*k     + iwy*(j-1) + iwz*i
   a1p  =  Vfc(ix1,jy1,kz1,Vx) +  U(ix1,jy1,kz1,Bx)*dsqinv
  da1p  = ( dQ(ix1,jy1,kz1,1 ) + dQ(ix1,jy1,kz1,3 )*dsqinv)*facm
   a1m  =  Vfc(ix1,jy1,kz1,Vx) -  U(ix1,jy1,kz1,Bx)*dsqinv
  da1m  = ( dQ(ix1,jy1,kz1,1 ) - dQ(ix1,jy1,kz1,3 )*dsqinv)*facm
   ix2   = iwx*(i-1) + iwy*k     + iwz*j
   jy2   = iwx*j     + iwy*(i-1) + iwz*k
   kz2   = iwx*k     + iwy*j     + iwz*(i-1)
   a2m  =  Vfc(ix2,jy2,kz2,Vy) -  U(ix2,jy2,kz2,By)*dsqinv
  da2m  = ( dQ(ix2,jy2,kz2,2 ) - dQ(ix2,jy2,kz2,4 )*dsqinv)*facm
   a2p  =  Vfc(ix2,jy2,kz2,Vy) +  U(ix2,jy2,kz2,By)*dsqinv
  da2p  = ( dQ(ix2,jy2,kz2,2 ) + dQ(ix2,jy2,kz2,4 )*dsqinv)*facm

  C1pLL = a1p+da1p*dye(j-1)*S1-da1p*(a2m+da2m*dxe(i-1)*S2)*dt
  C1pLL = C1pLL/(1.d0-da1p*da2m*dt*dt)
  C2mLL = a2m+da2m*dxe(i-1)*S2-da2m*(a1p+da1p*dye(j-1)*S1)*dt
  C2mLL = C2mLL/(1.d0-da1p*da2m*dt*dt)        
  C1mLL = a1m+da1m*dye(j-1)*S1-da1m*(a2p+da2p*dxe(i-1)*S2)*dt
  C1mLL = C1mLL/(1.d0-da1m*da2p*dt*dt)
  C2pLL = a2p+da2p*dxe(i-1)*S2-da2p*(a1m+da1m*dye(j-1)*S1)*dt
  C2pLL = C2pLL/(1.d0-da1m*da2p*dt*dt)

!!!!! for C1 > 0, C2 < 0
   S1   = -1.d0; S2   =  1.d0
   ix3   = iwx*i     + iwy*k     + iwz*j
   jy3   = iwx*j     + iwy*i     + iwz*k
   kz3   = iwx*k     + iwy*j     + iwz*i
   a1p  =  Vfc(ix3,jy3,kz3,Vx) +  U(ix3,jy3,kz3,Bx)*dsqinv
  da1p  = ( dQ(ix3,jy3,kz3,1 ) + dQ(ix3,jy3,kz3,3 )*dsqinv)*facm
   a1m  =  Vfc(ix3,jy3,kz3,Vx) -  U(ix3,jy3,kz3,Bx)*dsqinv
  da1m  = ( dQ(ix3,jy3,kz3,1 ) - dQ(ix3,jy3,kz3,3 )*dsqinv)*facm
   ix4   = iwx*(i-1) + iwy*k     + iwz*j
   jy4   = iwx*j     + iwy*(i-1) + iwz*k
   kz4   = iwx*k     + iwy*j     + iwz*(i-1)
   a2m  =  Vfc(ix4,jy4,kz4,Vy) -  U(ix4,jy4,kz4,By)*dsqinv
  da2m  = ( dQ(ix4,jy4,kz4,2 ) - dQ(ix4,jy4,kz4,4 )*dsqinv)*facm
   a2p  =  Vfc(ix4,jy4,kz4,Vy) +  U(ix4,jy4,kz4,By)*dsqinv
  da2p  = ( dQ(ix4,jy4,kz4,2 ) + dQ(ix4,jy4,kz4,4 )*dsqinv)*facm

  C1pRL = a1p+da1p*dye(j  )*S1-da1p*(a2m+da2m*dxe(i-1)*S2)*dt
  C1pRL = C1pRL/(1.d0-da1p*da2m*dt*dt)
  C2mRL = a2m+da2m*dxe(i-1)*S2-da2m*(a1p+da1p*dye(j  )*S1)*dt
  C2mRL = C2mRL/(1.d0-da1p*da2m*dt*dt)
  C1mRL = a1m+da1m*dye(j  )*S1-da1m*(a2p+da2p*dxe(i-1)*S2)*dt
  C1mRL = C1mRL/(1.d0-da1m*da2p*dt*dt)
  C2pRL = a2p+da2p*dxe(i-1)*S2-da2p*(a1m+da1m*dye(j  )*S1)*dt
  C2pRL = C2pRL/(1.d0-da1m*da2p*dt*dt)

!!!!! for C1 < 0, C2 > 0
   S1   =  1.d0; S2   = -1.d0
   ix5   = iwx*i     + iwy*k     + iwz*(j-1)
   jy5   = iwx*(j-1) + iwy*i     + iwz*k
   kz5   = iwx*k     + iwy*(j-1) + iwz*i
   a1p  =  Vfc(ix5,jy5,kz5,Vx) +  U(ix5,jy5,kz5,Bx)*dsqinv
  da1p  = ( dQ(ix5,jy5,kz5,1 ) + dQ(ix5,jy5,kz5,3 )*dsqinv)*facm
   a1m  =  Vfc(ix5,jy5,kz5,Vx) -  U(ix5,jy5,kz5,Bx)*dsqinv
  da1m  = ( dQ(ix5,jy5,kz5,1 ) - dQ(ix5,jy5,kz5,3 )*dsqinv)*facm
   ix6   = iwx*i     + iwy*k     + iwz*j
   jy6   = iwx*j     + iwy*i     + iwz*k
   kz6   = iwx*k     + iwy*j     + iwz*i
   a2m  =  Vfc(ix6,jy6,kz6,Vy) -  U(ix6,jy6,kz6,By)*dsqinv
  da2m  = ( dQ(ix6,jy6,kz6,2 ) - dQ(ix6,jy6,kz6,4 )*dsqinv)*facm
   a2p  =  Vfc(ix6,jy6,kz6,Vy) +  U(ix6,jy6,kz6,By)*dsqinv
  da2p  = ( dQ(ix6,jy6,kz6,2 ) + dQ(ix6,jy6,kz6,4 )*dsqinv)*facm

  C1pLR = a1p+da1p*dye(j-1)*S1-da1p*(a2m+da2m*dxe(i  )*S2)*dt
  C1pLR = C1pLR/(1.d0-da1p*da2m*dt*dt)
  C2mLR = a2m+da2m*dxe(i  )*S2-da2m*(a1p+da1p*dye(j-1)*S1)*dt
  C2mLR = C2mLR/(1.d0-da1p*da2m*dt*dt)
  C1mLR = a1m+da1m*dye(j-1)*S1-da1m*(a2p+da2p*dxe(i  )*S2)*dt
  C1mLR = C1mLR/(1.d0-da1m*da2p*dt*dt)
  C2pLR = a2p+da2p*dxe(i  )*S2-da2p*(a1m+da1m*dye(j-1)*S1)*dt
  C2pLR = C2pLR/(1.d0-da1m*da2p*dt*dt)

!!!!! for C1 < 0, C2 < 0
   S1   = -1.d0; S2   = -1.d0
   ix7   = iwx*i     + iwy*k     + iwz*j
   jy7   = iwx*j     + iwy*i     + iwz*k
   kz7   = iwx*k     + iwy*j     + iwz*i
   a1p  =  Vfc(ix7,jy7,kz7,Vx) +  U(ix7,jy7,kz7,Bx)*dsqinv
  da1p  = ( dQ(ix7,jy7,kz7,1 ) + dQ(ix7,jy7,kz7,3 )*dsqinv)*facm
   a1m  =  Vfc(ix7,jy7,kz7,Vx) -  U(ix7,jy7,kz7,Bx)*dsqinv
  da1m  = ( dQ(ix7,jy7,kz7,1 ) - dQ(ix7,jy7,kz7,3 )*dsqinv)*facm
   a2m  =  Vfc(ix7,jy7,kz7,Vy) -  U(ix7,jy7,kz7,By)*dsqinv
  da2m  = ( dQ(ix7,jy7,kz7,2 ) - dQ(ix7,jy7,kz7,4 )*dsqinv)*facm
   a2p  =  Vfc(ix7,jy7,kz7,Vy) +  U(ix7,jy7,kz7,By)*dsqinv
  da2p  = ( dQ(ix7,jy7,kz7,2 ) + dQ(ix7,jy7,kz7,4 )*dsqinv)*facm

  C1pRR = a1p+da1p*dye(j  )*S1-da1p*(a2m+da2m*dxe(i  )*S2)*dt
  C1pRR = C1pRR/(1.d0-da1p*da2m*dt*dt)
  C2mRR = a2m+da2m*dxe(i  )*S2-da2m*(a1p+da1p*dye(j  )*S1)*dt
  C2mRR = C2mRR/(1.d0-da1p*da2m*dt*dt)
  C1mRR = a1m+da1m*dye(j  )*S1-da1m*(a2p+da2p*dxe(i  )*S2)*dt
  C1mRR = C1mRR/(1.d0-da1m*da2p*dt*dt)
  C2pRR = a2p+da2p*dxe(i  )*S2-da2p*(a1m+da1m*dye(j  )*S1)*dt
  C2pRR = C2pRR/(1.d0-da1m*da2p*dt*dt)

!***** Weightng Factor for Upwinding *****
  ifW(1,1) = 0; ifW(2,1) = 0; ifW(3,1) = 0; ifW(4,1) = 0; ifW(1,2) = 0; ifW(2,2) = 0; ifW(3,2) = 0; ifW(4,2) = 0
  W(1,1) = 0.d0; W(2,1) = 0.d0; W(3,1) = 0.d0; W(4,1) = 0.d0; W(1,2) = 0.d0; W(2,2) = 0.d0; W(3,2) = 0.d0; W(4,2) = 0.d0;
  
  if((C1pLL.gt.0.d0).and.(C2mLL.gt.0.d0)) ifW(1,1)=1
!  ifW(1,1) = int( 0.51d0-dsign(0.25d0,-C1pLL)-dsign(0.25d0,-C2mLL) ) !C1pLL>0,C2mLL>0
  if((C1mLL.gt.0.d0).and.(C2pLL.gt.0.d0)) ifW(1,2)=1
!  ifW(1,2) = int( 0.51d0-dsign(0.25d0,-C1mLL)-dsign(0.25d0,-C2pLL) ) !C1mLL>0,C2pLL>0
  if((C1pRL.gt.0.d0).and.(C2mRL.lt.0.d0)) ifW(2,1)=1
!  ifW(2,1) = int( 0.51d0-dsign(0.25d0,-C1pRL)-dsign(0.25d0, C2mRL) ) !C1pRL>0,C2mRL<0
  if((C1mRL.gt.0.d0).and.(C2pRL.lt.0.d0)) ifW(2,2)=1
!  ifW(2,2) = int( 0.51d0-dsign(0.25d0,-C1mRL)-dsign(0.25d0, C2pRL) ) !C1mRL>0,C2pRL<0
  if((C1pLR.lt.0.d0).and.(C2mLR.gt.0.d0)) ifW(3,1)=1
!  ifW(3,1) = int( 0.51d0-dsign(0.25d0, C1pLR)-dsign(0.25d0,-C2mLR) ) !C1pLR<0,C2mLR>0
  if((C1mLR.lt.0.d0).and.(C2pLR.gt.0.d0)) ifW(3,2)=1
!  ifW(3,2) = int( 0.51d0-dsign(0.25d0, C1mLR)-dsign(0.25d0,-C2pLR) ) !C1mLR<0,C2pLR>0
  if((C1pRR.lt.0.d0).and.(C2mRR.lt.0.d0)) ifW(4,1)=1
!  ifW(4,1) = int( 0.51d0-dsign(0.25d0, C1pRR)-dsign(0.25d0, C2mRR) ) !C1pRR<0,C2mRR<0
  if((C1mRR.lt.0.d0).and.(C2pRR.lt.0.d0)) ifW(4,2)=1
!  ifW(4,2) = int( 0.51d0-dsign(0.25d0, C1mRR)-dsign(0.25d0, C2pRR) ) !C1mRR<0,C2pRR<0
  
  W(1,1) = dble(ifW(1,1)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(1,1)+0.51d0*ifW(1,2)))
  W(1,2) = dble(ifW(1,2)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(1,1)+0.51d0*ifW(1,2)))
  W(2,1) = dble(ifW(2,1)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(2,1)+0.51d0*ifW(2,2)))
  W(2,2) = dble(ifW(2,2)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(2,1)+0.51d0*ifW(2,2)))
  W(3,1) = dble(ifW(3,1)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(3,1)+0.51d0*ifW(3,2)))
  W(3,2) = dble(ifW(3,2)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(3,1)+0.51d0*ifW(3,2)))
  W(4,1) = dble(ifW(4,1)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(4,1)+0.51d0*ifW(4,2)))
  W(4,2) = dble(ifW(4,2)) + (1.d13-1.d0)*dble(int(0.51d0*ifW(4,1)+0.51d0*ifW(4,2)))
  
!***** Calculate Weighted Characteristics *****

  C1p = 0.d0; C2m = 0.d0
  if(ifW(1,1)+ifW(2,1)+ifW(3,1)+ifW(4,1).eq.0) goto 3333
  C1p = W(1,1)*dabs(C2mLL)*C1pLL+W(2,1)*dabs(C2mRL)*C1pRL+W(3,1)*dabs(C2mLR)*C1pLR+W(4,1)*dabs(C2mRR)*C1pRR
  C1p = C1p/(  W(1,1)*dabs(C2mLL)+W(2,1)*dabs(C2mRL)+W(3,1)*dabs(C2mLR)+W(4,1)*dabs(C2mRR)  )
  C2m = W(1,1)*dabs(C1pLL)*C2mLL+W(2,1)*dabs(C1pRL)*C2mRL+W(3,1)*dabs(C1pLR)*C2mLR+W(4,1)*dabs(C1pRR)*C2mRR
  C2m = C2m/(  W(1,1)*dabs(C1pLL)+W(2,1)*dabs(C1pRL)+W(3,1)*dabs(C1pLR)+W(4,1)*dabs(C1pRR)  )
  3333  continue
 
  C1m = 0.d0; C2p = 0.d0
  if(ifW(1,2)+ifW(2,2)+ifW(3,2)+ifW(4,2).eq.0) goto 3334
  C1m = W(1,2)*dabs(C2pLL)*C1mLL+W(2,2)*dabs(C2pRL)*C1mRL+W(3,2)*dabs(C2pLR)*C1mLR+W(4,2)*dabs(C2pRR)*C1mRR
  C1m = C1m/(  W(1,2)*dabs(C2pLL)+W(2,2)*dabs(C2pRL)+W(3,2)*dabs(C2pLR)+W(4,2)*dabs(C2pRR)  )
  C2p = W(1,2)*dabs(C1mLL)*C2pLL+W(2,2)*dabs(C1mRL)*C2pRL+W(3,2)*dabs(C1mLR)*C2pLR+W(4,2)*dabs(C1mRR)*C2pRR
  C2p = C2p/(  W(1,2)*dabs(C1mLL)+W(2,2)*dabs(C1mRL)+W(3,2)*dabs(C1mLR)+W(4,2)*dabs(C1mRR)  )
  3334  continue
  
!***** contain EMF *****
  v1a = C1p+C1m
  v2a = C2p+C2m
  B1a = C1p-C1m
  B2a = C2p-C2m
  
  EMF(ix0,jy0,kz0,Ez) = 0.25d0*( v1a*B2a-v2a*B1a )*dsqrt(dnc(ix0,jy0,kz0))
  
end do; end do; end do

END SUBROUTINE CCT


SUBROUTINE SOURCE(dt)
USE comvar
USE mpivar
USE chmvar
INCLUDE 'mpif.h'

double precision  dt
DOUBLE PRECISION :: ndpold,ndHold,ndH2old,ndHeold,ndHepold,ndCold,ndCpold,ndCOold,T
DOUBLE PRECISION :: zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
DOUBLE PRECISION :: temp1,temp2,temp3,omeps,eps
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Tn,Pn,Qx,Qy,Qz
double precision  :: mmean,rtTx,rtTy,rtTz,tcd,CooL

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
  nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
  ndtot(i,j,k) = ndp(i,j,k)+ndH(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
end do; end do; end do

if(ifrad.eq.2) then
  call SHIELD()
end if

!*******************************(A2)

!*******************************(B2)

!*******************************(B3)

!***** 2nd order explicit scheme ( with cooling function ) *****
if(ifthrm.eq.2) then

ALLOCATE(Tn(0:ndx,0:ndy,0:ndz),Pn(0:ndx,0:ndy,0:ndz),Qx(0:ndx,0:ndy,0:ndz),Qy(0:ndx,0:ndy,0:ndz),Qz(0:ndx,0:ndy,0:ndz))

N_MPI(20) = 2
N_MPI(1)  = 1
N_MPI(2)  = 5
iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)

N_MPI(10) = 10; CALL BC_MPI_OT(1,1); N_MPI(10) = 0 ! for chemical boundary

do k = 0, Ncellz+1; do j = 0, Ncelly+1; do i = 0, Ncellx+1
  Tn(i,j,k) = U(i,j,k,5)/( kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)) )
  Pn(i,j,k) = U(i,j,k,5)
end do; end do; end do

do k = 0, Ncellz; do j = 0, Ncelly; do i = 0, Ncellx
  rtTx    = dsqrt(  ( dx(i)*Tn(i+1,j,k)+dx(i+1)*Tn(i,j,k) ) /( dx(i) + dx(i+1) )  )
  rtTy    = dsqrt(  ( dy(j)*Tn(i,j+1,k)+dy(j+1)*Tn(i,j,k) ) /( dy(j) + dy(j+1) )  )
  rtTz    = dsqrt(  ( dz(k)*Tn(i,j,k+1)+dz(k+1)*Tn(i,j,k) ) /( dz(k) + dz(k+1) )  )
  Qx(i,j,k) = Kcond * dmin1(rtTx,3.873d0) * ( Tn(i+1,j,k)-Tn(i,j,k) )*2.d0 /( dx(i) + dx(i+1) )
  Qy(i,j,k) = Kcond * dmin1(rtTy,3.873d0) * ( Tn(i,j+1,k)-Tn(i,j,k) )*2.d0 /( dy(j) + dy(j+1) )
  Qz(i,j,k) = Kcond * dmin1(rtTz,3.873d0) * ( Tn(i,j,k+1)-Tn(i,j,k) )*2.d0 /( dz(k) + dz(k+1) )
end do; end do; end do

!if(NRANK==40) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),'point100'

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
!----- Cooling ---------------------------------------------------------
  Call Fcool( CooL,Tn(i,j,k),i,j,k)
  gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1
  if( dt .le. 0.2d0*Pn(i,j,k)/(gammi1*dabs(CooL)) ) then
    U(i,j,k,5) = U(i,j,k,5) - gammi1*CooL*dt*0.5d0 !explicit
  else
    Call IMC( U(i,j,k,5),ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)+ndH2(i,j,k),dt*0.5d0,i,j,k ) !implicit
  end if
!----- Conduction ------------------------------------------------------
  tcd = gammi1*Kcond*dsqrt( Tn(i,j,k) )
  tcd = 0.49d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)+ndH2(i,j,k))*kb/tcd*(dx(i)**2.d0)
  tcd = dsign(0.5d0,tcd-dt)
  U(i,j,k,5) = U(i,j,k,5) + (tcd+0.5d0)*gammi1*dt*0.5d0* &
  ( (Qx(i,j,k)-Qx(i-1,j,k))/dx(i) + (Qy(i,j,k)-Qy(i,j-1,k))/dy(j) + (Qz(i,j,k)-Qz(i,j,k-1))/dz(k) )
  U(i,j,k,5) = dmax1(U(i,j,k,5),pmin)
  U(i,j,k,5) = dmin1(U(i,j,k,5),pmax)
end do; end do; end do

!if(NRANK==40) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),'point101'

N_MPI(20) = 1
N_MPI(1)  = 5
iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)

do k = 0, Ncellz+1; do j = 0, Ncelly+1; do i = 0, Ncellx+1
  Tn(i,j,k) = U(i,j,k,5)/( kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)) )
end do; end do; end do

!if(NRANK==40) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),Tn(1,1,1),'point101a'

do k = 0, Ncellz; do j = 0, Ncelly; do i = 0, Ncellx
  rtTx    = dsqrt(  ( dx(i)*Tn(i+1,j,k)+dx(i+1)*Tn(i,j,k) ) /( dx(i) + dx(i+1) )  )
  rtTy    = dsqrt(  ( dy(j)*Tn(i,j+1,k)+dy(j+1)*Tn(i,j,k) ) /( dy(j) + dy(j+1) )  )
  rtTz    = dsqrt(  ( dz(k)*Tn(i,j,k+1)+dz(k+1)*Tn(i,j,k) ) /( dz(k) + dz(k+1) )  )
  Qx(i,j,k) = Kcond * dmin1(rtTx,3.873d0) * ( Tn(i+1,j,k)-Tn(i,j,k) )*2.d0 /( dx(i) + dx(i+1) )
  Qy(i,j,k) = Kcond * dmin1(rtTy,3.873d0) * ( Tn(i,j+1,k)-Tn(i,j,k) )*2.d0 /( dy(j) + dy(j+1) )
  Qz(i,j,k) = Kcond * dmin1(rtTz,3.873d0) * ( Tn(i,j,k+1)-Tn(i,j,k) )*2.d0 /( dz(k) + dz(k+1) )
end do; end do; end do

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
!----- Cooling ---------------------------------------------------------
  Call Fcool( CooL,Tn(i,j,k),i,j,k)
  gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1
  if( dt .le. 0.2d0*Pn(i,j,k)/(gammi1*dabs(CooL)) ) then
    U(i,j,k,5) = Pn(i,j,k) - gammi1*CooL*dt !explicit
  else
    Call IMC( Pn(i,j,k),ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)+ndH2(i,j,k),dt,i,j,k ) !implicit
    U(i,j,k,5) = Pn(i,j,k)
  end if
!if(NRANK==40 .and. k==1 .and. j==1 .and. i==1) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),Tn(1,1,1),'point101b'
!----- Conduction ------------------------------------------------------
  tcd = gammi1*Kcond*dsqrt( Tn(i,j,k) )
  tcd = 0.49d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)+ndH2(i,j,k))*kb/tcd*(dx(i)**2.d0)
  tcd = dsign(0.5d0,tcd-dt)
  U(i,j,k,5) = U(i,j,k,5) + (tcd+0.5d0)*gammi1*dt* &
  ( (Qx(i,j,k)-Qx(i-1,j,k))/dx(i) + (Qy(i,j,k)-Qy(i,j-1,k))/dy(j) + (Qz(i,j,k)-Qz(i,j,k-1))/dz(k) )
!if(NRANK==40 .and. k==1 .and. j==1 .and. i==1) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),Tn(1,1,1),'point101c'
  U(i,j,k,5) = dmax1(U(i,j,k,5),pmin)
!if(NRANK==40 .and. k==1 .and. j==1 .and. i==1) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),Tn(1,1,1),'point101d'
  U(i,j,k,5) = dmin1(U(i,j,k,5),pmax)
!if(NRANK==40 .and. k==1 .and. j==1 .and. i==1) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,5),sngl(U(33,33,33,1)),Tn(1,1,1),'point101e'
end do; end do; end do

DEALLOCATE(Tn,Pn,Qx,Qy,Qz)

end if
!*******************************(C3)
!if(NRANK==40) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)),'point102'
!*******************************(A1)
if(ifchem.eq.2) then


do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103am'

  ndpold=ndp(i,j,k); ndHold=ndH(i,j,k); ndH2old=ndH2(i,j,k); ndHeold=ndHe(i,j,k)
  ndHepold=ndHep(i,j,k); ndCold=ndC(i,j,k); ndCpold=ndCp(i,j,k); ndCOold=ndCO(i,j,k)
  T = U(i,j,k,5)/kb/( ndpold+ndHold+ndH2old+ndHeold+ndHepold )
  call RATES(i,j,k,T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)
  
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103al'

! H recombination & ionization by CR
  temp1 = kHrec*nde(i,j,k)
  temp2 = dexp(-dt*(zeta+temp1))
  temp3 = 1.d0/(zeta+temp1)
  call Omexp(omeps,dt*(zeta+temp1))
  ndH(i,j,k) = ( temp1*ndHold + zeta*ndHold*temp2 + temp1*ndpold*omeps )*temp3
  ndp(i,j,k) = (  zeta*ndpold +temp1*ndpold*temp2 +  zeta*ndHold*omeps )*temp3
  ndHold = ndH(i,j,k); ndpold = ndp(i,j,k); nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ak'

!He recombination & ionization by CR
  temp1 = kHerec*nde(i,j,k)
  temp2 = dexp(-dt*(zeta+temp1))
  temp3 = 1.d0/(zeta+temp1)
  call Omexp(omeps,dt*(zeta+temp1))
  ndHe(i,j,k)  = ( temp1*ndHeold + zeta*ndHeold*temp2 + temp1*ndHepold*omeps )*temp3
  ndHep(i,j,k) = ( zeta*ndHepold +temp1*ndHepold*temp2+  zeta*ndHeold*omeps )*temp3
  ndHeold = ndHe(i,j,k); ndHepold = ndHep(i,j,k); nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103aj'

!H2 formation & dissociation by UV & electron collision
  temp1 = kH2ph + kH2de*nde(i,j,k)
  temp2 = temp1 + 2.d0*kH2*ndtot(i,j,k)
  temp3 = dexp(-dt*temp2)
  call Omexp(omeps,dt*temp2)
  ndH(i,j,k) = ( ndHold+2.d0*ndH2old*omeps )*temp1 + 2.d0*kH2*ndHold*ndtot(i,j,k)*temp3
  ndH(i,j,k) = ndH(i,j,k)/temp2
  ndH2(i,j,k)= ( 2.d0*ndH2old+ndHold*omeps )*kH2*ndtot(i,j,k) + temp1*ndH2old*temp3
  ndH2(i,j,k)= ndH2(i,j,k)/temp2
  ndHold = ndH(i,j,k); ndH2old = ndH2(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ai'
  
!H2 dissociation by H collision
  temp1 = ndHold + 2.d0*ndH2old
  temp2 = dexp(-dt*temp1*kH2dH)
  temp3 = 1.d0/(ndHold + 2.d0*ndH2old*temp2)
  ndH(i,j,k) = ndHold*temp1*temp3
  ndH2(i,j,k)= ndH2old*temp1*temp2*temp3
  ndHold = ndH(i,j,k); ndH2old = ndH2(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ah'
  
!CO formation
  temp1 = dexp(-dt*kCO*ndtot(i,j,k))
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,temp1  
eps = dt*kCO*ndtot(i,j,k)
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,eps
 call Omexp(omeps,eps)
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,omeps
  ndCO(i,j,k)  = ndCOold + omeps*ndCpold
  ndCp(i,j,k)  = temp1*ndCpold
  ndCOold = ndCO(i,j,k); ndCpold = ndCp(i,j,k); nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
  
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ag'

!CO dissociation by UV & electron collision & H collision
  temp1 = dexp(-dt*(kCOph+kCOde*nde(i,j,k)+kCOdH*ndHold))
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,temp1,dt,kCOph,kCOde,nde(33,33,33),kCOdH,ndHold
  eps = dt*(kCOph+kCOde*nde(i,j,k)+kCOdH*ndHold)
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,eps
 call Omexp(omeps,eps)
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,omeps,eps
  ndC(i,j,k) = ndCold + omeps*ndCOold
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,ndCold,ndCOold,ndC(33,33,33)
  ndCO(i,j,k)= temp1*ndCOold
  ndCold = ndC(i,j,k); ndCOold = ndCO(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103af'
  
!C recombination & ionization by CR & UV
  temp1 = kCrec*nde(i,j,k)
  temp2 = dexp(-dt*(kCi+temp1))
  temp3 = 1.d0/(kCi+temp1)
  call Omexp(omeps,dt*(kCi+temp1))
  ndC(i,j,k) = ( temp1*ndCold +   kCi*ndCold*temp2 + temp1*ndCpold*omeps )*temp3
  ndCp(i,j,k)= (  kCi*ndCpold + temp1*ndCpold*temp2+   kCi*ndCold*omeps )*temp3
  ndCold = ndC(i,j,k); ndCpold = ndCp(i,j,k); nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ae'

!H, He, C ionization by e collision
   ndH(i,j,k) = ndHold *dexp(-dt*kHie *nde(i,j,k))
  ndHe(i,j,k) = ndHeold*dexp(-dt*kHeie*nde(i,j,k))
   ndC(i,j,k) = ndCold *dexp(-dt*kCie *nde(i,j,k))
  call Omexp(omeps,dt*kHie *nde(i,j,k));  ndp(i,j,k) = omeps*ndHold + ndpold
  call Omexp(omeps,dt*kHeie*nde(i,j,k)); ndHep(i,j,k)= omeps*ndHeold+ ndHepold
  call Omexp(omeps,dt*kCie *nde(i,j,k)); ndCp(i,j,k) = omeps*ndCold + ndCpold
  ndHold = ndH(i,j,k); ndHeold  =  ndHe(i,j,k);  ndCold =  ndC(i,j,k)
  ndpold = ndp(i,j,k); ndHepold = ndHep(i,j,k); ndCpold = ndCp(i,j,k) 
  nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ad'

!H, He, C ionization by H or H2 or p collision
   temp1 = 1.d0/( 1.d0+kHiH*dt*ndHold )
   ndH(i,j,k) = ndHold*temp1
   ndp(i,j,k) = ndpold+temp1*kHiH*dt*ndHold**2
   temp1 = ndH(i,j,k)+ndp(i,j,k)+ndH2old
  ndHe(i,j,k) = ndHeold*dexp(-dt*kHeiH*temp1)
   ndC(i,j,k) = ndCold *dexp(-dt*kCiH *temp1)
  call Omexp(omeps,dt*kHeiH*temp1); ndHep(i,j,k)= omeps*ndHeold+ ndHepold
  call Omexp(omeps,dt*kCiH *temp1); ndCp(i,j,k) = omeps*ndCold + ndCpold
  ndHold = ndH(i,j,k); ndHeold  =  ndHe(i,j,k);  ndCold =  ndC(i,j,k)
  ndpold = ndp(i,j,k); ndHepold = ndHep(i,j,k); ndCpold = ndCp(i,j,k) 
  nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ac'

!H2 dissiciation by Hep recombination
  temp1 = ndHepold-ndH2old
  if(temp1.ge.0.d0) then
    temp3 = dexp(-dt*kH2dHep*temp1)
    temp2 = 1.d0/( ndHepold-ndH2old*temp3 )
    ndHep(i,j,k) = ndHepold*temp1*temp2
     ndH2(i,j,k) =  ndH2old*temp1*temp2*temp3
    temp3 = ndHepold*( 1.d0-temp1*temp2 )
      ndp(i,j,k) =  ndpold + temp3
      ndH(i,j,k) =  ndHold + temp3
     ndHe(i,j,k) = ndHeold + temp3
  else
    temp3 = dexp(dt*kH2dHep*temp1)
    temp2 = 1.d0/( ndHepold*temp3-ndH2old )
    ndHep(i,j,k) = ndHepold*temp1*temp2*temp3
     ndH2(i,j,k) =  ndH2old*temp1*temp2
    temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
      ndp(i,j,k) =  ndpold + temp3
      ndH(i,j,k) =  ndHold + temp3
     ndHe(i,j,k) = ndHeold + temp3
  end if
   ndHold =  ndH(i,j,k);   ndpold =   ndp(i,j,k); ndH2old = ndH2(i,j,k)
  ndHeold = ndHe(i,j,k); ndHepold = ndHep(i,j,k) 

!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103ab'

!CO dissiciation by Hep recombination
  temp1 = ndHepold-ndCOold
  if(temp1.ge.0.d0) then
    temp3 = dexp(-dt*kCOdHep*temp1)
    temp2 = 1.d0/( ndHepold-ndCOold*temp3 )
    ndHep(i,j,k) = ndHepold*temp1*temp2
     ndCO(i,j,k) =  ndCOold*temp1*temp2*temp3
    temp3 = ndHepold*( 1.d0-temp1*temp2 )
     ndCp(i,j,k) = ndCpold + temp3
     ndHe(i,j,k) = ndHeold + temp3
  else
    temp3 = dexp(dt*kCOdHep*temp1)
    temp2 = 1.d0/( ndHepold*temp3-ndCOold )
    ndHep(i,j,k) = ndHepold*temp1*temp2*temp3
     ndCO(i,j,k) =  ndCOold*temp1*temp2
    temp3 = ndHepold*( 1.d0-temp1*temp2*temp3 )
     ndCp(i,j,k) = ndCpold + temp3
     ndHe(i,j,k) = ndHeold + temp3
  end if
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103aa'

    ndH(i,j,k) = dmax1( ndHmin  ,ndH(i,j,k)   )
    ndp(i,j,k) = dmax1( ndpmin  ,ndp(i,j,k)   )
   ndH2(i,j,k) = dmax1( ndH2min ,ndH2(i,j,k)  )
   ndHe(i,j,k) = dmax1( ndHemin ,ndHe(i,j,k)  )
  ndHep(i,j,k) = dmax1( ndHepmin,ndHep(i,j,k) )
    ndC(i,j,k) = dmax1( ndCmin  ,ndC(i,j,k)   )
   ndCp(i,j,k) = dmax1( ndCpmin ,ndCp(i,j,k)  )
   ndCO(i,j,k) = dmax1( ndCOmin ,ndCO(i,j,k)  )
!if(NRANK==40 .and. k==33 .and. j==33 .and. i==33) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)), &
!ndH(33,33,33),ndp(33,33,33),ndH2(33,33,33),ndHe(33,33,33),ndHep(33,33,33),ndC(33,33,33),ndCp(33,33,33),ndCO(33,33,33),'point103a'
  nde(i,j,k) = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
  ndtot(i,j,k) = ndp(i,j,k)+ndH(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
  U(i,j,k,1) = mH*ndp(i,j,k)+mH*ndH(i,j,k)+mH2*ndH2(i,j,k)+mHe*ndHe(i,j,k)+mHe*ndHep(i,j,k)
end do; end do; end do

end if
if(NRANK==40) write(*,*) NRANK,U(33,33,33,1),U(33,33,33,2),U(33,33,33,5),sngl(U(33,33,33,1)),'point103'


END SUBROUTINE SOURCE



SUBROUTINE Stblty(tLMT)
USE comvar
USE mpivar
USE chmvar

double precision  tLMT,alpha,tauC,Nn,Tn,dl!,cs
double precision  CooL

tLMT = tfinal
!cs =

do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
  Nn = ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
  Tn = U(i,j,k,5)/(kb*Nn)
!-----------------------------( Courant condition )
  dl    = dmin1(dx(i),dy(j),dz(k))
  gammi1 =   3.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+5.d0*ndH2(i,j,k)
  gammi1 = ( 2.d0*(ndH(i,j,k)+ndp(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))+2.d0*ndH2(i,j,k) )/gammi1
  alpha = gammi1*Kcond*dsqrt( dmin1(Tn,15.d0) )
  alpha = 0.4d0*Nn*kb/alpha*(dl**2.d0)
!-----------------------------( Avoid over cooling )
  Call Fcool(CooL,Tn,i,j,k)
  tauC  =  U(i,j,k,5)/gammi1/dabs(CooL)
  tauC  =  dmax1( 0.2d0*tauC , 2.5d-4 )


!if(NRANK==40) write(*,*) NRANK,tauC,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point201'

  tLMT =  dmin1( tLMT , tauC  )


!if(NRANK==40) write(*,*) NRANK,tLMT,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point202'

  tLMT =  dmin1( tLMT , alpha )

!if(NRANK==40) write(*,*) NRANK,tLMT,U(33,33,33,1),U(33,33,33,2),sngl(U(33,33,33,1)),'point203'
end do; end do; end do
if(tLMT.lt.0.d0) write(5,*) time,NRANK,'err at Stblty'

END SUBROUTINE Stblty


SUBROUTINE Fcool(CooL,T,i,j,k)
USE comvar
USE chmvar
double precision :: CooL,T,Av1,Av2,x1,x2
double precision :: Laml,Lamc,Lamo,Lamd,LCOr,LCOH,LCOH2,pha,ncr
double precision :: Gampe,Gamcr,Gampd
double precision :: ATN1,ATN2,SHLD1,SHLD2
double precision :: tau1,tau2,ct1,ct2,ym1,ym2,fes
double precision :: tC1,tC2,fesC1,fesC2,tO1,tO2,fesO1,fesO2
double precision :: n1,n2,b21

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

Av1  = 1.63542d-3*Ntot(i,j,k,1); x1 = 6.1714d3*NH2(i,j,k,1)
Av2  = 1.63542d-3*Ntot(i,j,k,2); x2 = 6.1714d3*NH2(i,j,k,2)
ATN1 = ( dexp(-2.5d0*Av1)     + dexp(-2.5d0*Av2) ) * 0.5d0
ATN2 = ( dexp(-3.77358d0*Av1) + dexp(-3.77358d0*Av2) ) * 0.5d0
pha  = G0 * dsqrt(1.d3*T) * ATN1 / nde(i,j,k)
!------------------------- Lya Cooling
Laml = ndH(i,j,k)*nde(i,j,k) * 1.44049d9*dexp( -1.184d2/T)
Laml = Laml + ndp(i,j,k)*nde(i,j,k) * 1.48d6*dexp( -8.d1/T) !ion
!_________________________ CII Cooling L
call fesc(tCII(i,j,k,1),fesC1); call fesc(tCII(i,j,k,2),fesC2); b21 = fesC1+fesC2
call LEVC2(T,b21,ndH(i,j,k),ndH2(i,j,k),nde(i,j,k),ndCp(i,j,k),n1,n2)
Lamc = 6.0157d7*n2*b21
!***------------------------- OI Cooling
tO1  = 3.40057d-6*Ntot(i,j,k,1)/xo ; tO2 = 3.40057d-6*Ntot(i,j,k,2)/xo 
call fesc(tO1,fesO1); call fesc(tO2,fesO2)
Lamo = (dmax1(ndtot(i,j,k)*xo-ndCO(i,j,k),0.d0)/xo) * (ndH(i,j,k)+0.5d0*ndH2(i,j,k)) * 1.23916d1 * (T**0.4d0) * dexp( -0.228d0/T )
Lamo = Lamo * (fesO1+fesO2)
!***------------------------- DustRec Cooling
Lamd = nde(i,j,k)*ndtot(i,j,k)*6.06236d0 * (T**0.94d0) * ( pha**( 0.462628d0/(T**6.8d-2) ) )
!***------------------------- CO Cooling
ncr  = 3.3d6*(T**0.75d0)/(ndH(i,j,k)+ndp(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k))
tau1 = 1.33194d1*NCO(i,j,k,1)/(T*dv); tau2 = 1.33194d1*NCO(i,j,k,2)/(T*dv)
ct1  = tau1*dsqrt( 6.283d0*dlog(2.13d0+(tau1*0.36788d0)**2) ); ct2 = tau2*dsqrt( 6.283d0*dlog(2.13d0+(tau2*0.36788d0)**2) )
ym1  = dlog( 1.d0+ct1/(1.d0+1.d1*ncr) ); ym2  = dlog( 1.d0+ct2/(1.d0+1.d1*ncr) )
fes  = (2.d0+ym1+0.6d0*ym1**2)/(1.d0+ct1+ncr+1.5d0*dsqrt(ncr)) + (2.d0+ym2+0.6d0*ym2**2)/(1.d0+ct2+ncr+1.5d0*dsqrt(ncr))
LCOr = ndCO(i,j,k) * 1.91505d10*(T**2)*fes
LCOH = ndH(i,j,k) *ndCO(i,j,k) * 7.96086d4*dsqrt(T)*dexp(-(2.d0/T)**3.43d0)*dexp(-3.08d0/T)
LCOH2= ndH2(i,j,k)*ndCO(i,j,k) * 3.60834d4*T*dexp(-(3.14d2/T)**0.333d0)*dexp(-3.08d0/T)
!***------------------------- Photo-electric Heating
Gampe = ndtot(i,j,k) * 2.56526d3 * G0*ATN1 * &
       ( 7.382d-3*(T**0.7d0)/(1.d0+2.d-4*pha) + 4.9d-2/(1.d0+4.d-3*(pha**0.73d0)) )
!------------------------- CR Heating
Gamcr = (ndH(i,j,k)+ndHe(i,j,k)+ndH2(i,j,k)) * 1.89435d0
!***------------------------- Photo-destruction Heating
SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))
Gampd = ndH2(i,j,k) * 4.16362d4 * G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0

CooL  = Laml + Lamc + Lamo + Lamd + LCOr + LCOH + LCOH2 - Gampe - Gamcr - Gampd

END SUBROUTINE Fcool

SUBROUTINE fesc(tau,fes)
double precision tau,fes
fes =   (0.5d0+dsign(0.5d0,0.1d0-1.d-16-tau))*(0.5d0-0.585d0*tau+0.4563d0*tau**2) &
      + (0.5d0+dsign(0.5d0,11.9025d0-(tau-3.55d0)**2))*(1.d0-dexp(-2.34d0*tau))/(4.68d0*tau+1.d-8) &
      + (0.5d0+dsign(0.5d0,tau-7.d0))/(4.d0*tau*dsqrt(dlog(dmax1(tau*0.56419d0,4.d0)))+1.d-8)
END SUBROUTINE fesc

SUBROUTINE LEVC2(T,b21,ndH,ndH2,nde,ndCp,n1,n2)
double precision :: T,b21,ndH,ndH2,nde,ndCp,n1,n2
double precision :: c21,c12,A21
c21 = 8.854d-8*nde/dsqrt(T) + 9.399d-10*(T**0.07d0)*(ndH+ndH2); c12 = 2.d0*c21*dexp(-0.092d0/T); A21 = 2.4d-6*b21
n2  = ndCp*c12/(c12+c21+A21); n1  = ndCp*(c21+A21)/(c12+c21+A21)
END SUBROUTINE LEVC2


SUBROUTINE IMC( P,n,dt,i,j,k )
USE comvar
double precision P,n,T,dt
double precision Pu,Pd,Pm,fev,iud
double precision CooL,Pmold,nkbi
Pu    = 1.d10
Pd    = 1.d-10
Pm    = 1.d1
Pmold = -1.d0
nkbi  = 1.d0/(kb*n)
do kkk = 1,30
  T = Pm*nkbi
  call Fcool(CooL,T,i,j,k)
  fev = Pm + CooL*dt*gammi1 - P
  iud = dsign(1.d0,fev)
  Pd  = dmax1(-iud*Pm,Pd)
  Pu  = dmin1(2.d0*Pm/(iud+1.d0+1.d-10),Pu)
  Pm  = 1.d1**(0.5d0*(dlog10(Pu)+dlog10(Pd)))
  if(dabs(Pm-Pmold).lt.1.d-5) goto 835
  Pmold = Pm
end do
835 continue
P = Pm
END SUBROUTINE IMC


SUBROUTINE SHIELD()
USE comvar
USE mpivar
USE chmvar
INCLUDE 'mpif.h'
double precision, dimension(:,:,:), allocatable :: Nttot,NtH2,NtC,NtCO,tau,temp
double precision :: tNtot,tNH2,tNC,tNCO,ttau,dxh
double precision :: T,fesC1,fesC2,n1,n2,b21,dvin
INTEGER MSTATUS(MPI_STATUS_SIZE)

ALLOCATE( Nttot(ndy,ndz,0:NSPLTx-1),NtH2(ndy,ndz,0:NSPLTx-1),NtC(ndy,ndz,0:NSPLTx-1),NtCO(ndy,ndz,0:NSPLTx-1) )
ALLOCATE(   tau(ndy,ndz,0:NSPLTx-1),temp(ndx,ndy,ndz) )

dxh = 0.5d0 * dx(1)

dvin = 1.d0/dv 
do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
  call fesc(tCII(i,j,k,1),fesC1); call fesc(tCII(i,j,k,2),fesC2); b21 = fesC1+fesC2
  T = U(i,j,k,5)/( kb*(ndp(i,j,k)+ndH(i,j,k)+ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)) )
  call LEVC2(T,b21,ndH(i,j,k),ndH2(i,j,k),nde(i,j,k),ndCp(i,j,k),n1,n2)
  temp(i,j,k) = 1.15563d1*(n1-n2)*dvin
end do; end do;end do

do k = 1, Ncellz; do j = 1, Ncelly
  tNtot=0.d0; tNH2=0.d0; tNC=0.d0; tNCO=0.d0; ttau=0.d0
  do i = 1, Ncellx
    tNtot = tNtot + dxh * ndtot(i,j,k)
    tNH2  = tNH2  + dxh *  ndH2(i,j,k)
    tNC   = tNC   + dxh *   ndC(i,j,k)
    tNCO  = tNCO  + dxh *  ndCO(i,j,k)
    ttau  = ttau  + dxh *  temp(i,j,k)
    Ntot(i,j,k,1) = tNtot
    NH2(i,j,k,1)  = tNH2
    NnC(i,j,k,1)  = tNC
    NCO(i,j,k,1)  = tNCO
    tCII(i,j,k,1) = ttau
    tNtot = tNtot + dxh * ndtot(i,j,k)
    tNH2  = tNH2  + dxh *  ndH2(i,j,k)
    tNC   = tNC   + dxh *   ndC(i,j,k)
    tNCO  = tNCO  + dxh *  ndCO(i,j,k)
    ttau  = ttau  + dxh *  temp(i,j,k)
  end do
  tNtot=0.d0; tNH2=0.d0; tNC=0.d0; tNCO=0.d0; ttau=0.d0
  do i = Ncellx, 1, -1
    tNtot = tNtot + dxh * ndtot(i,j,k)
    tNH2  = tNH2  + dxh *  ndH2(i,j,k)
    tNC   = tNC   + dxh *   ndC(i,j,k)
    tNCO  = tNCO  + dxh *  ndCO(i,j,k)
    ttau  = ttau  + dxh *  temp(i,j,k)
    Ntot(i,j,k,2) = tNtot
    NH2(i,j,k,2)  = tNH2
    NnC(i,j,k,2)  = tNC
    NCO(i,j,k,2)  = tNCO
    tCII(i,j,k,2) = ttau
    tNtot = tNtot + dxh * ndtot(i,j,k)
    tNH2  = tNH2  + dxh *  ndH2(i,j,k)
    tNC   = tNC   + dxh *   ndC(i,j,k)
    tNCO  = tNCO  + dxh *  ndCO(i,j,k)
    ttau  = ttau  + dxh *  temp(i,j,k)
  end do
  Nttot(j,k,IST)=tNtot; NtH2(j,k,IST)=tNH2; NtC(j,k,IST)=tNC; NtCO(j,k,IST)=tNCO; tau(j,k,IST)=ttau
end do; end do
      
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
DO L = 1, NSPLTx-1
  ITO    = NRANK + L
  ITOist = IST   + L
  IFR    = NRANK - L
  IFRist = IST   - L
  IF(ITOist.GE.NSPLTx) THEN; ITO = ITO - NSPLTx; END IF
  IF(IFRist.LT.0     ) THEN; IFR = IFR + NSPLTx; IFRist = IFRist + NSPLTx; END IF
  CALL MPI_SENDRECV(Nttot(1,1,IST   ),ndy*ndz,MPI_REAL8,ITO,1, &
                    Nttot(1,1,IFRist),ndy*ndz,MPI_REAL8,IFR,1,MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV( NtH2(1,1,IST   ),ndy*ndz,MPI_REAL8,ITO,1, &
                     NtH2(1,1,IFRist),ndy*ndz,MPI_REAL8,IFR,1,MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(  NtC(1,1,IST   ),ndy*ndz,MPI_REAL8,ITO,1, &
                      NtC(1,1,IFRist),ndy*ndz,MPI_REAL8,IFR,1,MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV( NtCO(1,1,IST   ),ndy*ndz,MPI_REAL8,ITO,1, &
                     NtCO(1,1,IFRist),ndy*ndz,MPI_REAL8,IFR,1,MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(  tau(1,1,IST   ),ndy*ndz,MPI_REAL8,ITO,1, &
                      tau(1,1,IFRist),ndy*ndz,MPI_REAL8,IFR,1,MPI_COMM_WORLD,MSTATUS,IERR)
END DO
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

DO L = 1, NSPLTx-1
  ITO   = IST + L
  IF(ITO.GE.NSPLTx) GOTO 177
  do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
    Ntot(i,j,k,2) = Ntot(i,j,k,2) + Nttot(j,k,ITO)
     NH2(i,j,k,2) =  NH2(i,j,k,2) +  NtH2(j,k,ITO)
     NnC(i,j,k,2) =  NnC(i,j,k,2) +   NtC(j,k,ITO)
     NCO(i,j,k,2) =  NCO(i,j,k,2) +  NtCO(j,k,ITO)
    tCII(i,j,k,2) = tCII(i,j,k,2) +   tau(j,k,ITO)
  end do; end do; end do
END DO
177 CONTINUE

DO L = 1, NSPLTx-1
  ITO   = IST - L
  IF(ITO.LT.0) GOTO 178
  do k = 1, Ncellz; do j = 1, Ncelly; do i = 1, Ncellx
    Ntot(i,j,k,1) = Ntot(i,j,k,1) + Nttot(j,k,ITO)
     NH2(i,j,k,1) =  NH2(i,j,k,1) +  NtH2(j,k,ITO)
     NnC(i,j,k,1) =  NnC(i,j,k,1) +   NtC(j,k,ITO)
     NCO(i,j,k,1) =  NCO(i,j,k,1) +  NtCO(j,k,ITO)
    tCII(i,j,k,1) = tCII(i,j,k,1) +   tau(j,k,ITO)
  end do; end do; end do
END DO
178 CONTINUE
DEALLOCATE( Nttot,NtH2,NtC,NtCO )
DEALLOCATE( tau,temp )

END SUBROUTINE SHIELD



SUBROUTINE RATES(i,j,k,T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep)
USE comvar
USE chmvar
DOUBLE PRECISION :: T,zeta,kHrec,kHerec,kH2,kH2ph,kH2dH,kH2de,kCO,kCOph,kCi,kCrec,kCOde,kCOdH,kHie,kHeie,kCie,kHiH,kHeiH,kCiH,kCOdHep,kH2dHep
DOUBLE PRECISION :: Av1,Av2,x1,x2,ATN2,ATN3,ATN4,ATN5,SHLD1,SHLD2,SHLC1,SHLC2

!( 1 pc * 5.3d-22 = 1.63542d-3 )
!( 1 pc * 2.d-15  = 6.1714d3 )
!( 1 pc * 1.d-17  = 3.0857d1 )
!( 1 pc * 1.405656457d-22  = 4.33743413d-4 )

Av1  = 1.63542d-3*Ntot(i,j,k,1); x1 = 6.1714d3*NH2(i,j,k,1)
Av2  = 1.63542d-3*Ntot(i,j,k,2); x2 = 6.1714d3*NH2(i,j,k,2)
ATN2 = ( dexp(-3.77358d0*Av1) +dexp(-3.77358d0*Av2) )*0.5d0
ATN3 = ( dexp(-2.3585d0*Av1)  +dexp(-2.3585d0*Av2) )*0.5d0

!( 1 pc * 5.e-22  = 1.54285d-3 )
!( 1 pc / 9.98337d16 = 3.09084d1 )
ATN4 = dmin1(1.d0, (3.0857d1*NH2(i,j,k,1))**(-3.d-2)*dexp(-1.54285d-3*NH2(i,j,k,1)) ) &
      +dmin1(1.d0, (3.0857d1*NH2(i,j,k,2))**(-3.d-2)*dexp(-1.54285d-3*NH2(i,j,k,2)) )
ATN4 = ATN4*0.5d0

!( 1 pc / 3.e15   = 1.02857d3 )
!( 1 pc / 8.2189020e14 = 3.75439d3 )
!( 1 pc / 4.2e16  = 7.3469d1 )
!( 4.35068d15 / 1 pc = 1.40995d-3 )
!( 8.97293d18 / 1 pc = 2.90791d0  )
!if( k==33 .and. j==33 .and. i==33) write(*,*) NCO(33,33,33,1),NCO(33,33,33,2),'11'

ATN5 = &
 (0.5d0-dsign(0.5d0,NCO(i,j,k,1)-1.40995d-3))*dexp(-1.d0*(1.02857d3*NCO(i,j,k,1))**0.6) &
+(0.5d0+dsign(0.5d0,NCO(i,j,k,1)-1.40995d-3))*dmin1( (3.75439d3*NCO(i,j,k,1)+1.d-100)**(-0.75d0),(7.3469d1*NCO(i,j,k,1)+1.d-100)**(-1.3d0) )
!if( k==33 .and. j==33 .and. i==33) write(*,*) ATN5,'12'
ATN5 = ATN5 + &
 (0.5d0-dsign(0.5d0,NCO(i,j,k,2)-1.40995d-3))*dexp(-1.d0*(1.02857d3*NCO(i,j,k,2))**0.6) &
+(0.5d0+dsign(0.5d0,NCO(i,j,k,2)-1.40995d-3))*dmin1( (3.75439d3*NCO(i,j,k,2)+1.d-100)**(-0.75d0),(7.3469d1*NCO(i,j,k,2)+1.d-100)**(-1.3d0) )
!if( k==33 .and. j==33 .and. i==33) write(*,*) ATN5,'13'
ATN5 = ATN5*0.5d0
!if( k==33 .and. j==33 .and. i==33) write(*,*) ATN5,'14'
!H Ionization
zeta  = 9.4671d-4
!H Recombination
kHrec = (0.5d0+dsign(0.5d0,15.78d0-T))*0.45d0*dlog(1.578d2/T) &
       +(0.5d0-dsign(0.5d0,15.78d0-T))*0.4d0*dsqrt(1.578d2/T)
kHrec = 2.05572d1*(T**(-0.5d0))*kHrec
!He Recombination
kHerec= 6.37623d1*(T**(-0.672d0))
!H2 formation
kH2   = 3.4569d-3*dsqrt(T)/(1.d0+1.26491d0*dsqrt(T+Tgr)+2.d0*T+8.d0*T**2)
!***H2 Photo-dissociation
SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))
kH2ph = 1.04138d3 * G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0
!H2 destruction by H collision
kH2dH = 1.07294d5*dexp(-4.39d1/T)
!H2 destruction by e collision
kH2de = 1.133d6*(T**0.35d0)*dexp(-1.02d2/T)
!***CO formation
kCO   = 1.57785d-2/(1.d0+G0*ATN3/(ndH2(i,j,k)*xo))
!***CO Photo-dissociation
kCOph = 3.155d3*G0*ATN3*ATN4*ATN5
!if( k==33 .and. j==33 .and. i==33) write(*,*) G0,ATN3,ATN4,ATN5,kCOph
!***C ionization
SHLC1 = dexp(-2.6d0*Av1-3.0857d1*NnC(i,j,k,1)-4.33743413d-4*NH2(i,j,k,1))/(1.d0+4.33743413d-4*NH2(i,j,k,1))
SHLC2 = dexp(-2.6d0*Av2-3.0857d1*NnC(i,j,k,2)-4.33743413d-4*NH2(i,j,k,2))/(1.d0+4.33743413d-4*NH2(i,j,k,2))
kCi   = 3.97d0*zeta + 6.62697d3*G0 * ( SHLC1+SHLC2 )*0.5d0
!C recombination
kCrec = 2.69267d1*(T**(-0.62d0))*(1.d0+1.d-4*dsqrt(nde(i,j,k))/(T**2))
!CO destruction by e collision
kCOde = 1.4d4*(T**(-0.5d0))*dexp(-1.14d2/T)
!CO destruction by H collision
kCOdH = 1.91858d3*(T**(-0.5d0))*dexp(-7.77d1/T)

!H  ionization by e collision
kHie  = 5.76d4*dsqrt(T)*dexp(-1.58d2/T)
!He ionization by e collision
kHeie = 2.76d4*(T**0.43d0)*dexp(-2.85d2/T)
!C  ionization by e collision
kCie  = 1.74d5*(T**0.4d0)*dexp(-1.31d2/T)

!H  ionization by H collision
kHiH  = 9.79d0*dsqrt(T)*dexp(-1.58d2/T)
!He ionization by H or p or H2 collision
kHeiH = 4.57d0*(T**0.43d0)*dexp(-2.85d2/T)
!C  ionization by H or p or H2 collision
kCiH  = 3.01d1*(T**0.4d0)*dexp(-1.31d2/T)
!H2 destruction by Hep recombination
kH2dHep = 1.1676d0*dexp(-0.035d0/T)
!CO destruction by Hep recombination
kCOdHep = 5.0491d4

END SUBROUTINE RATES


SUBROUTINE Omexp(omeps,eps)
DOUBLE PRECISION :: omeps,eps

omeps = ( 0.5d0+dsign(0.5d0,eps-1.d-4) )*( 1.d0-dexp(-eps) ) + &
( 0.5d0-dsign(0.5d0,eps-1.d-4) )*( eps-0.5d0*eps**2+0.16666666666666667d0*eps**3-4.16666666666666667d-2*eps**4 )

END SUBROUTINE Omexp


SUBROUTINE DISSIP()
USE comvar
USE mpivar

integer :: offset
double precision  :: divv(0:11,0:ndy,0:ndz)

N_MPI(20) = 3
N_MPI(1)  = 2
N_MPI(2)  = 3
N_MPI(3)  = 4
iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)

if((IST.eq.0).or.(IST.eq.NSPLTx-1)) then

  offset = 0; if(IST.eq.NSPLTx-1) offset = Ncellx-10
  do k=0,Ncellz+1; do j=0,Ncelly+1; do i=0,11
    ix = i+offset
    divv(i,j,k) = ( U(ix+1,j,k,2)-U(ix-1,j,k,2)+U(ix,j+1,k,3)-U(ix,j-1,k,3)+U(ix,j,k+1,4)-U(ix,j,k-1,4) )*0.25d0
  end do; end do; end do

  do k=1,Ncellz; do j=1,Ncelly; do i=1,10
    ix = i+offset
    U(ix,j,k,2) = U(ix,j,k,2) + 0.1d0*(divv(i+1,j,k)-divv(i-1,j,k))
    U(ix,j,k,3) = U(ix,j,k,3) + 0.1d0*(divv(i,j+1,k)-divv(i,j-1,k))
    U(ix,j,k,4) = U(ix,j,k,4) + 0.1d0*(divv(i,j,k+1)-divv(i,j,k-1))
  end do; end do; end do

end if

END SUBROUTINE DISSIP

