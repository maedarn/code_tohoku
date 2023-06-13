MODULE comvar
!INTEGER, parameter :: ndx=130, ndy=130, ndz=130, ndmax=130, Dim=3 !1024^3
INTEGER, parameter :: ndx=66, ndy=66, ndz=66, ndmax=66, Dim=3 !512^3
!INTEGER, parameter :: ndx=34, ndy=34, ndz=34, ndmax=34, Dim=3
!INTEGER, parameter :: ndx=18, ndy=18, ndz=18, ndmax=18, Dim=3
DOUBLE PRECISION, dimension(-1-1:ndx+1) :: x,dx
DOUBLE PRECISION, dimension(-1-1:ndy+1) :: y,dy
DOUBLE PRECISION, dimension(-1-1:ndz+1) :: z,dz
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: U, Bcc, Blg, Vfc, EMF
DOUBLE PRECISION, dimension(:,:,:),   allocatable :: dnc, xlag, dxlagM

DOUBLE PRECISION, parameter :: kb=8.63359d0, Kcond=1.6384d-2
DOUBLE PRECISION  :: gamma,gammi1,gammi2,gammi3,gampl1,gampl2,gampl3
DOUBLE PRECISION  :: CFL,facdep,tfinal,time,phr(-1:400)
DOUBLE PRECISION  :: pmin,pmax,rmin,rmax
double precision :: shusoku1=0.0d0 , phiratio=1.0d0!1.0d0/3.0d0
INTEGER :: Ncellx,Ncelly,Ncellz,iwx,iwy,iwz,maxstp,nitera
INTEGER :: ifchem,ifthrm,ifrad,ifgrv

!DOUBLE PRECISION :: cg=1.0d0,sourratio=0.5d0
DOUBLE PRECISION, parameter :: sourratio=0.5d0,adiff=0.25d0,rratio=0.20d0,rmove=0.d0!rmove=0.25d0!rmove=0.35d0!,adiff=0.375d0,cg=1.0d0,
double precision :: dx1,dy1,dz1,ddd,rrsph3,rrsph3x,rrsph3y,rrsph3z,tratio=0.5d0,cg=1.d0,Tdiff=0.2d0,rncn=0.d0,vmove=0.1d0!*cg/rmove
double precision :: kappa=1.d0/0.2d0,Msph1=0.d0,di_pos,tratio1
INTEGER, parameter :: mvstp=2!cfratio=5
INTEGER :: svci=50!cfratio=5
INTEGER ::ntdiv=5
!integer ifevogrv,ifevogrv2
!character(25) :: dir='/work/maedarn/3DMHD/test/' !samplecnv2
!character(43) :: dir='/work/maedarn/3DMHD/tel/tel-mesh128-cy-k20/'
!character(43) :: dir='/work/maedarn/3DMHD/tel/ts-paformance-time/'
character(43) :: dir='/data/group1/z40309n/telegraph/performance/'
character(17)  :: svdir
integer :: ntimeint=0,lsphmax=4,iwxts,iwyts,iwzts
DOUBLE PRECISION, dimension(:,:), allocatable :: time_pfm
END MODULE comvar

MODULE mpivar
INTEGER :: NPE,NRANK, NSPLTx,NSPLTy,NSPLTz,IST,JST,KST,LEFT,RIGT,BOTM,TOP,UP,DOWN
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
INTEGER :: point1(0:15),point2(0:15),NGL,NGcr,Nmem1,Nmem2,wvnum=1
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Phi ! , Phiexa
double precision, dimension(:,:,:), allocatable :: Phidt! , Phicgp , Phicgm
DOUBLE PRECISION :: Lbox,Lboxx,Lboxy,Lboxz
!double precision :: deltalength , cgcsratio= 1.0d0,cgratio1=0.2d0 !, shusoku1=0.0d0
double precision ::  cgcsratio= 1.0d0,cgratio1=0.2d0,rhomean !, shusoku1=0.0d0

!DOUBLE PRECISION , dimension(:,:,:,:), allocatable ::  Phicgm , Phi1step , Phi2step , Phicgp
DOUBLE PRECISION , dimension(:,:,:,:), allocatable ::  Phigrd
DOUBLE PRECISION , dimension(:,:,:), allocatable ::  Phiexa,Phiexab1,Phiexab2
DOUBLE PRECISION , dimension(:,:,:,:), allocatable ::  Phiwv, Phigrdwv!,Phiwvpre, Phigrdwvpre

INTEGER :: pointb1(0:15),pointb2(0:15)
DOUBLE PRECISION, dimension(:,:,:), allocatable :: bphil,bphir
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: bphigrdxl,bphigrdxr
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
!USE ISO_FORTRAN_ENV
!integer :: i_flow, i_flow_end=4000
integer :: i_flow_end=1
double precision :: dt
!character*5 :: NPENUM

CALL MPI_INIT(IERR)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPE  ,IERR)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,NRANK,IERR)
!NPE=NPE/40
!NRANK=NRANK/40
!write(*,*) NPE,NRANK
!write(*,*) 'OK'
!----- Prepare MPI SPLIT -----------------------------------------------!

!if(NPE.eq.4)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 1; end if
if(NPE.eq.4)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 1; end if
if(NPE.eq.8)    then; NSPLTx = 2; NSPLTy = 2; NSPLTz = 2; end if
if(NPE.eq.16)   then; NSPLTx = 4; NSPLTy = 2; NSPLTz = 2; end if
if(NPE.eq.32)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 2; end if
if(NPE.eq.40)   then; NSPLTx = 5; NSPLTy = 4; NSPLTz = 2; end if
if(NPE.eq.48)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 3; end if
if(NPE.eq.96)   then; NSPLTx = 6; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.64)   then; NSPLTx = 4; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.128)  then; NSPLTx = 8; NSPLTy = 4; NSPLTz = 4; end if
if(NPE.eq.192)  then; NSPLTx = 8; NSPLTy = 6; NSPLTz = 4; end if
if(NPE.eq.240)  then; NSPLTx = 8; NSPLTy = 6; NSPLTz = 5; end if
if(NPE.eq.256)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 4; end if
if(NPE.eq.384)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 6; end if
if(NPE.eq.480)  then; NSPLTx =10; NSPLTy = 8; NSPLTz = 6; end if
if(NPE.eq.512)  then; NSPLTx = 8; NSPLTy = 8; NSPLTz = 8; end if
if(NPE.eq.720)  then; NSPLTx =10; NSPLTy = 9; NSPLTz = 8; end if
if(NPE.eq.768)  then; NSPLTx =12; NSPLTy = 8; NSPLTz = 8; end if
if(NPE.eq.960)  then; NSPLTx =12; NSPLTy =10; NSPLTz = 8; end if
if(NPE.eq.1024) then; NSPLTx = 8; NSPLTy = 8; NSPLTz =16; end if
if(NPE.eq.1200) then; NSPLTx =12; NSPLTy =10; NSPLTz =10; end if
if(NPE.eq.1440) then; NSPLTx =12; NSPLTy =12; NSPLTz =10; end if
if(NPE.eq.1536) then; NSPLTx =16; NSPLTy =12; NSPLTz = 8; end if
if(NPE.eq.1920) then; NSPLTx =16; NSPLTy =12; NSPLTz =10; end if
if(NPE.eq.2400) then; NSPLTx =20; NSPLTy =12; NSPLTz =10; end if
if(NPE.eq.4800) then; NSPLTx =20; NSPLTy =20; NSPLTz =12; end if
if(NPE.eq.4800) then; NSPLTx =20; NSPLTy =20; NSPLTz =12; end if
if(NPE.eq.9600) then; NSPLTx =24; NSPLTy =20; NSPLTz =20; end if
if(NPE.eq.24000)then; NSPLTx =40; NSPLTy =30; NSPLTz =20; end if
if(NPE.eq.24000)then; NSPLTx =40; NSPLTy =30; NSPLTz =20; end if
if(NPE.eq.36480)then; NSPLTx =48; NSPLTy =40; NSPLTz =19; end if
if(NPE.eq.48000)then; NSPLTx =40; NSPLTy =40; NSPLTz =30; end if



!loopbc=NSPLTy*NSPLTz/(ndx-2)-1

!   write(*,*) 'OK1',NRANK
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
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'OK2'
ALLOCATE( U(-1:ndx,-1:ndy,-1:ndz,8) )
ALLOCATE(ndH(-1:ndx,-1:ndy,-1:ndz),ndp(-1:ndx,-1:ndy,-1:ndz),ndH2(-1:ndx,-1:ndy,-1:ndz),ndHe(-1:ndx,-1:ndy,-1:ndz), &
       ndHep(-1:ndx,-1:ndy,-1:ndz),ndC(-1:ndx,-1:ndy,-1:ndz),ndCp(-1:ndx,-1:ndy,-1:ndz),ndCO(-1:ndx,-1:ndy,-1:ndz), &
         nde(-1:ndx,-1:ndy,-1:ndz),ndtot(-1:ndx,-1:ndy,-1:ndz),Ntot(-1:ndx,-1:ndy,-1:ndz,2),                        &
         NH2(-1:ndx,-1:ndy,-1:ndz,2),NnC(-1:ndx,-1:ndy,-1:ndz,2),NCO(-1:ndx,-1:ndy,-1:ndz,2),tCII(-1:ndx,-1:ndy,-1:ndz,2) )
ALLOCATE(DTF(-1:(ndx-2)*NSPLTx+2,-1:ndy,-1:ndz))
ALLOCATE(Phi(-1:ndx,-1:ndy,-1:ndz))
ALLOCATE(time_pfm(0:NPE-1,1:12))

!*********grvwave*********
ALLOCATE(Phiexa(-1-1:ndx+1,-1-1:ndy+1,-1-1:ndz+1))
!ALLOCATE(Phiexab1(-1-1:ndx+1,-1-1:ndy+1,-1-1:ndz+1),Phiexab2(-1-1:ndx+1,-1-1:ndy+1,-1-1:ndz+1))
ALLOCATE(Phigrd(-1:ndx,-1:ndy,-1:ndz,1:wvnum))

ALLOCATE(Phiwv(-1:ndx,-1:ndy,-1:ndz,1:wvnum))!,Phiwvpre(-1:ndx,-1:ndy,-1:ndz,1:wvnum))
ALLOCATE(Phigrdwv(-1:ndx,-1:ndy,-1:ndz,1:wvnum))!,Phigrdwvpre(-1:ndx,-1:ndy,-1:ndz,1:wvnum))
ALLOCATE(bphil(-3:ndy+2,-3:ndz+2,-1:1     ))
ALLOCATE(bphir(-3:ndy+2,-3:ndz+2,ndx-2:ndx))
ALLOCATE(bphigrdxl(-1:ndy,-1:ndz,-1:1     ,1:wvnum))
ALLOCATE(bphigrdxr(-1:ndy,-1:ndz,ndx-2:ndx,1:wvnum))
!*********grvwave*********

!write(*,*) 'OK3'

call INITIA
call SELFGRAVWAVE(0.0d0,0)
!call SELFGRAVWAVE(0.0d0,4)
!call SELFGRAVWAVE(0.0d0,4)
!write(*,*)'dt',dx1/cg*tratio
!call fapp_start("loop1",1,0)
!time_pfm(:,:)=0.d0
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!time_pfm(NRANK,1)=MPI_WTIME()
!write(*,*)
!call fapp_start("loop1",1,0)
!call fipp_start
!do i_flow=1,i_flow_end
!do k=-1,ndz; do j=-1,ndy; do i=-1,ndx
!Phiwv(i,j,k,1)=1.d0
!Phigrdwv(i,j,k,1)=1.d0*dt
!enddo; enddo; enddo
dt=dx1/cg*tratio
call slvmuscle(dt)
!enddo

!call SELFGRAVWAVE(0.0d0,4)
!time_pfm(NRANK,2)=MPI_WTIME()
!time_pfm(NRANK,3)=(time_pfm(NRANK,2)-time_pfm(NRANK,1))/dble(i_flow_end)

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!do Nroot=0,NPE-1
!CALL MPI_BCAST(time_pfm(Nroot,1),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
!CALL MPI_BCAST(time_pfm(Nroot,2),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
!CALL MPI_BCAST(time_pfm(Nroot,3),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
!end do

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!WRITE(NPENUM,'(I5.5)') NRANK
 !if(NRANK==0) then
!open(521,file=dir//svdir//'/CPU_TIME'//NPENUM//'.DAT',access='stream',FORM='UNFORMATTED')
 !do Nroot=0,NPE-1
!write(521) time_pfm(NRANK,1),time_pfm(NRANK,2),time_pfm(NRANK,3)
 !end do
!close(521)
!endif

!call fipp_stop
!call fapp_stop("loop1",1,0)

!call EVOLVE

!write(*,*) 'OK'

DEALLOCATE(U)
DEALLOCATE(ndH,ndp,ndH2,ndHe,ndHep,ndC,ndCp,ndCO,nde,ndtot,Ntot,NH2,NnC,NCO,tCII)
DEALLOCATE(DTF)
DEALLOCATE(Phi)

!********gravwave**********
DEALLOCATE(Phiexa,Phigrd)
!DEALLOCATE(Phiexab1,phiexab2)
DEALLOCATE(Phiwv,Phigrdwv)
!DEALLOCATE(Phiwvpre,Phigrdwvpre)
DEALLOCATE(bphil,bphir,bphigrdxl,bphigrdxr)
DEALLOCATE(time_pfm)
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

integer :: Np1x, Np2x, Np1y, Np2y, Np1z, Np2z, nunit, ix, jy, kz
double precision ::  ql1x,ql2x,ql1y,ql2y,ql1z,ql2z,dinit1,dinit2,pinit1,pinit2, &
           vinitx1,vinitx2,vinity1,vinity2,vinitz1,vinitz2,           &
           binitx1,binitx2,binity1,binity2,binitz1,binitz2
double precision, dimension(:), allocatable :: x_i,y_i,z_i,dx_i,dy_i,dz_i
double precision :: pi!,amp!,xpi,ypi,zpi,phase1,phase2,phase3
double precision :: Hini,pini,H2ini,Heini,Hepini,Cini,COini,Cpini,dBC
!double precision :: ampn(2048),ampn0(2048)
character*3 :: NPENUM!,MPIname
!INTEGER :: MSTATUS(MPI_STATUS_SIZE)
!double precision, dimension(:,:), allocatable :: plane,rand
!integer i3,i4,i2y,i2z,rsph2,pls
!double precision cenx,ceny,cenz,rsph,rrsph,Hsheet,censh,minexa,rsph3

open(8,file=dir//'INPUT3D.DAT')
  read(8,*)  svdir
  read(8,*)  cg,Tdiff,rncn,ntdiv,svci,vmove
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
  read(8,*)  CFL,facdep,di_pos
  read(8,*)  maxstp,nitera,tfinal
  read(8,*)  BCx1,BCx2,BCy1,BCy2,BCz1,BCz2
  read(8,*)  ifchem,ifthrm,ifrad,ifgrv
  read(8,*)  iwxts,iwyts,iwzts,tratio1
close(8)
tratio=tratio*tratio1
vmove=vmove*cg/rmove
kappa=0.5d0/Tdiff


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
!write(*,*) NRANK,'INIT'
ALLOCATE(dx_i(-1-1:Ncellx*NSPLTx+2+1)); ALLOCATE(dy_i(-1-1:Ncelly*NSPLTy+2+1)); ALLOCATE(dz_i(-1-1:Ncellz*NSPLTz+2+1))
ALLOCATE( x_i(-1-1:Ncellx*NSPLTx+2+1)); ALLOCATE( y_i(-1-1:Ncelly*NSPLTy+2+1)); ALLOCATE( z_i(-1-1:Ncellz*NSPLTz+2+1))

do i = -1-1, Np1x
  dx_i(i) = ql1x/dble(Np1x)
end do
do i = Np1x+1, Ncellx*NSPLTx+2+1
  dx_i(i) = ql2x/dble(Np2x)
end do
do j = -1-1, Ncelly*NSPLTy+2+1
  dy_i(j) = ql1y/dble(Np1y)
end do
do k = -1-1, Ncellz*NSPLTz+2+1
  dz_i(k) = ql1z/dble(Np1z)
end do

!dx=dy=dz
dx1= dx_i(0)
dy1= dy_i(0)
dz1= dz_i(0)
!dx=dy=dz
dt=dx1/cg*tratio

x_i(-1) = -dx_i(0)
x_i(-2) = x_i(-1)-dx_i(0)
do i = 0, Ncellx*NSPLTx+2+1
   x_i(i) = x_i(i-1) + dx_i(i)
   !write(*,*) 'x', x_i(i)
end do
y_i(-1) = -dy_i(0)
y_i(-2) = y_i(-1)-dy_i(0)
do j = 0, Ncelly*NSPLTy+2+1
   y_i(j) = y_i(j-1) + dy_i(j)
   !write(*,*) 'y', y_i(j)
end do
z_i(-1) = -dz_i(0)
z_i(-2) = z_i(-1)-dz_i(0)
do k = 0, Ncellz*NSPLTz+2+1
   z_i(k) = z_i(k-1) + dz_i(k)
   !write(*,*) 'z', z_i(k)
end do

do i = -1-1, Ncellx+2+1
  ix    =  IST*Ncellx + i
  x(i)  =  x_i(ix)
  dx(i) =  dx_i(ix)
end do
do j = -1-1, Ncelly+2+1
  jy    =  JST*Ncelly + j
  y(j)  =  y_i(jy)
  dy(j) =  dy_i(jy)
end do
do k = -1-1, Ncellz+2+1
  kz    =  KST*Ncellz + k
  z(k)  =  z_i(kz)
  dz(k) =  dz_i(kz)
end do



dt=dx(1)/cg*tratio
DEALLOCATE(dx_i); DEALLOCATE(dy_i); DEALLOCATE(dz_i); DEALLOCATE(x_i); DEALLOCATE(y_i); DEALLOCATE(z_i)

nunit=1
!***** Read Initial Conditions *****!
if(nunit.eq.1) goto 120
  WRITE(NPENUM,'(I3.3)') NRANK
  open(unit=8,file=dir//'000'//NPENUM//'.dat',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
  Lboxx=ql1x+ql2x
  Lboxy=ql1y+ql2y
  Lboxz=ql1z+ql2z
  !call SELFGRAVWAVE(0.0d0,1)
  call SELFGRAVWAVE(0.0d0,0) !密度場の生成の時
  !call SELFGRAVWAVE(0.0d0,6) !calculate cg
end if
END SUBROUTINE INITIA

RECURSIVE subroutine SELFGRAVWAVE(dt,mode)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'
  integer :: mode,count=1,ndt1=0,svc1=0!,svci=50,rdnum
  DOUBLE PRECISION  :: dt, eps=1.0d-3
  character(3) NPENUM
  character(6) countcha
  integer nt
  double precision :: tdm,pi=3.14159265358979323846d0,amp,wvpt

  !**************** INITIALIZEATION **************
  if(mode==0) then
     Phi(:,:,:) = 0.0d0
     Phiwv(:,:,:,:)=0.d0
     Phigrdwv(:,:,:,:)=0.d0

     do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
         Phiwv(i,j,k,1)=0.d0!Phiexa(i,j,k)
         !Phigrdwv(i,j,k,1)=cg*Phigrd(i,j,k,1)*2.d0/2.d0/kappa+Phiwv(i,j,k,1)
         Phigrdwv(i,j,k,1)=0.d0!cg*Phigrd(i,j,k,1)+kappa*Phiexa(i,j,k)
     end do
     end do
     end do

     !test of wave propagation
     do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
     !xpi = 0.5d0*( x(i)+x(i-1) )
     amp = 1.d5
     wvpt = 4.d0
     U(i,j,k,1) = amp*dcos(dble(iwxts)*wvpt*2.d0*pi*x(i)/Lboxx)*dcos(dble(iwyts)*wvpt*2.d0*pi*y(j)/Lboxy)*dcos(dble(iwzts)*wvpt*2.d0*pi*z(k)/Lboxz)
     Phiexa(i,j,k) = -amp*G4pi/wvpt/wvpt/((2.d0*pi/Lboxx)**2.d0+(2.d0*pi/Lboxy)**2.d0+(2.d0*pi/Lboxz)**2.d0)*dcos(dble(iwxts)*wvpt*2.d0*pi*x(i)/Lboxx)*dcos(dble(iwyts)*wvpt*2.d0*pi*y(j)/Lboxy)*dcos(dble(iwzts)*wvpt*2.d0*pi*z(k)/Lboxz)
     Phiwv(i,j,k,1)   = 0.d0
!-amp*G4pi/wvpt/wvpt/((2.d0*pi/Lboxx)**2.d0+(2.d0*pi/Lboxy)**2.d0+(2.d0*pi/Lboxz)**2.d0)*dcos(dble(iwxts)*wvpt*2.d0*pi*x(i)/Lboxx)*dcos(dble(iwyts)*wvpt*2.d0*pi*y(j)/Lboxy)*dcos(dble(iwzts)*wvpt*2.d0*pi*z(k)/Lboxz)
     Phigrdwv(i,j,k,1)= 0.d0
!-amp*G4pi/wvpt/wvpt/((2.d0*pi/Lboxx)**2.d0+(2.d0*pi/Lboxy)**2.d0+(2.d0*pi/Lboxz)**2.d0)*dcos(dble(iwxts)*wvpt*2.d0*pi*x(i)/Lboxx)*dcos(dble(iwyts)*wvpt*2.d0*pi*y(j)/Lboxy)*dcos(dble(iwzts)*wvpt*2.d0*pi*z(k)/Lboxz)
     
     end do; end do; end do
     !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
     !write(*,*)'initial',Lbox
  end if


  !****************GRAVITY SOLVER*****************

  if(mode==2) then
     N_MPI(20)=1; N_MPI(1)=1
     iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)
     !---debug---
     !call  SELFGRAVWAVE(0.0d0,4)
     !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
     !write(*,*) '------pb1-------' ,Nrank

     !****calcurate bc****
     !call collect()
     !Call PB( 0)
     !Call PB(-1)
     !Call PB(-2)
     !call pbphigrd(dt)
     !****calcurate bc****
     !ndt1=ndt1+1
     !if(ndt1==maxstp*nitera/mvstp+1) then
     !call  SELFGRAVWAVE(0.0d0,4)
 !    call move()
 !    Phiwv(:,:,:,:)=0.d0
 !    Phigrdwv(:,:,:,:)=0.d0
     !ndt1=0
     !endif
 !    call movesph(dt)
     tdm=dt!/dble(ntdiv)
     do nt=1,ntdiv
    ! iwx=1;iwy=1;iwz=1
    ! call BCgrv(100,1,8)
    ! if(mod(svc1,svci)==0) then
    ! call SELFGRAVWAVE(0.0d0,4)
     !call movesph(dt)
    ! endif
     call slvmuscle(tdm)
    ! svc1=svc1+1
     enddo
  end if



  !***************SAVE PHI & PHIDT FOR DEBUG & INITIAL**************
  if(mode==4) then
     !write(*,*) 'save???'
     WRITE(NPENUM,'(I3.3)') NRANK
     WRITE(countcha,'(I6.6)') count
     write(*,*)'SAVE_Phi_pre',count,dir,svdir
     open(17,file=dir//svdir//'/PHI'//countcha//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     write(*,*)'SAVE_Phi_pre_op',count,dir,svdir,Ncellz,Ncelly,Ncellx
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI2step'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !write(*,*) 'save?????'

 
     !-------------------TEST---------------------
     !iwx = 1; iwy = 1; iwz = 1
     !call BCgrv(101)
     !call BCgrv(102)
     do k = 1, Ncellz
        !write(*,*) 'write',NRANK,k,sngl(Phiwv(1,1,k,1)),sngl(Phigrdwv(1,1,k,1)),sngl(Phiexa(1,1,k)),sngl(cg*Phigrd(1,1,k,1)+kappa*Phiexa(1,1,k)),sngl(U(1,1,k,1))
        do j = 1, Ncelly
           do i = 1, Ncellx
           !write(28) sngl(Phiwv(i,j,k,1)),sngl(Phigrdwv(i,j,k,1)),sngl(Phiexa(i,j,k)),sngl(Phigrd(i,j,k,1)),sngl(U(i,j,k,1))
           !write(*,*) 'write',NRANK,i,j,k,sngl(Phiwv(i,j,k,1)),sngl(Phigrdwv(i,j,k,1)),sngl(Phiexa(i,j,k)),sngl(cg*Phigrd(i,j,k,1)+kappa*Phiexa(i,j,k)),sngl(U(i,j,k,1))
           write(17) sngl(Phiwv(i,j,k,1)),sngl(Phigrdwv(i,j,k,1)),sngl(Phiexa(i,j,k)),sngl(cg*Phigrd(i,j,k,1)+kappa*Phiexa(i,j,k)),sngl(U(i,j,k,1))
           !write(28) Phiwv(i,j,k,1),Phigrdwv(i,j,k,1)
          enddo
        end do
        !write(*,*) sngl(Phiwv(8,8,k,1)),sngl(Phigrdwv(8,8,k,1))
     end do
     close(17)
     count=count+1
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !write(*,*)'SAVE_Phi',count,dir,svdir
end if
end subroutine SELFGRAVWAVE



subroutine slvmuscle(dt)
use comvar
use slfgrv
use mpivar
INCLUDE 'mpif.h'
double precision :: dt,dtratio=dsqrt(3.0d0),coeffx=0.d0,coeffy=0.d0,coeffz=0.d0!,rhomean
integer :: i=0!,n,m,l!,countn
double precision :: adiff2=0.5d0,dtration=0.5d0
double precision :: nu2,w=6.0d0,eps=1.0d-10! , deltap,deltam,deltalen !kappa -> comver  better?
integer :: cnt=0
DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phiu,Phiugrd!,Phi2dt,Phi2dtgrd!,Phigrad,Phipre,Phipregrd,Phi2dt,Phi2dtgrd
DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phiy,Phiygrd!,Phiprez,Phipregrdz!,Phiprey_swp,Phipregrdy_swp
!DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phivec,Phivecgrd
!character(5) name
integer :: Lnum,Mnum,is,ie,n_exp=13,N_ol=2,idm!,idm,hazi,Ncell,Ncm,Ncl
!DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,-1:ndz) :: ul,ur
DOUBLE PRECISION , dimension(-1:ndx) :: slop,slopgrd
!double precision :: rho(-1:ndx,-1:ndy,-1:ndz)
double precision  grdxy1,grdyz1,grdzx1
!double precision  grdxy1zp,grdxy1zm,grdxy1mn,grdyz1xp,grdyz1xm,grdyz1mn,grdzx1yp,grdzx1ym,grdzx1mn
!double precision :: Phiwvpre(-1:ndx,-1:ndy,-1:ndz,1:1)!,Phigrdwvpre(-1:ndx,-1:ndy,-1:ndz,1:1)
double precision :: expand_exp,expand_dx!,expand_trm,expand_dbi
double precision :: kp_i,exp_m,exp_p,exp_k
integer :: iswp1,iswp2,i_flow, i_flow_end=4000
double precision :: delp,delm,delpgrd,delmgrd,phiwv_d!,phigrdwv_d
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt
!INTEGER :: blki=4*1024/8,ii,blkii=96,blkjj=16,jj


!do i_flow=1,i_flow_end
!call fapp_start("loop1",1,0)
do i_flow=1,i_flow_end

call fapp_start("loop1",1,0)

N_ol=1
idm=1

CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
LEFTt = LEFT!; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
RIGTt = RIGT!; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
Phiwv(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(1            ,-1,-1,idm),1,VECU,LEFT,1, &
Phiwv(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt

CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
      Phiwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
      Phiwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
TOP = TOPt; BOTM = BOTMt

CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
Phiwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
 Phiwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
UP = UPt; DOWN = DOWNt

call fapp_stop("loop1",1,0)


call fapp_start("loop2",1,0)
expand_dx=-2.d0*kappa * 0.5d0 * dt
expand_exp=dexp(expand_dx)
kp_i=1.0/(2.d0*kappa+1.d-10)
exp_m=(1.d0-expand_exp)
exp_p=(1.d0+expand_exp)
exp_k=exp_m*kp_i


!call fapp_start("loop1",1,0)
!do k=-1,ndz; do j=-1,ndy; do i=-1,ndx
do k=0,ndz-1; do j=0,ndy-1
!do ii=1,ndx-1,blki
!do i=ii,min(ii+blki-1,ndx-1)
do i=0,ndx-1
phiwv_d=Phiwv(i,j,k,1)
!enddo
!do i=ii,min(ii+blki-1,ndx-1)
!phiwv_d=Phiwv(i,j,k,1)
Phiwv(i,j,k,1)    = 0.5d0*Phiwv(i,j,k,1)*exp_p+Phigrdwv(i,j,k,1)*exp_k
Phigrdwv(i,j,k,1) = 0.5d0*Phigrdwv(i,j,k,1)*exp_p+0.5d0*kappa*phiwv_d*exp_m
!enddo
enddo
enddo; enddo

call fapp_stop("loop2",1,0)



call fapp_start("loop3",1,0)
do k=1,ndz-2; do j=1,ndy-2; do i=1,ndx-2
grdxy1=adiff*Phiwv(i+1,j+1,k,1)+adiff*Phiwv(i-1,j-1,k,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,1)
grdyz1=adiff*Phiwv(i,j+1,k+1,1)+adiff*Phiwv(i,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,1)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,1)
grdzx1=adiff*Phiwv(i+1,j,k+1,1)+adiff*Phiwv(i-1,j,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,1)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,1)

Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
(-2.d0*cg*cg*grdxy1/dx1/dy1 &
 -2.d0*cg*cg*grdyz1/dy1/dz1 &
 -2.d0*cg*cg*grdzx1/dz1/dx1) *dt * dtration &
-G4pi*cg*cg*U(i,j,k,1)*dt * dtration
enddo; enddo; enddo
call fapp_stop("loop3",1,0)


call fapp_start("loop4",1,0)
N_ol=2
idm=1

CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
LEFTt = LEFT!; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
RIGTt = RIGT!; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
Phiwv(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(1            ,-1,-1,idm),1,VECU,LEFT,1, &
Phiwv(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phigrdwv(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
Phigrdwv(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phigrdwv(1            ,-1,-1,idm),1,VECU,LEFT,1, &
Phigrdwv(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt

CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
      Phiwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
      Phiwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phigrdwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
     Phigrdwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phigrdwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
     Phigrdwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
TOP = TOPt; BOTM = BOTMt

CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
Phiwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
Phiwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phigrdwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
Phigrdwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phigrdwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
Phigrdwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
UP = UPt; DOWN = DOWNt

call fapp_stop("loop4",1,0)

call fapp_start("loop5",1,0)

is = 1
ie = ndx-2
nu2 = cg * dt / dx1


DO Lnum = -1, ndz;DO Mnum = -1, ndy
do i = is-1 , ie+1
delp = Phiwv(i+1,Mnum,Lnum,1)-Phiwv(i  ,Mnum,Lnum,1)
delm = Phiwv(i  ,Mnum,Lnum,1)-Phiwv(i-1,Mnum,Lnum,1)
slop(i) = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )

delpgrd = Phigrdwv(i+1,Mnum,Lnum,1)-Phigrdwv(i,Mnum,Lnum,1)
delmgrd = Phigrdwv(i,Mnum,Lnum,1)-Phigrdwv(i-1,Mnum,Lnum,1)
slopgrd(i) = dmax1( 0.d0,(2.d0*delpgrd*delmgrd+eps)/(delpgrd**2+delmgrd**2+eps) )
end do
do i = is-1,ie+1
Phiu(i,Mnum,Lnum) = Phiwv(i,Mnum,Lnum,1)- 0.5d0 * nu2 * ( Phiwv(i,Mnum,Lnum,1) - Phiwv(i-1,Mnum,Lnum,1)) &
     + 0.25d0 * 1.d0 * slop(i) * ((1.0d0-slop(i)*1.d0/3.0d0)*(Phiwv(i,Mnum,Lnum,1)-Phiwv(i-1,Mnum,Lnum,1)) + &
     (1.0d0+slop(i)*1.d0/3.0d0)*(Phiwv(i+1,Mnum,Lnum,1) - Phiwv(i,Mnum,Lnum,1)))

Phiugrd(i,Mnum,Lnum) = Phigrdwv(i,Mnum,Lnum,1) + 0.5d0 * nu2 * ( Phigrdwv(i+1,Mnum,Lnum,1) - Phigrdwv(i,Mnum,Lnum,1)) &
     - 0.25d0 * 1.d0 * slopgrd(i) * ((1.0d0+slopgrd(i)*1.d0/3.0d0)*(Phigrdwv(i,Mnum,Lnum,1)-Phigrdwv(i-1,Mnum,Lnum,1)) + &
     (1.0d0-slopgrd(i)*1.d0/3.0d0)*(Phigrdwv(i+1,Mnum,Lnum,1) - Phigrdwv(i,Mnum,Lnum,1)))
end do
end DO;end DO

!do jj=-1,ndz,blkjj
!do ii=-1,ndy,blkii
!do Lnum=jj,min(jj+blkjj-1,ndz)
!do Mnum=ii,min(ii+blkii-1,ndy)
!do i = is,ie
!iswp1=Mnum
!iswp2=Lnum
!Phiy(i,iswp2,iswp1) = Phiwv(i,Mnum,Lnum,1) - nu2 * (Phiu(i,Mnum,Lnum) - Phiu(i-1,Mnum,Lnum))
!Phiygrd(i,iswp2,iswp1) = Phigrdwv(i,Mnum,Lnum,1) + nu2 * (Phiugrd(i+1,Mnum,Lnum) - Phiugrd(i,Mnum,Lnum))
!enddo
!enddo
!enddo
!enddo
!enddo



DO Lnum = -1, ndz;DO Mnum = -1, ndy;do i = is,ie
iswp1=Mnum
iswp2=Lnum
!Phiwv(i,Mnum,Lnum,1) = Phiwv(i,Mnum,Lnum,1) - nu2 * (Phiu(i,Mnum,Lnum) - Phiu(i-1,Mnum,Lnum))
Phiy(i,iswp2,iswp1) = Phiwv(i,Mnum,Lnum,1) - nu2 * (Phiu(i,Mnum,Lnum) - Phiu(i-1,Mnum,Lnum))
!Phigrdwv(i,Mnum,Lnum,1) = Phigrdwv(i,Mnum,Lnum,1) + nu2 * (Phiugrd(i+1,Mnum,Lnum) - Phiugrd(i,Mnum,Lnum))
Phiygrd(i,iswp2,iswp1) = Phigrdwv(i,Mnum,Lnum,1) + nu2 * (Phiugrd(i+1,Mnum,Lnum) - Phiugrd(i,Mnum,Lnum))
end do;end DO;end DO

call fapp_stop("loop5",1,0)
call fapp_start("loop6",1,0)

is = 1
ie = ndy-2
nu2 = cg * dt / dy1

DO Mnum = -1, ndz;DO Lnum = 1, ndx-2
do i = is-1 , ie+1
!delp = Phiwv(Lnum,i+1,Mnum,1)-Phiwv(Lnum,i,Mnum,1)
!delm = Phiwv(Lnum,i,Mnum,1)-Phiwv(Lnum,i-1,Mnum,1)
delp = Phiy(Lnum,Mnum,i+1)-Phiy(Lnum,Mnum,i)
delm = Phiy(Lnum,Mnum,i)-Phiy(Lnum,Mnum,i-1)
slop(i) = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )

delpgrd = Phiygrd(Lnum,Mnum,i+1)-Phiygrd(Lnum,Mnum,i)
delmgrd = Phiygrd(Lnum,Mnum,i)-Phiygrd(Lnum,Mnum,i-1)
slopgrd(i) = dmax1( 0.d0,(2.d0*delpgrd*delmgrd+eps)/(delpgrd**2+delmgrd**2+eps) )
end do
do i = is-1,ie+1
!Phiu(Lnum,i,Mnum) = Phiwv(Lnum,i,Mnum,1)- 0.5d0 * nu2 * ( Phiwv(Lnum,i,Mnum,1) - Phiwv(Lnum,i-1,Mnum,1)) &
!     + 0.25d0 * 1.d0 * slop(i) * ((1.0d0-slop(i)*1.d0/3.0d0)*(Phiwv(Lnum,i,Mnum,1)-Phiwv(Lnum,i-1,Mnum,1)) + &
!     (1.0d0+slop(i)*1.d0/3.0d0)*(Phiwv(Lnum,i+1,Mnum,1) - Phiwv(Lnum,i,Mnum,1)))
Phiu(Lnum,Mnum,i) = Phiy(Lnum,Mnum,i)- 0.5d0 * nu2 * ( Phiy(Lnum,Mnum,i) - Phiy(Lnum,Mnum,i-1)) &
     + 0.25d0 * 1.d0 * slop(i) * ((1.0d0-slop(i)*1.d0/3.0d0)*(Phiy(Lnum,Mnum,i)-Phiy(Lnum,Mnum,i-1)) + &
     (1.0d0+slop(i)*1.d0/3.0d0)*(Phiy(Lnum,Mnum,i+1) - Phiy(Lnum,Mnum,i)))

!Phiugrd(Lnum,i,Mnum) = Phigrdwv(Lnum,i,Mnum,1) + 0.5d0 * nu2 * ( Phigrdwv(Lnum,i+1,Mnum,1) - Phigrdwv(Lnum,i,Mnum,1)) &
!     - 0.25d0 * 1.d0 * slopgrd(i) * ((1.0d0+slopgrd(i)*1.d0/3.0d0)*(Phigrdwv(Lnum,i,Mnum,1)-Phigrdwv(Lnum,i-1,Mnum,1)) + &
!     (1.0d0-slopgrd(i)*1.d0/3.0d0)*(Phigrdwv(Lnum,i+1,Mnum,1) - Phigrdwv(Lnum,i,Mnum,1)))
Phiugrd(Lnum,Mnum,i) = Phiygrd(Lnum,Mnum,i) + 0.5d0 * nu2 * ( Phiygrd(Lnum,Mnum,i+1) - Phiygrd(Lnum,Mnum,i)) &
     - 0.25d0 * 1.d0 * slopgrd(i) * ((1.0d0+slopgrd(i)*1.d0/3.0d0)*(Phiygrd(Lnum,Mnum,i)-Phiygrd(Lnum,Mnum,i-1)) + &
     (1.0d0-slopgrd(i)*1.d0/3.0d0)*(Phiygrd(Lnum,Mnum,i+1) - Phiygrd(Lnum,Mnum,i)))
end do
end DO;end DO

!do i = is,ie
!do ii=-1,ndz,blkii
!do jj=-1,ndx-2,blkjj
!do Mnum=ii,min(ii+blkii-1,ndy)
!do Lnum=jj,min(jj+blkjj-1,ndz)
!iswp1=i
!iswp2=Mnum
!Phiwv(Lnum,iswp1,iswp2,1) = Phiy(Lnum,Mnum,i) - nu2 * (Phiu(Lnum,Mnum,i) - Phiu(Lnum,Mnum,i-1))
!Phigrdwv(Lnum,iswp1,iswp2,1) = Phiygrd(Lnum,Mnum,i) + nu2 * (Phiugrd(Lnum,Mnum,i+1) - Phiugrd(Lnum,Mnum,i))
!enddo
!enddo
!enddo
!enddo
!enddo

do i = is,ie;DO Mnum = -1, ndz;DO Lnum = 1, ndx-2
iswp1=i
iswp2=Mnum
!Phiwv(Lnum,i,Mnum,1) = Phiwv(Lnum,i,Mnum,1) - nu2 * (Phiu(Lnum,i,Mnum) - Phiu(Lnum,i-1,Mnum))
!Phigrdwv(Lnum,i,Mnum,1) = Phigrdwv(Lnum,i,Mnum,1) + nu2 * (Phiugrd(Lnum,i+1,Mnum) - Phiugrd(Lnum,i,Mnum))
Phiwv(Lnum,iswp1,iswp2,1) = Phiy(Lnum,Mnum,i) - nu2 * (Phiu(Lnum,Mnum,i) - Phiu(Lnum,Mnum,i-1))
Phigrdwv(Lnum,iswp1,iswp2,1) = Phiygrd(Lnum,Mnum,i) + nu2 * (Phiugrd(Lnum,Mnum,i+1) - Phiugrd(Lnum,Mnum,i))
end do;end DO;end DO

call fapp_stop("loop6",1,0)
call fapp_start("loop7",1,0)


is = 1
ie = ndz-2
nu2 = cg * dt / dz1
DO Lnum = 1, ndy-2;DO Mnum = 1, ndx-2
do i = is-1 , ie+1
delp = Phiwv(Mnum,Lnum,i+1,1)-Phiwv(Mnum,Lnum,i,1)
delm = Phiwv(Mnum,Lnum,i,1)  -Phiwv(Mnum,Lnum,i-1,1)
slop(i) = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )

delpgrd = Phigrdwv(Mnum,Lnum,i+1,1)-Phigrdwv(Mnum,Lnum,i,1)
delmgrd = Phigrdwv(Mnum,Lnum,i,1)  -Phigrdwv(Mnum,Lnum,i-1,1)
slopgrd(i) = dmax1( 0.d0,(2.d0*delpgrd*delmgrd+eps)/(delpgrd**2+delmgrd**2+eps) )
end do
do i = is-1,ie+1
Phiu(Mnum,Lnum,i) = Phiwv(Mnum,Lnum,i,1)- 0.5d0 * nu2 * ( Phiwv(Mnum,Lnum,i,1) - Phiwv(Mnum,Lnum,i-1,1))&
     + 0.25d0 * 1.d0 * slop(i) * ((1.0d0-slop(i)*1.d0/3.0d0)*(Phiwv(Mnum,Lnum,i,1)   -Phiwv(Mnum,Lnum,i-1,1)) + &
     (1.0d0+slop(i)*1.d0/3.0d0)   *(Phiwv(Mnum,Lnum,i+1,1) -Phiwv(Mnum,Lnum,i  ,1)))

Phiugrd(Mnum,Lnum,i) = Phigrdwv(Mnum,Lnum,i,1) + 0.5d0 * nu2 * ( Phigrdwv(Mnum,Lnum,i+1,1) - Phigrdwv(Mnum,Lnum,i,1))&
     - 0.25d0 * 1.d0 * slopgrd(i) * ((1.0d0+slopgrd(i)*1.d0/3.0d0)*(Phigrdwv(Mnum,Lnum,i  ,1) -Phigrdwv(Mnum,Lnum,i-1,1)) + &
     (1.0d0-slop(i)*1.d0/3.0d0)      *(Phigrdwv(Mnum,Lnum,i+1,1) -Phigrdwv(Mnum,Lnum,i  ,1)))
end do
end DO;end DO
DO Lnum = 1, ndy-2;DO Mnum = 1, ndx-2;do i = is,ie
Phiwv(Mnum,Lnum,i,1)    = Phiwv(Mnum,Lnum,i,1) - nu2 * (Phiu(Mnum,Lnum,i) - Phiu(Mnum,Lnum,i-1))
Phigrdwv(Mnum,Lnum,i,1) = Phigrdwv(Mnum,Lnum,i,1) + nu2 * (Phiugrd(Mnum,Lnum,i+1) - Phiugrd(Mnum,Lnum,i))
end do;end DO;end DO

call fapp_stop("loop7",1,0)

call fapp_start("loop8",1,0)
CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
LEFTt = LEFT!; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
RIGTt = RIGT!; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
Phiwv(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(1            ,-1,-1,idm),1,VECU,LEFT,1, &
Phiwv(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt

CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
      Phiwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
      Phiwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
TOP = TOPt; BOTM = BOTMt

CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
CALL MPI_SENDRECV(Phiwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
Phiwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_SENDRECV(Phiwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
 Phiwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
UP = UPt; DOWN = DOWNt


do k=0,ndz-1; do j=0,ndy-1; do i=0,ndx-1
grdxy1=adiff*Phiwv(i+1,j+1,k,1)+adiff*Phiwv(i-1,j-1,k,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,1)
grdyz1=adiff*Phiwv(i,j+1,k+1,1)+adiff*Phiwv(i,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,1)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,1)
grdzx1=adiff*Phiwv(i+1,j,k+1,1)+adiff*Phiwv(i-1,j,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,1)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,1)

Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
(-2.d0*cg*cg*grdxy1/dx1/dy1 &
 -2.d0*cg*cg*grdyz1/dy1/dz1 &
 -2.d0*cg*cg*grdzx1/dz1/dx1) *dt * dtration &
-G4pi*cg*cg*U(i,j,k,1)*dt * dtration
enddo; enddo; enddo


expand_dx=-2.d0*kappa * 0.5d0 * dt
expand_exp=dexp(expand_dx)
kp_i=1.0/(2.d0*kappa+1.d-10)
exp_m=(1.d0-expand_exp)
exp_p=(1.d0+expand_exp)
exp_k=exp_m*kp_i

!call fapp_start("loop1",1,0)
!do k=-1,ndz; do j=-1,ndy; do i=-1,ndx
do k=1,ndz-2; do j=1,ndy-2
!do ii=1,ndx-1,blki
!do i=ii,min(ii+blki-1,ndx-1)
do i=1,ndx-2
phiwv_d=Phiwv(i,j,k,1)
!enddo
!do i=ii,min(ii+blki-1,ndx-1)
!phiwv_d=Phiwv(i,j,k,1)
Phiwv(i,j,k,1)    = 0.5d0*Phiwv(i,j,k,1)*exp_p+Phigrdwv(i,j,k,1)*exp_k
Phigrdwv(i,j,k,1) = 0.5d0*Phigrdwv(i,j,k,1)*exp_p+0.5d0*kappa*phiwv_d*exp_m
!enddo
enddo
enddo; enddo
call fapp_stop("loop8",1,0)
enddo

!call fipp_stop
end subroutine slvmuscle
