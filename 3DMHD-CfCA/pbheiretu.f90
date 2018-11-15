MODULE comvar
!INTEGER, parameter :: ndx=130, ndy=130, ndz=130, ndmax=130, Dim=3 !1024^3
!INTEGER, parameter :: ndx=66, ndy=66, ndz=66, ndmax=66, Dim=3 !512^3
INTEGER, parameter :: ndx=34, ndy=34, ndz=34, ndmax=34, Dim=3
!INTEGER, parameter :: ndx=18, ndy=18, ndz=18, ndmax=18, Dim=3
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
INTEGER :: ifchem,ifthrm,ifrad,ifgrv,klrmax,loopbc=3
END MODULE comvar

MODULE mpivar
INTEGER :: NPE,NRANK, NSPLTx,NSPLTy,NSPLTz, IST,JST,KST, LEFT,RIGT,BOTM,TOP,UP,DOWN
INTEGER :: BCx1,BCx2,BCy1,BCy2,BCz1,BCz2, N_MPI(20)
!DOUBLE PRECISION  :: BBRV(10,2,2),BBRV_cm(8)
!REAL*4, dimension(:,:,:), allocatable :: DTF
END MODULE mpivar


MODULE slfgrv
DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
INTEGER :: point1(0:15),point2(0:15),NGL,NGcr,Nmem1,Nmem2
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Phi
DOUBLE PRECISION :: Lbox

INTEGER :: pointb1(0:15),pointb2(0:15)
DOUBLE PRECISION, dimension(:,:), allocatable :: bphi1,bphi2
END MODULE slfgrv

PROGRAM MAIN
USE comvar
USE mpivar
!USE chmvar
USE slfgrv
INCLUDE 'mpif.h'

!CALL MPI_INIT(IERR)
!CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPE  ,IERR)
!CALL MPI_COMM_RANK(MPI_COMM_WORLD,NRANK,IERR)
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

!double precision,  parameter :: pi = 3.14159265359d0
!DOUBLE PRECISION, dimension(:,:,:), allocatable :: temp1,temp2

DOUBLE PRECISION, dimension(:,:,:),  allocatable :: data
!complex*16, dimension(:,:), allocatable :: speq
!DOUBLE PRECISION :: rho

!DOUBLE PRECISION, dimension(:,:),  allocatable :: dat1,dat2
!complex*16, dimension(:), allocatable :: spe1,spe2
!double precision :: kap,temp1r,temp1i,temp2r,temp2i,facG,fac,dxx,dyy,dzz,zp1,zp2
!double precision, dimension(:,:,:), allocatable :: fint0,fint1

character*4 fnum
character*2 fn

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

!------input------
   Ncellx=32
   Ncelly=32
   Ncellz=32
!------input------
!write(*,*) 'OK1'

IST = mod(NRANK,NSPLTx); KST = NRANK/(NSPLTx*NSPLTy); JST = NRANK/NSPLTx-NSPLTy*KST
!write(*,*) 'OK11'
LEFT = NRANK - 1            ; if(IST.eq.0       ) LEFT = NRANK + (NSPLTx-1)
RIGT = NRANK + 1            ; if(IST.eq.NSPLTx-1) RIGT = NRANK - (NSPLTx-1)
BOTM = NRANK - NSPLTx       ; if(JST.eq.0       ) BOTM = NRANK + NSPLTx*(NSPLTy-1)
TOP  = NRANK + NSPLTx       ; if(JST.eq.NSPLTy-1) TOP  = NRANK - NSPLTx*(NSPLTy-1)
DOWN = NRANK - NSPLTx*NSPLTy; if(KST.eq.0       ) DOWN = NRANK + NSPLTx*NSPLTy*(NSPLTz-1)
UP   = NRANK + NSPLTx*NSPLTy; if(KST.eq.NSPLTz-1) UP   = NRANK - NSPLTx*NSPLTy*(NSPLTz-1)

ALLOCATE( U(-1:ndx,-1:ndy,-1:ndz,8) )
!ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1))
ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:Ncellx+2))
WRITE(fnum,'(I4.4)') NRANK
do k=1,Ncellz
do j=1,Ncelly
do i=-1,Ncellx+2
   !U(i,j,k,1) = dble(NRANK*Ncellx*Ncelly*Ncellz+i+j*Ncelly+k*Ncellz)
   !U(i,j,k,1) = dble(NRANK)
   U(i,j,k,1) = dble(NRANK + i*1000)
end do;end do;end do

open(10,FILE='/work/maedarn/3DMHD/test/init'//fnum//'.dat')
do i=-1,Ncellx+2
nccy = Ncelly; nccz = Ncellz
do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
!do i=-1,Ncellx+2
   data(jy,kz,i) = i + jy + kz !U(i,j,k,1)
   write(10,*) data(jy,kz,i)
end do;end do
write(10,*)
write(10,*)
end do
close(10)

!write(*,*) 'ok1'
                    !count,blocklength,stride
CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)

!open(18,FILE='/work/maedarn/3DMHD/test/snrs'//fnum//'.dat')
do nlp2 = 0 , loopbc , 1
!KST = (NPE-1)/(NSPLTx*NSPLTy); JST = NRANK/NSPLTx-NSPLTy*KST
   klrmax = ((NSPLTy) + NSPLTy * (NSPLTz-1)) * nlp2
   !write(*,*) klrmax
do Nlp = 1,NSPLTy*NSPLTz-1

  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

  Nis = JSs + NSPLTy*KSs
  !kls = Nis + 1
  kls = Nis - 1 + klrmax
  Nir = JST + NSPLTy*KST
  !klr = Nir + 1
  klr = Nir - 1 + klrmax

  if(kls.gt.Ncellx+2) then; isend = MPI_PROC_NULL; kls = Ncellx+3; end if
  if(klr.gt.Ncellx+2) then; irecv = MPI_PROC_NULL; klr = Ncellx+3; end if
     !write(*,*) kls,klr,isend,irecv,KSs,JSs,KSr,JSr,Nlp,NRANK!,&
     write(*,*) klr+1,NRANK!,&
!          JST*Ncelly+1,KST*Ncellz+1,kls,JSr*Ncelly+1,KSr*Ncellz+1,klr
  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls),1,VECU,isend,1, & !send    VECU is already created.
       data(JSr*Ncelly+1,KSr*Ncellz+1,klr),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv     !! no relation send buf !!
!  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls-1),1,VECU,isend,1, & !send    VECU is already created.
!       data(JSr*Ncelly+1,KSr*Ncellz+1,klr-1),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv     !! no relation send buf !!
!  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls-2),1,VECU,isend,1, & !send    VECU is already created.
!       data(JSr*Ncelly+1,KSr*Ncellz+1,klr-2),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv     !! no relation send buf !!
!  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls+1),1,VECU,isend,1, & !send    VECU is already created.
!       data(JSr*Ncelly+1,KSr*Ncellz+1,klr+1),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv     !! no relation send buf !!
!  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls+2),1,VECU,isend,1, & !send    VECU is already created.
!       data(JSr*Ncelly+1,KSr*Ncellz+1,klr+2),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv     !! no relation send buf !!
end do
!end do
!close(18)
if(klr.le.Ncellx+2) then
WRITE(fn,'(I2.2)') klr+1
open(19,FILE='/work/maedarn/3DMHD/test/final'//fn//fnum//'.dat')
end if
!do i=-1,Ncellx*NSPLTx+2
do k=1,Ncellz*NSPLTz
do j=1,Ncelly*NSPLTy
!do i=1,Ncellx*NSPLTx
   write(19,*)  data(j,k,klr)
end do;end do!;end do
if(klr.ne.Ncellx+3) then
close(19)
end if
end do
CALL MPI_TYPE_FREE(VECU,IERR)
!end do
DEALLOCATE(U,data)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL MPI_FINALIZE(IERR)
END program main
