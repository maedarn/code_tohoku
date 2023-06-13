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
