
subroutine SELFGRAVWAVE(dt,mode)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
integer :: mode,MRANK,count=0
DOUBLE PRECISION  :: dt,dxi
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU
character(3) NPENUM
character(6) countcha
double precision tfluid , cs
double precision dt_mpi_gr(0:NPE-1),dt_gat_gr(0:NPE-1),maxcs,tcool,cgtime
double precision :: ave1,ave1pre,ave2(0:NPE-1),ave,avepre,ave2_gather(0:NPE-1) , eps=1.0d-3

!**************** INITIALIZEATION **************
if(mode==0) then
   Phi(:,:,:)=0.0d0
   Phidt(:,:,:)=0.0d0
end if
!**************** INITIALIZEATION **************



!****************read INITIAL CONDITION**************
if(mode==1) then
      WRITE(NPENUM,'(I3.3)') NRANK
      open(unit=8,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
      open(unit=18,file='/work/maedarn/3DMHD/test/PHIDTINI/INIPHIDT'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
      do k = -1, Ncellz+2
         do j = -1, Ncelly+2
            !do i = -1, Ncellx+2
            read(8) (Phi(i,j,k),i=-1,Ncellx+2)
      !enddo
         end do
      end do
      close(8)
      do k = -1, Ncellz+2
         do j = -1, Ncelly+2
            !do i = -1, Ncellx+2
            read(8) (Phidt(i,j,k),i=-1,Ncellx+2)
            !enddo
         end do
      end do
      close(18)
end if
!****************read INITIAL CONDITION**************



!****************GRAVITY SOLVER*****************

if(mode==2) then
  N_MPI(20)=1; N_MPI(1)=1
  iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)
  !call PB()
  !call mglin(Nmem1,Nmem2,2,5,5)
  !write(*,*) NRANK,Phi(0,0,0),Phi(1,1,1),Phi(Ncellx,Ncelly,Ncellz),dt,'-------33--333---33---'
  call gravslv(dt)
  !DEALLOCATE(bphi1,bphi2)
  if(NRANK==0) then
     write(*,*) NRANK,Phi(0,0,0),Phi(1,1,1),Phi(Ncellx,Ncelly,Ncellz),dt,'-------33-----33---'
     write(*,*) NRANK,Phi(0,0,0),Phidt(1,1,1),Phidt(Ncellx,Ncelly,Ncellz),dt,'-------33-----33---'
  end if
  N_ol = 2
  !*************一応***************
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !*************一応***************

!!!!!!!!!!!!!!!!!!!!!!!!!!!PHIDT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!も

                      !count, blocklength, stride
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  !LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL !x exact
  !RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL !x exact
  LEFTt = LEFT!; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL 周期
  RIGTt = RIGT!; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL 周期
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

  !*********use phi exact**********
  !if(IST.eq.0       ) then; do k=1,Ncellz; do j=1,Ncelly
  !  Phi(0       ,j,k) = Phi(1     ,j,k); Phi(-1       ,j,k) = Phi(1     ,j,k) !grad=0
  !end do; end do; end if
  !if(IST.eq.NSPLTx-1) then; do k=1,Ncellz; do j=1,Ncelly
  !  Phi(Ncellx+1,j,k) = Phi(Ncellx,j,k); Phi(Ncellx+2,j,k) = Phi(Ncellx,j,k)
  !end do; end do; end if
  !*********use phi exact**********
  !write(*,*) NRANK,Phi(0,0,0),Phi(1,1,1),Phi(Ncellx,Ncelly,Ncellz),'-------33-----33--66666-'
end if
!****************GRAVITY SOLVER*****************


!**********acceraration because of gravity******
if(mode==3) then !acceraration because of gravity
  dxi = 1.d0/(12.d0*dx(0))
  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
    U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
    U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
    U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0
  end do;end do;end do
end if
!**********acceraration because of gravity******

!***************SAVE PHI & PHIDT FOR DEBUG & INITIAL**************
if(mode==4) then
   !write(*,*) 'save???'
   WRITE(NPENUM,'(I3.3)') NRANK
   WRITE(countcha,'(I6.6)') count
   open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
   open(unit=38,file='/work/maedarn/3DMHD/test/PHIDTINI/INIPHIDT'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
   !write(*,*) 'save?????'

   !-------------------INITIAL---------------------
!   do k = -1, Ncellz+2
!      do j = -1, Ncelly+2
!         !do i = -1, Ncellx+2
!         write(28) (Phi(i,j,k),i=-1,Ncellx+2)
!         !enddo
!      end do
!   end do
!   close(28)
!   do k = -1, Ncellz+2
!      do j = -1, Ncelly+2
!         !do i = -1, Ncellx+2
!         write(38) (Phidt(i,j,k),i=-1,Ncellx+2) 
!         !enddo
!      end do
!   end do
!   close(38)
   !-------------------INITIAL---------------------


   !-------------------TEST---------------------
   do k = -1, Ncellz+2
      do j = -1, Ncelly+2
         !do i = -1, Ncellx+2
         write(28) (sngl(Phi(i,j,k)),i=-1,Ncellx+2)
         !enddo
      end do
   end do
   close(28)
   do k = -1, Ncellz+2
      do j = -1, Ncelly+2
         !do i = -1, Ncellx+2
         write(38) (sngl(Phidt(i,j,k)),i=-1,Ncellx+2)
         !enddo
      end do
   end do
   close(38)
   !-------------------TEST---------------------
   count=count+1
   !write(*,*) 'save????????????'
end if

!***************SAVE PHI & PHIDT FOR DEBUG & INITIAL**************


!***************SABILITY1**************
if(mode==5) then
   tfluid = dt/CFL
   write(*,*) '------------------tfluid-------------' , tfluid
   cs = deltalength/tfluid
   if(cs > cg) then
      write(*,*)
      write(*,*) '---------------------------------'
      write(*,*) '  err -------- cs > cg ' ,cs,cg , NRANK
      write(*,*) '---------------------------------'
      write(*,*)
   end if
   if(cs  > cg * 0.2d0) then
      write(*,*)
      write(*,*) '---------------------------------'
      write(*,*) '  caution ------- cg ' ,cs,cg , NRANK
      write(*,*) '---------------------------------'
      write(*,*)
   end if
end if
!***************SABILITY1**************

!***************SABILITY2**************
if(mode==7) then
   tcool = dt
   !cgtime = deltalength/cg * CFL
   cgtime = deltalength/cg * 0.2d0 !3D
   if(dt > cgtime) then
      write(*,*)
      write(*,*) '---------------------------------'
      write(*,*) ' dt > cgtime  ' , dt , cgtime , NRANK
      write(*,*) '---------------------------------'
      write(*,*)
      dt = cgtime
   end if
end if
!***************SABILITY2**************

!*****************cg****************
if(mode==6) then
   !tfluid = dt/CFL
   !cs = deltalength/tfluid
   tfluid=tfinal
   call Couran(tfluid)
   cs = deltalength/tfluid
   dt_mpi_gr(NRANK) = cs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! for MPI
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   CALL MPI_GATHER(dt_mpi_gr(NRANK),1,MPI_REAL8,   &
        dt_gat_gr       ,1,MPI_REAL8,   &
        0            ,MPI_COMM_WORLD,IERR)
!   CALL MPI_GATHER(st_mpi_gr(NRANK),1,MPI_INTEGER, &
!        st_gat_gr       ,1,MPI_INTEGER, &
!        0            ,MPI_COMM_WORLD,IERR)
   IF(NRANK.EQ.0)  THEN

      !---------------debug-------------------
      write(*,*)
      write(*,*)
      write(*,*) '-------------NRANK==0-----------',NRANK
      write(*,*)
      write(*,*)
      !---------------debug-------------------
      !dt  = tfinal
      !dtt = tfinal
      maxcs = 0.0d0
      do i_t = 0, NPE-1
         maxcs  = dmax1( maxcs, dt_gat_gr(i_t) )
         write(*,*) maxcs , '-----------maxcs-------'
!         if(dt.lt.dtt) st = st_gat(i_t)
!         dtt = dt
      end do
      cg = cgcsratio * maxcs
   END IF
   CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
   CALL MPI_BCAST(cg,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !cg = cgcsratio * maxcs
   write(*,*) NRANK,cg,'---------------cg------------------'
end if
!*****************cg****************

!*****************shuusoku****************
if(mode==8) then

   ave1=1.0d2
   ave1pre=1.0d2

!do k=1,Ncellz
!   do j=1,Ncelly
!      do i=1,Ncellx
         !ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k) + 1.0d-10) + ave1
!         ave1 = dabs(Phi(i,j,k)-Phidt(i,j,k)) + ave1
!      end do
!   end do
!end do

!do k=1,Ncellz
!   do j=1,Ncelly
!      do i=1,Ncellx
!         ave1pre=ave1
!         ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k)) + ave1
         !ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k) + 1.0d-10) + ave1
         !ave1 = dabs(Phi(i,j,k)-Phidt(i,j,k))
!         ave1 = dmax1( ave1pre , ave1 )
!      end do
!   end do
!end do

do k=1,Ncellz
   do j=1,Ncelly
      do i=1,Ncellx
         if(Phidt(i,j,k) .ne. 0.0d0) then
         ave1pre=ave1
         ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k)) + ave1
         !ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k) + 1.0d-10) + ave1
         !ave1 = dabs(Phi(i,j,k)-Phidt(i,j,k))
         ave1 = dmax1( ave1pre , ave1 )
         end if
      end do
   end do
end do

!---------------debug-------------------
!write(*,*) '-------------1---------3-----------',NRANK , ave1
!---------------debug-------------------

!ave1=ave1/dble((Ncellx*Ncelly*Ncellz))
ave2(NRANK)=ave1

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CALL MPI_GATHER(ave2(NRANK),1,MPI_REAL8,   &
     ave2_gather       ,1,MPI_REAL8,   &
     0            ,MPI_COMM_WORLD,IERR)


IF(NRANK.EQ.0)  THEN
   ave=1.0d2
   avepre=1.0d2
   !---------------debug-------------------
   !write(*,*)
   !write(*,*)
   !write(*,*) '-------------NRANK==000-----------',NRANK
   !write(*,*)
   !write(*,*)
   !---------------debug-------------------


   do i_t = 0, NPE-1
      !ave  = ave2_gather(i_t) + ave
      ave=dmax1(ave,ave2_gather(i_t))
      write(*,*) ave , '-----ave-----',i_t
   end do


   !ave=ave/dble(NPE)
   !---------------debug-------------------
   !write(*,*) '-------------1-ok??-----------',NRANK
   !---------------debug-------------------

   if(ave < eps) then
      shusoku1=1.0d2
   ELSE
      shusoku1=0.0d0
      !else !if(ave > eps) then
   end if
END IF

!---------------debug-------------------
!write(*,*) '-------------1-ok??*****-----------',NRANK
!---------------debug-------------------

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!CALL MPI_BCAST(maxcs,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!CALL MPI_BCAST(ave,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)
CALL MPI_BCAST(shuskou1,1,MPI_REAL8,0,MPI_COMM_WORLD,IERR)

!if(ave < eps) then  aveの初期条件で入ってしまう?
!   shusoku1=1.0d2
!else !if(ave > eps) then
!end if

!---------------debug-------------------
!write(*,*) '-------------1---------44-----------',NRANK
!---------------debug-------------------
end if
!*****************shuusoku****************
end subroutine SELFGRAVWAVE

subroutine gravslv(dt)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Phidummy !, Phidtdummy
double precision :: nu2 , w=6.0d0 , dt2 , dt
character(3) rn

nu2 = cg * dt / deltalength
nu2 = nu2 * nu2
!write(*,*) nu2,cg,dt,deltalength,'-----------???????------------' , NRANK
dt2 = dt * dt
write(rn,'(i3.3)') NRANK

ALLOCATE(Phidummy(-1:ndx,-1:ndy,-1:ndz))
!ALLOCATE(Phidummy(-1:ndx,-1:ndy,-1:ndz))
open(122,file='kakuninpost'//rn//'.dat')
open(123,file='kakuninpre'//rn//'.dat')

do k = -1 , Ncellz+2
   do j = -1 , Ncelly+2
      do i = -1 , Ncellx+2
         Phidummy(i,j,k) = Phi(i,j,k)
         !Phidtdummy(i,j,k) = Phidt(i,j,k)
         write(123,*) Phi(i,j,k)
      end do
   end do
end do
close(123)
do k = 1 , Ncellz
   do j = 1 , Ncelly
      do i = 1 , Ncellx
         Phi(i,j,k) = 2.0d0*Phidummy(i,j,k) - Phidt(i,j,k) + nu2 * (Phidummy(i-1,j,k) + Phidummy(i+1,j,k) + &
              Phidummy(i,j-1,k) + Phidummy(i,j+1,k) + Phidummy(i,j,k-1) + Phidummy(i,j,k+1) - w * Phidummy(i,j,k)) + &
              U(i,j,k,1) * G4pi * dt2
         write(122,*) Phi(i,j,k)
      end do
   end do
end do
close(122)

do k = -1 , Ncellz+2
   do j = -1 , Ncelly+2
      do i = -1 , Ncellx+2
         Phidt(i,j,k) = Phidummy(i,j,k)
      end do
   end do
end do

DEALLOCATE(Phidummy)
!DEALLOCATE(Phidtdummy)
end subroutine gravslv

subroutine gravslvMUSCL3D(dt,direction,mode)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Phidummy , gradPhidt !, Phidtdummy
double precision :: nu2 , w=6.0d0 , dt2 , dt, kappa , deltap,deltam
character(3) rn
DOUBLE PRECISION, dimension(:), allocatable :: grad1DPhi , Phi1D
!DOUBLE PRECISION, dimension(:), allocatable :: grad1DPhix , Phi1Dx , grad1DPhiy , Phi1Dy , grad1DPhiz , Phi1Dz
double precision dimension(:,:,:), allocatable :: intgral !intgrl(0:NPE-1)
integer direction , mode , invdt , loopmode , dloop , Ncell , Ncell1 , Ncell2
DOUBLE PRECISION, dimension(:,:), allocatable :: intC , Phisurface
DOUBLE PRECISION, dimension(:,:,:), allocatable :: f,g

nu2 = cg * dt / deltalength
!nu2 = nu2 * nu2
dt2 = dt * dt
invdt = 1.0d0 / dtpre

!-------------1st order---------------
!------------- dPhi/dt ----------------
!allocate(gradPhidt(-1:Ncellx+2,-1:Ncelly+2,-1:Ncellz+2))
!do k=-1,Ncellz+2
!   do j=-1,Ncelly+2
!      do i=-1,Ncellx+2
!         gradPhidt(i,j,k) = invdt * (Phi(i,j,k) - Phidt(i,j,k)) !dtは過去の
!      end do
!   end do
!end do
!------------- dPhi/dt ----------------

!-------------2nd order---------------
!------------- dPhi/dt ----------------
allocate(gradPhidt(-1:Ncellx+2,-1:Ncelly+2,-1:Ncellz+2))
do k=-1,Ncellz+2
   do j=-1,Ncelly+2
      do i=-1,Ncellx+2
         gradPhidt(i,j,k) = invdt * ( 3.0d0 * Phi(i,j,k) - 4.0d0 * Phidt(i,j,k) + Phi2dt(i,j,k)) !dtは過去の,後退差分
      end do
   end do
end do
!------------- dPhi/dt ----------------

ALLOCATE(Phi1D(-1:Ncell+2))
ALLOCATE(grad1DPhi(0:Ncell+1))
ALLOCATE(integral(1:Ncelly,1:Ncellz,0:NPE-1))
allocate(intC(1:Ncelly,1:Ncellz))
!allocate(Phisurface(1:Ncelly,1:Ncellz))





do k=-1,Ncellz+2
   do j=-1,Ncelly+2
      do i=-1,Ncellx+2
         gradPhidt(i,j,k) =  gradPhidt(i,j,k) * deltalength ! dphi/dt = a(x,y,z) , a(x,y,z) * dx   !constant mesh
      end do
   end do
end do




do dloop = 1,3 !three-direction

!--------------integral----------------
!--------------?台形法?-----------------OK cell center

   if(doolp==1) then
      Ncell=Ncellx
      Ncell1=Ncelly
      Ncell2=Ncellz
      allocate(Phisurface(1:Ncelly,1:Ncellz))
      do i = 1 , Ncell
         !Phi1D(i) = gradPhidt(i,j,k) * deltalength ! dphi/dt = a(x,y,z) , a(x,y,z) * dx
         Phisurface(j,k)=gradPhidt(i,j,k)
      end do
   elseif(dloop==2) then
      Ncell=Ncelly
      Ncell1=Ncellx
      Ncell2=Ncellz
      allocate(Phisurface(1:Ncellx,1:Ncellz))
      do j = 1 , Ncell
         !Phi1D(i) = gradPhidt(i,j,k) * deltalength
         Phisurface(i,k)=gradPhidt(i,j,k)
      end do
   else
      Ncell=Ncellz
      Ncell1=Ncellx
      Ncell2=Ncelly
      allocate(Phisurface(1:Ncellx,1:Ncelly))
      do k = 1 , Ncell
         !Phi1D(i) = gradPhidt(i,j,k) * deltalength
         Phisurface(i,j)=gradPhidt(i,j,k)
      end do
   end if
!--------------integral----------------
!------------integral-mpi--------------

integral(j,k,NRANK)=0.0d0
!do i=1,Ncell
   integral(j,k,NRANK) = Phisurface(i,j)
!end do


do Nroot=0,NPE-1
  CALL MPI_BCAST(integral(1,1,Nroot),(Ncell1)*(Ncell2),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do


intC(:,:)=0.0d0
do Nroot=0,NPE-1
   ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
   if(KST==KSTt .and. JST==JSTt) then
      intC(:,:) = integral(:,:,Nroot)+intC(:,:) !C=0
   end if
end do
!------------integral-mpi--------------



!------------split f,g-----------------
allocate(f(-1:Ncell+2,-1:Ncell+2,-1:Ncell+2))
allocate(g(-1:Ncell+2,-1:Ncell+2,-1:Ncell+2))
f(i,j,k) = 0.5d0 * Phi(i,j,k) + intgG(j,k) !x+cg*t
g(i,j,k) = 0.5d0 * Phi(i,j,k) - intgG(j,k) !x-cg*t
!------------split f,g-----------------


!-------------MUSCL solver-------------
kappa=1.0d0/3.0d0
!ul(i,j,k) = u(i,j,k) + 0.25d0 * s * ((1-kappa*s)*gradum + (1+kappa*s)*gradup) !j-1
!ur(i,j,k) = u(i,j,k) - 0.25d0 * s * ((1-kappa*s)*gradup + (1+kappa*s)*gradum) !j+1
deltap = f(i+2,j,k) - f(i+1,j,k)
deltam = f(i+1,j,k) - f(i  ,j,k)
fluxf(i,j,k) = f(i+1,j,k) - gradf * 0.25d0 *( (1.0d0 - kappa * gradf) * deltap + (1.0d0 + kappa * gradf) * deltam)  !ur_{j+1/2}
!fluxf(i,j,k) = - cg * fluxf(i,j,k)

!fluxf = 0.5d0 * 

!f(i,j,k) = f(i,j,k) - nu2 * 0.5d0 * (fluxf(i,j,k) - fluxf(i-1,j,k)) !dt/2 timestep
fdt2(i,j,k) = f(i,j,k) - nu2 * 0.5d0 * (fluxf(i,j,k) - fluxf(i-1,j,k)) !dt/2 timestep


deltap = g(i+1,j,k) - g(i  ,j,k)
deltam = g(i  ,j,k) - g(i-1,j,k)
fluxg(i,j,k) = g(i  ,j,k) + gradg * 0.25d0 *( (1.0d0 - kappa * gradg) * deltam + (1.0d0 + kappa * gradg) * deltap)  !ul_{j+1/2}
fluxg(i,j,k) =   cg * fluxg(i,j,k)


!fluxr = 0.5d0 * cg * (ul(i,j,k)+ur(i,j,k) + ul(i,j,k) - ur(i,j,k)) !向きによって片方で良い
!fluxl = 0.5d0 * cg * (ul(i,j,k)+ur(i,j,k) - ul(i,j,k) + ur(i,j,k)) !向きによって片方で良い


!u(i,j,k) = u(i.j,k) - dt/dx * (fluxr-fluxl)
!-------------MUSCL solver-------------


!------------newPHI------------------
Phi(i,j,k) = f(i,j,k) + g(i,j,k)
!------------newPHI------------------


end subroutine gravslvMUSCL1D




subroutine slvgrvsource(dt,gradPhidt)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'
  double precision dt

  Phi(i,j,k) = 0.5d0 * G4pi * U(i,j,k,1) * dt * dt + dt * gradPhidt(i,j,k) + Phi(i,j,k)

end subroutine slvgrvsource


subroutine INTEGRAL1D(NCELL1,NCELL2,NCELL3,loop,gradPhidt,Phisurface)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'
  integer NCELL1,NCELL2,NCELL3,loop
  double precision gradPhidt(-1:Ncellx+2,-1:Ncelly+2,-1:Ncellz+2)
  double precision Phisurface(1:NCELL2,1:NCELL3)
  do i = 1 , Ncell
     Phisurface(j,k)=gradPhidt(i,j,k)
  end do

  integral(j,k,NRANK) = Phisurface(i,j)
  do Nroot=0,NPE-1
     CALL MPI_BCAST(integral(1,1,Nroot),(Ncell1)*(Ncell2),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
  end do


  intC(:,:)=0.0d0
  if(loop==1) then
     do Nroot=0,NPE-1
        ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
        if(KST==KSTt .and. JST==JSTt .and. ISTt < IST) then
           intC(:,:) = integral(:,:,Nroot)+intC(:,:) !C=0
        end if
     end do
     do i = 1 , Ncell
       gradPhidt(i,:,:) = gradPhidt(i,:,:) + intC(:,:)
     end do
  end if
  if(loop==2) then
     do Nroot=0,NPE-1
        ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
        if(KST==KSTt .and. IST==ISTt .and. JSTt < JST) then
           intC(:,:) = integral(:,:,Nroot)+intC(:,:) !C=0
        end if
     end do
     do j = 1 , Ncell
       gradPhidt(:,j,:) = gradPhidt(:,j,:) + intC(:,:)
     end do
  end if
  if(loop==1) then
     do Nroot=0,NPE-1
        ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
        if(IST==ISTt .and. JST==JSTt .and. KSTt < KST) then
           intC(:,:) = integral(:,:,Nroot)+intC(:,:) !C=0
        end if
     end do
  end if


!subroutine gradforintegral(Phi1D,grad1DPhi,len)
!USE comvar
!USE mpivar
!USE slfgrv
!INCLUDE 'mpif.h'
!  integer len
!  double precision grad1DPhi(0:len),Phi1D(-1,len+1)
  !double precision Phi1Dr(-1,len+1),Phi1Dl(-1,len+1)
!  double precision invdlen
!  invdlen = 1.0d0 / deltalength

!  do i = 0 , len
!     grad1DPhi(i) = 0.5d0 * (Phi1D(i-1) + Phi1D(i+1)) * invdlen
!  end do
!end subroutine gradforintegral


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
character(3) fn,deep
!character(3) lname
character(2) lcRANK
character(1) lRANK
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

nl = NGcr-1
!nl = NGcr  !===================????==================
pointb2(nl) = 1
nx=(2**NGcr)/NSPLTx+2; ny=(2**NGcr)/NSPLTy+2; nz=(2**NGcr)/NSPLTz+2
3 continue
pointb2(nl+1)=pointb2(nl)+max(nx*ny,ny*nz,nz*nx)
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

!***************fordebug*****************
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!***************fordebug*****************

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

  !***************fordebug*****************
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !***************fordebug*****************

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

  call rlft3(data(1,1,klr),speq(1,klr),nn1,nn2,1)

  kz = klr
  zp1 = x(kz)-0.5d0*dzz
  zp2 = Lbox - zp1
  temp1r = dat1(1,1) - data(1,1,klr) * 0.5d0*zp1 * facG
  temp1i = dat1(2,1) - data(2,1,klr) * 0.5d0*zp1 * facG
  temp2r = dat2(1,1) - data(1,1,klr) * 0.5d0*zp2 * facG
  temp2i = dat2(2,1) - data(2,1,klr) * 0.5d0*zp2 * facG

  do m=1,nn2/2+1; do l=1,nn1/2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)+1.0d-100
    !kap = sqrt(kap)
    dat1(2*l-1,m) = dat1(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat1(2*l  ,m) = dat1(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat2(2*l-1,m) = dat2(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    dat2(2*l  ,m) = dat2(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    !write(*,*) dat1(2*l-1,m),dat1(2*l  ,m),dat2(2*l-1,m),dat2(2*l  ,m),'PPPPPPPPPPPPPPPPP'
  end do;end do

  l=nn1/2+1
  do m=1,nn2/2+1
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)+1.0d-100
    !kap = sqrt(kap)
    spe1(m) = spe1(m) + speq(m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    spe2(m) = spe2(m) + speq(m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do

  do m2=nn2/2+2,nn2; m=nn2+2-m2; do l=1,nn1/2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)+1.0d-100
    !kap = sqrt(kap)
    dat1(2*l-1,m2) = dat1(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat1(2*l  ,m2) = dat1(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat2(2*l-1,m2) = dat2(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    dat2(2*l  ,m2) = dat2(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do;end do

  l=nn1/2+1
  do m2=nn2/2+2,nn2; m=nn2+2-m2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)+1.0d-100
    !kap = sqrt(kap)
    spe1(m2) = spe1(m2) + speq(m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    spe2(m2) = spe2(m2) + speq(m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
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

!write(*,*) 'PB1'

!write(fn,'(i3.3)') NRANK
!open(30,file='bpl'//fn//'.DAT')
!open(560,file='bpr'//fn//'.DAT')

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
!  write(30,*) bphi2(pointb2(NGL)+n,1)
!  write(560,*) bphi2(pointb2(NGL)+n,2)
end do
!write(30,*)
!write(560,*)
end do

!close(30)
!close(560)

!write(*,*) 'PB2'

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

ncx = Ncellx+1; ncy = Ncelly+1; ncz = Ncellz+1
w  = 0.125d0
do lc=NGL-1,NGcr-1,-1
  nfx = ncx; nfy = ncy; nfz = ncz
  ncx = ncx/2+1; ncy = ncy/2+1; ncz = ncz/2+1
!  write(lcRANK,'(i2.2)') lc
  do l = 1, 2
  if(l.eq.1) then; mfx=nfy; mfy=nfz; mcx=ncy; mcy=ncz; end if
  if(l.eq.2) then; mfx=nfy; mfy=nfz; mcx=ncy; mcy=ncz; end if

!     write(lRANK,'(i1.1)') l
!     open(211,file='bpl'//lRANK//lcRANK//fn//'.dat')

  do jc = 2,mcy-1; jf = 2*jc-1
  do ic = 2,mcx-1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 4.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1    ,l)+bphi2(pointb2(lc+1)+kf+1    ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
!    write(211,*)  bphi2(pointb2(lc)+kc,l) ,lc,ncx,'pos1'
  end do
  end do

  do jc = 2,mcy-1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
!     write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos2'
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
!     write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos3'
  end do
  do ic = 2,mcx-1; if = 2*ic-1
    jc = 1; jf = 2*jc-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
!     write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos4'
    jc = mcy; jf = 2*jc-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )
!    write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos5'
  end do

    jc = 1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
!    write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos6'
    jc = 1 ; jf = 2*jc-1
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
!    write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos7'
    jc = mcy; jf = 2*jc-1
    ic = 1 ; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )
!    write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos8'
    jc = mcy; jf = 2*jc-1
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1 ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )

!    write(211,*)  bphi2(pointb2(lc)+kc,l) ,'pos9'
!    close(211)
  end do
end do

goto 977
nx=(2**NGL)/NSPLTx+2
do lc=NGL-1,NGcr-1,-1
   write(lcRANK,'(i2.2)') lc
   do l=1,2
      if(l==1) then
         open(211,file='bpl'//lcRANK//fn//'.dat')
         do i=1,nx*nx
            write(211,*) bphi2(pointb2(lc)+i,l)
         end do
      end if
      if(l==2) then
         open(201,file='bpr'//lcRANK//fn//'.dat')
         do i=1,nx*nx
            write(201,*) bphi2(pointb2(lc)+i,l)
         end do
      end if
   end do
   nx=nx/2+1
end do
977 continue
!write(*,*) 'PB3'

nz=(2**NGcr)/NSPLTz+1; ny=(2**NGcr)/NSPLTy+1; nx=(2**NGcr)/NSPLTx+1; n1=2**NGcr+1
ALLOCATE( temp1(n1,n1,n1), temp2(0:nx,0:ny,0:nz) )

!write(fn,'(i3.3)') NRANK
!open(40,file='bp2l'//fn//'.DAT')
!open(570,file='bp2r'//fn//'.DAT')

!open(50,file='bpkaku2l'//fn//'.DAT')
!open(580,file='bpkaku2r'//fn//'.DAT')

!write(*,*)  pointb2(NGcr-1) , pointb2(NGcr) , pointb1(NGcr) , n1 , nz ,pointb2(NGcr)-pointb2(NGcr-1),pointb1(NGcr+1)-pointb1(NGcr) ,'PBPB'

do k=0,nz; kk = (ny+1)*k + pointb2(NGcr) !original + pointb2(NGcr)
do j=0,ny; n = j + kk
   temp2(1 ,j,k) = bphi2(n,1)
!   write(50,*) temp2(1 ,j,k)
end do;end do
call collect( temp1,temp2,n1,n1,n1,nx,ny,nz )
do k=1,n1; kk = n1*(k-1) + pointb1(NGcr)
do j=1,n1; n = j-1 + kk
   bphi1(n,1) = temp1(1 ,j,k)
!   write(40,*) bphi1(n,1)
end do;end do

!write(*,*) 'PB4'

do k=0,nz; kk = (ny+1)*k + pointb2(NGcr) !original + pointb2(NGcr)
do j=0,ny; n = j + kk
   temp2(nx,j,k) = bphi2(n,2)
!   write(580,*) temp2(nx,j,k)
end do;end do
call collect( temp1,temp2,n1,n1,n1,nx,ny,nz )
do k=1,n1; kk = n1*(k-1) + pointb1(NGcr)
do j=1,n1; n = j-1 + kk
   bphi1(n,2) = temp1(n1,j,k)
!   write(570,*) bphi1(n,2)
end do;end do


!close(40)
!close(570)
!close(50)
!close(580)

DEALLOCATE(temp1,temp2)

!write(*,*) 'PB5'

nc = (2**NGcr)+1
w  = 0.125d0
do lc=NGcr-1,1,-1
!   write(lcRANK,'(i2.2)') lc
  nf = nc
  nc = nc/2+1
  do l = 1, 2

!     write(lRANK,'(i1.1)') l
!     open(231,file='bpr'//lRANK//lcRANK//fn//'.dat')

  do jc = 2,nc-1; jf = 2*jc-1
  do ic = 2,nc-1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 4.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
!     write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos1'
  end do
  end do

  do jc = 2,nc-1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
!     write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos2'
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
!     write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos3'
  end do
  do ic = 2,nc-1; if = 2*ic-1
    jc = 1; jf = 2*jc-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
!    write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos4'
    jc = nc; jf = 2*jc-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )
!     write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos5'
  end do

    jc = 1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
!     write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos6'
    jc = 1 ; jf = 2*jc-1
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
!     write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos7'
    jc = nc; jf = 2*jc-1
    ic = 1 ; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )
!    write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos8'
    jc = nc; jf = 2*jc-1
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )

!    write(231,*)  bphi1(pointb1(lc)+kc,l) ,lc,nc,'pos9'
!    close(231)
  end do
end do




END SUBROUTINE PB


SUBROUTINE rlft3(data,speq,nn1,nn2,isign)
INTEGER isign,nn1,nn2
COMPLEX*16 data(nn1/2,nn2),speq(nn2)
INTEGER i1,i2,j1,j2,nn(2)
DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 
COMPLEX*16 c1,c2,h1,h2,w
c1=cmplx(0.5d0,0.0d0) 
c2=cmplx(0.0d0,-0.5d0*isign) 
theta=6.28318530717959d0/dble(isign*nn1) 
wpr=-2.0d0*sin(0.5d0*theta)**2 
wpi=sin(theta) 
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
  ntot=ntot*nn(idim)
enddo
nprev=1 
do idim=1,ndim
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
END SUBROUTINE fourn


SUBROUTINE collect(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2)
USE mpivar
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: u1(nx1,ny1,nz1),u2(0:nx2,0:ny2,0:nz2)
double precision :: tMPI(0:nx2,0:ny2,0:nz2,0:NPE-1)

!***************fordebug*****************
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!***************fordebug*****************

do k=0,nz2; do j=0,ny2; do i=0,nx2
  tMPI(i,j,k,NRANK)=u2(i,j,k)
end do;end do;end do
do Nroot=0,NPE-1
  CALL MPI_BCAST(tMPI(0,0,0,Nroot),(nx2+1)*(ny2+1)*(nz2+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
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
 !********debug**********
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!********debug**********
END SUBROUTINE collect
