subroutine SELFGRAVWAVE(dt,mode)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'
  integer :: mode,MRANK,count=0,rdnum,rddmy
  DOUBLE PRECISION  :: dt,dxi
  !INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt
  INTEGER :: MSTATUS(MPI_STATUS_SIZE)
  DOUBLE PRECISION  :: VECU
  real(4) :: dmy(1:28)
  character(3) NPENUM
  character(6) countcha
  double precision tfluid , cs
  double precision dt_mpi_gr(0:NPE-1),dt_gat_gr(0:NPE-1),maxcs,tcool,cgtime!,sourcedt
  double precision :: ave1,ave1pre,ave2(0:NPE-1),ave,avepre,ave2_gather(0:NPE-1) , eps=1.0d-3
  !double precision , dimension(:,:,:) , allocatable :: stbPhi
  !double precision , dimension(-1:Ncellx+2,-1:Ncelly,-1:Ncellz) :: Phipregrad,Phipregraddum
  !**************** INITIALIZEATION **************
  if(mode==0) then
     Phi(:,:,:) = 0.0d0
     Phiwv(:,:,:,:)=0.d0
     Phigrdwv(:,:,:,:)=0.d0
  end if
  !**************** INITIALIZEATION **************

  !****************read INITIAL CONDITION**************
  if(mode==1) then
     goto 5401
     WRITE(NPENUM,'(I3.3)') NRANK
     open(unit= 8,file=dir//'PHIINI/PHIINIwv'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=18,file=dir//'PHIINI/PHIGRDwv'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           do i = -1, Ncellx+2
           read(8) (Phiwv(i,j,k,rdnum),rdnum=1,wvnum)
           enddo
        end do
     end do
     close(8)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           do i = -1, Ncellx+2
           read(18) (Phigrdwv(i,j,k,rdnum),rdnum=1,wvnum)
           enddo
        end do
     end do
     close(18)
     5401 continue

     count=100
     WRITE(NPENUM,'(I3.3)') NRANK
     WRITE(countcha,'(I6.6)') count
     !open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIDTINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=8,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgp'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=18,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgm'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file=dir//'PHIINI/PHI'//countcha//NPENUM//'.DAT',FORM='UNFORMATTED')
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           do i = -1, Ncellx+2
           read(28) (dmy(rddmy),rddmy=1,28)

Phiwv(i,j,k,1) = dble(dmy(1))
Phiwv(i,j,k,2) = dble(dmy(2))
Phiwv(i,j,k,3) = dble(dmy(3))
Phiwv(i,j,k,4) = dble(dmy(4))
Phiwv(i,j,k,5) = dble(dmy(5))
Phiwv(i,j,k,6) = dble(dmy(6))
Phiwv(i,j,k,7) = dble(dmy(7))
Phiwv(i,j,k,8) = dble(dmy(8))
Phigrdwv(i,j,k,1) = dble(dmy(9))
Phigrdwv(i,j,k,2) = dble(dmy(10))
Phigrdwv(i,j,k,3) = dble(dmy(11))
Phigrdwv(i,j,k,4) = dble(dmy(12))
Phigrdwv(i,j,k,5) = dble(dmy(13))
Phigrdwv(i,j,k,6) = dble(dmy(14))
Phigrdwv(i,j,k,7) = dble(dmy(15))
Phigrdwv(i,j,k,8) = dble(dmy(16))
Phiexa(i,j,k) = dble(dmy(17))
Phigrd(i,j,k,1) = dble(dmy(18))
Phigrd(i,j,k,2) = dble(dmy(19))
Phigrd(i,j,k,3) = dble(dmy(20))
Phigrd(i,j,k,4) = dble(dmy(21))
Phigrd(i,j,k,5) = dble(dmy(22))
Phigrd(i,j,k,6) = dble(dmy(23))
Phigrd(i,j,k,7) = dble(dmy(24))
Phigrd(i,j,k,8) = dble(dmy(25))
Phiexab1(i,j,k) = dble(dmy(26))
Phiexab2(i,j,k) = dble(dmy(27))
U(i,j,k,1) = dble(dmy(28))
          enddo
        end do
     end do
     close(28)
  end if
  !****************read INITIAL CONDITION**************

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
     call slvmuscle(dt)
  end if

  !****************GRAVITY SOLVER*****************

  !**********acceraration because of gravity******
  if(mode==3) then !acceraration because of gravity
     iwx = 1; iwy = 1; iwz = 1
     !call BCgrv(101)
     !call BCgrv(102)
     dxi = 1.d0/(12.d0*dx(0))
     do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
        !U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
        !U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
        !U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0

        !U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
        !U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
        !U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0
     end do;end do;end do
  end if
  !**********acceraration because of gravity******

  !***************SAVE PHI & PHIDT FOR DEBUG & INITIAL**************
  if(mode==4) then
     !write(*,*) 'save???'
     WRITE(NPENUM,'(I3.3)') NRANK
     WRITE(countcha,'(I6.6)') count
     !open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIDTINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=8,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgp'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=18,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgm'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file=dir//'PHIINI/PHI'//countcha//NPENUM//'.DAT',access='stream',FORM='UNFORMATTED')
     !open(unit=29,file=dir//'PHIINI/FPHI'//countcha//NPENUM//'.DAT',FORM='FORMATTED') 
     !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI2step'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !write(*,*) 'save?????'

     !-------------------INITIAL---------------------
     !   do k = -1, Ncellz+2
     !      do j = -1, Ncelly+2
     !         do i = -1, Ncellx+2
     !         write(28) (Phiwv(i,j,k,rdnum),rdnum=1,wvnum)
     !         enddo
     !      end do
     !   end do
     !   close(28)
     !   do k = -1, Ncellz+2
     !      do j = -1, Ncelly+2
     !         do i = -1, Ncellx+2
     !         write(38) (Phigrdwv(i,j,k,rdnum),rdnum=1,wvnum)
     !         enddo
     !      end do
     !   end do
     !   close(38)
     !-------------------INITIAL---------------------

     !-------------------TEST---------------------
     !iwx = 1; iwy = 1; iwz = 1
     !call BCgrv(101)
     !call BCgrv(102)
     !Phiwv(:,:,:,:)=dble(NRANK)
     !do k = -2, Ncellz+3
     !   do j = -2, Ncelly+3
     !      do i = -2, Ncellx+3
     !      Phiwv(i,j,k,1)=dble(NRANK)
     !      Phiwv(i,j,k,2)=dble(NRANK)
     !      enddo
     !   end do
     !end do

     iwx=1;iwy=1;iwz=1
     call BCgrv(100,1,8)
     do k = -2, Ncellz+3
        do j = -2, Ncelly+3
           do i = -2, Ncellx+3
           write(28) sngl(Phiwv(i,j,k,1)),sngl(Phiwv(i,j,k,2)),&!sngl(Phiwv(i,j,k,3)),sngl(Phiwv(i,j,k,4)),&!sngl(Phiwv(i,j,k,5)),sngl(Phiwv(i,j,k,6)),&
                   !sngl(Phiwv(i,j,k,7)),sngl(Phiwv(i,j,k,8)),
                 sngl(Phigrdwv(i,j,k,1)),&!,sngl(Phigrdwv(i,j,k,2)),sngl(Phigrdwv(i,j,k,3)),sngl(Phigrdwv(i,j,k,4)),&
               !sngl(Phigrdwv(i,j,k,5)),sngl(Phigrdwv(i,j,k,6)),sngl(Phigrdwv(i,j,k,7)),sngl(Phigrdwv(i,j,k,8))
                      sngl(Phiexa(i,j,k)),sngl(Phigrd(i,j,k,1)),sngl(U(i,j,k,1))!,&
               !sngl(Phigrd(i,j,k,2)),sngl(Phigrd(i,j,k,3)),sngl(Phigrd(i,j,k,4)),sngl(Phigrd(i,j,k,5)),sngl(Phigrd(i,j,k,6)),sngl(Phigrd(i,j,k,7)),sngl(Phigrd(i,j,k,8)),&
               !sngl(Phiexab1(i,j,k)),sngl(Phiexab2(i,j,k))
     !      write(29,*) sngl(Phiwv(i,j,k,1)),sngl(Phiwv(i,j,k,2)),&!sngl(Phiwv(i,j,k,3)),sngl(Phiwv(i,j,k,4)),&!sngl(Phiwv(i,j,k,5)),sngl(Phiwv(i,j,k,6)),&
                   !sngl(Phiwv(i,j,k,7)),sngl(Phiwv(i,j,k,8)),
     !            sngl(Phigrdwv(i,j,k,1)),&!,sngl(Phigrdwv(i,j,k,2)),sngl(Phigrdwv(i,j,k,3)),sngl(Phigrdwv(i,j,k,4)),&
               !sngl(Phigrdwv(i,j,k,5)),sngl(Phigrdwv(i,j,k,6)),sngl(Phigrdwv(i,j,k,7)),sngl(Phigrdwv(i,j,k,8))
     !                 sngl(Phiexa(i,j,k)),sngl(Phigrd(i,j,k,1)),sngl(U(i,j,k,1))!,&
               !sngl(Phigrd(i,j,k,2)),sngl(Phigrd(i,j,k,3)),sngl(Phigrd(i,j,k,4)),sngl(Phigrd(i,j,k,5)),sngl(Phigrd(i,j,k,6)),sngl(Phigrd(i,j,k,7)),sngl(Phigrd(i,j,k,8)),&
               !sngl(Phiexab1(i,j,k)),sngl(Phiexab2(i,j,k))
          enddo
        end do
     end do
     close(28)
     !close(29)
!     do k = -1, Ncellz+2
!        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
!           write(38) (sngl(Phidt(i,j,k)),i=-1,Ncellx+2)
           !enddo
!        end do
!     end do
!     close(38)
     !-------------------TEST---------------------
     count=count+1
     !write(*,*) 'save????????????'
  end if

  !***************SAVE PHI & PHIDT FOR DEBUG & INITIAL**************

  !***************SABILITY1**************
if(mode==5) then
   tfluid = dt/CFL
   write(*,*) '------------------tfluid-------------' , tfluid
   !cs = deltalength/tfluid
   cs = dx(1)/tfluid
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
   !cgtime = deltalength/cg * 0.2d0 !3D
   cgtime = dx(1)/cg * cgratio1 !3D-isotropic
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
   !cs = deltalength/tfluid
   cs = dx(1)/tfluid !3D-isotropic
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
!      cg = cgcsratio * maxcs
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

   ave1=1.0d2  ! これはおかしい
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



!do k=1,Ncellz
!   do j=1,Ncelly
!      do i=1,Ncellx
!         if(Phidt(i,j,k) .ne. 0.0d0) then
!         ave1pre=ave1
!         ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k))! + ave1
         !ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k) + 1.0d-10) + ave1
         !ave1 = dabs(Phi(i,j,k)-Phidt(i,j,k))
!         ave1 = dmax1( ave1pre , ave1 )
!         end if
!      end do
!   end do
!end do

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



!***************SABILITY-exa**************
if(mode==9) then
   !sourcedt = dt:
   !cgtime = deltalength/cg * CFL
   !call STBLphi(sourcedt,Phipregrad)
   !dt = sourcedt
end if
!***************SABILITY-exa**************
end subroutine SELFGRAVWAVE


subroutine slvmuscle(dt)
  use comvar
  use slfgrv
  use mpivar
  INCLUDE 'mpif.h'
  double precision :: dt,dtratio=dsqrt(3.0d0),coeffx=0.d0,coeffy=0.d0,coeffz=0.d0!,rhomean
  integer :: i=0,n,m,l,countn
  double precision :: rho(-2:ndx+1,-2:ndy+1,-2:ndz+1),Phiwvtest(-2:ndx+1,-2:ndy+1,-2:ndz+1)
  double precision :: Phiwvdffxpyp,Phiwvdffxmyp,Phiwvdffypzp,Phiwvdffymzp,Phiwvdffzpxp,Phiwvdffzmxp, &
                      Phiwvdffxpym,Phiwvdffxmym,Phiwvdffypzm,Phiwvdffymzm,Phiwvdffzpxm,Phiwvdffzmxm
  double precision  grdxy1,grdyz1,grdzx1,grdxy2,grdyz2,grdzx2,grdxy3,grdyz3,grdzx3,grdxy4,grdyz4,grdzx4,&
                    grdxy5,grdyz5,grdzx5,grdxy7,grdyz7,grdzx7,grdxy8,grdyz8,grdzx8,grdxy6,grdyz6,grdzx6

  double precision  grdxy1zp,grdxy1zm,grdxy1mn,grdyz1xp,grdyz1xm,grdyz1mn,grdzx1yp,grdzx1ym,grdzx1mn
  double precision  grdxy2zp,grdxy2zm,grdxy2mn,grdyz2xp,grdyz2xm,grdyz2mn,grdzx2yp,grdzx2ym,grdzx2mn
  double precision  grdxy2zpp,grdxy2zmm,grdyz2xpp,grdyz2xmm,grdzx2ypp,grdzx2ymm
  double precision :: adiff2=0.5d0
  double precision dtt2

  dtt2=dt!*0.3d0

  !rhomean=0.d0
  do l=1,ndz-2
  do m=1,ndy-2
  do n=1,ndx-2
     !rho(n,m,l) = U(n,m,l,1)
     rho(n,m,l) = U(n,m,l,1)!-rhomean
  !   rhomean=rhomean+rho(i,j,k)
  end do;end do;end do
  

  !do l=-1,ndz
  !do m=-1,ndy
  !do n=-1,ndx
  !   rho(n,m,l) = U(n,m,l,1)-rhomean
  !   source(n,m,l,1)=-source(n,m,l,1)+G4pi*rho(n,m,l)
  !   source(n,m,l,2)=-source(n,m,l,2)+G4pi*rho(n,m,l)
  !   source(n,m,l,3)=-source(n,m,l,3)+G4pi*rho(n,m,l)
     !source(n,m,l,1)=source(n,m,l,1)-G4pi*rho(n,m,l)
     !source(n,m,l,2)=source(n,m,l,2)-G4pi*rho(n,m,l)
     !source(n,m,l,3)=source(n,m,l,3)-G4pi*rho(n,m,l)
  !end do;end do;end do
  !****************slv-wv****************
  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   !Phiwvtest(:,:,:)=NRANK
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwvtest,dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   do k=1,ndz-2
     do j=1,ndy-2
        do i=1,ndx-2
!     Phiwv(i,j,k,1) = Phiwv(i,j,k,1)+cg*dt*Phigrdwv(i,j,k,1)*0.5d0
!     Phiwv(i,j,k,2) = Phiwv(i,j,k,2)+cg*dt*Phigrdwv(i,j,k,2)*0.5d0
!     Phiwv(i,j,k,3) = Phiwv(i,j,k,3)+cg*dt*Phigrdwv(i,j,k,3)*0.5d0
!     Phiwv(i,j,k,4) = Phiwv(i,j,k,4)+cg*dt*Phigrdwv(i,j,k,4)*0.5d0
!     Phiwv(i,j,k,5) = Phiwv(i,j,k,5)+cg*dt*Phigrdwv(i,j,k,5)*0.5d0
!     Phiwv(i,j,k,6) = Phiwv(i,j,k,6)+cg*dt*Phigrdwv(i,j,k,6)*0.5d0
!     Phiwv(i,j,k,7) = Phiwv(i,j,k,7)+cg*dt*Phigrdwv(i,j,k,7)*0.5d0
!     Phiwv(i,j,k,8) = Phiwv(i,j,k,8)+cg*dt*Phigrdwv(i,j,k,8)*0.5d0

!Phiwv(i,j,k,1) = Phiwv(i,j,k,1)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,1)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
Phiwv(i,j,k,1) = 0.5d0*Phiwv(i,j,k,1)*(1.d0+dexp(-1.d0/Tdiff * 0.5d0 * dt))+0.5d0*Phigrdwv(i,j,k,1)*(1.d0-dexp(-1.0d0/Tdiff * 0.5d0 * dt))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(-1.d0+dexp(-1.0d0/Tdiff * 0.5d0 * dt))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0*(0.5d0 * dt)
Phiwv(i,j,k,2) = 0.5d0*Phiwv(i,j,k,2)*(1.d0+dexp(-1.d0/Tdiff * 0.5d0 * dtt2))+0.5d0*Phigrdwv(i,j,k,2)*(1.d0-dexp(-1.0d0/Tdiff * 0.5d0 * dtt2))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(-1.d0+dexp(-1.0d0/Tdiff * 0.5d0 * dtt2))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0*(0.5d0 * dtt2)

!Phiwv(i,j,k,2) = Phiwv(i,j,k,2)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,2)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!Phiwv(i,j,k,3) = Phiwv(i,j,k,3)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,3)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!Phiwv(i,j,k,4) = Phiwv(i,j,k,4)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,4)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!Phiwv(i,j,k,5) = Phiwv(i,j,k,5)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,5)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!Phiwv(i,j,k,6) = Phiwv(i,j,k,6)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,6)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!Phiwv(i,j,k,7) = Phiwv(i,j,k,7)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,7)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!Phiwv(i,j,k,8) = Phiwv(i,j,k,8)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,8)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))

!wp2(j,k,2) = wp2(j,k,2)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+wp1(j,k,2)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!wp2(j,k,3) = wp2(j,k,3)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+wp1(j,k,3)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
!wp2(j,k,4) = wp2(j,k,4)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+wp1(j,k,4)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
        enddo
     enddo
   enddo
   !iwx=1;iwy=1;iwz=1
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%


   !%%%%%%%%%%%%%%%%%phigrd(t+0.5*dt)%%%%%%%%%%%%%%%%%%
   iwx=1;iwy=0;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-2,-2,-2,1),dt*0.5d0,2)
   call muslcslv1D(Phigrdwv(-2,-2,-2,2),dtt2*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-2,-2,-2,1),dt*0.5d0,2)
   call muslcslv1D(Phigrdwv(-2,-2,-2,2),dtt2*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,2)
   !call BCgrv(110,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-2,-2,-2,1),dt*0.5d0,2)
   call muslcslv1D(Phigrdwv(-2,-2,-2,2),dtt2*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
  

   !call collectPhi()
      do k=1,ndz-2
        do j=1,ndy-2
          do i=1,ndx-2
Phigrdwv(i,j,k,1) = 0.5d0*Phiwv(i,j,k,1)*(1.d0-dexp(-1.d0/Tdiff *0.5d0* dt))+0.5d0*Phigrdwv(i,j,k,1)*(1.d0+dexp(-1.0d0/Tdiff *0.5d0* dt))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-1.0d0/Tdiff *0.5d0* dt))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0* (0.5d0*dt)
Phigrdwv(i,j,k,2) = 0.5d0*Phiwv(i,j,k,2)*(1.d0-dexp(-1.d0/Tdiff *0.5d0* dtt2))+0.5d0*Phigrdwv(i,j,k,2)*(1.d0+dexp(-1.0d0/Tdiff *0.5d0* dtt2))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-1.0d0/Tdiff *0.5d0* dtt2))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0* (0.5d0*dtt2)

!   Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,1) &
!   -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy1mn/dx1/dy1 &
!   -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz1mn/dy1/dz1 &
!   -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx1mn/dz1/dx1 &
!   -cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

          enddo
        enddo
      enddo

iwx=1;iwy=1;iwz=1
call BCgrv(100,1,8)

   do k=1,ndz-2
     do j=1,ndy-2
       do i=1,ndx-2

grdxy1=adiff*Phiwv(i+1,j+1,k,1)+adiff*Phiwv(i-1,j-1,k,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,1)
grdyz1=adiff*Phiwv(i,j+1,k+1,1)+adiff*Phiwv(i,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,1)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,1)
grdzx1=adiff*Phiwv(i+1,j,k+1,1)+adiff*Phiwv(i-1,j,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,1)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,1)

grdxy1zp=adiff*Phiwv(i+1,j+1,k+1,1)+adiff*Phiwv(i-1,j-1,k+1,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k+1,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k+1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k+1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k+1,1)
grdxy1zm=adiff*Phiwv(i+1,j+1,k-1,1)+adiff*Phiwv(i-1,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k-1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k-1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k-1,1)
grdxy1mn=(grdxy1+grdxy1zp+grdxy1zm)/3.d0

grdyz1xp=adiff*Phiwv(i+1,j+1,k+1,1)+adiff*Phiwv(i+1,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i+1,j+1,k-1,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k+1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i+1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j+1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k+1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i+1,j-1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k-1,1)
grdyz1xm=adiff*Phiwv(i-1,j+1,k+1,1)+adiff*Phiwv(i-1,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j-1,k+1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i-1,j,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j+1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k+1,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i-1,j-1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k-1,1)
grdyz1mn=(grdyz1+grdyz1xp+grdyz1xm)/3.d0

grdzx1yp=adiff*Phiwv(i+1,j+1,k+1,1)+adiff*Phiwv(i-1,j+1,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j+1,k+1,1)+(adiff-0.5d0)*Phiwv(i+1,j+1,k-1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j+1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j+1,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j+1,k,1)
grdzx1ym=adiff*Phiwv(i+1,j-1,k+1,1)+adiff*Phiwv(i-1,j-1,k-1,1)+(adiff-0.5d0)*Phiwv(i-1,j-1,k+1,1)+(adiff-0.5d0)*Phiwv(i+1,j-1,k-1,1) &
+(4.d0*adiff-1.d0)*Phiwv(i,j-1,k,1)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k+1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j-1,k,1)+&
(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k-1,1)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j-1,k,1)
grdzx1mn=(grdzx1+grdzx1yp+grdzx1ym)/3.d0


grdxy2=(-(-Phiwv(i+2,j+2,k,2)+8.d0*Phiwv(i+2,j+1,k,2)-8.d0*Phiwv(i+2,j-1,k,2)+Phiwv(i+2,j-2,k,2))/12.d0 &
   +8.d0*(-Phiwv(i+1,j+2,k,2)+8.d0*Phiwv(i+1,j+1,k,2)-8.d0*Phiwv(i+1,j-1,k,2)+Phiwv(i+1,j-2,k,2))/12.d0 &
   -8.d0*(-Phiwv(i-1,j+2,k,2)+8.d0*Phiwv(i-1,j+1,k,2)-8.d0*Phiwv(i-1,j-1,k,2)+Phiwv(i-1,j-2,k,2))/12.d0 &
        +(-Phiwv(i-2,j+2,k,2)+8.d0*Phiwv(i-2,j+1,k,2)-8.d0*Phiwv(i-2,j-1,k,2)+Phiwv(i-2,j-2,k,2))/12.d0 &
   )/12.d0
grdyz2=(-(-Phiwv(i,j+2,k+2,2)+8.d0*Phiwv(i,j+2,k+1,2)-8.d0*Phiwv(i,j+2,k-1,2)+Phiwv(i,j+2,k-2,2))/12.d0 &
   +8.d0*(-Phiwv(i,j+1,k+2,2)+8.d0*Phiwv(i,j+1,k+1,2)-8.d0*Phiwv(i,j+1,k-1,2)+Phiwv(i,j+1,k-2,2))/12.d0 &
   -8.d0*(-Phiwv(i,j-1,k+2,2)+8.d0*Phiwv(i,j-1,k+1,2)-8.d0*Phiwv(i,j-1,k-1,2)+Phiwv(i,j-1,k-2,2))/12.d0 &
        +(-Phiwv(i,j-2,k+2,2)+8.d0*Phiwv(i,j-2,k+1,2)-8.d0*Phiwv(i,j-2,k-1,2)+Phiwv(i,j-2,k-2,2))/12.d0 &
   )/12.d0
grdzx2=(-(-Phiwv(i+2,j,k+2,2)+8.d0*Phiwv(i+1,j,k+2,2)-8.d0*Phiwv(i-1,j,k+2,2)+Phiwv(i-2,j,k+2,2))/12.d0 &
   +8.d0*(-Phiwv(i+2,j,k+1,2)+8.d0*Phiwv(i+1,j,k+1,2)-8.d0*Phiwv(i-1,j,k+1,2)+Phiwv(i-2,j,k+1,2))/12.d0 &
   -8.d0*(-Phiwv(i+2,j,k-1,2)+8.d0*Phiwv(i+1,j,k-1,2)-8.d0*Phiwv(i-1,j,k-1,2)+Phiwv(i-2,j,k-1,2))/12.d0 &
        +(-Phiwv(i+2,j,k-2,2)+8.d0*Phiwv(i+1,j,k-2,2)-8.d0*Phiwv(i-1,j,k-2,2)+Phiwv(i-2,j,k-2,2))/12.d0 &
   )/12.d0

grdxy2zp =(-(-Phiwv(i+2,j+2,k+1,2)+8.d0*Phiwv(i+2,j+1,k+1,2)-8.d0*Phiwv(i+2,j-1,k+1,2)+Phiwv(i+2,j-2,k+1,2))/12.d0 &
      +8.d0*(-Phiwv(i+1,j+2,k+1,2)+8.d0*Phiwv(i+1,j+1,k+1,2)-8.d0*Phiwv(i+1,j-1,k+1,2)+Phiwv(i+1,j-2,k+1,2))/12.d0 &
      -8.d0*(-Phiwv(i-1,j+2,k+1,2)+8.d0*Phiwv(i-1,j+1,k+1,2)-8.d0*Phiwv(i-1,j-1,k+1,2)+Phiwv(i-1,j-2,k+1,2))/12.d0 &
           +(-Phiwv(i-2,j+2,k+1,2)+8.d0*Phiwv(i-2,j+1,k+1,2)-8.d0*Phiwv(i-2,j-1,k+1,2)+Phiwv(i-2,j-2,k+1,2))/12.d0 &
      )/12.d0
grdxy2zm=(-(-Phiwv(i+2,j+2,k-1,2)+8.d0*Phiwv(i+2,j+1,k-1,2)-8.d0*Phiwv(i+2,j-1,k-1,2)+Phiwv(i+2,j-2,k-1,2))/12.d0 &
     +8.d0*(-Phiwv(i+1,j+2,k-1,2)+8.d0*Phiwv(i+1,j+1,k-1,2)-8.d0*Phiwv(i+1,j-1,k-1,2)+Phiwv(i+1,j-2,k-1,2))/12.d0 &
     -8.d0*(-Phiwv(i-1,j+2,k-1,2)+8.d0*Phiwv(i-1,j+1,k-1,2)-8.d0*Phiwv(i-1,j-1,k-1,2)+Phiwv(i-1,j-2,k-1,2))/12.d0 &
          +(-Phiwv(i-2,j+2,k-1,2)+8.d0*Phiwv(i-2,j+1,k-1,2)-8.d0*Phiwv(i-2,j-1,k-1,2)+Phiwv(i-2,j-2,k-1,2))/12.d0 &
      )/12.d0
grdxy2zpp=(-(-Phiwv(i+2,j+2,k+2,2)+8.d0*Phiwv(i+2,j+1,k+2,2)-8.d0*Phiwv(i+2,j-1,k+2,2)+Phiwv(i+2,j-2,k+2,2))/12.d0 &
      +8.d0*(-Phiwv(i+1,j+2,k+2,2)+8.d0*Phiwv(i+1,j+1,k+2,2)-8.d0*Phiwv(i+1,j-1,k+2,2)+Phiwv(i+1,j-2,k+2,2))/12.d0 &
      -8.d0*(-Phiwv(i-1,j+2,k+2,2)+8.d0*Phiwv(i-1,j+1,k+2,2)-8.d0*Phiwv(i-1,j-1,k+2,2)+Phiwv(i-1,j-2,k+2,2))/12.d0 &
           +(-Phiwv(i-2,j+2,k+2,2)+8.d0*Phiwv(i-2,j+1,k+2,2)-8.d0*Phiwv(i-2,j-1,k+2,2)+Phiwv(i-2,j-2,k+2,2))/12.d0 &
      )/12.d0
grdxy2zmm=(-(-Phiwv(i+2,j+2,k+2,2)+8.d0*Phiwv(i+2,j+1,k+2,2)-8.d0*Phiwv(i+2,j-1,k+2,2)+Phiwv(i+2,j-2,k+2,2))/12.d0 &
      +8.d0*(-Phiwv(i+1,j+2,k+2,2)+8.d0*Phiwv(i+1,j+1,k+2,2)-8.d0*Phiwv(i+1,j-1,k+2,2)+Phiwv(i+1,j-2,k+2,2))/12.d0 &
      -8.d0*(-Phiwv(i-1,j+2,k+2,2)+8.d0*Phiwv(i-1,j+1,k+2,2)-8.d0*Phiwv(i-1,j-1,k+2,2)+Phiwv(i-1,j-2,k+2,2))/12.d0 &
           +(-Phiwv(i-2,j+2,k+2,2)+8.d0*Phiwv(i-2,j+1,k+2,2)-8.d0*Phiwv(i-2,j-1,k+2,2)+Phiwv(i-2,j-2,k+2,2))/12.d0 &
      )/12.d0
grdxy2mn=adiff2*grdxy2-(adiff2-1.d0)*2.d0/3.d0*grdxy2zp-(adiff2-1.d0)*2.d0/3.d0*grdxy2zm &
+(adiff2-1.d0)/2.d0/3.d0*grdxy2zpp+(adiff2-1.d0)/2.d0/3.d0*grdxy2zmm

grdyz2xp=(-(-Phiwv(i+1,j+2,k+2,2)+8.d0*Phiwv(i+1,j+2,k+1,2)-8.d0*Phiwv(i+1,j+2,k-1,2)+Phiwv(i+1,j+2,k-2,2))/12.d0 &
     +8.d0*(-Phiwv(i+1,j+1,k+2,2)+8.d0*Phiwv(i+1,j+1,k+1,2)-8.d0*Phiwv(i+1,j+1,k-1,2)+Phiwv(i+1,j+1,k-2,2))/12.d0 &
     -8.d0*(-Phiwv(i+1,j-1,k+2,2)+8.d0*Phiwv(i+1,j-1,k+1,2)-8.d0*Phiwv(i+1,j-1,k-1,2)+Phiwv(i+1,j-1,k-2,2))/12.d0 &
          +(-Phiwv(i+1,j-2,k+2,2)+8.d0*Phiwv(i+1,j-2,k+1,2)-8.d0*Phiwv(i+1,j-2,k-1,2)+Phiwv(i+1,j-2,k-2,2))/12.d0 &
     )/12.d0
grdyz2xm=(-(-Phiwv(i-1,j+2,k+2,2)+8.d0*Phiwv(i-1,j+2,k+1,2)-8.d0*Phiwv(i-1,j+2,k-1,2)+Phiwv(i-1,j+2,k-2,2))/12.d0 &
     +8.d0*(-Phiwv(i-1,j+1,k+2,2)+8.d0*Phiwv(i-1,j+1,k+1,2)-8.d0*Phiwv(i-1,j+1,k-1,2)+Phiwv(i-1,j+1,k-2,2))/12.d0 &
     -8.d0*(-Phiwv(i-1,j-1,k+2,2)+8.d0*Phiwv(i-1,j-1,k+1,2)-8.d0*Phiwv(i-1,j-1,k-1,2)+Phiwv(i-1,j-1,k-2,2))/12.d0 &
          +(-Phiwv(i-1,j-2,k+2,2)+8.d0*Phiwv(i-1,j-2,k+1,2)-8.d0*Phiwv(i-1,j-2,k-1,2)+Phiwv(i-1,j-2,k-2,2))/12.d0 &
     )/12.d0
grdyz2xpp=(-(-Phiwv(i+2,j+2,k+2,2)+8.d0*Phiwv(i+2,j+2,k+1,2)-8.d0*Phiwv(i+2,j+2,k-1,2)+Phiwv(i+2,j+2,k-2,2))/12.d0 &
      +8.d0*(-Phiwv(i+2,j+1,k+2,2)+8.d0*Phiwv(i+2,j+1,k+1,2)-8.d0*Phiwv(i+2,j+1,k-1,2)+Phiwv(i+2,j+1,k-2,2))/12.d0 &
      -8.d0*(-Phiwv(i+2,j-1,k+2,2)+8.d0*Phiwv(i+2,j-1,k+1,2)-8.d0*Phiwv(i+2,j-1,k-1,2)+Phiwv(i+2,j-1,k-2,2))/12.d0 &
           +(-Phiwv(i+2,j-2,k+2,2)+8.d0*Phiwv(i+2,j-2,k+1,2)-8.d0*Phiwv(i+2,j-2,k-1,2)+Phiwv(i+2,j-2,k-2,2))/12.d0 &
      )/12.d0
grdyz2xmm=(-(-Phiwv(i-2,j+2,k+2,2)+8.d0*Phiwv(i-2,j+2,k+1,2)-8.d0*Phiwv(i-2,j+2,k-1,2)+Phiwv(i-2,j+2,k-2,2))/12.d0 &
      +8.d0*(-Phiwv(i-2,j+1,k+2,2)+8.d0*Phiwv(i-2,j+1,k+1,2)-8.d0*Phiwv(i-2,j+1,k-1,2)+Phiwv(i-2,j+1,k-2,2))/12.d0 &
      -8.d0*(-Phiwv(i-2,j-1,k+2,2)+8.d0*Phiwv(i-2,j-1,k+1,2)-8.d0*Phiwv(i-2,j-1,k-1,2)+Phiwv(i-2,j-1,k-2,2))/12.d0 &
           +(-Phiwv(i-2,j-2,k+2,2)+8.d0*Phiwv(i-2,j-2,k+1,2)-8.d0*Phiwv(i-2,j-2,k-1,2)+Phiwv(i-2,j-2,k-2,2))/12.d0 &
      )/12.d0
grdyz2mn=adiff2*grdyz2-(adiff2-1.d0)*2.d0/3.d0*grdyz2xp-(adiff2-1.d0)*2.d0/3.d0*grdyz2xm &
+(adiff2-1.d0)/2.d0/3.d0*grdyz2xpp+(adiff2-1.d0)/2.d0/3.d0*grdyz2xmm

grdzx2yp=(-(-Phiwv(i+2,j+1,k+2,2)+8.d0*Phiwv(i+1,j+1,k+2,2)-8.d0*Phiwv(i-1,j+1,k+2,2)+Phiwv(i-2,j+1,k+2,2))/12.d0 &
     +8.d0*(-Phiwv(i+2,j+1,k+1,2)+8.d0*Phiwv(i+1,j+1,k+1,2)-8.d0*Phiwv(i-1,j+1,k+1,2)+Phiwv(i-2,j+1,k+1,2))/12.d0 &
     -8.d0*(-Phiwv(i+2,j+1,k-1,2)+8.d0*Phiwv(i+1,j+1,k-1,2)-8.d0*Phiwv(i-1,j+1,k-1,2)+Phiwv(i-2,j+1,k-1,2))/12.d0 &
          +(-Phiwv(i+2,j+1,k-2,2)+8.d0*Phiwv(i+1,j+1,k-2,2)-8.d0*Phiwv(i-1,j+1,k-2,2)+Phiwv(i-2,j+1,k-2,2))/12.d0 &
     )/12.d0
grdzx2ym=(-(-Phiwv(i+2,j-1,k+2,2)+8.d0*Phiwv(i+1,j-1,k+2,2)-8.d0*Phiwv(i-1,j-1,k+2,2)+Phiwv(i-2,j-1,k+2,2))/12.d0 &
     +8.d0*(-Phiwv(i+2,j-1,k+1,2)+8.d0*Phiwv(i+1,j-1,k+1,2)-8.d0*Phiwv(i-1,j-1,k+1,2)+Phiwv(i-2,j-1,k+1,2))/12.d0 &
     -8.d0*(-Phiwv(i+2,j-1,k-1,2)+8.d0*Phiwv(i+1,j-1,k-1,2)-8.d0*Phiwv(i-1,j-1,k-1,2)+Phiwv(i-2,j-1,k-1,2))/12.d0 &
          +(-Phiwv(i+2,j-1,k-2,2)+8.d0*Phiwv(i+1,j-1,k-2,2)-8.d0*Phiwv(i-1,j-1,k-2,2)+Phiwv(i-2,j-1,k-2,2))/12.d0 &
     )/12.d0
grdzx2ypp=(-(-Phiwv(i+2,j+2,k+2,2)+8.d0*Phiwv(i+1,j+2,k+2,2)-8.d0*Phiwv(i-1,j+2,k+2,2)+Phiwv(i-2,j+2,k+2,2))/12.d0 &
      +8.d0*(-Phiwv(i+2,j+2,k+1,2)+8.d0*Phiwv(i+1,j+2,k+1,2)-8.d0*Phiwv(i-1,j+2,k+1,2)+Phiwv(i-2,j+2,k+1,2))/12.d0 &
      -8.d0*(-Phiwv(i+2,j+2,k-1,2)+8.d0*Phiwv(i+1,j+2,k-1,2)-8.d0*Phiwv(i-1,j+2,k-1,2)+Phiwv(i-2,j+2,k-1,2))/12.d0 &
           +(-Phiwv(i+2,j+2,k-2,2)+8.d0*Phiwv(i+1,j+2,k-2,2)-8.d0*Phiwv(i-1,j+2,k-2,2)+Phiwv(i-2,j+2,k-2,2))/12.d0 &
      )/12.d0
grdzx2ymm=(-(-Phiwv(i+2,j-2,k+2,2)+8.d0*Phiwv(i+1,j-2,k+2,2)-8.d0*Phiwv(i-1,j-2,k+2,2)+Phiwv(i-2,j-2,k+2,2))/12.d0 &
      +8.d0*(-Phiwv(i+2,j-2,k+1,2)+8.d0*Phiwv(i+1,j-2,k+1,2)-8.d0*Phiwv(i-1,j-2,k+1,2)+Phiwv(i-2,j-2,k+1,2))/12.d0 &
      -8.d0*(-Phiwv(i+2,j-2,k-1,2)+8.d0*Phiwv(i+1,j-2,k-1,2)-8.d0*Phiwv(i-1,j-2,k-1,2)+Phiwv(i-2,j-2,k-1,2))/12.d0 &
           +(-Phiwv(i+2,j-2,k-2,2)+8.d0*Phiwv(i+1,j-2,k-2,2)-8.d0*Phiwv(i-1,j-2,k-2,2)+Phiwv(i-2,j-2,k-2,2))/12.d0 &
      )/12.d0
grdzx2mn=adiff2*grdzx2-(adiff2-1.d0)*2.d0/3.d0*grdzx2yp-(adiff2-1.d0)*2.d0/3.d0*grdzx2ym &
+(adiff2-1.d0)/2.d0/3.d0*grdzx2ypp+(adiff2-1.d0)/2.d0/3.d0*grdzx2ymm

!Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,1) &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy1/dx1/dy1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz1/dy1/dz1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx1/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
(-2.d0*2.d0*Tdiff*cg*cg*grdxy1mn/dx1/dy1 &
 -2.d0*2.d0*Tdiff*cg*cg*grdyz1mn/dy1/dz1 &
 -2.d0*2.d0*Tdiff*cg*cg*grdzx1mn/dz1/dx1) *dt

Phigrdwv(i,j,k,2) = Phigrdwv(i,j,k,2) +&
(-2.d0*2.d0*Tdiff*cg*cg*grdxy2mn/dx1/dy1 &
 +2.d0*2.d0*Tdiff*cg*cg*grdyz2mn/dy1/dz1 &
 +2.d0*2.d0*Tdiff*cg*cg*grdzx2mn/dz1/dx1) *dtt2

!Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,1) &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy1mn/dx1/dy1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz1mn/dy1/dz1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx1mn/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

       enddo
     enddo
   enddo

      do k=1,ndz-2
        do j=1,ndy-2
          do i=1,ndx-2
Phigrdwv(i,j,k,1) = 0.5d0*Phiwv(i,j,k,1)*(1.d0-dexp(-1.d0/Tdiff *0.5d0* dt))+0.5d0*Phigrdwv(i,j,k,1)*(1.d0+dexp(-1.0d0/Tdiff *0.5d0* dt))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-1.0d0/Tdiff *0.5d0* dt))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0* (0.5d0*dt)
Phigrdwv(i,j,k,2) = 0.5d0*Phiwv(i,j,k,2)*(1.d0-dexp(-1.d0/Tdiff *0.5d0* dtt2))+0.5d0*Phigrdwv(i,j,k,2)*(1.d0+dexp(-1.0d0/Tdiff *0.5d0* dtt2))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-1.0d0/Tdiff *0.5d0* dtt2))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0* (0.5d0*dtt2)

!   Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,1) &
!   -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy1mn/dx1/dy1 &
!   -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz1mn/dy1/dz1 &
!   -2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx1mn/dz1/dx1 &
!   -cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

          enddo
        enddo
      enddo
   !iwx=1;iwy=1;iwz=1
   !call BCgrv(110,1,8)

   iwx=0;iwy=0;iwz=1
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-2,-2,-2,1),dt*0.5d0,2)
   call muslcslv1D(Phigrdwv(-2,-2,-2,2),dtt2*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-2,-2,-2,1),dt*0.5d0,2)
   call muslcslv1D(Phigrdwv(-2,-2,-2,2),dtt2*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,2)
   !call BCgrv(110,1,8)
   iwx=1;iwy=0;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-2,-2,-2,1),dt*0.5d0,2)
   call muslcslv1D(Phigrdwv(-2,-2,-2,2),dtt2*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
   !%%%%%%%%%%%%%%%%%phigrd(t+0.5*dt)%%%%%%%%%%%%%%%%%%


  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
      do k=1,ndz-2
        do j=1,ndy-2
           do i=1,ndx-2
   !     Phiwv(i,j,k,1) = Phiwv(i,j,k,1)+cg*dt*Phigrdwv(i,j,k,1)*0.5d0
   !     Phiwv(i,j,k,2) = Phiwv(i,j,k,2)+cg*dt*Phigrdwv(i,j,k,2)*0.5d0
   !     Phiwv(i,j,k,3) = Phiwv(i,j,k,3)+cg*dt*Phigrdwv(i,j,k,3)*0.5d0
   !     Phiwv(i,j,k,4) = Phiwv(i,j,k,4)+cg*dt*Phigrdwv(i,j,k,4)*0.5d0
   !     Phiwv(i,j,k,5) = Phiwv(i,j,k,5)+cg*dt*Phigrdwv(i,j,k,5)*0.5d0
   !     Phiwv(i,j,k,6) = Phiwv(i,j,k,6)+cg*dt*Phigrdwv(i,j,k,6)*0.5d0
   !     Phiwv(i,j,k,7) = Phiwv(i,j,k,7)+cg*dt*Phigrdwv(i,j,k,7)*0.5d0
   !     Phiwv(i,j,k,8) = Phiwv(i,j,k,8)+cg*dt*Phigrdwv(i,j,k,8)*0.5d0
Phiwv(i,j,k,1) = 0.5d0*Phiwv(i,j,k,1)*(1.d0+dexp(-1.d0/Tdiff * 0.5d0 * dt))+0.5d0*Phigrdwv(i,j,k,1)*(1.d0-dexp(-1.0d0/Tdiff * 0.5d0 * dt))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(-1.d0+dexp(-1.0d0/Tdiff * 0.5d0 * dt))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0*(0.5d0 * dt)
Phiwv(i,j,k,2) = 0.5d0*Phiwv(i,j,k,2)*(1.d0+dexp(-1.d0/Tdiff * 0.5d0 * dtt2))+0.5d0*Phigrdwv(i,j,k,2)*(1.d0-dexp(-1.0d0/Tdiff * 0.5d0 * dtt2))+&
(-cg*cg*Tdiff*Tdiff*G4pi*rho(i,j,k))*(-1.d0+dexp(-1.0d0/Tdiff * 0.5d0 * dtt2))-cg*cg*2.d0*Tdiff*G4pi*rho(i,j,k)*0.5d0*(0.5d0 * dtt2)

   !Phiwv(i,j,k,1) = Phiwv(i,j,k,1)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,1)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,2) = Phiwv(i,j,k,2)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,2)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,3) = Phiwv(i,j,k,3)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,3)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,4) = Phiwv(i,j,k,4)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,4)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,5) = Phiwv(i,j,k,5)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,5)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,6) = Phiwv(i,j,k,6)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,6)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,7) = Phiwv(i,j,k,7)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,7)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !Phiwv(i,j,k,8) = Phiwv(i,j,k,8)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,8)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))

   !wp2(j,k,2) = wp2(j,k,2)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+wp1(j,k,2)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !wp2(j,k,3) = wp2(j,k,3)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+wp1(j,k,3)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
   !wp2(j,k,4) = wp2(j,k,4)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+wp1(j,k,4)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
           enddo
        enddo
      enddo
   !iwx=1;iwy=1;iwz=1
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-2,-2,-2,1),dt*0.25d0,1)
   call muslcslv1D(Phiwv(-2,-2,-2,2),dtt2*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  !****************slv-wv****************

  !iwx=1; iwy=1; iwz=1
  !call BCgrv(102,1)
  !call BCgrv(102,2)
  !call BCgrv(102,3)
end subroutine slvmuscle


subroutine BCgrv(mode,is,ie)
  use comvar
  use mpivar
  use slfgrv
  INCLUDE 'mpif.h'
  integer ::  N_ol=3,i,mode,j,k,idm,is,ie
  INTEGER :: MSTATUS(MPI_STATUS_SIZE)
  DOUBLE PRECISION  :: VECU
 
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!***************BC-for-Phigrd***********************
if(mode==100) then
  IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+3+1)*(Ncellz+6),N_ol,ndx+3+1,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT!; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT!; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
  do idm=is,ie
  CALL MPI_SENDRECV(Phiwv(Ncellx+1-N_ol,-2,-2,idm),1,VECU,RIGT,1, &
  Phiwv(       1-N_ol,-2,-2,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !IF((IST.eq.0).and.(mod(idm,2)==1)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
  !   Phiwv(IX,JY,KZ,idm)= bphil(JY,KZ,IX)
     !Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
  !   END DO;END DO;END DO
  !END IF
  !IF((IST.eq.0).and.(mod(idm,2)==0)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
  !   Phiwv(IX,JY,KZ,idm)=Phiwv(IX+1,JY,KZ,idm)
  !   END DO;END DO;END DO
  !END IF
!  IF(IST.eq.0) THEN
!     DO KZ = -2, Ncellz+3; DO JY = -2, Ncelly+3; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     !Phiwv(IX,JY,KZ,idm)= bphil(JY,KZ,IX)
!     Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
!     END DO;END DO;END DO
!  END IF
  enddo
  do idm=is,ie
  CALL MPI_SENDRECV(Phiwv(1            ,-2,-2,idm),1,VECU,LEFT,1, &
  Phiwv(Ncellx+1     ,-2,-2,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  !IF((IST.eq.NSPLTx-1).and.(mod(idm,2)==0)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
  !   Phiwv(IX,JY,KZ,idm)= bphir(JY,KZ,IX)
     !Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
  !   END DO;END DO;END DO
  !END IF
  !IF((IST.eq.NSPLTx-1).and.(mod(idm,2)==1)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX =Ncellx+N_ol, Ncellx+1, -1
  !   Phiwv(IX,JY,KZ,idm)=Phiwv(IX-1,JY,KZ,idm)
  !   END DO;END DO;END DO
  !END IF
!  IF(IST.eq.NSPLTx-1) THEN
!     DO KZ = -2, Ncellz+3; DO JY = -2, Ncelly+3; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     !Phiwv(IX,JY,KZ,idm)= bphir(JY,KZ,IX)
!     Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
!     END DO;END DO;END DO
!  END IF
  enddo
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF

IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+6,N_ol*(ndx+3+1),(ndx+3+1)*(ndy+3+1),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phiwv(-2,Ncelly+1-N_ol,-2,idm),1,VECU,TOP ,1, &
        Phiwv(-2,       1-N_ol,-2,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  ! IF(JST.eq.0) THEN
  !      DO KZ = -2, Ncellz+3; DO JY =  1-N_ol, 0; DO IX = -2, Ncellx+3
  !      Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
  !      END DO;END DO;END DO
  !   END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phiwv(-2,1            ,-2,idm),1,VECU,BOTM,1, &
        Phiwv(-2,Ncelly+1     ,-2,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  ! IF(JST.eq.NSPLTy-1) THEN
  !      DO KZ = -2, Ncellz+3; DO JY = Ncelly+1, Ncelly+N_ol; DO IX = -2, Ncellx+3
  !      Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
  !      END DO;END DO;END DO
  !   END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'BCinside-dy',iwx,iwy,iwz
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+3+1)*(ndy+3+1),N_ol*(ndx+3+1)*(ndy+3+1),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
  CALL MPI_SENDRECV(Phiwv(-2,-2,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
  Phiwv(-2,-2,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
 !  IF(KST.eq.0) THEN
 !       DO KZ = 1-N_ol, 0; DO JY = -2, Ncelly+3; DO IX = -2, Ncellx+3
 !       Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
 !       END DO;END DO;END DO
 !    END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phiwv(-2,-2,1            ,idm),1,VECU,DOWN,1, &
   Phiwv(-2,-2,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
 ! IF(KST.eq.NSPLTz-1) THEN
 !    DO KZ =  Ncellz+1, Ncellz+N_ol; DO JY = -2, Ncelly+3; DO IX = -2, Ncellx+3
 !    Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
 !    END DO;END DO;END DO
 ! END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt
END IF

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'BCinside-dz',iwx,iwy,iwz
endif
!***************BC-for-Phiwv***********************


!***************BC-for-Phiwvgrd***********************
if(mode==110) then
  IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+3+1)*(Ncellz+6),N_ol,ndx+3+1,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT!; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT!; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL

  do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(Ncellx+1-N_ol,-2,-2,idm),1,VECU,RIGT,1, &
  Phigrdwv(       1-N_ol,-2,-2,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)


  !IF((IST.eq.0).and.(mod(idm,2)==0)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
  !   Phigrdwv(IX,JY,KZ,idm)= bphigrdxl(JY,KZ,IX,idm)
     !Phigrdwv(IX,JY,KZ,idm) = Phigrd(IX,JY,KZ,idm)
  !   END DO;END DO;END DO
  !END IF
  !IF((IST.eq.0).and.(mod(idm,2)==1)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
  !   Phigrdwv(IX,JY,KZ,idm)= Phigrdwv(IX+1,JY,KZ,idm)
  !   END DO;END DO;END DO
  !END IF

 ! IF(IST.eq.0) THEN
 !    DO KZ = -2, Ncellz+3; DO JY = -2, Ncelly+3; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     !Phigrdwv(IX,JY,KZ,idm)= bphigrdxl(JY,KZ,IX,idm)
 !    Phigrdwv(IX,JY,KZ,idm) = cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
 !    END DO;END DO;END DO
 ! END IF
  enddo
  do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(1            ,-2,-2,idm),1,VECU,LEFT,1, &
  Phigrdwv(Ncellx+1     ,-2,-2,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  !IF((IST.eq.NSPLTx-1).and.(mod(idm,2)==1)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
  !   Phigrdwv(IX,JY,KZ,idm)= bphigrdxr(JY,KZ,IX,idm)
     !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
  !   END DO;END DO;END DO
  !END IF
  !IF((IST.eq.NSPLTx-1).and.(mod(idm,2)==0)) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX =Ncellx+N_ol, Ncellx+1, -1
  !   Phigrdwv(IX,JY,KZ,idm)= Phigrdwv(IX-1,JY,KZ,idm)
  !   END DO;END DO;END DO
  !END IF
 ! IF(IST.eq.NSPLTx-1) THEN
 !    DO KZ = -2, Ncellz+3; DO JY = -2, Ncelly+3; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     !Phigrdwv(IX,JY,KZ,idm)= bphigrdxr(JY,KZ,IX,idm)
 !    Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
 !    END DO;END DO;END DO
 ! END IF
  enddo

CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF


IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+6,N_ol*(ndx+3+1),(ndx+3+1)*(ndy+3+1),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-2,Ncelly+1-N_ol,-2,idm),1,VECU,TOP ,1, &
        Phigrdwv(-2,       1-N_ol,-2,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
 !  IF(JST.eq.0) THEN
 !       DO KZ = -2, Ncellz+3; DO JY = 1-N_ol, 0; DO IX = -2, Ncellx+3
 !       Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
 !       END DO;END DO;END DO
 !    END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-2,1            ,-2,idm),1,VECU,BOTM,1, &
        Phigrdwv(-2,Ncelly+1     ,-2,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
 !  IF(JST.eq.NSPLTy-1) THEN
 !       DO KZ = -2, Ncellz+3; DO JY = Ncelly+1, Ncelly+N_ol; DO IX =  -2, Ncellx+3
        !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
 !       Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
 !       END DO;END DO;END DO
 !    END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+3+1)*(ndy+3+1),N_ol*(ndx+3+1)*(ndy+3+1),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(-2,-2,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
  Phigrdwv(-2,-2,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
!   IF(KST.eq.0) THEN
!        DO KZ = 1-N_ol, 0; DO JY = -2, Ncelly+3; DO IX = -2, Ncellx+3
        !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
!        Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
!        END DO;END DO;END DO
!     END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-2,-2,1            ,idm),1,VECU,DOWN,1, &
   Phigrdwv(-2,-2,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  IF(KST.eq.NSPLTz-1) THEN
!       DO KZ =  Ncellz+1, Ncellz+N_ol; DO JY = -2, Ncelly+3; DO IX =  -2, Ncellx+3
       !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
!       Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
!       END DO;END DO;END DO
!    END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt
END IF
endif
!***************BC-for-Phiwvgrd***********************
end subroutine BCgrv


subroutine muslcslv1D(Phiv,dt,mode)
  use comvar
  !INCLUDE 'mpif.h'
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam ,deltalen !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  DOUBLE PRECISION, dimension(-2:ndx+1,-2:ndy+1,-2:ndz+1) :: Phigrad,Phipre,Phiv,Phi2dt,Phiu
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie,idm
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'INSIDE00',dt,mode,Phiv(-2,-2,-2),Phiv(ndx+1,ndy+1,ndz+1)
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) 'INSIDE0'


if(iwx.eq.1) then; Ncell = ndx+1; Ncm = ndy+1; Ncl = ndz+1; deltalen=dx1; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy+1; Ncm = ndz+1; Ncl = ndx+1; deltalen=dy1; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz+1; Ncm = ndx+1; Ncl = ndy+1; deltalen=dz1; endif! BT1 = 1; BT2 = 2; VN = 4; end if

  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !write(*,*) 'INSIDE0'
  !----kyoukai-----
   !if(hazi==1)then
   !   is = 2
   !   ie = Ncell-3
   !end if
   !if(hazi==2)then
      is = 1
      ie = Ncell-3
   !end if
  !----kyoukai-----
  nu2 = cg * dt / deltalen
  Phipre(:,:,:) = Phiv(:,:,:)
  !------------ul.solver.+cg-------------

  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     DO Lnum = 1-3, Ncl!-2+3
        DO Mnum = 1-3, Ncm!-2+3
           do i = is-2,ie+3
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
     !do i=ist-1,ndx-ien+1 !一次なので大丈夫
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz)- 0.5d0 * nu2 * ( Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
              !Phi2dt(ix,jy) = Phipre(ix,jy)- 0.5d0 * nu2 * ( Phiu(ix,jy) - Phiu(ixm,jym))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
      DO Lnum = 1-3, Ncl!-2+3
        DO Mnum = 1-3, Ncm!-2+3
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) - nu2 * (Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
           end do
        end DO
     end DO
  end if
  !------------ul.solver.+cg-------------



  !------------ul.solver.-cg-------------
  if(mode==2) then

     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,11,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,11)
     !------------calcurate dt/2------------
     DO Lnum = 1-3, Ncl!-2+3
        DO Mnum = 1-3, Ncm!-2+3
           do i = is-3,ie+2
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i=ist-1,ndx-ien+1
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz) + 0.5d0 * nu2 * ( Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
              !Phi2dt(ix,jy) = Phipre(ix,jy) + 0.5d0 * nu2 * ( Phiu(ixp,jyp) - Phiu(ix,jy))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4,is,ie)
     

     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
     DO Lnum = 1-3, Ncl!-2+3
        DO Mnum = 1-3, Ncm!-2+3
           do i = is,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) + nu2 * (Phiu(ix,jy,kz) - Phiu(ixm,jym,kzm))
              !Phiv(ix,jy) = Phipre(ix,jy) + nu2 * (Phiu(ix,jy) - Phiu(ixm,jym))
           end do
        end DO
     end DO

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do

  end if
  !------------ul.solver.-cg-------------

  cnt=cnt+2
end subroutine muslcslv1D


subroutine fluxcal(preuse,pre,uin,ep,kappa,mode,is,ie)
  use comvar
  !INCLUDE 'mpif.h'
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-2:ndx+1,-2:ndy+1,-2:ndz+1) :: ul,ur,pre,preuse,uin
  DOUBLE PRECISION , dimension(-2:ndx+1) :: slop  !------------- need allocation --------------
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie,ixpp,jypp,kzpp
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  DOUBLE PRECISION dqp12,dqm12,dqp32,ph,minmdp12,sgnp12,minmdm12,sgnm12,minmdp32,sgnp32 &
  ,d3qp12bar,d3qm12bar,d3qp12til,d3qp32til,minmd,sgn,ql,qr
    DOUBLE PRECISION , dimension(-2:ndx+1,-2:ndy+1,-2:ndz+1) :: d3qp12
  !uin(:)=0.0d0
  if(iwx.eq.1) then; Ncell = ndx+1; Ncm = ndy+1; Ncl = ndz+1;  end if
     if(iwy.eq.1) then; Ncell = ndy+1; Ncm = ndz+1; Ncl = ndx+1;  end if
        if(iwz.eq.1) then; Ncell = ndz+1; Ncm = ndx+1; Ncl = ndy+1;  end if



           !call vanalbada(pre,slop)
           if(mode==1) then
              DO Lnum = 1-3, Ncl!-2+3
              DO Mnum = 1-3, Ncm!-2+3
              !call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-2,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixpp= iwx*(i+2)+ iwy*Lnum + iwz*Mnum
              jypp= iwx*Mnum + iwy*(i+2)+ iwz*Lnum
              kzpp= iwx*Lnum + iwy*Mnum + iwz*(i+2)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              sgnm12=dsign(1.d0,pre(ix,jy,kz)-pre(ixm,jym,kzm))
              minmdm12=sgnm12*dmax1(0.d0,dmin1(dabs(pre(ix,jy,kz)-pre(ixm,jym,kzm)),sgnm12*2.d0*(pre(ixp,jyp,kzp)-pre(ix,jy,kz))&
              ,sgnm12*2.d0*(pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp))))

              sgnp12=dsign(1.d0,pre(ixp,jyp,kzp)-pre(ix,jy,kz))
              minmdp12=sgnp12*dmax1(0.d0,dmin1(dabs(pre(ixp,jyp,kzp)-pre(ix,jy,kz)),sgnp12*2.d0*(pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp))&
              ,sgnp12*2.d0*(pre(ix,jy,kz)-pre(ixm,jym,kzm))))

              sgnp32=dsign(1.d0,pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp))
              minmdp32=sgnp32*dmax1(0.d0,dmin1(dabs(pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp)),sgnp32*2.d0*(pre(ix,jy,kz)-pre(ixm,jym,kzm))&
              ,sgnp32*2.d0*(pre(ixp,jyp,kzp)-pre(ix,jy,kz))))

              d3qp12(ix,jy,kz)=pre(ixp,jyp,kzp)-pre(ix,jy,kz)-(minmdm12-2.d0*minmdp12+minmdp32)/6.d0

              !ul(ix,jy,kz) = preuse(ix,jy,kz) + 0.25d0 * ep * slop(i) &
              !     * ((1.0d0-slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
              !     (1.0d0+slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i+1/2
              !uin(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
              DO Lnum = 1-3, Ncl!-2+3
              DO Mnum = 1-3, Ncm!-2+3
              !call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixpp= iwx*(i+2)+ iwy*Lnum + iwz*Mnum
              jypp= iwx*Mnum + iwy*(i+2)+ iwz*Lnum
              kzpp= iwx*Lnum + iwy*Mnum + iwz*(i+2)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              sgn=dsign(1.d0,d3qp12(ixm,jym,kzm))
              d3qm12bar=sgn*dmax1(0.d0,dmin1(dabs(d3qp12(ixm,jym,kzm)),sgn*4.d0*(d3qp12(ix,jy,kz))))
              !sgn=dsign(1.d0,d3qp12(ix,jy))
              !d3qp12bar=sgn*dmax1(0,dmin1(dabs(d3qp12(ix,jy)),sgn*4.d0*(d3qp12(ixp,jyp))))

              sgn=dsign(1.d0,d3qp12(ix,jy,kz))
              d3qp12til=sgn*dmax1(0.d0,dmin1(dabs(d3qp12(ix,jy,kz)),sgn*4.d0*(d3qp12(ixm,jym,kzm))))

              !ul(ix,jy,kz) = preuse(ix,jy,kz) + 0.25d0 * ep * slop(i) &
              !     * ((1.0d0-slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
              !     (1.0d0+slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i+1/2
              !uin(ix,jy,kz)=ul(ix,jy,kz)
              uin(ix,jy,kz)=preuse(ix,jy,kz)+d3qm12bar/6.d0+d3qp12til/3.d0
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !uin(:)=ul(:)
           end if


           if(mode==4) then
              DO Lnum = 1-3, Ncl!-2+3
              DO Mnum = 1-3, Ncm!-2+3
              !call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-2,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixpp= iwx*(i+2)+ iwy*Lnum + iwz*Mnum
              jypp= iwx*Mnum + iwy*(i+2)+ iwz*Lnum
              kzpp= iwx*Lnum + iwy*Mnum + iwz*(i+2)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

              sgnm12=dsign(1.d0,pre(ix,jy,kz)-pre(ixm,jym,kzm))
              minmdm12=sgnm12*dmax1(0.d0,dmin1(dabs(pre(ix,jy,kz)-pre(ixm,jym,kzm)),sgnm12*2.d0*(pre(ixp,jyp,kzp)-pre(ix,jy,kz))&
              ,sgnm12*2.d0*(pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp))))

              sgnp12=dsign(1.d0,pre(ixp,jyp,kzp)-pre(ix,jy,kz))
              minmdp12=sgnp12*dmax1(0.d0,dmin1(dabs(pre(ixp,jyp,kzp)-pre(ix,jy,kz)),sgnp12*2.d0*(pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp))&
              ,sgnp12*2.d0*(pre(ix,jy,kz)-pre(ixm,jym,kzm))))

              sgnp32=dsign(1.d0,pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp))
              minmdp32=sgnp32*dmax1(0.d0,dmin1(dabs(pre(ixpp,jypp,kzpp)-pre(ixp,jyp,kzp)),sgnp32*2.d0*(pre(ix,jy,kz)-pre(ixm,jym,kzm))&
              ,sgnp32*2.d0*(pre(ixp,jyp,kzp)-pre(ix,jy,kz))))

              d3qp12(ix,jy,kz)=pre(ixp,jyp,kzp)-pre(ix,jy,kz)-(minmdm12-2.d0*minmdp12+minmdp32)/6.d0
              !do i = ist-1,ndx-ien+1
              !ur(ix,jy,kz) = preuse(ix,jy,kz) - 0.25d0 * ep * slop(i) &
              !     * ((1.0d0+slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
              !     (1.0d0-slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i-1/2
              !uin(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO

            
              DO Lnum = 1-3, Ncl!-2+3
              DO Mnum = 1-3, Ncm!-2+3
              !call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-2,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixpp= iwx*(i+2)+ iwy*Lnum + iwz*Mnum
              jypp= iwx*Mnum + iwy*(i+2)+ iwz*Lnum
              kzpp= iwx*Lnum + iwy*Mnum + iwz*(i+2)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

              sgn=dsign(1.d0,d3qp12(ix,jy,kz))
              d3qp12bar=sgn*dmax1(0.d0,dmin1(dabs(d3qp12(ix,jy,kz)),sgn*4.d0*(d3qp12(ixp,jyp,kzp))))
              !sgn=dsign(1.d0,d3qp12(ixm,jym,kzm))
              !d3qm12bar=sgn*dmax1(0.d0,dmin1(dabs(d3qp12(ixm,jym,kzm)),sgn*4.d0*(d3qp12(ix,jy,kz))))

              sgn=dsign(1.d0,d3qp12(ixp,jyp,kzp))
              d3qp32til=sgn*dmax1(0.d0,dmin1(dabs(d3qp12(ixp,jyp,kzp)),sgn*4.d0*(d3qp12(ix,jy,kz))))
              !sgn=dsign(1.d0,d3qp12(ix,jy,kz))
              !d3qp12til=sgn*dmax1(0.d0,dmin1(dabs(d3qp12(ix,jy,kz)),sgn*4.d0*(d3qp12(ixm,jym,kzm))))

              uin(ix,jy,kz)=preuse(ixp,jyp,kzp)-d3qp32til/6.d0-d3qp12bar/3.d0
              !uin(ix,jy,kz)=preuse(ix,jy,kz)+d3qm12bar/6.d0+d3qp12til/3.d0

              !do i = ist-1,ndx-ien+1
              !ur(ix,jy,kz) = preuse(ix,jy,kz) - 0.25d0 * ep * slop(i) &
              !     * ((1.0d0+slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
              !     (1.0d0-slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i-1/2
              !uin(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !uin(:)=ur(:)
           end if

           if(mode==10) then
              DO Lnum = 1-3, Ncl!-2+3
              DO Mnum = 1-3, Ncm!-2+3
              !call vanalbada(pre,slop)
              do i = is-3,ie+3
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-2,ndx-ien+2
              ul(ix,jy,kz) = preuse(ix,jy,kz)
              uin(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
           end if

           if(mode==11) then
              DO Lnum = 1-3, Ncl!-2+3
              DO Mnum = 1-3, Ncm!-2+3
              !call vanalbada(pre,slop)
              do i = is-3,ie+3
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = is,ie
              ur(ix,jy,kz) = preuse(ix,jy,kz)
              uin(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
           end if


end subroutine fluxcal


SUBROUTINE collect()
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: tMPI(1:Ncellx,1:Ncelly,1:Ncellz,0:NPE-1)

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

rhomean=0.d0
do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
  tMPI(i,j,k,NRANK)=U(i,j,k,1)
end do;end do;end do
do Nroot=0,NPE-1
  CALL MPI_BCAST(tMPI(1,1,1,Nroot),(Ncellx)*(Ncelly)*(Ncellz),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
do Nroot=0,NPE-1
 ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
do kk=1,Ncellz!; k=KSTt*Ncellz+kk
do jj=1,Ncelly!; j=JSTt*Ncelly+jj
do ii=1,Ncellz!; i=ISTt*Ncellx+ii
    !u1(i,j,k) = tMPI(ii,jj,kk,Nroot)
    rhomean = tMPI(ii,jj,kk,Nroot)+rhomean
end do;end do;end do;end do
rhomean=rhomean/dble(Ncellx*NSPLTx)/dble(Ncelly*NSPLTy)/dble(Ncellz*NSPLTz)

!rhomean=0.d0
END SUBROUTINE collect

SUBROUTINE collectrho()
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: meanrho(0:NPE-1)

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

meanrho(:)=0.d0
rhomean=0.d0
do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
  meanrho(NRANK)=U(i,j,k,1)+meanrho(NRANK)
end do;end do;end do
meanrho(NRANK)=meanrho(NRANK)/(dble(Ncellx*Ncelly*Ncellz))
do Nroot=0,NPE-1
  CALL MPI_BCAST(meanrho(Nroot),1,MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
do Nroot=0,NPE-1
    rhomean = meanrho(Nroot)+rhomean
end do
rhomean=rhomean/dble(NPE)

!rhomean=0.d0
END SUBROUTINE collectrho

SUBROUTINE collectPhi()
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: tMPI(1:Ncellx,1:Ncelly,1:Ncellz,0:NPE-1,wvnum),u1(NSPLTx*Ncellx,NSPLTy*Ncelly,NSPLTz*Ncellz,wvnum)
integer :: num1

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

u1(:,:,:,:)=0.d0

do num1=1,wvnum
do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
  tMPI(i,j,k,NRANK,num1)=Phiwv(i,j,k,num1)
end do;end do;end do
do Nroot=0,NPE-1
  CALL MPI_BCAST(tMPI(1,1,1,Nroot,num1),(Ncellx)*(Ncelly)*(Ncellz),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do

do Nroot=0,NPE-1
 ISTt = mod(Nroot,NSPLTx); KSTt = Nroot/(NSPLTx*NSPLTy); JSTt = Nroot/NSPLTx-NSPLTy*KSTt
do kk=1,Ncellz; k=KSTt*Ncellz+kk
do jj=1,Ncelly; j=JSTt*Ncelly+jj
do ii=1,Ncellz; i=ISTt*Ncellx+ii
    u1(i,j,k,num1) = tMPI(ii,jj,kk,Nroot,num1)
end do;end do;end do;end do


do kk=-2,Ncellz+3; k=KST*Ncellz+kk
do jj=-2,Ncelly+3; j=JST*Ncelly+jj
do ii=-2,Ncellz+3; i=IST*Ncellx+ii
    if((j.eq.NSPLTy*Ncelly+2).and.(JST.eq.NSPLTy-1)) j  = 2
    if((k.eq.NSPLTz*Ncellz+2).and.(KST.eq.NSPLTz-1)) k = 2
    if((j.eq.NSPLTy*Ncelly+3).and.(JST.eq.NSPLTy-1)) j  = 3
    if((k.eq.NSPLTz*Ncellz+3).and.(KST.eq.NSPLTz-1)) k = 3
    if((j.eq.NSPLTy*Ncelly+1).and.(JST.eq.NSPLTy-1)) j  = 1
    if((k.eq.NSPLTy*Ncelly+1).and.(KST.eq.NSPLTz-1)) k = 1
    if((j.eq.0  ).and.(JST.eq.0       )) j  = Ncelly*NSPLTy
    if((k.eq.0  ).and.(KST.eq.0       )) k = Ncellz*NSPLTz
    if((j.eq.-1  ).and.(JST.eq.0       )) j  = Ncelly*NSPLTy-1
    if((k.eq.-1  ).and.(KST.eq.0       )) k = Ncellz*NSPLTz-1
    if((j.eq.-2  ).and.(JST.eq.0       )) j  = Ncelly*NSPLTy-2
    if((k.eq.-2  ).and.(KST.eq.0       )) k = Ncellz*NSPLTz-2
    if((i.eq.NSPLTx*Ncellx+2).and.(IST.eq.NSPLTx-1)) i  = 2
    if((i.eq.NSPLTx*Ncellx+3).and.(IST.eq.NSPLTx-1)) i  = 3
    if((i.eq.NSPLTx*Ncellx+1).and.(IST.eq.NSPLTx-1)) i  = 1
    if((i.eq.0  ).and.(JST.eq.0       )) i  = Ncellx*NSPLTx
    if((i.eq.-1  ).and.(JST.eq.0       )) i  = Ncellx*NSPLTx-1
    if((i.eq.-2  ).and.(JST.eq.0       )) i  = Ncellx*NSPLTx-2

    Phiwv(ii,jj,kk,num1)=u1(i,j,k,num1)
end do;end do;end do


IF(IST.eq.0) THEN
   DO KZ = -3, Ncellz+3; DO JY = -3, Ncelly+3; DO IX = 1-3, 0
   Phiwv(IX,JY,KZ,num1)= bphil(JY,KZ,IX)
   END DO;END DO;END DO
END IF
IF(IST.eq.NSPLTx-1) THEN
   DO KZ = -3, Ncellz+3; DO JY = -3, Ncelly+3; DO IX = Ncellx+1, Ncellx+3
   Phiwv(IX,JY,KZ,num1)= bphir(JY,KZ,IX)
   END DO;END DO;END DO
END IF
enddo

END SUBROUTINE collectPhi
