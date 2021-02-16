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
     open(unit=28,file=dir//'PHIINI/PHI'//countcha//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           do i = -1, Ncellx+2
           write(28) sngl(Phiwv(i,j,k,1)),&!sngl(Phiwv(i,j,k,2)),sngl(Phiwv(i,j,k,3)),sngl(Phiwv(i,j,k,4)),&!sngl(Phiwv(i,j,k,5)),sngl(Phiwv(i,j,k,6)),&
                   !sngl(Phiwv(i,j,k,7)),sngl(Phiwv(i,j,k,8)),
                 sngl(Phigrdwv(i,j,k,1)),&!,sngl(Phigrdwv(i,j,k,2)),sngl(Phigrdwv(i,j,k,3)),sngl(Phigrdwv(i,j,k,4)),&
               !sngl(Phigrdwv(i,j,k,5)),sngl(Phigrdwv(i,j,k,6)),sngl(Phigrdwv(i,j,k,7)),sngl(Phigrdwv(i,j,k,8))
                      sngl(Phiexa(i,j,k)),sngl(Phigrd(i,j,k,1)),sngl(U(i,j,k,1))!,&
               !sngl(Phigrd(i,j,k,2)),sngl(Phigrd(i,j,k,3)),sngl(Phigrd(i,j,k,4)),sngl(Phigrd(i,j,k,5)),sngl(Phigrd(i,j,k,6)),sngl(Phigrd(i,j,k,7)),sngl(Phigrd(i,j,k,8)),&
               !sngl(Phiexab1(i,j,k)),sngl(Phiexab2(i,j,k))
          enddo
        end do
     end do
     close(28)
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
   !sourcedt = dt
   !cgtime = deltalength/cg * CFL
   !call STBLphi(sourcedt,Phipregrad)
   !dt = sourcedt
end if
!***************SABILITY-exa**************
end subroutine SELFGRAVWAVE


subroutine slvmuscle(dt)
  use comvar
  use slfgrv
  INCLUDE 'mpif.h'
  double precision :: dt,dtratio=dsqrt(3.0d0),coeffx=0.d0,coeffy=0.d0,coeffz=0.d0!,rhomean
  integer :: i=0,n,m,l,countn
  double precision :: rho(-1:ndx,-1:ndy,-1:ndz)
  double precision :: Phiwvdffxpyp,Phiwvdffxmyp,Phiwvdffypzp,Phiwvdffymzp,Phiwvdffzpxp,Phiwvdffzmxp, &
                      Phiwvdffxpym,Phiwvdffxmym,Phiwvdffypzm,Phiwvdffymzm,Phiwvdffzpxm,Phiwvdffzmxm
  double precision  grdxy1,grdyz1,grdzx1,grdxy2,grdyz2,grdzx2,grdxy3,grdyz3,grdzx3,grdxy4,grdyz4,grdzx4,&
                    grdxy5,grdyz5,grdzx5,grdxy7,grdyz7,grdzx7,grdxy8,grdyz8,grdzx8,grdxy6,grdyz6,grdzx6

  double precision  grdxy1zp,grdxy1zm,grdxy1mn,grdyz1xp,grdyz1xm,grdyz1mn,grdzx1yp,grdzx1ym,grdzx1mn

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
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
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

Phiwv(i,j,k,1) = Phiwv(i,j,k,1)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,1)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
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
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%


   !%%%%%%%%%%%%%%%%%phigrd(t+0.5*dt)%%%%%%%%%%%%%%%%%%
   iwx=1;iwy=0;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,4),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,7),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,3),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,8),dt*0.5d0,2)
   !call BCgrv(110,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,5),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
  

   !call collectPhi()
   iwx=1;iwy=1;iwz=1
   call BCgrv(100,1,8)
   do k=1,ndz-2
     do j=1,ndy-2
       do i=1,ndx-2
!       Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)-cg*G4pi*rho(i,j,k)*dt &
!       -0.5d0*dt*cg*(Phiwv(i+1,j+1,k,1)-Phiwv(i+1,j-1,k,1)-Phiwv(i-1,j+1,k,1)+Phiwv(i-1,j-1,k,1))/dx1/dy1 &
!       -0.5d0*dt*cg*(Phiwv(i,j+1,k+1,1)-Phiwv(i,j+1,k-1,1)-Phiwv(i,j-1,k+1,1)+Phiwv(i,j-1,k-1,1))/dy1/dz1 &
!       -0.5d0*dt*cg*(Phiwv(i+1,j,k+1,1)-Phiwv(i-1,j,k+1,1)-Phiwv(i+1,j,k-1,1)+Phiwv(i-1,j,k-1,1))/dz1/dx1
!       Phigrdwv(i,j,k,2) = Phigrdwv(i,j,k,2)-cg*G4pi*rho(i,j,k)*dt &
!       -0.5d0*dt*cg*(Phiwv(i+1,j+1,k,2)-Phiwv(i+1,j-1,k,2)-Phiwv(i-1,j+1,k,2)+Phiwv(i-1,j-1,k,2))/dx1/dy1 &
!       +0.5d0*dt*cg*(Phiwv(i,j+1,k+1,2)-Phiwv(i,j+1,k-1,2)-Phiwv(i,j-1,k+1,2)+Phiwv(i,j-1,k-1,2))/dy1/dz1 &
!       +0.5d0*dt*cg*(Phiwv(i+1,j,k+1,2)-Phiwv(i-1,j,k+1,2)-Phiwv(i+1,j,k-1,2)+Phiwv(i-1,j,k-1,2))/dz1/dx1
!       Phigrdwv(i,j,k,3) = Phigrdwv(i,j,k,3)-cg*G4pi*rho(i,j,k)*dt &
!       +0.5d0*dt*cg*(Phiwv(i+1,j+1,k,3)-Phiwv(i+1,j-1,k,3)-Phiwv(i-1,j+1,k,3)+Phiwv(i-1,j-1,k,3))/dx1/dy1 &
!       +0.5d0*dt*cg*(Phiwv(i,j+1,k+1,3)-Phiwv(i,j+1,k-1,3)-Phiwv(i,j-1,k+1,3)+Phiwv(i,j-1,k-1,3))/dy1/dz1 &
!       -0.5d0*dt*cg*(Phiwv(i+1,j,k+1,3)-Phiwv(i-1,j,k+1,3)-Phiwv(i+1,j,k-1,3)+Phiwv(i-1,j,k-1,3))/dz1/dx1
!       Phigrdwv(i,j,k,4) = Phigrdwv(i,j,k,4)-cg*G4pi*rho(i,j,k)*dt &
!       +0.5d0*dt*cg*(Phiwv(i+1,j+1,k,4)-Phiwv(i+1,j-1,k,4)-Phiwv(i-1,j+1,k,4)+Phiwv(i-1,j-1,k,4))/dx1/dy1 &
!       -0.5d0*dt*cg*(Phiwv(i,j+1,k+1,4)-Phiwv(i,j+1,k-1,4)-Phiwv(i,j-1,k+1,4)+Phiwv(i,j-1,k-1,4))/dy1/dz1 &
!       +0.5d0*dt*cg*(Phiwv(i+1,j,k+1,4)-Phiwv(i-1,j,k+1,4)-Phiwv(i+1,j,k-1,4)+Phiwv(i-1,j,k-1,4))/dz1/dx1
!       Phigrdwv(i,j,k,5) = Phigrdwv(i,j,k,5)-cg*G4pi*rho(i,j,k)*dt &
!       -0.5d0*dt*cg*(Phiwv(i+1,j+1,k,5)-Phiwv(i+1,j-1,k,5)-Phiwv(i-1,j+1,k,5)+Phiwv(i-1,j-1,k,5))/dx1/dy1 &
!       +0.5d0*dt*cg*(Phiwv(i,j+1,k+1,5)-Phiwv(i,j+1,k-1,5)-Phiwv(i,j-1,k+1,5)+Phiwv(i,j-1,k-1,5))/dy1/dz1 &
!       +0.5d0*dt*cg*(Phiwv(i+1,j,k+1,5)-Phiwv(i-1,j,k+1,5)-Phiwv(i+1,j,k-1,5)+Phiwv(i-1,j,k-1,5))/dz1/dx1
!       Phigrdwv(i,j,k,6) = Phigrdwv(i,j,k,6)-cg*G4pi*rho(i,j,k)*dt &
!       -0.5d0*dt*cg*(Phiwv(i+1,j+1,k,6)-Phiwv(i+1,j-1,k,6)-Phiwv(i-1,j+1,k,6)+Phiwv(i-1,j-1,k,6))/dx1/dy1 &
!       -0.5d0*dt*cg*(Phiwv(i,j+1,k+1,6)-Phiwv(i,j+1,k-1,6)-Phiwv(i,j-1,k+1,6)+Phiwv(i,j-1,k-1,6))/dy1/dz1 &
!       -0.5d0*dt*cg*(Phiwv(i+1,j,k+1,6)-Phiwv(i-1,j,k+1,6)-Phiwv(i+1,j,k-1,6)+Phiwv(i-1,j,k-1,6))/dz1/dx1
!       Phigrdwv(i,j,k,7) = Phigrdwv(i,j,k,7)-cg*G4pi*rho(i,j,k)*dt &
!       +0.5d0*dt*cg*(Phiwv(i+1,j+1,k,7)-Phiwv(i+1,j-1,k,7)-Phiwv(i-1,j+1,k,7)+Phiwv(i-1,j-1,k,7))/dx1/dy1 &
!       -0.5d0*dt*cg*(Phiwv(i,j+1,k+1,7)-Phiwv(i,j+1,k-1,7)-Phiwv(i,j-1,k+1,7)+Phiwv(i,j-1,k-1,7))/dy1/dz1 &
!       +0.5d0*dt*cg*(Phiwv(i+1,j,k+1,7)-Phiwv(i-1,j,k+1,7)-Phiwv(i+1,j,k-1,7)+Phiwv(i-1,j,k-1,7))/dz1/dx1
!       Phigrdwv(i,j,k,8) = Phigrdwv(i,j,k,8)-cg*G4pi*rho(i,j,k)*dt &
!       +0.5d0*dt*cg*(Phiwv(i+1,j+1,k,8)-Phiwv(i+1,j-1,k,8)-Phiwv(i-1,j+1,k,8)+Phiwv(i-1,j-1,k,8))/dx1/dy1 &
!       +0.5d0*dt*cg*(Phiwv(i,j+1,k+1,8)-Phiwv(i,j+1,k-1,8)-Phiwv(i,j-1,k+1,8)+Phiwv(i,j-1,k-1,8))/dy1/dz1 &
!       -0.5d0*dt*cg*(Phiwv(i+1,j,k+1,8)-Phiwv(i-1,j,k+1,8)-Phiwv(i+1,j,k-1,8)+Phiwv(i-1,j,k-1,8))/dz1/dx1

!grdxy1= +(Phiwv(i+1,j+1,k,1)/3.d0+Phiwv(i+1,j,k,1)/2.d0-Phiwv(i+1,j-1,k,1)+Phiwv(i+1,j-2,k,1)/6.d0)/3.d0 &
!        +(Phiwv(i  ,j+1,k,1)/3.d0+Phiwv(i  ,j,k,1)/2.d0-Phiwv(i  ,j-1,k,1)+Phiwv(i  ,j-2,k,1)/6.d0)/2.d0 &
!        -(Phiwv(i-1,j+1,k,1)/3.d0+Phiwv(i-1,j,k,1)/2.d0-Phiwv(i-1,j-1,k,1)+Phiwv(i-1,j-2,k,1)/6.d0) &
!        +(Phiwv(i-2,j+1,k,1)/3.d0+Phiwv(i-2,j,k,1)/2.d0-Phiwv(i-2,j-1,k,1)+Phiwv(i-2,j-2,k,1)/6.d0)/6.d0

!grdyz1= +(Phiwv(i,j+1,k+1,1)/3.d0+Phiwv(i,j+1,k,1)/2.d0-Phiwv(i,j+1,k-1,1)+Phiwv(i,j+1,k-2,1)/6.d0)/3.d0 &
!        +(Phiwv(i,j  ,k+1,1)/3.d0+Phiwv(i,j  ,k,1)/2.d0-Phiwv(i,j  ,k-1,1)+Phiwv(i,j  ,k-2,1)/6.d0)/2.d0 &
!        -(Phiwv(i,j-1,k+1,1)/3.d0+Phiwv(i,j-1,k,1)/2.d0-Phiwv(i,j-1,k-1,1)+Phiwv(i,j-1,k-2,1)/6.d0) &
!        +(Phiwv(i,j-2,k+1,1)/3.d0+Phiwv(i,j-2,k,1)/2.d0-Phiwv(i,j-2,k-1,1)+Phiwv(i,j-2,k-2,1)/6.d0)/6.d0

!grdzx1= +(Phiwv(i+1,j,k+1,1)/3.d0+Phiwv(i,j,k+1,1)/2.d0-Phiwv(i-1,j,k+1,1)+Phiwv(i-2,j,k+1,1)/6.d0)/3.d0 &
!        +(Phiwv(i+1,j,k  ,1)/3.d0+Phiwv(i,j,k  ,1)/2.d0-Phiwv(i-1,j,k  ,1)+Phiwv(i-2,j,k  ,1)/6.d0)/2.d0 &
!        -(Phiwv(i+1,j,k-1,1)/3.d0+Phiwv(i,j,k-1,1)/2.d0-Phiwv(i-1,j,k-1,1)+Phiwv(i-2,j,k-1,1)/6.d0) &
!        +(Phiwv(i+1,j,k-2,1)/3.d0+Phiwv(i,j,k-2,1)/2.d0-Phiwv(i-1,j,k-2,1)+Phiwv(i-2,j,k-2,1)/6.d0)/6.d0

!grdxy1=(-(-Phiwv(i+2,j+2,k,1)+8.d0*Phiwv(i+2,j+1,k,1)-8.d0*Phiwv(i+2,j-1,k,1)+Phiwv(i+2,j-2,k,1))/12.d0 &
!   +8.d0*(-Phiwv(i+1,j+2,k,1)+8.d0*Phiwv(i+1,j+1,k,1)-8.d0*Phiwv(i+1,j-1,k,1)+Phiwv(i+1,j-2,k,1))/12.d0 &
!   -8.d0*(-Phiwv(i-1,j+2,k,1)+8.d0*Phiwv(i-1,j+1,k,1)-8.d0*Phiwv(i-1,j-1,k,1)+Phiwv(i-1,j-2,k,1))/12.d0 &
!        +(-Phiwv(i-2,j+2,k,1)+8.d0*Phiwv(i-2,j+1,k,1)-8.d0*Phiwv(i-2,j-1,k,1)+Phiwv(i-2,j-2,k,1))/12.d0 &
!        )/12.d0
!grdyz1=(-(-Phiwv(i,j+2,k+2,1)+8.d0*Phiwv(i,j+2,k+1,1)-8.d0*Phiwv(i,j+2,k-1,1)+Phiwv(i,j+2,k-2,1))/12.d0 &
!   +8.d0*(-Phiwv(i,j+1,k+2,1)+8.d0*Phiwv(i,j+1,k+1,1)-8.d0*Phiwv(i,j+1,k-1,1)+Phiwv(i,j+1,k-2,1))/12.d0 &
!   -8.d0*(-Phiwv(i,j-1,k+2,1)+8.d0*Phiwv(i,j-1,k+1,1)-8.d0*Phiwv(i,j-1,k-1,1)+Phiwv(i,j-1,k-2,1))/12.d0 &
!        +(-Phiwv(i,j-2,k+2,1)+8.d0*Phiwv(i,j-2,k+1,1)-8.d0*Phiwv(i,j-2,k-1,1)+Phiwv(i,j-2,k-2,1))/12.d0 &
!        )/12.d0
!grdzx1=(-(-Phiwv(i+2,j,k+2,1)+8.d0*Phiwv(i+1,j,k+2,1)-8.d0*Phiwv(i-1,j,k+2,1)+Phiwv(i-2,j,k+2,1))/12.d0 &
!   +8.d0*(-Phiwv(i+2,j,k+1,1)+8.d0*Phiwv(i+1,j,k+1,1)-8.d0*Phiwv(i-1,j,k+1,1)+Phiwv(i-2,j,k+1,1))/12.d0 &
!   -8.d0*(-Phiwv(i+2,j,k-1,1)+8.d0*Phiwv(i+1,j,k-1,1)-8.d0*Phiwv(i-1,j,k-1,1)+Phiwv(i-2,j,k-1,1))/12.d0 &
!        +(-Phiwv(i+2,j,k-2,1)+8.d0*Phiwv(i+1,j,k-2,1)-8.d0*Phiwv(i-1,j,k-2,1)+Phiwv(i-2,j,k-2,1))/12.d0 &
!       )/12.d0

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

!grdxy2=adiff*Phiwv(i+1,j+1,k,2)+adiff*Phiwv(i-1,j-1,k,2)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,2)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,2) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,2)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,2)
!grdyz2=adiff*Phiwv(i,j+1,k+1,2)+adiff*Phiwv(i,j-1,k-1,2)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,2)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,2) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,2)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,2)
!grdzx2=adiff*Phiwv(i+1,j,k+1,2)+adiff*Phiwv(i-1,j,k-1,2)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,2)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,2) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,2)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,2)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,2)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,2)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,2)

!grdxy3=adiff*Phiwv(i+1,j+1,k,3)+adiff*Phiwv(i-1,j-1,k,3)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,3)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,3) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,3)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,3)
!grdyz3=adiff*Phiwv(i,j+1,k+1,3)+adiff*Phiwv(i,j-1,k-1,3)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,3)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,3) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,3)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,3)
!grdzx3=adiff*Phiwv(i+1,j,k+1,3)+adiff*Phiwv(i-1,j,k-1,3)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,3)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,3) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,3)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,3)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,3)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,3)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,3)

!grdxy4=adiff*Phiwv(i+1,j+1,k,4)+adiff*Phiwv(i-1,j-1,k,4)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,4)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,4) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,4)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,4)
!grdyz4=adiff*Phiwv(i,j+1,k+1,4)+adiff*Phiwv(i,j-1,k-1,4)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,4)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,4) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,4)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,4)
!grdzx4=adiff*Phiwv(i+1,j,k+1,4)+adiff*Phiwv(i-1,j,k-1,4)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,4)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,4) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,4)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,4)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,4)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,4)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,4)

!grdxy5=adiff*Phiwv(i+1,j+1,k,5)+adiff*Phiwv(i-1,j-1,k,5)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,5)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,5) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,5)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,5)
!grdyz5=adiff*Phiwv(i,j+1,k+1,5)+adiff*Phiwv(i,j-1,k-1,5)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,5)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,5) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,5)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,5)
!grdzx5=adiff*Phiwv(i+1,j,k+1,5)+adiff*Phiwv(i-1,j,k-1,5)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,5)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,5) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,5)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,5)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,5)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,5)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,5)

!grdxy6=adiff*Phiwv(i+1,j+1,k,6)+adiff*Phiwv(i-1,j-1,k,6)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,6)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,6) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,6)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,6)
!grdyz6=adiff*Phiwv(i,j+1,k+1,6)+adiff*Phiwv(i,j-1,k-1,6)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,6)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,6) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,6)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,6)
!grdzx6=adiff*Phiwv(i+1,j,k+1,6)+adiff*Phiwv(i-1,j,k-1,6)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,6)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,6) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,6)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,6)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,6)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,6)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,6)

!grdxy7=adiff*Phiwv(i+1,j+1,k,7)+adiff*Phiwv(i-1,j-1,k,7)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,7)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,7) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,7)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,7)
!grdyz7=adiff*Phiwv(i,j+1,k+1,7)+adiff*Phiwv(i,j-1,k-1,7)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,7)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,7) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,7)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,7)
!grdzx7=adiff*Phiwv(i+1,j,k+1,7)+adiff*Phiwv(i-1,j,k-1,7)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,7)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,7) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,7)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,7)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,7)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,7)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,7)

!grdxy8=adiff*Phiwv(i+1,j+1,k,8)+adiff*Phiwv(i-1,j-1,k,8)+(adiff-0.5d0)*Phiwv(i+1,j-1,k,8)+(adiff-0.5d0)*Phiwv(i-1,j+1,k,8) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,8)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,8)
!grdyz8=adiff*Phiwv(i,j+1,k+1,8)+adiff*Phiwv(i,j-1,k-1,8)+(adiff-0.5d0)*Phiwv(i,j+1,k-1,8)+(adiff-0.5d0)*Phiwv(i,j-1,k+1,8) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i,j+1,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,8)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j-1,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,8)
!grdzx8=adiff*Phiwv(i+1,j,k+1,8)+adiff*Phiwv(i-1,j,k-1,8)+(adiff-0.5d0)*Phiwv(i-1,j,k+1,8)+(adiff-0.5d0)*Phiwv(i+1,j,k-1,8) &
!+(4.d0*adiff-1.d0)*Phiwv(i,j,k,8)+(-2.d0*adiff+0.5d0)*Phiwv(i,j,k+1,8)+(-2.d0*adiff+0.5d0)*Phiwv(i+1,j,k,8)+&
!(-2.d0*adiff+0.5d0)*Phiwv(i,j,k-1,8)+(-2.d0*adiff+0.5d0)*Phiwv(i-1,j,k,8)
!grdch(i,j)=grdxy

!grdxy1=(-(-Phiwv(i+2,j+2,k,1)+8.d0*Phiwv(i+2,j+1,k,1)-8.d0*Phiwv(i+2,j-1,k,1)+Phiwv(i+2,j-2,k,1))/12.d0 &
!+8.d0*(-Phiwv(i+1,j+2,k,1)+8.d0*Phiwv(i+1,j+1,k,1)-8.d0*Phiwv(i+1,j-1,k,1)+Phiwv(i+1,j-2,k,1))/12.d0 &
!-8.d0*(-Phiwv(i-1,j+2,k,1)+8.d0*Phiwv(i-1,j+1,k,1)-8.d0*Phiwv(i-1,j-1,k,1)+Phiwv(i-1,j-2,k,1))/12.d0 &
!+(-Phiwv(i-2,j+2,k,1)+8.d0*Phiwv(i-2,j+1,k,1)-8.d0*Phiwv(i-2,j-1,k,1)+Phiwv(i-2,j-2,k,1))/12.d0 &
!)/12.d0
!grdyz1=(-(-Phiwv(i,j+2,k+2,1)+8.d0*Phiwv(i,j+2,k+1,1)-8.d0*Phiwv(i,j+2,k-1,1)+Phiwv(i,j+2,k-2,1))/12.d0 &
!+8.d0*(-Phiwv(i,j+1,k+2,1)+8.d0*Phiwv(i,j+1,k+1,1)-8.d0*Phiwv(i,j+1,k-1,1)+Phiwv(i,j+1,k-2,1))/12.d0 &
!-8.d0*(-Phiwv(i,j-1,k+2,1)+8.d0*Phiwv(i,j-1,k+1,1)-8.d0*Phiwv(i,j-1,k-1,1)+Phiwv(i,j-1,k-2,1))/12.d0 &
!+(-Phiwv(i,j-2,k+2,1)+8.d0*Phiwv(i,j-2,k+1,1)-8.d0*Phiwv(i,j-2,k-1,1)+Phiwv(i,j-2,k-2,1))/12.d0 &
!)/12.d0
!grdzx1=(-(-Phiwv(i+2,j,k+2,1)+8.d0*Phiwv(i+1,j,k+2,1)-8.d0*Phiwv(i-1,j,k+2,1)+Phiwv(i-2,j,k+2,1))/12.d0 &
!+8.d0*(-Phiwv(i+2,j,k+1,1)+8.d0*Phiwv(i+1,j,k+1,1)-8.d0*Phiwv(i-1,j,k+1,1)+Phiwv(i-2,j,k+1,1))/12.d0 &
!-8.d0*(-Phiwv(i+2,j,k-1,1)+8.d0*Phiwv(i+1,j,k-1,1)-8.d0*Phiwv(i-1,j,k-1,1)+Phiwv(i-2,j,k-1,1))/12.d0 &
!+(-Phiwv(i+2,j,k-2,1)+8.d0*Phiwv(i+1,j,k-2,1)-8.d0*Phiwv(i-1,j,k-2,1)+Phiwv(i-2,j,k-2,1))/12.d0 &
!)/12.d0

!Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,1) &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy1/dx1/dy1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz1/dy1/dz1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx1/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,1) &
-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy1mn/dx1/dy1 &
-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz1mn/dy1/dz1 &
-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx1mn/dz1/dx1 &
-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,2) = Phigrdwv(i,j,k,2)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,2) &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy2/dx1/dy1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz2/dy1/dz1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx2/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,3) = Phigrdwv(i,j,k,3)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,3) &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy3/dx1/dy1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz3/dy1/dz1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx3/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,4) = Phigrdwv(i,j,k,4)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,4) &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy4/dx1/dy1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz4/dy1/dz1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx4/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,5) = Phigrdwv(i,j,k,5)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,5) &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy5/dx1/dy1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz5/dy1/dz1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx5/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,6) = Phigrdwv(i,j,k,6)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,6) &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy6/dx1/dy1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz6/dy1/dz1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx6/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,7) = Phigrdwv(i,j,k,7)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,7) &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy7/dx1/dy1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz7/dy1/dz1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx7/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))

!Phigrdwv(i,j,k,8) = Phigrdwv(i,j,k,8)*dexp(-0.5d0/Tdiff * dt)+(Phiwv(i,j,k,8) &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdxy8/dx1/dy1 &
!+2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdyz8/dy1/dz1 &
!-2.d0*2.d0*Tdiff*2.d0*Tdiff*cg*cg*grdzx8/dz1/dx1 &
!-cg*cg*4.d0*Tdiff*Tdiff*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0/Tdiff * dt))
       enddo
     enddo
   enddo
   !iwx=1;iwy=1;iwz=1
   !call BCgrv(110,1,8)

   iwx=0;iwy=0;iwz=1
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,5),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,3),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,4),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,7),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,8),dt*0.5d0,2)
   !call BCgrv(110,1,8)
   iwx=1;iwy=0;iwz=0
   call BCgrv(110,1,8)
   call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,3),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,4),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,5),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,6),dt*0.5d0,1)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,7),dt*0.5d0,2)
   !call muslcslv1D(Phigrdwv(-1,-1,-1,8),dt*0.5d0,1)
   !call BCgrv(110,1,8)
   !%%%%%%%%%%%%%%%%%phigrd(t+0.5*dt)%%%%%%%%%%%%%%%%%%


  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=0;iwy=0;iwz=1
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
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

   Phiwv(i,j,k,1) = Phiwv(i,j,k,1)*dexp(-0.5d0/Tdiff * 0.5d0 * dt)+Phigrdwv(i,j,k,1)*(1.d0-dexp(-0.5d0/Tdiff * 0.5d0 * dt))
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
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
   !call BCgrv(100,1,8)
   iwx=0;iwy=1;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,1)
   !call BCgrv(100,1,8)
   iwx=1;iwy=0;iwz=0
   call BCgrv(100,1,8)
   call muslcslv1D(Phiwv(-1,-1,-1,1),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,2),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,3),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,4),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,5),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,6),dt*0.25d0,2)
   !call muslcslv1D(Phiwv(-1,-1,-1,7),dt*0.25d0,1)
   !call muslcslv1D(Phiwv(-1,-1,-1,8),dt*0.25d0,2)
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
  integer ::  N_ol=2,i,mode,j,k,idm,is,ie
  INTEGER :: MSTATUS(MPI_STATUS_SIZE)
  DOUBLE PRECISION  :: VECU
 
  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!***************BC-for-Phigrd***********************
if(mode==100) then
  IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
  do idm=is,ie
  CALL MPI_SENDRECV(Phiwv(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
  Phiwv(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
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
  IF(IST.eq.0) THEN
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     !Phiwv(IX,JY,KZ,idm)= bphil(JY,KZ,IX)
     Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
     END DO;END DO;END DO
  END IF
  enddo
  do idm=is,ie
  CALL MPI_SENDRECV(Phiwv(1            ,-1,-1,idm),1,VECU,LEFT,1, &
  Phiwv(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

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

  !IF(IST.eq.NSPLTx-1) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     !Phiwv(IX,JY,KZ,idm)= bphir(JY,KZ,IX)
  !   Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
  !   END DO;END DO;END DO
  !END IF
  enddo
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF

IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM ; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  ; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phiwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
        Phiwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
   IF(JST.eq.0) THEN
        DO KZ = -1, Ncellz+2; DO JY =  1-N_ol, 0; DO IX = -1, Ncellx+2
        Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
        END DO;END DO;END DO
     END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phiwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
        Phiwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
   !IF(JST.eq.NSPLTy-1) THEN
   !     DO KZ = -1, Ncellz+2; DO JY = Ncelly+1, Ncelly+N_ol; DO IX = -1, Ncellx+2
   !     Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
   !     END DO;END DO;END DO
   !  END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN ; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   ; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
  CALL MPI_SENDRECV(Phiwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
  Phiwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
   IF(KST.eq.0) THEN
        DO KZ = 1-N_ol, 0; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
        Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
        END DO;END DO;END DO
     END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phiwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
   Phiwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !IF(KST.eq.NSPLTz-1) THEN
  !   DO KZ =  Ncellz+1, Ncellz+N_ol; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
  !   Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
  !   END DO;END DO;END DO
  !END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt
END IF
endif
!***************BC-for-Phiwv***********************


!***************BC-for-Phiwvgrd***********************
if(mode==110) then
  IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL

  do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
  Phigrdwv(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)


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

  !IF(IST.eq.0) THEN
  !   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     !Phigrdwv(IX,JY,KZ,idm)= bphigrdxl(JY,KZ,IX,idm)
  !   Phigrdwv(IX,JY,KZ,idm) = cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
  !   END DO;END DO;END DO
  !END IF
  enddo
  do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(1            ,-1,-1,idm),1,VECU,LEFT,1, &
  Phigrdwv(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

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
  IF(IST.eq.NSPLTx-1) THEN
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     !Phigrdwv(IX,JY,KZ,idm)= bphigrdxr(JY,KZ,IX,idm)
     Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
     END DO;END DO;END DO
  END IF
  enddo

CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF


IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM ; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  ; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
        Phigrdwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
   !IF(JST.eq.0) THEN
   !     DO KZ = -1, Ncellz+2; DO JY = 1-N_ol, 0; DO IX = -1, Ncellx+2
   !     Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
   !     END DO;END DO;END DO
   !  END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
        Phigrdwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
   IF(JST.eq.NSPLTy-1) THEN
        DO KZ = -1, Ncellz+2; DO JY = Ncelly+1, Ncelly+N_ol; DO IX =  -1, Ncellx+2
        !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
        Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
        END DO;END DO;END DO
     END IF
   enddo
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN ; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   ; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
  Phigrdwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
   !IF(KST.eq.0) THEN
   !     DO KZ = 1-N_ol, 0; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
   !     !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
   !     Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
   !     END DO;END DO;END DO
   !  END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
   Phigrdwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  IF(KST.eq.NSPLTz-1) THEN
       DO KZ =  Ncellz+1, Ncellz+N_ol; DO JY = -1, Ncelly+2; DO IX =  -1, Ncellx+2
       !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
       Phigrdwv(IX,JY,KZ,idm)= cg*2.d0*Tdiff*Phigrd(IX,JY,KZ,idm)+Phiexa(IX,JY,KZ)
       END DO;END DO;END DO
    END IF
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
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam ,deltalen !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre,Phiv,Phi2dt,Phiu
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie,idm
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G


if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz; deltalen=dx1; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx; deltalen=dy1; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy; deltalen=dz1; endif! BT1 = 1; BT2 = 2; VN = 4; end if

  !----kyoukai-----
   !if(hazi==1)then
   !   is = 2
   !   ie = Ncell-3
   !end if
   !if(hazi==2)then
      is = 1
      ie = Ncell-2
   !end if
  !----kyoukai-----
  nu2 = cg * dt / deltalen
  Phipre(:,:,:) = Phiv(:,:,:)
  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     DO Lnum = 1-2, Ncl-2+2
        DO Mnum = 1-2, Ncm-2+2
           do i = is-1,ie+1+1
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
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,1,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1)
     !write(*,*) Phiu(127),'127-2'
     !do i = ist , ndx-ien
      DO Lnum = 1-2, Ncl-2+2
         DO Mnum = 1-2, Ncm-2+2
           do i = is,ie+1
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
     DO Lnum = 1-2, Ncl-2+2
        DO Mnum = 1-2, Ncm-2+2
           do i = is-1-1,ie+1
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
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
     DO Lnum = 1-2, Ncl-2+2
        DO Mnum = 1-2, Ncm-2+2
           do i = is-1,ie
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) + nu2 * (Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
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

!subroutine vanalbada(fg,gradfg,iwx,iwy,iwz)
subroutine vanlimter(inp,in,inm,flmt)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  double precision :: inp,in,inm

  delp = inp-in
  delm = in-inm
  flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
end subroutine vanlimter


subroutine vanalbada(Mnum,Lnum,Phipre,Phigrad,i_sta,i_end,dmein)
  use comvar
  double precision :: delp , delm ,flmt,eps=1.0d-10
  !integer :: i , ip , im , flmt ,eps=1.0d-10
  integer :: Mnum,Lnum,Ncell,i_sta,i_end,k,dmein
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm
  integer :: i , ip , im
  !DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phipre
  DOUBLE PRECISION, dimension(-1:dmein) :: Phigrad


  !if(iwx.eq.1) Ncell = ndx
  !if(iwy.eq.1) Ncell = ndy
  !if(iwz.eq.1) Ncell = ndz

  do i = i_sta-1 , i_end+1
     ix  = iwx*i    + iwy*Lnum + iwz*Mnum
     jy  = iwx*Mnum + iwy*i    + iwz*Lnum
     kz  = iwx*Lnum + iwy*Mnum + iwz*i
     ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
     jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
     kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
     ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
     jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
     kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)

     delp = Phipre(ixp,jyp,kzp)-Phipre(ix,jy,kz)
     delm = Phipre(ix,jy,kz)-Phipre(ixm,jym,kzm)
     flmt = dmax1( 0.d0,(2.d0*delp*delm+eps)/(delp**2+delm**2+eps) )
     !Phigrad(ix,jy,kz) = flmt
     Phigrad(i) = flmt
  end do

end subroutine vanalbada


subroutine fluxcal(preuse,pre,uin,ep,kappa,mode,is,ie)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,-1:ndz) :: ul,ur,pre,preuse,uin
  DOUBLE PRECISION , dimension(-1:ndx) :: slop  !------------- need allocation --------------
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  !uin(:)=0.0d0
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
              DO Lnum = 1-2, Ncl-2+2
              DO Mnum = 1-2, Ncm-2+2
              call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !call vanalbada(pre,slop)
              !do i = is,ie
              ul(ix,jy,kz) = preuse(ix,jy,kz) + 0.25d0 * ep * slop(i) &
                   * ((1.0d0-slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
                   (1.0d0+slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i+1/2
              uin(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !uin(:)=ul(:)
           end if


           if(mode==4) then
              DO Lnum = 1-2, Ncl-2+2
              DO Mnum = 1-2, Ncm-2+2
              call vanalbada(Mnum,Lnum,pre,slop,is,ie,Ncell)
              do i = is-1,ie+1
              ix  = iwx*i    + iwy*Lnum + iwz*Mnum
              jy  = iwx*Mnum + iwy*i    + iwz*Lnum
              kz  = iwx*Lnum + iwy*Mnum + iwz*i
              ixp = iwx*(i+1)+ iwy*Lnum + iwz*Mnum
              jyp = iwx*Mnum + iwy*(i+1)+ iwz*Lnum
              kzp = iwx*Lnum + iwy*Mnum + iwz*(i+1)
              ixm = iwx*(i-1)+ iwy*Lnum + iwz*Mnum
              jym = iwx*Mnum + iwy*(i-1)+ iwz*Lnum
              kzm = iwx*Lnum + iwy*Mnum + iwz*(i-1)
              !do i = ist-1,ndx-ien+1
              ur(ix,jy,kz) = preuse(ix,jy,kz) - 0.25d0 * ep * slop(i) &
                   * ((1.0d0+slop(i)*kappa)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
                   (1.0d0-slop(i)*kappa)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i-1/2
              uin(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !uin(:)=ur(:)
           end if

           if(mode==10) then
              DO Lnum = 1-2, Ncl-2+2
              DO Mnum = 1-2, Ncm-2+2
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
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
              DO Lnum = 1-2, Ncl-2+2
              DO Mnum = 1-2, Ncm-2+2
              !call vanalbada(pre,slop)
              do i = is-2,ie+2
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

subroutine pbphigrd(dt)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
integer :: i,j,k,c=0
double precision :: dxx,dt,prdt=0.0d0,ratio=dsqrt(3.0d0),dddt,rx,ry,rz,rx2,ry2,rz2
character(3) fn
character(5) ciii
character(1) lRANK
!DOUBLE PRECISION , dimension(-1:ndx,1:Dim) :: slop
!DOUBLE PRECISION flx1,fly1,flz1,flx2,fly2,flz2,flx3,fly3,flz3,bc11,bc12,bc21,bc22
!DOUBLE PRECISION fllx1,flly1,fllz1,fllx2,flly2,fllz2,fllx3,flly3,fllz3

!iwx = 0; iwy = 1; iwz = 1
!call BCgrv(101)
!call BCgrv(102)
!bphi1l(j,k,1)
!write(fn,'(i3.3)') NRANK
!write(ciii,'(i5.5)') c
!open(3,file='/work/maedarn/3DMHD/test/bcstep'//fn//'.DAT')
!open(4,file='/work/maedarn/3DMHD/test/bcstepscnd'//fn//'.DAT')
!open(5,file='/work/maedarn/3DMHD/test/source'//fn//ciii//'.DAT')


!******x-BC********
do k = -1,Ncellz+2
do j = -1,Ncelly+2
bphigrdxl(j,k,1,1) = (3.0d0*bphil(j,k,1)-4.0d0*bphil(j,k,0)+bphil(j,k,-1))*0.5d0/dx1
bphigrdxl(j,k,0,1) = (-bphil(j,k,-1)+bphil(j,k,1))*0.5d0/dx1
bphigrdxl(j,k,-1,1) = -(3.0d0*bphil(j,k,-1)-4.0d0*bphil(j,k,0)+bphil(j,k,1))*0.5d0/dx1
bphigrdxr(j,k,Ncellx+2,1) = (3.0d0*bphir(j,k,Ncellx+2)-4.0d0*bphir(j,k,Ncellx+1)+bphir(j,k,Ncellx))*0.5d0/dx1
bphigrdxr(j,k,Ncellx+1,1) = (-bphir(j,k,Ncellx)+bphir(j,k,Ncellx+2))*0.5d0/dx1
bphigrdxr(j,k,Ncellx,1) = -(3.0d0*bphir(j,k,Ncellx)-4.0d0*bphir(j,k,Ncellx+1)+bphir(j,k,Ncellx+2))*0.5d0/dx1
end do
end do

bphigrdxl(:,:,:,2)=bphigrdxl(:,:,:,1)
bphigrdxl(:,:,:,3)=bphigrdxl(:,:,:,1)
bphigrdxl(:,:,:,4)=bphigrdxl(:,:,:,1)
bphigrdxl(:,:,:,5)=bphigrdxl(:,:,:,1)
bphigrdxl(:,:,:,6)=bphigrdxl(:,:,:,1)
bphigrdxl(:,:,:,7)=bphigrdxl(:,:,:,1)
bphigrdxl(:,:,:,8)=bphigrdxl(:,:,:,1)

bphigrdxr(:,:,:,2)=bphigrdxr(:,:,:,1)
bphigrdxr(:,:,:,3)=bphigrdxr(:,:,:,1)
bphigrdxr(:,:,:,4)=bphigrdxr(:,:,:,1)
bphigrdxr(:,:,:,5)=bphigrdxr(:,:,:,1)
bphigrdxr(:,:,:,6)=bphigrdxr(:,:,:,1)
bphigrdxr(:,:,:,7)=bphigrdxr(:,:,:,1)
bphigrdxr(:,:,:,8)=bphigrdxr(:,:,:,1)


!if(IST==0) then
!do i=-1,0
!write(fn,'(i3.3)') NRANK/NSPLTx
!write(lRANK,'(i1.1)') i+1
!open(12,file=dir//'bcphispcsave'//lRANK//fn//'.DAT',FORM='UNFORMATTED')
do k=-1,ndz
do j=-1,ndy
do i=-1,0
   bphigrdxl(j,k,i,1)= cg*2.d0*Tdiff*( bphigrdxl(j,k,i,1) &
                    +(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1+(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,2)= cg*2.d0*Tdiff*(-bphigrdxl(j,k,i,2) &
                    -(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1+(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,3)= cg*2.d0*Tdiff*( bphigrdxl(j,k,i,3) &
                    -(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1+(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,4)= cg*2.d0*Tdiff*(-bphigrdxl(j,k,i,4) &
                    +(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1+(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,5)= cg*2.d0*Tdiff*( bphigrdxl(j,k,i,5) &
                    +(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1-(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,6)= cg*2.d0*Tdiff*(-bphigrdxl(j,k,i,6) &
                    -(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1-(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,7)= cg*2.d0*Tdiff*( bphigrdxl(j,k,i,7) &
                    -(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1-(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)
   bphigrdxl(j,k,i,8)= cg*2.d0*Tdiff*(-bphigrdxl(j,k,i,8) &
                    +(-bphil(j-1,k,i)+bphil(j+1,k,i))*0.5d0/dy1-(-bphil(j,k-1,i)+bphil(j,k+1,i))*0.5d0/dz1)+bphil(j,k,i)

   bphigrdxr(j,k,i+Ncellx+2,1)= cg*2.d0*Tdiff*( bphigrdxr(j,k,i+Ncellx+2,1) &
                    +(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1+(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,2)= cg*2.d0*Tdiff*(-bphigrdxr(j,k,i+Ncellx+2,2) &
                    -(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1+(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,3)= cg*2.d0*Tdiff*( bphigrdxr(j,k,i+Ncellx+2,3) &
                    -(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1+(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,4)= cg*2.d0*Tdiff*(-bphigrdxr(j,k,i+Ncellx+2,4) &
                    +(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1+(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,5)= cg*2.d0*Tdiff*( bphigrdxr(j,k,i+Ncellx+2,5) &
                    +(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1-(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,6)= cg*2.d0*Tdiff*(-bphigrdxr(j,k,i+Ncellx+2,6) &
                    -(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1-(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,7)= cg*2.d0*Tdiff*( bphigrdxr(j,k,i+Ncellx+2,7) &
                    -(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1-(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)
   bphigrdxr(j,k,i+Ncellx+2,8)= cg*2.d0*Tdiff*( -bphigrdxr(j,k,i+Ncellx+2,8) &
                    +(-bphir(j-1,k,i+Ncellx+2)+bphir(j+1,k,i+Ncellx+2))*0.5d0/dy1-(-bphir(j,k-1,i+Ncellx+2)+bphir(j,k+1,i+Ncellx+2))*0.5d0/dz1)&
                    +bphir(j,k,i+Ncellx+2)

!write(12) bphigrdxr(j,k,i+Ncellx+2,1),bphigrdxl(j,k,i,1)
end do
end do
!close(12)
end do
!endif
!******x-BC********


end subroutine pbphigrd


SUBROUTINE PB(pls)
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

DOUBLE PRECISION, dimension(:,:),  allocatable :: dat1,dat2!,bcl1,bcr2
DOUBLE PRECISION, dimension(:,:,:),  allocatable :: bcl1,bcr2
complex*16, dimension(:), allocatable :: spe1,spe2!,bcspel1,bcspel2
complex*16, dimension(:,:), allocatable :: bcspel1,bcspel2
double precision :: kap,temp1r,temp1i,temp2r,temp2i,facG,fac,dxx,dyy,dzz,zp1,zp2,LLl,LLr
double precision, dimension(:,:,:), allocatable :: fint0,fint1
double precision, dimension(:,:,:,:), allocatable :: bcsend
character*4 fnum
character(3) fn,deep
!character(3) lname
character(2) lcRANK
character(1) lRANK
integer pls,mis,mode,klrmax,npl2

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
iwx=1;iwy=1;iwz=1;N_MPI(20)=1;N_MPI(1)=1;CALL BC_MPI(2,1)


!MIYAMA method ---------------------------------------------------

!ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1),speq(Ncellz*NSPLTz,Ncellx))
!ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
!     dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))

ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:Ncellx+3),speq(Ncellz*NSPLTz,-1:Ncellx+2))
ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
     dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))
!allocate(bcsend(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:1,0:NPE-1))

allocate(bcl1(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:loopbc),bcr2(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:loopbc))
allocate(bcspel1(Ncellz*NSPLTz,-1:loopbc),bcspel2(Ncellz*NSPLTz,-1:loopbc))

!if(mode==1) then
!   LL=0.d0
!end if

!if(mode==2) then
!   LL=Lbox
!end if

bcl1(:,:,:)=0.d0; bcr2(:,:,:)=0.d0; bcspel1(:,:)=(0.d0,0.d0); bcspel2(:,:)=(0.d0,0.d0)
!nccy = Ncelly/NSPLTy; nccz = Ncellz/NSPLTz
nccy = Ncelly; nccz = Ncellz
do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
do i=-1,Ncellx+2
  !data(jy,kz,i) = U(i,j,k,1)!-rhomean
  data(jy,kz,i) = U(i,j,k,1)-rhomean
end do;end do;end do

                    !count, blocklength, stride
!CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
!CALL MPI_TYPE_COMMIT(VECU,IERR)

do nlp2 = 0 , loopbc , 1
   klrmax = ((NSPLTy) + NSPLTy * (NSPLTz-1)) * nlp2

CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
do Nlp = 1,NSPLTy*NSPLTz-1

  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

  Nis = JSs + NSPLTy*KSs
  !kls = Nis + 1 !+ pls
  kls = Nis - 1 + klrmax
  Nir = JST + NSPLTy*KST
  !klr = Nir + 1 !+ pls
  klr = Nir - 1 + klrmax

  !***************fordebug*****************
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !***************fordebug*****************

  if(kls.gt.Ncellx+2) then; isend = MPI_PROC_NULL; kls = Ncellx+3; end if
  if(klr.gt.Ncellx+2) then; irecv = MPI_PROC_NULL; klr = Ncellx+3; end if
  CALL MPI_SENDRECV(data(JST*Ncelly+1,KST*Ncellz+1,kls),1,VECU,isend,1, & !send
                    data(JSr*Ncelly+1,KSr*Ncellz+1,klr),1,VECU,irecv,1, MPI_COMM_WORLD,MSTATUS,IERR) !recv
end do
CALL MPI_TYPE_FREE(VECU,IERR)

dxx = dy(1); dyy = dz(1); dzz = dx(1)
facG = -G4pi*dzz

dat1(:,:)=0.d0; dat2(:,:)=0.d0; spe1(:)=(0.d0,0.d0); spe2(:)=(0.d0,0.d0)

!Nir = JST + NSPLTy*KST
!klr = Nir + 1

LLl = dzz*0.5d0 + dzz*dble(pls)
LLr = Lbox-dzz*0.5d0 - dzz*dble(pls)
!LLl = dzz*0.5d0 + dzz*dble(pls)


nn1 = Ncelly*NSPLTy; nn2 = Ncellz*NSPLTz

if(klr.le.Ncellx+2) then

  call rlft3(data(1,1,klr),speq(1,klr),nn1,nn2,1)

  kz = klr
  !zp1 = x(kz)-0.5d0*dzz
  !zp2 = Lbox - zp1
  if((x(kz)-0.5d0*dzz > 0.d0) .and. (klr<1))then
     goto 4269
  end if
  if((x(kz)-0.5d0*dzz < Lbox) .and. (klr>Ncellx))then
     goto 4269
  end if
  zp1 = (x(kz)-0.5d0*dzz)-LLl
  !zp1 = x(kz)-LLl
  zp2 = (x(kz)-0.5d0*dzz)-LLr !おそらくtot .ne. Lboxが原因？(PBini との違い)
  !zp2 = Lbox - zp1
  !zp2 = Lbox+4.0*dzz - zp1
  !zp1 = (x(kz) - 2.0d0*dzz )-0.5d0*dzz + dzz*dble(pls)
  !zp2 = Lbox - zp1
  zp1 = dabs(zp1)
  zp2 = dabs(zp2)
  !zp2 = Lbox/dble(NSPLTx) - zp1
  !write(*,*) zp1,x(kz),kz,'write-zp',NRANK,pls
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
4269 continue
do k=1,Ncellz*NSPLTz
   do j=1,Ncelly*NSPLTy
      dat1(j,k)=dat1(j,k)+bcl1(j,k,nlp2-1)
      dat2(j,k)=dat2(j,k)+bcr2(j,k,nlp2-1)
   end do
   spe1(k)=spe1(k)+bcspel1(k,nlp2-1)
   spe2(k)=spe2(k)+bcspel2(k,nlp2-1)
end do

do k=1,Ncellz*NSPLTz
   do j=1,Ncelly*NSPLTy
      bcl1(j,k,nlp2)=dat1(j,k)
      bcr2(j,k,nlp2)=dat2(j,k)
   end do
   bcspel1(k,nlp2)=spe1(k)
   bcspel2(k,nlp2)=spe2(k)
end do

end do

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


ncy=Ncelly+2
ncz=Ncellz+2

!if(IST==0) then
!write(fn,'(i3.3)') NRANK/NSPLTx
!write(lRANK,'(i1.1)') pls+2
!open(12,file=dir//'bcsave'//lRANK//fn//'.DAT',FORM='UNFORMATTED')
!do k=-1,ncz!; kk= (ncy+1)*k
!do j=-1,ncy!; n = j+kk
do k=-1-2,ncz+2!; kk= (ncy+1)*k
do j=-1-2,ncy+2
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 2
  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 2
  if((j.eq.ncy-1).and.(JST.eq.NSPLTy-1)) jb  = 1
  if((k.eq.ncz-1).and.(KST.eq.NSPLTz-1)) kbb = 1
  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz
  if((j.eq.-1  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy-1
  if((k.eq.-1  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz-1

  if((j.eq.ncy+2).and.(JST.eq.NSPLTy-1)) jb  = 4
  if((k.eq.ncz+2).and.(KST.eq.NSPLTz-1)) kbb = 4
  if((j.eq.ncy+1).and.(JST.eq.NSPLTy-1)) jb  = 3
  if((k.eq.ncz+1).and.(KST.eq.NSPLTz-1)) kbb = 3
  if((j.eq.-2  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy-2
  if((k.eq.-2  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz-2
  if((j.eq.-3  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy-3
  if((k.eq.-3  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz-3

  bphil(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  bphir(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
  !bphi2l(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  !bphi2r(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
!  write(12) bphil(j,k,1-abs(pls)),bphir(j,k,Ncellx+abs(pls))!,bphi2l(j,k,1-abs(pls)), bphi2r(j,k,Ncellx+abs(pls))
end do
end do
!close(12)
!endif

DEALLOCATE(data,speq)
DEALLOCATE(dat1,spe1,dat2,spe2)
DEALLOCATE(bcspel1,bcspel2,bcl1,bcr2)
END SUBROUTINE PB


SUBROUTINE PBini(pls)
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

DOUBLE PRECISION, dimension(:,:),  allocatable :: dat1,dat2!,bcl1,bcr2
DOUBLE PRECISION, dimension(:,:,:),  allocatable :: bcl1,bcr2
complex*16, dimension(:), allocatable :: spe1,spe2!,bcspel1,bcspel2
complex*16, dimension(:,:), allocatable :: bcspel1,bcspel2
double precision :: kap,temp1r,temp1i,temp2r,temp2i,facG,fac,dxx,dyy,dzz,zp1,zp2,LLl,LLr
double precision, dimension(:,:,:), allocatable :: fint0,fint1
double precision, dimension(:,:,:,:), allocatable :: bcsend
character*4 fnum
character(3) fn,deep
!character(3) lname
character(2) lcRANK
character(1) lRANK
integer pls,mis,mode,klrmax,npl2

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
iwx=1;iwy=1;iwz=1;N_MPI(20)=1;N_MPI(1)=1;CALL BC_MPI(2,1)


!MIYAMA method ---------------------------------------------------

!ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1),speq(Ncellz*NSPLTz,Ncellx))
!ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
!     dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))

ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1),speq(Ncellz*NSPLTz,Ncellx+1))
ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
     dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))
!allocate(bcsend(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:1,0:NPE-1))

allocate(bcl1(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:loopbc),bcr2(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:loopbc))
allocate(bcspel1(Ncellz*NSPLTz,-1:loopbc),bcspel2(Ncellz*NSPLTz,-1:loopbc))

!if(mode==1) then
!   LL=0.d0
!end if

!if(mode==2) then
!   LL=Lbox
!end if

bcl1(:,:,:)=0.d0; bcr2(:,:,:)=0.d0; bcspel1(:,:)=(0.d0,0.d0); bcspel2(:,:)=(0.d0,0.d0)
!nccy = Ncelly/NSPLTy; nccz = Ncellz/NSPLTz
nccy = Ncelly; nccz = Ncellz
do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
do i=1,Ncellx
  !data(jy,kz,i) = U(i,j,k,1)!-rhomean
  data(jy,kz,i) = U(i,j,k,1)-rhomean
end do;end do;end do

                    !count, blocklength, stride
!CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
!CALL MPI_TYPE_COMMIT(VECU,IERR)

do nlp2 = 0 , loopbc , 1
   klrmax = ((NSPLTy) + NSPLTy * (NSPLTz-1)) * nlp2

CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)
do Nlp = 1,NSPLTy*NSPLTz-1 ! sent same IST core

  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

  Nis = JSs + NSPLTy*KSs
  !kls = Nis + 1 !+ pls
  kls = Nis + 1 + klrmax
  !kls = Nis - 1 + klrmax
  Nir = JST + NSPLTy*KST
  !klr = Nir + 1 !+ pls
  klr = Nir + 1  + klrmax
  !klr = Nir - 1 + klrmax

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

LLl = dzz*0.5d0 + dzz*dble(pls)   !z'
LLr = Lbox+dzz*0.5d0-dzz*dble(NSPLTx*Ncellx-(IST*Ncellx+mod(pls,Ncellx)))!-dzz*dble(((IST+1)*Ncellx-mod(pls,Ncellx))) !- dzz*dble(NSPLTx*Ncellx-((IST+1)*Ncellx-mod(pls,Ncellx))) !z'



nn1 = Ncelly*NSPLTy; nn2 = Ncellz*NSPLTz

if(klr.le.Ncellx) then

  call rlft3(data(1,1,klr),speq(1,klr),nn1,nn2,1) ! deta=\hat{\rho}(z)

  kz = klr
  !zp1 = x(kz)-0.5d0*dzz
  !zp2 = Lbox - zp1
  if((x(kz)-0.5d0*dzz > 0.d0) .and. (klr<1))then !only use boundary x<0,x>Lbox
     goto 4269
  end if
  if((x(kz)-0.5d0*dzz < Lbox) .and. (klr>Ncellx))then
     goto 4269
  end if
  zp1 = (x(kz)-0.5d0*dzz)-LLl
  !zp1 = x(kz)-LLl
  !zp2 = (x(Ncellx+1-kz)-0.5d0*dzz)!-LLr
  !zp2 = -(x(kz)-0.5d0*dzz)+LLr
  zp2 = Lbox - zp1
  !zp1 = (x(kz) - 2.0d0*dzz )-0.5d0*dzz + dzz*dble(pls)
  !zp2 = Lbox - zp1
  zp1 = dabs(zp1)
  zp2 = dabs(zp2)
  !zp2 = Lbox/dble(NSPLTx) - zp1
  !write(*,*) zp1,x(kz),kz,'write-zp',NRANK,pls
  !------MIYAMA A.10-----------
  !======MIYAMA A.8 (2)========
  temp1r = dat1(1,1) - data(1,1,klr) * 0.5d0*zp1 * facG
  temp1i = dat1(2,1) - data(2,1,klr) * 0.5d0*zp1 * facG
  temp2r = dat2(1,1) - data(1,1,klr) * 0.5d0*zp2 * facG
  temp2i = dat2(2,1) - data(2,1,klr) * 0.5d0*zp2 * facG
  !======MIYMA A.8 (2)========
  do m=1,nn2/2+1; do l=1,nn1/2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)+1.0d-100
    !kap = sqrt(kap)
    dat1(2*l-1,m) = dat1(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG       !\hat{\rho}*G(m,n,z)
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
  !------MIYAMA A.10-----------

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
4269 continue
do k=1,Ncellz*NSPLTz
   do j=1,Ncelly*NSPLTy
      dat1(j,k)=dat1(j,k)+bcl1(j,k,nlp2-1)
      dat2(j,k)=dat2(j,k)+bcr2(j,k,nlp2-1)
   end do
   spe1(k)=spe1(k)+bcspel1(k,nlp2-1)
   spe2(k)=spe2(k)+bcspel2(k,nlp2-1)
end do

do k=1,Ncellz*NSPLTz
   do j=1,Ncelly*NSPLTy
      bcl1(j,k,nlp2)=dat1(j,k)
      bcr2(j,k,nlp2)=dat2(j,k)
   end do
   bcspel1(k,nlp2)=spe1(k)
   bcspel2(k,nlp2)=spe2(k)
end do

end do

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


ncy=Ncelly+2
ncz=Ncellz+2

!if(IST==0) then
!write(fn,'(i3.3)') NRANK/NSPLTx
!write(lRANK,'(i1.1)') pls+2
!open(12,file=dir//'bcsave'//lRANK//fn//'.DAT',FORM='UNFORMATTED')

if((pls/Ncellx)==IST) then
do k=-1,ncz!; kk= (ncy+1)*k
do j=-1,ncy!; n = j+kk
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 2
  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 2
  if((j.eq.ncy-1).and.(JST.eq.NSPLTy-1)) jb  = 1
  if((k.eq.ncz-1).and.(KST.eq.NSPLTz-1)) kbb = 1
  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz
  if((j.eq.-1  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy-1
  if((k.eq.-1  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz-1

  Phiexab1(mod(pls,Ncellx)+1,j,k) = dble(data(jb,kbb,1)) !dble(pls)!dble(data(jb,kbb,1))
  Phiexab2(mod(pls,Ncellx)+1,j,k) = dble(data(jb,kbb,2)) !dble((IST+1)*Ncellx-mod(pls,Ncellx))!dble(data(jb,kbb,2)) dble(NSPLTx*Ncellx-((IST+1)*Ncellx-mod(pls,Ncellx)))!dble(data(jb,kbb,2))
  !bphi2l(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  !bphi2r(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
!  write(12) bphil(j,k,1-abs(pls)),bphir(j,k,Ncellx+abs(pls))!,bphi2l(j,k,1-abs(pls)), bphi2r(j,k,Ncellx+abs(pls))
end do
!write(*,*) Phiexa(1+abs(pls),1,k),Phiexa(Ncellx-abs(pls),1,k)
end do
endif

!if((NSPLTx-1-pls/Ncellx)==IST) then
!do k=-1,ncz!; kk= (ncy+1)*k
!do j=-1,ncy!; n = j+kk
!  jb  = JST*Ncelly + j
!  kbb = KST*Ncellz + k
!  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 2
!  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 2
!  if((j.eq.ncy-1).and.(JST.eq.NSPLTy-1)) jb  = 1
!  if((k.eq.ncz-1).and.(KST.eq.NSPLTz-1)) kbb = 1
!  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
!  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz
!  if((j.eq.-1  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy-1
!  if((k.eq.-1  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz-1

  !Phiexab1(mod(pls,Ncellx)+1,j,k) = dble(data(jb,kbb,1)) !dble(pls)!dble(data(jb,kbb,1))
!  Phiexab2(mod(pls,Ncellx)+1,j,k) =dble(data(jb,kbb,2))!Lbox+dzz*0.5d0-dzz*dble(NSPLTx*Ncellx-(IST*Ncellx+mod(pls,Ncellx)))
!dble((IST+1)*Ncellx-mod(pls,Ncellx))!dble(data(jb,kbb,2)) dble(NSPLTx*Ncellx-((IST+1)*Ncellx-mod(pls,Ncellx)))!dble(data(jb,kbb,2))
  !bphi2l(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  !bphi2r(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
!  write(12) bphil(j,k,1-abs(pls)),bphir(j,k,Ncellx+abs(pls))!,bphi2l(j,k,1-abs(pls)), bphi2r(j,k,Ncellx+abs(pls))
!end do
 !write(*,*) Phiexa(1+abs(pls),1,k),Phiexa(Ncellx-abs(pls),1,k)
!end do
!endif
!close(12)
!endif

DEALLOCATE(data,speq)
DEALLOCATE(dat1,spe1,dat2,spe2)
DEALLOCATE(bcspel1,bcspel2,bcl1,bcr2)
END SUBROUTINE PBini



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

subroutine timesource(Phiv,source,dt,mode)
  use comvar
  !use grvver
  integer i,mode,j,k
  double precision :: dt,sdt,mindt,maxdt , epsl = 1.0d-4
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phiv,source
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G

  !mindt=1000.0d0
  maxdt=0.0d0

  if(mode==1) then
     do k=1,ndx-2
        do j=1,ndx-2
           do i=1,ndx-2
              if((source(i,j,k) .ne. 0.0d0) .and. (Phiv(i,j,k) .ne. 0.0d0))then
                 sdt = sourratio*dabs(Phiv(i,j,k) / (cg * G4pi * source(i,j,k) ))
                 !sdt = 0.2d0*dabs(Phiv(i)) / (cg * G4pi * source(i) )
                 !mindt=dmin1(mindt,sdt)
                 maxdt=dmax1(maxdt,sdt)
              end if
           end do
        end do
     end do
     if( (maxdt < dt) .and. (maxdt .gt. 0.0d0)) then
        dt = sdt
     end if
  end if


  if(mode==2) then
     do k=1,ndx-2
        do j=1,ndx-2
           do i=1,ndx-2
              if((source(i,j,k) .ne. 0.0d0) .and. (Phiv(i,j,k) .ne. 0.0d0))then
                 sdt = sourratio*dabs(Phiv(i,j,k) / ( cg * source(i,j,k) ))
                 !sdt = 0.05d0*dabs(Phiv(i)) / ( cg * source(i) )
                 !mindt=dmin1(mindt,sdt)
                 maxdt=dmax1(maxdt,sdt)
              end if
           end do
        end do
     end do
     write(*,*) maxdt,'maxdt'
     if( (maxdt < dt) .and. (maxdt .gt. 0.0d0)) then
        dt = sdt
     end if
  end if


  write(*,*) 'time source' , dt
end subroutine timesource

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


do kk=-1,Ncellz+2; k=KST*Ncellz+kk
do jj=-1,Ncelly+2; j=JST*Ncelly+jj
do ii=-1,Ncellz+2; i=IST*Ncellx+ii
    if((j.eq.NSPLTy*Ncelly+2).and.(JST.eq.NSPLTy-1)) j  = 2
    if((k.eq.NSPLTz*Ncellz+2).and.(KST.eq.NSPLTz-1)) k = 2
    if((j.eq.NSPLTy*Ncelly+1).and.(JST.eq.NSPLTy-1)) j  = 1
    if((k.eq.NSPLTy*Ncelly+1).and.(KST.eq.NSPLTz-1)) k = 1
    if((j.eq.0  ).and.(JST.eq.0       )) j  = Ncelly*NSPLTy
    if((k.eq.0  ).and.(KST.eq.0       )) k = Ncellz*NSPLTz
    if((j.eq.-1  ).and.(JST.eq.0       )) j  = Ncelly*NSPLTy-1
    if((k.eq.-1  ).and.(KST.eq.0       )) k = Ncellz*NSPLTz-1
    if((i.eq.NSPLTx*Ncellx+2).and.(IST.eq.NSPLTx-1)) i  = 2
    if((i.eq.NSPLTx*Ncellx+1).and.(IST.eq.NSPLTx-1)) i  = 1
    if((i.eq.0  ).and.(JST.eq.0       )) i  = Ncellx*NSPLTx
    if((i.eq.-1  ).and.(JST.eq.0       )) i  = Ncellx*NSPLTx-1

    Phiwv(ii,jj,kk,num1)=u1(i,j,k,num1)
end do;end do;end do


IF(IST.eq.0) THEN
   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-2, 0
   Phiwv(IX,JY,KZ,num1)= bphil(JY,KZ,IX)
   END DO;END DO;END DO
END IF
IF(IST.eq.NSPLTx-1) THEN
   DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+2
   Phiwv(IX,JY,KZ,num1)= bphir(JY,KZ,IX)
   END DO;END DO;END DO
END IF
enddo

END SUBROUTINE collectPhi
