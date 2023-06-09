subroutine SELFGRAVWAVE(dt,mode)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'
  integer :: mode,MRANK,count=0
  DOUBLE PRECISION  :: dt,dxi
  !INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt
  INTEGER :: MSTATUS(MPI_STATUS_SIZE)
  DOUBLE PRECISION  :: VECU
  character(3) NPENUM
  character(6) countcha
  double precision tfluid , cs
  double precision dt_mpi_gr(0:NPE-1),dt_gat_gr(0:NPE-1),maxcs,tcool,cgtime!,sourcedt
  double precision :: ave1,ave1pre,ave2(0:NPE-1),ave,avepre,ave2_gather(0:NPE-1) , eps=1.0d-3
  !double precision , dimension(:,:,:) , allocatable :: stbPhi
  !double precision , dimension(-1:Ncellx+2,-1:Ncelly,-1:Ncellz) :: Phipregrad,Phipregraddum
  !**************** INITIALIZEATION **************
  if(mode==0) then
     Phicgp(:,:,:,:)=0.0d0
     Phicgp(:,:,:,:)=0.0d0
     Phi1step(:,:,:,:)=0.0d0
     Phi2step(:,:,:,:)=0.0d0
     Phi(:,:,:) = 0.0d0
     source(:,:,:,:)=0.0d0
     sourcedt(:,:,:,:)=0.0d0
     sourcedt2(:,:,:,:)=0.0d0
  end if
  !**************** INITIALIZEATION **************

  !****************read INITIAL CONDITION**************
  if(mode==1) then
     WRITE(NPENUM,'(I3.3)') NRANK
     open(unit=8,file=dir//'PHIINI/INIPHIcgp'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=18,file=dir//'PHIINI/INIPHIcgm'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file=dir//'PHIINI/INIPHI1step'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=38,file=dir//'PHIINI/INIPHI2step'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(8) (Phicgp(i,j,k,1),Phicgp(i,j,k,2),Phicgp(i,j,k,3),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(8)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(18) (Phicgm(i,j,k,1),Phicgm(i,j,k,2),Phicgm(i,j,k,3),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(18)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(28) (Phi1step(i,j,k,1),Phi1step(i,j,k,2),Phi1step(i,j,k,3),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(28)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(38) (Phi2step(i,j,k,1),Phi2step(i,j,k,2),Phi2step(i,j,k,3),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(38)
  end if
  !****************read INITIAL CONDITION**************

  !****************GRAVITY SOLVER*****************

  if(mode==2) then
     N_MPI(20)=1; N_MPI(1)=1
     iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(2,1)
     !---debug---
     !call  SELFGRAVWAVE(0.0d0,4)
     !write(*,*) '------pb1-------' ,Nrank

     !calcurate Phi-bc
     Call PB( 0)
     Call PB(-1)
     Call PB(-2)
     !call BCgrv(101,1)
     !call BCgrv(101,2)
     !call BCgrv(101,3)
     call pbstep(dt)
     !call pbstep()
     !calcurate Phi-bc
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
     !open(unit=28,file=dir//'PHIINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file=dir//'PHIDTINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=8,file=dir//'PHIINI/INIPHIcgp'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=18,file=dir//'PHIINI/INIPHIcgm'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file=dir//'PHIINI/INIPHI1step'//NPENUM//countcha//'.DAT')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file=dir//'PHIINI/INIPHI2step'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
     iwx = 1; iwy = 1; iwz = 1
     !call BCgrv(101)
     !call BCgrv(102)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           do i = -1, Ncellx+2
  write(28,*) i,j,k, sngl(Phicgp(i,j,k,1)),sngl(Phicgp(i,j,k,2)),sngl(Phicgp(i,j,k,3)),&
    sngl(Phi1step(i,j,k,1)),sngl(Phi1step(i,j,k,2)),sngl(Phi1step(i,j,k,3)),&
                sngl(((Phicgp(i,j,k,1)+Phicgp(i,j,k,2)+Phicgp(i,j,k,3))))
                !sngl(Phicgp(i,j,k)),sngl(Phicgm(i,j,k)),sngl(Phi1step(i,j,k)),sngl(Phi2step(i,j,k)),&
                !sngl((Phicgp(i,j,k)+Phicgm(i,j,k))*0.5d0)!,i=-1,Ncellx+2)
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


!***************pointerforPB**************
if(mode.eq.10) then
  !call pinter(Nmem1,Nmem2,Ncellx,Ncelly,Ncellz)
end if
!***************pointerforPB**************

if(mode==11) then
   call cllsub(2,dt,0,0)
end if
end subroutine SELFGRAVWAVE


subroutine slvmuscle(dt)
  use comvar
  use slfgrv
  INCLUDE 'mpif.h'
  double precision :: dt,dtratio=dsqrt(3.0d0)
  integer :: i=0,n,m,l
  double precision :: rho(-1:ndx,-1:ndy,-1:ndz)
  integer ifEVOgrv,ifEVOgrv2

  ifEVOgrv  = 1
  ifEVOgrv2 = 1

  !do l=-1,ndz
  !do m=-1,ndy
  !do n=-1,ndx
  !   rho(n,m,l) = U(n,m,l,1)
  !end do;end do;end do

  do l=-1,ndz
  do m=-1,ndy
  do n=-1,ndx
     rho(n,m,l) = U(n,m,l,1)
     source(n,m,l,1)=-source(n,m,l,1)+G4pi*rho(n,m,l)
     source(n,m,l,2)=-source(n,m,l,2)+G4pi*rho(n,m,l)
     source(n,m,l,3)=-source(n,m,l,3)+G4pi*rho(n,m,l)
     !source(n,m,l,1)=source(n,m,l,1)-G4pi*rho(n,m,l)
     !source(n,m,l,2)=source(n,m,l,2)-G4pi*rho(n,m,l)
     !source(n,m,l,3)=source(n,m,l,3)-G4pi*rho(n,m,l)
  end do;end do;end do

  !write(*,*) 'OK or NOT 0'
  !call muslcslv1D(Phi1step(-1,-1,-1,1),rho,dt*0.5d0*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,2),rho,dt*0.5d0*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,3),rho,dt*0.5d0*dtratio,3,2)
  call muslcslv1D(Phi1step(-1,-1,-1,1),source(-1,-1,-1,1),dt*0.5d0*dtratio,4,2)
  call muslcslv1D(Phi1step(-1,-1,-1,2),source(-1,-1,-1,2),dt*0.5d0*dtratio,4,2)
  call muslcslv1D(Phi1step(-1,-1,-1,3),source(-1,-1,-1,3),dt*0.5d0*dtratio,4,2)
  !write(*,*) 'OK or NOT 1'
  iwx=1; iwy=1; iwz=1
  call BCgrv(102,1)
  call BCgrv(102,2)
  call BCgrv(102,3)
  !write(*,*) 'OK or NOT 2'

  call cllsub(4,dt*0.5d0*dtratio,ifEVOgrv,ifEVOgrv2)

  !write(*,*) 'OK or NOT 3'
  call muslcslv1D(Phicgp(-1,-1,-1,1),Phi1step(-1,-1,-1,1),dt*dtratio,3,2)
  call muslcslv1D(Phicgp(-1,-1,-1,2),Phi1step(-1,-1,-1,2),dt*dtratio,3,2)
  call muslcslv1D(Phicgp(-1,-1,-1,3),Phi1step(-1,-1,-1,3),dt*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,1),source(-1,-1,-1,1),dt*0.5d0*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,2),source(-1,-1,-1,2),dt*0.5d0*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,3),source(-1,-1,-1,3),dt*0.5d0*dtratio,3,2)
  iwx=1; iwy=1; iwz=1
  call BCgrv(101,1)
  call BCgrv(101,2)
  call BCgrv(101,3)

  !write(*,*) 'OK or NOT 4'
  !call cllsub(5,dt*dtratio/3.0d0,dtifEVOgrv,ifEVOgrv2)
  call cllsub(5,dt*dtratio,dtifEVOgrv,ifEVOgrv2)

  !write(*,*) 'OK or NOT 5'
  !call muslcslv1D(Phi1step(-1,-1,-1,1),rho,dt*0.5d0*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,2),rho,dt*0.5d0*dtratio,3,2)
  !call muslcslv1D(Phi1step(-1,-1,-1,3),rho,dt*0.5d0*dtratio,3,2)
  call muslcslv1D(Phi1step(-1,-1,-1,1),source(-1,-1,-1,1),dt*0.5d0*dtratio,4,2)
  call muslcslv1D(Phi1step(-1,-1,-1,2),source(-1,-1,-1,2),dt*0.5d0*dtratio,4,2)
  call muslcslv1D(Phi1step(-1,-1,-1,3),source(-1,-1,-1,3),dt*0.5d0*dtratio,4,2)
  iwx=1; iwy=1; iwz=1
  call BCgrv(102,1)
  call BCgrv(102,2)
  call BCgrv(102,3)

  call cllsub(4,dt*0.5d0*dtratio,ifEVOgrv,ifEVOgrv2)
  !write(*,*) 'OK or NOT 6'
end subroutine slvmuscle


subroutine cllsub(mode,dt,ifEVOgrv,ifEVOgrv2)
  use comvar
  use slfgrv
  INCLUDE 'mpif.h'
  integer :: mode !,ifEVOgrv=1,ifEVOgrv2=1
  double precision dt
  integer :: n,m,l
  double precision rho(-1:ndx,-1:ndy,-1:ndz)
  integer ifEVOgrv,ifEVOgrv2

  !write(*,*) 'cllsub',mode,dt,ifEVOgrv,ifEVOgrv2

  if(mode==4) then
     if(ifEVOgrv.eq.1) then
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phi1step(-1,-1,-1,1),rho,dt,2,1)
        call BCgrv(102,1)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phi1step(-1,-1,-1,2),rho,dt,2,2)
        call BCgrv(102,2)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phi1step(-1,-1,-1,3),rho,dt,2,2)
        call BCgrv(102,3)
        goto 1470
     end if
1470 continue
  end if


  if(mode==5) then
     if(ifEVOgrv2.eq.1) then
        iwx=1; iwy=0; iwz=0
        call muslcslv1D(Phicgp(-1,-1,-1,1),Phi1step(-1,-1,-1,1),dt,1,1)
        call muslcslv1D(Phicgp(-1,-1,-1,2),Phi1step(-1,-1,-1,2),dt,1,1)
        call muslcslv1D(Phicgp(-1,-1,-1,3),Phi1step(-1,-1,-1,3),dt,1,1)
        call BCgrv(101,1)
        call BCgrv(101,2)
        call BCgrv(101,3)
        iwx=0; iwy=1; iwz=0
        call muslcslv1D(Phicgp(-1,-1,-1,1),Phi1step(-1,-1,-1,1),dt,1,1)
        call muslcslv1D(Phicgp(-1,-1,-1,2),Phi1step(-1,-1,-1,2),dt,1,1)
        call muslcslv1D(Phicgp(-1,-1,-1,3),Phi1step(-1,-1,-1,3),dt,1,1)
        call BCgrv(101,1)
        call BCgrv(101,2)
        call BCgrv(101,3)
        iwx=0; iwy=0; iwz=1
        call muslcslv1D(Phicgp(-1,-1,-1,1),Phi1step(-1,-1,-1,1),dt,1,1)
        call muslcslv1D(Phicgp(-1,-1,-1,2),Phi1step(-1,-1,-1,2),dt,1,1)
        call muslcslv1D(Phicgp(-1,-1,-1,3),Phi1step(-1,-1,-1,3),dt,1,1)
        call BCgrv(101,1)
        call BCgrv(101,2)
        call BCgrv(101,3)
        goto 1788
     end if
1788 continue
  end if
end subroutine cllsub

subroutine BCgrv(mode,idm)
  use comvar
  use mpivar
  use slfgrv
  INCLUDE 'mpif.h'
  integer ::  N_ol=2,i,mode,j,k,idm
  INTEGER :: MSTATUS(MPI_STATUS_SIZE)
  DOUBLE PRECISION  :: VECU
  double precision , dimension(1:2) :: pl,pr

  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

  if(mode==100) then

IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL


  CALL MPI_SENDRECV(Phi(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
                    Phi(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  CALL MPI_SENDRECV(Phicgm(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
!       Phicgm(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.0) THEN
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 1,1
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phi(IX,JY,KZ)= bphil(JY,KZ,IX)
!     Phicgm(IX,JY,KZ)= bphi2l(JY,KZ,IX)
     END DO;END DO;END DO
  END IF

  CALL MPI_SENDRECV(Phi(1            ,-1,-1),1,VECU,LEFT,1, &
       Phi(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  CALL MPI_SENDRECV(Phicgm(1            ,-1,-1),1,VECU,LEFT,1, &
!       Phicgm(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.NSPLTx-1) THEN
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx, Ncellx+N_ol,1
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phi(IX,JY,KZ)= bphir(JY,KZ,IX)
!     Phicgm(IX,JY,KZ)= bphi2r(JY,KZ,IX)
     END DO;END DO;END DO
  END IF

CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF


IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
       Phi(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  CALL MPI_SENDRECV(Phicgm(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
!       Phicgm(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi(-1,1            ,-1),1,VECU,BOTM,1, &
       Phi(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  CALL MPI_SENDRECV(Phicgm(-1,1            ,-1),1,VECU,BOTM,1, &
!       Phicgm(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
       Phi(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  CALL MPI_SENDRECV(Phicgm(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
!       Phicgm(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi(-1,-1,1            ),1,VECU,DOWN,1, &
       Phi(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
!  CALL MPI_SENDRECV(Phicgm(-1,-1,1            ),1,VECU,DOWN,1, &
!       Phicgm(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt
END IF
endif

if(mode==101) then

IF(iwx.EQ.1) THEN
      !dim=1
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL


  CALL MPI_SENDRECV(Phicgp(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
                    Phicgp(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
       Phicgm(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.0) THEN
     if(idm==1)then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phicgp(IX,JY,KZ,idm)= bphixl(JY,KZ,IX,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
     END DO;END DO;END DO
     endif
     if(idm==2 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phicgp(IX,JY,KZ,idm)=  Phicgp(1,JY,KZ,idm)!Phicgp(IX+1,JY,KZ,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2r(JY,KZ,IX)
     END DO;END DO;END DO
     endif
  END IF

  CALL MPI_SENDRECV(Phicgp(1            ,-1,-1,idm),1,VECU,LEFT,1, &
       Phicgp(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(1            ,-1,-1,idm),1,VECU,LEFT,1, &
       Phicgm(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.NSPLTx-1) THEN
     if(idm==1)then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phicgp(IX,JY,KZ,idm)= bphixr(JY,KZ,IX,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2r(JY,KZ,IX)
     END DO;END DO;END DO
     write(*,*) Phicgp(Ncellx+1,1,1,idm),idm,'bcbc-kakunin',bphixr(1,1,Ncellx+1,idm)
     endif
     if(idm==2 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+N_ol, Ncellx+1,-1
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phicgp(IX,JY,KZ,idm)=  Phicgp(Ncellx,JY,KZ,idm)!Phicgp(IX-1,JY,KZ,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2r(JY,KZ,IX)
     END DO;END DO;END DO
     endif
     !write(*,*) NRANK,Phicgp(-1,1,1),Phicgm(-1,1,1),Phicgp(Ncellx+2,1,1),Phicgm(Ncellx+2,1,1),'---bcc1---'
  END IF

CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF


IF(iwy.EQ.1) THEN
   !idm=2
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
       Phicgp(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
       Phicgm(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
  IF(JST.eq.0) THEN
     if(idm==1 .or. idm==3) then
     !DO KZ = -1, Ncellz+2; DO IX = -1, Ncellx+2; DO JY = 1-N_ol, 0
      DO KZ = -1, Ncellz+2; DO JY = 1-N_ol, 0 ; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phicgp(IX,JY,KZ,idm)= Phicgp(IX,1,KZ,idm)!Phicgp(IX,JY+1,KZ,idm)
    !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
     END DO;END DO;END DO
     endif
!     if(idm==2) then
     !DO KZ = -1, Ncellz+2; DO IX = -1, Ncellx+2; DO JY = 1-N_ol, 0
!      DO KZ = -1, Ncellz+2; DO JY = 1-N_ol, 0 ; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
!     Phicgp(IX,JY,KZ,idm)= Phi(IX,JY,KZ)-Phicgp(IX,JY,KZ,1)-Phicgp(IX,JY,KZ,3)!Phicgp(IX,JY+1,KZ,idm)
    !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
!     END DO;END DO;END DO
!     endif
     !if(idm==2) then
     !!DO KZ = -1, Ncellz+2; DO IX = -1, Ncellx+2; DO JY = 1-N_ol, 0
     ! DO KZ = -1, Ncellz+2; DO JY = 1-N_ol, 0 ; DO IX = -1, Ncellx+2
     !!DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     !Phicgp(IX,JY,KZ,idm)= Phicgp(IX,JY+1,KZ,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
     !END DO;END DO;END DO
     !endif
  END IF
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,1            ,-1,idm),1,VECU,BOTM,1, &
       Phicgp(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,1            ,-1,idm),1,VECU,BOTM,1, &
       Phicgm(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
  IF(JST.eq.NSPLTy-1) THEN
     !DO KZ = -1, Ncellz+2; DO IX = -1, Ncellx+2; DO JY = Ncelly+1, Ncelly+N_ol,1
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
!     if(idm==2) then
!     DO KZ = -1, Ncellz+2; DO JY = Ncelly+N_ol, Ncelly+1,-1 ; DO IX = -1, Ncellx+2
!     Phicgp(IX,JY,KZ,idm)= Phi(IX,JY,KZ)-Phicgp(IX,JY,KZ,1)-Phicgp(IX,JY,KZ,3)!Phicgp(IX,JY-1,KZ,idm)!bphiyl(KZ,IX,JY,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
!     END DO;END DO;END DO
!     endif
     if(idm==1 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = Ncelly+N_ol, Ncelly+1,-1 ; DO IX = -1, Ncellx+2
     Phicgp(IX,JY,KZ,idm)= Phicgp(IX,Ncelly,KZ,idm)!Phicgp(IX,JY-1,KZ,idm)!bphiyl(KZ,IX,JY,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
     END DO;END DO;END DO
     endif
  END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt
END IF


IF(iwz.EQ.1) THEN
   !idm=3
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
       Phicgp(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
       Phicgm(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
  IF(KST.eq.0) THEN
!     if(idm==3) then
     !DO IX = -1, Ncellx+2; DO JY = -1, Ncelly+2; DO KZ = 1-N_ol, 0
!     DO KZ = 1-N_ol, 0; DO JY = -1, Ncelly+2 ; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
!     Phicgp(IX,JY,KZ,idm)= Phi(IX,JY,KZ)-Phicgp(IX,JY,KZ,1)-Phicgp(IX,JY,KZ,2)!Phicgp(IX,JY,KZ+1,idm)!bphizl(IX,JY,KZ,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
!     END DO;END DO;END DO
!     end if
     if(idm==1 .or. idm==2) then
     !DO IX = -1, Ncellx+2; DO JY = -1, Ncelly+2; DO KZ = 1-N_ol, 0
     DO KZ = 1-N_ol, 0; DO JY = -1, Ncelly+2 ; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phicgp(IX,JY,KZ,idm)=  Phicgp(IX,JY,1,idm)!Phicgp(IX,JY,KZ+1,idm)!bphizl(IX,JY,KZ,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
     END DO;END DO;END DO
     end if
  END IF
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,-1,1            ,idm),1,VECU,DOWN,1, &
       Phicgp(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,-1,1            ,idm),1,VECU,DOWN,1, &
       Phicgm(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
  IF(KST.eq.NSPLTz-1) THEN
!     if(idm==3) then
!     DO KZ = Ncellz+N_ol, Ncellz+1,-1; DO JY = -1, Ncelly+2 ; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
!     Phicgp(IX,JY,KZ,idm)= Phi(IX,JY,KZ)-Phicgp(IX,JY,KZ,1)-Phicgp(IX,JY,KZ,2)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
!     END DO;END DO;END DO
!     endif
     !DO IX = -1, Ncellx+2; DO JY = -1, Ncelly+2; DO KZ = Ncellz+1, Ncellz+N_ol,1
     if(idm==1 .or. idm==2)then
     DO KZ = Ncellz+N_ol, Ncellz+1,-1; DO JY = -1, Ncelly+2 ; DO IX = -1, Ncellx+2
    !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phicgp(IX,JY,KZ,idm)= Phicgp(IX,JY,Ncellz,idm)!Phicgp(IX,JY,KZ-1,idm)!bphizl(IX,JY,KZ,idm)
     !Phicgm(IX,JY,KZ,idm)= bphi2l(JY,KZ,IX)
     END DO;END DO;END DO
     endif
  END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt
END IF
endif



if(mode==102) then

IF(iwx.EQ.1) THEN
      !idm=1
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
!********************************  BC for the leftsides of domains  *****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
       Phi1step(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(Ncellx+1-N_ol,-1,-1,idm),1,VECU,RIGT,1, &
       Phi2step(       1-N_ol,-1,-1,idm),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.0) THEN
     if(idm==1) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0,1
 ! DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phi1step(IX,JY,KZ,idm)= bstepxl(JY,KZ,IX,idm)
     !Phi2step(IX,JY,KZ,idm)= bstepxl(JY,KZ,IX)
  END DO;END DO;END DO
  endif
  if(idm==2 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0,1
        ! DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
        Phi1step(IX,JY,KZ,idm)= Phi1step(1,JY,KZ,idm)
        !Phi2step(IX,JY,KZ,idm)= bstepxl(JY,KZ,IX)
     END DO;END DO;END DO
  end if
  END IF
  CALL MPI_SENDRECV(Phi1step(1            ,-1,-1,idm),1,VECU,LEFT,1, &
       Phi1step(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(1            ,-1,-1,idm),1,VECU,LEFT,1, &
       Phi2step(Ncellx+1     ,-1,-1,idm),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)


  IF(IST.eq.NSPLTx-1) THEN
     if(idm==1) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phi1step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX,idm)
     !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
  END DO;END DO;END DO
  endif
  if(idm==2 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
        !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
        Phi1step(IX,JY,KZ,idm)=Phi1step(Ncellx,JY,KZ,idm)
        !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
     END DO;END DO;END DO
  endif
END IF

!write(*,*) NRANK,'--bc--5--',Phi1step(0,1,1),Phi2step(0,1,1),Phi1step(Ncellx+1,1,1),Phi2step(Ncellx,1,1)

!************************************************************************
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF


IF(iwy.EQ.1) THEN
   !idm=2
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
  CALL MPI_SENDRECV(Phi1step(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
       Phi1step(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
       Phi2step(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  IF(JST.eq.0) THEN
     if(idm==2) then
     DO KZ = -1, Ncellz+2; DO JY = 1-N_ol, 0,1; DO IX = -1, Ncellx+2
 ! DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
        Phi1step(IX,JY,KZ,idm)= bstepyl(KZ,IX,JY,idm)
        !Phi1step(IX,JY,KZ,idm)= 0.0d0
     !Phi2step(IX,JY,KZ,idm)= bstepxl(JY,KZ,IX)
  END DO;END DO;END DO
endif
if(idm==1 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
        !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
        Phi1step(IX,JY,KZ,idm)=Phi1step(IX,1,KZ,idm)
        !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
     END DO;END DO;END DO
     END if
  END IF
!**************************************  BC for the upsides of domains  ****
  CALL MPI_SENDRECV(Phi1step(-1,1            ,-1,idm),1,VECU,BOTM,1, &
       Phi1step(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,1            ,-1,idm),1,VECU,BOTM,1, &
       Phi2step(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  IF(JST.eq.NSPLTy-1) THEN
     if(idm==2) then
    !DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
     DO KZ = -1, Ncellz+2; DO JY = Ncelly+1, Ncelly+N_ol,1; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phi1step(IX,JY,KZ,idm)= bstepyr(KZ,IX,JY,idm)
     !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
  END DO;END DO;END DO
  endif
  if(idm==1 .or. idm==3) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
        !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
        Phi1step(IX,JY,KZ,idm)=Phi1step(IX,Ncelly,KZ,idm)
        !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
     END DO;END DO;END DO
     END if
  END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  TOP = TOPt; BOTM = BOTMt
END IF

IF(iwz.EQ.1) THEN
   !idm=3
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
  CALL MPI_SENDRECV(Phi1step(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
       Phi1step(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
       Phi2step(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  IF(KST.eq.0) THEN
     if(idm==3) then
     DO KZ = 1-N_ol, 0,1; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
 ! DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phi1step(IX,JY,KZ,idm)= bstepzl(IX,JY,KZ,idm)
     !Phi2step(IX,JY,KZ,idm)= bstepxl(JY,KZ,IX)
  END DO;END DO;END DO
       endif
if(idm==1 .or. idm==2) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
        !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
        Phi1step(IX,JY,KZ,idm)=Phi1step(IX,JY,1,idm)
        !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
     END DO;END DO;END DO
     END if
  END IF
!**************************************  BC for the upsides of domains  ****
  CALL MPI_SENDRECV(Phi1step(-1,-1,1            ,idm),1,VECU,DOWN,1, &
       Phi1step(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,-1,1            ,idm),1,VECU,DOWN,1, &
       Phi2step(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  IF(KST.eq.NSPLTz-1) THEN
     if(idm==3) then
     !DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
     DO KZ =Ncellz+1, Ncellz+N_ol,1 ; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phi1step(IX,JY,KZ,idm)= bstepzr(IX,JY,KZ,idm)
     !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
  END DO;END DO;END DO
endif
if(idm==1 .or. idm==2) then
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol,1
        !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
        Phi1step(IX,JY,KZ,idm)=Phi1step(IX,JY,Ncellz,idm)
        !Phi2step(IX,JY,KZ,idm)= bstepxr(JY,KZ,IX)
     END DO;END DO;END DO
     END if
 END IF
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
  UP = UPt; DOWN = DOWNt
END IF
endif
!----debug----
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
end subroutine BCgrv


subroutine muslcslv1D(Phiv,source,dt,mode,hazi)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam  !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  !DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre,fluxphi&
       ,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie,idm
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G


  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy; endif! BT1 = 1; BT2 = 2; VN = 4; end if

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
  nu2 = cg * dt / deltalength
  Phipre(:,:,:) = Phiv(:,:,:)
  !write(name,'(i5.5)') cnt
  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
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
      DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
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
     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
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
              !do i=ist-1,ndx-ien+1
              Phi2dt(ix,jy,kz) = Phipre(ix,jy,kz) + 0.5d0 * nu2 * ( Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
           end do
        end DO
     end DO
     !------------calcurate dt/2------------
     call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,1.d0/3.0d0,4,is,ie)
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4)

     !do i = ist , ndx-ien
     DO Lnum = 1, Ncl-2
        DO Mnum = 1, Ncm-2
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
              Phiv(ix,jy,kz) = Phipre(ix,jy,kz) + nu2 * (Phiu(ixp,jyp,kzp) - Phiu(ix,jy,kz))
           end do
        end DO
     end DO

     !do i=-1,ndx
     !   write(202,*) i, Phiv(i)
     !end do

  end if
  !------------ul.solver.-cg-------------


  !--------------source------------------
  if(mode==3) then
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i = 1,ndx-2
              !Phiv(i,j,k) =  cg * G4pi * source(i,j,k) * dt * phiratio + Phipre(i,j,k)
              !Phiv(i,j,k) =  cg * source(i,j,k) * dt * phiratio + Phipre(i,j,k)
              Phiv(i,j,k) =  cg * source(i,j,k) * dt + Phipre(i,j,k)
              !write(*,*) source(i,j,k),Phipre(i,j,k),Phiv(i,j,k),'-mode3-'
           end do
        end do
     end do
  end if

  if(mode==4) then
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i = 1,ndx-2
              Phiv(i,j,k) = -cg * source(i,j,k) * dt + Phipre(i,j,k)
           end do
        end do
     end do
  end if
  !--------------source------------------

!  close(201)
!  close(202)
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
              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
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
              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
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
              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
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
              DO Lnum = 1, Ncl-2
              DO Mnum = 1, Ncm-2
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

subroutine pbstep(dt)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
integer :: i,j,k,c=0
double precision :: dxx,dt,prdt=0.0d0,ratio=dsqrt(3.0d0),dddt,rx,ry,rz,rx2,ry2,rz2
character(3) fn
character(5) ciii
!DOUBLE PRECISION , dimension(-1:ndx,1:Dim) :: slop
DOUBLE PRECISION flx1,fly1,flz1,flx2,fly2,flz2,flx3,fly3,flz3,bc11,bc12,bc21,bc22
DOUBLE PRECISION fllx1,flly1,fllz1,fllx2,flly2,fllz2,fllx3,flly3,fllz3

!dxx=dx(1)
!dxx=-dx(1)
!dxx=-deltalength
dxx=deltalength
!iwx = 0; iwy = 1; iwz = 1
!call BCgrv(101)
!call BCgrv(102)
!bphi1l(j,k,1)
write(fn,'(i3.3)') NRANK
write(ciii,'(i5.5)') c
!open(3,file=dir//'bcstep'//fn//'.DAT')
!open(4,file=dir//'bcstepscnd'//fn//'.DAT')
!open(5,file=dir//'source'//fn//ciii//'.DAT')


sourcedt2(:,:,:,:)=sourcedt(:,:,:,:)
sourcedt(:,:,:,:)=Phi1step(:,:,:,:)
source(:,:,:,:)=0.0d0

dddt=1.0d0/dt/ratio/cg
write(*,*) 'dddt',dddt

!do
!call vanalbada(Mnum,Lnum,Phipre,Phigrad,i_sta,i_end,dmein)

do k = 1,ndz-2
do j = 1,ndy-2
do i = 1,ndx-2

!call vanlimter(sourcedt(i+1,j,k,1),sourcedt(i,j,k,1),sourcedt(i-1,j,k,1),flx1)
!call vanlimter(sourcedt(i+1,j,k,2),sourcedt(i,j,k,2),sourcedt(i-1,j,k,2),flx2)
!call vanlimter(sourcedt(i+1,j,k,3),sourcedt(i,j,k,3),sourcedt(i-1,j,k,3),flx3)
!call vanlimter(sourcedt(i,j+1,k,1),sourcedt(i,j,k,1),sourcedt(i,j-1,k,1),fly1)
!call vanlimter(sourcedt(i,j+1,k,2),sourcedt(i,j,k,2),sourcedt(i,j-1,k,2),fly2)
!call vanlimter(sourcedt(i,j+1,k,3),sourcedt(i,j,k,3),sourcedt(i,j-1,k,3),fly3)
!call vanlimter(sourcedt(i,j,k+1,1),sourcedt(i,j,k,1),sourcedt(i,j,k-1,1),flz1)
!call vanlimter(sourcedt(i,j,k+1,2),sourcedt(i,j,k,2),sourcedt(i,j,k-1,2),flz2)
!call vanlimter(sourcedt(i,j,k+1,3),sourcedt(i,j,k,3),sourcedt(i,j,k-1,3),flz3)

!flx1 = dsign(1.d0,sourcedt(i,j,k,1))*dmax1(0.0d0,dmin1(dabs(sourcedt(i,j,k,1)),dsign(1.0d0,sourcedt(i,j,k,1))&
!     *sourcedt(i+1,j,k,1)))!/dabs(sourcedt(i,j,k,1))
!fly2 = dsign(1.d0,sourcedt(i,j,k,2))*dmax1(0.0d0,dmin1(dabs(sourcedt(i,j,k,2)),dsign(1.0d0,sourcedt(i,j,k,2))&
!     *sourcedt(i,j+1,k,2)))!/dabs(sourcedt(i,j,k,2))
!flz3 = dsign(1.d0,sourcedt(i,j,k,3))*dmax1(0.0d0,dmin1(dabs(sourcedt(i,j,k,3)),dsign(1.0d0,sourcedt(i,j,k,3))&
!     *sourcedt(i,j,k+1,3)))!/dabs(sourcedt(i,j,k,3))
!fllx1 = dsign(1.d0,sourcedt(i-1,j,k,1))*dmax1(0.0d0,dmin1(dabs(sourcedt(i-1,j,k,1)),dsign(1.0d0,sourcedt(i-1,j,k,1))&
!     *sourcedt(i,j,k,1)))!/dabs(sourcedt(i-1,j,k,1))
!flly2 = dsign(1.d0,sourcedt(i,j-1,k,2))*dmax1(0.0d0,dmin1(dabs(sourcedt(i,j-1,k,2)),dsign(1.0d0,sourcedt(i,j-1,k,2))&
!     *sourcedt(i,j,k,2)))!/dabs(sourcedt(i,j-1,k,2))
!fllz3 = dsign(1.d0,sourcedt(i,j,k-1,3))*dmax1(0.0d0,dmin1(dabs(sourcedt(i,j,k-1,3)),dsign(1.0d0,sourcedt(i,j,k-1,3))&
!     *sourcedt(i,j,k,3)))!/dabs(sourcedt(i,j,k-1,3))

rx = (sourcedt(i,j,k,1)-sourcedt(i-1,j,k,1)+1.d-20)/(sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1)+1.d-20)
ry = (sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2)+1.d-20)/(sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2)+1.d-20)
rz = (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3)+1.d-20)/(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3)+1.d-20)
rx2 = (sourcedt(i-1,j,k,1)-sourcedt(i-2,j,k,1)+1.d-20)/(sourcedt(i,j,k,1)-sourcedt(i-1,j,k,1)+1.d-20)
ry2 = (sourcedt(i,j-1,k,2)-sourcedt(i,j-2,k,2)+1.d-20)/(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2)+1.d-20)
rz2 = (sourcedt(i,j,k-1,3)-sourcedt(i,j,k-2,3)+1.d-20)/(sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3)+1.d-20)

flx1=dmax1(0.0d0,dmin1(1.0d0,rx))
fly2=dmax1(0.0d0,dmin1(1.0d0,ry))
flz3=dmax1(0.0d0,dmin1(1.0d0,rz))
fllx1=dmax1(0.0d0,dmin1(1.0d0,rx2))
flly2=dmax1(0.0d0,dmin1(1.0d0,ry2))
fllz3=dmax1(0.0d0,dmin1(1.0d0,rz2))

!call vanlimter(inp,in,inm,flx)
!call vanlimter(inp,in,inm,flx)
!call vanalbada(,Lnum,Phipre,Phigrad,i_sta,i_end,dmein)
!sourece(i,j,k,1)= sourcedt(i,j,k,1)-(sourcedt(i,j,k,2)+sourece2dt(i,j,k,2))-(sourcedt(i,j,k,3)+sourece2dt(i,j,k,3))&
!     +ratio*dt/deltalength*0.5d0*(sourcedt(i,j+1,k,2)-sourcedt(i,j-1,k,2) + sourcedt(i,j,k+1,3)-sourcedt(i,j,k-1,3))
!sourece(i,j,k,2)= sourcedt(i,j,k,2)-(sourcedt(i,j,k,1)+sourece2dt(i,j,k,1))-(sourcedt(i,j,k,3)+sourece2dt(i,j,k,3))&
!     +ratio*dt/deltalength*0.5d0*(sourcedt(i+1,j,k,1)-sourcedt(i-1,j-1,k,1) + sourcedt(i,j,k+1,2)-sourcedt(i,j,k-1,2))
!sourece(i,j,k,3)= sourcedt(i,j,k,3)-(sourcedt(i,j,k,2)+sourece2dt(i,j,k,2))-(sourcedt(i,j,k,1)+sourece2dt(i,j,k,1))&
!     +ratio*dt/deltalength*0.5d0*(sourcedt(i,j+1,k,2)-sourcedt(i,j-1,k,2) + sourcedt(i-1,j,k,1)-sourcedt(i-1,j,k,1))
!source(i,j,k,1)= -(sourcedt(i,j,k,2)+sourcedt2(i,j,k,2))*dddt-(sourcedt(i,j,k,3)+sourcedt2(i,j,k,3))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j-1,k,2))*fly2 + (sourcedt(i,j,k+1,3)-sourcedt(i,j,k-1,3))*flz3)
!source(i,j,k,2)= -(sourcedt(i,j,k,1)+sourcedt2(i,j,k,1))*dddt-(sourcedt(i,j,k,3)+sourcedt2(i,j,k,3))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i+1,j,k,1)-sourcedt(i-1,j,k,1))*flx1 + (sourcedt(i,j,k+1,3)-sourcedt(i,j,k-1,3))*flz3)
!source(i,j,k,3)= -(sourcedt(i,j,k,2)+sourcedt2(i,j,k,2))*dddt-(sourcedt(i,j,k,1)+sourcedt2(i,j,k,1))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j-1,k,2))*fly2 + (sourcedt(i-1,j,k,1)-sourcedt(i-1,j,k,1))*flx1)

!source!
!source(i,j,k,1)= -(sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt-(sourcedt(i,j,k,3)-sourcedt2(i,j,k,3))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2-(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
!     +(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3))*flz3 - (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3))*fllz3)
!source(i,j,k,2)= -(sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt-(sourcedt(i,j,k,3)-sourcedt2(i,j,k,3))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1-(sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 &
!     +(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3))*flz3 - (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3))*fllz3)
!source(i,j,k,3)= -(sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt-(sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt&
!     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2-(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
!     +(sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1 - (sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 )
!source!


  !source(i,j,k,1)= (sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt+(sourcedt(i,j,k,3)-sourcedt2(i,j,k,3))*dddt&
  !   +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2-(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
  !   +(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3))*flz3 - (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3))*fllz3)
  !source(i,j,k,2)= (sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt+(sourcedt(i,j,k,3)-sourcedt2(i,j,k,3))*dddt&
  !   +1.0d0/deltalength*0.5d0*((sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1-(sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 &
  !   +(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3))*flz3 - (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3))*fllz3)
  !source(i,j,k,3)= (sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt+(sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt&
  !   +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2-(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
  !   +(sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1 - (sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 )

source(i,j,k,1)= (sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt+(sourcedt(i,j,k,3)-sourcedt2(i,j,k,3))*dddt&
     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2+(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
     +(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3))*flz3 + (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3))*fllz3)
source(i,j,k,2)= (sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt+(sourcedt(i,j,k,3)-sourcedt2(i,j,k,3))*dddt&
     +1.0d0/deltalength*0.5d0*((sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1+(sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 &
     +(sourcedt(i,j,k+1,3)-sourcedt(i,j,k,3))*flz3 + (sourcedt(i,j,k,3)-sourcedt(i,j,k-1,3))*fllz3)
source(i,j,k,3)= (sourcedt(i,j,k,2)-sourcedt2(i,j,k,2))*dddt+(sourcedt(i,j,k,1)-sourcedt2(i,j,k,1))*dddt&
     +1.0d0/deltalength*0.5d0*((sourcedt(i,j+1,k,2)-sourcedt(i,j,k,2))*fly2+(sourcedt(i,j,k,2)-sourcedt(i,j-1,k,2))*flly2 &
     +(sourcedt(i+1,j,k,1)-sourcedt(i,j,k,1))*flx1 + (sourcedt(i,j,k,1)-sourcedt(i+1,j,k,1))*fllx1 )
!write(5,*) source(i,j,k,1),source(i,j,k,2),source(i,j,k,3),sourcedt(i,j,k,1),sourcedt2(i,j,k,1),sourcedt(i,j,k,2),&
!     sourcedt2(i,j,k,2),sourcedt(i,j,k,3),sourcedt2(i,j,k,3),flx1,fly2,flz3,fllx1,flly2,fllz3
end do
end do
end do

!close(5)
c=c+1

bphixl(:,:,:,:)=0.d0
bphixr(:,:,:,:)=0.d0
bphiyl(:,:,:,:)=0.d0
bphiyr(:,:,:,:)=0.d0
bphizl(:,:,:,:)=0.d0
bphizr(:,:,:,:)=0.d0

bstepxl(:,:,:,:)=0.d0
bstepxr(:,:,:,:)=0.d0
bstepyl(:,:,:,:)=0.d0
bstepyr(:,:,:,:)=0.d0
bstepzl(:,:,:,:)=0.d0
bstepzr(:,:,:,:)=0.d0

do k = -1,ndz
do j = -1,ndy
do i = -1,ndx
Phi(i,j,k)=Phicgp(i,j,k,1)+Phicgp(i,j,k,2)+Phicgp(i,j,k,3)
end do
end do
end do

!iwx=1;iwy=1;iwz=1
!call BCgrv(100,0)

do k = -1,Ncellz+2
do j = -1,Ncelly+2

bstepxl(j,k,1,1) = (3.0d0*bphil(j,k,1)-4.0d0*bphil(j,k,0)+bphil(j,k,-1))*0.5d0/dxx
bstepxl(j,k,0,1) = (-bphil(j,k,-1)+bphil(j,k,1))*0.5d0/dxx
bstepxl(j,k,-1,1) = -(3.0d0*bphil(j,k,-1)-4.0d0*bphil(j,k,0)+bphil(j,k,1))*0.5d0/dxx
bstepxr(j,k,Ncellx+2,1) = (3.0d0*bphir(j,k,Ncellx+2)-4.0d0*bphir(j,k,Ncellx+1)+bphir(j,k,Ncellx))*0.5d0/dxx
bstepxr(j,k,Ncellx+1,1) = (-bphir(j,k,Ncellx)+bphir(j,k,Ncellx+2))*0.5d0/dxx
bstepxr(j,k,Ncellx,1) = -(3.0d0*bphir(j,k,Ncellx)-4.0d0*bphir(j,k,Ncellx+1)+bphir(j,k,Ncellx+2))*0.5d0/dxx

!bstepxl(j,k,1,1) = 3.0d0*(3.0d0*bphil(j,k,1)-4.0d0*bphil(j,k,0)+bphil(j,k,-1))*0.5d0/dxx
!bstepxl(j,k,0,1) = 3.0d0*(-bphil(j,k,-1)+bphil(j,k,1))*0.5d0/dxx
!bstepxl(j,k,-1,1) = -3.0d0*(3.0d0*bphil(j,k,-1)-4.0d0*bphil(j,k,0)+bphil(j,k,1))*0.5d0/dxx
!bstepxr(j,k,Ncellx+2,1) = 3.0d0*(3.0d0*bphir(j,k,Ncellx+2)-4.0d0*bphir(j,k,Ncellx+1)+bphir(j,k,Ncellx))*0.5d0/dxx
!bstepxr(j,k,Ncellx+1,1) = 3.0d0*(-bphir(j,k,Ncellx)+bphir(j,k,Ncellx+2))*0.5d0/dxx
!bstepxr(j,k,Ncellx,1) = -3.0d0*(3.0d0*bphir(j,k,Ncellx)-4.0d0*bphir(j,k,Ncellx+1)+bphir(j,k,Ncellx+2))*0.5d0/dxx
!write(3,*) bstepxl(j,k,1,1),bstepxl(j,k,0,1),bstepxl(j,k,-1,1),bstepxr(j,k,Ncellx+2,1),bstepxr(j,k,Ncellx+1,1),bstepxr(j,k,Ncellx,1),&
!     bstepxl(j,k,1,1),bstepxl(j,k,0,1),bstepxl(j,k,-1,1),bstepxr(j,k,Ncellx+2,1),bstepxr(j,k,Ncellx+1,1),bstepxr(j,k,Ncellx,1),'1'
end do
end do

!close(3)

Phi1step(     0,      -1,:,1)= -(2.d0*Phi1step(1,-1,:,1)-Phi1step(2,-1,:,1))
Phi1step(Ncellx,      -1,:,1)= -(Phi1step(Ncellx-1,-1,:,1)-2.d0*Phi1step(Ncellx,-1,:,1))
Phi1step(     0,       0,:,1)= -(2.d0*Phi1step(1,0,:,1)-Phi1step(2,0,:,1))
Phi1step(Ncellx,       0,:,1)= -(Phi1step(Ncellx-1,0,:,1)-2.d0*Phi1step(Ncellx,0,:,1))
Phi1step(     0,Ncelly+2,:,1)= -(2.d0*Phi1step(1,Ncelly+2,:,1)-Phi1step(2,Ncelly+2,:,1))
Phi1step(Ncellx,Ncelly+2,:,1)= -(Phi1step(Ncellx-1,Ncelly+2,:,1)-2.d0*Phi1step(Ncellx,Ncelly+2,:,1))
Phi1step(     0,Ncelly+1,:,1)= -(2.d0*Phi1step(1,Ncelly+1,:,1)-Phi1step(2,Ncelly+1,:,1))
Phi1step(Ncellx,Ncelly+1,:,1)= -(Phi1step(Ncellx-1,Ncelly+1,:,1)-2.d0*Phi1step(Ncellx,Ncelly+1,:,1))

Phi1step(:,      -1,     0,3)= -(2.d0*Phi1step(:,-1,1,3)-Phi1step(:,-1,2,3))
Phi1step(:,      -1,Ncellz,3)= -(Phi1step(:,-1,Ncellz-1,3)-2.d0*Phi1step(:,-1,Ncellz,3))
Phi1step(:,       0,     0,3)= -(2.d0*Phi1step(:,0,1,3)-Phi1step(:,0,2,3))
Phi1step(:,       0,Ncellz,3)= -(Phi1step(:,0,Ncellz-1,3)-2.d0*Phi1step(:,0,Ncellz,3))
Phi1step(:,Ncelly+2,     0,3)= -(2.d0*Phi1step(:,Ncelly+2,1,3)-Phi1step(:,Ncelly+2,2,3))
Phi1step(:,Ncelly+2,Ncellz,3)= -(Phi1step(:,Ncelly+2,Ncellz-1,3)-2.d0*Phi1step(:,Ncelly+2,Ncellz,3))
Phi1step(:,Ncelly+1,     0,3)= -(2.d0*Phi1step(:,Ncelly+1,1,3)-Phi1step(:,Ncelly+1,2,3))
Phi1step(:,Ncelly+1,Ncellz,3)= -(Phi1step(:,Ncelly+1,Ncellz-1,3)-2.d0*Phi1step(:,Ncelly+1,Ncellz,3))

dxx=-deltalength

do k = -1,Ncellx+2
do j = -1,Ncellz+2

!bstepyl(j,k,1,2) = (3.0d0*Phi(k,1,j)-4.0d0*Phi(k,0,j)+Phi(k,-1,j))*0.5d0/dxx
!bstepyl(j,k,0,2) = (-Phi(k,-1,j)+Phi(k,1,j))*0.5d0/dxx
!bstepyl(j,k,-1,2) = -(3.0d0*Phi(k,-1,j)-4.0d0*Phi(k,0,j)+Phi(k,1,j))*0.5d0/dxx
!bstepyr(j,k,Ncelly+2,2) = (3.0d0*Phi(k,Ncelly+2,j)-4.0d0*Phi(k,Ncelly+1,j)+Phi(k,Ncelly,j))*0.5d0/dxx
!bstepyr(j,k,Ncelly+1,2) = (-Phi(k,Ncelly,j)+Phi(k,Ncelly+2,j))*0.5d0/dxx
!bstepyr(j,k,Ncelly,2) = -(3.0d0*Phi(k,Ncelly,j)-4.0d0*Phi(k,Ncelly+1,j)+Phi(k,Ncelly+2,j))*0.5d0/dxx

!bstepyl(j,k,1,2) = (3.0d0*Phicgp(k,1,j,2)-4.0d0*Phicgp(k,0,j,2)+Phicgp(k,-1,j,2))*0.5d0/dxx
!bstepyl(j,k,0,2) = (-Phicgp(k,-1,j,2)+Phicgp(k,1,j,2))*0.5d0/dxx
!bstepyl(j,k,-1,2) = -(3.0d0*Phicgp(k,-1,j,2)-4.0d0*Phicgp(k,0,j,2)+Phicgp(k,1,j,2))*0.5d0/dxx
!bstepyr(j,k,Ncelly+2,2) = (3.0d0*Phicgp(k,Ncelly+2,j,2)-4.0d0*Phicgp(k,Ncelly+1,j,2)+Phicgp(k,Ncelly,j,2))*0.5d0/dxx
!bstepyr(j,k,Ncelly+1,2) = (-Phicgp(k,Ncelly,j,2)+Phicgp(k,Ncelly+2,j,2))*0.5d0/dxx
!bstepyr(j,k,Ncelly,2) = -(3.0d0*Phicgp(k,Ncelly,j,2)-4.0d0*Phicgp(k,Ncelly+1,j,2)+Phicgp(k,Ncelly+2,j,2))*0.5d0/dxx
!bc21
!bc22

!bstepyl(j,k, 1,2) = G4pi*U(k,1,j,1)-(Phi1step(k+1,1,j,1)-Phi1step(k-1,1,j,1))*0.5d0/dxx-(Phi1step(k,1,j+1,3)-Phi1step(k,1,j-1,3))*0.5d0/dxx
!bstepyl(j,k,-1,2) = Phi1step(k,1,j,2)-(G4pi*U(k,-1,j,1)-(Phi1step(k+1,-1,j,1)-Phi1step(k-1,-1,j,1))*0.5d0/dxx&
!        -(Phi1step(k,-1,j+1,3)-Phi1step(k,-1,j-1,3))*0.5d0/dxx)*2.d0*dxx
!bstepyl(j,k, 0,2) = (3.d0*Phi1step(k,1,j,2)+bstepyl(j,k,-1,2)-(G4pi*U(k,0,j,1)-(Phi1step(k+1,0,j,1)-Phi1step(k-1,0,j,1))*0.5d0/dxx&
!        -(Phi1step(k,0,j+1,3)-Phi1step(k,0,j-1,3))*0.5d0/dxx)*2.d0*dxx)*0.25d0
!bstepyl(j,k,-1,2) = G4pi*U(k,-1,j,1)-(Phi1step(k+1,-1,j,1)-Phi1step(k-1,-1,j,1))*0.5d0/dxx-(Phi1step(k,-1,j+1,3)-Phi1step(k,-1,j-1,3))*0.5d0/dxx
!bstepyr(j,k,Ncelly+2,2) = -Phi1step(k,Ncelly,j,2)+(G4pi*U(k,Ncelly+2,j,1)-(Phi1step(k+1,Ncelly+2,j,1)-Phi1step(k-1,Ncelly+2,j,1))*0.5d0/dxx-&
!     (Phi1step(k,Ncelly+2,j+1,3)-Phi1step(k,Ncelly+2,j-1,3))*0.5d0/dxx)*2.d0*dxx
!bstepyr(j,k,Ncelly+1,2) =(Phi1step(k,Ncelly,j,2)+3.d0*bstepyr(j,k,Ncelly+2,2) -(G4pi*U(k,Ncelly+1,j,1)-(Phi1step(k+1,Ncelly+1,j,1)-&
!     Phi1step(k-1,Ncelly+1,j,1))*0.5d0/dxx-&
!     (Phi1step(k,Ncelly+1,j+1,3)-Phi1step(k,Ncelly+1,j-1,3))*0.5d0/dxx)*2.d0*dxx)*0.25d0
!bstepyr(j,k,Ncelly  ,2) = G4pi*U(k,1,j,1)-(Phi1step(k+1,1,j,1)-Phi1step(k-1,1,j,1))*0.5d0/dxx-(Phi1step(k,1,j+1,3)-Phi1step(k,1,j-1,3))*0.5d0/dxx


bstepyl(j,k,-1,2) = Phi1step(k,1,j,2)-(G4pi*U(k,-1,j,1)-(Phi1step(k+1,-1,j,1)-Phi1step(k-1,-1,j,1))*0.5d0/dxx&
        -(Phi1step(k,-1,j+1,3)-Phi1step(k,-1,j-1,3))*0.5d0/dxx)*2.d0*dxx
bstepyl(j,k, 0,2) = (3.d0*Phi1step(k,1,j,2)+bstepyl(j,k,-1,2)-(G4pi*U(k,0,j,1)-(Phi1step(k+1,0,j,1)-Phi1step(k-1,0,j,1))*0.5d0/dxx&
        -(Phi1step(k,0,j+1,3)-Phi1step(k,0,j-1,3))*0.5d0/dxx)*2.d0*dxx)*0.25d0
bstepyr(j,k,Ncelly+2,2) = -Phi1step(k,Ncelly,j,2)+(G4pi*U(k,Ncelly+2,j,1)-(Phi1step(k+1,Ncelly+2,j,1)-Phi1step(k-1,Ncelly+2,j,1))*0.5d0/dxx-&
     (Phi1step(k,Ncelly+2,j+1,3)-Phi1step(k,Ncelly+2,j-1,3))*0.5d0/dxx)*2.d0*dxx
bstepyr(j,k,Ncelly+1,2) =(Phi1step(k,Ncelly,j,2)+3.d0*bstepyr(j,k,Ncelly+2,2) -(G4pi*U(k,Ncelly+1,j,1)-(Phi1step(k+1,Ncelly+1,j,1)-&
     Phi1step(k-1,Ncelly+1,j,1))*0.5d0/dxx-&
     (Phi1step(k,Ncelly+1,j+1,3)-Phi1step(k,Ncelly+1,j-1,3))*0.5d0/dxx)*2.d0*dxx)*0.25d0

!write(3,*) bstepyl(j,k,1,2),bstepyl(j,k,0,2),bstepyl(j,k,-1,2),bstepyr(j,k,Ncellx+2,2),bstepyr(j,k,Ncellx+1,2),bstepyr(j,k,Ncellx,2),&
!     bstepyl(j,k,1,2),bstepyl(j,k,0,2),bstepyl(j,k,-1,2),bstepyr(j,k,Ncellx+2,2),bstepyr(j,k,Ncellx+1,2),bstepyr(j,k,Ncellx,2),'2'
end do
end do
!bstepyl(:,:,:,2)=0.0d0
!bstepyr(:,:,:,2)=0.0d0

!Phi1step(0,:,-1,1)=-(2.d0*Phi1step(1,:,-1,1)-Phi1step(2,:,-1,1))
!Phi1step(Ncellx,:,-1,1)=-(Phi1step(Ncellx-1,:,-1,1)-2.d0*Phi1step(Ncellx,:,-1,1))
!Phi1step(0,:,0,1)=-(2.d0*Phi1step(1,:,0,1)-Phi1step(2,:,0,1))
!Phi1step(Ncellx,:,0,1)=-(Phi1step(Ncellx-1,:,0,1)-2.d0*Phi1step(Ncellx,:,0,1))

!Phi1step(:,-1,     0,2)=2.d0*Phi1step(:,-1,1,3)-Phi1step(:,-1,2,3)
!Phi1step(:,-1,Ncellz,2)=Phi1step(:,-1,Ncellz-1,3)-2.d0*Phi1step(:,-1,Ncellz,3)
!Phi1step(:, 0,     0,2)=2.d0*Phi1step(:,0,1,3)-Phi1step(:,0,2,3)
!Phi1step(:, 0,Ncellz,2)=Phi1step(:,0,Ncellz-1,3)-2.d0*Phi1step(:,0,Ncellz,3)

Phi1step(     0,:,      -1,1)= -(2.d0*Phi1step(1,:,-1,1)-Phi1step(2,:,-1,1))
Phi1step(Ncellx,:,      -1,1)= -(     Phi1step(Ncellx-1,:,-1,1)-2.d0*Phi1step(Ncellx,:,-1,1))
Phi1step(     0,:,       0,1)= -(2.d0*Phi1step(1,:,0,1)-Phi1step(2,:,0,1))
Phi1step(Ncellx,:,       0,1)= -(     Phi1step(Ncellx-1,:,0,1)-2.d0*Phi1step(Ncellx,:,0,1))
Phi1step(     0,:,Ncellz+2,1)= -(2.d0*Phi1step(1,:,Ncellz+2,1)-Phi1step(2,:,Ncellz+2,1))
Phi1step(Ncellx,:,Ncellz+2,1)= -(     Phi1step(Ncellx-1,:,Ncellz+2,1)-2.d0*Phi1step(Ncellx,:,Ncellz+2,1))
Phi1step(     0,:,Ncellz+1,1)= -(2.d0*Phi1step(1,:,Ncellz+1,1)-Phi1step(2,:,Ncellz+1,1))
Phi1step(Ncellx,:,Ncellz+1,1)= -(     Phi1step(Ncellx-1,:,Ncellz+1,1)-2.d0*Phi1step(Ncellx,:,Ncellz+1,1))

Phi1step(:,     0,      -1,2)= -(2.d0*Phi1step(:,1,-1,2)-Phi1step(:,2,-1,2))
Phi1step(:,Ncelly,      -1,2)= -(Phi1step(:,Ncelly-1,-1,2)-2.d0*Phi1step(:,Ncelly,-1,2))
Phi1step(:,     0,       0,2)= -(2.d0*Phi1step(:,1,0,2)-Phi1step(:,2,0,2))
Phi1step(:,Ncelly,       0,2)= -(Phi1step(:,Ncelly-1,0,2)-2.d0*Phi1step(:,Ncelly,0,2))
Phi1step(:,     0,Ncellz+2,2)= -(2.d0*Phi1step(:,1,Ncellz+2,2)-Phi1step(:,2,Ncellz+2,2))
Phi1step(:,Ncelly,Ncellz+2,2)= -(Phi1step(:,Ncelly-1,Ncellz+2,2)-2.d0*Phi1step(:,Ncelly,Ncellz+2,2))
Phi1step(:,     0,Ncellz+1,2)= -(2.d0*Phi1step(:,1,Ncellz+1,2)-Phi1step(:,2,Ncellz+1,2))
Phi1step(:,Ncelly,Ncellz+1,2)= -(Phi1step(:,Ncelly-1,Ncellz+1,2)-2.d0*Phi1step(:,Ncelly,Ncellz+1,2))

do k = -1,Ncelly+2
do j = -1,Ncellx+2

!bstepzl(j,k,1,3) = (3.0d0*Phi(j,k,1)-4.0d0*Phi(j,k,0)+Phi(j,k,-1))*0.5d0/dxx
!bstepzl(j,k,0,3) = (-Phi(j,k,-1)+Phi(j,k,1))*0.5d0/dxx
!bstepzl(j,k,-1,3) = -(3.0d0*Phi(j,k,-1)-4.0d0*Phi(j,k,0)+Phi(j,k,1))*0.5d0/dxx
!bstepzr(j,k,Ncellz+2,3) = (3.0d0*Phi(j,k,Ncellz+2)-4.0d0*Phi(j,k,Ncellz+1)+Phi(j,k,Ncellz))*0.5d0/dxx
!bstepzr(j,k,Ncellz+1,3) = (-Phi(j,k,Ncellz)+Phi(j,k,Ncellz+2))*0.5d0/dxx
!bstepzr(j,k,Ncellz,3) = -(3.0d0*Phi(j,k,Ncellz)-4.0d0*Phi(j,k,Ncellx+1)+Phi(j,k,Ncellz+2))*0.5d0/dxx

!bstepzl(j,k,1,3) = (3.0d0*Phicgp(j,k,1,3)-4.0d0*Phicgp(j,k,0,3)+Phicgp(j,k,-1,3))*0.5d0/dxx
!bstepzl(j,k,0,3) = (-Phicgp(j,k,-1,3)+Phicgp(j,k,1,3))*0.5d0/dxx
!bstepzl(j,k,-1,3) = -(3.0d0*Phicgp(j,k,-1,3)-4.0d0*Phicgp(j,k,0,3)+Phicgp(j,k,1,3))*0.5d0/dxx
!bstepzr(j,k,Ncellz+2,3) = (3.0d0*Phicgp(j,k,Ncellz+2,3)-4.0d0*Phicgp(j,k,Ncellz+1,3)+Phicgp(j,k,Ncellz,3))*0.5d0/dxx
!bstepzr(j,k,Ncellz+1,3) = (-Phicgp(j,k,Ncellz,3)+Phicgp(j,k,Ncellz+2,3))*0.5d0/dxx
!bstepzr(j,k,Ncellz,3) = -(3.0d0*Phicgp(j,k,Ncellz,3)-4.0d0*Phicgp(j,k,Ncellx+1,3)+Phicgp(j,k,Ncellz+2,3))*0.5d0/dxx

bstepzl(j,k,-1,3) = Phi1step(j,k,1,3)-(G4pi*U(j,k,-1,1)-(Phi1step(j+1,k,-1,1)-Phi1step(j-1,k,-1,1))*0.5d0/dxx&
        -(Phi1step(j,k+1,-1,2)-Phi1step(j,k-1,-1,2))*0.5d0/dxx)*2.d0*dxx
bstepzl(j,k, 0,3) = (3.d0*Phi1step(j,k,1,3)+bstepzl(j,k,-1,3)-(G4pi*U(j,k,0,1)-(Phi1step(j+1,k,0,1)-Phi1step(j-1,k,0,1))*0.5d0/dxx&
        -(Phi1step(j,k+1,0,2)-Phi1step(j,k-1,0,2))*0.5d0/dxx)*2.d0*dxx)*0.25d0

bstepzr(j,k,Ncellz+2,3) = -Phi1step(j,k,Ncellz,3)+(G4pi*U(j,k,Ncellz+2,1)-(Phi1step(j+1,k,Ncellz+2,1)-Phi1step(j-1,k,Ncellz+2,1))*0.5d0/dxx-&
     (Phi1step(j,k+1,Ncellz+2,2)-Phi1step(j,k-1,Ncellz+2,2))*0.5d0/dxx)*2.d0*dxx
bstepzr(j,k,Ncellz+1,3) =(Phi1step(j,k,Ncellz,3)+3.d0*bstepzr(j,k,Ncellz+2,3) -(G4pi*U(j,k,Ncellz+1,1)-(Phi1step(j+1,k,Ncellz+1,1)-&
     Phi1step(j-1,k,Ncellz+1,1))*0.5d0/dxx-&
     (Phi1step(j,k+1,Ncellz+1,2)-Phi1step(j,k-1,Ncellz+1,2))*0.5d0/dxx)*2.d0*dxx)*0.25d0

!write(3,*) bstepzl(j,k,1,3),bstepzl(j,k,0,3),bstepzl(j,k,-1,3),bstepzr(j,k,Ncellx+2,3),bstepzr(j,k,Ncellx+1,3),bstepzr(j,k,Ncellx,3),&
!     bstepzl(j,k,1,3),bstepzl(j,k,0,3),bstepzl(j,k,-1,3),bstepzr(j,k,Ncellx+2,3),bstepzr(j,k,Ncellx+1,3),bstepzr(j,k,Ncellx,3),'3'
end do
end do
!close(3)
!bstepzl(:,:,:,3)=0.0d0
!bstepzr(:,:,:,3)=0.0d0

do k = -1,Ncellz+2
do j = -1,Ncelly+2

!bphixl(j,k,-1,2) = Phicgp(0,j,k,2)
bphixl(j,k,-1,2) = Phicgp(1,j,k,2)
bphixl(j,k, 0,2) = Phicgp(1,j,k,2)
!bphixr(j,k,Ncellx+2,2) = Phicgp(Ncellx+1,j,k,2)
bphixr(j,k,Ncellx+2,2) = Phicgp(Ncellx,j,k,2)
bphixr(j,k,Ncellx+1,2) = Phicgp(Ncellx  ,j,k,2)

bphixl(j,k,-1,3) = Phicgp(0,j,k,3)
bphixl(j,k, 0,3) = Phicgp(1,j,k,3)
bphixr(j,k,Ncellx+2,3) = Phicgp(Ncellx+1,j,k,3)
bphixr(j,k,Ncellx+1,3) = Phicgp(Ncellx  ,j,k,3)

!bphixl(j,k,-1,1) = bphil(j,k,-1)-bphixl(j,k,-1,2)-bphixl(j,k,-1,3)
!bphixl(j,k, 0,1) = bphil(j,k, 0)-bphixl(j,k, 0,2)-bphixl(j,k, 0,3)
!bphixr(j,k,Ncellx+2,1) = bphir(j,k,Ncellx+2)-bphixr(j,k,Ncellx+2,2)-bphixr(j,k,Ncellx+2,3)
!bphixr(j,k,Ncellx+1,1) = bphir(j,k,Ncellx+1)-bphixr(j,k,Ncellx+1,2)-bphixr(j,k,Ncellx+1,3)

bphixl(j,k,-1,1) = bphil(j,k,-1)-Phicgp(-1,j,k,2)-Phicgp(-1,j,k,3)
bphixl(j,k, 0,1) = bphil(j,k, 0)-Phicgp( 0,j,k,2)-Phicgp( 0,j,k,3)
bphixr(j,k,Ncellx+2,1) = bphir(j,k,Ncellx+2)-Phicgp(Ncellx+2,j,k,2)-Phicgp(Ncellx+2,j,k,3)
bphixr(j,k,Ncellx+1,1) = bphir(j,k,Ncellx+1)-Phicgp(Ncellx+1,j,k,2)-Phicgp(Ncellx+1,j,k,3)

!bphixl(j,k,-1,1) = bphil(j,k,-1)+(-Phicgp(-1,j,k,2)-Phicgp(-1,j,k,3))/3.0d0
!bphixl(j,k, 0,1) = bphil(j,k, 0)+(-Phicgp( 0,j,k,2)-Phicgp( 0,j,k,3))/3.0d0
!bphixr(j,k,Ncellx+2,1) = bphir(j,k,Ncellx+2)+(-Phicgp(Ncellx+2,j,k,2)-Phicgp(Ncellx+2,j,k,3))/3.0d0
!bphixr(j,k,Ncellx+1,1) = bphir(j,k,Ncellx+1)+(-Phicgp(Ncellx+1,j,k,2)-Phicgp(Ncellx+1,j,k,3))/3.0d0

!bphixl(j,k,-1,1) = 3.0d0*bphil(j,k,-1)-Phicgp(-1,j,k,2)-Phicgp(-1,j,k,3)
!bphixl(j,k, 0,1) = 3.0d0*bphil(j,k, 0)-Phicgp( 0,j,k,2)-Phicgp( 0,j,k,3)
!bphixr(j,k,Ncellx+2,1) = 3.0d0*bphir(j,k,Ncellx+2)-Phicgp(Ncellx+2,j,k,2)-Phicgp(Ncellx+2,j,k,3)
!bphixr(j,k,Ncellx+1,1) = 3.0d0*bphir(j,k,Ncellx+1)-Phicgp(Ncellx+1,j,k,2)-Phicgp(Ncellx+1,j,k,3)

!write(4,*) bphixl(j,k,-1,1),bphixl(j,k,-1,2),bphixl(j,k,-1,3),bphixr(j,k,Ncellx+1,1),bphixl(j,k,Ncellx+1,2),bphixl(j,k,Ncellx+1,3)
end do
end do

do k = -1,Ncellx+2
do j = -1,Ncellz+2

bphiyl(j,k,-1,1) = Phicgp(k,0,j,1)
bphiyl(j,k, 0,1) = Phicgp(k,1,j,1)
bphiyr(j,k,Ncelly+2,1) = Phicgp(k,Ncelly+1,j,1)
bphiyr(j,k,Ncelly+1,1) = Phicgp(k,Ncelly  ,j,1)

bphiyl(j,k,-1,3) = Phicgp(k,0,j,3)
bphiyl(j,k, 0,3) = Phicgp(k,1,j,3)
bphiyr(j,k,Ncelly+2,3) = Phicgp(k,Ncelly+1,j,3)
bphiyr(j,k,Ncelly+1,3) = Phicgp(k,Ncelly  ,j,3)

bphiyl(j,k,-1,2) = Phi(k,-1,j)-bphiyl(j,k,-1,1)-bphiyl(j,k,-1,3)
bphiyl(j,k, 0,2) = Phi(k, 0,j)-bphiyl(j,k, 0,1)-bphiyl(j,k, 0,3)
bphiyr(j,k,Ncelly+2,2) = Phi(k,Ncelly+2,j)-bphiyr(j,k,Ncelly+2,1)-bphiyr(j,k,Ncelly+2,3)!-bphiyr(k,Ncelly+2,j,1)-bphiyr(k,Ncelly+2,j,3)
bphiyr(j,k,Ncelly+1,2) = Phi(k,Ncelly+1,j)-bphiyr(j,k,Ncelly+1,1)-bphiyr(j,k,Ncelly+1,3)!-bphiyr(k,Ncelly+1,j,1)-bphiyr(k,Ncelly+1,j,3)

!write(4,*) bphiyl(j,k,-1,1),bphiyl(j,k,-1,2),bphiyl(j,k,-1,3),bphiyr(j,k,Ncelly+1,1),bphiyl(j,k,Ncelly+1,2),bphiyl(j,k,Ncelly+1,3)
end do
end do

do k = -1,Ncelly+2
do j = -1,Ncellx+2

bphizl(j,k,-1,1) = Phicgp(j,k,0,1)
bphizl(j,k, 0,1) = Phicgp(j,k,1,1)
bphizr(j,k,Ncellz+2,1) = Phicgp(j,k,Ncellz+1,1)
bphizr(j,k,Ncellz+1,1) = Phicgp(j,k,Ncellz  ,1)

bphizl(j,k,-1,2) = Phicgp(j,k,0,2)
bphizl(j,k, 0,2) = Phicgp(j,k,1,2)
bphizr(j,k,Ncellz+2,2) = Phicgp(j,k,Ncellz+1,2)
bphizr(j,k,Ncellz+1,2) = Phicgp(j,k,Ncellz  ,2)

bphizl(j,k,-1,3) = Phi(j,k,-1)-bphizl(j,k,-1,1)-bphizl(j,k,-1,2)
bphizl(j,k, 0,3) = Phi(j,k, 0)-bphizl(j,k, 0,1)-bphizl(j,k, 0,2)
bphizr(j,k,Ncellz+2,3) = Phi(j,k,Ncellz+2)-bphizr(j,k,Ncellz+2,1)-bphizr(j,k,Ncellz+2,2)!-bphizr(k,Ncellz+2,j,1)-bphizr(k,Ncellz+2,j,2)
bphizr(j,k,Ncellz+1,3) = Phi(j,k,Ncellz+1)-bphizr(j,k,Ncellz+1,1)-bphizr(j,k,Ncellz+1,2)!-bphizr(k,Ncellz+1,j,1)-bphizr(k,Ncellz+1,j,2)

!write(4,*) bphizl(j,k,-1,1),bphizl(j,k,-1,2),bphizl(j,k,-1,3),bphizr(j,k,Ncellz+1,1),bphizl(j,k,Ncellz+1,2),bphizl(j,k,Ncellz+1,3)
end do
end do
!close(4)

end subroutine pbstep


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
  data(jy,kz,i) = U(i,j,k,1)
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
  zp2 = (x(kz)-0.5d0*dzz)-LLr
!  zp2 = Lbox - zp1
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

write(fn,'(i3.3)') NRANK
write(lRANK,'(i1.1)') pls+2
open(12,file=dir//'bcsave'//lRANK//fn//'.DAT')
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

  bphil(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  bphir(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
  !bphi2l(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  !bphi2r(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
  write(12,*) bphil(j,k,1-abs(pls)),bphir(j,k,Ncellx+abs(pls))!,bphi2l(j,k,1-abs(pls)), bphi2r(j,k,Ncellx+abs(pls))
end do
end do
close(12)

DEALLOCATE(data,speq)
DEALLOCATE(dat1,spe1,dat2,spe2)
DEALLOCATE(bcspel1,bcspel2,bcl1,bcr2)
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
