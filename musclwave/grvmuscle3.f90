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
  double precision dt_mpi_gr(0:NPE-1),dt_gat_gr(0:NPE-1),maxcs,tcool,cgtime,sourcedt
  double precision :: ave1,ave1pre,ave2(0:NPE-1),ave,avepre,ave2_gather(0:NPE-1) , eps=1.0d-3
  !double precision , dimension(:,:,:) , allocatable :: stbPhi
  !double precision , dimension(-1:Ncellx+2,-1:Ncelly,-1:Ncellz) :: Phipregrad,Phipregraddum
  !**************** INITIALIZEATION **************
  if(mode==0) then
     Phicgp(:,:,:)=0.0d0
     Phicgp(:,:,:)=0.0d0
     Phi1step(:,:,:)=0.0d0
     Phi2step(:,:,:)=0.0d0
  end if
  !**************** INITIALIZEATION **************



  !****************read INITIAL CONDITION**************
  if(mode==1) then
     WRITE(NPENUM,'(I3.3)') NRANK
     open(unit=8,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgp'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=18,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgm'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI1step'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=38,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI2step'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(8) (Phicgp(i,j,k),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(8)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(18) (Phicgm(i,j,k),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(18)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(28) (Phi1step(i,j,k),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(28)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           read(38) (Phi2step(i,j,k),i=-1,Ncellx+2)
           !enddo
        end do
     end do
     close(38)
  end if
  !****************read INITIAL CONDITION**************

  !****************GRAVITY SOLVER*****************

  if(mode==2) then
     N_MPI(20)=1; N_MPI(1)=1
     iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)
     !---debug---
     !call  SELFGRAVWAVE(0.0d0,4)
     write(*,*) '------pb1-------' ,Nrank
     Call PB(0)
     Call PB(-1)
     Call PB(-2)
     call pbstep()
     write(*,*) '------pb2-------' ,Nrank
     call slvmuscle(dt)
     write(*,*) '------slv-------' ,Nrank
  end if

  !****************GRAVITY SOLVER*****************




  !**********acceraration because of gravity******
  if(mode==3) then !acceraration because of gravity
     iwx = 1; iwy = 1; iwz = 1
     call BCgrv(101)
     call BCgrv(102)
     dxi = 1.d0/(12.d0*dx(0))
     do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
        !U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
        !U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
        !U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0
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
     !open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIDTINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=8,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgp'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=18,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgm'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI1step'//NPENUM//countcha//'.DAT')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI2step'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
     call BCgrv(101)
     call BCgrv(102)
     do k = -1, Ncellz+2
        do j = -1, Ncelly+2
           !do i = -1, Ncellx+2
           write(28,*) (sngl(Phicgp(i,j,k)),sngl(Phicgm(i,j,k)),sngl(Phi1step(i,j,k)),sngl(Phi2step(i,j,k)),&
                sngl((Phicgp(i,j,k)+Phicgm(i,j,k))*0.5d0),i=-1,Ncellx+2)
           !enddo
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
   cgtime = dx(1)/cg * 0.2d0 !3D-isotropic
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

do k=1,Ncellz
   do j=1,Ncelly
      do i=1,Ncellx
         if(Phidt(i,j,k) .ne. 0.0d0) then
         ave1pre=ave1
         ave1 = dabs((Phi(i,j,k)-Phidt(i,j,k))/Phidt(i,j,k))! + ave1
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



!***************SABILITY-exa**************
if(mode==9) then
   sourcedt = dt
   !cgtime = deltalength/cg * CFL
   !call STBLphi(sourcedt,Phipregrad)
   dt = sourcedt
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
  double precision dt
  integer :: i=0,n,m,l
  double precision rho(-1:ndx,-1:ndy,-1:ndz)
  integer ifEVOgrv,ifEVOgrv2

  ifEVOgrv=1+mod(i,6)
  ifEVOgrv2=1+mod(i,6)

  do l=-1,ndz
  do m=-1,ndy
  do n=-1,ndx
     rho(n,m,l) = U(n,m,l,1)
  end do;end do;end do
  !write(*,*) '--bc--111--',Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1),&
  !     rho(1,1,1)
  !write(*,*) ifEVOgrv,ifEVOgrv2
  !call cllsub(2,dt)
  !call cllsub(3,dt)
  iwx=1; iwy=1; iwz=1
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  !write(*,*) 'BC------1'
  call BCgrv(102)
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !write(*,*) '--bc--112--',Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1)
  !write(*,*) 'BC------2'
  call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !write(*,*) 'BC------2.2'
  call muslcslv1D(Phi2step,rho,dt*0.5d0,3,2)
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  ! write(*,*) 'BC------3'
  !iwx=1; iwy=1; iwz=1
  !call BCgrv(102)
  !call cllsub(3,dt)
  !write(*,*) '--bc--113--',Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1)
  call cllsub(4,dt,ifEVOgrv,ifEVOgrv2)
  !call cllsub(3,dt)
  !call cllsub(1,dt)
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  ! write(*,*) 'BC------4'
  iwx=1; iwy=1; iwz=1
  call BCgrv(101)
  call BCgrv(102)
  call muslcslv1D(Phicgp,Phi1step,dt,4,2)
  call muslcslv1D(Phicgm,Phi2step,dt,4,2)
  ! write(*,*) 'BC------5'
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  !call cllsub(1,dt)
  !write(*,*) '--bc--114--',Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1)
  call cllsub(5,dt,dtifEVOgrv,ifEVOgrv2)
  !call cllsub(3,dt)
  ! write(*,*) 'BC------6'
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  iwx=1; iwy=1; iwz=1
  call BCgrv(102)
  call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
  call muslcslv1D(Phi2step,rho,0.5d0*dt,3,2)
  ! write(*,*) 'BC------7'
  !---debug---
  !call  SELFGRAVWAVE(0.0d0,4)
  !call cllsub(3,dt)
  !write(*,*) '--bc--115--',Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1)
   call cllsub(4,dt,dtifEVOgrv,ifEVOgrv2)
  !  write(*,*) 'BC------8'
  !ifEVOgrv=1+mod(i,6)
   !ifEVOgrv2=1+mod(i,6)
   !---debug---
   !call  SELFGRAVWAVE(0.0d0,4)
   !write(*,*) '--bc--116--',Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1)
  i = i+1
end subroutine slvmuscle


subroutine cllsub(mode,dt,ifEVOgrv,ifEVOgrv2)
  use comvar
  !use grvvar
  use slfgrv
  INCLUDE 'mpif.h'
  integer :: mode !,ifEVOgrv=1,ifEVOgrv2=1
  double precision dt
  integer :: n,m,l
  double precision rho(-1:ndx,-1:ndy,-1:ndz)
  integer ifEVOgrv,ifEVOgrv2
  !do l=-1,ndz
  !do m=-1,ndy
  !do n=-1,ndx
  !   rho(n,m,l) = U(n,m,l,1)
  !end do;end do;end do
  if(mode==1) then
     call BCgrv(1)
     call BCgrv(8)
     call BCgrv(28)
     call BCgrv(9)
     call BCgrv(29)
  end if

  if(mode==2) then
     !call time(dt)
     do l=-1,ndz
     do m=-1,ndy
     do n=-1,ndx
        rho(n,m,l) = U(n,m,l,1)
     end do;end do;end do
     iwx=1; iwy=1; iwz=1
     call BCgrv(101)
     call BCgrv(102)
     call timesource(Phicgp,Phi1step,dt,2)
     call timesource(Phi1step,rho,dt,1)
     call timesource(Phicgm,Phi2step,dt,2)
     call timesource(Phi2step,rho,dt,1)
  end if

  if(mode==3) then
     call BCgrv(4)
     call BCgrv(3)
     call BCgrv(18)
     call BCgrv(38)
     call BCgrv(19)
     call BCgrv(39)
  end if

  if(mode==4) then
     if(ifEVOgrv.eq.1) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 2
        goto 1470
     end if
     if(ifEVOgrv.eq.2) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !ifEVOgrv = 3
        goto 1470
     end if
     if(ifEVOgrv.eq.3) then
        !write(*,*)'in3-1'
        !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        iwx=0; iwy=0; iwz=1
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
        !write(*,*)'in3-2'
        iwx=1; iwy=0; iwz=0
        call BCgrv(102)
        !write(*,*)'in3-2-1'
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        !write(*,*)'in3-2-2'
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !write(*,*)'in3-3'
        iwx=0; iwy=1; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 4
        goto 1470
     end if
     if(ifEVOgrv.eq.4) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 5
        goto 1470
     end if
     if(ifEVOgrv.eq.5) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 6
        goto 1470
     end if
     if(ifEVOgrv.eq.6) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(102)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !ifEVOgrv = 1
        goto 1470
     end if
1470 continue
  end if


  if(mode==5) then
     if(ifEVOgrv2.eq.1) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 2
        goto 1788
     end if
     if(ifEVOgrv2.eq.2) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        !ifEVOgrv2 = 3
        goto 1788
     end if
     if(ifEVOgrv2.eq.3) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 4
        goto 1788
     end if
     if(ifEVOgrv2.eq.4) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 5
        goto 1788
     end if
     if(ifEVOgrv2.eq.5) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 6
        goto 1788
     end if
     if(ifEVOgrv2.eq.6) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(101)
        call BCgrv(102)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        !ifEVOgrv2 = 1
        goto 1788
     end if
1788 continue
  end if
end subroutine cllsub

subroutine BCgrv(mode)
  use comvar
  use mpivar
  !use grvvar
  use slfgrv
  INCLUDE 'mpif.h'
  integer ::  N_ol=2,i,mode,j,k
  INTEGER :: MSTATUS(MPI_STATUS_SIZE)
  DOUBLE PRECISION  :: VECU
  double precision , dimension(1:2) :: pl,pr

CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

if(mode==101) then

IF(iwx.EQ.1) THEN
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!********************************  BC for the leftsides of domains  *****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
       Phicgp(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
       Phicgm(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.0) THEN
     !DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 1
     DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phicgp(IX,JY,KZ)= bphi1l(JY,KZ,IX)
     Phicgm(IX,JY,KZ)= bphi2l(JY,KZ,IX)
  END DO;END DO;END DO
  END IF
 !END DO
 !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(1            ,-1,-1),1,VECU,LEFT,1, &
       Phicgp(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(1            ,-1,-1),1,VECU,LEFT,1, &
       Phicgm(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

  IF(IST.eq.NSPLTx-1) THEN
     !DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx, Ncellx+N_ol
     DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phicgp(IX,JY,KZ)= bphi1r(JY,KZ,IX)
     Phicgm(IX,JY,KZ)= bphi2r(JY,KZ,IX)
  END DO;END DO;END DO
END IF
  !END DO
!************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF


IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
       Phicgp(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
       Phicgm(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,1            ,-1),1,VECU,BOTM,1, &
       Phicgp(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,1            ,-1),1,VECU,BOTM,1, &
       Phicgm(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF


IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
       Phicgp(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
       Phicgm(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phicgp(-1,-1,1            ),1,VECU,DOWN,1, &
       Phicgp(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phicgm(-1,-1,1            ),1,VECU,DOWN,1, &
       Phicgm(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF
endif





if(mode==102) then

IF(iwx.EQ.1) THEN
!  write(*,*) '--bc--1--',Phi1step(1,-1,-1),Phi1step(Ncellx+1,-1,-1),Phi2step(1,-1,-1),Phi2step(Ncellx+1,-1,-1),NRANK
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!********************************  BC for the leftsides of domains  *****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
       Phi1step(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(Ncellx+1-N_ol,-1,-1),1,VECU,RIGT,1, &
       Phi2step(       1-N_ol,-1,-1),1,VECU,LEFT,1, MPI_COMM_WORLD,MSTATUS,IERR)

!  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!   write(*,*) '--bc--2--',Phi1step(1,-1,-1),Phi1step(Ncellx+1,-1,-1),Phi2step(1,-1,-1),Phi2step(Ncellx+1,-1,-1)
  IF(IST.eq.0) THEN
     !DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 1
  DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     Phi1step(IX,JY,KZ)= bstep1l(JY,KZ,IX)
     Phi2step(IX,JY,KZ)= bstep2l(JY,KZ,IX)
  END DO;END DO;END DO
END IF
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!   write(*,*) '--bc--3--',Phi1step(1,-1,-1),Phi1step(Ncellx+1,-1,-1),Phi2step(1,-1,-1),Phi2step(Ncellx+1,-1,-1)
 !END DO
 !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(1            ,-1,-1),1,VECU,LEFT,1, &
       Phi1step(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(1            ,-1,-1),1,VECU,LEFT,1, &
       Phi2step(Ncellx+1     ,-1,-1),1,VECU,RIGT,1, MPI_COMM_WORLD,MSTATUS,IERR)

!  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!  write(*,*) '--bc--4--'
  IF(IST.eq.NSPLTx-1) THEN
     !DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx, Ncellx+N_ol
     DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     Phi1step(IX,JY,KZ)= bstep1r(JY,KZ,IX)
     Phi2step(IX,JY,KZ)= bstep2r(JY,KZ,IX)
  END DO;END DO;END DO
END IF
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!write(*,*) '--bc--5--'
  !END DO
!************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF

!write(*,*) '--bc--6--',NRANK,Phi1step(-1,1,-1),Phi1step(-1,Ncelly+1,-1),Phi2step(-1,1,-1),Phi2step(-1,Ncelly+1,-1)
IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
       Phi1step(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,Ncelly+1-N_ol,-1),1,VECU,TOP ,1, &
       Phi2step(-1,       1-N_ol,-1),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(-1,1            ,-1),1,VECU,BOTM,1, &
       Phi1step(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,1            ,-1),1,VECU,BOTM,1, &
       Phi2step(-1,Ncelly+1     ,-1),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF

!write(*,*) '--bc--7--',NRANK,Phi1step(-1,-1,1),Phi1step(-1,-1,Ncellz+1),Phi2step(-1,-1,1),Phi2step(-1,-1,Ncellz+1)
IF(iwz.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(1,N_ol*(ndx+2)*(ndy+2),N_ol*(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
!*************************************  BC for the downsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
       Phi1step(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,-1,Ncellz+1-N_ol),1,VECU,UP  ,1, &
       Phi2step(-1,-1,       1-N_ol),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!**************************************  BC for the upsides of domains  ****
  !DO K = 1, N_MPI(20)
  CALL MPI_SENDRECV(Phi1step(-1,-1,1            ),1,VECU,DOWN,1, &
       Phi1step(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  CALL MPI_SENDRECV(Phi2step(-1,-1,1            ),1,VECU,DOWN,1, &
       Phi2step(-1,-1,Ncellz+1     ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  !END DO
!***************************************************************************
  CALL MPI_TYPE_FREE(VECU,IERR)
END IF
endif
end subroutine BCgrv


subroutine muslcslv1D(Phiv,source,dt,mode,hazi)
  use comvar
  double precision :: nu2 , w=6.0d0 , dt2 , dt , deltap,deltam !kappa -> comver  better?
  integer :: direction , mode , invdt , loopmode , dloop,cnt=0
  !DOUBLE PRECISION :: fluxf(-1:ndx,-1:ndy,-1:ndz),fluxg(-1:ndx,-1:ndy,-1:ndz)
  !DOUBLE PRECISION, dimension(-1:ndx) :: Phigrad,Phipre,fluxphi,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  DOUBLE PRECISION, dimension(-1:ndx,-1:ndy,-1:ndz) :: Phigrad,Phipre,fluxphi&
       ,Phiv,source,Phi2dt,Phiu,sourcepre,sourcepri
  character(5) name
  integer Ncell,Ncm,Ncl,ix,jy,kz,Lnum,Mnum,hazi,is,ie
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G


  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy; endif! BT1 = 1; BT2 = 2; VN = 4; end if

  !if(iwx.eq.1) then
  !   Ncell = ndx
  !   Ncm = ndy
  !   Ncl = ndz
  !endif
  !if(iwy.eq.1) then
  !   Ncell = ndy
  !   Ncm = ndz
  !   Ncl = ndx
  !endif
  !if(iwz.eq.1) then
  !   Ncell = ndz
  !   Ncm = ndx
  !   Ncl = ndy
  !endif

  !----kyoukai-----
   if(hazi==1)then
      is = 2
      ie = Ncell-3
   end if
   if(hazi==2)then
      is = 1
      ie = Ncell-2
   end if
   !----kyoukai-----
  nu2 = cg * dt / dx(1)
  Phipre(:,:,:) = Phiv(:,:,:)
  !write(name,'(i5.5)') cnt
  !------------ul.solver.+cg-------------
  if(mode==1) then
     call fluxcal(Phipre,Phipre,Phiu,0.0d0,1.d0/3.0d0,10,is,ie)
     !call fluxcal(Phipre,Phipre,Phiu,0.0d0,0.0d0,10)
     !------------calcurate dt/2------------
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
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
      DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
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
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
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
     DO Lnum = 1, Ncl
        DO Mnum = 1, Ncm
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
           do i=is,ie
              Phiv(i,j,k) =  -cg * G4pi * source(i,j,k) * dt + Phipre(i,j,k)
           end do
        end do
     end do
  end if

  if(mode==4) then
     do k = 1,ndz-2
        do j = 1,ndy-2
           do i=is,ie
              Phiv(i,j,k) = cg * source(i,j,k) * dt + Phipre(i,j,k)
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


  if(iwx.eq.1) Ncell = ndx
  if(iwy.eq.1) Ncell = ndy
  if(iwz.eq.1) Ncell = ndz

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
  DOUBLE PRECISION , dimension(-1:ndx) :: slop
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  DOUBLE PRECISION, parameter :: G=1.11142d-4, G4pi=12.56637d0*G
  !uin(:)=0.0d0
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
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
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
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
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
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
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
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

subroutine pbstep()
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
integer i,j,k
double precision dxx
character(3) fn
dxx=dx(1)
!iwx = 0; iwy = 1; iwz = 1
!call BCgrv(101)
!call BCgrv(102)
!bphi1l(j,k,1)
!write(fn,'(i3.3)') NRANK
!open(3,file='/work/maedarn/3DMHD/test/bcstep'//fn//'.DAT')

do k=1,Ncellz
do j=1,Ncelly

bstep1l(j,k,1)=(3.0d0*bphi1l(j,k,1)-4.0d0*bphi1l(j,k,0)+bphi1l(j,k,-1))*0.5d0/dxx
bstep1l(j,k,0)=(-bphi1l(j,k,-1)+bphi1l(j,k,1))*0.5d0/dxx
bstep1l(j,k,-1)=-(3.0d0*bphi1l(j,k,-1)-4.0d0*bphi1l(j,k,0)+bphi1l(j,k,1))*0.5d0/dxx
bstep1r(j,k,Ncellx+2)=(3.0d0*bphi1r(j,k,Ncellx+2)-4.0d0*bphi1r(j,k,Ncellx+1)+bphi1r(j,k,Ncellx))*0.5d0/dxx
bstep1r(j,k,Ncellx+1)=(-bphi1r(j,k,Ncellx)+bphi1r(j,k,Ncellx+2))*0.5d0/dxx
bstep1r(j,k,Ncellx)=-(3.0d0*bphi1r(j,k,Ncellx)-4.0d0*bphi1r(j,k,Ncellx+1)+bphi1r(j,k,Ncellx+2))*0.5d0/dxx

bstep2l(j,k,1)=-(3.0d0*bphi2l(j,k,1)-4.0d0*bphi2l(j,k,0)+bphi2l(j,k,-1))*0.5d0/dxx
bstep2l(j,k,0)=-(-bphi2l(j,k,-1)+bphi2l(j,k,1))*0.5d0/dxx
bstep2l(j,k,-1)=(3.0d0*bphi2l(j,k,-1)-4.0d0*bphi2l(j,k,0)+bphi2l(j,k,1))*0.5d0/dxx
bstep2r(j,k,Ncellx+2)=-(3.0d0*bphi2r(j,k,Ncellx+2)-4.0d0*bphi2r(j,k,Ncellx+1)+bphi2r(j,k,Ncellx))*0.5d0/dxx
bstep2r(j,k,Ncellx+1)=-(-bphi2r(j,k,Ncellx)+bphi2r(j,k,Ncellx+2))*0.5d0/dxx
bstep2r(j,k,Ncellx)=(3.0d0*bphi2r(j,k,Ncellx)-4.0d0*bphi2r(j,k,Ncellx+1)+bphi2r(j,k,Ncellx+2))*0.5d0/dxx

!write(3,*) bstep1l(j,k,1),bstep1l(j,k,0),bstep1l(j,k,-1),bstep1r(j,k,Ncellx+2),bstep1r(j,k,Ncellx+1),bstep1r(j,k,Ncellx),&
!     bstep2l(j,k,1),bstep2l(j,k,0),bstep2l(j,k,-1),bstep2r(j,k,Ncellx+2),bstep2r(j,k,Ncellx+1),bstep2r(j,k,Ncellx)
end do
end do
!close(3)
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

DOUBLE PRECISION, dimension(:,:),  allocatable :: dat1,dat2
complex*16, dimension(:), allocatable :: spe1,spe2
double precision :: kap,temp1r,temp1i,temp2r,temp2i,facG,fac,dxx,dyy,dzz,zp1,zp2
double precision, dimension(:,:,:), allocatable :: fint0,fint1
double precision, dimension(:,:,:,:), allocatable :: bcsend
character*4 fnum
character(3) fn,deep
!character(3) lname
character(2) lcRANK
character(1) lRANK
integer pls,mis

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

!nccy = Ncelly/NSPLTy; nccz = Ncellz/NSPLTz
nccy = Ncelly; nccz = Ncellz
do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
do i=-1,Ncellx+2
  data(jy,kz,i) = U(i,j,k,1)
end do;end do;end do

                    !count, blocklength, stride
CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)

do Nlp = 1,NSPLTy*NSPLTz-1

  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE
  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs
  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE
  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr

  Nis = JSs + NSPLTy*KSs
  kls = Nis + 1 !+ pls
  Nir = JST + NSPLTy*KST
  klr = Nir + 1 !+ pls

  !***************fordebug*****************
  !CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
  !***************fordebug*****************

  if(kls.gt.Ncellx) then; isend = MPI_PROC_NULL; kls = Ncellx+3; end if
  if(klr.gt.Ncellx) then; irecv = MPI_PROC_NULL; klr = Ncellx+3; end if
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

if(klr.le.Ncellx+2) then

  call rlft3(data(1,1,klr),speq(1,klr),nn1,nn2,1)

  kz = klr
  zp1 = x(kz)-0.5d0*dzz
  zp2 = Lbox - zp1
  !zp2 = Lbox/dble(NSPLTx) - zp1
  write(*,*) zp1,zp2,100.d0-zp2,x(kz),kz,'write-zp',NRANK,pls
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

!write(fn,'(i3.3)') NRANK
!write(lRANK,'(i1.1)') pls+2
!open(12,file='/work/maedarn/3DMHD/test/bcsave'//lRANK//fn//'.DAT')
!ncx=Ncellx+1; ncy=Ncelly+1; ncz=Ncellz+1
!do k=0,ncz; kk= (ncy+1)*k
!do j=0,ncy; n = j+kk
!ncx=Ncellx-2
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!if(NRANK==0) then
ncy=Ncelly; ncz=Ncellz
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!do Nroot=0,NPE-1
do k=1,ncz!; kk= (ncy+1)*k
do j=1,ncy!; n = j+kk
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  !if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 1
  !if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 1
  !if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  !if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz

  bphi1l(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  bphi1r(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
  bphi2l(j,k,1-abs(pls)) = dble(data(jb,kbb,1))
  bphi2r(j,k,Ncellx+abs(pls)) = dble(data(jb,kbb,2))
  !bcsend(j,k,1-abs(pls),Nroot)=dble(data(jb,kbb,1))
  !  write(12,*) bphi1l(j,k,1-abs(pls)),bphi1r(j,k,Ncellx+abs(pls)),bphi2l(j,k,1-abs(pls)), bphi2r(j,k,Ncellx+abs(pls))
  !CALL MPI_BCAST(bcsend(1,1,1,Nroot),(Ncelly)*(Ncellz)*(1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
end do
end do
!close(12)
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!do Nroot=0,NPE-1
!CALL MPI_BCAST(bcsend(1,1,1-abs(pls),Nroot),(Ncelly)*(Ncellz)*(1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
!  CALL MPI_BCAST(tMPI(0,0,0,Nroot),(nx2+1)*(ny2+1)*(nz2+1),MPI_REAL8,Nroot,MPI_COMM_WORLD,IERR)
!end do
!end if

DEALLOCATE(data,speq)
!DEALLOCATE(bcsend)
DEALLOCATE(dat1,spe1,dat2,spe2)
!-----------------------------------------------------------------

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
                 sdt = sourratio*dabs(Phiv(i,j,k)) / (cg * G4pi * source(i,j,k) )
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
                 sdt = sourratio*dabs(Phiv(i,j,k)) / ( cg * source(i,j,k) )
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
