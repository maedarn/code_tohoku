subroutine SELFGRAVWAVE(dt,mode)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'

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


     !write(*,*) '------pb1-------' ,Nrank
     !Call PB()
     !write(*,*) '------pb2-------' ,Nrank


     !*********use phi exact**********
     if(IST.eq.0       ) then
        !write(*,*) NRANK,Phi(1     ,1,1),'boundary'
        do k=1,Ncellz; do j=1,Ncelly
           Phi(1     ,j,k)=PBexatest
           Phi(0       ,j,k) = Phi(1     ,j,k); Phi(-1       ,j,k) = Phi(1     ,j,k) !grad=0
     end do; end do; end if
     if(IST.eq.NSPLTx-1) then
        !write(*,*) NRANK,Phi(Ncellx     ,1,1),'boundary'
        do k=1,Ncellz; do j=1,Ncelly
           Phi(Ncellx,j,k)=PBexatest
           Phi(Ncellx+1,j,k) = Phi(Ncellx,j,k); Phi(Ncellx+2,j,k) = Phi(Ncellx,j,k)
     end do; end do; end if
     !*********use phi exact**********



     !call mglin(Nmem1,Nmem2,2,5,5)
     !write(*,*) NRANK,Phi(0,0,0),Phi(1,1,1),Phi(Ncellx,Ncelly,Ncellz),dt,'-------33--333---33---'
     write(*,*) '------gr1-------' ,Nrank
     call         end do
  end do
end do
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
     !open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     !open(unit=38,file='/work/maedarn/3DMHD/test/PHIDTINI/INIPHI'//NPENUM//countcha//'.DAT',FORM='UNFORMATTED')!,FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=8,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgp'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=18,file='/work/maedarn/3DMHD/test/PHIINI/INIPHIcgm'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=28,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI1step'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
     open(unit=38,file='/work/maedarn/3DMHD/test/PHIINI/INIPHI2step'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
   call STBLphi(sourcedt,Phipregrad)
   dt = sourcedt
end if
!***************SABILITY-exa**************


!***************pointerforPB**************
if(mode.eq.10) then
  call pinter(Nmem1,Nmem2,Ncellx,Ncelly,Ncellz)
end if
!***************pointerforPB**************

end subroutine SELFGRAVWAVE


subroutine slvmuscle


  call cllsub(2,dt)
  call cllsub(3,dt)
  call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
  call muslcslv1D(Phi2step,rho,dt*0.5d0,3,2)
  call cllsub(3,dt)
  call cllsub(4,dt)
  call cllsub(3,dt)
  call cllsub(1,dt)
  call muslcslv1D(Phicgp,Phi1step,dt,4,2)
  call muslcslv1D(Phicgm,Phi2step,dt,4,2)
  call cllsub(1,dt)
  call cllsub(5,dt)
  call cllsub(3,dt)
  call muslcslv1D(Phi1step,rho,dt*0.5d0,3,2)
  call muslcslv1D(Phi2step,rho,0.5d0*dt,3,2)
  call cllsub(3,dt)
  call cllsub(4,dt)
  ifEVOgrv=1+mod(i,6)
  ifEVOgrv2=1+mod(i,6)

end subroutine slvmuscle


subroutine cllsub(mode,dt)
  use comvar
  use grvvar
  integer :: mode !,ifEVOgrv=1,ifEVOgrv2=1
  double precision dt
  if(mode==1) then
     call BCgrv(1)
     call BCgrv(8)
     call BCgrv(28)
     call BCgrv(9)
     call BCgrv(29)
  end if

  if(mode==2) then
     call time(dt)
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
        call BCgrv(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BCgrv(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 2
        goto 1470
     end if
     if(ifEVOgrv.eq.2) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BCgrv(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !ifEVOgrv = 3
        goto 1470
     end if
     if(ifEVOgrv.eq.3) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 4
        goto 1470
     end if
     if(ifEVOgrv.eq.4) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BCgrv(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 5
        goto 1470
     end if
     if(ifEVOgrv.eq.5) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BCgrv(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        !ifEVOgrv = 6
        goto 1470
     end if
     if(ifEVOgrv.eq.6) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(19)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(39)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(18)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,2)
        call BCgrv(38)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(3)
        call muslcslv1D(Phi1step,rho,dt*0.5d0,2,1)
        call BCgrv(4)
        call muslcslv1D(Phi2step,rho,dt*0.5d0,1,1)
        !ifEVOgrv = 1
        goto 1470
     end if
1470 continue
  end if


  if(mode==5) then
     if(ifEVOgrv2.eq.1) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BCgrv(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 2
        goto 1788
     end if
     if(ifEVOgrv2.eq.2) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=0; iwz=1
        call BCgrv(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BCgrv(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        !ifEVOgrv2 = 3
        goto 1788
     end if
     if(ifEVOgrv2.eq.3) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BCgrv(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=1; iwz=0
        call BCgrv(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 4
        goto 1788
     end if
     if(ifEVOgrv2.eq.4) then
        iwx=1; iwy=0; iwz=0
        call BCgrv(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BCgrv(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 5
        goto 1788
     end if
     if(ifEVOgrv2.eq.5) then
        iwx=0; iwy=1; iwz=0
        call BCgrv(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BCgrv(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        iwx=0; iwy=0; iwz=1
        call BCgrv(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        !ifEVOgrv2 = 6
        goto 1788
     end if
     if(ifEVOgrv2.eq.6) then
        iwx=0; iwy=0; iwz=1
        call BCgrv(9)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(29)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=0; iwy=1; iwz=0
        call BCgrv(8)
        call muslcslv1D(Phicgp,Phi1step,dt,1,2)
        call BCgrv(28)
        call muslcslv1D(Phicgm,Phi2step,dt,2,2)
        iwx=1; iwy=0; iwz=0
        call BCgrv(1)
        call muslcslv1D(Phicgp,Phi1step,dt,1,1)
        call BCgrv(1)
        call muslcslv1D(Phicgm,Phi2step,dt,2,1)
        !ifEVOgrv2 = 1
        goto 1788
     end if
1788 continue
  end if
end subroutine cllsub

subroutine BCgrv(mode)
  use comvar
  use grvvar
  integer :: i,mode,j,k
  double precision , dimension(1:2) :: pl,pr

  if(mode==1) then
     !---------kotei-x-----------
     !goto 100
     !---------Phi-------------
     Phicgp(1,:,:)= bcphi1(1,:,:)
     Phicgp(0,:,:)= bcphi1(2,:,:)
     Phicgp(-1,:,:)= bcphi1(3,:,:)
     Phicgp(ndx-2,:,:)= bcphi2(1,:,:)
     Phicgp(ndx-1,:,:)= bcphi2(2,:,:)
     Phicgp(ndx,:,:)= bcphi2(3,:,:)

     Phicgm(1,:,:)= bcphi1(1,:,:)
     Phicgm(0,:,:)= bcphi1(2,:,:)
     Phicgm(-1,:,:)= bcphi1(3,:,:)
     Phicgm(ndx-2,:,:)= bcphi2(1,:,:)
     Phicgm(ndx-1,:,:)= bcphi2(2,:,:)
     Phicgm(ndx,:,:)= bcphi2(3,:,:)

     !---------Phi-------------

  end if

  if(mode==3)then
     !---------kotei-x-----------
     !-------Phi1step+cg-----------
     !goto 700
     Phi1step(1,:,:)= Phigrd(1,:,:)
     Phi1step(0,:,:)= Phigrd(0,:,:)
     Phi1step(-1,:,:)=Phigrd(-1,:,:)
     Phi1step(ndx-2,:,:)= Phigrd(ndx-2,:,:)
     Phi1step(ndx-1,:,:)= Phigrd(ndx-1,:,:)
     Phi1step(ndx,:,:)= Phigrd(ndx,:,:)
     !700 continue
     !-------Phi1step-----------
     !---------kotei-x-----------
  end if

  if(mode==4) then
     !---------kotei-x-----------
     !-------Phi1step-cg-----------
     !goto 701
     Phi2step(1,:,:)= -Phigrd(1,:,:)
     Phi2step(0,:,:)= -Phigrd(2,:,:)
     Phi2step(-1,:,:)=-Phigrd(3,:,:)
     Phi2step(ndx-2,:,:)= -Phigrd(ndx-2,:,:)
     Phi2step(ndx-1,:,:)= -Phigrd(ndx-1,:,:)
     Phi2step(ndx,:,:)= -Phigrd(ndx,:,:)
     !701 continue
     !-------Phi1step-----------
     !---------kotei-x-----------
  end if


  !----- y-periodic ------
  if(mode==8) then
     !-------period2-----------
     !     goto 102
     do k=-1,ndz
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgp(i,ndx-2,k)
           pr(1)= Phicgp(i,ndx-3,k)
           pl(1)= Phicgp(i,1,k)
           pl(2)= Phicgp(i,2,k)
           Phicgp(i,-1,k)=pr(1)
           Phicgp(i,0,k)=pr(2)
           Phicgp(i,ndx-1,k)=pl(1)
           Phicgp(i,ndx,k)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==18) then
     do k=-1,ndz
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi1step(i,ndx-2,k)
           pr(1)= Phi1step(i,ndx-3,k)
           pl(1)= Phi1step(i,1,k)
           pl(2)= Phi1step(i,2,k)
           Phi1step(i,-1,k)=pr(1)
           Phi1step(i,0,k)=pr(2)
           Phi1step(i,ndx-1,k)=pl(1)
           Phi1step(i,ndx,k)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if




  if(mode==28) then
     !-------period2-----------
     !     goto 102
     do k=-1,ndz
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgm(i,ndx-2,k)
           pr(1)= Phicgm(i,ndx-3,k)
           pl(1)= Phicgm(i,1,k)
           pl(2)= Phicgm(i,2,k)
           Phicgm(i,-1,k)=pr(1)
           Phicgm(i,0,k)=pr(2)
           Phicgm(i,ndx-1,k)=pl(1)
           Phicgm(i,ndx,k)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==38) then
     do k=-1,ndz
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi2step(i,ndx-2,k)
           pr(1)= Phi2step(i,ndx-3,k)
           pl(1)= Phi2step(i,1,k)
           pl(2)= Phi2step(i,2,k)
           Phi2step(i,-1,k)=pr(1)
           Phi2step(i,0,k)=pr(2)
           Phi2step(i,ndx-1,k)=pl(1)
           Phi2step(i,ndx,k)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if


  !----- y-periodic ------


  !----- z-periodic ------
  if(mode==9) then
     !-------period2-----------
     !     goto 102
     do j=-1,ndy
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgp(i,j,ndx-2)
           pr(1)= Phicgp(i,j,ndx-3)
           pl(1)= Phicgp(i,j,1)
           pl(2)= Phicgp(i,j,2)
           Phicgp(i,j,-1)=pr(1)
           Phicgp(i,j,0)=pr(2)
           Phicgp(i,j,ndx-1)=pl(1)
           Phicgp(i,j,ndx)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==19) then
     do j=-1,ndy
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi1step(i,j,ndx-2)
           pr(1)= Phi1step(i,j,ndx-3)
           pl(1)= Phi1step(i,j,1)
           pl(2)= Phi1step(i,j,2)
           Phi1step(i,j,-1)=pr(1)
           Phi1step(i,j,0)=pr(2)
           Phi1step(i,j,ndx-1)=pl(1)
           Phi1step(i,j,ndx)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if




  if(mode==29) then
     !-------period2-----------
     !     goto 102
     do j=-1,ndy
        do i=-1,ndx
           !---------Phi-------------
           pr(2)= Phicgm(i,j,ndx-2)
           pr(1)= Phicgm(i,j,ndx-3)
           pl(1)= Phicgm(i,j,1)
           pl(2)= Phicgm(i,j,2)
           Phicgm(i,j,-1)=pr(1)
           Phicgm(i,j,0)=pr(2)
           Phicgm(i,j,ndx-1)=pl(1)
           Phicgm(i,j,ndx)=pl(2)
           !---------Phi-------------
        end do
     end do
  end if
  if(mode==39) then
     do j=-1,ndy
        do i=-1,ndx
           !-------Phi1step-----------
           pr(2)= Phi2step(i,j,ndx-2)
           pr(1)= Phi2step(i,j,ndx-3)
           pl(1)= Phi2step(i,j,1)
           pl(2)= Phi2step(i,j,2)
           Phi2step(i,j,-1)=pr(1)
           Phi2step(i,j,0)=pr(2)
           Phi2step(i,j,ndx-1)=pl(1)
           Phi2step(i,j,ndx)=pl(2)
           !-------Phi1step-----------
!102  continue
           !-------period2-----------
        end do
     end do
  end if


  !----- z-periodic ------
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


  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz; endif!  BT1 = 2; BT2 = 3; VN = 2; end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx; endif! BT1 = 3; BT2 = 1; VN = 3; end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy; endif! BT1 = 1; BT2 = 2; VN = 4; end if


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
  nu2 = cg * dt / dx
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

subroutine fluxcal(preuse,pre,u,ep,kappa,mode,is,ie)
  use comvar
  double precision :: ep , kappa
  DOUBLE PRECISION , dimension(-1:ndx,-1:ndy,-1:ndz) :: ul,ur,pre,preuse,u
  DOUBLE PRECISION , dimension(-1:ndx) :: slop
  integer :: i,mode,Ncell,Ncl,Ncm,j,k,Lnum,Mnum
  integer ix,jy,kz,ixp,jyp,kzp,ixm,jym,kzm,is,ie
  !u(:)=0.0d0
  if(iwx.eq.1) then; Ncell = ndx; Ncm = ndy; Ncl = ndz;  end if
     if(iwy.eq.1) then; Ncell = ndy; Ncm = ndz; Ncl = ndx;  end if
        if(iwz.eq.1) then; Ncell = ndz; Ncm = ndx; Ncl = ndy;  end if

           !call vanalbada(pre,slop)
           if(mode==1) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(Ncm,Ncl,pre,slop,is,ie,Ncell)
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
              u(ix,jy,kz)=ul(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !u(:)=ul(:)
           end if


           if(mode==4) then
              DO Lnum = 1, Ncl
              DO Mnum = 1, Ncm
              call vanalbada(Ncm,Ncl,pre,slop,is,ie,Ncell)
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
              u(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
              !write(*,*) slop(127),'127slop'
              !write(*,*) slop(ndx-ien),ndx-ien,slop(ndx-ien+1)
              !write(*,*) u(2)
              !u(:)=ur(:)
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
              u(ix,jy,kz)=ul(ix,jy,kz)
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
              u(ix,jy,kz)=ur(ix,jy,kz)
              end do
              end DO
              end DO
           end if


end subroutine fluxcal
