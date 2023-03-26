RECURSIVE subroutine SELFGRAVWAVE(dt,mode)
  USE comvar
  USE mpivar
  USE slfgrv
  INCLUDE 'mpif.h'
  integer :: mode,MRANK,count=0,rdnum,rddmy,ndt1=0,svc1=0!,svci=50
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
  integer nt
  double precision :: tdm
  !double precision , dimension(:,:,:) , allocatable :: stbPhi
  !double precision , dimension(-1:Ncellx+2,-1:Ncelly,-1:Ncellz) :: Phipregrad,Phipregraddum
  !**************** INITIALIZEATION **************
  if(mode==0) then
     Phi(:,:,:) = 0.0d0
     Phiwv(:,:,:,:)=0.d0
     Phigrdwv(:,:,:,:)=0.d0

     do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
         Phiwv(i,j,k,1)=Phiexa(i,j,k)
         !Phigrdwv(i,j,k,1)=cg*Phigrd(i,j,k,1)*2.d0/2.d0/kappa+Phiwv(i,j,k,1)
         Phigrdwv(i,j,k,1)=cg*Phigrd(i,j,k,1)+kappa*Phiexa(i,j,k)
     end do
     end do
     end do
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
     !ndt1=ndt1+1
     !if(ndt1==maxstp*nitera/mvstp+1) then
     !call  SELFGRAVWAVE(0.0d0,4)
!call move()
     Phiwv(:,:,:,:)=0.d0
     Phigrdwv(:,:,:,:)=0.d0
     !ndt1=0
     !endif
     !call movesph(dt)
     tdm=dt!/dble(ntdiv)
     do nt=1,ntdiv
     iwx=1;iwy=1;iwz=1
     call BCgrv(100,1,8)
     if(mod(svc1,svci)==0) then
     call SELFGRAVWAVE(0.0d0,4)
     endif
     call slvmuscle(tdm)
     svc1=svc1+1
     enddo
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
     open(unit=28,file=dir//svdir//'/PHI'//countcha//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
           !write(28) sngl(Phiwv(i,j,k,1)),sngl(Phigrdwv(i,j,k,1)),sngl(Phiexa(i,j,k)),sngl(Phigrd(i,j,k,1)),sngl(U(i,j,k,1))
           write(28) sngl(Phiwv(i,j,k,1)),sngl(Phigrdwv(i,j,k,1)),sngl(Phiexa(i,j,k)),sngl(cg*Phigrd(i,j,k,1)+kappa*Phiexa(i,j,k)),sngl(U(i,j,k,1))
          enddo
        end do
        !write(*,*) sngl(Phiwv(8,8,k,1)),sngl(Phigrdwv(8,8,k,1))
     end do
     close(28)
     count=count+1
  end if

  if(mode==40) then
     !write(*,*) 'save???'
     WRITE(NPENUM,'(I3.3)') NRANK
     !WRITE(countcha,'(I6.6)') count
     open(unit=28,file=dir//svdir//'/PHIF'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
           write(28) Phiwv(i,j,k,1),Phigrdwv(i,j,k,1),Phiexa(i,j,k),Phigrd(i,j,k,1),U(i,j,k,1)
          enddo
        end do
        !write(*,*) sngl(Phiwv(8,8,k,1)),sngl(Phigrdwv(8,8,k,1))
     end do
     close(28)
     !count=count+1
  end if

    if(mode==41) then
     !write(*,*) 'save???'
     WRITE(NPENUM,'(I3.3)') NRANK
     !WRITE(countcha,'(I6.6)') count
     open(unit=28,file=dir//svdir//'/PHIF'//NPENUM//'.DAT',FORM='UNFORMATTED') !,CONVERT='LITTLE_ENDIAN')
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
           read(28) Phiwv(i,j,k,1),Phigrdwv(i,j,k,1),Phiexa(i,j,k),Phigrd(i,j,k,1),U(i,j,k,1)
          enddo
        end do
        !write(*,*) sngl(Phiwv(8,8,k,1)),sngl(Phigrdwv(8,8,k,1))
     end do
     close(28)
     !count=count+1
  end if

  !***************SAVE PHI & PHIDT FOR DEBUG & INITIAL**************
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
  double precision  grdxy2zp,grdxy2zm,grdxy2mn,grdyz2xp,grdyz2xm,grdyz2mn,grdzx2yp,grdzx2ym,grdzx2mn
  double precision  grdxy2zpp,grdxy2zmm,grdyz2xpp,grdyz2xmm,grdzx2ypp,grdzx2ymm
  double precision :: adiff2=0.5d0
  double precision dtt2
  double precision :: Phiwvpre(-1:ndx,-1:ndy,-1:ndz,1:2),Phigrdwvpre(-1:ndx,-1:ndy,-1:ndz,1:2)

  dtt2=dt!*0.3d0
  !rhomean=0.d0
  do l=1,ndz-2
  do m=1,ndy-2
  do n=1,ndx-2
     !rho(n,m,l) = U(n,m,l,1)
     rho(n,m,l) = U(n,m,l,1)!-rhomean
  !   rhomean=rhomean+rho(i,j,k)
  end do;end do;end do
  


  !****************slv-wv****************
  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=1;iwy=0;iwz=0
  call BCgrv(100,1,8)
  !Phiwvtest(:,:,:)=NRANK
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)
  !call muslcslv1D(Phiwvtest,dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-1,-1,-1,2),dtt2*0.5d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)
  !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dtt2*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
  !call BCgrv(100,1,8)

  iwx=0;iwy=1;iwz=0
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)
  !call muslcslv1D(Phiwv(-1,-1,-1,2),dtt2*0.5d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,1)
  iwx=0;iwy=1;iwz=0
  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)
  !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dtt2*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,2)
  !call BCgrv(110,1,8)
  !call BCgrv(100,1,8)

  iwx=0;iwy=0;iwz=1
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)
  !call muslcslv1D(Phiwv(-1,-1,-1,2),dtt2*0.5d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
  !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
  !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)
  !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dtt2*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
  !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
  !call BCgrv(110,1,8)
  !call BCgrv(100,1,8)

do k=-1,ndz
do j=-1,ndy
do i=-1,ndx
  Phiwvpre(i,j,k,1)=Phiwv(i,j,k,1)
  !Phiwvpre(i,j,k,2)=Phiwv(i,j,k,2)
  Phigrdwvpre(i,j,k,1)=Phigrdwv(i,j,k,1)
  !Phigrdwvpre(i,j,k,2)=Phigrdwv(i,j,k,2)
enddo
enddo
enddo



  do k=1,ndz-2
     do j=1,ndy-2
        do i=1,ndx-2
Phiwv(i,j,k,1) = 0.5d0*Phiwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa * 0.5d0 * dt))+Phigrdwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa * 0.5d0 * dt))/(2.d0*kappa+1.d-10)

Phigrdwv(i,j,k,1) = 0.5d0*kappa*Phiwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa *0.5d0* dt))+0.5d0*Phigrdwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa *0.5d0* dt))
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


     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0*2.d0*kappa * dt)+(Phiwv(i,j,k,1) &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdxy1/dx1/dy1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdyz1/dy1/dz1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdzx1/dz1/dx1 &
     !-cg*cg*4.d0/2.d0/kappa/2.d0/kappa*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0*2.d0*kappa * dt))

     Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
     (-2.d0*cg*cg*grdxy1/dx1/dy1 &
      -2.d0*cg*cg*grdyz1/dy1/dz1 &
      -2.d0*cg*cg*grdzx1/dz1/dx1) *dt &
     -G4pi*cg*cg*rho(i,j,k)*dt


     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
     !(-2.d0*2.d0/2.d0/kappa*cg*cg*grdxy1/dx1/dy1 &
     ! -2.d0*2.d0/2.d0/kappa*cg*cg*grdyz1/dy1/dz1 &
     ! -2.d0*2.d0/2.d0/kappa*cg*cg*grdzx1/dz1/dx1) *dt

     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0*2.d0*kappa * dt)+(Phiwv(i,j,k,1) &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdxy1mn/dx1/dy1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdyz1mn/dy1/dz1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdzx1mn/dz1/dx1 &
     !-cg*cg*4.d0/2.d0/kappa/2.d0/kappa*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0*2.d0*kappa * dt))

            enddo
          enddo
        enddo

do k=-1,ndz
do j=-1,ndy
do i=-1,ndx
  Phiwvpre(i,j,k,1)=Phiwv(i,j,k,1)
  !Phiwvpre(i,j,k,2)=Phiwv(i,j,k,2)
  Phigrdwvpre(i,j,k,1)=Phigrdwv(i,j,k,1)
  !Phigrdwvpre(i,j,k,2)=Phigrdwv(i,j,k,2)
enddo
enddo
enddo


        do k=1,ndz-2
           do j=1,ndy-2
              do i=1,ndx-2
Phiwv(i,j,k,1) = 0.5d0*Phiwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa * 0.5d0 * dt))+Phigrdwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa * 0.5d0 * dt))/(2.d0*kappa+1.d-10)

Phigrdwv(i,j,k,1) = 0.5d0*kappa*Phiwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa *0.5d0* dt))+0.5d0*Phigrdwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa *0.5d0* dt))
              enddo
           enddo
         enddo



 !iwx=1;iwy=1;iwz=1
 !call BCgrv(100,1,8)
 iwx=0;iwy=0;iwz=1
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)
 !call muslcslv1D(Phiwv(-1,-1,-1,2),dtt2*0.5d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)
 !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dtt2*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
 !call BCgrv(100,1,8)

 iwx=0;iwy=1;iwz=0
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)
 !call muslcslv1D(Phiwv(-1,-1,-1,2),dtt2*0.5d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,1)
 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)
 !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dtt2*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,2)
 !call BCgrv(110,1,8)
 !call BCgrv(100,1,8)

 iwx=1;iwy=0;iwz=0
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)
 !call muslcslv1D(Phiwv(-1,-1,-1,2),dtt2*0.5d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,3),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,4),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,5),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,6),dt*0.25d0,2)
 !call muslcslv1D(Phiwv(-2,-2,-2,7),dt*0.25d0,1)
 !call muslcslv1D(Phiwv(-2,-2,-2,8),dt*0.25d0,2)
 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)
 !call muslcslv1D(Phigrdwv(-1,-1,-1,2),dtt2*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,3),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,4),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,5),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,6),dt*0.5d0,1)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,7),dt*0.5d0,2)
 !call muslcslv1D(Phigrdwv(-2,-2,-2,8),dt*0.5d0,1)
 !call BCgrv(100,1,8)
 !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%



  !****************slv-wv****************
  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=0;iwy=1;iwz=0
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

  iwx=0;iwy=0;iwz=1
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

  iwx=1;iwy=0;iwz=0
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

do k=-1,ndz
do j=-1,ndy
do i=-1,ndx
  Phiwvpre(i,j,k,1)=Phiwv(i,j,k,1)
  !Phiwvpre(i,j,k,2)=Phiwv(i,j,k,2)
  Phigrdwvpre(i,j,k,1)=Phigrdwv(i,j,k,1)
  !Phigrdwvpre(i,j,k,2)=Phigrdwv(i,j,k,2)
enddo
enddo
enddo



  do k=1,ndz-2
     do j=1,ndy-2
        do i=1,ndx-2
Phiwv(i,j,k,1) = 0.5d0*Phiwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa * 0.5d0 * dt))+Phigrdwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa * 0.5d0 * dt))/(2.d0*kappa+1.d-10)

Phigrdwv(i,j,k,1) = 0.5d0*kappa*Phiwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa *0.5d0* dt))+0.5d0*Phigrdwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa *0.5d0* dt))
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


     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0*2.d0*kappa * dt)+(Phiwv(i,j,k,1) &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdxy1/dx1/dy1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdyz1/dy1/dz1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdzx1/dz1/dx1 &
     !-cg*cg*4.d0/2.d0/kappa/2.d0/kappa*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0*2.d0*kappa * dt))

     Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
     (-2.d0*cg*cg*grdxy1/dx1/dy1 &
      -2.d0*cg*cg*grdyz1/dy1/dz1 &
      -2.d0*cg*cg*grdzx1/dz1/dx1) *dt &
     -G4pi*cg*cg*rho(i,j,k)*dt


     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
     !(-2.d0*2.d0/2.d0/kappa*cg*cg*grdxy1/dx1/dy1 &
     ! -2.d0*2.d0/2.d0/kappa*cg*cg*grdyz1/dy1/dz1 &
     ! -2.d0*2.d0/2.d0/kappa*cg*cg*grdzx1/dz1/dx1) *dt

     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0*2.d0*kappa * dt)+(Phiwv(i,j,k,1) &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdxy1mn/dx1/dy1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdyz1mn/dy1/dz1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdzx1mn/dz1/dx1 &
     !-cg*cg*4.d0/2.d0/kappa/2.d0/kappa*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0*2.d0*kappa * dt))

            enddo
          enddo
        enddo

do k=-1,ndz
do j=-1,ndy
do i=-1,ndx
  Phiwvpre(i,j,k,1)=Phiwv(i,j,k,1)
  !Phiwvpre(i,j,k,2)=Phiwv(i,j,k,2)
  Phigrdwvpre(i,j,k,1)=Phigrdwv(i,j,k,1)
  !Phigrdwvpre(i,j,k,2)=Phigrdwv(i,j,k,2)
enddo
enddo
enddo


        do k=1,ndz-2
           do j=1,ndy-2
              do i=1,ndx-2
Phiwv(i,j,k,1) = 0.5d0*Phiwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa * 0.5d0 * dt))+Phigrdwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa * 0.5d0 * dt))/(2.d0*kappa+1.d-10)

Phigrdwv(i,j,k,1) = 0.5d0*kappa*Phiwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa *0.5d0* dt))+0.5d0*Phigrdwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa *0.5d0* dt))
              enddo
           enddo
         enddo



 !iwx=1;iwy=1;iwz=1
 !call BCgrv(100,1,8)
 iwx=1;iwy=0;iwz=0
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)


 iwx=0;iwy=0;iwz=1
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

 iwx=0;iwy=1;iwz=0
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

 !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%






  !****************slv-wv****************
  !%%%%%%%%%%%%%%%%%phi(t+0.5*dt)%%%%%%%%%%%%%%%%%%
  iwx=0;iwy=0;iwz=1
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

  iwx=1;iwy=0;iwz=0
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

  iwx=0;iwy=1;iwz=0
  call BCgrv(100,1,8)
  call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

  call BCgrv(110,1,8)
  call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

do k=-1,ndz
do j=-1,ndy
do i=-1,ndx
  Phiwvpre(i,j,k,1)=Phiwv(i,j,k,1)
  !Phiwvpre(i,j,k,2)=Phiwv(i,j,k,2)
  Phigrdwvpre(i,j,k,1)=Phigrdwv(i,j,k,1)
  !Phigrdwvpre(i,j,k,2)=Phigrdwv(i,j,k,2)
enddo
enddo
enddo



  do k=1,ndz-2
     do j=1,ndy-2
        do i=1,ndx-2
Phiwv(i,j,k,1) = 0.5d0*Phiwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa * 0.5d0 * dt))+Phigrdwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa * 0.5d0 * dt))/(2.d0*kappa+1.d-10)

Phigrdwv(i,j,k,1) = 0.5d0*kappa*Phiwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa *0.5d0* dt))+0.5d0*Phigrdwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa *0.5d0* dt))
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


     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0*2.d0*kappa * dt)+(Phiwv(i,j,k,1) &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdxy1/dx1/dy1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdyz1/dy1/dz1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdzx1/dz1/dx1 &
     !-cg*cg*4.d0/2.d0/kappa/2.d0/kappa*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0*2.d0*kappa * dt))

     Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
     (-2.d0*cg*cg*grdxy1/dx1/dy1 &
      -2.d0*cg*cg*grdyz1/dy1/dz1 &
      -2.d0*cg*cg*grdzx1/dz1/dx1) *dt &
     -G4pi*cg*cg*rho(i,j,k)*dt


     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1) +&
     !(-2.d0*2.d0/2.d0/kappa*cg*cg*grdxy1/dx1/dy1 &
     ! -2.d0*2.d0/2.d0/kappa*cg*cg*grdyz1/dy1/dz1 &
     ! -2.d0*2.d0/2.d0/kappa*cg*cg*grdzx1/dz1/dx1) *dt

     !Phigrdwv(i,j,k,1) = Phigrdwv(i,j,k,1)*dexp(-0.5d0*2.d0*kappa * dt)+(Phiwv(i,j,k,1) &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdxy1mn/dx1/dy1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdyz1mn/dy1/dz1 &
     !-2.d0*2.d0/2.d0/kappa*2.d0/2.d0/kappa*cg*cg*grdzx1mn/dz1/dx1 &
     !-cg*cg*4.d0/2.d0/kappa/2.d0/kappa*G4pi*rho(i,j,k))*(1.d0-dexp(-0.5d0*2.d0*kappa * dt))

            enddo
          enddo
        enddo

do k=-1,ndz
do j=-1,ndy
do i=-1,ndx
  Phiwvpre(i,j,k,1)=Phiwv(i,j,k,1)
  !Phiwvpre(i,j,k,2)=Phiwv(i,j,k,2)
  Phigrdwvpre(i,j,k,1)=Phigrdwv(i,j,k,1)
  !Phigrdwvpre(i,j,k,2)=Phigrdwv(i,j,k,2)
enddo
enddo
enddo


        do k=1,ndz-2
           do j=1,ndy-2
              do i=1,ndx-2
Phiwv(i,j,k,1) = 0.5d0*Phiwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa * 0.5d0 * dt))+Phigrdwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa * 0.5d0 * dt))/(2.d0*kappa+1.d-10)

Phigrdwv(i,j,k,1) = 0.5d0*kappa*Phiwvpre(i,j,k,1)*(1.d0-dexp(-2.d0*kappa *0.5d0* dt))+0.5d0*Phigrdwvpre(i,j,k,1)*(1.d0+dexp(-2.d0*kappa *0.5d0* dt))
              enddo
           enddo
         enddo



 !iwx=1;iwy=1;iwz=1
 !call BCgrv(100,1,8)
 iwx=0;iwy=1;iwz=0
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)


 iwx=1;iwy=0;iwz=0
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)

 iwx=0;iwy=0;iwz=1
 call BCgrv(100,1,8)
 call muslcslv1D(Phiwv(-1,-1,-1,1),dt  *0.5d0,1)

 call BCgrv(110,1,8)
 call muslcslv1D(Phigrdwv(-1,-1,-1,1),dt  *0.5d0,2)



     iwx=1;iwy=1;iwz=1
     call BCgrv(100,1,8)
     call BCgrv(110,1,8)
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
  IF(IST.eq.NSPLTx-1) THEN
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = Ncellx+1, Ncellx+N_ol
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = Ncellx, Ncellx+N_ol
     !Phiwv(IX,JY,KZ,idm)= bphir(JY,KZ,IX)
     Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
     END DO;END DO;END DO
  END IF
  enddo
CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF

IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
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
   IF(JST.eq.NSPLTy-1) THEN
        DO KZ = -1, Ncellz+2; DO JY = Ncelly+1, Ncelly+N_ol; DO IX = -1, Ncellx+2
        Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
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
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
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
  IF(KST.eq.NSPLTz-1) THEN
     DO KZ =  Ncellz+1, Ncellz+N_ol; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
     Phiwv(IX,JY,KZ,idm)= Phiexa(IX,JY,KZ)
     END DO;END DO;END DO
  END IF
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

  IF(IST.eq.0) THEN
     DO KZ = -1, Ncellz+2; DO JY = -1, Ncelly+2; DO IX = 1-N_ol, 0
     !DO KZ = 1, Ncellz; DO JY = 1, Ncelly; DO IX = 1-N_ol, 1
     !Phigrdwv(IX,JY,KZ,idm)= bphigrdxl(JY,KZ,IX,idm)
     Phigrdwv(IX,JY,KZ,idm) = cg*Phigrd(IX,JY,KZ,idm)+kappa*Phiexa(IX,JY,KZ)
     END DO;END DO;END DO
  END IF
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
     Phigrdwv(IX,JY,KZ,idm)= cg*Phigrd(IX,JY,KZ,idm)+kappa*Phiexa(IX,JY,KZ)
     END DO;END DO;END DO
  END IF
  enddo

CALL MPI_TYPE_FREE(VECU,IERR)
LEFT = LEFTt; RIGT = RIGTt
END IF


IF(iwy.EQ.1) THEN
  CALL MPI_TYPE_VECTOR(Ncellz+4,N_ol*(ndx+2),(ndx+2)*(ndy+2),MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  BOTMt = BOTM !; IF(JST.eq.0       ) BOTM = MPI_PROC_NULL
  TOPt  = TOP  !; IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-1,Ncelly+1-N_ol,-1,idm),1,VECU,TOP ,1, &
        Phigrdwv(-1,       1-N_ol,-1,idm),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
   IF(JST.eq.0) THEN
        DO KZ = -1, Ncellz+2; DO JY = 1-N_ol, 0; DO IX = -1, Ncellx+2
        Phigrdwv(IX,JY,KZ,idm)= cg*Phigrd(IX,JY,KZ,idm)+kappa*Phiexa(IX,JY,KZ)
        END DO;END DO;END DO
     END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-1,1            ,-1,idm),1,VECU,BOTM,1, &
        Phigrdwv(-1,Ncelly+1     ,-1,idm),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
   IF(JST.eq.NSPLTy-1) THEN
        DO KZ = -1, Ncellz+2; DO JY = Ncelly+1, Ncelly+N_ol; DO IX =  -1, Ncellx+2
        !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
        Phigrdwv(IX,JY,KZ,idm)= cg*Phigrd(IX,JY,KZ,idm)+kappa*Phiexa(IX,JY,KZ)
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
  DOWNt = DOWN !; IF(KST.eq.0       ) DOWN = MPI_PROC_NULL
  UPt   = UP   !; IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
!*************************************  BC for the downsides of domains  ****
   do idm=is,ie
  CALL MPI_SENDRECV(Phigrdwv(-1,-1,Ncellz+1-N_ol,idm),1,VECU,UP  ,1, &
  Phigrdwv(-1,-1,       1-N_ol,idm),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
   IF(KST.eq.0) THEN
        DO KZ = 1-N_ol, 0; DO JY = -1, Ncelly+2; DO IX = -1, Ncellx+2
        !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
        Phigrdwv(IX,JY,KZ,idm)= cg*Phigrd(IX,JY,KZ,idm)+kappa*Phiexa(IX,JY,KZ)
        END DO;END DO;END DO
     END IF
   enddo
!**************************************  BC for the upsides of domains  ****
   do idm=is,ie
   CALL MPI_SENDRECV(Phigrdwv(-1,-1,1            ,idm),1,VECU,DOWN,1, &
   Phigrdwv(-1,-1,Ncellz+1     ,idm),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
  IF(KST.eq.NSPLTz-1) THEN
       DO KZ =  Ncellz+1, Ncellz+N_ol; DO JY = -1, Ncelly+2; DO IX =  -1, Ncellx+2
       !Phigrdwv(IX,JY,KZ,idm)= Phigrd(IX,JY,KZ,idm)
       Phigrdwv(IX,JY,KZ,idm)= cg*Phigrd(IX,JY,KZ,idm)+kappa*Phiexa(IX,JY,KZ)
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
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,1,is,ie)
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
     !call fluxcal(Phi2dt,Phipre,Phiu,1.0d0,0.0d0,4,is,ie)

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


subroutine fluxcal(preuse,pre,uin,ep,kappa1,mode,is,ie)
  use comvar
  double precision :: ep , kappa1
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
                   * ((1.0d0-slop(i)*kappa1)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
                   (1.0d0+slop(i)*kappa1)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i+1/2
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
                   * ((1.0d0+slop(i)*kappa1)*(pre(ix,jy,kz)-pre(ixm,jym,kzm)) + &
                   (1.0d0-slop(i)*kappa1)*(pre(ixp,jyp,kzp) - pre(ix,jy,kz))) !i-1/2
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


subroutine move()
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'
double precision :: rsph3,tint=0.d0,dt,dinit1


   !********************sphere***********************
     dinit1=ddd
     write(*,*) ddd

     rrsph3x=rrsph3x+dx1
     !rrsph3y=rrsph3y+dy1
     !rrsph3z=rrsph3z+dz1

    do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2

      rsph3 =dsqrt((rrsph3x-x(i))**2+(rrsph3y-y(j))**2+(rrsph3z-z(k))**2)

     if(rsph3 .le. rrsph3 ) then
        !write(*,*) 'inside', rsph3,rrsph3,dinit1
        U(i,j,k,1) = dinit1
        !U(i,j,k,1) = dinit1*(-6.d0/rrsph3**2.d0 + 4.d0*rsph3**2.d0/rrsph3**4.d0)*dexp(-rsph3**2.d0/rrsph3**2.d0)
     else
        !write(*,*) 'outside', rsph3,rrsph3,dinit1
        U(i,j,k,1) = 0.d0
        !U(i,j,k,1) = dinit1*(-6.d0/rrsph3**2.d0 + 4.d0*rsph3**2.d0/rrsph3**4.d0)*dexp(-rsph3**2.d0/rrsph3**2.d0)
     end if
  end do
  end do
  end do

  !do i=1,Ncellx/2,-1
  !  call PBini(i)
  !enddo

  do k = -1-1, Ncellz+2+1; do j = -1-1, Ncelly+2+1; do i = -1-1, Ncellx+2+1
    !rsph3 =dsqrt( (ql1x-x_i(i2)-dx1/2.d0)**2 + (ql1y-y_i(i2y)-dy1/2.d0)**2 + (ql1z-z_i(i2z)-dz1/2.d0)**2 )
    rsph3 =dsqrt((rrsph3x-x(i))**2+(rrsph3y-y(j))**2+(rrsph3z-z(k))**2)
    if(rsph3 .le. rrsph3 ) then
        !Phiexa(i,j,k)=G4pi/6.d0*dinit1*(rsph3*dx1)**2
        Phiexa(i,j,k)=G4pi/6.d0*dinit1*(rsph3)**2
      !  Phiexa(i,j,k)=dinit1*dexp(-rsph3**2.d0/rrsph3**2.d0)*G4pi
        !U(i,j,k,1) = dinit1
        !write(*,*) 'in'
     else
     !Phiexa(i,j,k)=-G4pi/rsph3/dx1/3.d0*dinit1*(rrsph3*dx1)**3+G4pi/2.d0*dinit1*(rrsph3*dx1)**2
     Phiexa(i,j,k)=-G4pi/rsph3/3.d0*dinit1*(rrsph3)**3+G4pi/2.d0*dinit1*(rrsph3)**2
      !Phiexa(i,j,k)=dinit1*dexp(-rsph3**2.d0/rrsph3**2.d0)*G4pi
        !U(i,j,k,1) = 0.d0
     end if
  end do
  end do
  end do


  !call collect()

  !do i=0,-(Ncellx/2-1),-1
  !do pls=0,Ncellx*NSPLTx-1
  !  call PBini(pls)
  !enddo
  !write(*,*) Phiexa(1,1,1)

  do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
     Phigrd(i,j,k,1)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,2)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,3)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,4)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,5)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,6)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,7)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
     Phigrd(i,j,k,8)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                      +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
  end do
  end do
  end do
end subroutine move

subroutine movesph(dt)
USE comvar
USE mpivar
USE chmvar
USE slfgrv
INCLUDE 'mpif.h'
double precision :: rsph3,tint=0.d0,dt

  tint=tint+dt
  !rrsph3 = ql1x * rratio
  !dinit1=dinit1*((ql1x*0.2d0)**3.d0)/(rrsph3**3.d0)
  !rrsph3 = ql1x*rratio
  !dinit1=dinit1*((ql1x*0.2d0)**3.d0)/(rrsph3**3.d0)
  !dinit1=dinit1*4.d0*3.1415926536d0/3.d0*rrsph3**3.d0/(ql1x+ql2x)**3.d0
  rrsph3x=0.5d0+rmove*dsin(0.d0+vmove*tint)
  rrsph3y=0.5d0+rmove*dcos(0.d0+vmove*tint)
  rrsph3z=0.5d0

  write(*,*)'MOVE',rrsph3x,rrsph3y,rrsph3z,tint

  do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   i2 =  IST*Ncellx+i
   i2y = JST*Ncelly+j
   i2z = KST*Ncellz+k
   !rsph=dsqrt( (cenx-dble(i2))**2 + (ceny-dble(i2y))**2 + (cenz-dble(i2z))**2 )
   !rsph=dsqrt( (cenx-dble(i2))**2 + (ceny-dble(i2y))**2 + (cenz-dble(i2z))**2 )
   !rsph3 =dsqrt( (rrsph3x-x_i(i2))**2 + (rrsph3y-y_i(i2y))**2 + (rrsph3z-z_i(i2z))**2 )
   rsph3 =dsqrt((rrsph3x-x(i))**2+(rrsph3y-y(j))**2+(rrsph3z-z(k))**2)
   !rrsph3 = ql1x*rratio
   

  !Mcnst=4.d0*3.14159265358979d0/3.d0*dinit1*(ql1x*0.2d0)**3.d0
  !dinit1=dinit1*((ql1x*0.2d0)**3.d0)/(rrsph3**3.d0)

   if(rsph3 .le. rrsph3 ) then
      !write(*,*) 'inside', rsph3,rrsph3,dinit1
      U(i,j,k,1) = ddd
!      U(i,j,k,2) = 0.0d0
!      U(i,j,k,3) = 0.0d0
!      U(i,j,k,4) = 0.0d0
!      U(i,j,k,5) = pinit1
!      U(i,j,k,6) = 0.0d0
!      U(i,j,k,7) = 0.0d0
!      U(i,j,k,8) = 0.0d0
!      ndH(i,j,k)   = Hini
!      ndp(i,j,k)   = pini
!      ndH2(i,j,k)  = H2ini
!      ndHe(i,j,k)  = Heini
!      ndHep(i,j,k) = Hepini
!      ndC(i,j,k)   = Cini
!      ndCO(i,j,k)  = COini
!      ndCp(i,j,k)  = Cpini
!      nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
!      ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
!      Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
!      Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
   else
      !write(*,*) 'outside', rsph3,rrsph3,dinit1
      U(i,j,k,1) = 0.0d0
!      U(i,j,k,2) = 0.0d0
!      U(i,j,k,3) = 0.0d0
!      U(i,j,k,4) = 0.0d0
!      U(i,j,k,5) = 0.0d0
!      U(i,j,k,6) = 0.0d0
!      U(i,j,k,7) = 0.0d0
!      U(i,j,k,8) = 0.0d0
!      ndH(i,j,k)   = 0.0d0
!      ndp(i,j,k)   = 0.0d0
!      ndH2(i,j,k)  = 0.0d0
!      ndHe(i,j,k)  = 0.0d0
!      ndHep(i,j,k) = 0.0d0
!      ndC(i,j,k)   = 0.0d0
!      ndCO(i,j,k)  = 0.0d0
!      ndCp(i,j,k)  = 0.0d0
!      nde(i,j,k)   = ndp(i,j,k)+ndHep(i,j,k)+ndCp(i,j,k)
!      ndtot(i,j,k) = ndH(i,j,k)+ndp(i,j,k)+2.d0*ndH2(i,j,k)+ndHe(i,j,k)+ndHep(i,j,k)
!      Ntot(i,j,k,1)=0.d0; NH2(i,j,k,1)=0.d0; NnC(i,j,k,1)=0.d0; tCII(i,j,k,1)=0.d0
!      Ntot(i,j,k,2)=0.d0; NH2(i,j,k,2)=0.d0; NnC(i,j,k,2)=0.d0; tCII(i,j,k,2)=0.d0
   end if
end do
end do
end do

!do i=1,Ncellx/2,-1
!  call PBini(i)
!enddo

do k = -1-1, Ncellz+2+1; do j = -1-1, Ncelly+2+1; do i = -1-1, Ncellx+2+1
   !i2 = IST*Ncellx+i
   !i2y = JST*Ncelly+j
   !i2z = KST*Ncellz+k
   !rrsph3x=50.d0+rmove*dsin(0.d0+tint)
   !rrsph3y=50.d0+rmove*dcos(0.d0+tint)
   !rrsph3z=50.d0
   !rsph=dsqrt( (cenx-dble(i2))**2 + (ceny-dble(i2y))**2 + (cenz-dble(i2z))**2 )
   rsph3 =dsqrt((rrsph3x-x(i))**2+(rrsph3y-y(j))**2+(rrsph3z-z(k))**2)
   !rsph3 =dsqrt( (ql1x-x_i(i2))**2 + (ql1y-y_i(i2y))**2 + (ql1z-z_i(i2z))**2 )
   !rsph3 =dsqrt( (ql1x+0.5d0*dx1-x_i(i2))**2 + (ql1y+0.5d0*dy1-y_i(i2y))**2 + (ql1z+0.5d0*dz1-z_i(i2z))**2 )
   !rrsph3 = ql1x*0.4d0!+dx1*0.5d0
   !rrsph3 = ql1x*rratio

   !write(*,*)rsph3,rrsph3,x_i(i2),y_i(i2y),z_i(i2z)
   if(rsph3 .le. rrsph3 ) then
      !Phiexa(i,j,k)=G4pi/6.d0*dinit1*(rsph3*dx1)**2
   Phiexa(i,j,k)=G4pi/6.d0*ddd*(rsph3)**2.d0
      !U(i,j,k,1) = dinit1
      !write(*,*) 'in'
   else
   !Phiexa(i,j,k)=-G4pi/rsph3/dx1/3.d0*dinit1*(rrsph3*dx1)**3+G4pi/2.d0*dinit1*(rrsph3*dx1)**2
    Phiexa(i,j,k)=-G4pi/rsph3/3.d0*ddd*(rrsph3)**3+G4pi/2.d0*ddd*(rrsph3)**2
    !U(i,j,k,1) = 0.d0
   end if
end do
end do
end do


!call collect()

!do i=0,-(Ncellx/2-1),-1
!do pls=0,Ncellx*NSPLTx-1
!  call PBini(pls)
!enddo
!write(*,*) Phiexa(1,1,1)

do k = -1, Ncellz+2; do j = -1, Ncelly+2; do i = -1, Ncellx+2
   Phigrd(i,j,k,1)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,2)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,3)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,4)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1+(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,5)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,6)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,7)= (-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    -(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
   Phigrd(i,j,k,8)=-(-Phiexa(i-1,j,k)+Phiexa(i+1,j,k))*0.5d0/dx1 &
                    +(-Phiexa(i,j-1,k)+Phiexa(i,j+1,k))*0.5d0/dy1-(-Phiexa(i,j,k-1)+Phiexa(i,j,k+1))*0.5d0/dz1
end do
end do
end do


!dinit1=0.0d0
 !6001 continue
end subroutine movesph
