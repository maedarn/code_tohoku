SUBROUTINE GRAVTY(dt,mode)
USE comvar
USE mpivar
USE slfgrv
!implicit none
INCLUDE 'mpif.h'
DOUBLE PRECISION  :: dt,dxi
INTEGER :: LEFTt,RIGTt,TOPt,BOTMt,UPt,DOWNt,mode,N_ol,i,j,k
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
DOUBLE PRECISION  :: VECU

if(mode.eq.1) then
  call pinter(Nmem1,Nmem2,Ncellx,Ncelly,Ncellz)
end if

if(mode.eq.2) then
  N_MPI(20)=1; N_MPI(1)=1
  iwx = 1; iwy = 1; iwz = 1; CALL BC_MPI(1,1)
  call PB()
  call mglin(Nmem1,Nmem2,2,5,5)
  DEALLOCATE(bphi1,bphi2)

  N_ol = 2
                      !count, blocklength, stride
  CALL MPI_TYPE_VECTOR((ndy+2)*(Ncellz+4),N_ol,ndx+2,MPI_REAL8,VECU,IERR)
  CALL MPI_TYPE_COMMIT(VECU,IERR)
  LEFTt = LEFT; IF(IST.eq.0       ) LEFT = MPI_PROC_NULL
  RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGT = MPI_PROC_NULL
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

  if(IST.eq.0       ) then; do k=1,Ncellz; do j=1,Ncelly
    Phi(0       ,j,k) = Phi(1     ,j,k); Phi(-1       ,j,k) = Phi(1     ,j,k) !grad=0
  end do; end do; end if
  if(IST.eq.NSPLTx-1) then; do k=1,Ncellz; do j=1,Ncelly
    Phi(Ncellx+1,j,k) = Phi(Ncellx,j,k); Phi(Ncellx+2,j,k) = Phi(Ncellx,j,k)
  end do; end do; end if
end if

if(mode.eq.3) then !acceraration because of gravity
  dxi = 1.d0/(12.d0*dx(0))
  do k=1,Ncellz; do j=1,Ncelly; do i=1,Ncellx
    U(i,j,k,2) = U(i,j,k,2) - dt * ( -Phi(i+2,j,k)+8.d0*Phi(i+1,j,k)-8.d0*Phi(i-1,j,k)+Phi(i-2,j,k) ) * dxi *0.5d0
    U(i,j,k,3) = U(i,j,k,3) - dt * ( -Phi(i,j+2,k)+8.d0*Phi(i,j+1,k)-8.d0*Phi(i,j-1,k)+Phi(i,j-2,k) ) * dxi *0.5d0
    U(i,j,k,4) = U(i,j,k,4) - dt * ( -Phi(i,j,k+2)+8.d0*Phi(i,j,k+1)-8.d0*Phi(i,j,k-1)+Phi(i,j,k-2) ) * dxi *0.5d0
  end do;end do;end do
end if

END SUBROUTINE GRAVTY

SUBROUTINE pinter(Need1,Need2,Ncellx,Ncelly,Ncellz)
USE slfgrv
USE mpivar
integer nl,nx,ny,nz,Ncellx,Ncelly,Ncellz,Need1,Need2

!***  finest grid pointer : point(NGL) ***
!*** coasest grid pointer : point(0)   ***

NGL = min0(Ncellx*NSPLTx,Ncelly*NSPLTy,Ncellz*NSPLTz)
NGL = int(dlog(dble(NGL))/dlog(2.d0)+1.d-3)


!*** Unsplit Pointer ***!
NGcr = max0(NSPLTx,NSPLTy,NSPLTz)
NGcr = int(dlog(dble(NGcr))/dlog(2.d0)+1.d-3) + 1

point1(1) = 1
nl = 1
nx=3; ny=3; nz=3
2 continue
point1(nl+1)=point1(nl)+(nx)*(ny)*(nz)
nx=nx*2-1; ny=ny*2-1; nz=nz*2-1
nl = nl+1
if(nl.ne.NGcr+1) goto 2

Need1 = point1(NGcr+1)

!*** MPI Split Pointer ***!
nl = NGcr-1
point2(nl) = 1
nx=(2**NGcr)/NSPLTx+2; ny=(2**NGcr)/NSPLTy+2; nz=(2**NGcr)/NSPLTz+2
3 continue
point2(nl+1)=point2(nl)+(nx)*(ny)*(nz)
nx=nx*2-2; ny=ny*2-2; nz=nz*2-2
write(*,*) point2(nl+1)
nl = nl+1
if(nl.ne.NGL+1) goto 3

Need2 = point2(NGL+1)

if(NRANK.eq.0) write(*,*) 'NGL=',NGL
if(NRANK.eq.0) write(*,*) 'NGcr=',NGcr
if(NRANK.eq.0) write(*,*) 'need1=',Need1
if(NRANK.eq.0) write(*,*) 'need2=',Need2

END SUBROUTINE pinter

!***********************************
! set BC in interp (except in addint) and the first slvsml when fixed boundary except phi=0
!
SUBROUTINE mglin(Need1,Need2,ncycle,NPRE,NPOST)
USE comvar
USE mpivar
USE slfgrv
INCLUDE 'mpif.h'
DOUBLE PRECISION :: cphi1(Need1),crho1(Need1),cres1(Need1),crhs1(Need1) !Unsplit
DOUBLE PRECISION :: cphi2(Need2),crho2(Need2),cres2(Need2),crhs2(Need2) !MPI Split
DOUBLE PRECISION :: tMPI(Need1,0:NPE-1)
integer :: debug1=0,ngrid,ncx,ncy,ncz,ncxcr,ncycr,nczcr,ncycle,jcycle,nfx,nfy,nfz,NPRE,NPOST
integer kk,ii,jj,jpre,jpost,Need1,Need2,mode,i,j,k,testi,testj,testk
!Pre-BC for rho is necessary
!if(debug1==0) then
!   open(502,file='inirhs.dat')
!endif

do k=0,Ncellz+1; kk=(Ncellx+2)*(Ncelly+2)*k+point2(NGL)
do j=0,Ncelly+1; jj=(Ncellx+2)*j+kk
do i=0,Ncellx+1; nc = i+jj
   crhs2(nc) = U(i,j,k,1)*G4pi
  ! write(502,*) crhs2(nc)
end do;end do;end do

!if(debug1==0) then
!   close(502)
!   debug=5
!endif

ncx=(Ncellx+1)/2+1; ncy=(Ncelly+1)/2+1; ncz=(Ncellz+1)/2+1
ngrid=NGL-1

call rstrctMPI(crho2(point2(ngrid)),crhs2(point2(ngrid+1)),ncx,ncy,ncz,0)

do while(ngrid.ne.NGcr)
  ncx=ncx/2+1; ncy=ncy/2+1; ncz=ncz/2+1
  ngrid=ngrid-1
  call rstrctMPI(crho2(point2(ngrid)),crho2(point2(ngrid+1)),ncx,ncy,ncz,0)
end do
ncxcr=ncx; ncycr=ncy; nczcr=ncz

ncx=2**NGcr+1;ncy=ncx;ncz=ncx
!write(*,*) ncx,ncxcr,'nnnnnnnnccccccccccccxxxxxxxxx'
call collect( crho1(point1(NGcr)),crho2(point2(NGcr)),ncx,ncy,ncz,ncxcr,ncycr,nczcr )

do while(ngrid.ne.1)
  ncx=ncx/2+1; ncy=ncy/2+1; ncz=ncz/2+1
  ngrid=ngrid-1
  call rstrct(crho1(point1(ngrid)),crho1(point1(ngrid+1)),ncx,ncy,ncz,0)
end do

call slvsmlb(cphi1(point1(1)),crho1(point1(1))) !BC set is necessary

!Here nc=3
ngrid = NGL
do j=2,ngrid

  IF(j.le.NGcr) THEN !*** generate candidate sol. from j-1 to j (upward) *** 
     ncx=ncx*2-1; ncy=ncy*2-1; ncz=ncz*2-1
    ! write(*,*) NRANK,j,cphi1(point1(j)),cphi1(point1(j-1)),jcycle,jj ,'pointphi2-1'
     call interp(cphi1(point1(j)),cphi1(point1(j-1)),ncx,ncy,ncz,pointb1(j),1)  !BC set is necessary
     ! write(*,*) NRANK,j,cphi1(point1(j)),cphi1(point1(j-1)),crhs1(point1(j)),crho1(point1(j)),jcycle,jj ,'pointphi2-2'
      call copy(crhs1(point1(j)),crho1(point1(j)),ncx,ncy,ncz)
     !  write(*,*) NRANK,j,crhs1(point1(j)),crho1(point1(j)),jcycle,jj ,'pointphi2-3'
    if(j.eq.NGcr) then
       ncx=ncxcr; ncy=ncycr; ncz=nczcr
      !  write(*,*) NRANK,j,cphi2(point2(NGcr)),cphi1(point1(NGcr)),jcycle,jj ,'pointphi2-4'
        call divide( cphi2(point2(NGcr)),cphi1(point1(NGcr)),ncx,ncy,ncz,2**NGcr+1,2**NGcr+1,2**NGcr+1 )
       ! write(*,*) NRANK,j,cphi2(point2(NGcr)),cphi1(point1(NGcr)), cphi2(point2(NGcr)),cphi1(point1(NGcr)),jcycle,jj ,'pointphi2-5'
        call copyMPI(crhs2(point2(NGcr)),crho2(point2(NGcr)),ncx,ncy,ncz)
        ! write(*,*) NRANK,j,cphi2(point2(NGcr)),cphi1(point1(NGcr)),jcycle,jj ,'pointphi2-6'
    end if
  ELSE
     ncx=ncx*2-1; ncy=ncy*2-1; ncz=ncz*2-1
     ! write(*,*) NRANK,j,cphi2(point2(j)),cphi2(point2(j-1)),pointb2(j),jcycle,jj ,'pointphi2-7'
      call interpMPI(cphi2(point2(j)),cphi2(point2(j-1)),ncx,ncy,ncz,pointb2(j),1)  !BC set is necessary
      ! write(*,*) NRANK,j,cphi2(point2(j)),cphi2(point2(j-1)),pointb2(j),jcycle,jj ,'pointphi2-8'
       if(j.ne.ngrid) call copyMPI(crhs2(point2(j)),crho2(point2(j)),ncx,ncy,ncz)
      !  write(*,*) NRANK,j,jcycle,crhs2(point2(j)),crho2(point2(j)),jj ,'pointphi2-9'
  END IF
  !write(*,*) NRANK,j,Phi(33,33,33) ,'pointphi1'



  !*********************dbug***********************
  goto 10111
  if(NRANK==40) then
     if(j.le.NGcr) then
        do testi=1,ncx**3
           write(*,*) cphi1(point1(j)+testi),cphi1(point1(j-1)+testi/8+1),crhs1(point1(j)+testi),crho1(point1(j)+testi),ncx,j,testi,'phipint1-1'
        end do
        if(j.eq.NGcr) then
           do testi=1,ncx**3
              write(*,*) cphi2(point2(NGcr)+testi),cphi1(point1(NGcr)+testi), cphi2(point2(NGcr)+testi),cphi1(point1(NGcr)+testi),ncx,j,testi,'phipint1-2'
           end do
        end if
     else
        do testi=1,ncx**3
           write(*,*) cphi2(point2(j)+testi),cphi2(point2(j-1)+testi/8+1),crhs2(point2(j)+testi),crho2(point2(j)+testi),ncx,j,testi,'phipint1-3'
        end do
     end if
  end if
10111 continue
  !*********************dbug***********************





  do jcycle=1,ncycle !V-cycle

    nfx=ncx; nfy=ncy; nfz=ncz
    do jj=j,2,-1          !*** DOWNWARD *****************************
                          !    phi + rhs --> res -->      (level N  )
                          !                          rhs  (level N-1)
       IF(jj.lt.NGcr) THEN !*** generate residual from jj to jj-1 ****
      !    write(*,*) NRANK,j,jcycle,cphi1(point1(jj)),crhs1(point1(jj)),jj ,'pointphi2-10'
        do jpre=1,NPRE
          mode=2; if((jj.ne.j).and.(jpre.eq.1)) mode=1
          call relax(cphi1(point1(jj)),crhs1(point1(jj)),nfx,nfy,nfz,mode)
       end do
 !       write(*,*) NRANK,j,jcycle,cphi1(point1(jj)),crhs1(point1(jj)),cres1(point1(jj)),jj ,'pointphi2-11'
        call resid(cres1(point1(jj)),cphi1(point1(jj)),crhs1(point1(jj)),nfx,nfy,nfz)
  !      write(*,*) NRANK,j,jcycle,cphi1(point1(jj)),crhs1(point1(jj)),cres1(point1(jj)),jj ,'pointphi2-12'
        nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
        call rstrct(crhs1(point1(jj-1)),cres1(point1(jj)),nfx,nfy,nfz,1)  !fill0 at BC below this subroutine is necessary
   !     write(*,*) NRANK,j,jcycle,cphi1(point1(jj)),crhs1(point1(jj)),cres1(point1(jj)),jj ,'pointphi2-13'
        call  fill0(cphi1(point1(jj-1)),nfx,nfy,nfz)
    !    write(*,*) NRANK,j,jcycle,cphi1(point1(jj)),crhs1(point1(jj)),cres1(point1(jj)),jj ,'pointphi2-14'
     ELSE
         NPRE1 = NPRE; if(j.ge.NGL) NPRE1 = 2
     !    write(*,*) NRANK,j,jcycle,cphi2(point2(jj)),crhs2(point2(jj)),jj ,'pointphi2-15'
        do jpre=1,NPRE1
          mode=2; if((jj.ne.j).and.(jpre.eq.1)) mode=1
          call relaxMPI(cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz,mode)
       end do
      ! write(*,*) NRANK,j,cres2(point2(jj)),cphi2(point2(jj)),crhs2(point2(jj)),jcycle,jj ,'pointphi2**************'
       call residMPI(cres2(point2(jj)),cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz)
       !write(*,*) NRANK,j,cres2(point2(jj)),cphi2(point2(jj)),crhs2(point2(jj)),jcycle,jj ,'pointphi3'
        if(jj.eq.NGcr) then
           nfx=2**NGcr+1;nfy=nfx;nfz=nfx
        !    write(*,*) NRANK,j,cphi1(point1(NGcr)),cphi2(point2(NGcr)),cres1(point1(NGcr)),cres2(point2(NGcr)),jcycle,jj ,'pointphi3-1'
            call collect( cphi1(point1(NGcr)),cphi2(point2(NGcr)),nfx,nfy,nfz,ncxcr,ncycr,nczcr ) !necessary at upward loop below
         !    write(*,*) NRANK,j,cphi1(point1(NGcr)),cphi2(point2(NGcr)),cres1(point1(NGcr)),cres2(point2(NGcr)),jcycle,jj ,'pointphi3-2'
             call collect( cres1(point1(NGcr)),cres2(point2(NGcr)),nfx,nfy,nfz,ncxcr,ncycr,nczcr )
          !    write(*,*) NRANK,j,cphi1(point1(NGcr)),cphi2(point2(NGcr)),cres1(point1(NGcr)),cres2(point2(NGcr)),jcycle,jj ,'pointphi3-3'
              nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
           !    write(*,*) NRANK,j,crhs1(point1(jj-1)),cres1(point1(jj)),jcycle,jj ,'pointphi3-4'
               call rstrct(crhs1(point1(jj-1)),cres1(point1(jj)),nfx,nfy,nfz,1)  !fill0 at BC below this subroutine is necessary
            !    write(*,*) NRANK,j,cphi1(point1(jj-1)),crhs1(point1(jj-1)),cres1(point1(jj)),jcycle,jj ,'pointphi3-5'
                call  fill0(cphi1(point1(jj-1)),nfx,nfy,nfz)
             !    write(*,*) NRANK,j,cphi1(point1(jj-1)),crhs1(point1(jj-1)),cres1(point1(jj)),jcycle,jj ,'pointphi3-6'
        else
           nfx=nfx/2+1; nfy=nfy/2+1; nfz=nfz/2+1
       !     write(*,*) NRANK,j,crhs2(point2(jj-1)),cres2(point2(jj)),cphi2(point2(jj-1)),jcycle,jj ,'pointphi3-7'
            call rstrctMPI(crhs2(point2(jj-1)),cres2(point2(jj)),nfx,nfy,nfz,1)  !fill0 at BC below this subroutine is necessary
        !     write(*,*) NRANK,j,crhs2(point2(jj-1)),cres2(point2(jj)),cphi2(point2(jj-1)),jcycle,jj ,'pointphi3-8'
             call  fill0MPI(cphi2(point2(jj-1)),nfx,nfy,nfz)
         !     write(*,*) NRANK,j,crhs2(point2(jj-1)),cres2(point2(jj)),cphi2(point2(jj-1)),jcycle,jj ,'pointphi3-9'
        end if
     END IF

     !*********************dbug***********************
 goto    10112
  if(NRANK==40) then
     if(jj.lt.NGcr) then
        do testi=1,(nfx*2-1)**3
           write(*,*) cphi1(point1(jj)+testi),cphi1(point1(jj-1)+testi/8+1),crhs1(point1(jj)+testi),crhs1(point1(jj)-1+testi/8),cres1(point1(jj)+testi),nfx,jj,jcycle,testi,'phipint2-1'
        end do
        write(*,*)
     else
        if(jj.eq.NGcr) then
           do testi=1,(nfx*2-1)**3
              write(*,*) cphi1(point1(jj)+testi),cphi2(point1(jj)+testi),crhs1(point1(jj-1)+testi/8+1),cres1(point1(jj)+testi),cres2(point1(jj)+testi),nfx,jj,jcycle,testi,'phipint2-2'
           end do
            write(*,*)
        end if
        do testi=1,(nfx*2-1)**3
           write(*,*) cphi2(point2(jj)+testi),cphi2(point2(jj-1)+testi/8+1),crhs2(point2(jj)+testi),crhs2(point1(jj)-1+testi/8),cres2(point2(jj)+testi),nfx,jj,jcycle,testi,'phipint2-3'
        end do
         write(*,*)
     end if
  end if
10112 continue
  !*********************dbug***********************

end do
 !write(*,*) NRANK,j,cphi1(point1(1)),crhs1(point1(1)),jcycle,jj ,'pointphi3-10'
    call slvsml(cphi1(point1(1)),crhs1(point1(1)))  !BC set is unnecessary
!write(*,*) NRANK,j,cphi1(point1(1)),crhs1(point1(1)),jcycle,jj ,'pointphi3-11'
    nfx=3; nfy=3; nfz=3
    do jj=2,j             !*** UPWARD **********************************
                          !    phi --> phi +rhs --> phi      (level N  )
                          !  + phi                           (level N-1)
      IF(jj.le.NGcr) THEN !*** generate new solution from jj-1 to jj ***
         nfx=2*nfx-1; nfy=2*nfy-1; nfz=2*nfz-1
 !        write(*,*) NRANK,j,cphi1(point1(jj)),cphi1(point1(jj-1)),cres1(point1(jj)),jcycle,jj ,'pointphi3-12'
         call addint(cphi1(point1(jj)),cphi1(point1(jj-1)),cres1(point1(jj)),nfx,nfy,nfz)  !BC set is unnecessary
  !        write(*,*) NRANK,j,cphi1(point1(jj)),cphi1(point1(jj-1)),cres1(point1(jj)),jcycle,jj ,'pointphi3-13'
        if(jj.eq.NGcr) then
           nfx=ncxcr; nfy=ncycr; nfz=nczcr
   !         write(*,*) NRANK,j,cphi2(point2(NGcr)),cphi1(point1(NGcr)),jcycle,jj,nfx,nfy,nfz,2**NGcr+1,2**NGcr+1,2**NGcr+1 ,'pointphi3-14'
            call divide( cphi2(point2(NGcr)),cphi1(point1(NGcr)),nfx,nfy,nfz,2**NGcr+1,2**NGcr+1,2**NGcr+1 )
    !         write(*,*) NRANK,j,cphi2(point2(NGcr)),cphi1(point1(NGcr)),jcycle,cphi2(point2(jj)),crhs2(point2(jj)),jj ,'pointphi3-15'
          do jpost=1,NPOST
            call relaxMPI(cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz,2)
          end do
 !write(*,*) NRANK,j,cphi2(point2(NGcr)),cphi1(point1(NGcr)),jcycle,cphi2(point2(jj)),crhs2(point2(jj)),jj ,'pointphi3-15'
       else
  !        write(*,*) NRANK,j,cphi1(point1(jj)),crhs1(point1(jj)),jcycle,jj ,'pointphi4'
          do jpost=1,NPOST
            call relax(cphi1(point1(jj)),crhs1(point1(jj)),nfx,nfy,nfz,2)
         end do
   !      write(*,*) NRANK,j,cphi1(point1(jj)),crhs1(point1(jj)),jcycle,jj ,'pointphi5'
        end if
      ELSE
         nfx=2*nfx-1; nfy=2*nfy-1; nfz=2*nfz-1
    !      write(*,*) NRANK,j,cphi2(point2(jj)),cphi2(point2(jj-1)),cres2(point2(jj)),jcycle,jj ,'pointphi6'
          call addintMPI(cphi2(point2(jj)),cphi2(point2(jj-1)),cres2(point2(jj)),nfx,nfy,nfz)  !BC set is unnecessary
     !     write(*,*) NRANK,j,cphi2(point2(jj)),cphi2(point2(jj-1)),cres2(point2(jj)),jcycle,jj ,'pointphi7'
          NPOST1 = NPOST; if(j.ge.NGL) NPOST1 = 2
      !    write(*,*) NRANK,j,cphi2(point2(jj)),cphi2(point2(jj)),cres2(point2(jj)),jcycle,jj,nfx,nfy,nfz ,'pointphi7-1'
          do jpost=1,NPOST1
          call relaxMPI(cphi2(point2(jj)),crhs2(point2(jj)),nfx,nfy,nfz,2)
       end do
   !     write(*,*) NRANK,j,cphi2(point2(jj)),cphi2(point2(jj)),cres2(point2(jj)),jcycle,jj ,'pointphi7-2'
    END IF



    !*********************dbug***********************
    goto 10113
  if(NRANK==40) then
     if(jj.le.NGcr) then
        do testi=1,nfx**3
           write(*,*) cphi1(point1(jj)+testi),cphi1(point1(jj-1)+testi/8+1),crhs1(point1(jj)+testi),crhs1(point1(jj)-1+testi/8),cres1(point1(jj)+testi),nfx,jj,jcycle,testi,'phipint3-1'
        end do
        write(*,*)
         if(jj.eq.NGcr) then
           do testi=1,nfx**3
              write(*,*) cphi2(point1(jj)+testi),cphi1(point1(jj)+testi),cphi1(point1(jj-1)+testi/8+1),crhs2(point1(jj)+testi),cres1(point1(jj)+testi),nfx,jj,jcycle,testi,'phipint3-2'
           end do
           write(*,*)
         end if
     else
        do testi=1,nfx**3
           write(*,*) cphi2(point2(jj)+testi),cphi2(point2(jj-1)+testi/8+1),crhs2(point2(jj)+testi),crhs2(point1(jj)-1+testi/8),cres2(point2(jj)+testi),nfx,jj,jcycle,testi,'phipint3-3'
        end do
        write(*,*)
     end if
  end if
10113 continue
  !*********************dbug***********************


  end do

end do
end do

!ncx = Ncellx+1

do k=1,ncz; kk=(ncx+1)*(ncy+1)*k+point2(NGL)
do j=1,ncy; jj=(ncx+1)*j+kk
do i=1,ncx; ii = i+jj
      Phi(i,j,k) = cphi2(ii)
end do; end do; end do

END SUBROUTINE mglin


SUBROUTINE BCsgr_MPI(u,nx,ny,nz,lx1,lx2,ly1,ly2,lz1,lz2)
USE mpivar
INCLUDE 'mpif.h'
integer nx,ny,nz,lx1,lx2,ly1,ly2,lz1,lz2
DOUBLE PRECISION :: u(0:nx,0:ny,0:nz)
DOUBLE PRECISION  :: VECU
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
INTEGER :: LEFTt,RIGTt

!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!*** for X ***!
LEFTt = LEFT; IF(IST.eq.0)        LEFTt = MPI_PROC_NULL
RIGTt = RIGT; IF(IST.eq.NSPLTx-1) RIGTt = MPI_PROC_NULL
CALL MPI_TYPE_VECTOR((ny+1)*(nz+1),1,nx+1,MPI_REAL8,VECU,IERR); CALL MPI_TYPE_COMMIT(VECU,IERR)
if(lx1.eq.1) CALL MPI_SENDRECV(u(nx-1,0,0),1,VECU,RIGTt,1, &
                               u(   0,0,0),1,VECU,LEFTt,1, MPI_COMM_WORLD,MSTATUS,IERR)
if(lx2.eq.1) CALL MPI_SENDRECV(u(   1,0,0),1,VECU,LEFTt,1, &
                               u(nx  ,0,0),1,VECU,RIGTt,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
!*** for Y ***!
!IF(JST.eq.0)        BOTM = MPI_PROC_NULL
!IF(JST.eq.NSPLTy-1) TOP  = MPI_PROC_NULL
CALL MPI_TYPE_VECTOR(nz+1,nx+1,(nx+1)*(ny+1),MPI_REAL8,VECU,IERR); CALL MPI_TYPE_COMMIT(VECU,IERR)
if(ly1.eq.1) CALL MPI_SENDRECV(u(0,ny-1,0),1,VECU,TOP ,1, &
                               u(0,   0,0),1,VECU,BOTM,1, MPI_COMM_WORLD,MSTATUS,IERR)
if(ly2.eq.1) CALL MPI_SENDRECV(u(0,   1,0),1,VECU,BOTM,1, &
                               u(0,ny  ,0),1,VECU,TOP ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
!*** for z ***!
!IF(KST.eq.0)        DOWN = MPI_PROC_NULL
!IF(KST.eq.NSPLTz-1) UP   = MPI_PROC_NULL
CALL MPI_TYPE_VECTOR(1,(nx+1)*(ny+1),(nx+1)*(ny+1),MPI_REAL8,VECU,IERR); CALL MPI_TYPE_COMMIT(VECU,IERR)
if(lz1.eq.1) CALL MPI_SENDRECV(u(0,0,nz-1),1,VECU,UP  ,1, &
                               u(0,0,   0),1,VECU,DOWN,1, MPI_COMM_WORLD,MSTATUS,IERR)
if(lz2.eq.1) CALL MPI_SENDRECV(u(0,0,   1),1,VECU,DOWN,1, &
                               u(0,0,nz  ),1,VECU,UP  ,1, MPI_COMM_WORLD,MSTATUS,IERR)
CALL MPI_TYPE_FREE(VECU,IERR)
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

END SUBROUTINE BCsgr_MPI


SUBROUTINE collect(u1,u2,nx1,ny1,nz1,nx2,ny2,nz2)
USE mpivar
INCLUDE 'mpif.h'
integer nx1,ny1,nz1,nx2,ny2,nz2,nxed,nyed,nzed,kk,jj,ii,Nroot
INTEGER :: MSTATUS(MPI_STATUS_SIZE)
double precision :: u1(nx1,ny1,nz1),u2(0:nx2,0:ny2,0:nz2)
double precision :: tMPI(0:nx2,0:ny2,0:nz2,0:NPE-1)

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
END SUBROUTINE collect


SUBROUTINE divide(u2,u1,nx2,ny2,nz2,nx1,ny1,nz1)
  USE mpivar
  integer nx1,ny1,nz1,nx2,ny2,nz2,iist,jjst,kkst,ii,jj,kk,i,j,k
double precision :: u2(0:nx2,0:ny2,0:nz2),u1(nx1,ny1,nz1)

iist=0; IF(IST.eq.0) iist=1
jjst=0; IF(JST.eq.0) jjst=1
kkst=0; IF(KST.eq.0) kkst=1
do kk=kkst,nz2;k=KST*(nz2-1)+kk
  do jj=jjst,ny2;j=JST*(ny2-1)+jj
    do ii=iist,nx2;i=IST*(nx2-1)+ii
       u2(ii,jj,kk) = u1(i,j,k)
       !if(NRANK==40) then
       !   write(*,*) u2(ii,jj,kk) , u1(i,j,k) ,'divide'
       !end if
end do;end do;end do

END SUBROUTINE divide


SUBROUTINE rstrctMPI(uc,uf,nx,ny,nz,mode)
  USE mpivar
  integer nx,ny,nz,kc,kf,jc,jf,ic,if,izst,ixst,iyst
  integer ized,ixed,iyed
double precision :: uc(0:nx,0:ny,0:nz),uf(0:2*nx-1,0:2*ny-1,0:2*nz-1)
double precision, parameter :: w = 1.d0/12.d0


ixst=1; ixed=nx-1; iyst=1; iyed=ny-1; izst=1; ized=nz-1
IF(IST.eq.0) ixst = 2; IF(JST.eq.0) iyst = 2; IF(KST.eq.0) izst = 2



do kc=izst,ized; kf=2*kc-1
  do jc=iyst,iyed; jf=2*jc-1
    do ic=ixst,ixed; if=2*ic-1
      uc(ic,jc,kc)=0.5d0*uf(if,jf,kf)+ &
      w*(uf(if+1,jf,kf)+uf(if-1,jf,kf)+uf(if,jf+1,kf)+uf(if,jf-1,kf)+uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do


IF(IST.eq.0) THEN
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(1 ,jc,kc)=uf(1 ,jf,kf)
  end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
  nf=2*nx-1
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(nx,jc,kc)=uf(nf,jf,kf)
  end do; end do
END IF
IF(JST.eq.0) THEN
  do kc=1,nz; kf=2*kc-1; do ic=1,nx; if=2*ic-1
    uc(ic,1 ,kc)=uf(if,1 ,kf)
  end do; end do
END IF
IF(JST.eq.NSPLTy-1) THEN
  nf=2*ny-1
  do kc=1,nz; kf=2*kc-1; do ic=1,nx; if=2*ic-1
    uc(ic,ny,kc)=uf(if,nf,kf)
  end do; end do
END IF
IF(KST.eq.0) THEN
  do jc=1,ny; jf=2*jc-1; do ic=1,nx; if=2*ic-1
    uc(ic,jc,1 )=uf(if,jf,1 )
  end do; end do
END IF
IF(KST.eq.NSPLTz-1) THEN
  nf=2*nz-1
  do jc=1,ny; jf=2*jc-1; do ic=1,nx; if=2*ic-1
    uc(ic,jc,nz)=uf(if,jf,nf)
  end do; end do
END IF

IF(mode.eq.1) THEN

IF(IST.eq.0) THEN
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(1 ,jc,kc)=0.d0
  end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
  nf=2*nx-1
  do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
    uc(nx,jc,kc)=0.d0
  end do; end do
END IF

END IF

CALL BCsgr_MPI(uc,nx,ny,nz,1,1,1,1,1,1)

END SUBROUTINE rstrctMPI


SUBROUTINE rstrct(uc,uf,nx,ny,nz,mode)
   integer nx,ny,nz,kc,kf,jc,jf,ic,if,mode
double precision uc(nx,ny,nz),uf(2*nx-1,2*ny-1,2*nz-1)
double precision, parameter :: w = 1.d0/12.d0

do kc=2,nz-1; kf=2*kc-1
  do jc=2,ny-1; jf=2*jc-1
    do ic=2,nx-1; if=2*ic-1
      uc(ic,jc,kc)=0.5d0*uf(if,jf,kf)+ &
      w*(uf(if+1,jf,kf)+uf(if-1,jf,kf)+uf(if,jf+1,kf)+uf(if,jf-1,kf)+uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do

nf=2*nx-1
do kc=1,nz; kf=2*kc-1; do jc=1,ny; jf=2*jc-1
  uc(1 ,jc,kc)=uf(1 ,jf,kf) 
  uc(nx,jc,kc)=uf(nf,jf,kf)
end do; end do
nf=2*nz-1
do jc=1,ny; jf=2*jc-1; do ic=1,nx; if=2*ic-1
  uc(ic,jc,1 )=uf(if,jf,1 )
  uc(ic,jc,nz)=uf(if,jf,nf)
end do; end do
nf=2*ny-1
do kc=1,nz; kf=2*kc-1; do ic=1,nx; if=2*ic-1
  uc(ic,1 ,kc)=uf(if,1 ,kf)
  uc(ic,ny,kc)=uf(if,nf,kf)
end do; end do

IF(mode.eq.1) THEN

do kc=1,nz; do jc=1,ny
  uc(1 ,jc,kc) = 0.d0
  uc(nx,jc,kc) = 0.d0
end do; end do

END IF

END SUBROUTINE rstrct


SUBROUTINE interpMPI(uf,uc,nx,ny,nz,np,mode)
USE mpivar
USE slfgrv
 integer nx,ny,nz,kc,kf,jc,jf,ic,if,mode,np,nn,kk,i,j,k
double precision uf(0:nx,0:ny,0:nz),uc(0:nx/2+1,0:ny/2+1,0:nz/2+1)

do kc=1,nz/2+1; kf=2*kc-1
  do jc=1,ny/2+1; jf=2*jc-1
    do ic=1,nx/2+1
      uf(2*ic-1,jf,kf)=uc(ic,jc,kc)
    end do
  end do
end do

do kf=1,nz,2
  do jf=1,ny,2
    do if=2,nx,2
      uf(if,jf,kf)=.5d0*(uf(if+1,jf,kf)+uf(if-1,jf,kf))
    end do
  end do
end do
do kf=1,nz,2
  do jf=2,ny,2
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf+1,kf)+uf(if,jf-1,kf))
    end do
  end do
end do
do kf=2,nz,2
  do jf=1,ny
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do

IF(mode.eq.1) THEN

IF(IST.eq.0) THEN
do k=0,nz; kk=(ny+1)*k
do j=0,ny; nn=j+kk
  uf(1 ,j,k)=bphi2(np+nn,1)
end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
do k=0,nz; kk=(ny+1)*k
do j=0,ny; nn=j+kk
  uf(nx,j,k)=bphi2(np+nn,2)
end do; end do
END IF

END IF

CALL BCsgr_MPI(uf,nx,ny,nz,1,0,1,0,1,0)

END SUBROUTINE interpMPI


SUBROUTINE interp(uf,uc,nx,ny,nz,np,mode)
  USE slfgrv
   integer nx,ny,nz,kc,kf,jc,jf,ic,if,mode,np,nn,kk
double precision uf(nx,ny,nz),uc(nx/2+1,ny/2+1,nz/2+1)

do kc=1,nz/2+1; kf=2*kc-1
  do jc=1,ny/2+1; jf=2*jc-1
    do ic=1,nx/2+1
      uf(2*ic-1,jf,kf)=uc(ic,jc,kc)
    end do
  end do
end do

do kf=1,nz,2
  do jf=1,ny,2
    do if=2,nx,2
      uf(if,jf,kf)=.5d0*(uf(if+1,jf,kf)+uf(if-1,jf,kf))
    end do
  end do
end do
do kf=1,nz,2
  do jf=2,ny,2
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf+1,kf)+uf(if,jf-1,kf))
    end do
  end do
end do
do kf=2,nz,2
  do jf=1,ny
    do if=1,nx
      uf(if,jf,kf)=.5d0*(uf(if,jf,kf+1)+uf(if,jf,kf-1))
    end do
  end do
end do

IF(mode.eq.1) THEN !x面
  do kf=1,nz; kk=ny*(kf-1)
  do jf=1,ny; nn=jf-1+kk
    uf(1 ,jf,kf) = bphi1(np+nn,1)
    uf(nx,jf,kf) = bphi1(np+nn,2)
  end do; end do
END IF

END SUBROUTINE interp


SUBROUTINE addintMPI(uf,uc,res,nx,ny,nz)
  use mpivar
   integer nx,ny,nz,i,j,k,isw
   double precision res(0:nx,0:ny,0:nz),uc(0:nx/2+1,0:ny/2+1,0:nz/2+1),uf(0:nx,0:ny,0:nz)


!isw=2 !for speed up
do k=0,nz
  do j=0,ny
    do i=1,nx
      !if(NRANK==40) then
      !    write(*,*) uf(i,j,k),res(i,j,k),k,j,i,'addint1'
      ! end if
    end do
 !   isw=3-isw
 end do
 !isw=3-isw
end do


call interpMPI(res,uc,nx,ny,nz,0,0)
!do k=0,nz; do j=0,ny; do i=0,nx
!  uf(i,j,k)=uf(i,j,k)+res(i,j,k)
!end do; end do; end do

isw=2 !for speed up
do k=0,nz
  do j=0,ny
    do i=isw-1,nx,2
       uf(i,j,k)=uf(i,j,k)+res(i,j,k)
      ! if(NRANK==40) then
      !    write(*,*) uf(i,j,k),res(i,j,k),k,j,i,'addint2'
      ! end if
    end do
    isw=3-isw
  end do
  isw=3-isw
end do

END SUBROUTINE addintMPI


SUBROUTINE addint(uf,uc,res,nx,ny,nz)
  use mpivar
   integer nx,ny,nz,i,j,k,isw
double precision res(nx,ny,nz),uc(nx/2+1,ny/2+1,nz/2+1),uf(nx,ny,nz)
call interp(res,uc,nx,ny,nz,0,0)
!do k=1,nz; do j=1,ny; do i=1,nx
!uf(i,j,k)=uf(i,j,k)+res(i,j,k)
!end do; end do; end do

!do k=0,nz
!  do j=0,ny
!    do i=1,nx
      !if(NRANK==40) then
      !    write(*,*) uf(i,j,k),res(i,j,k),k,j,i,'addint2.5'
      ! end if
!    end do
 !   isw=3-isw
! end do
 !isw=3-isw
!end do

isw=1 ! for speed up
do k=1,nz
  do j=1,ny
    do i=isw,nx,2
       uf(i,j,k)=uf(i,j,k)+res(i,j,k)
       !if(NRANK==40) then
       !   write(*,*) uf(i,j,k),res(i,j,k),k,j,i,'addint3'
       !end if
    end do
    isw=3-isw
  end do
end do

END SUBROUTINE addint


SUBROUTINE slvsml(u,rhs)
  USE slfgrv
  use mpivar
double precision  u(3,3,3),rhs(3,3,3)
double precision h
integer :: i,j,k,count=0
h=0.5d0*Lbox

u(1,1,1)=0.d0;u(3,1,1)=0.d0
u(1,2,1)=0.d0;u(3,2,1)=0.d0
u(1,3,1)=0.d0;u(3,3,1)=0.d0

u(1,1,2)=0.d0;u(3,1,2)=0.d0
u(1,2,2)=0.d0;u(3,2,2)=0.d0
u(1,3,2)=0.d0;u(3,3,2)=0.d0

u(1,1,3)=0.d0;u(3,1,3)=0.d0
u(1,2,3)=0.d0;u(3,2,3)=0.d0
u(1,3,3)=0.d0;u(3,3,3)=0.d0

u(2,1,1) = 0.025d0*h*h*(-8.d0*rhs(2,1,1)-2.d0*rhs(2,1,2)-2.d0*rhs(2,1,3)-2.d0*rhs(2,2,1)-rhs(2,2,2)-rhs(2,2,3)-2.d0*rhs(2,3,1)-rhs(2,3,2)-rhs(2,3,3))
u(2,2,1) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-rhs(2,1,2)-rhs(2,1,3)-8.d0*rhs(2,2,1)-2.d0*rhs(2,2,2)-2.d0*rhs(2,2,3)-2.d0*rhs(2,3,1)-rhs(2,3,2)-rhs(2,3,3))
u(2,3,1) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-rhs(2,1,2)-rhs(2,1,3)-2.d0*rhs(2,2,1)-rhs(2,2,2)-rhs(2,2,3)-8.d0*rhs(2,3,1)-2.d0*rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,1,2) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-8.d0*rhs(2,1,2)-2.d0*rhs(2,1,3)-rhs(2,2,1)-2.d0*rhs(2,2,2)-rhs(2,2,3)-rhs(2,3,1)-2.d0*rhs(2,3,2)-rhs(2,3,3))
u(2,2,2) = 0.025d0*h*h*(-rhs(2,1,1)-2.d0*rhs(2,1,2)-rhs(2,1,3)-2.d0*rhs(2,2,1)-8.d0*rhs(2,2,2)-2.d0*rhs(2,2,3)-rhs(2,3,1)-2.d0*rhs(2,3,2)-rhs(2,3,3))
u(2,3,2) = 0.025d0*h*h*(-rhs(2,1,1)-2.d0*rhs(2,1,2)-rhs(2,1,3)-rhs(2,2,1)-2.d0*rhs(2,2,2)-rhs(2,2,3)-2.d0*rhs(2,3,1)-8.d0*rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,1,3) = 0.025d0*h*h*(-2.d0*rhs(2,1,1)-2.d0*rhs(2,1,2)-8.d0*rhs(2,1,3)-rhs(2,2,1)-rhs(2,2,2)-2.d0*rhs(2,2,3)-rhs(2,3,1)-rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,2,3) = 0.025d0*h*h*(-rhs(2,1,1)-rhs(2,1,2)-2.d0*rhs(2,1,3)-2.d0*rhs(2,2,1)-2.d0*rhs(2,2,2)-8.d0*rhs(2,2,3)-rhs(2,3,1)-rhs(2,3,2)-2.d0*rhs(2,3,3))
u(2,3,3) = 0.025d0*h*h*(-rhs(2,1,1)-rhs(2,1,2)-2.d0*rhs(2,1,3)-rhs(2,2,1)-rhs(2,2,2)-2.d0*rhs(2,2,3)-2.d0*rhs(2,3,1)-2.d0*rhs(2,3,2)-8.d0*rhs(2,3,3))

!if(NRANK==40) then
!do i=1,3
!   do j=1,3
!      do k=1,3
!         write(*,*) u(i,j,k), 'slvsml'
!      end do
!   end do
!end do
!count=2
!end if

END SUBROUTINE slvsml


SUBROUTINE slvsmlb(u,rhs)
USE slfgrv
double precision  u(3,3,3),rhs(3,3,3)
double precision h
integer i,j,k
!integer :: aa,bb,cc,dd=0
!character(14) ff
h=0.5d0*Lbox

u(1,1,1)=bphi1(1,1);  u(3,1,1)=bphi1(1,2)
u(1,2,1)=bphi1(2,1);  u(3,2,1)=bphi1(2,2)
u(1,3,1)=bphi1(3,1);  u(3,3,1)=bphi1(3,2)

u(1,1,2)=bphi1(4,1);  u(3,1,2)=bphi1(4,2)
u(1,2,2)=bphi1(5,1);  u(3,2,2)=bphi1(5,2)
u(1,3,2)=bphi1(6,1);  u(3,3,2)=bphi1(6,2)

u(1,1,3)=bphi1(7,1);  u(3,1,3)=bphi1(7,2)
u(1,2,3)=bphi1(8,1);  u(3,2,3)=bphi1(8,2)
u(1,3,3)=bphi1(9,1);  u(3,3,3)=bphi1(9,2)

u(2,1,1) = 0.025d0* (8.d0*u(1,1,1) + 2.d0*u(1,1,2) + 2.d0*u(1,1,3) + 2.d0*u(1,2,1) + u(1,2,2) + u(1,2,3) + 2.d0*u(1,3,1) + &
u(1,3,2) + u(1,3,3) + 8.d0*u(3,1,1) + 2.d0*u(3,1,2) + 2.d0*u(3,1,3) + 2.d0*u(3,2,1) + u(3,2,2) + u(3,2,3) + 2.d0*u(3,3,1) + u(3,3,2) + u(3,3,3) &
 - 8.d0*h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,2,1) = 0.025d0* (2.d0*u(1,1,1) + u(1,1,2) + u(1,1,3) + 8.d0*u(1,2,1) + 2.d0*u(1,2,2) + 2.d0*u(1,2,3) + 2.d0*u(1,3,1) + &
u(1,3,2) + u(1,3,3) + 2.d0*u(3,1,1) + u(3,1,2) + u(3,1,3) + 8.d0*u(3,2,1) + 2.d0*u(3,2,2) + 2.d0*u(3,2,3) + 2.d0*u(3,3,1) + u(3,3,2) + u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - 8.d0*h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,3,1) = 0.025d0* (2.d0*u(1,1,1) + u(1,1,2) + u(1,1,3) + 2.d0*u(1,2,1) + u(1,2,2) + u(1,2,3) + 8.d0*u(1,3,1) + 2.d0*u(1,3,2) + &
2.d0*u(1,3,3) + 2.d0*u(3,1,1) + u(3,1,2) + u(3,1,3) + 2.d0*u(3,2,1) + u(3,2,2) + u(3,2,3) + 8.d0*u(3,3,1) + 2.d0*u(3,3,2) + 2.d0*u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
8.d0*h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,1,2) = 0.025d0* (2.d0*u(1,1,1) + 8.d0*u(1,1,2) + 2.d0*u(1,1,3) + u(1,2,1) + 2.d0*u(1,2,2) + u(1,2,3) + u(1,3,1) + 2.d0*u(1,3,2) + &
u(1,3,3) + 2.d0*u(3,1,1) + 8.d0*u(3,1,2) + 2.d0*u(3,1,3) + u(3,2,1) + 2.d0*u(3,2,2) + u(3,2,3) + u(3,3,1) + 2.d0*u(3,3,2) + u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - 8.d0*h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,2,2) = 0.025d0* ( u(1,1,1) + 2.d0*u(1,1,2) + u(1,1,3) + 2.d0*u(1,2,1) + 8.d0*u(1,2,2) + 2.d0*u(1,2,3) + u(1,3,1) + 2.d0*u(1,3,2) + &
u(1,3,3) + u(3,1,1) + 2.d0*u(3,1,2) + u(3,1,3) + 2.d0*u(3,2,1) + 8.d0*u(3,2,2) + 2.d0*u(3,2,3) + u(3,3,1) + 2.d0*u(3,3,2) + u(3,3,3) &
- h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - 8.d0*h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - h*h*rhs(2,3,3))
u(2,3,2) = 0.025d0* ( u(1,1,1) + 2.d0*u(1,1,2) + u(1,1,3) + u(1,2,1) + 2.d0*u(1,2,2) + u(1,2,3) + 2.d0*u(1,3,1) + 8.d0*u(1,3,2) + &
2.d0*u(1,3,3) + u(3,1,1) + 2.d0*u(3,1,2) + u(3,1,3) + u(3,2,1) + 2.d0*u(3,2,2) + u(3,2,3) + 2.d0*u(3,3,1) + 8.d0*u(3,3,2) + 2.d0*u(3,3,3) &
- h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - 8.d0*h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,1,3) = 0.025d0* (2.d0*u(1,1,1) + 2.d0*u(1,1,2) + 8.d0*u(1,1,3) + u(1,2,1) + u(1,2,2) + 2.d0*u(1,2,3) + u(1,3,1) + u(1,3,2) + &
2.d0*u(1,3,3) + 2.d0*u(3,1,1) + 2.d0*u(3,1,2) + 8.d0*u(3,1,3) + u(3,2,1) + u(3,2,2) + 2.d0*u(3,2,3) + u(3,3,1) + u(3,3,2) + 2.d0*u(3,3,3) &
- 2.d0*h*h*rhs(2,1,1) - 2.d0*h*h*rhs(2,1,2) - 8.d0*h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,2,3) = 0.025d0* ( u(1,1,1) + u(1,1,2) + 2.d0*u(1,1,3) + 2.d0*u(1,2,1) + 2.d0*u(1,2,2) + 8.d0*u(1,2,3) + u(1,3,1) + u(1,3,2) + &
2.d0*u(1,3,3) + u(3,1,1) + u(3,1,2) + 2.d0*u(3,1,3) + 2.d0*u(3,2,1) + 2.d0*u(3,2,2) + 8.d0*u(3,2,3) + u(3,3,1) + u(3,3,2) + 2.d0*u(3,3,3) &
- h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - 2.d0*h*h*rhs(2,2,1) - 2.d0*h*h*rhs(2,2,2) - 8.d0*h*h*rhs(2,2,3) - &
h*h*rhs(2,3,1) - h*h*rhs(2,3,2) - 2.d0*h*h*rhs(2,3,3))
u(2,3,3) = 0.025d0* ( u(1,1,1) + u(1,1,2) + 2.d0*u(1,1,3) + u(1,2,1) + u(1,2,2) + 2.d0*u(1,2,3) + 2.d0*u(1,3,1) + 2.d0*u(1,3,2) + &
8.d0*u(1,3,3) + u(3,1,1) + u(3,1,2) + 2.d0*u(3,1,3) + u(3,2,1) + u(3,2,2) + 2.d0*u(3,2,3) + 2.d0*u(3,3,1) + 2.d0*u(3,3,2) + 8.d0*u(3,3,3) &
- h*h*rhs(2,1,1) - h*h*rhs(2,1,2) - 2.d0*h*h*rhs(2,1,3) - h*h*rhs(2,2,1) - h*h*rhs(2,2,2) - 2.d0*h*h*rhs(2,2,3) - &
2.d0*h*h*rhs(2,3,1) - 2.d0*h*h*rhs(2,3,2) - 8.d0*h*h*rhs(2,3,3))


!do k=1,3
!   do j=1,3
!      do i=1,3
!         write(*,*) u(i,j,k) , '========slvb========'
!      end do
!   end do
!end do
!write (ff,'("slvsmlb",i3.3,".dat")') dd
!open(505,file=ff)
!do aa=1,3
!   do bb=1,3
!      do cc=1,3
!         write(505,*) u(aa,bb,cc)
!      end do
!   end do
!end do
!close(505)
!dd=dd+1
END SUBROUTINE slvsmlb


SUBROUTINE relaxMPI(u,rhs,nx,ny,nz,mode)
USE mpivar
USE slfgrv
integer :: i,j,k,mode,li,lj,lk,nx,ny,nz,isw,jsw,ifl,shu,count=0
double precision, parameter :: w=1.d0/6.d0
double precision u(0:nx,0:ny,0:nz),rhs(0:nx,0:ny,0:nz)
double precision h,h2


if(NRANK==40) then
write(*,*) mode,nx,ny,nz
do i=0,nx
   do j=0,ny
      do k=0,nz
         write(*,*) u(i,j,k),rhs(i,j,k),count ,'rMPI=========rMPI'
      end do
   end do
end do
count = 1+count
end if
!integer :: bugnum=0
!character(14) bugn

h=Lbox/dble(((nx-1)*NSPLTx))
h2=h*h


li=0; IF(IST.eq.0) li=1
lj=0; IF(JST.eq.0) lj=1
lk=0; IF(KST.eq.0) lk=1
do jsw=2,1,-1 !Red-Black Gauss-Seidel
  isw=jsw
  if(mode.eq.1) then
     ! write (bugn,'("rMPI",i3.3,".dat")') bugnum
     ! open(506,file=bugn)

    do k=1,nz-1
      do j=1,ny-1
         do i=isw,nx-1,2
            ifl=li*int(1/i)
           ! IF(NRANK==40) then
       ! write(*,*) u(i,j,k),rhs(i,j,k),w,h2,ifl,0.5d0+dsign(0.5d0,dble(ifl)-0.5d0),i,j,k,NRANK,isw,nx,ny,nz,'point2-mpirelax-2'
           ! end IF
          u(i,j,k)=-w*h2*rhs(i,j,k) &
          *dble((1-ifl)) + (0.5d0+dsign(0.5d0,dble(ifl)-0.5d0))*u(i,j,k)
!		  if(NRANK.eq.0) write(*,*) 'i=',i,'ifl=',ifl
!		  if(NRANK.eq.0) write(*,*) (1-ifl),(0.5d0+dsign(0.5d0,ifl-0.5d0))
          !write(506,*) u(i,j,k),rhs(i,j,k),w,h2,NRANK,'point1'

          !IF(NRANK==40) then
        !   write(*,*) u(i,j,k),rhs(i,j,k),w,h2,ifl,0.5d0+dsign(0.5d0,dble(ifl)-0.5d0),i,j,k,NRANK,nx,isw,'point2-mpirelax-1'
         ! end IF
       end do
        isw=3-isw
      end do
      isw=3-isw
   end do
   !close(506)
   !bugnum=bugnum+1
   mode=2
else
   !write (bugn,'("rMPI",i3.3,".dat")') bugnum
   ! open(507,file=bugn)
   do k=1,nz-1
      do j=1,ny-1
         do i=isw,nx-1,2
            ifl=li*int(1/i)
             IF(NRANK==40) then
                !      write(*,*) u(i,j,k),rhs(i,j,k),w,h2,ifl,0.5d0+dsign(0.5d0,dble(ifl)-0.5d0),isw,i,j,k,nx,NRANK,'point2-mpirelax1'
                !      write(*,*) u(i+1,j,k),u(i-1,j,k),u(i,j+1,k),u(i,j-1,k),u(i,j,k+1),u(i,j,k-1),h2*rhs(i,j,k),isw,i,j,k,nx,NRANK,'point2-mpirelax11'
            end IF
          u(i,j,k)=w*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-h2*rhs(i,j,k)) &
               *dble((1-ifl)) + (0.5d0+dsign(0.5d0,dble(ifl)-0.5d0))*u(i,j,k) !==============double=======================
          !  write(507,*) u(i,j,k),rhs(i,j,k),w,h2,NRANK,'point2-mpirelax'
          !  IF(NRANK==40) then
          !      write(*,*) u(i,j,k),rhs(i,j,k),w,h2,ifl,0.5d0+dsign(0.5d0,dble(ifl)-0.5d0),i,j,k,NRANK,'point2-mpirelax2'
          !  end IF

          if(NRANK==40.and.i==nx/2.and.j==ny/2.and.k==nz/2) then
              write(*,*) u(i,j,k),shu, 'shusoku-MPI????'
              shu=shu+1
          end if
        end do
        isw=3-isw
      end do
      isw=3-isw
   end do
   !close(507)
   !bugnum=bugnum+1
end if
if(NRANK==40) then
write(*,*) mode,nx,ny,nz
do i=0,nx
   do j=0,ny
      do k=0,nz
         write(*,*) u(i,j,k),rhs(i,j,k),count ,'rMPI=========rMPI----2'
      end do
   end do
end do
count = 1+count
end if
  CALL BCsgr_MPI(u,nx,ny,nz,1,1,1,1,1,1)
end do

END SUBROUTINE relaxMPI


SUBROUTINE relax(u,rhs,nx,ny,nz,mode)
  USE slfgrv
  use mpivar
  integer :: i,j,k,mode,nx,ny,nz,isw,jsw,km,kp,jm,jp,shu=0
double precision, parameter :: w=1.d0/6.d0
double precision u(nx,ny,nz),rhs(nx,ny,nz)
double precision h,h2
h=Lbox/dble((nx-1))
h2=h*h

do jsw=1,2 !Red-Black Gauss-Seidel
  isw=jsw
  if(mode.eq.1) then
    do k=1,nz
      do j=1,ny
        do i=isw+1,nx-1,2
          u(i,j,k)=-w*h2*rhs(i,j,k)
        end do
        isw=3-isw
      end do
    end do
    mode=2
  else
    do k=1,nz
      km = k-1; if(k.eq.1) km = nz
!      km = (k-1)*(1+isign(1,k-2)   )/2 + nz*(1-isign(1,k-2)   )/2 !original
      kp = k+1; if(k.eq.nz) kp = 1
!      kp = (k+1)*(1+isign(1,nz-1-k))/2 +  1*(1-isign(1,nz-1-k))/2 !original
      do j=1,ny
        jm = j-1; if(j.eq.1) jm = ny
!        jm = (j-1)*(1+isign(1,j-2)   )/2 + ny*(1-isign(1,j-2)   )/2 !original
        jp = j+1; if(j.eq.ny) jp = 1
!        jp = (j+1)*(1+isign(1,ny-1-j))/2 +  1*(1-isign(1,ny-1-j))/2 !original
        do i=isw+1,nx-1,2
           u(i,j,k)=w*(u(i+1,j,k)+u(i-1,j,k)+u(i,jp,k)+u(i,jm,k)+u(i,j,kp)+u(i,j,km)-h2*rhs(i,j,k))
           !if(NRANK==40.and.i==nx/2.and.j==ny/2.and.k==nz/2) then
           !   write(*,*) u(i,j,k),shu, 'shusoku????'
           !   shu=shu+1
           !end if
       end do
       isw=3-isw
    end do
    end do
  end if
end do

END SUBROUTINE relax


SUBROUTINE residMPI(res,u,rhs,nx,ny,nz)
USE mpivar
USE slfgrv
integer i,j,k,nx,ny,nz
double precision res(0:nx,0:ny,0:nz),rhs(0:nx,0:ny,0:nz),u(0:nx,0:ny,0:nz)
double precision h,h2i
h=Lbox/dble(((nx-1)*NSPLTx))
h2i=1.d0/(h*h)

!is=1; IF(IST.eq.0) is=2
!js=1; IF(JST.eq.0) js=2
!ks=1; IF(KST.eq.0) ks=2
!do k=ks,nz-1; do j=js,ny-1; do i=is,nx-1
!  res(i,j,k)=-h2i*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-6.d0*u(i,j,k))+rhs(i,j,k) 
!end do; end do; end do

!write(*,*) 'residMPI'

do k=1,nz-1; do j=1,ny-1; do i=1,nx-1
 !  if(NRANK==40) then
 !     write(*,*) NRANK , i,j,k , h,h2i , u(i,j,k),rhs(i,j,k),res(i,j,k),'rMPI--1'
 !  end if
   res(i,j,k)=-h2i*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-6.d0*u(i,j,k))+rhs(i,j,k)
 !  if(NRANK==40) then
 !     write(*,*) NRANK , i,j,k , h,h2i , u(i,j,k),rhs(i,j,k),res(i,j,k),'rMPI--2'
 !  end if
end do; end do; end do

IF(IST.eq.0) THEN
do k=1,nz,2; do j=1,ny,2 !for speed up
!do k=1,nz; do j=1,ny
  res(1 ,j,k)=0.d0
end do; end do
END IF
IF(IST.eq.NSPLTx-1) THEN
do k=1,nz,2; do j=1,ny,2 !for speed up
!do k=1,nz; do j=1,ny
  res(nx,j,k)=0.d0
end do; end do
END IF

CALL BCsgr_MPI(res,nx,ny,nz,1,1,1,1,1,1)

END SUBROUTINE residMPI


SUBROUTINE resid(res,u,rhs,nx,ny,nz)
  USE slfgrv
  integer i,j,k,nx,ny,nz
double precision res(nx,ny,nz),rhs(nx,ny,nz),u(nx,ny,nz)
double precision h,h2i
h=Lbox/dble((nx-1))
h2i=1.d0/(h*h)
do k=1,nz; do j=1,ny; do i=2,nx-1
  res(i,j,k)=-h2i*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)+u(i,j,k+1)+u(i,j,k-1)-6.d0*u(i,j,k))+rhs(i,j,k)
end do; end do; end do

do k=1,nz,2; do j=1,ny,2 !for speed up
!do k=1,nz; do j=1,ny
  res(1 ,j,k)=0.d0
  res(nx,j,k)=0.d0
end do; end do

END SUBROUTINE resid


SUBROUTINE copy(aout,ain,nx,ny,nz)
   integer i,j,k,nx,ny,nz
double precision ain(nx,ny,nz),aout(nx,ny,nz)
do k=1,nz; do j=1,ny; do i=1,nx; aout(i,j,k)=ain(i,j,k); end do; end do; end do
END SUBROUTINE copy

SUBROUTINE copyMPI(aout,ain,nx,ny,nz)
   integer i,j,k,nx,ny,nz
double precision ain(0:nx,0:ny,0:nz),aout(0:nx,0:ny,0:nz)
do k=0,nz; do j=0,ny; do i=0,nx; aout(i,j,k)=ain(i,j,k); end do; end do; end do
END SUBROUTINE copyMPI

SUBROUTINE fill0(u,nx,ny,nz)
   integer i,j,k,nx,ny,nz
double precision u(nx,ny,nz)
do j=1,ny; do i=1,nx
  u(i,j,1 )=0.d0
  u(i,j,nz)=0.d0
end do; end do
do k=1,nz; do i=1,nx
  u(i,1 ,k)=0.d0
  u(i,ny,k)=0.d0
end do; end do
do k=1,nz; do j=1,ny
  u(1 ,j,k)=0.d0
  u(nx,j,k)=0.d0
end do; end do
END SUBROUTINE fill0

SUBROUTINE fill0MPI(u,nx,ny,nz)
   integer i,j,k,nx,ny,nz
double precision u(0:nx,0:ny,0:nz)
do j=0,ny; do i=0,nx
  u(i,j,0 )=0.d0
  u(i,j,1 )=0.d0
  u(i,j,nz)=0.d0
end do; end do
do k=0,nz; do i=0,nx
  u(i,0 ,k)=0.d0
  u(i,1 ,k)=0.d0
  u(i,ny,k)=0.d0
end do; end do
do k=0,nz; do j=0,ny
  u(0 ,j,k)=0.d0
  u(1 ,j,k)=0.d0
  u(nx,j,k)=0.d0
end do; end do
END SUBROUTINE fill0MPI


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

integer nl,nc,nx,ny,nz,i,j,k,Nlp,isend,KSs,JSs,kz,jy,irecv,KSr,JSr
integer Nis,kls,Nir,klr,nn1,nn2,m,m2,ncx,ncy,ncz,lc,nfx,nfy,nfz,mcy,mcx
integer kf,kc,jf,jc,if,ic,kk,jb,kbb

double precision w

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
    kap = sqrt(kap)
    dat1(2*l-1,m) = dat1(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat1(2*l  ,m) = dat1(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat2(2*l-1,m) = dat2(2*l-1,m) + data(2*l-1,m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    dat2(2*l  ,m) = dat2(2*l  ,m) + data(2*l  ,m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do;end do

  l=nn1/2+1
  do m=1,nn2/2+1
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)
    spe1(m) = spe1(m) + speq(m,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    spe2(m) = spe2(m) + speq(m,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do

  do m2=nn2/2+2,nn2; m=nn2+2-m2; do l=1,nn1/2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)
    dat1(2*l-1,m2) = dat1(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat1(2*l  ,m2) = dat1(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*exp(-zp1*kap)/kap *facG
    dat2(2*l-1,m2) = dat2(2*l-1,m2) + data(2*l-1,m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
    dat2(2*l  ,m2) = dat2(2*l  ,m2) + data(2*l  ,m2,klr)* 0.5d0*exp(-zp2*kap)/kap *facG
  end do;end do

  l=nn1/2+1
  do m2=nn2/2+2,nn2; m=nn2+2-m2
    kap = 4.d0*( sin(pi*(l-1)/nn1)**2/dxx**2 + sin(pi*(m-1)/nn2)**2/dyy**2 )
    kap = sqrt(kap)
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
end do; end do


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
  do l = 1, 2
  if(l.eq.1) then; mfx=nfy; mfy=nfz; mcx=ncy; mcy=ncz; end if
  if(l.eq.2) then; mfx=nfy; mfy=nfz; mcx=ncy; mcy=ncz; end if

  do jc = 2,mcy-1; jf = 2*jc-1
  do ic = 2,mcx-1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 4.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1    ,l)+bphi2(pointb2(lc+1)+kf+1    ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
  end do
  end do

  do jc = 2,mcy-1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)+bphi2(pointb2(lc+1)+kf+mfx+1,l) )
  end do
  do ic = 2,mcx-1; if = 2*ic-1
    jc = 1; jf = 2*jc-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    jc = mcy; jf = 2*jc-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 5.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )
  end do

    jc = 1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    jc = 1 ; jf = 2*jc-1
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1   ,l)+ &
                                                                                                    bphi2(pointb2(lc+1)+kf+mfx+1,l) )
    jc = mcy; jf = 2*jc-1
    ic = 1 ; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*(                                bphi2(pointb2(lc+1)+kf+1   ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )
    jc = mcy; jf = 2*jc-1
    ic = mcx; if = 2*ic-1
    kf = if+(mcx*2)*jf
    kc = ic+(mcx+1)*jc
    bphi2(pointb2(lc)+kc,l) = 6.d0*w*bphi2(pointb2(lc+1)+kf,l) + w*( bphi2(pointb2(lc+1)+kf-1 ,l)+ &
                                                                     bphi2(pointb2(lc+1)+kf-mfx-1,l)                                )

  end do
end do

nz=(2**NGcr)/NSPLTz+1; ny=(2**NGcr)/NSPLTy+1; nx=(2**NGcr)/NSPLTx+1; n1=2**NGcr+1
ALLOCATE( temp1(n1,n1,n1), temp2(0:nx,0:ny,0:nz) )


do k=0,nz; kk = (ny+1)*k + pointb2(NGcr)
do j=0,ny; n = j + kk
  temp2(1 ,j,k) = bphi2(n,1)
end do;end do
call collect( temp1,temp2,n1,n1,n1,nx,ny,nz )
do k=1,n1; kk = n1*(k-1) + pointb1(NGcr)
do j=1,n1; n = j-1 + kk
  bphi1(n,1) = temp1(1 ,j,k)
end do;end do

do k=0,nz; kk = (ny+1)*k + pointb2(NGcr)
do j=0,ny; n = j + kk
  temp2(nx,j,k) = bphi2(n,2)
end do;end do
call collect( temp1,temp2,n1,n1,n1,nx,ny,nz )
do k=1,n1; kk = n1*(k-1) + pointb1(NGcr)
do j=1,n1; n = j-1 + kk
  bphi1(n,2) = temp1(n1,j,k)
end do;end do

DEALLOCATE(temp1,temp2)


nc = (2**NGcr)+1
w  = 0.125d0
do lc=NGcr-1,1,-1
  nf = nc
  nc = nc/2+1
  do l = 1, 2

  do jc = 2,nc-1; jf = 2*jc-1
  do ic = 2,nc-1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 4.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
  end do
  end do

  do jc = 2,nc-1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)+bphi1(pointb1(lc+1)+kf+nf,l) )
  end do
  do ic = 2,nc-1; if = 2*ic-1
    jc = 1; jf = 2*jc-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
    jc = nc; jf = 2*jc-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 5.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )
  end do

    jc = 1; jf = 2*jc-1
    ic = 1; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
    jc = 1 ; jf = 2*jc-1
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                                                  bphi1(pointb1(lc+1)+kf+nf,l) )
    jc = nc; jf = 2*jc-1
    ic = 1 ; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*(                              bphi1(pointb1(lc+1)+kf+1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )
    jc = nc; jf = 2*jc-1
    ic = nc; if = 2*ic-1
    kf = if-1+(nc*2-1)*(jf-1)
    kc = ic-1+(nc)    *(jc-1)
    bphi1(pointb1(lc)+kc,l) = 6.d0*w*bphi1(pointb1(lc+1)+kf,l) + w*( bphi1(pointb1(lc+1)+kf-1 ,l)+ &
                                                                     bphi1(pointb1(lc+1)+kf-nf,l)                            )

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
END SUBROUTINE
