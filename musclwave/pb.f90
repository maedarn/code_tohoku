  mudule aavar
  double precision , dimension(:,:,:) , allocatable :: phibc
  end module

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

!MIYAMA method ---------------------------------------------------

!ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,Ncellx+1),speq(Ncellz*NSPLTz,Ncellx))
ALLOCATE(data(Ncelly*NSPLTy,Ncellz*NSPLTz,-1:Ncellx+2),speq(Ncellz*NSPLTz,-1:Ncellx+2))
ALLOCATE(dat1(Ncelly*NSPLTy,Ncellz*NSPLTz),spe1(Ncellz*NSPLTz), &
         dat2(Ncelly*NSPLTy,Ncellz*NSPLTz),spe2(Ncellz*NSPLTz))

!nccy = Ncelly/NSPLTy; nccz = Ncellz/NSPLTz
nccy = Ncelly; nccz = Ncellz
do k=1,Ncellz; kz=KST*Ncellz+k
do j=1,Ncelly; jy=JST*Ncelly+j
do i=-1,Ncellx+2
  data(jy,kz,i) = U(i,j,k,1)
end do;end do;end do

!***************fordebug*****************
!CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
!***************fordebug*****************

                    !count, blocklength, stride
CALL MPI_TYPE_VECTOR(Ncellz,Ncelly,Ncelly*NSPLTy,MPI_REAL8,VECU,IERR)
CALL MPI_TYPE_COMMIT(VECU,IERR)

do Nlp = 1,NSPLTy*NSPLTz-1

  isend = NRANK + NSPLTx*Nlp; if(isend.ge.NPE) isend = isend - NPE !x=const面に送る  : ifは0に行く時
  KSs = isend/(NSPLTx*NSPLTy); JSs = isend/NSPLTx-NSPLTy*KSs !その面の y,z 位置
  irecv = NRANK - NSPLTx*Nlp; if(irecv.lt.0  ) irecv = irecv + NPE !x=const面 逆から
  KSr = irecv/(NSPLTx*NSPLTy); JSr = irecv/NSPLTx-NSPLTy*KSr !その面の y,z 位置

  Nis = JSs + NSPLTy*KSs !JST = NRANK/NSPLTx-NSPLTy*KST x位置の逆とき
  kls = Nis + 1
  !Nir = JSr + NSPLTy*KSr
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



!-----for debug------
write(*,*) '--PB3--'
!-----for debug------





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
  !IF(IST==0) then
  !   Phi(1,j,k)= bphi2(pointb2(NGL)+n,1)
  !end IF
!  write(30,*) bphi2(pointb2(NGL)+n,1)
!  write(560,*) bphi2(pointb2(NGL)+n,2)
end do
!write(30,*)
!write(560,*)
end do

!IF(IST==0) then

!   Phi(i,j,k) =






IF(IST==0) then
do k=0,ncz; kk= (ncy+1)*k
do j=0,ncy; n = j+kk
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 1
  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 1
  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz

  Phi(1,j,k)= bphi2(pointb2(NGL)+n,1)

end do
end do
end IF


IF(IST==NSPLTx-1) then
do k=0,ncz; kk= (ncy+1)*k
do j=0,ncy; n = j+kk
  jb  = JST*Ncelly + j
  kbb = KST*Ncellz + k
  if((j.eq.ncy).and.(JST.eq.NSPLTy-1)) jb  = 1
  if((k.eq.ncz).and.(KST.eq.NSPLTz-1)) kbb = 1
  if((j.eq.0  ).and.(JST.eq.0       )) jb  = Ncelly*NSPLTy
  if((k.eq.0  ).and.(KST.eq.0       )) kbb = Ncellz*NSPLTz

  Phi(Ncellx,j,k)= bphi2(pointb2(NGL)+n,2)

end do
end do
end IF


!-----for debug------
write(*,*) '--PB4--'
!-----for debug------


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
!goto 4589

DEALLOCATE(data,speq)
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
