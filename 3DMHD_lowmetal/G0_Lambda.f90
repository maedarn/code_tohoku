PROGRAM main
double precision, parameter :: Zsolar=1.d-2,zeta_cr=1.d0
DOUBLE PRECISION, parameter :: mH=1.d0, mHe=4.d0, mH2=2.d0, mC=12.d0, mCO=28.d0, TCMB=3.d-3
DOUBLE PRECISION, parameter :: G0=1.d0, xc=1.4d-4*Zsolar, xo=3.2d-4*Zsolar, dv=2.d0, Tgr=5.d-3, fgr=1.d0*Zsolar
REAL(4), ALLOCATABLE, DIMENSION(:,:,:,:) :: U,NH2,Nfgr
REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: G03D
!CHARACTER(49) :: dir='/glv0/maedarn/low_metal/samplecnv-10m0-2-v100-tp/'
character(39) :: dir='/glv0/maedarn/low_metal/samplecnv-10m4/'
CHARACTER(3) timech
REAL(4) :: tNH2,tNfgr,dxh,Lbox
integer :: val,nx,ny,nz,initime,timeloop,timestep
DOUBLE PRECISION :: Av1,Av2,x1,x2,ATN2,SHLD1,SHLD2

U(:,:,:,:)=0.e0

!----parameter----
nx=512
ny=512
nz=512
val=20 !18
Lbox=100.e0
!Msun=1.473d-2 !1pc * 1pc * 1pc * 1m_p/cc
!Msun=2.4d-2 !1pc * 1pc * 1pc * 1m_p/cc
initime=1
timeloop=96
timestep=1
jump=1
time1=130
dxh=100.e0/512.e0*0.5e0
!----parameter----

allocate(U(1:nx,1:ny,1:nz,val),G03D(1:nx,1:ny,1:nz),NH2(1:nx,1:ny,1:nz,val),Nfgr(1:nx,1:ny,1:nz,val))

write(timech,'(i3.3)') time1
open(110,file=dir//'All/WCNALL'//timech//'.DAT',access='stream',FORM='UNFORMATTED')

do k=1,nz
   do j=1,ny
      do i=1,nx
         read(110) (U(i,j,k,n),n=1,val)
      end do
   end do
end do
close(110)

do k = 1, nz; do j = 1, ny
  tNH2=0.d0; tNfgr=0.d0
  do i = 1, nx
    tNH2  = tNH2  + dxh *  U(i,j,k,11) !11=nH2
    tNfgr   = tNfgr   + dxh *  U(i,j,k,19) !* ndtot(i,j,k) 19=nfgr
    NH2(i,j,k,1)  = tNH2
    Nfgr(i,j,k,1)   = tNfgr
    tNH2  = tNH2  + dxh *  U(i,j,k,11)
    tNfgr   = tNfgr   + dxh *  U(i,j,k,19) !* ndtot(i,j,k)
  end do
  tNH2=0.d0; tNfgr=0.d0
  do i = nx, 1, -1
    tNH2  = tNH2  + dxh *  U(i,j,k,11)
    tNfgr   = tNfgr   + dxh *  U(i,j,k,19) !* ndtot(i,j,k)
    NH2(i,j,k,2)  = tNH2
    Nfgr(i,j,k,2)   = tNfgr
    tNH2  = tNH2  + dxh *  U(i,j,k,11)
    tNfgr   = tNfgr   + dxh *  U(i,j,k,19) !* ndtot(i,j,k)
  end do
  !NtH2(j,k,IST)=tNH2; Ntfgr(j,k,IST)=tNfgr
end do; end do


do k = 1, nz; do j = 1, ny; do i = 1, nx
 Av1  = 1.63542d-3*dble(Nfgr(i,j,k,1)); x1 = 6.1714d3*dble(NH2(i,j,k,1))
 Av2  = 1.63542d-3*dble(Nfgr(i,j,k,2)); x2 = 6.1714d3*dble(NH2(i,j,k,2))
 ATN2 = ( dexp(-3.77358d0*Av1) +dexp(-3.77358d0*Av2) )*0.5d0

 SHLD1 = 0.965d0/(1.d0+x1/dv)**2 + 0.035d0/dsqrt(1.d0+x1)*dexp(-8.5d-4*dsqrt(1.d0+x1))
 SHLD2 = 0.965d0/(1.d0+x2/dv)**2 + 0.035d0/dsqrt(1.d0+x2)*dexp(-8.5d-4*dsqrt(1.d0+x2))

 G03D(i,j,k)=sngl(G0*ATN2 * ( SHLD1+SHLD2 )*0.5d0)
end do; end do; end do


open(120,file=dir//'All/G03D'//timech//'.DAT',access='stream',FORM='UNFORMATTED')

do k=1,nz
   do j=1,ny
      do i=1,nx
         write(120) G03D(i,j,k)
      end do
   end do
end do
close(120)

deallocate(U,G03D,NH2,Nfgr)
END PROGRAM main