!******************************************************************!
!multiglid method    - numerical recipe inF77 -
!mglin rstrct interp addint slvsml relax resid copy fill0 maloc
!------------------------------------------------------------------!
!arrange some subroutine
!main initia position save comvar BC converge
!******************************************************************!

module comvar
  implicit none
  !integer(8),PARAMETER :: NG=5,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-int(100/3)
  integer,PARAMETER :: NG=5
  integer(8),PARAMETER :: MEMLEN=5000
  integer, PARAMETER :: NPRE=5,NPOST=1 !ガウスサイデル反復,smoosing
  integer mem
  DOUBLE PRECISION z(MEMLEN)
end module comvar

program main
  implicit none
  INTEGER :: n=32
  integer ncycle
  double precision u(1:32,1:32) , Px(1:32),Py(1:32)
  call INITIA(u,n)
  call position(Px,Py,n)
  call mglin(u,n,ncycle)
  call save(u,Px,Py,n)
end program main

subroutine INITIA(u,n)
  integer n
  double precision u(n,n)
  integer i,j
  do i = 1,n
     do j=1,n
        u(i,j) = 0.0d0
     end do
  end do
  !u(16,16) = 100.0d0
end subroutine INITIA

subroutine position(Px,Py,n)
  integer n
  integer i
  double precision hx,hy
  double precision Px(n),Py(n)
  hx=1.0d0/dble(n)
  hy=1.0d0/dble(n)
  do i=1,n
     Px(i) = hx*i
  end do
  do i=1,n
     Py(i) = hy*i
  end do
end subroutine position


subroutine save(u,Px,Py,n)
  integer n
  double precision u(n,n)
  double precision Px(n),Py(n)
  integer i,j
  !CHARACTER(22) :: dir='/Users/ryunosukemaeda/'
  CHARACTER(32) :: dir='/Users/maeda/Desktop/code/multi/'
  character(4) :: filenm='ml01'
100 format(E19.10e3,E19.10e3,E19.10e3)
  open(10,file=dir//filenm//'.dat')
  do i = 1,n
     do j=1,n
        write(10,100) Px(i) , Py(j) , u(i,j)
     end do
     write(10,*)
  end do
end subroutine save


SUBROUTINE mglin(u,n,ncycle)
  use comvar
  INTEGER n,ncycle
  !integer,PARAMETER :: NG=8,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3
  !integer, PARAMETER :: NPRE=1,NPOST=1
  INTEGER j,jcycle,jj,jpost,jpre,nf,ngrid,nn,ires(NG), irho(NG),irhs(NG),iu(NG),maloc
  !DOUBLE PRECISION z
  double precision u(n,n)
  !COMMON /memory/ z(MEMLEN),mem  !Storage for grid functions is allocated by maloc
  mem=0 !from array z.
  nn=n/2+1     !nはx,yのメッシュ数,nnは半分の地点
  ngrid=NG-1 !レベル2**Nみたいなもの（何回イタレーションするか）
  irho(ngrid) = maloc(int(nn**2))   ! Allocate storage for r.h.s. on grid NG − 1,
  call rstrct(z(irho(ngrid)),u,nn)   !and ll it by restricting from the ne grid. 初期のuにはrhs(密度)が収められている
1 if (nn.gt.3) then   !Similarly allocate storage and ll r.h.s. on all
     nn=nn/2+1 !coarse grids.
     ngrid=ngrid-1
     irho(ngrid)=maloc(int(nn**2))
     call rstrct(z(irho(ngrid)),z(irho(ngrid+1)),nn) !粗いグリッドと細かいグリッドを引数に取る
     goto 1
  endif
  nn=3
  iu(1)=maloc(int(nn**2)) !1=ngrid
  irhs(1)=maloc(int(nn**2))
  call slvs  double precision u(n,n) , cnv(0:10) , sum
  integer :: i=1
  integer j,k

  cnv(0)=0.0d0
  cnv(i)=0.0d0
  sum =0.0d0
  do j=1,n
     do k=1,n
        sum = u(i,k) + sum
     end do
  end do
  cnv(i)=sum
  write(*,*) cnv(i)-cnv(i-1)
  i=i+1
end subroutine converge

