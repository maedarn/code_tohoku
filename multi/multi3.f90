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
  integer :: ncycle=1
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
  call slvsml(z(iu(1)),z(irho(1))) !Initial solution on coarsest grid. 初期値uの
  ngrid=NG
  do j=2,ngrid! Nested iteration loop. 粗い位置から始まる (前のv-loopnoの一番下から)
     nn=2*nn-1 !ふやしていく
     iu(j)=maloc(int(nn**2))
     irhs(j)=maloc(int(nn**2))
     ires(j)=maloc(int(nn**2))
     call interp(z(iu(j)),z(iu(j-1)),nn) !Interpolate from coarse grid to next ner grid. zに格納（挿入）細かく

     if (j.ne.ngrid) then !.ne. = not equal
        call copy(z(irhs(j)),z(irho(j)),nn) !Set up r.h.s. 初期の密度を右辺の項に代入  nn=3の時２など？ zを使う(rstrctで用いた値を使う)
     else
        call copy(z(irhs(j)),u,nn) !j最大(ngrid)では元の密度分布を使う
     endif

     do jcycle=1,ncycle !V-cycle loop.
        nf=nn
        do jj=j,2,-1 !Downward stoke of the V.
           do jpre=1,NPRE !Pre-smoothing.ガウスサイデル法の回数
              call relax(z(iu(jj)),z(irhs(jj)),nf) !収束させる。
           enddo

           call resid(z(ires(jj)),z(iu(jj)),z(irhs(jj)),nf) !残差
           nf=nf/2+1
           call rstrct(z(irhs(jj-1)),z(ires(jj)),nf) !粗いメッシュの残差からのrhsを計測

           !Restriction of the residual is the next r.h.s.

           call fill0(z(iu(jj-1)),nf) !Zero for initial guess in next relaxation. 初めに求めたポテンシャルを初期化(iu(j)を求めるために作った)
        enddo
        call slvsml(z(iu(1)),z(irhs(1))) !Bottom of V: solve on coarsest grid. iu(1)が初期化されたのでもう一度呼び出す。
        nf=3
        do jj=2,j !Upward stroke of V.
           nf=2*nf-1
           call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),nf)

           !Use res for temporary storage inside addint.

           do jpost=1,NPOST !Post-smoothing.
              call relax(z(iu(jj)),z(irhs(jj)),nf) !さらに収束させる(残差を小さく)
           enddo
        enddo
     enddo
  enddo
  call copy(u,z(iu(ngrid)),n) !Return solution in u.
  return
END SUBROUTINE mglin

SUBROUTINE rstrct(uc,uf,nc)
  INTEGER nc
  DOUBLE PRECISION uc(nc,nc),uf(2*nc-1,2*nc-1)
  !Half-weighting restriction. nc is the coarse-grid dimension. The ne-grid solution is input
  !in uf(1:2*nc-1,1:2*nc-1), the coarse-grid solution is returned in uc(1:nc,1:nc).
  INTEGER ic,iff,jc,jf
  do  jc=2,nc-1 !Interior points.
     jf=2*jc-1
     do  ic=2,nc-1
        iff=2*ic-1
        uc(ic,jc)=0.5d0*uf(iff,jf)+0.125d0*(uf(iff+1,jf)+ uf(iff-1,jf)+uf(iff,jf+1)+uf(iff,jf-1)) !平均
     enddo
  enddo
  do  ic=1,nc !Boundary points.
     uc(ic,1)=uf(2*ic-1,1)
     uc(ic,nc)=uf(2*ic-1,2*nc-1)
  enddo
  do  jc=1,nc !Boundary points.
     uc(1,jc)=uf(1,2*jc-1)
     uc(nc,jc)=uf(2*nc-1,2*jc-1)
  enddo

  uc((nc+1)/2,(nc+1)/2) = uf((nc+1),(nc+1)) !boundarypoint(center)

  return
END SUBROUTINE rstrct

SUBROUTINE interp(uf,uc,nf)
  INTEGER nf
  DOUBLE PRECISION uc(nf/2+1,nf/2+1),uf(nf,nf)
  INTEGER ic,iff,jc,jf,nc
  !Coarse-to-ne prolongation by bilinear interpolation. nf is the ne-grid dimension. The
  !coarse-grid solution is input as uc(1:nc,1:nc), where nc = nf=2 + 1. The ne-grid
  !solution is returned in uf(1:nf,1:nf).
  nc=nf/2+1
  do  jc=1,nc !Do elements that are copies. 増やしたグリッドに元の値（雑な）を代入
     jf=2*jc-1
     do  ic=1,nc
        uf(2*ic-1,jf)=uc(ic,jc)
     enddo
  enddo

  do jf=1,nf,2 !Do odd-numbered columns, interpolating verdo 内挿
     do iff=2,nf-1,2 !tically.
        uf(iff,jf)=0.5d0*(uf(iff  double precision u(n,n) , cnv(0:10) , sum
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

