!********************************************************!
!multiglid method     numerical recipe inF77
!********************************************************!
module comvar
  implicit none
  !integer(8),PARAMETER :: NG=5,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-int(100/3)
  integer,PARAMETER :: NG=5,MEMLEN=4437
  integer, PARAMETER :: NPRE=1,NPOST=1
  integer mem
  DOUBLE PRECISION z(MEMLEN)
end module comvar

program main
  implicit none
  INTEGER :: n=32
  integer ncycle
  double precision u(1:32,1:32)
  call INITIA(u,n)
  call mglin(u,n,ncycle)
  call save(u,n)
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
  u(16,16) = 5.0d0      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
end subroutine INITIA

subroutine save(u,n)
  integer n
  double precision u(n,n)
  integer i,j
  !CHARACTER(22) :: dir='/Users/ryunosukemaeda/'
  CHARACTER(32) :: dir='/Users/maeda/Desktop/code/multi/'
  character(4) :: filenm='ml01'
100 format(E19.10e3)
  open(10,file=dir//filenm//'.dat')
  do i = 1,n
     do j=1,n
        write(10,100) u(i,j)
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
  irho(ngrid) = maloc(nn**2)   ! Allocate storage for r.h.s. on grid NG − 1,
  call rstrct(z(irho(ngrid)),u,nn)   !and ll it by restricting from the ne grid. 初期のuにはrhs(密度)が収められている
1 if (nn.gt.3) then   !Similarly allocate storage and ll r.h.s. on all
     nn=nn/2+1 !coarse grids.
     ngrid=ngrid-1
     irho(ngrid)=maloc(nn**2)
     call rstrct(z(irho(ngrid)),z(irho(ngrid+1)),nn) !粗いグリッドと細かいグリッドを引数に取る
     goto 1
  endif
  nn=3
  iu(1)=maloc(nn**2) !1=ngrid
  irhs(1)=maloc(nn**2)
  call slvsml(z(iu(1)),z(irho(1))) !Initial solution on coarsest grid. 初期値uの
  ngrid=NG
  do j=2,ngrid! Nested iteration loop. 粗い位置から始まる (前のv-loopnoの一番下から)
     nn=2*nn-1 !ふやしていく
     iu(j)=maloc(nn**2)
     irhs(j)=maloc(nn**2)
     ires(j)=maloc(nn**2)
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
        uf(iff,jf)=0.5d0*(uf(if+1,jf)+uf(iff-1,jf))
     enddo
  enddo
  do jf=2,nf-1,2 !Do even-numbered columns, interpolating hordo
     do iff=1,nf !izontally.
        uf(iff,jf)=0.5d0*(uf(iff,jf+1)+uf(iff,jf-1))
     enddo
  enddo
  return
END SUBROUTINE interp

SUBROUTINE addint(uf,uc,res,nf)
  INTEGER nf
  DOUBLE PRECISION res(nf,nf),uc(nf/2+1,nf/2+1),uf(nf,nf)
  !C USES interp
  !Does coarse-to-ne interpolation and adds result to uf. nf is the ne-grid dimension. The
  !coarse-grid solution is input as uc(1:nc,1:nc), where nc = nf=2 + 1. The ne-grid
  !solution is returned in uf(1:nf,1:nf). res(1:nf,1:nf) is used for temporary storage.
  INTEGER i,j
  call interp(res,uc,nf)
  do  j=1,nf
     do  i=1,nf
        uf(i,j)=uf(i,j)+res(i,j) !新しいu(ポテンシャル)を作成
     enddo
  enddo
  return
END SUBROUTINE addint

SUBROUTINE slvsml(u,rhs)
  DOUBLE PRECISION rhs(3,3),u(3,3)
  !C USES fill0
  !Solution of the model problem on the coarsest grid, where h = 1
  !2 . The right-hand side is
  !input in rhs(1:3,1:3) and the solution is returned in u(1:3,1:3).
  DOUBLE PRECISION h
  call fill0(u,3) !uの初期化 ここで密度から切り替わる(ポテンシャルに) このサブルーチン内では
  h=0.5d0
  u(2,2)= -h*h*rhs(2,2)/4.d0 !rhsは元の配列(例えば密度) , 初期の長さを１としているのでh=0.5d0  rhs=f:銀本P40 逆行列解ける
  return  !ただし全ての境界のuを0としている
END SUBROUTINE slvsml

SUBROUTINE relax(u,rhs,n) !ガウスサイデル
  INTEGER n
  DOUBLE PRECISION rhs(n,n),u(n,n)
  !Red-black Gauss-Seidel relaxation for model problem. The current value of the solution
  !u(1:n,1:n) is updated, using the right-hand side function rhs(1:n,1:n).
  INTEGER i,ipass,isw,j,jsw
  DOUBLE PRECISION h,h2
  h=1.d0/(n-1) !メッシュ間隔
  h2=h*h
  jsw=1
  do ipass=1,2 ! Red and black sweeps.
     isw=jsw
     do j=2,n-1
        do i=isw+1,n-1,2 !Gauss-Seidel formula.
           u(i,j) = 0.25d0*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-h2*rhs(i,j)) !銀本P40,P49
        enddo
        isw=3-isw
     enddo
     jsw=3-jsw
  enddo
  return
END SUBROUTINE relax

SUBROUTINE resid(res,u,rhs,n)
  INTEGER n
  DOUBLE PRECISION res(n,n),rhs(n,n),u(n,n)
  !Returns minus the residual for the model problem. Input quantities are u(1:n,1:n) and
  !rhs(1:n,1:n), while res(1:n,1:n) is returned.
  INTEGER i,j
  DOUBLE PRECISION h,h2i
  h=1.d0/(n-1)
  h2i=1.d0/(h*h)
  do j=2,n-1 !Interior points.
     do i=2,n-1
        res(i,j)=-h2i*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4.d0*u(i,j))+rhs(i,j) !残差（ガウスサイデル法による）
     enddo
  enddo
  do i=1,n !Boundary points.
     res(i,1)=0.d0
     res(i,n)=0.d0
     res(1,i)=0.d0
     res(n,i)=0.d0
  enddo
  return
END SUBROUTINE resid

SUBROUTINE copy(aout,ain,n)
  INTEGER n
  DOUBLE PRECISION ain(n,n),aout(n,n)
  !Copies ain(1:n,1:n) to aout(1:n,1:n).
  INTEGER i,j
  do i=1,n
     do j=1,n
        aout(j,i)=ain(j,i)
     enddo
  enddo
  return
END SUBROUTINE copy

SUBROUTINE fill0(u,n)
  INTEGER n
  DOUBLE PRECISION u(n,n)
  !Fills u(1:n,1:n) with zeros.
  INTEGER i,j
  do j=1,n
     do i=1,n
        u(i,j)=0.d0
     enddo
  enddo
  return
END SUBROUTINE fill0

FUNCTION maloc(len) !len=格子の個数
  use comvar
  INTEGER maloc,len
  !integer , PARAMETER :: NG=8,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3 !for mglin
  !integer , PARAMETER :: NG=8,MEMLEN=17*2**(2*NG)/3+18*2**NG+10*NG-86/3 !for mgfas, N.B.!
  !INTEGER mem
  !DOUBLE PRECISION z
  !COMMON /memory/ z(MEMLEN),mem
  !Dynamical storage allocation. Returns integer pointer to the starting position for len array
  !elements in the array z. The preceding array element is lled with the value of len, and
  !the variable mem is updated to point to the last element of z that has been used.
  if (mem+len+1.gt.MEMLEN) then
     write(*,*) 'insufficient memory in maloc'
     stop
  else
  z(mem+1)=len
  maloc=mem+2
  mem=mem+len+1
  endif
  return
END FUNCTION maloc
