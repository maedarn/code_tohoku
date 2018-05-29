
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
  !integer,PARAMETER :: NG=9
  integer,PARAMETER :: NG=5
  !integer(8),PARAMETER :: MEMLEN=1500000 !周回を増やす時もきをつけて
  integer(8),PARAMETER :: MEMLEN=5000 !周回を増やす時もきをつけて
  integer mem
  DOUBLE PRECISION z(MEMLEN)
end module comvar

program main
  implicit none
  !INTEGER :: n=513,tsave=0,itr
  INTEGER :: n=33,tsave=0,itr,itr2,itr3,itr4
  integer :: ncycle=1
  integer :: NPRE=1,NPOST=1 !ガウスサイデル反復,smoosing
  !double precision u(1:513,1:513),Px(1:513),Py(1:513),uBCx1(1:513),uBCy1(1:513),uBCxn(1:513),uBCyn(1:513)
  double precision u(1:33,1:33),Px(1:33),Py(1:33),uBCx1(1:33),uBCy1(1:33),uBCxn(1:33),uBCyn(1:33)
  character(4) chr
  !write(*,*) 'aaaa'
  call position(Px,Py,n)
  do itr4=1,5
     do itr3=1,5
        do itr2=1,5
           call INITIA(u,n)
           !call BChazi(uBCx1,uBCy1,uBCxn,uBCyn,n)
           do itr=1,10
              call save(u,Px,Py,n,tsave,itr2,itr3,itr4)
              !write(*,*) 'aaaa'
              call mglin(u,n,ncycle,uBCx1,uBCy1,uBCxn,uBCyn,NPRE,NPOST)
              !write(*,*) 'aaaa'
           end do
           !call position(Px,Py,n)
           call save(u,Px,Py,n,tsave,itr2,itr3,itr4)
           ncycle=ncycle+1
           tsave = 0
        end do
        NPRE=1+NPRE
     end do
     NPOST=1+NPOST
  end do
end program

subroutine INITIA(u,n)
  integer n
  double precision u(n,n)
  integer i,j
  !double precision,intent(out) :: rnd
  integer :: seedsize
  integer, allocatable :: seed(:,:)

  seedsize=n

  do i = 1,n
     do j=1,n
        u(i,j) = 0.0d0
     end do
  end do
  !*****************************************************
  !call random_seed(size=n**2) !初期値のサイズを取得
  !allocate(seed(seedsize,seedsize)) !配列の割り当て
  !do i = 1, seedsize
  !   do j=1,seedsize
  !      call system_clock(count=seed(i,j)) !時間を取得
  !   end do
  !end do
  !call random_seed(put=seed(:,:)) !初期値を与える

  !do i = 1,n
  !   do j=1,n
  !      call random_number(u)
  !      u(i,j)=u*dble(n)*dble(n)
  !   end do
  !end do

  call random_number(u(:,:))
  !u(:,:)=(u(:,:)-0.5d0)*dble(n)*dble(n)
  u(:,:)=(u(:,:)-0.5d0)*10.0d0
  !*****************************************************
  !u(16,16) = 100.0d0
  !u(16,17) = 100.0d0
  !u(17,16) = 100.0d0
  !u(257,257) = 100.0d0
  !u(17,17) = 100.0d0
  !write(*,*) '**********',u(16,16),'**********'
  !write(*,*) '--------------INT----------------'
end subroutine INITIA

subroutine position(Px,Py,n)
  integer n
  integer i
  double precision hx,hy
  double precision Px(n),Py(n)
  hx=1.0d0/dble(n)
  hy=1.0d0/dble(n)
  !write(*,*) n,'99999999999'
  do i=1,n
     Px(i) = hx*i
  end do
  do i=1,n
     Py(i) = hy*i
  end do
  ! write(*,*) '--------------pos----------------'
end subroutine position


subroutine save(u,Px,Py,n,tsave,itr2,itr3,itr4)
  integer n,tsave,itr2,itr3,itr4
  double precision  u(n,n)
  double precision :: gosa(-1:100)=0.0d0 !itrの回数
  double precision Px(n),Py(n)
  integer i,j
  !CHARACTER(22) :: dir='/Users/ryunosukemaeda/'
  CHARACTER(32) :: dir='/Users/maeda/Desktop/code/multi/'
  character(4) :: filenm='ml-1',chr
  character(4) :: filenm2='mlg1'
  character(2) sv
  character(1) i2,i3,i4
  write(sv,'(I2.2)') tsave
  write(i2,'(I1.1)') itr2
  write(i3,'(I1.1)') itr3
  write(i4,'(I1.1)') itr4
  if(tsave==0)then
     gosa(:)0.0d0
  end if

100 format(E19.10e3,E19.10e3,E19.10e3)
  open(10,file=dir//filenm//i4//i3//i2//sv//'.dat')
  do i = 1,n
     do j=1,n
        write(10,100) Px(i) , Py(j) , u(i,j)
        gosa(tsave)=u(i,j)**2+gosa(tsave)
     end do
     write(10,*)
  end do
  ! write(*,*) '--------------sav----------------'
  open(11,file=dir//filenm2//i4//i3//i2//'.dat')
  write(11,100)  -gosa(tsave)+gosa(tsave-1)
  tsave=tsave+1
end subroutine save


SUBROUTINE mglin(u,n,ncycle,uBCx1,uBCy1,uBCxn,uBCyn,NPRE,NPOST)
  use comvar
  INTEGER n,ncycle,NPRE,NPOST
  !integer,PARAMETER :: NG=8,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3
  !integer, PARAMETER :: NPRE=1,NPOST=1
  INTEGER j,jcycle,jj,jpost,jpre,nf,ngrid,nn,ires(NG), irho(NG),irhs(NG),iu(NG),maloc
  !DOUBLE PRECISION z
  double precision u(n,n),uBCx1(n),uBCy1(n),uBCxn(n),uBCyn(n)
  !COMMON /memory/ z(MEMLEN),mem  !Storage for grid functions is allocated by maloc
  mem=0 !from array z.
  nn=n/2+1     !nはx,yのメッシュ数,nnは半分の地点
  ngrid=NG-1 !レベル2**Nみたいなもの（何回イタレーションするか）
  irho(ngrid) = maloc(int(nn**2))   ! Allocate storage for r.h.s. on grid NG − 1,
  !write(*,*) '**********',u(8,8),'**********'
  !z(irho(ngrid)+1)=1.0d0
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
     nn=2*nn-1 !ふやしていく 細かく
     write(*,*) j
     iu(j)=maloc(int(nn**2))
     irhs(j)=maloc(int(nn**2))
     ires(j)=maloc(int(nn**2))
     ! write(*,*) j
     call interp(z(iu(j)),z(iu(j-1)),nn) !Interpolate from coarse grid to next ner grid. zに格納（挿入）細かく

     if (j.ne.ngrid) then !.ne. = not equal
        call copy(z(irhs(j)),z(irho(j)),nn) !Set up r.h.s. 初期の密度を右辺の項に代入  zを使う(rstrctで用いた値を使う)
     else
        call copy(z(irhs(j)),u,nn) !j最大(ngrid)では元の密度分布を使う
     endif

     do jcycle=1,ncycle !V-cycle loop.
        nf=nn
        write(*,*) 'nf=' , nf
        do jj=j,2,-1 !Downward stoke of the V.  jはグリットの荒さ(どんどん細かく)
           do jpre=1,NPRE !Pre-smoothing.ガウスサイデル法の回数
              call relax(z(iu(jj)),z(irhs(jj)),nf) !収束させる。
           enddo

           call resid(z(ires(jj)),z(iu(jj)),z(irhs(jj)),nf) !残差
           nf=nf/2+1
           call rstrct(z(irhs(jj-1)),z(ires(jj)),nf) !粗いメッシュの残差からのrhsを計測(残差による湧き出しの項)

           !Restriction of the residual is the next r.h.s.

           call fill0(z(iu(jj-1)),nf) !Zero for initial guess in next relaxation. 初めに求めたポテンシャルを初期化(iu(j)を求めるために作った)
        enddo
        !write(*,*) j
        call slvsml(z(iu(1)),z(irhs(1))) !Bottom of V: solve on coarsest grid. iu(1)が初期化されたのでもう一度呼び出す。
        !write(*,*) j
        nf=3
        do jj=2,j !Upward stroke of V.
           nf=2*nf-1 !グリッドの目
           write(*,*) j
           call addint(z(iu(jj)),z(iu(jj-1)),z(ires(jj)),nf) !粗いメッシュの残差が作ったポテンシャルを細かいグリッドの新たなzに足す。ただし初期化しているのはu(NG-1)までなので最終的なuは正しい値になる。
           write(*,*) j
           !Use res for temporary storage inside addint.

           do jpost=1,NPOST !Post-smoothing.
              !write(*,*) z(iu(jj)+nf**2-2*nf) ,'!!!!!!!!!!!!!!!'
              call relax(z(iu(jj)),z(irhs(jj)),nf) !さらに収束させる(残差を小さく)これは残差の作るポテンシャルを残差の作る密度で収束させている。
              !write(*,*) z(iu(jj)+nf**2-2*nf) ,'!!!!!!!!!!!!!!!'
           enddo
           !call converge(z(iu(j)),nf)
        enddo
        !write(*,*) j
     enddo
     !call converge(z(iu(j)),nf)
  enddo
  !write(*,*) u(n,n/2),'????????????????????????'
  !write(*,*) z(iu(ngrid)) ,'!!!!!!!!!!!!!!!'
  !write(*,*) z(iu(ngrid)+n**2-n/2),'!!!!!!!!!!!!!!!'
  call copy(u,z(iu(ngrid)),n) !Return solution in u.
  !write(*,*) u(n,n/2),'????????????????????????'
  return
END SUBROUTINE mglin

SUBROUTINE rstrct(uc,uf,nc)
  INTEGER nc
  DOUBLE PRECISION uc(nc,nc),uf(2*nc-1,2*nc-1)
  !Half-weighting restriction. nc is the coarse-grid dimension. The ne-grid solution is input
  !in uf(1:2*nc-1,1:2*nc-1), the coarse-grid solution is returned in uc(1:nc,1:nc).
  INTEGER :: ic,iff,jc,jf,count=1
  !write(*,*) uc(1,2),uc(2,1) ,'888888888888888888888888888' 後ろから格納
  !uc(1,2)=0.d0
  !uc(2,1)=0.d0
  if(count==1) call INITIA(uf,2*nc-1)
  count=count+1
  do  jc=2,nc-1 !Interior points.
     jf=2*jc-1
     do  ic=2,nc-1
        iff=2*ic-1
        uc(ic,jc)=0.5d0*uf(iff,jf)+0.125d0*(uf(iff+1,jf)+ uf(iff-1,jf)+uf(iff,jf+1)+uf(iff,jf-1)) !平均
     enddo
  enddo
  do  ic=1,nc !Boundary points. 境界では残差はもともと０
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
  DOUBLE PRECISION uc(int(nf/2)+1,int(nf/2)+1),uf(nf,nf)
  INTEGER ic,iff,jc,jf,nc
  !Coarse-to-ne prolongation by bilinear interpolation. nf is the ne-grid dimension. The
  !coarse-grid solution is input as uc(1:nc,1:nc), where nc = nf=2 + 1. The ne-grid
  !solution is returned in uf(1:nf,1:nf).
  nc=int(nf/2)+1
  !write(*,*) '----------OK----------'
  !***********************BC**************
  uc(1,:) = 0.0d0
  uc(nc,:) = 0.0d0
  uc(:,1) = 0.0d0
  uc(:,nc) = 0.0d0
  !***********************BC**************
  ! write(*,*) '----------OK----------'
  do  jc=1,nc !Do elements that are copies. 増やしたグリッドに元の値（雑な）を代入
     jf=2*jc-1
     do  ic=1,nc
        uf(2*ic-1,jf)=uc(ic,jc)
     enddo
  enddo

  do jf=1,nf,2 !Do odd-numbered columns, interpolating verdo 内挿
     do iff=2,nf-1,2 !tically.
        uf(iff,jf)=0.5d0*(uf(iff+1,jf)+uf(iff-1,jf))
     enddo
  enddo
  do jf=2,nf-1,2 !Do even-numbered columns, interpolating hordo
     do iff=1,nf !izontally.
        uf(iff,jf)=0.5d0*(uf(iff,jf+1)+uf(iff,jf-1))
     enddo
  enddo
  !write(*,*) '----------OK----------'
  !***********************BC**************
  uf(1,:) = 0.0d0
  uf(nf,:) = 0.0d0
  uf(:,1) = 0.0d0
  uf(:,nf) = 0.0d0
  !***********************BC**************
  !write(*,*) '----------OK----------'
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
  !write(*,*) j
  !write(*,*) '--------------------',uf(nf/2,nf/2),uf(nf,nf),'--------------------'
  call interp(res,uc,nf)
  ! write(*,*) '--------------------',uf(nf/2,nf/2),uf(nf,nf),'--------------------'
  do  j=1,nf
     do  i=1,nf
        uf(i,j)=uf(i,j)+res(i,j) !新しいu(ポテンシャル)を作成
     enddo
  enddo
  !write(*,*) '--------------------',uf(nf/2,nf/2),uf(nf,nf),'--------------------'
  return
END SUBROUTINE addint

SUBROUTINE slvsml(u,rhs) !式19.0.6
  DOUBLE PRECISION rhs(3,3),u(3,3)
  !C USES fill0
  !Solution of the model problem on the coarsest grid, where h = 1
  !2 . The right-hand side is
  !input in rhs(1:3,1:3) and the solution is returned in u(1:3,1:3).
  DOUBLE PRECISION h
  call fill0(u,3) !uの初期化 ここで密度から切り替わる(ポテンシャルに) このサブルーチン内では
  h=0.5d0
  !write(*,*) '**********',u(2,2),'**********'
  u(2,2)= -h*h*rhs(2,2)*0.25d0 !rhsは元の配列(例えば密度) , 初期の長さを１としているのでh=0.5d0  rhs=f:銀本P40 逆行列解ける
  !write(*,*) '**********',u(2,2),'**********'
  !u(2,2) = 5.0d0 !境界条件(potential center)
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
        isw=3-isw !偶奇交互に
     enddo
     jsw=3-jsw
  enddo
  ! write(*,*) u(n,n/2),'????????????????????????'
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
  !write(*,*) ain(1,1),'77777777777777777777777'
  !write(*,*) ain(n,n/2),'77777777777777777777777'
  do i=1,n
     do j=1,n
        aout(j,i)=ain(j,i)
     enddo
  enddo
  !write(*,*) aout(n,n/2),'77777777777777777777777'
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
  write(*,*) 'len=',len
  if (mem+len+1.gt.MEMLEN) then
     write(*,*) 'insufficient memory in maloc'
     stop
  else
     z(mem+1)=len
     maloc=mem+2
     mem=mem+len+1
     write(*,*)'mem=', mem
  endif
  return
END FUNCTION maloc

subroutine BC(u,n) !初めの境界の値って覚えさせておいたほうがいいのか？
  integer n
  integer i
  double precision u(n,n)
  integer :: a=1,b=0
  !BC 端
  if(a==1) then !固定端
     do i=1,n
        u(1,i) = 0.0d0
        u(n,i) = 0.0d0
     enddo
     do i=1,n
        u(i,1) = 0.0d0
        u(i,n) = 0.0d0
     enddo
  else if(a==3) then !ノイマン
     do i=1,n
        u(1,i) = u(2,i)
        u(n,i) = u(n-1,i)
     enddo
     do i=1,n
        u(i,1) = u(i,2)
        u(i,n) = u(i,n-1)
     enddo
  else !周期境界
     do i=1,n
        u(1,i) = u(n,i)
        u(n,i) = u(1,i)
     enddo
     do i=1,n
        u(i,1) = u(i,n)
        u(i,n) = u(i,1)
     enddo
  end if
  !BC 中心
  if(b>1) then
     if(b==2) then !固定端
        u(n/2,n/2) = 0.0d0
        u(n/2+1,n/2) = 0.0d0
        u(n/2,n/2+1) = 0.0d0
        u(n/2+1,n/2+1) = 0.0d0
     else if(b==3) then !ノイマン 中心が平坦(傾きを持たすようにしてもいい)
        u(n/2+1,n/2) = u(n/2,n/2)
        u(n/2,n/2+1) =  u(n/2,n/2)
        u(n/2+1,n/2+1) =  u(n/2,n/2)
     else !フリー？
        u(n/2,n/2) = u(n/2-1,n/2-1)
        u(n/2,n/2+1) = u(n/2-1,n/2+2)
        u(n/2+1,n/2) =  u(n/2+2,n/2-1)
        u(n/2+1,n/2+1) =  u(n/2+2,n/2+2)
     end if
  end if
end subroutine BC

subroutine converge(u,n)
  integer n
  double precision u(n,n) , cnv(0:50) , sum
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

subroutine BChazi(uBCx1,uBCy1,uBCxn,uBCyn,n)
  integer n
  integer i
  double precision uBCx1(n),uBCy1(n),uBCxn(n),uBCyn(n)  !yとxが重なるところではyを0としておく。
  do i=1,n
     uBCx1(i)=0.0d0
  end do
  do i=1,n
     uBCxn(i)=0.0d0
  end do
  do i=1,n
     uBCy1(i)=0.0d0
  end do
  do i=1,n
     uBCyn(i)=0.0d0
  end do
end subroutine BChazi

