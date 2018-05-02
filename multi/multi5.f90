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
  integer(8),PARAMETER :: MEMLEN=5000 !周回を増やす時もきをつけて
  integer, PARAMETER :: NPRE=50,NPOST=1 !ガウスサイデル反復,smoosing
  integer mem
  DOUBLE PRECISION z(MEMLEN)
end module comvar

program main
  implicit none
  INTEGER :: n=32
  integer :: ncycle=5
  double precision u(1:32,1:32),Px(1:32),Py(1:32),uBCx1(1:32),uBCy1(1:32),uBCxn(1:32),uBCyn(1:32)
  call INITIA(u,n)
  call BChazi(uBCx1,uBCy1,uBCxn,uBCyn,n)
  call position(Px,Py,n)
  call mglin(u,n,ncycle,uBCx1,uBCy1,uBCxn,uBCyn)
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
  ! u(16,16) = 100.0d0
  ! u(16,17) = 100.0d0
  ! u(17,16) = 100.0d0
  ! u(17,17) = 100.0d0
   !write(*,*) '**********',u(16,16),'**********'
 write(*,*) '--------------INT----------------'
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
  write(*,*) '--------------pos----------------'
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
  write(*,*) '--------------sav----------------'
end subroutine save


SUBROUTINE mglin(u,n,ncycle,uBCx1,uBCy1,uBCxn,uBCyn)
  use comvar
  INTEGER n,ncycle
  !integer,PARAMETER :: NG=8,MEMLEN=13*2**(2*NG)/3+14*2**NG+8*NG-100/3
  !integer, PARAMETER :: NPRE=1,NPOST=1
  INTEGER j,jcycle,jj,jpost,jpre,nf,ngrid,nn,ires(NG), irho(NG),irhs(NG),iu(NG),maloc
  !DOUBLE PRECISION z
  double precision u(n,n),uBCx1(n),uBCy1(n),uBCxn(n),uBCyn(n)
  !COMMON /memory/ z(MEMLEN),mem  !Storage for grid functions is allocated by maloc
  mem=0 !from array z.
  nn=n/2+1     !nはx,yのメッシュ数,nnは半分の地点
 write(*,*) '----------OK----------'
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
  write(*,*) j
  call interp(res,uc,nf)
  do  j=1,nf
     do  i=1,nf
        uf(i,j)=uf(i,j)+res(i,j) !新しいu(ポテンシャル)を作成
     enddo
  enddo
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
  write(*,*) '**********',u(2,2),'**********'
  u(2,2)= -h*h*rhs(2,2)*0.25d0 !rhsは元の配列(例えば密度) , 初期の長さを１としているのでh=0.5d0  rhs=f:銀本P40 逆行列解ける
 write(*,*) '**********',u(2,2),'**********'
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

