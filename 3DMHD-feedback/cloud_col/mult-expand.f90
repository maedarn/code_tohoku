program multipole_expansion
  implicit none
  
  ! 必要なライブラリをインクルードする
  include 'mpif.h'
  include 'math.h'
  
  ! その他の変数を宣言する
  integer :: rank, size, ierr
  integer :: n, m, l, i, j, k
  real(kind=8) :: r, theta, phi
  real(kind=8), dimension(3) :: x, multipole
  
  ! MPIの初期化
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  
  ! メインの計算ルーチン
  do i = rank + 1, n, size
    do j = 1, m
      do k = 1, l
        r=dsqrt(x(1)**2+x(2)**2+x(3)**2)
        theta=arccos(x(3)/r)
        phi=dsign(1.d0,x(2))*arccos(x(1)/dsqrt(x(1)**2+x(2)**2))
        ! x座標、y座標、z座標を計算する
        x(1) = r * sin(theta) * cos(phi)
        x(2) = r * sin(theta) * sin(phi)
        x(3) = r * cos(theta)
        
        ! 多重極展開を計算する
        call multipole_expansion(x, multipole)
      end do
    end do
  end do
  
  ! MPIの終了処理
  call MPI_Finalize(ierr)
end program multipole_expansion

subroutine multipole_expansion(x, multipole)
    implicit none
    double precision Ylm
    complex ephi
    
    ! 入力引数
    real(kind=8), dimension(3), intent(in) :: x
    
    ! 出力引数
    real(kind=8), dimension(3), intent(out) :: multipole

    Ylm=dsqrt(dble(2*l+1)/(4.d0*pi)*dble(factorial(l-m))/dble(factorial(l+m)))*plgndr(l,m,x)
    ephi1=DCMPLX(Ylm*dcos(m*phi),Ylm*dsin(m*phi))
    ehpi2=DCONJG(ephi1)
    
    
    do l=0,lmax
    multipole=-G*4.d0*pi/(2.d0*dble(l)+1
    enddo
    
end subroutine multipole_expansion


FUNCTION plgndr(l,m,x)
INTEGER l,m
double precision plgndr,x
INTEGER i,ll
double precision fact,pll,pmm,pmmp1,somx2

if((m.lt.0).or.(m.gt.l).or.(abs(x).gt.1.d0)) goto 295 !'bad arguments in plgndr'

pmm=1.d0
 if(m.gt.0) then
 somx2=sqrt((1.d0-x)*(1.d0+x))
 fact=1.d0
  do  i=1,m
   pmm=-pmm*fact*somx2 !numerical receipt (6.8.8)
   fact=fact+2.d0
  enddo
 endif
 if(l.eq.m) then
 plgndr=pmm
 else
 pmmp1=x*dble(2*m+1)*pmm !numerical receipt (6.8.9)
 m+1
  if(l.eq.m+1) then
  plgndr=pmmp1
  else
   do ll=m+2,l
    pll=(x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm)/dble(ll-m) !numerical receipt (6.8.7)
    pmm=pmmp1
    pmmp1=pll
   enddo
  plgndr=pll
  endif
 endif

295 continue
return
END FUNCTION plgndr

recursive function factorial(n) result(factorial_n)
    implicit none
    integer(int32), intent(in) :: n
    integer(int64) :: factorial_n

    if (n > 0) then
        factorial_n = n*factorial(n - 1)
        return
   end if

   factorial_n = 1
end function


subroutine centroid()

do i=1,nx; do j=1,ny; do k=1,nz
xg=rho(i,j,k)*x(i)+xg
yg=rho(i,j,k)*y(j)+yg
zg=rho(i,j,k)*z(k)+zg

mass=rho(i,j,k)+mass
enddo; enddo; enddo
xg=xg/mass
yg=yg/mass
zg=zg/mass
end subroutine centroid
