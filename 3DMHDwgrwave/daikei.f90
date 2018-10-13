program integral_comp
  implicit none
  !c local
  integer m,n ! 刻み数n=2**m
  real h,integral ! 刻み幅h, 数値積分の中間変数
  real trapezoid,simpson ! 台形法とSimpson 法の解
  real x ! 積分変数
  real exact ! 解析解
  integer i ! ループ用変数
  !c function:
  real y ! 被積分関数
  !c begin:
  exact = exp(1.0)-1
  do m=0,6
     n=2**m ! 刻み数
     h=1.0/n ! 刻み幅
     !c 台形法
     integral=y(0.0)+y(1.0) ! 積分端での値
     do i=1,n-1
        x = h*i
        integral=integral+2.0*y(x)
     end do
     trapezoid=integral*h/2.0 ! 台形法の解
     !c Simpson 法
     integral=y(0.0)+y(1.0) ! 積分端での値
     do i=1,n/2
        x = h*(2*i-1)
        integral=integral+4.0*y(x)
     end do
     do i=1,n/2-1
        x = h*(2*i)
        integral=integral+2.0*y(x)
     end do
     simpson=integral*h/3.0 ! Simpson 法の解
     !c 数値積分の結果と解析解の差の出力
     write(*,’(I3,F8.5,4F12.7)’) n,h,
     & trapezoid,trapezoid-exact,
     & simpson,simpson-exact
  end do
  stop
end program integral_comp
!c
!c 被積分関数の定義
!c
real function y(x)
  implicit none
  !c input:
  real x
  !c begin:
  y = exp(x)
  return
end function y
