set term png
set out 'riemannv5.png'
plot 'riemanndata0025.dat' using 1:3 lw 3 w l
set output
set term x11
!open riemannv5.png









