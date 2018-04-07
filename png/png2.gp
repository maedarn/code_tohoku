set term png
set out '001.png'
plot '001.DAT' using 1:9 lw 3 w l
set output
set term x11
!open 001.png









