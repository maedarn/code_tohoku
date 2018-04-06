set terminal png
set output "a001.png"
set xlabel "position [x]"
set ylabel "2-rho , 3.4.5-v , 6-p , 7.8.9-B"
set title "alfven"
set xrange [ 0 : 10 ]
set yrange [0.9999990:1.0000020]
set mxtics 5
set mytics 5
set xtics 1
set ytics 0.50
plot 'a001.DAT' using 1:2 lw 3  w l