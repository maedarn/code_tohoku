set term x11
 unset key
  set grid x
  set grid y
  set tics font ",15"
  set xlabel  font ',20'
  set ylabel  font ',20'
  set logscale
  set title  font ',20'
  set xrange [0.01:1000]
  set yrange [10:100000]
  file(n) = sprintf("flwb%04d.dat",n)
  plot for [i=1:500] file(i) every ::513:0:513:0 u ($2/1.27):($6*115.8) pt 6 ps 1 lc 3