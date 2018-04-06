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
  set yrange [0.01:100000]
plot "b001.dat" every ::500:0:500:0 u 2:6 pt 6 ps 1 lc 3
st = 2
to = 999
offset = 50
do for[i = st : to : offset]{
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
  set yrange [0.01:100000]
    input = sprintf("b%03d.dat", i)
   plot input every ::500:0:500:0 u 2:6 pt 6 ps 1 lc 3
   }
