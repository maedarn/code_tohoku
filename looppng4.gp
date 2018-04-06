reset
set term  x11

    set style line 1 lc rgb "red" lw 3
    set style line 1 pointtype 6 ps 5
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
plot "b001.dat" every ::512:0:512:0  u 2:6 pt 6 ps 1 lc 3
st = 2
to = 999
offset = 10
do for[i = st : to : offset]{
    input = sprintf("b%03d.dat", i)
   replot input every ::512:0:512:0 u 2:6 pt 6 ps 1 lc 3
   }

reset

