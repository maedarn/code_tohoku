reset
st = 1
to = 999
offset = 1
do for[i = st : to : offset]{
    input = sprintf("HPba%03d.dat", i)
    set term  png
    set out sprintf("HPbap%03d.png", i)
    set style line 1 lc rgb "red" lw 3
    set style line 1 pointtype 6 ps 5
unset key
set grid x
  set grid y
  set tics font ",15"
  set xlabel  font ',20'
  set ylabel  font ',20'
  set title  font ',20'
  set xrange [0:600]
  set yrange [0:100]
    plot input using 1:6 w l lw 3
    set output
}
reset