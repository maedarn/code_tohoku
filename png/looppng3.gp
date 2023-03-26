n=1
if(n==1)\
set term  png ;\
 set out 'cool.png' ;\
    set style line 1 lc rgb "red" lw 3 ;\
    set style line 1 pointtype 6 ps 5 ;\
  unset key;\
  set grid x;\
  set grid y;\
  set tics font ",15";\
  set xlabel  font ',20';\
  set ylabel  font ',20';\
  set logscale;\
  set title  font ',20';\
  set xrange [0.01:1000];\
  set yrange [0.01:100000];\
  m=2;\
  n=1+1;\
  i=0;\
  st=1;\
plot "b001.dat"  u 2:6 pt 6 ps 1 lc 3 ;\
else;\
if(i <= 999)
i = st + 1;\
    input = sprintf("b%03d.dat", i);\
   replot input every ::500:0:500:0 u 2:6 pt 6 ps 1 lc 3;\
m = m+1;\
reread;\
if(m==1000) set output;
