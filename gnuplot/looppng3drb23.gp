reset
st = 2
to = 16
offset = 1
do for[i = st : to : offset]{
set term  png
set out sprintf("time-lp%03d.png", i-1)
set logscale y
set logscale x
set grid x
set grid y
set xlabel  font 'Times,30'
set ylabel  font 'Times,30'
set tics font "Times,20"
set format y "10^{%L}"
set format x "10^{%L}"
set title  font ',30'
set ylabel offset -2,0
set xlabel offset 0,-6
set lmargin 8
set tmargin 2
set bmargin 6
set size square
unset key
set xrange [100:10000000]
set yrange [100:10000000]
j=i+16
k=j+16
plot 'tffbind.DAT' using i:j w lp, 'tffbind.DAT' using i:k w lp,x
set output
}
reset