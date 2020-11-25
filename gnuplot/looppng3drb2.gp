reset
st = 1
to = 1000
offset = 1
do for[i = st : to : offset]{
input = sprintf("phi%05d.dat", i)
set term  png
set out sprintf("phi%05d.png", i)
unset key
set grid x
set grid y
set tics font ",15"
set xlabel  font ',20'
set ylabel  font ',20'
set title  font ',20'
set xrange [0:100]
set yrange [-0.00000005:0.00000005]
plot input using 1:3 w l
set output
}
reset