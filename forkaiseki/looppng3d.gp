reset
st = 0
to = 30
offset = 1
do for[i = st : to : offset]{
input = sprintf("2Dall%03d.dat", i)
set term  png
set out sprintf("2Dall%03d.png", i)
unset key
set grid x
set grid y
set size square
set tics font ",15"
set xlabel  font ',20'
set ylabel  font ',20'
set title  font ',20'
set xrange [0:20]
set yrange [0:20]
set ticslevel 0
set pm3d map
set palette defined (-3 "blue", 0 "white", 1 "red")
splot input using 1:2:3
set output
}
reset