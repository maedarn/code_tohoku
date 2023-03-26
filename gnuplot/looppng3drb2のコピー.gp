reset
st = 0
to = 100
offset = 1
do for[i = st : to : offset]{
input = sprintf("phi%05d.dat", i)
set term  png
set out sprintf("phi-equi%05d.png", i)
unset key
set grid x
set grid y
set size square
set tics font ",15"
set xlabel  font ',20'
set ylabel  font ',20'
set title  font ',20'
set xrange [0:100]
set yrange [0:100]
set ticslevel 0
set pm3d map
set cbrange[-0.1:1.0]
set palette defined(1"#ff0000",2"#ff8000",3"#ffff00",4"#80ff00",5"#00ff00",6"#00ff80",7"#00ffff",8"#0080ff",9"#0000ff")
splot input using 1:2:11
set output
}
reset
