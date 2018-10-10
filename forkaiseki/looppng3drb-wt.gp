reset
st = 0
to = 200
offset = 1
do for[i = st : to : offset]{
n=i+1
timefile = "timein.DAT"
time = system("cat " . timefile . " | awk \'NR==" . n . "{printf(\"%f\", $1" . ")}\'")
time = time*1
input = sprintf("2Dall%03d.dat", i)
set term  png
set out sprintf("2Dall%03d.png", i)
unset key
#set grid x
#set grid y
set size square
set tics font ",15"
set xlabel  font ',20'
set ylabel  font ',20'
set title  font ',20'
set title sprintf("time = %.5f Myr",time)
set xrange [0:10]
set yrange [0:10]
set ticslevel 0
set pm3d map
set logscale cb
set cbrange[1:1000]
set palette defined(1"#ff0000",2"#ff8000",3"#ffff00",4"#80ff00",5"#00ff00",6"#00ff80",7"#00ffff",8"#0080ff",9"#0000ff")
splot input using 1:2:3
set output
}
reset