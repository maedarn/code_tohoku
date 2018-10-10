unset key
set size square
set tics font "Times New Roman,30"
set xlabel  font 'Times New Roman,30'
set ylabel  font 'Times New Roman,30'
set title  font 'Times New Roman,30'
set xlabel offset 0,-0.5
set ylabel offset -3,0
set cblabel offset 6,0
set xrange [0:100]
set yrange [0:100]
set ticslevel 0
set pm3d
set pm3d map
set logscale cb
set format cb "10^{%L}"
set cbrange[1:1000]
set palette defined(1"#ff0000",2"#ff8000",3"#ffff00",4"#80ff00",5"#00ff00",6"#00ff80",7"#00ffff",8"#0080ff",9"#0000ff")
