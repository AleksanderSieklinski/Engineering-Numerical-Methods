set terminal pngcairo size 600,500 enhanced font 'Verdana,10'
set lmargin at screen 0.15
set rmargin at screen 0.75
set output 'res/temperature_0100.png'
set title 'temperature at it = 100'
set xlabel 'x'
set ylabel 'y'
set pm3d map
splot 'out/temperature_0100.dat' using 1:2:3 with pm3d

set output 'res/temperature_0200.png'
set title 'temperature at it = 200'
splot 'out/temperature_0200.dat' using 1:2:3 with pm3d

set output 'res/temperature_0500.png'
set title 'temperature at it = 500'
splot 'out/temperature_0500.dat' using 1:2:3 with pm3d

set output 'res/temperature_1000.png'
set title 'temperature at it = 1000'
splot 'out/temperature_1000.dat' using 1:2:3 with pm3d

set output 'res/temperature_2000.png'
set title 'temperature at it = 2000'
splot 'out/temperature_2000.dat' using 1:2:3 with pm3d