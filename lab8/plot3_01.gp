# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------
reset
set term gif size 800,300 animate delay 10
set output "res/anim01.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]
set xlabel "X"
set ylabel "Y"

do for [i=0:n] {
  file = sprintf("out/zad5_it=%i_0.1.txt",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 