set xlabel 'X'
set ylabel 'Y'
set title 'Solution'
set key off
set term png
set output "plot.png"
splot 'solution.dat' with lines
pause -1