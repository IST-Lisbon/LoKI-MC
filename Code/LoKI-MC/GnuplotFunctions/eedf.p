# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/eedf.dat"

set size 0.6,0.9
set origin 0.2,0.05

set log y
set format y "%2.tx10^{%T}";
set title "EEDF" 
set xlabel "Energy (eV)"
set ylabel "EEDF (eV^{-3/2})"
plot "LoKI-MC/GnuplotFunctions/TempData/eedf.dat" using 1:2 title 'EEDF' with lines
