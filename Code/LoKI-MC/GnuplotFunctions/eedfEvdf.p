# Gnuplot script file for plotting data in the files "LoKI-MC/GnuplotFunctions/TempData/eedf.dat" and "LoKI-MC/GnuplotFunctions/TempData/evdf.dat"

set multiplot                    # get into multiplot mode

set size 0.48,0.8

set origin 0.02,0.1
set log y
set format y "%2.tx10^{%T}";
set title "EEDF"
set xlabel "Energy (eV)"
set ylabel "EEDF (eV^{-3/2})"
plot "LoKI-MC/GnuplotFunctions/TempData/eedf.dat" using 1:2 title 'EEDF' with lines enhanced 

set origin 0.5,0.1
unset log x
unset log y
unset format
set title "EVDF" 
set xlabel "v_r(m/s)" offset 0,graph -0.05
set ylabel "v_z (m/s)" 
#set pm3d
#set palette rgbformulae 22,13,-31
set xtics offset 0,graph -0.04 out nomirror
set mxtics
show mxtics
set ytics out nomirror
set mytics
show mytics
plot "LoKI-MC/GnuplotFunctions/TempData/evdf.dat" using 1:2:3
clear
plot [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX*2] [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX] "LoKI-MC/GnuplotFunctions/TempData/evdf.dat" using 1:2:3 title 'EVDF' with image enhanced 

unset multiplot