# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat"

set multiplot                       # get into multiplot mode
set key outside above

set size 0.46,0.5;

set origin 0.02,0.5
set   autoscale                        # scale axes automatically
unset log                               # remove any log-scaling
if (!(xLabel eq "E-angle (º)")){
	set log x
}
unset label                            # remove any previous labels
unset yrange
set xlabel xLabel
set ylabel "Reduced diffusion (1/(ms))"
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:2 title 'Transverse, Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:6 title 'Transverse, Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:3 title 'Longitudinal, Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:7 title 'Longitudinal, Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:13 title 'Eedf' with linespoints

set origin 0.52,0.5
unset yrange
set log y
set ylabel "Reduced mobility (1/(msV))"
stats "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" u 14 nooutput
set yrange [0.8*STATS_min : 1.2*STATS_max]
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:4 title 'Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:8 title 'Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:14 title 'Eedf' with linespoints

set origin 0.02,0.0
set ylabel "Electron temperature (eV)" 
unset yrange
stats "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" u 15 nooutput
min_y = STATS_min;
max_y = STATS_max;
stats "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" u (0.666666667*$10) nooutput
if (STATS_min < min_y) {min_y = STATS_min}
if (STATS_max > max_y) {max_y = STATS_max}
set yrange [0.8*min_y : 1.2*max_y]
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:5 title 'Charac, Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:9 title 'Charac, Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:15 title 'Charac, Eedf' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:(0.666666667*$10) title 'T_e' with linespoints

set origin 0.52,0.0
unset yrange
# check if all elements in the arrays are zero
stats 'LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat' using 11 nooutput
ionizationSum = STATS_sum
stats 'LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat' using 12 nooutput
attachmentSum = STATS_sum

if (ionizationSum != 0 && attachmentSum != 0){
	set log y
	set log y2
	set y2tics
	set ytics nomirror
	set ylabel "Ion. Coeff. (m^{-3})" 
	set y2label "Att. Coeff. (m^{-3})" 
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:11 title 'Ionization' with linespoints axis x1y1, \
		 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:12 title 'Attachment' with linespoints axis x1y2
}
if (attachmentSum == 0){
	set log y
	set ylabel "Ion. Coeff. (m^{-3})" 
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:11 title 'Ionization' with linespoints
}
if (ionizationSum == 0){
	set log y
	set ylabel "Att. Coeff. (m^{-3})" 
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:12 title 'Attachment' with linespoints
}	

unset multiplot