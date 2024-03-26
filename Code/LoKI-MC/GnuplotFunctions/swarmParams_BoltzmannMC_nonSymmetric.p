# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat"

set multiplot                       # get into multiplot mode
set key outside above samplen 5

set size 0.46,0.5;

set origin 0.02,0.5
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
if (!(xLabel eq "E-angle (º)")){
	set log x
}
unset label                            # remove any previous labels
set xlabel xLabel
set ylabel "Reduced diffusion (1/(ms))"
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:2 title 'xx, Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:5 title 'xx, Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:3 title 'yy, Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:6 title 'yy, Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:4 title 'zz, Flux' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:7 title 'zz, Bulk' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:17 title 'Eedf' with linespoints


set origin 0.52,0.5
unset yrange
set log y
set ylabel "Absol. Drift Velocity (m/s)"
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:8 title 'Flux, x' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:9 title 'Flux, y' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:10 title 'Flux, z' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:11 title 'Bulk, x' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:12 title 'Bulk, y' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:13 title 'Bulk, z' with linespoints	 

set origin 0.02,0.0
set ylabel "Electron Temperature (eV)" 
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:19 title 'Charac, Eedf' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:(0.666666667*$14) title 'T_e' with linespoints

set origin 0.52,0.0
# check if all elements in the arrays are zero
stats 'LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat' using 15 nooutput
ionizationSum = STATS_sum
stats 'LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat' using 16 nooutput
attachmentSum = STATS_sum

if (ionizationSum != 0 && attachmentSum != 0){
	set log y
	set log y2
	set y2tics
	set ytics nomirror
	set ylabel "Ion. Coeff. (m^{-3})" 
	set y2label "Att. Coeff. (m^{-3})" 
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:15 title 'Ionization' with linespoints axis x1y1, \
		 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:16 title 'Attachment' with linespoints axis x1y2
}
if (attachmentSum == 0){
	set log y
	set ylabel "Ion. Coeff. (m^{-3})" 
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:15 title 'Ionization' with linespoints
}
if (ionizationSum == 0){
	set log y
	set ylabel "Att. Coeff. (m^{-3})" 
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:16 title 'Attachment' with linespoints
}

unset multiplot