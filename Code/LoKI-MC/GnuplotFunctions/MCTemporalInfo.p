# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat"
set multiplot                       # get into multiplot mode

set size 0.46,0.5;

set origin 0.27,0.5
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels                   
set xlabel "Time (s)" 
set ylabel "Mean Energy (eV)"
set key outside above
if (steadyStateTime != -1){
	stats "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" u 1:2 nooutput
	set parametric
	set trange [0:STATS_max_y]
	plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:2 title 'Mean Energy' with lines, \
		 steadyStateTime,t title 'Steady-State time'
}
else{
	plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:2 title 'Mean Energy' with lines
}

set origin 0.02,0.0
unset yrange
set xlabel "Time (s)" 
set ylabel "Mean position (m)" 
plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:3 title 'x' with lines, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:4 title 'y' with lines, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:5 title 'z' with lines

set origin 0.52,0.0
set xlabel "Time (s)" 
set ylabel "Squared swarm width (m^2)"
plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:6 title 'x' with lines, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:7 title 'y' with lines, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo.dat" using 1:8 title 'z' with lines


unset multiplot