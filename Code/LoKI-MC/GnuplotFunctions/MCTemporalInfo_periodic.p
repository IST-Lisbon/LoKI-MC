# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat"
set multiplot                       # get into multiplot mode

set size 0.46,0.5;

set origin 0.27,0.5
set   autoscale                        # scale axes automatically
unset log                              # remove any log-scaling
unset label                            # remove any previous labels                   
set key outside above
set y2tics
set ytics nomirror
set xlabel "Phase (rad)" 
set ylabel "Mean Energy (eV)"
set y2label "E/N (Td)"
set xrange [0:2*pi]
plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:4 title 'Mean Energy' with lines axis x1y1, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:3 title 'E/N' with lines axis x1y2
unset y2tics
unset y2label

set origin 0.02,0.0
set ylabel "Velocity (m/s)" 
plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:5 title 'Flux, x' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:6 title 'Flux, y' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:7 title 'Flux, z' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:8 title 'Bulk, x' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:9 title 'Bulk, y' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:10 title 'Bulk, z' with linespoints

set origin 0.52,0.0
set ylabel "N D (1/(ms)" 
plot "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:11 title 'Flux, xx' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:15 title 'Flux, yy' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:19 title 'Flux, zz' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:20 title 'Bulk, xx' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:24 title 'Bulk, yy' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/MCTemporalInfo_periodic.dat" using 1:28 title 'Bulk, zz' with linespoints


unset multiplot