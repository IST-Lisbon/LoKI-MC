# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat"

set multiplot                       # get into multiplot mode
set key outside above

set size 0.46,0.5;

set origin 0.02,0.5
set   autoscale                                      # scale axes automatically
unset log y                                          # remove any log-scaling
set log x
unset label                                          # remove any previous labels
set xlabel label_x 
set ylabel "Reduced transverse diffusion (1/(ms))"
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:2 notitle with linespoints

set origin 0.52,0.5
unset yrange
set ylabel "Reduced mobility (1/(msV))"
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:3 notitle with linespoints

set origin 0.02,0.0
set ylabel "Electron energy (eV)" 
plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:7 title 'Characteristic' with linespoints, \
	 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:6 title 'Mean' with linespoints


set origin 0.52,0.0
set origin 0.52,0.0
# check if all elements in the arrays are zero
stats 'LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat' using 4 nooutput
townsendSum = STATS_sum
stats 'LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat' using 5 nooutput
attachmentSum = STATS_sum

if (townsendSum != 0 && attachmentSum != 0){
	set log y
	set log y2
	set y2tics
	set ytics nomirror
	set ylabel "Red. Townsend Coeff. (m^2)" #offset -1.0,-5
	set y2label "Red. Attachment Coeff. (m^2)" #offset 3.5,-5
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:4 title 'Townsend' with linespoints axis x1y1, \
		 "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:5 title 'Attachment' with linespoints axis x1y2
}
if (attachmentSum == 0){
	set log y
	set ylabel "Red. Townsend Coeff. (m^2)" #offset -1.0,-5
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:4 title 'Townsend' with linespoints
}
if (townsendSum == 0){
	set log y 
	set ylabel "Red. Attachment Coeff. (m^2)" #offset -1.0,-5
	plot "LoKI-MC/GnuplotFunctions/TempData/swarmParams.dat" using 1:5 title 'Attachment' with linespoints
}	 

unset multiplot