# Gnuplot script file for plotting data in file "LoKI-MC/GnuplotFunctions/TempData/power.dat"

set key outside above

set size 0.92,0.95
set origin 0.04,0.04

unset log y
if (!(xLabel eq "E-angle (º)")){
	set log x
}
set xlabel xLabel
set ylabel "Normalized power"
if (!(xLabel eq "Time (s)")){
	set yrange [-1:1]
}
plot "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:2 title 'Field' with linespoints linetype rgb "black" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:3 title '  Elastic' with linespoints linetype rgb "red" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:4 notitle with linespoints linetype rgb "red" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:7 title '  Rotations' with linespoints linetype rgb "forest-green" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:8 notitle with linespoints linetype rgb "forest-green" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:9 title '  Vibrations' with linespoints linetype rgb "blue" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:10 notitle with linespoints linetype rgb "blue" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:11 title ' Electronic' with linespoints linetype rgb "mediumpurple3" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:12 notitle with linespoints linetype rgb "mediumpurple3" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:13 title ' Ionization' with linespoints linetype rgb "brown" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:14 title ' Attachment' with linespoints linetype rgb "aquamarine" linewidth 2,\
	 "LoKI-MC/GnuplotFunctions/TempData/power.dat" using 1:15 title ' Growth' with linespoints linetype rgb "gold" linewidth 2

