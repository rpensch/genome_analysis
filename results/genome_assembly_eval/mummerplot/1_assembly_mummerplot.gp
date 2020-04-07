set terminal x11 font "Courier,8"
set xtics rotate ( \
 "tig00000037" 1.0, \
 "tig00000041" 222791.0, \
 "tig00000042" 262700.0, \
 "tig00000044" 290263.0, \
 "tig00000046" 306313.0, \
 "tig00000047" 319437.0, \
 "tig00000051" 335917.0, \
 "tig00000100" 345716.0, \
 "tig00000141" 348174.0, \
 "tig00000151" 351936.0, \
 "tig00000173" 355367.0, \
 "tig00000174" 368855.0, \
 "" 3143732 \
)
set size 1,1
set grid
unset key
set border 5
set tics scale 0
set xlabel "REF"
set ylabel "%SIM"
set format "%.0f"
set mouse format "%.0f"
set mouse mouseformat "[%.0f, %.0f]"
if(GPVAL_VERSION < 5) { set mouse clipboardformat "[%.0f, %.0f]" }
set xrange [1:3143732]
set yrange [1:110]
set style line 1  lt 2 lw 3
set style line 2  lt 2 lw 3
set style line 3  lt 1 lw 3 pt 6 ps 1
plot \
 "1_assembly_mummerplot.fplot" title "FWD" w l ls 1, \
 "1_assembly_mummerplot.rplot" title "REV" w l ls 2, \
 "1_assembly_mummerplot.hplot" title "HLT" w lp ls 3
print "-- INTERACTIVE MODE --"
print "consult gnuplot docs for command list"
print "mouse 1: coords to clipboard"
print "mouse 2: mark on plot"
print "mouse 3: zoom box"
print "'h' for help in plot window"
print "enter to exit"
pause -1
