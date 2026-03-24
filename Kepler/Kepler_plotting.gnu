
set size ratio -1
set xrange [-11:11]
set yrange [-11:11]

do for [k=1:30000:100] {
    plot "Kepler_ENERGY_output.txt" every ::1::k using 6:7 with lines notitle, \
         "Kepler_ENERGY_output.txt" every ::k::k using 6:7 with points pt 7 ps 1.5 notitle
    pause 0.01
}