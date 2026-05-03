set terminal gif animate delay 1 size 700,700
set output "Kepler_orbit_VEL.gif"

set size ratio -1
set xrange [-11:11]
set yrange [-11:11]

do for [k=1:30000:100] {
    plot "Kepler_VEL_output.txt" every ::1::k using 6:7 with lines notitle, \
         "Kepler_VEL_output.txt" every ::k::k using 6:7 with points pt 7 ps 1.5 notitle
}

unset output