set terminal pdfcairo enhanced color font ',16'
#set output "cutoff_-12.pdf"
set output "rk4_SHO.pdf"
set border linewidth 1
#set size square (1,1)

set title "RK4 Solution for Simple Harmonic Oscillator"
#set title "Magnification for M_{lens} = 10^{-12}M_{solar} r=0.001"
#set title "iic"
set xlabel r"x"
set ylabel "y(x)

#set xtics (1082.9, 1082.95, 1083, 1083.05, 1083.1)
#set logscale yx
#set yrange [0.94:1.005]
#set xtics
set yrange [-3:3]
set xrange [0:6.28]
#set style line 2 lt 1 lc rgb "dark-green" lw 4 pt 6 ps 1 

set key bottom left

set grid
set tics scale 1.5

plot "rk4_SHO.txt" u 1:2 w l lw 2 lc rgb "black" dt 1 title "position",\
"rk4_SHO.txt" u 1:3 w l lw 2 lc rgb "black" dt 2 title "velocity",\
"rk4_SHO.txt" u 1:4 w l lw 2 lc rgb "red" dt 3 title "analytic position"
