set term pngcairo enhanced font "Arial,24" size 800,600
set output "volume-presen.png"
set gri
set xlabel "Time (s)"
set ylabel "Vapor volume (m^3)"
set format y "%.1t x 10^{%T}"
set xtic 0.01
set xran [0:0.04]
plot "volume.out" u 1:3 w l t ""

