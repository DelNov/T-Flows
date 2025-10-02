set term pngcairo enhanced font "Arial,24" size 800,600
set output "massTransfer-presen.png"
set gri
set xlabel "Time (s)"
set ylabel "Mass transfer (kg/s)"
set format y "%.1t x 10^{%T}"
set xtic 0.01
set yran [-3e-5:1e-5]
set xran [0:0.04]
set key bottom
plot "massTransfer.out" u 1:2 every ::::7300 w l t "Vaporization","" u 1:3 every ::::7300 w l t "Condensation"

