set term png
set output "massTransfer.png"
set gri
set xlabel "Time (s)"
set ylabel "Mass transfer (kg/s)"
plot "massTransfer.out" u 1:2 w l t "Vaporization","" u 1:3 w l t "Condensation"

