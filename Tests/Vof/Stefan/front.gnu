set gri
set xlabel "Time (s)"
set ylabel "Front position (m)"
plot "stefans_solution.dat" u 1:2 w l t "T-Flows","interface_position.dat" u 1:2 w l t "Analytical solution"
