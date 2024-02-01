set term png
set output "courant.png"
set xlabel "Time step"
set ylabel "Courant number"
set gri
plot "courant.out" u 4 w l

