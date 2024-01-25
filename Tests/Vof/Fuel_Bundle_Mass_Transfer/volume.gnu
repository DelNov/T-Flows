set term png
set output "volume.png"
set gri
set format y "%4.2E"
#set format y "%E"
set xlabel "Time (s)"
set ylabel "Vapor volume (m^3)"
plot "volume.out" u 1:3 w l t ""

