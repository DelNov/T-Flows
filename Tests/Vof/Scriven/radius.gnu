set term png
set output "radius.png"
beta = 2.57864
alpha = 0.677/(958.4*4216.0)
r0 = 5e-5
t0 = r0**2.0/(4.0*beta**2.0*alpha)
set xran [0:0.002]
set gri
set xlabel "Time (s)"
set ylabel "Radius (m)"
set title "dT = 0.75 K"
plot 2.0*beta*sqrt(alpha*x),"bench-data.dat" u ($1+t0):(3.0*$3/(4.0*pi))**(1.0/3.0) w l 

