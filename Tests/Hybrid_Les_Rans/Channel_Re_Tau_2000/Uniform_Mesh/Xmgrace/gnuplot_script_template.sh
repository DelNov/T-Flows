#!/bin/gnuplot
# this gnuplot template is used to plot data from .dat file produced by
# Save_Results.f90
# it is executed in ./Xmgrace directory
# variables DAT_FILE_WITH_RESULTS_MACRO and TEX_FILE_WITH_PLOT_MACRO 
# must be substituted

#-------------------------------------------------------------------------------
#-------- composition
#-------------------------------------------------------------------------------
A4w=29.7 #cm
A4h=21   #cm

# offsets
set lmargin 6
set rmargin 2
set tmargin 0.5
set bmargin 2.5

#-------------------------------------------------------------------------------
#-------- functions
#-------------------------------------------------------------------------------
# if 0 < x <= 14 then func = x, NaN otherwise
linear_law(x) = (x > 0 && x <= 14) ? x : 1/0
k = 0.41
log_law(x) = (x >= 14) ? 1/k * log(x) + 5.2 : 1/0
#-------------------------------------------------------------------------------
#-------- output settings
#-------------------------------------------------------------------------------
set terminal epslatex size A4w cm,A4h*2/3 cm standalone color colortext 20
set output 'TEX_FILE_WITH_PLOT_MACRO.tex'
set multiplot layout 1,2 columnsfirst
#-------------------------------------------------------------------------------
#-------- axes ranges and ticks
#-------------------------------------------------------------------------------
set border linewidth 2

set logscale x 10

set xtics nomirror
set ytics nomirror

set format x "$10^{%L}$"
#-------------------------------------------------------------------------------
#-------- line/point style
#-------------------------------------------------------------------------------
# --- blue
set style line 1 lc rgb '#0060ad' lt 1 lw 4 pt 7 ps 1.5
# --- red
set style line 2 lc rgb '#dd181f' lt 1 lw 4 pt 5 ps 1.5
# --- green
set style line 3 lc rgb '#1bb532' lt 1 lw 4 pt 9 ps 1.5
# --- black dashed
set style line 4 lc rgb '#000000' lt 2 lw 4 pt 1 ps 1.5
# --- black
set style line 5 lc rgb '#000000' lt 1 lw 4 pt 4 ps 1.5
# --- black
set style line 6 lc rgb '#ff8c00' lt 1 lw 4 pt 3 ps 1.5
# --- purple
set style line 7 lc rgb '#CC00CC' lt 1 lw 4 pt 4 ps 1.5
#-------------------------------------------------------------------------------
#-------- legeng style
#-------------------------------------------------------------------------------
set key box opaque
set key samplen 1.
set key spacing 1.2
set key box ls 5 lw 4
set key top left # box position
set key Left      # text in box alignment
#-------------------------------------------------------------------------------
#-------- main block
#-------------------------------------------------------------------------------
set xlabel '$y^+$' offset 0., 0.5
set xrange [1e-1:5e3]
set xtics offset 0,0.0 border 0,10,5e3 scale 5
set mxtics 5

# ----- U_plus
set ylabel "$U^+$" offset 0.5, 0.0
set yrange [1e-1:35]
set ytics offset 0,0.0 border 0,10,35 scale 5
set mytics 5

plot \
1/0 with l ls 1 t '$y^+$', \
(linear_law(x)) with l ls 1 not, \
1/0 with l ls 4 t '$\frac{1}{k} log(y^+) + 5.2$', \
(log_law(x)) with l ls 4 not, \
1/0 with lp ls 2 ps 5 t 'Current', \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($2) with lp ls 2 not

# ----- k_plus
unset logscale x
set format x "$%2.1t\\\cdot10^{%L}$"
set xlabel '$y^+$' offset 0., 0.5
set xrange [1e-1:2e3]
set xtics offset 0,0.0 border 0,500,1500 scale 5
set mxtics 5

set ylabel "$k^+$" offset 0.5, 0.0
set yrange [0:7]
set ytics offset 0,0.0 border 0,1,7 scale 5
set mytics 2

plot \
1/0 with lp ls 1 ps 5 t '$k_{res}$', \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($3) with lp ls 1 not, \
1/0 with lp ls 2 ps 5 t '$k_{tot}$', \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($5) with lp ls 2 not

# necessary line
unset multiplot