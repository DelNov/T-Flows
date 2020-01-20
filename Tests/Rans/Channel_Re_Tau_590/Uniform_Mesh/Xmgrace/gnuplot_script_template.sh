#!/bin/gnuplot
# this gnuplot template is used to plot data from .dat file produced by
# Save_Results.f90
# it is executed in ./Xmgrace directory
# variables DAT_FILE_WITH_RESULTS_MACRO and TEX_FILE_WITH_PLOT_MACRO 
# must be substituted

#-------------------------------------------------------------------------------
#-------- input data
#-------------------------------------------------------------------------------
U_PLUS   = sprintf('"u_plus_re_tau_590.dat"')
K_PLUS   = sprintf('"k_plus_re_tau_590.dat"')
EPS_PLUS = sprintf('"eps_plus_re_tau_590.dat"')
UV_PLUS  = sprintf('"uv_plus_re_tau_590.dat"')
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
#-------- output settings
#-------------------------------------------------------------------------------
set terminal epslatex size A4w cm,A4h*2/3 cm standalone color colortext 20
set output 'TEX_FILE_WITH_PLOT_MACRO.tex'
set multiplot layout 2,2 columnsfirst
#-------------------------------------------------------------------------------
#-------- axes ranges and ticks
#-------------------------------------------------------------------------------
set border linewidth 2

set xtics nomirror
set ytics nomirror

set format '$%g$'
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
set key spacing 1.1
set key box ls 5 lw 4
set key top right # box position
set key Left      # text in box alignment
#-------------------------------------------------------------------------------
#-------- main block
#-------------------------------------------------------------------------------
set xlabel '$y^+$' offset 0., 0.5
set xrange [0:300]
set xtics offset 0,0.0 border 0,50,290 scale 5
set mxtics 5

# ----- U_plus
set ylabel "$U^+$" offset 0.5, 0.0
set yrange [0:30]
set ytics offset 0,0.0 border 0,10,30 scale 5
set mytics 5

plot \
1/0 with p ls 1 ps 5 t 'DNS $Re_\tau=590$', \
@U_PLUS usi ($1):($2) with p ls 1 not, \
1/0 with lp ls 2 ps 5 t 'Current', \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($2) with lp ls 2 not

# remove legend
unset key

# ----- k_plus
set ylabel "$k^+$" offset 0.5, 0.0
set yrange [0:5]
set ytics offset 0,0.0 border 0,2,5 scale 5
set mytics 2

plot \
@K_PLUS usi ($1):($2) with p ls 1, \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($3) with lp ls 2

# ----- eps_plus
set ylabel "$\\varepsilon^+$" offset 0.5, 0.0
set yrange [0:0.3]
#set autoscale y
set ytics offset 0,0.0 border 0,0.1,0.3 scale 5
set mytics 2

plot \
@EPS_PLUS usi ($1):($2) with p ls 1, \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($4) with lp ls 2

# ----- uv_plus
set ylabel "$uv^+$" offset 0.5, 0.0
set yrange [0:1]
set ytics offset 0,0.0 border 0,0.5,1 scale 5
set mytics 5

plot \
@UV_PLUS usi ($1):($2) with p ls 1, \
'DAT_FILE_WITH_RESULTS_MACRO' usi ($1):($5) with lp ls 2

# necessary line
unset multiplot
