#!/bin/python2
# this python template is used to plot data from .dat file produced by
# Save_Results.f90
# it is executed in ./Xmgrace directory
# As argument you should supply a valid ../file.dat and 5 columns to read
# Launch as python2 python_matplotlib_script.py ../file.dat 1 2 3 4 5

# Header
import numpy              as np
import matplotlib.pyplot  as plt
import matplotlib.lines   as lines
import sys

# Update default matplotlib settings
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '24'
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.major.size'] = 24
plt.rcParams['xtick.minor.size'] = 16
plt.rcParams['ytick.major.size'] = 24
plt.rcParams['ytick.minor.size'] = 16
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['xtick.minor.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.0
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.figsize'] = [3508./300, 2480./300] # A4 at 300 dpi
plt.rcParams['image.cmap'] = 'jet'
plt.rcParams["legend.edgecolor"] = 'black'
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.borderaxespad"] = 0.
plt.rcParams["legend.labelspacing"] = 0.

#-Reference data
U_plus_ref_data   = 'u_plus_re_tau_590.dat'
k_plus_ref_data   = 'k_plus_re_tau_590.dat'
eps_plus_ref_data = 'eps_plus_re_tau_590.dat'
uv_plus_ref_data  = 'uv_plus_re_tau_590.dat'

y1, u1   = np.loadtxt(U_plus_ref_data,   unpack=True)
y2, k1   = np.loadtxt(k_plus_ref_data,   unpack=True)
y3, eps1 = np.loadtxt(eps_plus_ref_data, unpack=True)
y4, uv1  = np.loadtxt(uv_plus_ref_data,  unpack=True)

#-Input data
input_name = sys.argv[1]
y5_c   = int(sys.argv[2])-1
u2_c   = int(sys.argv[3])-1
k2_c   = int(sys.argv[4])-1
eps2_c = int(sys.argv[5])-1
uv2_c  = int(sys.argv[6])-1
y5, u2, k2, eps2, uv2 = np.loadtxt(input_name, \
  usecols=(y5_c,u2_c,k2_c,eps2_c,uv2_c), unpack=True)

# Assemble all data in one array ([x1, y1, x2, y2])
data = np.array([[y1, u1,   y5, u2  ], \
                 [y2, k1,   y5, k2  ], \
                 [y3, eps1, y5, eps2], \
                 [y4, uv1,  y5, uv2 ]] )

# Limits for each plot ([x_min, x_max, x_label_dx, y_min, y_max, y_label_dy])
axis = np.array([[0., 300., 50., 0., 30., 10.], \
                 [0., 300., 50., 0., 5. ,  2.], \
                 [0., 300., 50., 0., 0.3, 0.1], \
                 [0., 300., 50., 0., 1. , 0.5]] )

# Labels for each plot ([xlabel, ylabel])
label = np.array([[r'$y^+$', r'$U^+$'          ], \
                  [r'$y^+$', r'$k^+$'          ], \
                  [r'$y^+$', r'$\varepsilon^+$'], \
                  [r'$y^+$', r'$uv^+$'         ]] )

# Layout
fig, ax = plt.subplots(2, 2)
plt.subplots_adjust(left=0.05, right=0.97, top=0.97, bottom=0.05, \
                    hspace=0.15, wspace=0.15)

# Plot id [i,j] or [p]
p = -1
for i in range(2):
  for j in range(2):
    p = p + 1

    # Add 1d plots
    ax[i,j].plot(data[p,0], data[p,1], 'b-', \
                 data[p,2], data[p,3], 'r-', \
                 ms=30, markerfacecolor='none', markeredgewidth=3)
    ax[i,j].plot(data[p,0], data[p,1], 'b.', \
                 data[p,2], data[p,3], 'r.', \
                 ms=5, markeredgewidth=3)

    # Set axes
    ax[i,j].axis([axis[p,0], axis[p,1], axis[p,3], axis[p,4]])

    # Set scale x:y = 1 by tweaking axes -> no control over margins
    #ax[i,j].set_aspect('equal', adjustable='box')

    # Ticks and labels
    x_ticks = np.arange(axis[p,0], axis[p,1]+0.001,axis[p,2])
    y_ticks = np.arange(axis[p,3], axis[p,4]+0.001,axis[p,5])
    ax[i,j].xaxis.set_ticks(x_ticks)
    ax[i,j].yaxis.set_ticks(y_ticks)
    # Apply ticks on canvas
    fig.canvas.draw()
    # Replace last tick label
    xlabel = [k.get_text() for k in ax[i,j].get_xticklabels()]
    xlabel[-1] = label[p,0]
    ax[i,j].set_xticklabels(xlabel)
    ylabel = [k.get_text() for k in ax[i,j].get_yticklabels()]
    ylabel[-1] = label[p,1]
    ax[i,j].set_yticklabels(ylabel)

# Legend
leg1 = lines.Line2D([], [], color='b', marker='.',
                    ms=30, label=r'DNS $Re_\tau=590$')
leg2 = lines.Line2D([], [], color='r', marker='.',
                    ms=30, label=r'Present')
lines  = [leg1, leg2]
labels = [line.get_label() for line in lines]
plt.legend(lines, labels)

# Show
#plt.show()

# Save
plt.savefig(input_name.replace('.dat','.png').replace('../','./'))
