#!/bin/python
# this python template is used to plot data from .dat file produced by
# Save_Results.f90
# it is executed automatically by test_build.sh script
#
# You can also launch it manually according to readme_python_matplotlib_script

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
plt.rcParams['font.size'] = '22.325'
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.major.pad'] = 15  # distance from axis to Mticks label
plt.rcParams['ytick.major.pad'] = 15  # distance from axis to Mticks label
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
plt.rcParams['figure.figsize'] = [3508./300, 2480./300]  # A4 at 300 dpi
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 100
plt.rcParams['image.cmap'] = 'jet'
plt.rcParams['legend.frameon'] = True
plt.rcParams['legend.fancybox'] = False
plt.rcParams['legend.framealpha'] = 1
plt.rcParams['legend.edgecolor'] = 'k'
plt.rcParams['legend.labelspacing'] = 0
plt.rcParams['legend.handlelength'] = 1.5
plt.rcParams['legend.handletextpad'] = 0.5
plt.rcParams['legend.columnspacing'] = 0.1
plt.rcParams['legend.borderpad'] = 0.1
plt.rcParams['lines.linewidth'] = 2.
plt.rcParams['lines.markeredgewidth'] = 2.0
plt.rcParams['lines.markersize'] = 15
plt.rcParams['legend.numpoints'] = 2

#-Reference data
U_plus_ref_data   = 'u_plus_re_tau_2000.dat'
kin_plus_ref_data = 'kin_tol_plus_les.dat'

y1, u1 = np.loadtxt(U_plus_ref_data, unpack=True)
y2, k1 = np.loadtxt(kin_plus_ref_data, unpack=True)

#-Input data
input_name  = sys.argv[1] # str type
output_name = sys.argv[2] # str type
y3_c    = int(sys.argv[3]) - 1
u2_c    = int(sys.argv[4]) - 1
kres2_c = int(sys.argv[5]) - 1
kmod2_c = int(sys.argv[6]) - 1
y3, u2, kres2, kmod2 = np.loadtxt(input_name, \
  usecols=(y3_c,u2_c,kres2_c,kmod2_c), unpack=True)
k2 = kres2 + kmod2

# Functions
xlin = np.arange(start=1e-1, stop=14+0.5, step=0.5)
ylin = xlin

k = 0.41
xlog = np.arange(start=14, stop=2e3+0.5, step=0.5)
ylog = 1./k*np.log(xlog) + 5.2

# Assemble all data in one array ([x1, y1, x2, y2])
data = np.array([[y1,    u1,   y3,   u2  ], \
                 [y2,    k1,   y3,   k2  ]] )

dataf = np.array([[xlin, ylin, xlog, ylog]] )

# Limits for each plot ([x_min, x_max, x_label_dx, y_min, y_max, y_label_dy])
axis = np.array([[-1, 4, 6, 0., 40., 10.  ], \
                 [-1, 4, 6, 0., 7.5, 2.5  ]] )

# Labels for each plot ([xlabel, ylabel])
label = np.array([[r'$y^+$', r'$U^+$'          ], \
                  [r'$y^+$', r'$k^+$'          ]] )

# Layout
fig, ax = plt.subplots(1, 2)
plt.subplots_adjust(left=0.10, right=1.0, top=0.97, bottom=0.09, \
                    hspace=0.30, wspace=0.30)

# Plot id [i] or [p]
p = -1
for i in range(2):
  p = p + 1

  # Add 1d plots
  ax[i].semilogx(data[p,0], data[p,1], 'b-', \
                 data[p,2], data[p,3], 'r-'  )

  ax[i].semilogx(data[p,0], data[p,1], 'b.', \
                 data[p,2], data[p,3], 'r.', \
                 ms=5, markeredgewidth=3)
  # Functions
  if p == 0:
    ax[i].semilogx(dataf[p,0], dataf[p,1], 'k--', \
                   dataf[p,2], dataf[p,3], 'k--'  )

  # Set axes
  ax[i].axis([10.**axis[p,0], 10.**axis[p,1], axis[p,3], axis[p,4]])

  # Ticks and labels
  #x_ticks = np.arange(axis[p,0], axis[p,1]+0.001, axis[p,2])
  x_ticks = np.logspace(axis[p,0], axis[p,1], axis[p,2], endpoint=True)
  y_ticks = np.arange(axis[p,3], axis[p,4]+0.001, axis[p,5])
  ax[i].xaxis.set_ticks(x_ticks)
  ax[i].yaxis.set_ticks(y_ticks)
  # Apply ticks on canvas
  fig.canvas.draw()
  # Replace last tick label
  xlabel = [k.get_text() for k in ax[i].get_xticklabels()]
  xlabel[-1] = label[p,0]
  ax[i].set_xticklabels(xlabel)
  ylabel = [k.get_text() for k in ax[i].get_yticklabels()]
  ylabel[-1] = label[p,1]
  ax[i].set_yticklabels(ylabel)

# Legend
leg1 = lines.Line2D([], [], color='b', marker='.',
                    ms=10, label=r'DNS $Re_\tau=2000$')
leg2 = lines.Line2D([], [], color='r', marker='.',
                    ms=10, label=r'Present')
leg3 = lines.Line2D([], [], color='k', linestyle='dashed', \
                    label=r'$\frac{1}{k}log(x + 5.2)$')
leg4 = lines.Line2D([], [], color='k', linestyle='dashed', \
                    label=r'$x$')
lines  = [leg1, leg2, leg3, leg4]
labels = [line.get_label() for line in lines]
ax[0].legend(lines, labels)

# Show
# plt.show()

# Save
plt.savefig(output_name, bbox_inches='tight', pad_inches=0)
