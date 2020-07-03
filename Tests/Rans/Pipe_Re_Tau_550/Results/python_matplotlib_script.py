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
R_data = 're_tau_550_ref_data.dat'

# r, Uz, r+ Uz+ uzur^+, k_rr, k_tt, k_zz, eps_rr, eps_tt, eps_zz
r1, uz1, rp1, uzp1, uzurp1, kp1_1, kp1_2, kp1_3, epsp1_1, epsp1_2, epsp1_3 = \
np.loadtxt(R_data, usecols=(0,3,1,2,7,4,5,6,9,10,11), comments='%', unpack=True)
kp1 = (kp1_1**2 + kp1_2**2 + kp1_3**2)/2
epsp1 = np.abs((epsp1_1 + epsp1_2 + epsp1_3)/2)

#-Input data
input_name_1 = sys.argv[1]      # 1
rp2_c      = int(sys.argv[2])-1 # 2
uzp2_c     = int(sys.argv[3])-1 # 3
kp2_c      = int(sys.argv[4])-1 # 4
epsp2_c    = int(sys.argv[5])-1 # 5
uvp2_c     = int(sys.argv[6])-1 # 6
input_name_2 = sys.argv[7]      # 7
r2_c       = int(sys.argv[8])-1 # 8
uz2_c      = int(sys.argv[9])-1 # 9

rp2, uzp2, kp2, epsp2, uvp2 = np.loadtxt(input_name_1, \
  usecols=(rp2_c, uzp2_c, kp2_c, epsp2_c, uvp2_c), unpack=True)
r2, uz2 = np.loadtxt(input_name_2, \
  usecols=(r2_c,uz2_c), unpack=True)

uvp2 = np.abs(uvp2)

# extract U bulk from res file
with open(input_name_1) as f:
  content = f.readlines()
content = [x.strip() for x in content]
Ub = float(content[0].split()[3]) #
uz2 = uz2/Ub

# Assemble all data in one array ([x1, y1, x2, y2])
data = np.array([[rp1, uzp1,  rp2, uzp2 ], \
                 [rp1, uzurp1,rp2, uvp2 ], \
                 [rp1, epsp1, rp2, epsp2], \
                 [rp1, kp1,   rp2, kp2  ], \
                 [r1,  uz1,   r2,  uz2  ], \
                 [r1,  uz1-10,   r2,  uz2-10  ]] )

# Limits for each plot ([x_min, x_max, x_label_dx, y_min, y_max, y_label_dy])
axis = np.array([[0., 100., 50.,  0., 25.,  10.], \
                 [0., 100., 50.,  0., 1.0 , 0.5], \
                 [0., 100., 50.,  0., 0.25, 0.1], \
                 [0., 100., 50.,  0., 5.00, 2.0], \
                 [0., 1.,   0.5,  0., 1.5 , 1. ], \
                 [0., 1.,   0.5,  0., 1.5 , 1. ]] )

# Labels for each plot ([xlabel, ylabel])
label = np.array([[r'$r^+$', r'$U_z^+$'        ], \
                  [r'$r^+$', r'$uv^+$'         ], \
                  [r'$r^+$', r'$\varepsilon^+$'], \
                  [r'$r^+$', r'$k^+$'          ], \
                  [r'$r$',   r'$U_z$'          ], \
                  [r'$$',   r'$$'              ]] )

# Layout
fig, ax = plt.subplots(2, 3)
plt.subplots_adjust(left=0.05, right=0.97, top=0.97, bottom=0.05, \
                    hspace=0.15, wspace=0.15)

# Plot id [i,j] or [p]
p = -1
for i in range(2):
  for j in range(3):
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
    if (p==5):
      plt.axis('off')

# Legend
leg1 = lines.Line2D([], [], color='b', marker='.',
                    ms=30, label=r'DNS $Re_\tau=550$')
leg2 = lines.Line2D([], [], color='r', marker='.',
                    ms=30, label=r'Present')
lines  = [leg1, leg2]
labels = [line.get_label() for line in lines]
plt.legend(lines, labels)

# Show
plt.show()

# Save
plt.savefig(input_name_1.replace('.dat','.png').replace('../','./'))
