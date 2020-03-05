#!/bin/python2
# this python template is used to plot data from ??????.dat files
# produced by test_build.sh(Process accuracy test)
# it is executed automatically by test_build.sh script
#
# You can also launch it manually according to readme_python_matplotlib_script

# Header
import numpy                as np
import matplotlib.pyplot    as plt
import matplotlib.lines     as lines
import matplotlib.offsetbox as offsetbox
import sys
import glob

# Update default matplotlib settings
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.size'] = '48'
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.major.pad'] = 15 # distance from axis to Mticks label
plt.rcParams['ytick.major.pad'] = 15 # distance from axis to Mticks label
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
plt.rcParams['figure.figsize'] = [3508./300, 2480./300] # A4 at 300 dpi
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.dpi'] = 100
plt.rcParams['image.cmap'] = 'jet'
plt.rcParams["legend.frameon"] = True
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.edgecolor"]  = 'k'
plt.rcParams["legend.labelspacing"] = 0
plt.rcParams["legend.handlelength"] = 1.5
plt.rcParams["legend.handletextpad"] = 0.5
plt.rcParams["legend.columnspacing"] = 0.1
plt.rcParams["legend.borderpad"] = 0.1
plt.rcParams['lines.linewidth'] = 2.
plt.rcParams['lines.markeredgewidth'] = 2.0
plt.rcParams['lines.markersize'] = 15
plt.rcParams["legend.numpoints"] = 2

# arg1 -> y, arg2 -> H, returns (3.0*(y/(0.5*H))*(1.0-y/H)
y_profile = lambda arg1, arg2: np.array(3.0*(arg1/(0.5*arg2))*(1.0-arg1/arg2));

files = np.sort(glob.glob('??????.dat'))

# Layout
fig, ax = plt.subplots(1, 1)
plt.subplots_adjust(left=0.10, right=1.0, top=0.97, bottom=0.09, \
                    hspace=0.15, wspace=0.15)

x_data = np.array([], dtype=np.int64)
y_data = np.array([], dtype=np.double)
z_data = np.array([], dtype=np.double)

print '{0:<1}{1:>17}{2:>18}{3:>18}'.format('#', 'Ny', 'L2(u)', 'LInf(u)')
for fi in range(len(files)):

  #-Input data
  y, u = np.loadtxt(files[fi], unpack=True)

  N_y = np.shape(y)[0]
  l2_error = np.sqrt(np.sum((u - y_profile(y, 1.0))**2)/N_y)
  linf_error = np.max(np.abs(u - y_profile(y, 1.0)))
  print '{0: 18d}{1: 18.8E}{2: 18.8E}'.format(N_y, l2_error, linf_error)

  x_data = np.append(x_data, N_y)
  y_data = np.append(y_data, l2_error)
  z_data = np.append(z_data, linf_error)

fi = open('results.dat','w')
data_to_write = np.array([x_data, y_data, z_data])

np.savetxt(fi, data_to_write.T)
fi.close()

# Add 1d plots
first_order_y_data  = np.ones(len(files)) * 1.2 * max(y_data[0], z_data[0])
second_order_y_data = np.ones(len(files)) * 0.8 * min(y_data[0], z_data[0])
for i in range(1,len(files)):
  first_order_y_data[i]  = first_order_y_data[i-1]/2.
  second_order_y_data[i] = second_order_y_data[i-1]/4.

# 1st order precision line
ax.loglog(x_data, first_order_y_data, \
  linestyle='dashed', color='r')

# 2nd order precision line
ax.loglog(x_data, second_order_y_data, \
  linestyle='dashed', color='g')

ax.loglog(x_data, y_data, \
  linestyle='solid', color='k', marker='o', \
  markerfacecolor='none')

ax.loglog(x_data, z_data, \
  linestyle='solid', color='b', marker='s', \
  markerfacecolor='none')

# Legend
leg1 = lines.Line2D([], [], color='r', linestyle='dashed', \
                    label=r'1st order')
leg2 = lines.Line2D([], [], color='g', linestyle='dashed', \
                    label=r'2nd order')
leg3 = lines.Line2D([], [], color='k', marker='o', markerfacecolor='none', \
                    label=r'$L_2$-error')
leg4 = lines.Line2D([], [], color='b', marker='s', markerfacecolor='none', \
                    label=r'$L_\infty$-error')
lines  = [leg1, leg2, leg3, leg4]
labels = [line.get_label() for line in lines]
plt.legend(lines, labels)

# Plot x label in a box
atext = offsetbox.AnchoredText(r'$N_y$', \
  loc='lower right', pad=0.05, borderpad=0.1, frameon=True)
atext.patch.set_edgecolor('w')
ax.add_artist(atext)

#-Input data
output_name = ""
if len(sys.argv) >= 2:
  output_name = sys.argv[1] # str type

# Correct if no input was given
if output_name == "" or output_name == "PNG_FILE_WITH_RESULTS_MACRO":
  output_name = 'results'

plt.savefig(output_name + '.png', bbox_inches='tight', pad_inches=0)

# Show
#plt.show()