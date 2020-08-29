import Const
from Xfig.plot_object     import plot_object
from Xfig.plot_grid       import plot_grid
from Xfig.plot_legend     import plot_legend
from Xfig.plot_spline     import plot_spline
from Xfig.write_comment   import write_comment

#===============================================================================
# Plot everything (the entire graph) from object list
#
# Parameters:
#   - file:       Xfig file's handle
#   - obj_list:   list of all objects representing modules or subroutines
#   - box_margins box margins which can be changed by command line option
# Returns:
#   - nothing
#-------------------------------------------------------------------------------
def plot_all(file, obj_list, spl_list, box_margins):

  # Plot boxes
  write_comment(file, "Objects", 5)
  for i in range(len(obj_list)):
    plot_object(file, obj_list[i])

  # Plot all connections
  write_comment(file, "Connections", 5)
  for s in range(len(spl_list)):
    plot_spline(file, spl_list[s])

  # Plot grid
  plot_grid(file, obj_list)

  # Plot legend
  plot_legend(file, obj_list, 0, -8 * Const.UNIT_BOX_HEIGHT)

