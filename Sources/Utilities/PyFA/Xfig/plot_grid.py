import grid
import Const
from Xfig.plot_line     import plot_line
from Xfig.plot_text     import plot_text
from Xfig.write_comment import write_comment

#===============================================================================
# Function to plot grid with coordinates for ROW BASED
#
# Parameters:
#   - file:         Xfig file's handle
#   - obj_list:     list of all objects
# Returns:
#   - nothing
# Used by:
#   - function for plotting everything (the entire graph)
#-------------------------------------------------------------------------------
def plot_grid(file, obj_list):

  # Print a comment
  write_comment(file, "Grid", 5)

  min_v = 0
  max_v = max(grid.x)
  min_h = 0
  max_h = max(grid.y)

  # Plot grid

  # Plot vertical lines
  for i in range(len(grid.x)):
    plot_line(file, grid.x[i], min_h, grid.x[i], max_h, Const.COLOR_GRID, 500)

  # Plot horizontal lines
  for j in range(len(grid.y)):
    plot_line(file, min_v, grid.y[j], max_v, grid.y[j], Const.COLOR_GRID, 500)

  # Plot coordinates of each spot in grid
  for i in range(len(grid.x)-1):
    for j in range(len(grid.y)-1):
      plot_text(file, "Right",                     \
                grid.x[i+1]-0.5, grid.y[j]+0.5,  \
                "({}, {})".format(i, j),         \
                Const.FONT_NORMAL, Const.COLOR_GRID, 500)

