import grid
from Xfig.find_width  import find_width
from Xfig.find_height import find_height

#===============================================================================
# Find grid and object coordinates
#
# Parameters:
#   - obj_list:     list of all objects
# Returns:
#   - nothing
#-------------------------------------------------------------------------------
def find_coordinates(obj_list, box_margins):

  # Grid coordinates
  n_row = 0
  n_col = 0
  for o in range(len(obj_list)):
    n_row = max(n_row, obj_list[o].row)
    n_col = max(n_col, obj_list[o].column)

  widths  = [0] * (n_col + 1)
  heights = [0] * (n_row + 1)
  grid.x  = [0] * (n_col + 2)
  grid.y  = [0] * (n_row + 2)

  for o in range(len(obj_list)):
    row = obj_list[o].row
    col = obj_list[o].column
    widths[col]  = max(widths [col], find_width (obj_list[o])  \
                 + box_margins * 2.0)
    heights[row] = max(heights[row], find_height(obj_list[o])  \
                 + box_margins * 2.0)

  for o in range(len(obj_list)):
    row = obj_list[o].row
    col = obj_list[o].column
    grid.x[col] = sum(widths [0:col])
    grid.y[row] = sum(heights[0:row])
  grid.x[n_col+1] = sum(widths)
  grid.y[n_row+1] = sum(heights)

  # Find object coordinates
  for o in range(len(obj_list)):
    row = obj_list[o].row
    col = obj_list[o].column
    xc = (grid.x[col] + grid.x[col+1]) * 0.5
    yc = (grid.y[row] + grid.y[row+1]) * 0.5
    obj_list[o].x0 = xc - obj_list[o].w * 0.5
    obj_list[o].y0 = yc - obj_list[o].h * 0.5
    obj_list[o].x1 = xc + obj_list[o].w * 0.5
    obj_list[o].y1 = yc + obj_list[o].h * 0.5

