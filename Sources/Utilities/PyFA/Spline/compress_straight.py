import math
import Const

#===============================================================================
# Walk from one object to another, avoiding all objects in the graph
#
# Parameters:
#   - x1, y1, ... x6, y6:  coordinates the way Ivan introduced them
#   - obj_list:            list of all objects
# Returns:
#   - x, y:                coordinates with all steps from one object to another
#-------------------------------------------------------------------------------
def compress_straight(spline):

  x = spline.x
  y = spline.y
  keep = []
  for i in range(spline.N_Points()):
    keep.append(True)

  #------------------------------------------------
  #
  # Eliminate the points in-between straight lines
  #
  #------------------------------------------------

  # Mark points in between straight lines for deletion
  for i in range(1, len(x)-1):
    dx_p = x[i+1] - x[i]
    dy_p = y[i+1] - y[i]
    dx_m = x[i]   - x[i-1]
    dy_m = y[i]   - y[i-1]
    if abs(dx_p - dx_m) < 0.4 and abs(dy_p - dy_m) < 0.4:  # GHOST NUMBERS
      keep[i] = False

  # Yet, keep the points next to ones which are kept (to preserve curves)
  keep_2 = keep[:]
  for i in range(1, len(x)-1):
    if not keep[i]:
      if keep[i-1] or keep[i+1]:
        keep_2[i] = True

  # Make a compressed list of x and y coordinates
  x_c = []
  y_c = []
  for i in range(0, len(x)):
    if keep_2[i]:
      x_c.append(x[i])
      y_c.append(y[i])

  spline.x = x_c
  spline.y = y_c

