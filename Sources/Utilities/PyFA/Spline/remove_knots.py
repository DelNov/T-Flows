import math
import Const

#===============================================================================
# Eliminates knots and hairpins in splines
#
# Parameters:
#   - x1, y1, ... x6, y6:  coordinates the way Ivan introduced them
#   - obj_list:            list of all objects
# Returns:
#   - x, y:                coordinates with all steps from one object to another
#-------------------------------------------------------------------------------
def remove_knots(spline):

  x = spline.x
  y = spline.y
  keep = []
  for i in range(spline.N_Points()):
    keep.append(True)

  # Mark points in between straight lines for deletion
  for i in range(2, len(x)-2):
    for j in range(len(x)-3, i+1, -1):
      if abs(x[i] - x[j]) < 0.4 and abs(y[i] - y[j]) < 0.4:  # GHOST NUMBERS
        for k in range(i,j):
          keep[k] = False
        break

  # Make a compressed list of x and y coordinates
  x_c = []
  y_c = []
  for i in range(0, len(x)):
    if keep[i]:
      x_c.append(x[i])
      y_c.append(y[i])

  spline.x = x_c
  spline.y = y_c

