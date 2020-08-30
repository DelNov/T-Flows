import Const
from Spline.Spline import Spline
from Spline.walk   import walk

#===============================================================================
# Function to plot spline (with 6 coordinates)
#
# Parameters:
#   - file:     Xfig file's handle
#   - object1:  starting object (spline starts at the rigth side of this object)
#   - object2:  ending object   (spline ends at the left side of this object)
#   - depth:    depth of plotted spline
# Returns:
#   - nothing
# Used by:
#   - function for plotting spline connections
#-------------------------------------------------------------------------------
def create(obj_list, spl_list, object1, object2, line_type, depth, offset, stride):

  # print("Connecting ", object1.name, "and", object2.name)
  # print("offset =   ", offset)
  # print("stride =   ", stride)

  xc1 = object1.x0 + object1.w * 0.5
  xc2 = object2.x0 + object2.w * 0.5

  # 0.7 in lines below is to avoid coinciding lines
  if abs(xc1 - xc2) <= offset:
    x1 = object1.x0              # start at the lhs of object1
    x2 = x1 - offset * 0.7       # continue to the left
    x6 = object2.x0              # end on the lhs of object1
    x5 = x6 - offset * 0.7       # come from left side

  elif xc1 < xc2 - offset:
    x1 = object1.x0 + object1.w  # start at the rhs of object1
    x2 = x1 + offset * 0.7       # continue to the right
    x6 = object2.x0              # end on the lhs of object1
    x5 = x6 - offset * 0.7       # come from left side

  else:
    x1 = object1.x0              # start at the lhs of the object1
    x2 = x1 - offset * 0.7       # continue to the left
    x6 = object2.x0 + object2.w  # end on the rhs of object2
    x5 = x6 + offset * 0.7       # come from right side

  # First height depends on line_type
  if line_type == "Continuous":
    # y1 = object1.y0 + object1.h * 0.5  # starts in the middle of object1
    y1 = object1.y0  \
       + Const.UNIT_BOX_HEIGHT * 0.5  # starts from the middle of header
  elif line_type == "Dashed":
    y1 = object1.y0  \
       + Const.UNIT_BOX_HEIGHT * 0.5  # starts from the middle of header

  # Second coordinate should be the same as first
  y2 = y1

  # Last coordinate for continous lines (use statements)
  if line_type == "Continuous":
    ind = object2.uses.index("use " + object1.name)
    y6 = object2.y0 + Const.UNIT_BOX_HEIGHT        \
                    + object2.N_Types()            \
                    + ind * Const.UNIT_BOX_HEIGHT  \
                    + 0.5 * Const.UNIT_BOX_HEIGHT

  # Last coordinate for dashed lines (call statements)
  elif line_type == "Dashed":
    y6 = object2.y0   \
       + Const.UNIT_BOX_HEIGHT * 0.5  # hits in the middle of the header

  # Penultimate coordinate should be the same as last
  y5 = y6

  # Create a new spline object and walk!
  spline = Spline(object1.name, object2.name, line_type, depth)

  spline.x,  \
  spline.y = walk(x1, y1, x2, y2, x5, y5, x6, y6, obj_list, spl_list, stride)

  return spline
