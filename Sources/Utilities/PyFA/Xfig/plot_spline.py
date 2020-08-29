import Const
from Spline.Spline      import Spline
from Xfig.write_comment import write_comment

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
def plot_spline(file, spline):

  # Print a comment
  write_comment(file, spline.object1 + " ---> " + spline.object2, 3)

  # Take aliases
  depth     = spline.depth
  line_type = spline.line_type
  x         = spline.x
  y         = spline.y

  # Start writing a spline
  if line_type == "Continuous":
    file.write("3 2 0 2 0 7 ")
    file.write("%5d" % (depth))
    file.write(" -1 -1 0.000 0 1 1 %6d" % len(x))
  elif line_type == "Dashed":
    file.write("3 2 1 2 0 7 ")
    file.write("%5d" % (depth))
    file.write(" -1 -1 8.000 0 1 1 %6d" % len(x))  # 8.000 is dash length

  # Arrow settings (135 and 180 are width and length, divide by 15)
  if line_type == "Continuous":
    file.write("\n 1 1 2.00 135.00 180.00")
    file.write("\n 6 1 2.00 135.00 180.00")
  elif line_type == "Dashed":
    file.write("\n 1 0 2.00 135.00 180.00")
    file.write("\n 6 0 2.00 135.00 180.00")

  cnt = 0
  for i in range(len(x)):
    if cnt % 4 == 0:
      file.write("\n       ")
    file.write(" %9d %9d" % ( x[i] * Const.XFIG_SCALE,  \
                              y[i] * Const.XFIG_SCALE))
    cnt = cnt + 1

  cnt = 0
  for i in range(len(x)):
    if cnt % 4 == 0:
      file.write("\n       ")
    if i == 0 or i == len(x)-1:
      file.write(" 0.000")
    else:
      file.write(" 1.000")
    cnt = cnt + 1

  file.write("\n")

