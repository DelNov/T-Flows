import Const

#===============================================================================
# Function to plot spline (with 2 coordinates)
#
# Parameters:
#   - file:       Xfig file's handle
#   - x0:         first coordinate on x axis
#   - y0:         first coordinate on y axis
# Returns:
#   - nothing
# Used by:
#   - function for plotting spline connections for legend
#-------------------------------------------------------------------------------
def plot_spline_legend(file, x0, y0, width, line_type):

  x1 = x0 + width

  if line_type == "Continuous":
    file.write("3 0 0 2 0 7 ")                     # 2 --> thickness
    file.write("%5d" % (50))                       # depth
    file.write(" -1 -1 0.000 0 1 1 2")             # 2 --> number of points
  else:
    file.write("3 0 1 2 0 7 ")                     # 2 --> thickness
    file.write("%5d" % (50))                       # depth
    file.write(" -1 -1 4.000 0 1 1 2")             # 2 --> number of points

  # Arrow settings (135 and 180 are width and length, divide by 15)
  if line_type == "Continuous":
    file.write("\n 1 1 2.00 135.00 180.00")        # arrow settings
    file.write("\n 6 1 2.00 135.00 180.00")        # arrow settings
  else:
    file.write("\n 1 0 2.00 135.00 180.00")        # arrow settings
    file.write("\n 6 0 2.00 135.00 180.00")        # arrow settings

  file.write("\n%9d %9d" % ( (x0) * Const.XFIG_SCALE,  \
                             (y0) * Const.XFIG_SCALE))
  file.write("%9d %9d" %   ( (x1) * Const.XFIG_SCALE,  \
                             (y0) * Const.XFIG_SCALE))
  file.write("\n 0.000 0.000\n")

