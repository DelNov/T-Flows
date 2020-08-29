import Const
from Xfig.code_color import code_color

#===============================================================================
# Function to plot an empty module frame
#
# Parameters:
#   - file:            Xfig file's handle
#   - object:          object to plot title
#   - box_width:       box width in centimeters
#   - box_height:      box height in centimeters
# Returns:
#   - nothing
# Used by:
#   - function for plotting module name box (header box)
#-------------------------------------------------------------------------------
def plot_title_frame(file, object, box_width, box_height):

  file.write("2 2 0 ")
  file.write("%3d "     % Const.THICKNESS)
  file.write("0")
  if object.Type()   == "Subroutine":
    file.write("%3d " % code_color(Const.COLOR_HEADER_SUBROUTINE))
  elif object.Type() == "Program":
    file.write("%3d " % code_color(Const.COLOR_HEADER_PROGRAM))
  elif object.Type() == "Module":
    file.write("%3d " % code_color(Const.COLOR_HEADER_MODULE))
  elif object.Type() == "Function":
    file.write("%3d " % code_color(Const.COLOR_HEADER_FUNCTION))
  file.write("15 -1 20 0.000 0 0 -1 0 0 5\n")
  file.write("%9d %9d"  % ( object.x0            *Const.XFIG_SCALE,   \
                            object.y0            *Const.XFIG_SCALE))
  file.write("%9d %9d"  % ((object.x0+box_width) *Const.XFIG_SCALE,   \
                            object.y0            *Const.XFIG_SCALE))
  file.write("%9d %9d"  % ((object.x0+box_width) *Const.XFIG_SCALE,   \
                           (object.y0+box_height)*Const.XFIG_SCALE))
  file.write("%9d %9d"  % ( object.x0            *Const.XFIG_SCALE,   \
                           (object.y0+box_height)*Const.XFIG_SCALE))
  file.write("%9d %9d\n"% ( object.x0            *Const.XFIG_SCALE,   \
                            object.y0            *Const.XFIG_SCALE))

