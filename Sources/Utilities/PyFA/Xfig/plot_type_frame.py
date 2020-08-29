import Const
from Xfig.code_color import code_color

#===============================================================================
# Function to plot an empty type statement box (frame without text)
#
# Parameters:
#   - file:            Xfig file's handle
#   - box_width:       box width in centimeters
#   - object:          object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting type statement box
#-------------------------------------------------------------------------------
def plot_type_frame(file, box_width, object):

  if object.Type() == "Module":
    color = Const.COLOR_HEADER_MODULE
  if object.Type() == "Subroutine":
    color = Const.COLOR_HEADER_SUBROUTINE
  if object.Type() == "Function":
    color = Const.COLOR_HEADER_FUNCTION
  if object.Type() == "Program":
    color = Const.COLOR_HEADER_PROGRAM

  x0 =  object.x0              * Const.XFIG_SCALE
  x1 = (object.x0 + box_width) * Const.XFIG_SCALE
  y0 = (object.y0 + ( 1 )      * Const.UNIT_BOX_HEIGHT ) * Const.XFIG_SCALE
  y1 = y0 + object.N_Types()   * Const.UNIT_BOX_HEIGHT * Const.XFIG_SCALE

  file.write("2 2 0 ")
  file.write("%3d "     % Const.THICKNESS)
  file.write("0")
  file.write("%3d "     % code_color(color))
  file.write("11 -1 30 0.000 0 0 -1 0 0 5\n")         # 30*5 = 150% intensity
  file.write("%9d %9d"  % ( x0, y0) )
  file.write("%9d %9d"  % ( x1, y0) )
  file.write("%9d %9d"  % ( x1, y1) )
  file.write("%9d %9d"  % ( x0, y1) )
  file.write("%9d %9d\n"% ( x0, y0) )

