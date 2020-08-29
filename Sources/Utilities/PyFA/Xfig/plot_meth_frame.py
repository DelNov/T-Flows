import Const
from Xfig.code_color import code_color

#===============================================================================
# Function to plot an empty method box (frame without text)
#
# Parameters:
#   - file:            Xfig file's handle
#   - box_width:       box width in centimeters
#   - object:          object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting methods box
#-------------------------------------------------------------------------------
def plot_meth_frame(file, box_width, object):

  x0 =  object.x0              * Const.XFIG_SCALE
  x1 = (object.x0 + box_width) * Const.XFIG_SCALE
  y0 = (object.y0 + ( 1                             \
           + object.N_Types()                       \
           + object.N_Uses()                        \
           + object.H_Vars() ) * Const.UNIT_BOX_HEIGHT ) * Const.XFIG_SCALE
  y1 = y0 + object.H_Methods() * Const.UNIT_BOX_HEIGHT * Const.XFIG_SCALE

  file.write("2 2 0 ")
  file.write("%3d "       % Const.THICKNESS)
  file.write("0")
  file.write("%3d "       % code_color(Const.COLOR_BOX))
  file.write("13 -1 20 0.000 0 0 -1 0 0 5\n")
  file.write("%9d %9d"  % ( x0, y0) )
  file.write("%9d %9d"  % ( x1, y0) )
  file.write("%9d %9d"  % ( x1, y1) )
  file.write("%9d %9d"  % ( x0, y1) )
  file.write("%9d %9d\n"% ( x0, y0) )

