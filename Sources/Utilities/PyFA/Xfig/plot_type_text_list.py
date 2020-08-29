import Const
from Xfig.find_width import find_width
from Xfig.plot_text  import plot_text

#===============================================================================
# Function to plot type statements (text)
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting type statements box
#-------------------------------------------------------------------------------
def plot_type_text_list(file, object):

  if object.N_Types() != 0:
    for i in range(object.N_Types()):
      plot_text(file, "Left",                                    \
                object.x0 + Const.UNIT_BOX_HEIGHT*0.333,         \
                object.y0 + Const.UNIT_BOX_HEIGHT*(i+1)          \
                          + Const.UNIT_BOX_HEIGHT*0.75,          \
                object.types[i], Const.FONT_HEADER, "Black", 10)

