import Const
from Xfig.find_width import find_width
from Xfig.plot_text  import plot_text

#===============================================================================
# Function to plot use statements (text)
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting use statements box
#-------------------------------------------------------------------------------
def plot_use_text_list(file, object):

  for i in range(object.N_Uses()):
    plot_text(file, "Left",                                          \
              object.x0 + Const.UNIT_BOX_HEIGHT*0.333,               \
              object.y0 + Const.UNIT_BOX_HEIGHT*object.N_Types()     \
                        + Const.UNIT_BOX_HEIGHT*(i+1)                \
                        + Const.UNIT_BOX_HEIGHT*0.75,                \
              object.uses[i], Const.FONT_NORMAL, "Black", 10)

