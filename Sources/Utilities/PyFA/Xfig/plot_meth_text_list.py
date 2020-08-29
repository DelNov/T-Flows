import Const
from Xfig.plot_text import plot_text

#===============================================================================
# Function to plot methods (text)
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting methods box
#-------------------------------------------------------------------------------
def plot_meth_text_list(file, object):

  if not object.methods_hidden:
    for i in range(object.H_Methods()):
      plot_text(file, "Left",                                         \
                object.x0 + Const.UNIT_BOX_HEIGHT*0.333,              \
                object.y0 + Const.UNIT_BOX_HEIGHT*object.N_Types()    \
                          + Const.UNIT_BOX_HEIGHT*object.H_Vars()     \
                          + Const.UNIT_BOX_HEIGHT*object.N_Uses()     \
                          + Const.UNIT_BOX_HEIGHT*(i+1)               \
                          + Const.UNIT_BOX_HEIGHT*0.75,               \
                object.methods[i], Const.FONT_NORMAL, "Black", 10)

