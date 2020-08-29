import Const
from Xfig.find_width import find_width
from Xfig.plot_text  import plot_text

#===============================================================================
# Function to plot variables (text)
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting variable box
#-------------------------------------------------------------------------------
def plot_var_text_list(file, object):

  if not object.vars_hidden:
    for i in range(object.H_Vars()):
      plot_text(file, "Left",                                         \
                object.x0 + Const.UNIT_BOX_HEIGHT*0.333,              \
                object.y0 + Const.UNIT_BOX_HEIGHT*object.N_Types()    \
                          + Const.UNIT_BOX_HEIGHT*object.N_Uses()     \
                          + Const.UNIT_BOX_HEIGHT*(i+1)               \
                          + Const.UNIT_BOX_HEIGHT*0.75,               \
                object.vars[i], Const.FONT_NORMAL, "Black", 10)

