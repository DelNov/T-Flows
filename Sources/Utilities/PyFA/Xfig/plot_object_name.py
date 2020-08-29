import Const
from Xfig.find_width       import find_width
from Xfig.plot_title_frame import plot_title_frame
from Xfig.plot_text        import plot_text

#===============================================================================
# Function to plot object name box (function header box)
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot (function)
# Returns:
#   - nothing
# Used by:
#   - function for plotting function box
#-------------------------------------------------------------------------------
def plot_object_name(file, object):

  box_width = find_width(object)

  # Plot module framing box first
  plot_title_frame(file, object, box_width, Const.UNIT_BOX_HEIGHT)

  # Plot name
  full_name = " " + object.name
  if object.Type() == "Function":
    full_name = " " + object.fun_type + " " + object.name

  plot_text(file, "Center",                          \
            object.x0 + box_width*0.5,               \
            object.y0 + Const.UNIT_BOX_HEIGHT*0.75,  \
            full_name, Const.FONT_HEADER, "Black", 10)

