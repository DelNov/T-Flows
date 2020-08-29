import Const
from Xfig.find_height import find_height
from Xfig.find_width  import find_width

#===============================================================================
# Function to start definition of a compound in xfig
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot (module)
# Returns:
#   - nothing
# Used by:
#   - plot_module
#-------------------------------------------------------------------------------
def plot_object_start_compound(file, object):

  box_width  = find_width(object)
  box_height = find_height(object)

  file.write("6 %9d %9d %9d %9d\n" % (                                   \
     object.x0             * Const.XFIG_SCALE - Const.COMPOUND_MARGIN,   \
     object.y0             * Const.XFIG_SCALE - Const.COMPOUND_MARGIN,   \
    (object.x0+box_width)  * Const.XFIG_SCALE + Const.COMPOUND_MARGIN,   \
    (object.y0+box_height) * Const.XFIG_SCALE + Const.COMPOUND_MARGIN))

