#===============================================================================
# Function to end a compound in xfig
#
# Parameters:
#   - file:         Xfig file's handle
#   - x0:           object position on x axis in centimeters
#   - y0:           object position on y axis in centimeters
#   - text:         module name)
#   - object:       object to plot (module)
# Returns:
#   - nothing
# Used by:
#   - plot_module
#-------------------------------------------------------------------------------
def plot_object_end_compound(file, object):
  file.write("-6\n")

