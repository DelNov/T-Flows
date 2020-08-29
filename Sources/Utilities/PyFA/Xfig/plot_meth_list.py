from Xfig.find_width          import find_width
from Xfig.plot_meth_frame     import plot_meth_frame
from Xfig.plot_meth_text_list import plot_meth_text_list

#===============================================================================
# Function to plot methods box with text
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - functions for plotting modules, subroutines and functions
#-------------------------------------------------------------------------------
def plot_meth_list(file, object):

  box_width = find_width(object)

  # Plot methods framing box first
  plot_meth_frame(file, box_width, object)

  # Plot text
  plot_meth_text_list(file, object)

