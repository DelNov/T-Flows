from Xfig.find_width         import find_width
from Xfig.plot_use_frame     import plot_use_frame
from Xfig.plot_use_text_list import plot_use_text_list

#===============================================================================
# Function to plot use statements box with text
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - functions for plotting modules, subroutines and functions
#-------------------------------------------------------------------------------
def plot_use_list(file, object):

  box_width = find_width(object)

  # Plot use statements framing box first
  plot_use_frame(file, box_width, object)

  # Plot text
  plot_use_text_list(file, object)

