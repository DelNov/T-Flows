from Xfig.find_width          import find_width
from Xfig.plot_type_frame     import plot_type_frame
from Xfig.plot_type_text_list import plot_type_text_list

#===============================================================================
# Function to plot type statements box with text
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - functions for plotting modules, subroutines and functions
#-------------------------------------------------------------------------------
def plot_type_list(file, object):

  box_width = find_width(object)

  # Plot type statement framing box first
  plot_type_frame(file, box_width, object)

  # Plot text
  plot_type_text_list(file, object)

