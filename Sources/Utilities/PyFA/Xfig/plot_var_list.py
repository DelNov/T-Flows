from Xfig.find_width         import find_width
from Xfig.plot_var_frame     import plot_var_frame
from Xfig.plot_var_text_list import plot_var_text_list

#===============================================================================
# Function to plot variable box with text
#
# Parameters:
#   - file:         Xfig file's handle
#   - object:       object to plot
# Returns:
#   - nothing
# Used by:
#   - functions for plotting modules, subroutines and functions
#-------------------------------------------------------------------------------
def plot_var_list(file, object):

  box_width = find_width(object)

  # Plot variable framing box first
  plot_var_frame(file, box_width, object)

  # Plot text
  plot_var_text_list(file, object)

