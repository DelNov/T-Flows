from Xfig.plot_object_name           import plot_object_name
from Xfig.plot_use_list              import plot_use_list
from Xfig.plot_meth_list             import plot_meth_list
from Xfig.plot_var_list              import plot_var_list
from Xfig.plot_type_list             import plot_type_list
from Xfig.plot_object_end_compound   import plot_object_end_compound
from Xfig.plot_object_start_compound import plot_object_start_compound
from Xfig.write_comment              import write_comment

#===============================================================================
# Function to plot object box
#
# Parameters:
#   - file:            Xfig file's handle
#   - object:          object to plot
# Returns:
#   - nothing
# Used by:
#   - function for plotting module/subroutine/function/program
#-------------------------------------------------------------------------------
def plot_object(file, object):

  # Print a comment
  write_comment(file, object.name, 3)

  # Start a compound around the module
  plot_object_start_compound(file, object)

  # Plot a header text box
  plot_object_name(file, object)

  # Plot a type statements box
  if object.N_Types() > 0:
    plot_type_list(file, object)

  # If use statements have been found, plot use text box
  if object.N_Uses() > 0:
    plot_use_list(file, object)

  # If variables have been found, plot variables text box
  if object.H_Vars() > 0:
    plot_var_list(file, object)

  # If methods have been found, plot methods text box
  if object.H_Methods() > 0:
    plot_meth_list(file, object)

  # End the compound around the module
  plot_object_end_compound(file, object)

