import Const

#===============================================================================
# Choose box width depending on longest string
#
# Parameters:
#   - object:   name of the Fortran file being read (.f90)
# Returns:
#   - var_width:  box width in Xfig drawing units
# Used by:
#   - functions which plot boxes
# Warning:
#   - Uses ghost parameter 0.4 to convert width from characters to Xfig units
#-------------------------------------------------------------------------------
def find_width(object):

  # Find variable with maximum length
  longest_variable_name = ""
  if object.H_Vars() > 0 and not object.vars_hidden:
    longest_variable_name = max(object.vars, key=len)

  # Find max length of methods
  longest_method_name = ""
  if object.H_Methods() > 0 and not object.methods_hidden:
    longest_method_name = max(object.methods, key=len)

  # Find longest type name
  longest_type_name = ""
  if object.N_Types() > 0: longest_type_name = max(object.types, key=len)

  # Find longest use statement
  longest_use_name = ""
  if object.N_Uses() > 0: longest_use_name = max(object.uses,  key=len)

  lengths = [len(object.name) + len(object.fun_type),  \
             len(longest_variable_name),               \
             len(longest_method_name),                 \
             len(longest_use_name),                    \
             len(longest_type_name)]

  box_width = max(lengths)
  box_width = box_width  \
            * Const.UNIT_BOX_HEIGHT * 0.5  # gives the best ratio for width

  return box_width

