import Finder
from Objects.Module import Module

#===============================================================================
# Import attributes from fortran files to module object
#
# Parameters:
#   - file_path:    path to .f90 file
# Returns:
#   - module:       object of type "Module" with assigned attributes
# Used by:
#   - Function for appending modules (all module objects) into a list
#-------------------------------------------------------------------------------
def module_class(file_path):

  module_name  = Finder.get_mod(file_path)
  use_list     = Finder.get_use(file_path)
  var_list     = Finder.get_var(file_path)
  meth_list    = Finder.get_meth(file_path)
  call_list    = Finder.get_call(file_path)
  type_stat    = Finder.get_type(file_path)
  # Function type could come here, but it is set to "" in mother
  path         = file_path

  module = Module(module_name,  \
                  path,         \
                  use_list,     \
                  var_list,     \
                  meth_list,    \
                  call_list,    \
                  type_stat)

  return module

