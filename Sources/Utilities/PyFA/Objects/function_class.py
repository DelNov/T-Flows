import Finder
from Objects.Function import Function

#===============================================================================
# Import attributes from fortran files to function object
#
# Parameters:
#   - file_path:    path to .f90 file
# Returns:
#   - function:     object of type "Function" with assigned attributes
# Used by:
#   - Function for appending functions (all function objects) into a list
#-------------------------------------------------------------------------------
def function_class(file_path):

  fun_name   = Finder.get_fun(file_path)
  use_list   = Finder.get_use(file_path)
  var_list   = Finder.get_var(file_path)
  meth_list  = []
  call_list  = Finder.get_call(file_path)
  type_stat  = Finder.get_type(file_path)
  fun_type   = Finder.get_fun_type(file_path)
  path       = file_path

  function = Function(fun_name,   \
                      path,       \
                      use_list,   \
                      var_list,   \
                      meth_list,  \
                      call_list,  \
                      type_stat,  \
                      fun_type)

  return function

