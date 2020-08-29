import Finder
from Objects.Program import Program

#===============================================================================
# Import attributes from fortran files to program object
#
# Parameters:
#   - file_path:    path to .f90 file
# Returns:
#   - program:      object of type "Program" with assigned attributes
# Used by:
#   - Function for appending programs (all program objects) into a list
#-------------------------------------------------------------------------------
def program_class(file_path):

  prog_name  = Finder.get_prog(file_path)
  use_list   = Finder.get_use(file_path)
  var_list   = Finder.get_var(file_path)
  meth_list  = []
  call_list  = Finder.get_call(file_path)
  type_stat  = Finder.get_type(file_path)
  # Function type could come here, but it is set to "" in mother
  path       = file_path

  program = Program(prog_name,   \
                    path,        \
                    use_list,    \
                    var_list,    \
                    meth_list,   \
                    call_list,   \
                    type_stat)
  return program

