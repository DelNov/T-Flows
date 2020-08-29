import Finder
from Objects.program_class    import program_class
from Objects.object_use_level import object_use_level

#===============================================================================
# Function for creating programs and appending into list
#
# Parameters:
#   - file_paths:     list of all paths to .f90 files
# Returns:
#   - program_list:   list with only program objects
# Used by:
#   - Function for creating complete and updated object list
#-------------------------------------------------------------------------------
def program_objects(file_paths, mod_list):

  program_list = []

  for i in range(len(file_paths)):
    program_name = Finder.get_prog(file_paths[i]) # find program from file paths
    if program_name != 0:                 # if it is program then append to list
      program_list.append(program_class(file_paths[i]))

  return program_list

