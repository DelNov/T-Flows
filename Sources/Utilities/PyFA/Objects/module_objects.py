import Finder
from Objects.module_class     import module_class
from Objects.object_use_level import object_use_level

#===============================================================================
# Function for creating modules and appending into list
#
# Parameters:
#   - file_paths:     list of all paths to .f90 files
# Returns:
#   - mod_list:       list with only module objects
# Used by:
#   - Function for creating complete and updated object list
#-------------------------------------------------------------------------------
def module_objects(file_paths, mod_list):

  modules_list = []

  for i in range(len(file_paths)):
    module_name = Finder.get_mod(file_paths[i]) # find modules from file paths

    if module_name != []:                       # if it is module append to list
      modules_list.append(module_class(file_paths[i]))

  return modules_list

