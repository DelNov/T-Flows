import Finder
from Objects.subroutine_class import subroutine_class
from Objects.object_level     import object_level

#===============================================================================
# Function for creating subroutines and appending into list
#
# Parameters:
#   - file_paths:     list of all paths to .f90 files
# Returns:
#   - sub_list:       list with only subrorutine objects
# Used by:
#   - Function for creating complete and updated object list
#-------------------------------------------------------------------------------
def subroutine_objects(file_paths, mod_list):

  subroutines_list = []

  for i in range(len(file_paths)):
    sub_name = Finder.get_sub(file_paths[i])  # find subroutines from file paths
    if sub_name != 0:                         # if it is sub then append to list
      subroutines_list.append(subroutine_class(file_paths[i]))

  sub_list = object_level(subroutines_list, mod_list)

  return sub_list

