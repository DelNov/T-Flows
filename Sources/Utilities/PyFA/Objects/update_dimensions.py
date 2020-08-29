import Xfig

#===============================================================================
# Function for updating (importing) coordinates of objects
#
# Parameters:
#   - obj_list:     list of objects
# Returns:
#   - obj_list:     list of objects with updated dimensions
# Used by:
#   - Function for creating complete and updated file list
#-------------------------------------------------------------------------------
def update_dimensions(obj_list):

  for o in range(len(obj_list)):
    obj_list[o].w = Xfig.find_width(obj_list[o])
    obj_list[o].h = Xfig.find_height(obj_list[o])

  return obj_list

