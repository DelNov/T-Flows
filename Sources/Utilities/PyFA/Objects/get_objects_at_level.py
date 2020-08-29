#===============================================================================
# Function for creating lists of objects at specific level
#  and updating x coordinates
#
# Parameters:
#   - obj_list:     list of all objects
#   - level:        level
# Returns:
#   - list:         list of objects at level level with updated x0 coordinate
# Used by:
#   - Function for creating lists of classes at specific level
#-------------------------------------------------------------------------------
def get_objects_at_level(obj_list, level):

  list = []

  for i in range(len(obj_list)):
    if obj_list[i].level == level:
      list.append(obj_list[i])

  return list

