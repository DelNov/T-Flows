#===============================================================================
# Function for creating lists of objects at specific row
#  and updating x coordinates
#
# Parameters:
#   - obj_list:     list of objects
#   - row:          row
# Returns:
#   - list:         list of objects with updated x0 coordinate
# Used by:
#   - Function for creating lists of classes at specific row
#-------------------------------------------------------------------------------
def row_list(obj_list,row):

  list = []
  for i in range(len(obj_list)):
    if obj_list[i].row == row:
      list.append(obj_list[i])

  return list

