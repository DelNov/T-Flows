This seems to be unused

#===============================================================================
# Function for finding maximum width of all objects
#
# Parameters:
#   - obj_list:    list of objects
# Returns:
#   - max_width:   maximum width of all objets (boxes)
# Used by:
#   - functions for creating and updating grid
#-------------------------------------------------------------------------------
def max_width(obj_list):

  widths_list = []

  for i in range(len(obj_list)):
    widths = obj_list[i].w
    widths_list.append(widths)

  max_width = max(widths_list)

  return max_width

