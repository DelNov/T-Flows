This seems to be unused

#===============================================================================
# Function for finding max height of all objects
#
# Parameters:
#   - obj_list:     list of objects
# Returns:
#   - max_height:   maximum height of all objets (boxes)
# Used by:
#   - functions for creating and updating grid
#-------------------------------------------------------------------------------
def max_height(obj_list):

  heights_list = []

  for i in range(len(obj_list)):
    heights = obj_list[i].y1 - obj_list[i].y0
    heights_list.append(heights)

  max_height = max(heights_list)

  return max_height

