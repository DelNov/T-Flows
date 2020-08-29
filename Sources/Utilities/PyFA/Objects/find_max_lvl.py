#===============================================================================
# Function to find maximum level of objects from list
#
# Parameters:
#   - obj_list:    list of objects
# Returns:
#   - lvl:         max level of objects from list
# Used by:
#   - functions for creating grid and updating coordinates of objects
#-------------------------------------------------------------------------------
def find_max_lvl(obj_list):

  lvls = []

  for i in range(len(obj_list)):
    lvl = obj_list[i].level
    lvls.append(lvl)

  lvl = max(lvls)

  return lvl

