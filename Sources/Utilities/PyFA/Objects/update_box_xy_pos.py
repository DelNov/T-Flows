#===============================================================================
# Function for updating only 1 box by row and column (change placement in grid)
#
# Parameters:
#   - obj_list:     list of objects
# Returns:
#   - list of objects
# Used by:
#   - Finder.py - Function for searching coordinates in file and updating them
#-------------------------------------------------------------------------------
def update_box_xy_pos(obj_list, name, x0, y0):

  # Strip name of whitespaces and new lines
  strip_name = name.replace(" ", "").replace("\n", "")

  # Assign new coordinates
  for i in range(len(obj_list)):

    # Strip object name of whitespaces and new lines
    strip_obj_name = obj_list[i].name.replace(" ", "").replace("\n", "")

    if strip_obj_name == strip_name:
      obj_list[i].x0 = x0
      obj_list[i].y0 = y0

  return obj_list

