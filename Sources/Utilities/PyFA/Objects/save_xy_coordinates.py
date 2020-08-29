#===============================================================================
# Function for saving real coordinates into .xy file
#
# Parameters:
#   - obj_list:     list of objects
#   - file_name:    name of saved .txt file
# Returns:
#   - nothing
# Note:
#   - creates a .txt file in PyFA folder
# Used by:
#   - main program (simple.py)
#-------------------------------------------------------------------------------
def save_xy_coordinates(obj_list, file_name, object_hierarchy):

  # Write list of all names into a .txt file
  text_file = open(file_name,"w")

  text_file.write("#========================================================\n")
  text_file.write("#  object_hierarchy\n")
  text_file.write("#--------------------------------------------------------\n")
  text_file.write("%s\n" % object_hierarchy)

  text_file.write("#========================================================\n")
  text_file.write("#       x;       y;  name\n")
  text_file.write("#--------------------------------------------------------\n")

  for o in range(len(obj_list)):

    # Write a separator between levels
    if o > 0:
      if object_hierarchy == "Column-Based":
        if obj_list[o].column > obj_list[o-1].column:
          text_file.write("#------------------------------------------------\n")
      if object_hierarchy == "Row-Based":
        if obj_list[o].row > obj_list[o-1].row:
          text_file.write("#------------------------------------------------\n")

    # Write coordinates
    text_file.write(" %8.3f;%8.3f;  %s\n" % (obj_list[o].x0,      \
                                             obj_list[o].y0,      \
                                             obj_list[o].name))

  text_file.close()
  print("File", file_name, \
        "with (x,y) coordinates has been created!")

