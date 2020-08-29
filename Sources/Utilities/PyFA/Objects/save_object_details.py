#===============================================================================
# Function for saving logical coordinates into .ij file
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
def save_object_details(obj_list, file_name, object_hierarchy):

  # Write list of all names into a .txt file
  text_file = open(file_name,"w")

  text_file.write("#========================================================\n")
  text_file.write("#  name                                     detail\n")
  text_file.write("#--------------------------------------------------------\n")

  for o in range(len(obj_list)):

    # Write a separator between levels
    if o > 0:
      if object_hierarchy == "Column-Based":
        if obj_list[o].column > obj_list[o-1].column:
          text_file.write("#----------------------------" +  \
                          "----------------------------\n")
      if object_hierarchy == "Row-Based":
        if obj_list[o].row > obj_list[o-1].row:
          text_file.write("#----------------------------" +  \
                          "----------------------------\n")

    # Write coordinates
    if obj_list[o].vars_hidden and obj_list[o].methods_hidden:
      text_file.write("   %-40s %s\n" % (obj_list[o].name + ";",  \
                                         "Minimal"))
    elif obj_list[o].vars_hidden:
      text_file.write("   %-40s %s\n" % (obj_list[o].name + ";",  \
                                         "Reduced"))
    else:
      text_file.write("   %-40s %s\n" % (obj_list[o].name + ";",  \
                                         "Normal"))

  text_file.close()
  print("File", file_name, \
        "with object details has been created!")

