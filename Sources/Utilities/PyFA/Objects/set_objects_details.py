#===============================================================================
# Function for creating complete and updated object list
#
# Parameters:
#   - file_paths:   paths to .f90 files
# Returns:
#   - obj_list:     list with all created and updated objects
# Used by:
#   - main program (simple.py)
#-------------------------------------------------------------------------------
def set_objects_details(obj_list, object_details):

  for o in range(len(obj_list)):

    if object_details == "Normal":
      obj_list[o].vars_hidden    = False
      obj_list[o].methods_hidden = False

    elif object_details == "Reduced":
      obj_list[o].vars_hidden    = True
      obj_list[o].methods_hidden = False

    elif object_details == "Minimal":
      obj_list[o].vars_hidden    = True
      obj_list[o].methods_hidden = True

  return obj_list

