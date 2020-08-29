#===============================================================================
# Function for separating member methods from all other types of objects
#
# Parameters:
#   - obj_list:     list of objects
# Returns:
#   - obj_list:     list of objects without subroutine objects already printed
#                   in module methods (functions)
# Used by:
#   - Function for creating complete and updated file list
#-------------------------------------------------------------------------------
def classify_objects(obj_list):

  obj_used = []
  obj_memb = []  # list of member methods

  for o in range(len(obj_list)):
    if "_Mod_" in obj_list[o].name:
      obj_memb.append(obj_list[o])
    else:
      obj_used.append(obj_list[o])

  for om in range(len(obj_memb)):
    i = obj_memb[om].name.find("_Mod_")
    mod_name = obj_memb[om].name[0:i+4]
    for ou in range(len(obj_used)):
      if obj_used[ou].name == mod_name:
        obj_memb[om].in_module = obj_used[ou]

  return obj_used, obj_memb

