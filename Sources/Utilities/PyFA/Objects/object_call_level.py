#===============================================================================
# Function for determining and importing levels of modules (iterate 8 times)
#
# Parameters:
#   - object_list:    list of module objects
#
# Returns:
#   - object_list:    list of module objects with imported correct levels
# Used by:
#   - Function for appending module objects into a list
#-------------------------------------------------------------------------------
def object_call_level(object_list, called_list):

  n = 0
  while n < 12:
    n += 1
    for i in range(len(object_list)):
      if object_list[i].N_Calls() > 0:        # if there are call statements
        obj_call_list = object_list[i].calls  # get use list of objects

        for k in range(len(obj_call_list)):     # for every called object
          for z in range(len(called_list)):
            bare_name = called_list[z].name.split("(")[0]
            if bare_name == obj_call_list[k]: # if module matches use
              if object_list[i].level <= called_list[z].level:
                object_list[i].level = called_list[z].level + 1

  return object_list
