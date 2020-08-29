#===============================================================================
# Function for determining and importing levels of modules (iterate 12 times)
#
# Parameters:
#   - object_list:    list of module objects
#
# Returns:
#   - object_list:    list of module objects with imported correct levels
# Used by:
#   - Function for appending module objects into a list
#-------------------------------------------------------------------------------
def object_use_level(object_list, modules_list):

  # If looking for its own self
  if modules_list == []: modules_list = object_list

  n = 0
  while n < 12:
    n += 1
    for i in range(len(object_list)):
      if object_list[i].N_Uses() > 0:        # if there are use statements
        obj_use_list = object_list[i].uses   # get use list of objects
        obj_use_list = [i.split()[1] for i in obj_use_list]    # only take name
        obj_use_list = ([s.strip(",") for s in obj_use_list])  # modules without
                                                               # other info
        obj_lvl      = []
        for k in range(len(obj_use_list)):        # for every use in object
          for z in range(len(modules_list)):
            if modules_list[z].name == obj_use_list[k]: # if module matches use
              lvl = modules_list[z].level
              obj_lvl.append(lvl)   # add level
        if obj_lvl == []:           # if obj_lvl is empty, level is 0
          obj_lvl = [0]
        else:
          obj_lvl = obj_lvl

        obj_lvl = max(obj_lvl)       # take the biggest used obj level from list
        object_list[i].level = obj_lvl + 1   # add 1 level to max level

  return object_list

