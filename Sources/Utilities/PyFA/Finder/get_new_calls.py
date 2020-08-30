import re

#===============================================================================
#  Function for updating use statements
#
# Parameters:
#   - file_paths:    fortran files with full paths in front
#   - obj_list:      list of all objects
# Returns:
#   - obj_list:      list of all objects updated
# Used by:
#   - Main program
#===============================================================================
def get_new_calls(file_paths, obj_list, obj_memb):

  print("Searching for call statements ...")

  # Get all functions names from obj_list into a list
  fun_list_names = []
  for o in range(len(obj_list)):
    if obj_list[o].Type() == "Function":
      name = obj_list[o].name
      if "(" in name:
        name = name.split("(")
        fun_list_names.append(name[0])
      else:
        fun_list_names.append(name[0])

  # Get all subroutine names from obj_list into a list
  sub_list_names = []
  for o in range(len(obj_list)):
    if obj_list[o].Type() == "Subroutine":
      name = obj_list[o].name
      if "(" in name:
        name = name.split("(")
        sub_list_names.append(name[0])
      else:
        sub_list_names.append(name[0])

  # Put all subroutine and function names in list
  calling_names = [*fun_list_names, *sub_list_names]

  # Search through all non-memer objects
  for o in range(len( (obj_list) )):               # through objects
    obj_list[o].calls = []
    file_path = obj_list[o].path
    for c in range(len(calling_names)):            # c - calling counter
      if calling_names[c] not in file_path:
        with open(file_path) as file:              # open each file
          for line in file:
            line = line.split('!',          1)[0]  # remove comment
            line = line.split('function',   1)[0]  # remove lines defining ...
            line = line.split('subroutine', 1)[0]  # ... functions or subs
            if calling_names[c] in line and "_Mod_" not in line:
              obj_list[o].calls.append(calling_names[c])

  # Search through all member objects
  for o in range(len( (obj_memb) )):             # through objects
    mod = obj_memb[o].in_module
    file_path = obj_memb[o].path
    for c in range(len(calling_names)):          # c - calling counter
      with open(file_path) as file:              # open each file
        for line in file:
          line = line.split('!',          1)[0]  # remove comment
          line = line.split('function',   1)[0]  # remove lines defining ...
          line = line.split('subroutine', 1)[0]  # ... functions or subroutines
          if calling_names[c] in line and "_Mod_" not in line:
            mod.calls.append(calling_names[c])

  # Remove duplicate entries from the list
  for o in range(len( (obj_list) )):
    obj_list[o].calls = list(dict.fromkeys(obj_list[o].calls))

  return obj_list

