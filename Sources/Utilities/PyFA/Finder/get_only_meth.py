import re

#===============================================================================
# Function to get only names of methods without module names as "prefix"
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - meth_list:     list of methods without module names (e.g. Fetch_Profile)
# Used by:
#   - Function in browse.py for checking directories
#-------------------------------------------------------------------------------
def get_only_meth(file_name_with_path):
  module_name = get_mod(file_name_with_path)

  methods = []

  with open(file_name_with_path) as file:                # open file
    for line in file:                                    # read line by line
      meths = re.findall("(?<=_Mod/)(.*)(?=.f90)", line) # search _Mod and .f90
      methods.append(meths)                              # add lines to list
  methods2 = [x for x in methods if x != []]             # remove empty lists

  flat_meth_list = []                        # create only one list with strings
  for sublist in methods2:                   # instead of list with lists
    for item in sublist:
      flat_meth_list.append(item)

  # Check if any method is found
  if len(module_name) == 0:
    meth_list = [""]
  elif flat_meth_list == []:
    meth_list = ["No methods defined"]
  elif flat_meth_list != []:
    meth_list = [i.split()[0] for i in flat_meth_list]

  return meth_list

