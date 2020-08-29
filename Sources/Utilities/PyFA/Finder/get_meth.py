import re
from Finder.get_mod           import get_mod

#===============================================================================
# Function to search through .f90 file and returns all methods(module functions)
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - meth_list:               list of all methods(module functions)
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_meth(file_name_with_path):
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

  meth_list = []

  if flat_meth_list != []:
    meth_list = [i.split()[0] for i in flat_meth_list]
    mod = get_mod(file_name_with_path)                     # get module name
    meth_list = [mod + "_" + x for x in meth_list]         # add module name

  return meth_list

