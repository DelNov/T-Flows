import re
from Finder.get_all_var       import get_all_var
from Finder.get_sub           import get_sub

#===============================================================================
# Function to decide number of printing variables (to print only global)
#
# Parameters:
#   - file_name_with_path:        fortran file with full path in front
# Returns:
#   - sub_var_list:               list of global variables in subroutines
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_var(file_name_with_path):

  var_list = get_all_var(file_name_with_path)
  sub_name = get_sub(file_name_with_path)
  sub_var_list = []

  if sub_name != 0:                                # if it is subroutine

    result = re.search("\((.*)\)", sub_name)

    if result:
      sub_var_list = result.group(0)

    if isinstance(sub_var_list, list):
      sub_var_list = var_list
    else:
      sub_var_list = sub_var_list.split(",")       # split by "," and remove ","
      sub_var_list = var_list[0:len(sub_var_list)] # return only sub variables

  else:                                            # if it is not subroutine
    sub_var_list = var_list                        # return all variables

  return sub_var_list

