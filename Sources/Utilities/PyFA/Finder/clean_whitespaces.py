import re

#===============================================================================
# Function to delete spaces in all strings in a list
#
# Parameters:
#   - list_item:     list that needs to bo "cleaned" of spaces in strings
# Returns:
#   - list_item:     list without spaces in all strings
# Used by:
#   - Function "get_all_var"
#-------------------------------------------------------------------------------
def clean_whitespaces(list_item):
  if isinstance(list_item, list):
    for index in range(len(list_item)):
      if isinstance(list_item[index], list):
        list_item[index] = clean_whitespaces(list_item[index])
      if not isinstance(list_item[index], (int, tuple, float, list)):
        list_item[index] = list_item[index].strip()

  return list_item

