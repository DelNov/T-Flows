import re

#===============================================================================
# Function to search through .f90 file and returns type statements
#
# Parameters:
#   - file_name_with_path:   fortran file with full path in front
# Returns:
#   - type_list:             list of type statements, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_type(file_name_with_path):

  type_name = []
  pattern   = re.compile("^\s+(?=type\s+)", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          type_name.append(( line.rstrip("\n")))   # add line with patt. to list

  type_name = [s.strip() for s in type_name if s.strip()] # remove whitespace

  # If you only want to take name of type statement without "type" or "only"
  type_name_list = [i.split()[1] for i in type_name]           # take type name
  type_name_list = ([s.strip("(") for s in type_name_list])    # remove ","
  type_name_list = [i.rsplit("(",1)[0] for i in type_name_list]

  # Solve problem with having "!" in strings
  type_list = []
  for i in range(len(type_name)):
    string = type_name[i]
    string = string.split('!')[0]
    type_list.append(string)

  if type_list != []:
    type_list = type_list
    type_list = list(set(type_list))

  return type_list

