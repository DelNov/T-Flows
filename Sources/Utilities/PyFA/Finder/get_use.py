import re

#===============================================================================
# Function to search through .f90 file and returns use statements
#
# Parameters:
#   - file_name_with_path:    fortran file with full path in front
# Returns:
#   - use_list:               list of use statements, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_use(file_name_with_path):

  use_name = []

  pattern    = re.compile("(use)\s", re.IGNORECASE)
  pattern2   = re.compile("!.*(use)\s", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          if not line.split(maxsplit=1)[0] == "!": # skip line start with "!"
            if not pattern2.search(line) != None:
              use_name.append(( line.rstrip("\n")))# add line with patt. to list

  use_name = [s.strip() for s in use_name if s.strip()] # remove whitespace

  # If you only want to take name of use statement without "type" or "only"
  use_name_list = [i.split()[1] for i in use_name]           # take use name
  use_name_list = ([s.strip(",") for s in use_name_list])    # remove ","
  use_name_list = ["use " + x for x in use_name_list]

  # Solve problem with having "!" in strings
  use_list = []
  for i in range(len(use_name_list)):
    string = use_name_list[i]
    string = string.split('!')[0]
    use_list.append(string)

  if use_list != []:
    use_list = use_list

  return use_list

