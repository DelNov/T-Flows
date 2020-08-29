import re

#===============================================================================
# Function to search through .f90 file and returns call statements
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - call_list:               list of call statements, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_call(file_name_with_path):

  call_name = []
  pattern   = re.compile("(call)\s", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          call_name.append(( line.rstrip("\n")))   # add line with patt. to list

  call_name = [s.strip() for s in call_name if s.strip()] # remove whitespace

  # If you only want to take name of call statement without "type" or "only"
  call_name_list = [i.split()[1] for i in call_name]           # take call name
  call_name_list = ([s.strip("(") for s in call_name_list])    # remove ","
  call_name_list = [i.rsplit("(")[0] for i in call_name_list]

  if call_name_list != []:                # call_name for whole line
    call_list = call_name_list            # call_name_list - take only name
    call_list = list(set(call_list))

  else:
    call_list = 0

  return call_list

