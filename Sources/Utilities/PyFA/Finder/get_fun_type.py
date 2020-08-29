import re

#===============================================================================
# Function to search through .f90 file and returns function type
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - fun_type:                type of the function, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_fun_type(file_name_with_path):

  function = []                                 # initialize
  pattern  = re.compile(".+?(?=function)", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          function.append(( line.rstrip("\n")))    # add line with patt. to list
  function = [s.strip() for s in function if s.strip()] # remove whitespaces

  if len(function) != 0:                      # if function is not empty
    fun_string = function[0]                  # take the first string

    fun_type = fun_string.split(" function")[0] # split by "function" and
                                                # take only part before phrase
    if fun_type.endswith("&"):
      fun_type = fun_type + ")"

    if fun_type.startswith("!"):
      fun_type = ""

  elif len(function) == 0:
    fun_type = ""

  return fun_type

