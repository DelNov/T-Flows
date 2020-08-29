import re

#===============================================================================
# Function to search through .f90 file and returns function name
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - fun_name:                name of the function, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_fun(file_name_with_path):

  function   = []
  pattern    = re.compile(".+?(?=function)", re.IGNORECASE)
  pattern2   = re.compile("^(.*[^&])\&$", re.IGNORECASE)
  pattern3   = re.compile(".+?(?=end)", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          if "!" in line:
            line = line.split("!")[0]
          if not "print" in line:
            function.append(( line.rstrip("\n")))    # add line with patt. to list

          if pattern2.search(line) != None:           # if "&" is found
            for line in myfile:
              if not pattern3.search(line) != None:   # if "end" is not found

                new_function = function
                new_function.append(( line.rstrip("\n"))),
                new_function = list(new_function)
                new_function = ''.join(new_function)
                new_function = new_function.replace(" ", "")
                new_function = re.sub("&", "", new_function)
                function[0]  = new_function

                if not pattern2.search(line) != None: # if "&" is not found
                  break                               # stop this inner for loop
                                                      # outer loop continues

  function = [s.strip() for s in function if s.strip()] # remove whitespaces

  if len(function) != 0:                      # if function is not empty
    fun_name = function[0]                    # take the first string

    if "function" in fun_name:
      fun_name = fun_name.split("function",1)[1]
      fun_name = fun_name.lstrip()
      fun_name = fun_name.rstrip()

    else:
      fun_name = 0

    if fun_name != 0:
      if fun_name.endswith("&"):
        fun_name = fun_name + ")"

      if fun_name.startswith("!"):
        fun_name = 0

  elif len(function) == 0:
    fun_name = 0

  return fun_name

