import re

#===============================================================================
# Function to search through .f90 file and returns module name
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - module_name:             name of the module, return [] if none
# Used by:
#   - functions for assigning attributes to objects ! Variables in Grid_Mod
#-------------------------------------------------------------------------------
def get_mod(file_name_with_path):

  module   = []                                        # initialize module list
  pattern  = re.compile(".+?(?=_Mod$)", re.IGNORECASE)
  pattern2 = re.compile("^\s*\S*use.*", re.IGNORECASE) # avoiding uses(for prog)

  with open (file_name_with_path, 'rt') as myfile:    # open file
    for line in myfile:                               # read line by line
      if pattern.search(line) != None:                # search for pattern
        if not line.startswith("!"):                  # skip line start with "!"
          if not line.split(maxsplit=1)[0] == "!":    # skip line start with "!"
            if not pattern2.search(line) != None:     # skip line with use
              if not "::" in line:
                module.append(( line.rstrip("\n")))   # add lines to list

  module = [s.strip() for s in module if s.strip()]   # remove whitespaces

  if len(module) != 0:                                # if module is not empty
    mod_string  = module[0]                           # take the first string
    module_name = re.sub("module ", "", mod_string)   # return module name

  elif len(module) == 0:
    module_name = []

  if "!" in module_name:
    module_name = []

  return module_name

