import re

#===============================================================================
# Function to search through .f90 file and returns program name
#
# Parameters:
#   - file_name_with_path:     fortran file with full path in front
# Returns:
#   - prog_name:               name of the function, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_prog(file_name_with_path):

  program = []                                 # initialize
  pattern    = re.compile("^\s*.program", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          program.append(( line.rstrip("\n")))     # add line with patt. to list
  program = [s.strip() for s in program if s.strip()] # remove whitespaces

  if len(program) != 0:                      # if program is not empty
    fun_string = program[0]                  # take the first string

    prog_name   = re.sub("program ", "", fun_string) # return program

    if prog_name.endswith("&"):
      prog_name = prog_name + ")"

    if prog_name.startswith("!"):
      prog_name = 0

  elif len(program) == 0:
    prog_name = 0

  return prog_name

