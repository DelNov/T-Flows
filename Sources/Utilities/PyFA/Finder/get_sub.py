import re

#===============================================================================
# Function to search through .f90 file and returns subroutine name
#
# Parameters:
#   - file_name_with_path:  fortran file with full path in front
# Returns:
#   - sub_name:             name of the subroutine, return 0 if none
# Used by:
#   - functions for assigning attributes to objects
#-------------------------------------------------------------------------------
def get_sub(file_name_with_path):

  subroutine = []
  pattern    = re.compile(".+?(?=subroutine)", re.IGNORECASE)
  pattern2   = re.compile("(.*)[&]\s*$", re.IGNORECASE)
  pattern3   = re.compile(".+?(?=end)", re.IGNORECASE)

  with open (file_name_with_path, 'rt') as myfile: # open file
    for line in myfile:                            # read line by line
      if pattern.search(line) != None:             # search for pattern
        if not line.startswith("!"):               # skip line starting with "!"
          if "!" in line:
            line = line.split("!")[0]
          subroutine.append(( line.rstrip("\n")))  # add line with patt. to list

          if pattern2.search(line) != None:           # if "&" is found
            for line in myfile:
              if "!" in line:
                line = line.split("!")[0]

              if not pattern3.search(line) != None:   # if "end" is not found
                new_subroutine = subroutine[0]
                line = line.rstrip("\n")
                line = ''.join(line)
                new_subroutine = new_subroutine + line

                # Editing of string
                new_subroutine = new_subroutine.replace(" ", "")
                new_subroutine = re.sub("subroutine", "", new_subroutine)
                new_subroutine = re.sub("&", "", new_subroutine)
                subroutine[0]  = new_subroutine

                if not pattern2.search(line) != None: # if "&" is not found
                  break                               # stop this inner for loop
                                                      # outer loop continues

  subroutine = [s.strip() for s in subroutine if s.strip()] # remove whitespaces

  if len(subroutine) != 0:                      # if subroutine is not empty
    sub_string = subroutine[0]                  # take the first string
    sub_name   = re.sub("subroutine ", "", sub_string)     # return subroutine

    if sub_name.endswith("&"):
      sub_name = sub_name + ")"

    if "recursive" in sub_name:
      sub_name = re.sub("recursive","",sub_name)
      sub_name = sub_name.lstrip()
      sub_name = sub_name.rstrip()

  elif len(subroutine) == 0:
    sub_name = 0

  return sub_name

