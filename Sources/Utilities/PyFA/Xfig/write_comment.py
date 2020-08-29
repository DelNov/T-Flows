#===============================================================================
# Function to write a comment in Xfig file
#
# Parameters:
#   - file:            Xfig file's handle
#   - text:            text to be written
#   - sep:             single character separator
# Returns:
#   - nothing
# Used by:
#   - function for plotting module/subroutine/function/program
#-------------------------------------------------------------------------------
def write_comment(file, text, height):

  sep = "-"
  if height == 5: sep = "="

  # Print a comment
  line1 = "#" + sep * 3 + sep * len(text) + sep * 3 + "#"
  line2 = "#   " + " " * len(text) + "   #"
  line3 = "#   " +         text    + "   #"
  if height > 1: file.write("%s\n" % line1)
  if height > 3: file.write("%s\n" % line2)
  file.write("%s\n" % line3)
  if height > 3: file.write("%s\n" % line2)
  if height > 1: file.write("%s\n" % line1)

