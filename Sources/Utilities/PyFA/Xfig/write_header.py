#===============================================================================
# Function to write xfig header
#
# Parameters:
#   - file:  xfig file's handle
# Returns:
#   - nothing
#-------------------------------------------------------------------------------
def write_header(file, command_line):

  file.write("#FIG 3.2\n")
  file.write("# Produced by command: ")
  file.write(command_line)
  file.write("\n")
  file.write("Landscape\n")
  file.write("Center\n")
  file.write("Metric\n")
  file.write("B0\n")
  file.write("100.00\n")
  file.write("Single\n")
  file.write("-2\n")
  file.write("1200 2\n")

