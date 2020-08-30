#===============================================================================
# Function to print help screen and exit the program
#
# Parameters:
#   - none
# Returns:
#   - nothing
# Used by:
#   - main program
#-------------------------------------------------------------------------------
def print_help_and_exit():
  print("\nUsage: pyfa.py [OPTIONS]")
  print("\nValid options are:\n")
  print("   -a [SWITCH]\
  Plot by specified object alignment: ")
  print("    \
   'straight' for straight alignment")
  print("    \
   'diagonal' for diagonal alignment")
  print("   -g [SWITCH]\
  Plot by specified object detail: ")
  print("    \
   'normal'  for normal representation")
  print("    \
   'reduced' for reduced representation")
  print("    \
   'minimal' for minimal representation")
  print("   -h [SWITCH]\
  Plot by specified object hierarchy: ")
  print("    \
   'row_dec'    for decreasing row based hierarchy")
  print("    \
   'column_dec' for decreasing column based hierarchy")
  print("    \
   'row_inc'    for increasing row based hierarchy")
  print("    \
   'column_inc' for increasing column based hierarchy")
  print("  -ij [FILE]  \
  Read (i,j) object coordinates from the file")
  print("   -m [MARGIN]\
  Set margin in cm for individual boxes")
  print("   -o [FILE]  \
  Read object details from the file. ")
  print("   -r [DIR]   \
  Root directory for browsing sources")
  print("   -s [FILE]  \
  Choose source list with paths for plotting")
  print("  -xy [FILE]  \
  Read (x,y) object coordinates from the file")
  print("\nExample1: pyfa.py -s source.list -a straight")
  print("Example2: pyfa.py -a straight\n")

  exit()
