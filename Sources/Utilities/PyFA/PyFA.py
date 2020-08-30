#!/usr/bin/python3

#===============================================================================
# Import libraries
#-------------------------------------------------------------------------------
import time
import os            # needed for getcwd
import sys           # needed to get command line arguments
import Const
import Xfig
import Finder
import Browse
import Objects
import Spline
from splash_screen       import splash_screen
from print_help_and_exit import print_help_and_exit

# Start measuring time
start = time.time()

# Initialize
file_paths = []
obj_list   = []

print_args       = "No"
object_hierarchy = "Row-Based-Decreasing"
object_details   = "Minimal"
box_margins      = Const.BOX_MARGINS

src_file     = "None"
r_specified  = "None"  # root directory
s_specified  = "None"  # list of sources
ij_specified = "None"
xy_specified = "None"
d_specified  = "None"

# If no command line arguments were specified, print help and exit
if len(sys.argv) == 1:
  splash_screen()
  print_help_and_exit()

#---------------------------------------
#
# Browse through command line arguments
#
#---------------------------------------
for j in range(1,len(sys.argv),2):

  if(len(sys.argv) > j+1):

    if str(sys.argv[j]) == "-a":

      if str(sys.argv[j+1]) == "yes":
        print_args  = "Yes"
      elif str(sys.argv[j+1]) == "no":
        print_args  = "No"

    # Check if object hierarchy was specified
    elif str(sys.argv[j]) == "-h":

      if str(sys.argv[j+1]) == "column_dec":
        object_hierarchy  = "Column-Based-Decreasing"
      elif str(sys.argv[j+1]) == "row_dec":
        object_hierarchy  = "Row-Based-Decreasing"
      elif str(sys.argv[j+1]) == "column_inc":
        object_hierarchy  = "Column-Based-Increasing"
      elif str(sys.argv[j+1]) == "row_inc":
        object_hierarchy  = "Row-Based-Increasing"
      else:
        print("Incorrect switch:", sys.argv[j+1], "after", sys.argv[j])
        print("Allowed switches are 'row_dec', 'column_dec',",  \
                                   "'row_inc' or 'column_inc'")
        print("Exiting the program")
        sys.exit()

    # Check if object representation was specified
    elif str(sys.argv[j]) == "-g":

      if str(sys.argv[j+1]) == "normal":
        object_details = "Normal"
      elif str(sys.argv[j+1]) == "reduced":
        object_details = "Reduced"
      elif str(sys.argv[j+1]) == "minimal":
        object_details = "Minimal"
      else:
        print("Incorrect switch:", sys.argv[j+1], "after", sys.argv[j])
        print("Allowed switches are: 'normal', 'reduced' or 'minimal'")
        print("Exiting the program")
        sys.exit()

    # Check if margins were specified
    elif str(sys.argv[j]) == "-m":

      box_margins = float(sys.argv[j+1])
      print("Box margins are set to:", str(sys.argv[j+1]))

    # Check if root for browsing was specified
    elif str(sys.argv[j]) == "-r":

      r_specified = sys.argv[j+1]
      print("Root directory for sources is:", str(sys.argv[j+1]))

    # Check if list of sources were specified
    elif str(sys.argv[j]) == "-o":

      d_specified = sys.argv[j+1]
      print("Object details are specified in:", str(sys.argv[j+1]))

    # Check if list of sources were specified
    elif str(sys.argv[j]) == "-s":

      s_specified = sys.argv[j+1]
      print("List of files is specified in:", str(sys.argv[j+1]))

    # Check if file with object (i,j) coordinates was specified
    elif str(sys.argv[j]) == "-ij":

      ij_specified = sys.argv[j+1]
      print("Object (i,j) coordinates are specified in:", str(sys.argv[j+1]))

    # Check if file with object (x,y) coordinates was specified
    elif str(sys.argv[j]) == "-xy":

      xy_specified = sys.argv[j+1]
      print("Object (x,y) coordinates are specified in:", str(sys.argv[j+1]))

    else:
      print("Unknow option:", sys.argv[j])
      print("Exiting the program")
      sys.exit()

#----------------
#
# Main algorithm
#
#----------------

#---------------------------------------------------------------
# List of sources was not specified, browse the whole directory
#---------------------------------------------------------------
if r_specified != "None":

  print("Analyzing Fortan sources in: " + r_specified)
  file_paths = Browse.source_paths(r_specified)

#---------------------------------------------
# List of sources was specified, read from it
#---------------------------------------------
else:

  # Get list of objects from source.list
  try: src_file = open(s_specified, 'rt')
  except:
    print("File", s_specified, "can't be found!  Exiting")
    sys.exit()

  with src_file:
    for line in src_file:
      if not line.startswith("#"):
        file_paths.append(line.rstrip("\n"))

  file_paths = list(filter(None, file_paths))

  for i in range(len(file_paths)):
    file_paths[i] = os.getcwd() + "/" + file_paths[i]

#------------------------------------------
# Take object list from file paths and ...
# ... find use level for each object
#------------------------------------------
obj_list, mod_list, obj_memb = Objects.get_obj_lists(file_paths, print_args)
obj_list = Objects.object_use_level(obj_list, mod_list)

#------------------------------------
# Find calls between objects and ...
# ... correct call levels for them
#------------------------------------
obj_list = Finder.get_new_calls(file_paths, obj_list, obj_memb)
obj_list = Objects.object_call_level(obj_list, obj_list)

if d_specified == "None":
  obj_list = Objects.set_objects_details(obj_list, object_details)
else:
  obj_list = Objects.load_object_details(d_specified, obj_list)

#--------------------------------------------------
# The following line updates object.h and object.w
# and it depends on object detail we want to have
#--------------------------------------------------
obj_list = Objects.update_dimensions(obj_list)

#----------------------------------------------
# Initial object placement, based on hierarchy
#----------------------------------------------
obj_list = Objects.place_objects(obj_list, object_hierarchy)

#-------------------------------------------------
# If logical coordinates specified, load them now
# (These should over-write those specified above)
#-------------------------------------------------
if ij_specified != "None":
  obj_list = Objects.load_ij_coordinates(ij_specified,     \
                                         obj_list,         \
                                         object_hierarchy)

#------------------------------------------
#
# Obviously the main function for plotting
#
#------------------------------------------

# Try to restore the full command line
command_line = '';
for arg in sys.argv:
  if ' ' in arg:
    command_line += '"{}" '.format(arg);
  else:
    command_line+="{} ".format(arg);

#---------------------------------------
# Find object coordinates in Xfig units
#---------------------------------------
Xfig.find_coordinates(obj_list, box_margins)
if xy_specified != "None":
  obj_list = Objects.load_xy_coordinates(xy_specified,      \
                                         obj_list,          \
                                         object_hierarchy)

#------------------------------------
# Create connections between objects
#------------------------------------          offset          stride
spl_list = Spline.connect_objects(obj_list, box_margins, box_margins * 0.5)

#----------------
# Open Xfig file
#----------------
file = open(Const.FIG_FILE_NAME, "w")

#------------------
# Write header out
#------------------
Xfig.write_header(file, command_line)

#-------------------------------------------
# Plot all fortran files starting from root
#-------------------------------------------
Xfig.plot_all(file, obj_list, spl_list, box_margins)

#-------------------------------------------------------------------------
# If neither (i,j) or (x,y) coordinates were not specified, save them now
#-------------------------------------------------------------------------
if ij_specified == "None" and xy_specified == "None":
  Objects.save_ij_coordinates(obj_list, Const.IJ_FILE_NAME, object_hierarchy)
  Objects.save_xy_coordinates(obj_list, Const.XY_FILE_NAME, object_hierarchy)

#---------------------------------------------------------
# If only (x,y) coordinates were not specified, save them
#---------------------------------------------------------
elif xy_specified == "None":
  Objects.save_xy_coordinates(obj_list, Const.XY_FILE_NAME, object_hierarchy)

if d_specified == "None":
  Objects.save_object_details(obj_list, Const.D_FILE_NAME, object_hierarchy)

# End
file.close()

# Print out execution time
end = time.time()
print("File", Const.FIG_FILE_NAME, "has been created!")
print("Execution time:", end - start, "seconds")
