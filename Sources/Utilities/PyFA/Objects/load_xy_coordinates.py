from Objects.update_box_xy_pos import update_box_xy_pos

#===============================================================================
# Function for searching coordinates in .ij file and update them
#
# Parameters:
#   - file_with_names:  file with names and coordinates
#   - obj_list:         list of objects
# Returns:
#   - obj_list:         list of objects with updated placements in grid
# Used by:
#   - main program (function for changing object placement in grid)
#===============================================================================
def load_xy_coordinates(file_with_names, obj_list, object_hierarchy):

  list = obj_list
  try: myfile = open(file_with_names, 'rt')
  except:
    print("File", file_with_names, "can't be found!  Exiting")
    sys.exit()

  with myfile:
    for line in myfile:
      if not line.startswith("#"):
        oh = line.strip("\n").strip(" ")
        if oh != object_hierarchy:
          print("Object hierarchy specified in file \"%s\" (%s) is not" %  \
                (file_with_names, oh))
          print("compatible with the one specified from command line (%s)" % \
                object_hierarchy)
          print("Exiting!")
          sys.exit(0)
        break

    for line in myfile:
      if not line.startswith("#"):
        data = line.split(";")
        obj_list = update_box_xy_pos(list,            \
                                     data[2],         \
                                     float(data[0]),  \
                                     float(data[1]))
  myfile.close()

  return obj_list

