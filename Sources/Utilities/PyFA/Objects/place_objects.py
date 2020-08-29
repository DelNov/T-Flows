from Objects.find_max_lvl         import find_max_lvl
from Objects.get_objects_at_level import get_objects_at_level
from Objects.evenize_list         import evenize_list

#===============================================================================
# Function for creating complete and updated object list
#
# Parameters:
#   - file_paths:   paths to .f90 files
# Returns:
#   - obj_list:     list with all created and updated objects
# Used by:
#   - main program (simple.py)
#-------------------------------------------------------------------------------
def place_objects(obj_list, object_hierarchy):

  max_lvl = find_max_lvl(obj_list)   # get max level

  levels_objects = []
  obj_ordered    = []

  # List of lists of levels
  max_obj = 0  # max number of object over all levels
  for l in range(max_lvl + 1):
    levels_objects.append( get_objects_at_level(obj_list, l) )
    max_obj = max(max_obj, len(levels_objects[-1]))

  # Assign values to coordinates
  for l in range(len(levels_objects)):
    objs_at_level = levels_objects[l]

    obj_pos = evenize_list(0, max_obj, len(objs_at_level))

    # Browse through all objects at level "l"
    for o in range(len(objs_at_level)):
      if object_hierarchy == "Column-Based-Decreasing":
        objs_at_level[o].column = max_lvl - l # decreasing level
        objs_at_level[o].row    = obj_pos[o]  # object

      elif object_hierarchy == "Row-Based-Decreasing":
        objs_at_level[o].row    = max_lvl - l # decreasing level
        objs_at_level[o].column = obj_pos[o]  # object

      if object_hierarchy == "Column-Based-Increasing":
        objs_at_level[o].column = l           # increasing level
        objs_at_level[o].row    = obj_pos[o]  # object

      elif object_hierarchy == "Row-Based-Increasing":
        objs_at_level[o].row    = l           # increasing level
        objs_at_level[o].column = obj_pos[o]  # object


      obj_ordered.append(objs_at_level[o])

  return obj_ordered

