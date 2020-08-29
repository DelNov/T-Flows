from Objects.classify_objects     import classify_objects
from Objects.module_objects       import module_objects
from Objects.program_objects      import program_objects
from Objects.subroutine_objects   import subroutine_objects
from Objects.function_objects     import function_objects
from Objects.update_dimensions    import update_dimensions
from Objects.object_use_level     import object_use_level

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
def get_obj_lists(file_paths):

  mod_list  = module_objects    (file_paths, [])        # list of module
  sub_list  = subroutine_objects(file_paths, mod_list)  # list of subroutine
  fun_list  = function_objects  (file_paths, mod_list)  # list of functions
  prog_list = program_objects   (file_paths, mod_list)  # list of programs

  obj_list  = [*mod_list,   \
               *sub_list,   \
               *fun_list,   \
               *prog_list]               # list of all classes(mod+sub+fun+prog)

  obj_list, obj_memb  = classify_objects(obj_list)

  return obj_list, mod_list, obj_memb

