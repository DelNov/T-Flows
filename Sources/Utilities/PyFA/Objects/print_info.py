This seems to be unused

#===============================================================================
# Function for printing out all object information
#
# Parameters:
#   - obj_list:   list of objects
# Returns:
#   - nothing
# Used by:
#   - main program, only for printing information
#-------------------------------------------------------------------------------
def print_info(obj_list):

  for i in range(len(obj_list)):
    print("\nName: ",            obj_list[i].name,        \
          "\nType: ",            obj_list[i].type,        \
          "\nModules used: ",    obj_list[i].use,         \
          "\nVariables: ",       obj_list[i].var,         \
          "\nMethods: ",         obj_list[i].meth,        \
          "\nCalls: ",           obj_list[i].call,        \
          "\nType statements: ", obj_list[i].type_stat,   \
          "\nLevel: ",           obj_list[i].level,       \
          "\nWidth: ",           obj_list[i].w,           \
          "\nHeight: ",          obj_list[i].h,           \
          "\nX0:",               obj_list[i].x0,          \
          "Y0:",                 obj_list[i].y0,          \
          "\nFile path:",        obj_list[i].path)

