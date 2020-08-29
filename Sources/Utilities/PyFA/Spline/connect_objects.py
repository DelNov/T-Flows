import Const
from Spline.Create   import Create
from Spline.Compress import Compress

#===============================================================================
# Function for plotting all spline connections
#
# Parameters:
#   - obj_list:  list of all objects representing modules or subroutines
# Returns:
#   - nothing
# Used by:
#   - function for plotting everything (the entire graph) from object list
#-------------------------------------------------------------------------------
def connect_objects(obj_list, offset, stride):

  use_objects  = []
  mod_objects  = []
  call_objects = []

  # Getting list with only modules
  for i in range(len(obj_list)):
    if obj_list[i].Type() == "Module":
      mod_objects.append(obj_list[i])

  # Getting list with objects that have use statements
  for i in range(len(obj_list)):
    if obj_list[i].uses != []:
      use_objects.append(obj_list[i])

  # Getting list with objects that have call statements
  for i in range(len(obj_list)):
    if obj_list[i].call != 0:
      call_objects.append(obj_list[i])

  splines = []

  # Creating connections for use statements
  # (inc_, max_ = 2.0, 5.0 was OK)
  # (inc_, max_ = 3.0, 8.0 was still OK)
  counter = 0.0
  inc_cnt = 2.0
  max_cnt = 5.0
  for i in range(len(use_objects)):
    use = use_objects[i].uses
    for k in range(len(use)):
      used = use[k]
      used = used.strip("use ")
      for m in range(len(mod_objects)):
        if used == mod_objects[m].name:
          splines.append(                                                 \
            Create(obj_list, splines, mod_objects[m], use_objects[i],     \
                   "Continuous", 101+len(splines),                        \
                   offset * (1.0 + counter / max_cnt), stride))
          counter += inc_cnt
          if counter > max_cnt+0.5: counter -= (max_cnt)

  # Creating connections for call statements
  for i in range(len(call_objects)):
    call = call_objects[i].call
    for k in range(len(call)):
      called = call[k]
      for m in range(len(obj_list)):
        if called in obj_list[m].name:
          splines.append(                                                 \
            Create(obj_list, splines, call_objects[i], obj_list[m],       \
                   "Dashed", 201+len(splines),                            \
                   offset * (1.0 + counter / max_cnt), stride))
          counter += inc_cnt
          if counter > max_cnt+0.5: counter -= (max_cnt)

  for i in range(len(splines)):
    Compress(splines[i])

  return splines
