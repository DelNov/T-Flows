#===============================================================================
# Removes path (can be a list of paths) from path list
#-------------------------------------------------------------------------------
def remove_path(file_paths,path):

  # Can be a list of paths or only 1 path <===============¬
  #                                                      ||
  new_list = [x for x in file_paths if not x.startswith(path)]

  return new_list

