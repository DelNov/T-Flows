#===============================================================================
# List of all _Mod directories in root
#-------------------------------------------------------------------------------
def source_mod_dirs(root):

  source_mod_dir = []

  for f_name in os.listdir(root):
    if f_name.endswith("_Mod"):                # looks for all _Mod directories
      source_mod_dir.append(f_name)

  return source_mod_dir

