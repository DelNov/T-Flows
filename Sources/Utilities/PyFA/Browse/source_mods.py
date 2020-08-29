#===============================================================================
# List of all fortran modules in root
#-------------------------------------------------------------------------------
def source_mods(root):

  source_mod = []

  for f_name in os.listdir(root):
    if f_name.endswith("_Mod.f90"):               # looks for all _Mod.f90 files
      source_mod.append(f_name)

  return source_mod

