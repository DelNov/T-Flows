#===============================================================================
# List of all fortran files except modules in root
#-------------------------------------------------------------------------------
def source_subs(root):

  source_sub = []
  source_mod = source_mods(root)

  for f_name in os.listdir(root):
    if f_name.endswith(".f90"):                     # looks for all .f90 files
      source_sub.append(f_name)

  source_sub = list(set(source_sub) - set(source_mod))  # removes _Mod.f90

  return source_sub

