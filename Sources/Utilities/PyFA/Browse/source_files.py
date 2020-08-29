#===============================================================================
# List of all fortran files in root
#-------------------------------------------------------------------------------
def source_files(root):

  fortran_files = []
  for f_name in os.listdir(root):                 # looks for all .f90 files
    if f_name.endswith(".f90"):
      fortran_files.append(f_name)

  return fortran_files

