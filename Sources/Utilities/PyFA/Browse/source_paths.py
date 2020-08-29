import os

#===============================================================================
# List of all fortran files in root and subdirectories
#-------------------------------------------------------------------------------
def source_paths(root):

  fortran_files = []
  for root, dirs, files in os.walk(root):
    for file in files:
      if file.endswith(".f90"):
        paths = os.path.join(root, file)
        fortran_files.append(paths)

  return fortran_files

