#===============================================================================
# Function to check all directories and unused files
#-------------------------------------------------------------------------------
def check_all(root):

  # Check directories for errors
  check_directories(root)

  # Print all unused files and subdirectories
  source_unused(root)
