#===============================================================================
# Browse through directories and subdirectories
#-------------------------------------------------------------------------------
def check_directories(root):

  # Get module files
  mod_files     = sorted(source_mods(root))
  mod_files     = ["/home/simcic/Development/Synthetic-Eddies/"    \
                   + s for s in mod_files]
  mod_names     = []
  mod_dirs      = sorted(source_mod_dirs(root))
  mod_dirs      = ["/home/simcic/Development/Synthetic-Eddies/"    \
                  + s for s in mod_dirs]
  mod_file_name = sorted([re.sub(".f90$", '', i) for i in mod_files])

  # Check if module names are same as their subdirectory names
  if mod_file_name == mod_dirs:
    print("\nNo errors. Modules have their corresponding directories.\n")
    for i in range(len(mod_files)):
      mod_name = finder.get_mod(mod_files[i])
      mod_names.append(mod_name)
    mod_names = sorted(mod_names)

    print("Modules: ", mod_names)

  # Check if files in subdirectories match methods of their modules
    for d in range(len(mod_names)):
      print("\nEntering directory: %s" % mod_names[d])

      files = listdir(join(root, mod_names[d]))
      files = sorted([re.sub(".f90$", '', i) for i in files])
      print("Files in module directory: ",files)

      methods = sorted(finder.get_only_meth(mod_files[d]))
      print("Methods of module:         ",methods)

      if files == methods:
        print("\n", mod_names[d], "is looking fine.")

      else:
        print("\n Found ERROR in module directory: ",mod_names[d])

  # If module names are not the same as their subdirectory names
  else:
    print("Module files or directories missing!")

