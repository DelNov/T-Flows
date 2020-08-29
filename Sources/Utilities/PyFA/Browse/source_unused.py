#===============================================================================
# Prints all unused files and directories in root
#-------------------------------------------------------------------------------
def source_unused(root):

  source_mod_dir = source_mod_dirs(root)
  source_mod     = source_mods(root)
  source_sub     = source_subs(root)


  directories       = [f for f in listdir(root) if isdir (join(root, f))]
  source_unused_dir = list(set(directories) - set(source_mod_dir))

  source_unused_file = []
  for f_name in os.listdir(root):
    source_unused_file.append(f_name)

  source_unused_file = list(set(source_unused_file)    \
                          - set(source_mod)            \
                          - set(source_mod_dir)        \
                          - set(source_sub)            \
                          - set(source_unused_dir))

  print("\nUnused directories:\n",sorted(source_unused_dir), \
        "\n\nUnused files:\n",sorted(source_unused_file),"\n")

  return source_unused_file

