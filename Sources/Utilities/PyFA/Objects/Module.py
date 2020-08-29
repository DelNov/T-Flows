from Objects.Object import Object  # mother

#===============================================================================
# Defining Module class
#
# Parameters:
#   - see explanation for Object class
#   - type_stat:  type statements of module
# Returns:
#   - nothing
# Used by:
#   - module_class (object)
#-------------------------------------------------------------------------------
class Module(Object):
  def __init__(self, name, path,              \
               use, var, meth, call, types):

    # Call mother's init
    Object.__init__(self, name, path,             \
                    use, var, meth, call, types, "")  # last is function  type

