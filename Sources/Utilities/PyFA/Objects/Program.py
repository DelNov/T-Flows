from Objects.Object import Object  # mother

#===============================================================================
# Defining Program class
#
# Parameters:
#   - see explanation for Object class
#   - type_stat:  type statements of subroutine
# Returns:
#   - nothing
# Used by:
#   - program_class (object)
#-------------------------------------------------------------------------------
class Program(Object):
  def __init__(self, name, path,              \
               use, var, meth, call, types):

    # Call mother's init
    Object.__init__(self, name, path,             \
                    use, var, meth, call, types, "")  # last is function type

