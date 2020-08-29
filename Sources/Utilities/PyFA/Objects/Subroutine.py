from Objects.Object import Object  # mother

#===============================================================================
# Defining Subroutine class
#
# Parameters:
#   - see explanation for Object class
#   - type_stat:  type statements of subroutine
#   - in_module:  holds the module name it is in
# Returns:
#   - nothing
# Used by:
#   - subroutine_class (object)
#-------------------------------------------------------------------------------
class Subroutine(Object):
  def __init__(self, name, path,              \
               use, var, meth, call, types):

    # Call mother's init
    Object.__init__(self, name, path,      \
                    use, var, meth, call, types, "")  # last is function type

    self.in_module = []

