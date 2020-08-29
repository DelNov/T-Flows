#===============================================================================
# Mother for differet programming constructs:
#
# Function, Module, Program and Subroutine
#
# Parameters:
#
#   - self:       its own self, reserved parameter in Python
#   - name:       name of the object
#   - path:       path to file where the object is defined
#   - use:        list of object's use statements
#   - var:        list of object's variables
#   - meth:       list of object's methods
#   - call:       list of object's call statements
# Returns:
#   - nothing, it is a class definition
# Used by:
#   - Function, Module, Program and Subroutine classes
#-------------------------------------------------------------------------------
class Object():

  def __init__(self, name, path,      \
               l_uses, l_vars, l_meth, call, types, f_type):

    self.name      = name
    self.path      = path

    # Logical properties
    self.uses     = l_uses  # list of use statements
    self.vars     = l_vars  # list of local variables
    self.methods  = l_meth  # list of member functions
    self.call     = call    # list of calls
    self.types    = types   # list of defined types
    self.fun_type = f_type  # function type (makes sense for Function only)

    self.level          = 0
    self.vars_hidden    = False
    self.methods_hidden = False

    # Geometrical properties
    self.x0        = 0.0
    self.y0        = 0.0
    self.row       = 0
    self.column    = 0
    self.w         = 0.0
    self.h         = 0.0

#===============================================================================
# Returns type of the object
#-------------------------------------------------------------------------------
  def Type(self):
    return self.__class__.__name__

#===============================================================================
# Returns number of type statements (new types)
#
# (It would be more elegant if empty list was still a list, like this: [])
#-------------------------------------------------------------------------------
  def N_Types(self):
    return len(self.types)

  def N_Uses(self):
    return len(self.uses)

  def H_Vars(self):
    if len(self.vars) == 0:
      return 0
    else:
      if self.vars_hidden:
        return 0.2
      else:
        return len(self.vars)

  def H_Methods(self):
    if len(self.methods) == 0:
      return 0
    else:
      if self.methods_hidden:
        return 0.2
      else:
        return len(self.methods)


