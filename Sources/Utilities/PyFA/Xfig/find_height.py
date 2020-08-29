import Const

#===============================================================================
# Function for calculating height of a box
#
# Parameters:
#   - object:    object for calculating height
# Returns:
#   - height     height of the box
# Used by:
#   - Function for updating object attributes
#-------------------------------------------------------------------------------
def find_height(object):

  height = (  1                      \
            + object.H_Vars()        \
            + object.H_Methods()     \
            + object.N_Uses()        \
            + object.N_Types() ) * Const.UNIT_BOX_HEIGHT

  return height

