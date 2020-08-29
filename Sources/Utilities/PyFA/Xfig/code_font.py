#===============================================================================
# Return the code value of a Xfig font
#
# Parameters:
#   - name:  font name, the same name as in Xfig
# Returns:
#   - number corresponding to font code, as defined in Xfig format
# Used by:
#   - functions which plot text
#-------------------------------------------------------------------------------
def code_font(name):

  if name   == "AvantGarde-Book":
    return  4
  elif name == "AvantGarde-Demi":
    return  6
  elif name == "Courier-Bold":
    return 14
  if name   == "Courier":
    return 12
  elif name == "Courier-Bold":
    return 14
  elif name == "Helvetica":
    return 16
  elif name == "Helvetica-Bold":
    return 18
  elif name == "Helvetica-Narrow":
    return 20
  elif name == "Helvetica-Narrow-Bold":
    return 22

