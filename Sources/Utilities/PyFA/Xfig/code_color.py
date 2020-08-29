#===============================================================================
# Return the code value of a Xfig box color
#
# Parameters:
#   - name:  color name, the same name as in Xfig
# Returns:
#   - number corresponding to color code, as defined in Xfig format
# Used by:
#   - functions which plot frames
#-------------------------------------------------------------------------------
def code_color(name):

  if name   == "Black":
    return  0
  elif name == "Blue":
    return  1
  elif name == "Green":
    return  2
  elif name == "Cyan":
    return  3
  elif name == "Red":
    return  4
  elif name == "Magenta":
    return  5
  elif name == "Yellow":
    return  6
  elif name == "White":
    return  7
  elif name == "LtBlue":
    return 11
  elif name == "Green4":
    return 12
  elif name == "Green3":
    return 13
  elif name == "Green2":
    return 14
  elif name == "Cyan4":
    return 15
  elif name == "Cyan3":
    return 16
  elif name == "Cyan2":
    return 17
  elif name == "Red4":
    return 18
  elif name == "Red3":
    return 19
  elif name == "Red2":
    return 20
  elif name == "Magenta4":
    return 21
  elif name == "Magenta3":
    return 22
  elif name == "Magenta2":
    return 23
  elif name == "Brown4":
    return 24
  elif name == "Brown3":
    return 25
  elif name == "Brown2":
    return 26
  elif name == "Pink4":
    return 27
  elif name == "Pink3":
    return 28
  elif name == "Pink2":
    return 29
  elif name == "Pink":
    return 30
  elif name == "Gold":
    return 31

