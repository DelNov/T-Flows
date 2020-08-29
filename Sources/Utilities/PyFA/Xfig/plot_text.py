import Const
from Xfig.code_font  import code_font
from Xfig.code_color import code_color

#===============================================================================
# Function to print centered frameless text
#
# Parameters:
#   - file:            Xfig file's handle
#   - alignment:       can be "Left", "Center" or "Right"
#   - x0:              object position on x axis in centimeters
#   - y0:              object position on y axis in centimeters
#   - text:            text to plot (header name)
#   - color:           in which color to print
#   - depth:           layer depth
# Returns:
#   - nothing
# Used by:
#   - function for plotting text in Xfig format
#-------------------------------------------------------------------------------
def plot_text(file, alignment, x0, y0, text, font, color, depth):

  # Alignments can be different
  if   alignment == "Left":
    file.write("4 0 %3d %3d -1 " % (code_color(color), depth))
  elif alignment == "Center":
    file.write("4 1 %3d %3d -1 " % (code_color(color), depth))
  elif alignment == "Right":
    file.write("4 2 %3d %3d -1 " % (code_color(color), depth))

  file.write("%5d" % code_font(font))             # font code
  file.write("%3d" % (Const.FONT_SIZE * 36))      # font size
  file.write(" 0.0000 4 ")
  text_width  = 3                                 # could be any value
  text_height = 3                                 # could be any value
  file.write("%9d" % (text_height * Const.XFIG_SCALE))  # text height xfig units
  file.write("%9d" % (text_width  * Const.XFIG_SCALE))  # text width xfig units
  file.write("%9d %9d" % ((x0)*Const.XFIG_SCALE, (y0)*Const.XFIG_SCALE))

  if alignment == "Center":
    file.write("%s\\001\n" % text)
  else:
    file.write("%s%s\\001\n" % (" ", text))

