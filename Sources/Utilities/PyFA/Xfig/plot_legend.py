import Const
from Objects import Function, Module, Program, Subroutine
from Xfig.plot_spline_legend import plot_spline_legend
from Xfig.plot_text          import plot_text
from Xfig.find_width         import find_width
from Xfig.plot_object_name   import plot_object_name
from Xfig.write_comment      import write_comment

#===============================================================================
# Function to plot legend
#
# Parameters:
#   - file:         Xfig file's handle
#   - obj_list:     list of all objects
# Returns:
#   - nothing
#-------------------------------------------------------------------------------
def plot_legend(file, obj_list, x0, y0):

  # Print a comment
  write_comment(file, "Legend", 5)

  p = ""  # path
  u = []  # use
  v = []  # var
  m = []  # methods
  c = []  # calls
  mobject = Module    ("       Module       ", p, u, v, m, c, [])
  mobject.x0 = x0
  mobject.y0 = y0
  sobject = Subroutine("     Subroutine     ", p, u, v, m, c, [])
  sobject.x0 = x0
  sobject.y0 = y0 + Const.UNIT_BOX_HEIGHT
  fobject = Function  ("      Function      ", p, u, v, m, c, [], "")
  fobject.x0 = x0
  fobject.y0 = y0 + Const.UNIT_BOX_HEIGHT * 2
  pobject = Program   ("      Program       ", p, u, v, m, c, [])
  pobject.x0 = x0
  pobject.y0 = y0 + Const.UNIT_BOX_HEIGHT * 3
  text_height = Const.UNIT_BOX_HEIGHT
  text_width  = find_width(sobject)

  # Boxes
  plot_object_name(file, mobject)
  plot_object_name(file, sobject)
  plot_object_name(file, fobject)
  plot_object_name(file, pobject)

  # Splines
  plot_spline_legend(file,                                           \
                     x0,  y0+Const.UNIT_BOX_HEIGHT*5.5, text_width,  \
                     "Continuous")
  plot_spline_legend(file,                                           \
                     x0, y0+Const.UNIT_BOX_HEIGHT*7.0, text_width,   \
                     "Dashed")
  plot_text(file, "Center",                   \
            x0 + text_width*0.5,              \
            y0 + Const.UNIT_BOX_HEIGHT*4.5    \
               + Const.UNIT_BOX_HEIGHT*0.75,  \
            " Use statements ", Const.FONT_HEADER, "Black", 10)
  plot_text(file, "Center",                   \
            x0 + text_width*0.5,              \
            y0 + Const.UNIT_BOX_HEIGHT*6.0    \
               + Const.UNIT_BOX_HEIGHT*0.75,  \
            " Call statements ", Const.FONT_HEADER, "Black", 10)

