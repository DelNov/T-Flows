class Spline():

  def __init__(self, obj1, obj2, line, deep):

    self.object1   = obj1
    self.object2   = obj2
    self.line_type = line
    self.depth     = deep

    self.x = []
    self.y = []

  def N_Points(self):
    return len(self.x)

