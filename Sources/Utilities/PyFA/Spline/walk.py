import math
import Const

#===============================================================================
# Walk from one object to another, avoiding all objects in the graph
#
# Parameters:
#   - x1, y1, ... x6, y6:  coordinates the way Ivan introduced them
#   - obj_list:            list of all objects
# Returns:
#   - x, y:                coordinates with all steps from one object to another
#-------------------------------------------------------------------------------
def walk(x1, y1, x2, y2, x5, y5, x6, y6, obj_list, spl_list, stride):

  # Walk
  x    = []
  y    = []
  dist = []
  keep = []

  x.append(x1)
  y.append(y1)

  x.append(x2)
  y.append(y2)

  #-----------
  #
  # Main loop
  #
  #-----------
  for i in range(0, 256):

    #--------------------------   3 2 1
    # Set eight possible direc    4 c 0
    #--------------------------   5 6 7
    step_x    = []
    step_y    = []
    step_dist = []

    n_direc = Const.WALK_DIRECTIONS
    for j in range(n_direc):
      a = j * 2 * math.pi / (n_direc)  # angle
      cos_a = math.cos(a)
      sin_a = math.sin(a)
      scaling = 1.0 / max(abs(cos_a), abs(sin_a))
      dx = stride * cos_a * scaling
      dy = stride * sin_a * scaling

      step_x.append(x[-1] + dx)
      step_y.append(y[-1] + dy)

    tol = Const.UNIT_BOX_HEIGHT / 2.0

    #-----------------------------------------------
    # Eliminate steps too close to existing splines
    #-----------------------------------------------
    eliminate_steps = []
    for s in range(len(step_x)):
      for o in range(len(spl_list)):
        for n in range(2, spl_list[o].N_Points() - 2):
          if step_x[s] >= (spl_list[o].x[n] - tol * 0.5) and \
             step_x[s] <= (spl_list[o].x[n] + tol * 0.5) and \
             step_y[s] >= (spl_list[o].y[n] - tol * 0.5) and \
             step_y[s] <= (spl_list[o].y[n] + tol * 0.5):
            eliminate_steps.append(s)
    eliminate_steps = list(set(eliminate_steps))
    eliminate_steps.sort(reverse = True)
    for e in range(len(eliminate_steps)):
      step_x.pop(eliminate_steps[e])
      step_y.pop(eliminate_steps[e])

    #---------------------------------------------
    # Eliminate steps which would fall in objects
    #---------------------------------------------
    eliminate_steps = []
    for o in range(len(obj_list)):
      for s in range(len(step_x)):
        if step_x[s] >= (obj_list[o].x0                 - tol) and \
           step_x[s] <= (obj_list[o].x0 + obj_list[o].w + tol) and \
           step_y[s] >= (obj_list[o].y0                 - tol) and \
           step_y[s] <= (obj_list[o].y0 + obj_list[o].h + tol):
          eliminate_steps.append(s)
    eliminate_steps = list(set(eliminate_steps))
    eliminate_steps.sort(reverse = True)
    for e in range(len(eliminate_steps)):
      step_x.pop(eliminate_steps[e])
      step_y.pop(eliminate_steps[e])

    #-------------------------------------
    # Eliminate steps which would go back
    #  (too close to its own old steps)
    #-------------------------------------
    eliminate_steps = []
    for s in range(len(step_x)):
      if step_x[s] >= (x[-2] - tol) and   \
         step_x[s] <= (x[-2] + tol) and   \
         step_y[s] >= (y[-2] - tol) and   \
         step_y[s] <= (y[-2] + tol):
        eliminate_steps.append(s)
    eliminate_steps = list(set(eliminate_steps))
    eliminate_steps.sort(reverse = True)
    for e in range(len(eliminate_steps)):
      step_x.pop(eliminate_steps[e])
      step_y.pop(eliminate_steps[e])

    #-----------------------------------------
    # From the remaining (possible) steps, do
    # find the one closest to the destination
    #-----------------------------------------
    for s in range(len(step_x)):
      dx = step_x[s] - x5
      dy = step_y[s] - y5
      step_dist.append(math.sqrt(dx*dx + dy*dy))

    # Index of direction with minimum distance
    if len(step_dist) > 0:
      min_dist = step_dist.index(min(step_dist))

      x.   append(step_x[min_dist])
      y.   append(step_y[min_dist])
      dist.append(min(step_dist))
      keep.append(True)

      # Check if converged
      if dist[-1] < (stride * 0.25):
        x = x[:-2]
        y = y[:-2]
        break

      # Check if it wobbles (only if you are close)
      if len(dist) > 2:
        if dist[-1] < (stride):
          if dist[-1] > dist[-2]:
            x = x[:-3]
            y = y[:-3]
            break

  x.   append(x5)
  y.   append(y5)
  keep.append(True)

  x.   append(x6)
  y.   append(y6)
  keep.append(True)

  return x, y

