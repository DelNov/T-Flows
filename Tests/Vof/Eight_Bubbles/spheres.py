#===============================================================================
#   This is to create eight spheres in Blender
#
#   You can run it interactivelly with:
#   > blender --background --python ./spheres.py
#-------------------------------------------------------------------------------

import bpy

# Delete initial cube in Blender
bpy.ops.object.delete(use_global=False, confirm=False)

# Set some parameters ...
S  = 3
P  = 0.2
R1 = 0.19
R2 = 0.17
numb = [-1, 1]

# ... and run!
for i in numb:
  for j in numb:
    for k in numb:
      x = i * P
      y = j * P
      z = k * P
      if i * j * k > 0:
        R = R1
      else:
        R = R2
      bpy.ops.mesh.primitive_ico_sphere_add(subdivisions   = S,         \
                                            radius         = R,         \
                                            enter_editmode = False,     \
                                            align          ='WORLD',    \
                                            location       =(x, y, z),  \
                                            scale          =(1, 1, 1))

# Finally export what you got
bpy.ops.export_mesh.stl(filepath='spheres.stl')


