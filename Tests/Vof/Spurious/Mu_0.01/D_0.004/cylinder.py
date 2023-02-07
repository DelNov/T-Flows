#===============================================================================
#   This is to create a cylinder in Blender
#
#   You can run it interactivelly with:
#   > blender --background --python ./cylinder.py
#
#   or:
#   > blender -b -P ./cylinder.py
#-------------------------------------------------------------------------------

import bpy

while bpy.data.objects:
  bpy.data.objects.remove(bpy.data.objects[0], do_unlink=True)

R  =  0.002
D  =  0.05
X  =  0.005
N  = 60
PI = 3.14159265359
A1 = PI / 2.0
A2 = 2.0*PI / N / 2.0
A2 = 0.0
bpy.ops.mesh.primitive_cylinder_add(vertices=N,             \
                                    radius=R,               \
                                    depth=D,                \
                                    enter_editmode=False,   \
                                    align='WORLD',          \
                                    location=(X, 0.0, X),   \
                                    rotation=(A1,  A2, 0),  \
                                    scale=(1, 1, 1))

# Finally export what you got
bpy.ops.export_mesh.stl(filepath='cylinder.stl')

