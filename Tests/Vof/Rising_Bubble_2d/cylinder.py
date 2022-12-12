#===============================================================================
#   This is to create a cylinder in Blender
#
#   You can run it interactivelly with:
#   > blender --background --python ./cylinder.py
#-------------------------------------------------------------------------------

import bpy

bpy.ops.object.delete(use_global=False, confirm=False)

N  = 60
PI = 3.14159265359
A1 = PI / 2.0
A2 = 2.0*PI / N / 2.0
A2 = 0.0
bpy.ops.mesh.primitive_cylinder_add(vertices=60,            \
                                    radius=0.199,           \
                                    depth=2,                \
                                    enter_editmode=False,   \
                                    align='WORLD',          \
                                    location=(0.5, 0, 0),   \
                                    rotation=(A1,  A2, 0),  \
                                    scale=(1, 1, 1))

# Finally export what you got
bpy.ops.export_mesh.stl(filepath='cylinder.stl')

