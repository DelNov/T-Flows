#===============================================================================
#   This is to create a single sphere in Blender
#
#   You can run it interactivelly with:
#   > blender --background --python ./sphere.py
#
#   or:
#   > blender -b -P ./sphere.py
#-------------------------------------------------------------------------------

import bpy

# Delete initial cube in Blender
while bpy.data.objects:
  bpy.data.objects.remove(bpy.data.objects[0], do_unlink=True)

# Set some parameters ...
S  = 3
R  = 0.25

bpy.ops.mesh.primitive_ico_sphere_add(subdivisions   = S,               \
                                      radius         = R,               \
                                      enter_editmode = False,           \
                                      align          ='WORLD',          \
                                      location       =(0.5, 0.5, 0.5),  \
                                      scale          =(1, 1, 1))

# Finally export what you got
bpy.ops.export_mesh.stl(filepath='sphere.stl')

