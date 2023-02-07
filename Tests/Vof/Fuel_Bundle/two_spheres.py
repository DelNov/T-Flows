#===============================================================================
#   This is to create a single sphere in Blender
#
#   You can run it interactivelly with:
#   > blender --background --python ./two_spheres.py
#
#   or:
#   > blender -b -P ./two_spheres.py
#-------------------------------------------------------------------------------

import bpy

# Delete initial cube in Blender
while bpy.data.objects:
  bpy.data.objects.remove(bpy.data.objects[0], do_unlink=True)

# Set some parameters ...
S  = 3
R  = 0.2

bpy.ops.mesh.primitive_ico_sphere_add(subdivisions   = S,                   \
                                      radius         = R,                   \
                                      enter_editmode = False,               \
                                      align          ='WORLD',              \
                                      location       =(0.625, 0.625, 0.4),  \
                                      scale          =(1, 1, 1))

bpy.ops.mesh.primitive_ico_sphere_add(subdivisions   = S,                   \
                                      radius         = R,                   \
                                      enter_editmode = False,               \
                                      align          ='WORLD',              \
                                      location       =(0.675, 0.675, 1.0),  \
                                      scale          =(1, 1, 1))

# Finally export what you got
bpy.ops.export_mesh.stl(filepath='two_spheres.stl')

