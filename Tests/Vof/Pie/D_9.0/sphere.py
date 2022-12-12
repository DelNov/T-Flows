#===============================================================================
#   This is to create a single sphere in Blender
#
#   You can run it interactivelly with:
#   > blender --background --python ./sphere.py
#-------------------------------------------------------------------------------
# Original parameters from Mijail:
#
# Number of ellipsoids
#   1
#
# Ellipsoid principal axes
#   0.0045   0.0045   0.0045
#
# Ellipsoid centre
#   0.0      0.009   0.0
#-------------------------------------------------------------------------------

import bpy

# Delete initial cube in Blender
bpy.ops.object.delete(use_global=False, confirm=False)

# Set some parameters ...
D  = 0.009
S  = 3
R  = D / 2.0
X  = 0.0
Y  = D
Z  = 0.0

bpy.ops.mesh.primitive_ico_sphere_add(subdivisions   = S,         \
                                      radius         = R,         \
                                      enter_editmode = False,     \
                                      align          ='WORLD',    \
                                      location       =(X, Y, Z),  \
                                      scale          =(1, 1, 1))

# Finally export what you got
bpy.ops.export_mesh.stl(filepath='sphere.stl')

