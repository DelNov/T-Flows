# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
stefants000 = FindSource('stefan-ts000*')

# create a new 'Transform'
transform1 = Transform(Input=stefants000)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=transform1.Transform)

# set active source
SetActiveSource(transform1)

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1149, 920]

# show data in view
transform1Display = Show(transform1, renderView1)

# trace defaults for the display properties.
transform1Display.Representation = 'Surface'

# hide data in view
Hide(stefants000, renderView1)

# Properties modified on transform1.Transform
transform1.Transform.Scale = [1.0, 0.01, 0.01]

# show data in view
transform1Display = Show(transform1, renderView1)

# reset view to fit data
renderView1.ResetCamera()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# hide data in view
Hide(stefants000, renderView1)

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [0.001, 0.0, 0.0038733505194375175]
renderView1.CameraFocalPoint = [0.001, 0.0, 0.0]
renderView1.CameraParallelScale = 0.001002496882788171

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).