# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
stefanfrontts000 = FindSource('stefan-front-ts000*')

# create a new 'Transform'
transform2 = Transform(Input=stefanfrontts000)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=transform2.Transform)

# find source
stefants000 = FindSource('stefan-ts000*')

# find source
transform1 = FindSource('Transform1')

# Properties modified on transform2.Transform
transform2.Transform.Scale = [1.0, 0.011, 0.011]

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1149, 920]

# show data in view
transform2Display = Show(transform2, renderView1)

# trace defaults for the display properties.
transform2Display.Representation = 'Surface'

# hide data in view
Hide(stefanfrontts000, renderView1)

# update the view to ensure updated data information
renderView1.Update()

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).