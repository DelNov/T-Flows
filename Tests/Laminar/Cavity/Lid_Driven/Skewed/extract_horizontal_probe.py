# trace generated using paraview version 5.10.0
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#--------------------------------------
# definitions for this particular case
#--------------------------------------
from math import *

# result file name
name_in  = 'skewed-ts012000.pvtu'
name_out = 'vel_horizontal_probe.csv'

# domain characteristics
angle = 45  # should be the same as the angle in skewed.geo file
L = 1.0
W = 0.2
S = L / 10000.0

# point positions
x1 = cos(angle/180.0 * pi) * (L*0.5)
x2 = x1 + L
y  = W*0.5
z1 = sin(angle/180.0 * pi) * (L*0.5)
z2 = z1


#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Partitioned Unstructured Grid Reader'
myResults = XMLPartitionedUnstructuredGridReader(registrationName=name_in, FileName=[name_in])

# Properties modified on myResults
myResults.CellArrayStatus = ['Velocity [m/s]']
myResults.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
myResultsDisplay = Show(myResults, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
myResultsDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clean to Grid'
cleantoGrid1 = CleantoGrid(registrationName='CleantoGrid1', Input=myResults)

# show data in view
cleantoGrid1Display = Show(cleantoGrid1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cleantoGrid1Display.Representation = 'Surface'

# hide data in view
Hide(myResults, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Data to Point Data'
cellDatatoPointData1 = CellDatatoPointData(registrationName='CellDatatoPointData1', Input=cleantoGrid1)

# show data in view
cellDatatoPointData1Display = Show(cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellDatatoPointData1Display.Representation = 'Surface'

# hide data in view
Hide(cleantoGrid1, renderView1)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Plot Over Line'
plotOverLine1 = PlotOverLine(registrationName='PlotOverLine1', Input=cellDatatoPointData1)

# Properties modified on plotOverLine1
plotOverLine1.Point1 = [x1+S, y, z1]
plotOverLine1.Point2 = [x2-S, y, z2]

# show data in view
plotOverLine1Display = Show(plotOverLine1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
plotOverLine1Display.Representation = 'Surface'

# Create a new 'Line Chart View'
lineChartView1 = CreateView('XYChartView')

# show data in view
plotOverLine1Display_1 = Show(plotOverLine1, lineChartView1, 'XYChartRepresentation')

# get layout
layout1 = GetLayoutByName("Layout #1")

# add view to a layout so it's visible in UI
AssignViewToLayout(view=lineChartView1, layout=layout1, hint=0)

# Properties modified on plotOverLine1Display_1
plotOverLine1Display_1.SeriesPlotCorner = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Velocity [m/s]_Magnitude', '0', 'Velocity [m/s]_X', '0', 'Velocity [m/s]_Y', '0', 'Velocity [m/s]_Z', '0', 'arc_length', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesLineStyle = ['Points_Magnitude', '1', 'Points_X', '1', 'Points_Y', '1', 'Points_Z', '1', 'Velocity [m/s]_Magnitude', '1', 'Velocity [m/s]_X', '1', 'Velocity [m/s]_Y', '1', 'Velocity [m/s]_Z', '1', 'arc_length', '1', 'vtkValidPointMask', '1']
plotOverLine1Display_1.SeriesLineThickness = ['Points_Magnitude', '2', 'Points_X', '2', 'Points_Y', '2', 'Points_Z', '2', 'Velocity [m/s]_Magnitude', '2', 'Velocity [m/s]_X', '2', 'Velocity [m/s]_Y', '2', 'Velocity [m/s]_Z', '2', 'arc_length', '2', 'vtkValidPointMask', '2']
plotOverLine1Display_1.SeriesMarkerStyle = ['Points_Magnitude', '0', 'Points_X', '0', 'Points_Y', '0', 'Points_Z', '0', 'Velocity [m/s]_Magnitude', '0', 'Velocity [m/s]_X', '0', 'Velocity [m/s]_Y', '0', 'Velocity [m/s]_Z', '0', 'arc_length', '0', 'vtkValidPointMask', '0']
plotOverLine1Display_1.SeriesMarkerSize = ['Points_Magnitude', '4', 'Points_X', '4', 'Points_Y', '4', 'Points_Z', '4', 'Velocity [m/s]_Magnitude', '4', 'Velocity [m/s]_X', '4', 'Velocity [m/s]_Y', '4', 'Velocity [m/s]_Z', '4', 'arc_length', '4', 'vtkValidPointMask', '4']

# update the view to ensure updated data information
lineChartView1.Update()

# save data
SaveData(name_out, proxy=plotOverLine1, PointDataArrays=['Velocity [m/s]', 'arc_length', 'vtkValidPointMask'],
    AddMetaData=0)

# layout/tab size in pixels
layout1.SetSize(1047, 793)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [0.8535533905932735, 0.1, 3.94399911709714]
renderView1.CameraFocalPoint = [0.8535533905932735, 0.1, 0.35355339059327373]
renderView1.CameraParallelScale = 0.9292757344261569
