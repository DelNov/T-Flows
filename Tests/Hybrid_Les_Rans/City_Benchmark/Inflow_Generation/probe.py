import os
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Current directory and file name input
current_directory = os.getcwd()
file_name = input(f"Enter the file name (in {current_directory}): ")
file_path = os.path.join(current_directory, file_name)

if not file_path.endswith('.vtu'):
    print("Error: The file must be a .vtu file.")
    exit()

# Use the selected file path
dataReader = XMLUnstructuredGridReader(registrationName=file_path, FileName=[file_path])

# Properties modified on dataReader
dataReader.CellArrayStatus = ['Turbulent Dissipation [m^2/s^3]',
                              'Turbulent Kinetic Energy [m^2/s^2]',
                              'Turbulent Quantity F22 [1]',
                              'Turbulent Quantity Zeta [1]',
                              'Velocity [m/s]']
dataReader.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
dataReaderDisplay = Show(dataReader, renderView1, 'UnstructuredGridRepresentation')

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# update the view to ensure updated data information
renderView1.Update()

# Get the bounding box of the domain
bounds = dataReader.GetDataInformation().GetBounds()

# Calculate the mean y-coordinate
y_mean = (bounds[2] + bounds[3]) / 2  # bounds[2] is ymin, bounds[3] is ymax
z_min  =  bounds[4]
z_max  =  bounds[5]

# Create a new 'Plot Over Line' for each coordinate
plotOverLine = PlotOverLine(registrationName=f'PlotOverLine', Input=dataReader)
plotOverLine.Point1 = [0.501, y_mean+0.001, z_min]
plotOverLine.Point2 = [0.501, y_mean+0.001, z_max]

# Properties modified on plotOverLine1
plotOverLine.SamplingPattern = 'Sample At Segment Centers'

# Save data to a file
output_file_path = file_path.replace('.vtu', f'.csv')
SaveData(output_file_path, proxy=plotOverLine, ChooseArraysToWrite=1,
         PointDataArrays=['Turbulent Dissipation [m^2/s^3]',
                          'Turbulent Kinetic Energy [m^2/s^2]', 
                          'Turbulent Quantity F22 [1]',
                          'Turbulent Quantity Zeta [1]', 
                          'Velocity [m/s]',
                          'arc_length',
                          'vtkValidPointMask'])

# Print success message for each file
print(f"Data successfully saved to {output_file_path}")

