import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'ADIOS2VTXReader'
ubp = ADIOS2VTXReader(registrationName='u.bp', FileName='../results/u.bp')

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
ubpDisplay = Show(ubp, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
ubpDisplay.Representation = 'Surface'
ubpDisplay.ColorArrayName = [None, '']
ubpDisplay.SelectTCoordArray = 'None'
ubpDisplay.SelectNormalArray = 'None'
ubpDisplay.SelectTangentArray = 'None'
ubpDisplay.OSPRayScaleArray = 'u_n'
ubpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
ubpDisplay.SelectOrientationVectors = 'None'
ubpDisplay.ScaleFactor = 0.10000000000000003
ubpDisplay.SelectScaleArray = 'None'
ubpDisplay.GlyphType = 'Arrow'
ubpDisplay.GlyphTableIndexArray = 'None'
ubpDisplay.GaussianRadius = 0.005000000000000001
ubpDisplay.SetScaleArray = ['POINTS', 'u_n']
ubpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
ubpDisplay.OpacityArray = ['POINTS', 'u_n']
ubpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
ubpDisplay.DataAxesGrid = 'GridAxesRepresentation'
ubpDisplay.PolarAxes = 'PolarAxesRepresentation'
ubpDisplay.ScalarOpacityUnitDistance = 0.09531842929969367
ubpDisplay.OpacityArrayName = ['POINTS', 'u_n']
ubpDisplay.SelectInputVectors = ['POINTS', 'u_n']
ubpDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
ubpDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
ubpDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(ubpDisplay, ('POINTS', 'u_n', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
ubpDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
ubpDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'u_n'
u_nLUT = GetColorTransferFunction('u_n')

# get opacity transfer function/opacity map for 'u_n'
u_nPWF = GetOpacityTransferFunction('u_n')

# get 2D transfer function for 'u_n'
u_nTF2D = GetTransferFunction2D('u_n')

# change representation type
ubpDisplay.SetRepresentationType('Outline')

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=ubp,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.ScaleFactor = 0.10000000000000003
glyph1.GlyphTransform = 'Transform2'

# Properties modified on glyph1
glyph1.OrientationArray = ['POINTS', 'u_n']
glyph1.ScaleArray = ['POINTS', 'u_n']

# show data in view
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.ColorArrayName = ['POINTS', 'u_n']
glyph1Display.LookupTable = u_nLUT
glyph1Display.SelectTCoordArray = 'None'
glyph1Display.SelectNormalArray = 'None'
glyph1Display.SelectTangentArray = 'None'
glyph1Display.OSPRayScaleArray = 'u_n'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.ScaleFactor = 0.10000000000100001
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.GaussianRadius = 0.00500000000005
glyph1Display.SetScaleArray = ['POINTS', 'u_n']
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityArray = ['POINTS', 'u_n']
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.PolarAxes = 'PolarAxesRepresentation'
glyph1Display.SelectInputVectors = ['POINTS', 'u_n']
glyph1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# show color bar/color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(774, 785)

# current camera placement for renderView1
renderView1.CameraPosition = [2.684267964989148, 2.0190645613395812, -1.6015702913275676]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
renderView1.CameraViewUp = [-0.278813078701974, 0.8922579367136417, 0.35516058553313873]
renderView1.CameraParallelScale = 0.8783332583601864

# save animation
SaveAnimation('../results/Flow_Ani.avi', renderView1, ImageResolution=[772, 784],
    FrameRate=30,
    FrameWindow=[0, 1])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(774, 785)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [2.684267964989148, 2.0190645613395812, -1.6015702913275676]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
renderView1.CameraViewUp = [-0.278813078701974, 0.8922579367136417, 0.35516058553313873]
renderView1.CameraParallelScale = 0.8783332583601864

