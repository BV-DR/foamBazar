#!/usr/bin/python
from __future__ import print_function
import vtk
import os

try :  from System.Chrono import Chrono
except :  print( 'Chrono module not found , please add "shared-code/python" in your PYTHON_PATH' )


def is_number(s):
   try: float(s); return True
   except ValueError: return False

def getAvailableTimeData(caseDir , parallel ) :
   if parallel : d = os.path.join(caseDir , "processor0")
   else : d = caseDir
   return sorted([float(f) for f in os.listdir(d) if is_number(f) and os.path.isdir(os.path.join(d , f)) ])
   


def getMeshPicture(meshFile,  camPosition , targetPosition , pictureFile, timeList = [0], startInteractive = False, mag = 4 , parallel = "auto", zoom = 1.0 , scale = (1,1,1)):
   """
   Function to generate picture out of openFoam results

   meshFile : mesh file or case directory
   """

   c = Chrono(start = True)

   baseFile, ext = os.path.splitext(pictureFile)[0:2]
   pathPic  = os.path.abspath(os.path.dirname( pictureFile ))
   if not os.path.exists(pathPic) : os.makedirs(pathPic)


   #Add a file name to the path. "p.foam" does not need to be created
   if os.path.isdir(meshFile) :
      caseDir = meshFile
      meshFile = os.path.join(caseDir , "p.foam")
   else :
      caseDir = os.path.dirname(meshFile)

   #Automatically choose "SetCaseType" ( parallel vs (serial/notRun))
   if parallel == "auto" :
      if os.path.exists( os.path.join(caseDir , "processor0") ) : parallel = True
      else : parallel = False

   if timeList == "all" :
      timeList = getAvailableTimeData( caseDir , parallel )
   elif timeList == "latest" :
      timeList = [getAvailableTimeData( caseDir , parallel )[-1]]
   elif timeList == "latest-1" :
      timeList = [getAvailableTimeData( caseDir , parallel )[-2]]


   print (timeList)



   #--- Read the openFoam file
   vtk_r = vtk.vtkPOpenFOAMReader()
   vtk_r.SetFileName(meshFile)

   if parallel : vtk_r.SetCaseType(0)
   else :        vtk_r.SetCaseType(1)    #0 = decomposed case, 1 = reconstructed case
   #vtk_r.ReadZonesOn()   # ?
   cdp = vtk.vtkCompositeDataPipeline()
   vtk_r.SetDefaultExecutivePrototype( cdp )
   vtk_r.SetDecomposePolyhedra(0)
   vtk_r.CreateCellToPointOn()
   vtk_r.DisableAllPatchArrays()
   vtk_r.SetPatchArrayStatus( "internalMesh", 1 )
   vtk_r.SetPatchArrayStatus( "structure", 1 )

   exe = vtk_r.GetExecutive()

   #--------------------------------- Free-surface
   aa = vtk.vtkAssignAttribute()
   aa.SetInputConnection( vtk_r.GetOutputPort() )
   aa.Assign( 'alpha.water', "SCALARS", "POINT_DATA" )

   isoContour = vtk.vtkContourFilter()
   isoContour.SetInputConnection( aa.GetOutputPort() )
   isoContour.SetValue(0, 0.5)
   isoContour.SetGenerateTriangles(0)
   isoContour.SetNumberOfContours(1)

   isoContourDataFilter = vtk.vtkCompositeDataGeometryFilter()
   isoContourDataFilter.SetInputConnection( isoContour.GetOutputPort() )


   transform = vtk.vtkTransform()
   transform.Scale( *scale )
   tf = vtk.vtkTransformPolyDataFilter()
   tf.SetInputConnection(isoContourDataFilter.GetOutputPort())
   tf.SetTransform(transform)


   # Use an ArrayCalculator to get the iso-contour elevation
   calc = vtk.vtkArrayCalculator()
   calc.SetInputConnection(tf.GetOutputPort())
   calc.SetAttributeModeToUsePointData()
   calc.AddCoordinateScalarVariable( "coordsZ", 2 )
   calc.SetFunction("coordsZ")
   calc.SetResultArrayName("coordsZ")

   #--- Map the data
   isoContourMapper = vtk.vtkDataSetMapper()
   isoContourMapper.SetInputConnection(calc.GetOutputPort())

   #---Select node property to display
   lutBlueRed = vtk.vtkLookupTable()
   lutBlueRed.SetHueRange( 2. / 3., 0. )
   lutBlueRed.SetVectorModeToMagnitude()
   lutBlueRed.Build()

   scalarBar = vtk.vtkScalarBarActor()
   scalarBar.SetLookupTable( lutBlueRed )
   scalarBar.SetWidth( 0.05 )
   scalarBar.SetHeight( 0.25 )
   scalarBar.SetTitle( "Z(m)" )
   scalarBar.GetTitleTextProperty().SetColor( 0., 0., 0. )
   scalarBar.GetLabelTextProperty().SetColor( 0., 0., 0. )

   isoContourMapper.SelectColorArray("coordsZ")
   #isoContourMapper.SetScalarModeToUsePointData()
   isoContourMapper.SetScalarModeToUsePointFieldData()
   isoContourMapper.SetLookupTable( lutBlueRed )
   isoContourMapper.SetScalarRange(-1.2 , 1.2)
   isoContourMapper.SetUseLookupTableScalarRange(0)

   #---Add actor
   fsActor = vtk.vtkActor()
   fsActor.GetProperty().SetEdgeVisibility(0)
   fsActor.SetMapper(isoContourMapper)
   fsActor.GetProperty().SetOpacity( 0.7 )


   #--------------------------------- Structure
   structureOnly  = vtk.vtkExtractBlock()
   structureOnly.SetInputConnection(vtk_r.GetOutputPort())
   structureOnly.AddIndex(2)

   structureDataFilter = vtk.vtkCompositeDataGeometryFilter()
   structureDataFilter.SetInputConnection( structureOnly.GetOutputPort() )
   structureMapper = vtk.vtkDataSetMapper()
   structureMapper.SetInputConnection(structureDataFilter.GetOutputPort())
   structureMapper.SelectColorArray("p_rgh")
   structureMapper.SetScalarModeToUsePointData()
   #structureMapper.SetScalarModeToUsePointFieldData()
   structureMapper.SetUseLookupTableScalarRange(0)
   structureMapper.SetScalarRange(-10000 , 10000)

   structureActor = vtk.vtkActor()
   structureActor.GetProperty().SetEdgeVisibility(0)
   structureActor.SetMapper(structureMapper)

   #--- Renderer
   renderer = vtk.vtkRenderer()
   renderer.AddActor(fsActor)        #Add the mesh to the view
   renderer.AddActor(structureActor) #Add the mesh to the view
   renderer.AddActor( scalarBar )
   renderer.SetBackground(1, 1, 1)   #White background
   renderer.GetActiveCamera().SetParallelProjection(0)

   #--- Rendering windows
   renWin = vtk.vtkRenderWindow()
   renWin.AddRenderer(renderer)

   #Avoid displaying interactive window
   if not startInteractive :
      renWin.SetOffScreenRendering(1)

   renWin.SetSize(1200, 1000)

   #Set view point
   camera = renderer.GetActiveCamera()
   camera.SetPosition( camPosition )
   camera.SetFocalPoint( targetPosition )
   camera.SetViewUp( 1. , 1. , 1.)
   camera.Zoom(zoom)
   
   c.printLap ('init')
   vtk_r.Update()
   c.printLap ('Update 1')
   exe.SetUpdateTimeStep( 0, timeList[0] )
   vtk_r.Update()  # not mandatory, but just postpone the waiting time
   vtk_r.Modified()
   c.printLap ('Update {}'.format(timeList[0]) )

   #To get interactive windows
   if startInteractive :
      iren = vtk.vtkRenderWindowInteractor()
      iren.SetRenderWindow(renWin)
      iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
      iren.Start()

   else :

      if ext == ".png" : writer = vtk.vtkPNGWriter()
      elif ext == ".jpg" : writer = vtk.vtkJPEGWriter()
      elif ext == ".bmp" : writer = vtk.vtkBMPWriter()
      elif ext == ".eps" : writer = vtk.vtkPostScriptWriter()
      elif ext == ".tiff" : writer = vtk.vtkTIFFWriter()
      else : print ("Picture extension not recognized")

      for itime, time in enumerate(timeList):
         print ("time" , time)
         if itime != 0:
            exe.SetUpdateTimeStep( 0, time )
            vtk_r.Update()
            vtk_r.Modified()   #Require to update the time step data
   
         renWin.Render()
         w2if = vtk.vtkWindowToImageFilter()
         w2if.SetMagnification(mag)   # => Resoulition of the picture
   
         w2if.SetInput(renWin)
         w2if.Update()
         pictureFile =  "{:}_{:03}{:}".format(baseFile, itime , ext)
   
         writer.SetFileName(pictureFile)
         writer.SetInputConnection(w2if.GetOutputPort())
         writer.Write()




def convert2avi( baseName , clean = False) :
   import subprocess
   subprocess.call( [os.path.join( r"C:\Program Files (x86)\ImageMagick-6.8.1-Q16","ffmpeg") ,  "-framerate",  "10",  "-i" ,  baseName+"_%03d.png" ,"-y", "-vcodec",  "mpeg4" ,  baseName+".avi" ] , shell = False   )
   if clean :
      #TODO : remove all png files
      pass
   return


if __name__ == "__main__" :

   """
      Example of use
      seaKeepingVtkPicture -case . -output anim.png -time 0,1,2
   """
   import argparse
   parser = argparse.ArgumentParser(description='openFoam vtk picture')
   parser.add_argument( '-case' ,   help='Quantity to plot' , type = str,  default = ".")
   parser.add_argument( '-time'  ,  help='Time to plot (use "," as separator)' , type = str,  default = "0")
   parser.add_argument( '-output' ,  help='Picture to generate' , type = str,  default = "none.png")
   parser.add_argument( '-interactive' ,  help='Picture to generate' , action="store_true")
   args = parser.parse_args()
   
   if args.time not in ['all' , "latest", "latest-1"] : args.time = [float(i) for i in args.time.split(',')]
   
   getMeshPicture( args.case, startInteractive = args.interactive,
                   camPosition = (600 , 600 , 300) , targetPosition =  (0.0 ,0.0 , 0.0),
                   timeList = args.time,  pictureFile = args.output , scale = [1.,1.,1.])

