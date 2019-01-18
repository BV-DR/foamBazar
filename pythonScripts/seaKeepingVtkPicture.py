#!/usr/bin/env python
import os
import numpy as np
import vtk

try:
    from Pluto.System import Chrono
except:
    print('Chrono module not found , please add "shared-code/python" in your PYTHON_PATH')


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def getAvailableTimeData(caseDir, parallel):
    if parallel:
        d = os.path.join(caseDir, "processor0")
    else:
        d = caseDir
    return sorted([float(f) for f in os.listdir(d) if is_number(f) and os.path.isdir(os.path.join(d, f))])


def getFreeSurfaceActor(vtk_r, scale = [1,1,1], fsRange = None):

    aa = vtk.vtkAssignAttribute()
    aa.SetInputConnection(vtk_r.GetOutputPort())
    aa.Assign('alpha.water', "SCALARS", "POINT_DATA")

    isoContour = vtk.vtkContourFilter()
    isoContour.SetInputConnection(aa.GetOutputPort())
    isoContour.SetValue(0, 0.5)
    # isoContour.SetGenerateTriangles(0)
    isoContour.SetNumberOfContours(1)

    isoContourDataFilter = vtk.vtkCompositeDataGeometryFilter()
    isoContourDataFilter.SetInputConnection(isoContour.GetOutputPort())

    transform = vtk.vtkTransform()
    transform.Scale(*scale)
    tf = vtk.vtkTransformPolyDataFilter()
    tf.SetInputConnection(isoContourDataFilter.GetOutputPort())
    tf.SetTransform(transform)

    # Use an ArrayCalculator to get the iso-contour elevation
    calc = vtk.vtkArrayCalculator()
    calc.SetInputConnection(tf.GetOutputPort())
    calc.SetAttributeModeToUsePointData()
    calc.AddCoordinateScalarVariable("coordsZ", 2)
    calc.SetFunction("coordsZ")
    calc.SetResultArrayName("coordsZ")

    #--- Map the data
    isoContourMapper = vtk.vtkDataSetMapper()
    isoContourMapper.SetInputConnection(calc.GetOutputPort())

    #---Select node property to display
    lutBlueRed = vtk.vtkLookupTable()
    lutBlueRed.SetHueRange(2. / 3., 0.)
    lutBlueRed.SetVectorModeToMagnitude()
    lutBlueRed.Build()

    scalarBar = vtk.vtkScalarBarActor()
    scalarBar.SetLookupTable(lutBlueRed)
    scalarBar.SetWidth(0.05)
    scalarBar.SetHeight(0.25)
    scalarBar.SetTitle("Z(m)")
    scalarBar.GetTitleTextProperty().SetColor(0., 0., 0.)
    scalarBar.GetLabelTextProperty().SetColor(0., 0., 0.)

    isoContourMapper.SelectColorArray("coordsZ")
    # isoContourMapper.SetScalarModeToUsePointData()
    isoContourMapper.SetScalarModeToUsePointFieldData()
    isoContourMapper.SetLookupTable(lutBlueRed)
    if fsRange is not None :
        isoContourMapper.SetScalarRange(*fsRange)
    isoContourMapper.SetUseLookupTableScalarRange(0)

    #---Add actor
    fsActor = vtk.vtkActor()
    fsActor.GetProperty().SetEdgeVisibility(0)
    fsActor.SetMapper(isoContourMapper)
    fsActor.GetProperty().SetOpacity(0.7)
    return fsActor, scalarBar



def getSymPlaneVtkActor(vtk_r, blockIndex, scalarField = "alpha.water", scalarRange = None):
    """return symmetry plane colored by field
    """
    #SymmetryPlane
    symPlane = vtk.vtkExtractBlock()
    symPlane.SetInputConnection(vtk_r.GetOutputPort())
    symPlane.AddIndex(blockIndex)
    symPlaneDataFilter = vtk.vtkCompositeDataGeometryFilter()
    symPlaneDataFilter.SetInputConnection(symPlane.GetOutputPort())
    symPlaneMapper = vtk.vtkDataSetMapper()
    symPlaneMapper.SetInputConnection(symPlaneDataFilter.GetOutputPort())
    symPlaneMapper.SelectColorArray(scalarField)
    #symPlaneMapper.SetScalarModeToUsePointData()
    symPlaneMapper.SetScalarModeToUseCellFieldData()
    symPlaneMapper.SetUseLookupTableScalarRange(0)
    if scalarRange is not None :
        symPlaneMapper.SetScalarRange(*scalarRange)

    symPlaneActor = vtk.vtkActor()
    symPlaneActor.GetProperty().SetEdgeVisibility(1)
    symPlaneActor.SetMapper(symPlaneMapper)
    return symPlaneActor

def getStuctureActor(vtk_r, blockIndex, scalarField = "p_rgh"):
    """Return an actor with the ship structure
    """
    #--------------------------------- Structure
    structureOnly = vtk.vtkExtractBlock()
    structureOnly.SetInputConnection(vtk_r.GetOutputPort())
    structureOnly.AddIndex(blockIndex)

    structureDataFilter = vtk.vtkCompositeDataGeometryFilter()
    structureDataFilter.SetInputConnection(structureOnly.GetOutputPort())

    structureMapper = vtk.vtkDataSetMapper()
    structureMapper.SetInputConnection(structureDataFilter.GetOutputPort())
    structureMapper.SelectColorArray(scalarField)
    #structureMapper.SetScalarModeToUsePointData()
    structureMapper.SetScalarModeToUsePointFieldData()
    structureMapper.SetUseLookupTableScalarRange(0)
    structureMapper.SetScalarRange(-10000, 10000)

    structureActor = vtk.vtkActor()
    structureActor.GetProperty().SetEdgeVisibility(0)
    structureActor.SetMapper(structureMapper)
    return structureActor


def writerFromExt(ext) :
    """Pick correct writter class based on extension
    """
    if ext == ".png":
        writer = vtk.vtkPNGWriter
    elif ext == ".jpg":
        writer = vtk.vtkJPEGWriter
    elif ext == ".bmp":
        writer = vtk.vtkBMPWriter
    elif ext == ".eps":
        writer = vtk.vtkPostScriptWriter
    elif ext == ".tiff":
        writer = vtk.vtkTIFFWriter
    else:
        print("Picture extension not recognized")
    return writer


def printCamera(cam):
    print("Position" , cam.GetPosition())
    print("Focal point" , cam.GetFocalPoint())
    print("Parallel" , cam.GetParallelProjection())
    print("ParallelScale" , cam.GetParallelScale())
    print("ViewAngle" , cam.GetViewAngle())

def setCamera( renderer, camPosition=None, targetPosition=None, viewAngle=30., scale=None, fitView = True, viewUp = [0,0,1], reduceFittedDistance = None ):

    camera = renderer.GetActiveCamera()

    if camPosition is not None : camera.SetPosition(camPosition)
    if targetPosition is not None : camera.SetFocalPoint(targetPosition)
    camera.SetViewUp(viewUp)
    if viewAngle is not None :  camera.SetViewAngle(viewAngle)

    if fitView or camPosition is None or targetPosition is None:
        renderer.ResetCamera()

    if viewAngle is None :
        if fitView :
            d = camera.GetDistance()
            a = camera.GetViewAngle()
            camera.SetParallelScale(d*np.tan(0.5*(a*np.pi/180)))
            camera.SetParallelProjection(True)
        else :
            camera.SetParallelProjection(True)
            camera.SetParallelScale(scale)

    printCamera(camera)

    return




def getMeshPicture( meshFile,
                    pictureFile,
                    cameraArgs = { "camPosition" : None, "viewAngle" : 30, "scale" : 1.0, "targetPosition" : None, "viewUp" : [0,0,1], "fitView" : True },
                    timeList=[0],
                    startInteractive=False, mag=4, parallel="auto",
                    fsArgs = {"scale" : (1, 1, 1), "fsRange" : None },
                    y0Args = {"scalarField" : "alpha.water", "scalarRange" : [0,1] },
                    structArgs = {"scalarField" : "p_rgh", },
                    hullPatch="ship",
                    ):
    """
    Function to generate picture out of openFoam results

    meshFile : mesh file or case directory
    """
    from tqdm import tqdm

    c = Chrono(start=True)
    baseFile, ext = os.path.splitext(pictureFile)[0:2]
    Writer = writerFromExt(ext)
    pathPic = os.path.abspath(os.path.dirname(pictureFile))
    if not os.path.exists(pathPic):
        os.makedirs(pathPic)

    # Add a file name to the path. "p.foam" does not need to be created
    if os.path.isdir(meshFile):
        caseDir = meshFile
        meshFile = os.path.join(caseDir, "p.foam")
    else:
        caseDir = os.path.dirname(meshFile)

    # Automatically choose "SetCaseType" ( parallel vs (serial/notRun))
    if parallel == "auto":
        if os.path.exists(os.path.join(caseDir, "processor0")):
            parallel = True
        else:
            parallel = False

    if type(timeList) == str :
        if timeList == "all":
            timeList = getAvailableTimeData(caseDir, parallel)
        elif timeList == "latest":
            timeList = [getAvailableTimeData(caseDir, parallel)[-1]]
        elif timeList == "latest-1":
            timeList = [getAvailableTimeData(caseDir, parallel)[-2]]
        else :
            raise(Exception("time should be 'latest', 'all', or a float iterable"))

    print("Generating picture for : ", ["{:.2f}".format(i) for i in timeList])

    #--- Read the openFoam file
    vtk_r = vtk.vtkPOpenFOAMReader()
    vtk_r.SetFileName(meshFile)

    if parallel:
        vtk_r.SetCaseType(0)
        print ("Using decomposed case")
    else:
        vtk_r.SetCaseType(1)  # 0 = decomposed case, 1 = reconstructed case
        print ("Using reconstructed case")
    #vtk_r.ReadZonesOn()
    cdp = vtk.vtkCompositeDataPipeline()
    vtk_r.SetDefaultExecutivePrototype(cdp)
    vtk_r.SetDecomposePolyhedra(0)
    vtk_r.CreateCellToPointOn()
    vtk_r.DisableAllPatchArrays()

    vtk_r.SetPatchArrayStatus("internalMesh", 1)
    vtk_r.SetPatchArrayStatus(hullPatch, 1)
    vtk_r.SetPatchArrayStatus("domainY0", 1)

    vtk_r.SetTimeValue(timeList[0])

    print ("Read first")
    vtk_r.Update()  # not mandatory, but just postpone the waiting time
    print ("Read done")

    iter = vtk_r.GetOutput().NewIterator()
    blockDict = {}
    while not iter.IsDoneWithTraversal():
        blockDict[ iter.GetCurrentMetaData().Get(vtk.vtkCompositeDataSet.NAME()) ] = iter.GetCurrentFlatIndex()
        iter.GoToNextItem()

    print (blockDict)

    #--- Renderer
    renderer = vtk.vtkRenderer()

    #--------------------------------- Free-surface (ISO alpha = 0.5)
    if fsArgs is not None :
        fsActor, scalarBar = getFreeSurfaceActor(vtk_r, **fsArgs)
        renderer.AddActor(fsActor)  # Add the mesh to the view
        renderer.AddActor(scalarBar)

    #--------------------------------- Ship surface
    if structArgs is not None :
        structureActor = getStuctureActor(vtk_r, blockIndex = blockDict[hullPatch], **structArgs )
        renderer.AddActor(structureActor)  # Add the mesh to the view

    if y0Args is not None :
        symActor = getSymPlaneVtkActor(vtk_r, blockIndex = blockDict["domainY0"], **y0Args )
        renderer.AddActor(symActor)  # Add the mesh to the view

    renderer.SetBackground(1, 1, 1)  # White background

    #--- Rendering windows
    renWin = vtk.vtkRenderWindow()

    # Avoid displaying interactive window
    if not startInteractive :
        renWin.SetOffScreenRendering(1)

    renWin.AddRenderer(renderer)
    renWin.SetSize(1650, 1050)

    # Set view point
    setCamera( renderer, **cameraArgs )

    # To get interactive windows
    if startInteractive:
        vtk_r.Modified()
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        iren.Start()

    else:
        for itime, time in tqdm(enumerate(timeList)):
            #exe.SetUpdateTimeStep( 0, time )
            #vtk_r.SetTimeValue( time )
            # vtk_r.Update()
            vtk_r.UpdateTimeStep(time)
            vtk_r.Modified()  # Require to update the time step data

            renWin.Render()
            w2if = vtk.vtkRenderLargeImage()
            w2if.SetMagnification(mag)   # => Resoulition of the picture

            w2if.SetInput(renderer)
            w2if.Update()
            pictureFile = "{:}_{:03}{:}".format(baseFile, itime, ext)

            writer = Writer()
            writer.SetFileName(pictureFile)
            writer.SetInputConnection(w2if.GetOutputPort())
            writer.Write()
            #c.printLap('Chrono t={} {}'.format(time, time))


if __name__ == "__main__":

    """
       Example of use
       seaKeepingVtkPicture -case . -output anim.png -time 0,1,2
    """
    import argparse
    parser = argparse.ArgumentParser(description='openFoam vtk picture')
    parser.add_argument('-case', '-c', help='Quantity to plot', type=str,  default=".")
    parser.add_argument('-time', '-t', help='Time to plot (use "," as separator)', type=str,  default="0")
    parser.add_argument('-output', '-o', help='Picture to generate', type=str,  default="none.png")
    parser.add_argument('-interactive', '-i', help='Interactive 3D view', action="store_true")
    args = parser.parse_args()

    if args.time not in ['all', "latest", "latest-1"]:
        args.time = [float(i) for i in args.time.split(',')]

    getMeshPicture( args.case,
                    pictureFile = args.output,
                    cameraArgs = {
                                   "camPosition" : (20., 40, 20),
                                   "targetPosition" : (2, 2, 2),
                                   "viewAngle" : 30.,
                                   "fitView" : True,
                                  },
                    timeList = "all", # np.arange(0.5, 2, 0.2),
                    mag = 6,
                    fsArgs = { 'fsRange' : [-0.15 , 0.15],  } ,
                    y0Args = { "scalarField" : "alpha.water", },
                    startInteractive = args.interactive)
