from os.path import join
from inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import DictProxy


class SnappyHexMeshDict(ReadWriteFile):
    """
       SnappyHexMeshDict dictionary
    """
    
    @classmethod
    def Build(cls , case, castellatedMesh=True, snap=True, addLayers=False, stlname="body",
                        patchName='ship', stlPatches=None, locationInMesh=[1e-2,5.0,0.5],
                        nCellsBetweenLevels=4, edgeLvl=1, hullLvl=[1,1], resolveFeatureAngle=45,
                        refinementLength=None, allowFreeStandingZoneFaces=True, snapTol=4.,
                        nSolveIter=30, relativeSizes=False, nSurfaceLayers=3, expansionRatio=1.3,
                        finalLayerThickness=0.1, minThickness=0.02, featureAngle=89, maxNonOrtho=55,
                        minTwist=0.05, nSmoothScale=5, errorReduction=0.85, noLayers=None, OFversion=5):
        """Build SnappyHexMesh file from a few parameters.
        
        Parameters
        ----------
        case : str
            Name of case
            
        castellatedMesh : bool, default True
            Logical defining if castellated mesh should be created
        snap : bool, default True
            Logical defining if snapping should be performed
        addLayers : bool, default True
            Logical defining if boundary layers should be added
          
        stlname : str, default 'body'
            Name of STL file used for snapping
        patchName : str, default 'ship'
            Name of patch created by snapping
        stlPatches : list of str
            TODO
        locationInMesh, list of floats, dimension(3), default [1e-2,5.0,0.5]
            Position of snapping in global mesh
            
        nCellsBetweenLevels : int, default 4
            TODO (option nCellsBetweenLevels)
        edgeLvl : int, default 1
            TODO (option features { file : "stl.eMesh", level : edgeLvl})
        hullLvl : list of int, dimension(2), default [1,1]
            TODO (option refinementSurfaces {patchName : {"level" : hulllvl}})
        resolveFeatureAngle : float, default 45.
            TODO (option resolveFeatureAngle)            
        refinementLength : list of float
            List defining proximity refinements around the body. The length of the list defines the number of proximity refinements.
        allowFreeStandingZoneFaces : 
            TODO (option allowFreeStandingZoneFaces)
            
        snapTol : float, defeult 4.
            Tolerence used for snapping.
        nSolveIter : int, default 30
            TODO (option nSolveIter)
            
        relativeSizes :
            TODO (option relativeSizes)
        nSurfaceLayers : int, default 3
            Number of cells in boundary layer
        expansionRatio : float, default 1.3
            Expansion ratio for boundary layer cells
        finalLayerThickness : flaot, default 0.1
            Thickness of boundary layer
        minThickness : float, default 0.02
            Minimal thickness of boundary layer
        featureAngle : float, default 89.
            TODO (option featureAngle for layers)
        noLayers: list of str
            Disable layers creation on selected patches (e.g. deck)
            
        maxNonOrtho : int, default 55
            TODO (option maxNonOrtho)
        minTwist : float, default 0.05
            TODO (option minTwist)
        nSmoothScale : int, default 5
            TODO (option nSmoothScale)
        errorReduction : float, default 0.85
            TODO (option errorReduction)
        OFversion: int or str, default 5
            OpenFOAM version used
        
        """
        
        res = cls(  name = join(case, getFilePath("snappyHexMeshDict") ), read = False )
        
        stlname = stlname.split('.stl')[0] #remove .stl extension
        
        res["#inputMode"] = "overwrite"
      
        res["castellatedMesh"] = castellatedMesh
        res["snap"] = snap
        res["addLayers"] = addLayers

        geometry = DictProxy()
        body = DictProxy() 
        body["type"]  = "triSurfaceMesh"
        body["name"]  = patchName
        body["patchInfo"]  = { "type" : "wall" }
        geometry['"'+stlname+'.stl"'] = body
        res["geometry"] = geometry
        
        castel = DictProxy()
        castel["maxLocalCells"] = 1000000
        castel["maxGlobalCells"] = 10000000
        castel["minRefinementCells"] = 0
        castel["nCellsBetweenLevels"] = nCellsBetweenLevels
        
        castel["locationInMesh"] = "({} {} {})".format(*locationInMesh)
        castel["features"] = [{'file': '"'+stlname+'.eMesh"', 'level': edgeLvl}]
        castel["refinementSurfaces"] = { patchName : { "level" : "({} {})".format(*hullLvl) } }
        
        castel["resolveFeatureAngle"] = resolveFeatureAngle
        
        if refinementLength is not None:
            nRefL = len(refinementLength)
            lvl = ''
            for i, rl in enumerate(refinementLength):
                lvl += ' ({:.3f} {:d})'.format(rl,nRefL-i)
            castel["refinementRegions"] = { patchName : { "mode" : "distance", "levels" : "("+lvl+")" } }
        else:
            castel["refinementRegions"] = {}

        castel["allowFreeStandingZoneFaces"] = allowFreeStandingZoneFaces
        res["castellatedMeshControls"] = castel
        
        snapControl = DictProxy()
        snapControl["nSmoothPatch"] = 3
        snapControl["tolerance"] = snapTol
        snapControl["nSolveIter"] = nSolveIter
        snapControl["nRelaxIter"] = 5
        snapControl["nFeatureSnapIter"] = 10
        snapControl["implicitFeatureSnap"] = False
        snapControl["explicitFeatureSnap"] = True
        snapControl["multiRegionFeatureSnap"] = False 
        res["snapControls"] = snapControl
            
        

        addLayersControl = DictProxy()
        addLayersControl["relativeSizes"] = relativeSizes
        
        patch = DictProxy()
        if stlPatches is not None:
            if  len(stlPatches)>1:
                if noLayers==None: noLayers = []
                elif len(noLayers)>0: print("    disable layers on patch(es):", noLayers)
                for sp in stlPatches:
                    if sp in noLayers: continue
                    patch[patchName+"_"+str(sp)] = { "nSurfaceLayers" : nSurfaceLayers }
            else:
                patch[patchName] = { "nSurfaceLayers" : nSurfaceLayers }
        else:
            patch[patchName] = { "nSurfaceLayers" : nSurfaceLayers }
        addLayersControl["layers"] = patch
        
        addLayersControl["nMedialAxisIter"] = 10
        addLayersControl["expansionRatio"] = expansionRatio
        addLayersControl["finalLayerThickness"] = finalLayerThickness
        addLayersControl["minThickness"] = minThickness
        addLayersControl["nGrow"] = 0
        addLayersControl["featureAngle"] = featureAngle
        addLayersControl["nRelaxIter"] = 5
        addLayersControl["nSmoothSurfaceNormals"] = 1
        addLayersControl["nSmoothNormals"] = 3
        addLayersControl["nSmoothThickness"] = 10
        addLayersControl["maxFaceThicknessRatio"] = 0.5
        addLayersControl["maxThicknessToMedialRatio"] = 0.3
        addLayersControl["minMedianAxisAngle"] = 90
        addLayersControl["nBufferCellsNoExtrude"] = 0
        addLayersControl["nLayerIter"] = 50
        addLayersControl["nRelaxedIter"] = 20
        if 'p' in str(OFversion).lower():
            addLayersControl["meshShrinker"] = "displacementMotionSolver"
            addLayersControl["solver"] = "displacementLaplacian"
            dlc = DictProxy()
            dlc["diffusivity"] = "quadratic inverseDistance 1("+patchName+")"
            addLayersControl["displacementLaplacianCoeffs"] = dlc
        res["addLayersControls"] = addLayersControl
            
        mesh = DictProxy()
        mesh["maxNonOrtho"] = maxNonOrtho
        mesh["maxBoundarySkewness"] = 20
        mesh["maxInternalSkewness"] = 4
        mesh["maxConcave"] = 80
        mesh["minVol"] = 1e-13
        mesh["minTetQuality"] = 1e-15
        mesh["minArea"] = -1
        mesh["minTwist"] = minTwist
        mesh["minDeterminant"] = 0.001
        mesh["minFaceWeight"] = 0.05
        mesh["minVolRatio"] = 0.01
        mesh["minTriangleTwist"] = -1
        mesh["nSmoothScale"] = nSmoothScale
        mesh["errorReduction"] = errorReduction
        mesh["minVolCollapseRatio"] = 0.1
        mesh["relaxed"] = { "maxNonOrtho" : 75 }
        res["meshQualityControls"] = mesh

        res["mergeTolerance"] = 1e-6
        res["debug"] = 0
        return res

if __name__ == "__main__" : 
   print(SnappyHexMeshDict.Build("test"))



