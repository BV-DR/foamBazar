from inputFiles.ofDictionary import ofDictionary
import sys
from PyFoam.Basics.DataStructures import DictProxy

"""
  Convenience class to simply write DecomposeParDict
"""

class snappyHexMeshDict(ofDictionary):
    """
       snappyHexMeshDict dictionary
    """
    def __init__(self , root, fdir, fid="snappyHexMeshDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)
        if self.exists: return
        self["#inputMode"] = "overwrite"
        self["castellatedMesh"] = True
        self["snap"] = True
        self["addLayers"] = True

        geometry = DictProxy()
        self["geometry"] = geometry
        
        castel = DictProxy()
        castel["maxLocalCells"] = 1000000
        castel["maxGlobalCells"] = 10000000
        castel["minRefinementCells"] = 0
        castel["nCellsBetweenLevels"] = 4        
        castel["resolveFeatureAngle"] = 30
        castel["allowFreeStandingZoneFaces"] = True
        castel["locationInMesh"] = "LOCATION_IN_MESH"
        castel["features"] = {}
        castel["refinementSurfaces"] = {}
        castel["refinementRegions"] = {}
        self["castellatedMeshControls"] = castel
        
        snapControl = DictProxy()
        snapControl["nSmoothPatch"] = 3
        snapControl["tolerance"] = 1
        snapControl["nSolveIter"] = 100
        snapControl["nRelaxIter"] = 5
        snapControl["nFeatureSnapIter"] = 10
        snapControl["implicitFeatureSnap"] = False
        snapControl["explicitFeatureSnap"] = True
        snapControl["multiRegionFeatureSnap"] = False 
        self["snapControls"] = snapControl

        addLayersControl = DictProxy()
        addLayersControl["relativeSizes"] = True
        addLayersControl["nMedialAxisIter"] = 10
        addLayersControl["expansionRatio"] = 1.0
        addLayersControl["finalLayerThickness"] = 0.3
        addLayersControl["minThickness"] = 0.01
        addLayersControl["nGrow"] = 0
        addLayersControl["featureAngle"] = 30
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
        patch = DictProxy()
        addLayersControl["layers"] = patch
        #if ofp:
        #    addLayersControl["meshShrinker"] = "displacementMotionSolver"
        #    addLayersControl["solver"] = "displacementLaplacian"
        #    dlc = DictProxy()
        #    dlc["diffusivity"] = "quadratic inverseDistance 1("+stlname+")"
        #    addLayersControl["displacementLaplacianCoeffs"] = dlc
        self["addLayersControls"] = addLayersControl
            
        meshQuality = DictProxy()
        meshQuality["maxNonOrtho"] = 65
        meshQuality["maxBoundarySkewness"] = 20
        meshQuality["maxInternalSkewness"] = 4
        meshQuality["maxConcave"] = 80
        meshQuality["minVol"] = 1e-13
        meshQuality["minTetQuality"] = 1e-15
        meshQuality["minArea"] = -1
        meshQuality["minTwist"] = 0.05
        meshQuality["minDeterminant"] = 0.001
        meshQuality["minFaceWeight"] = 0.05
        meshQuality["minVolRatio"] = 0.01
        meshQuality["minTriangleTwist"] = -1
        meshQuality["nSmoothScale"] = 4
        meshQuality["errorReduction"] = 0.75
        meshQuality["minVolCollapseRatio"] = 0.1
        meshQuality["relaxed"] = { "maxNonOrtho" : 75 }
        self["meshQualityControls"] = meshQuality
        self["mergeTolerance"] = 1e-6
        self["debug"] = 0

        # update according to user input
        self.update(**kwargs)
        
    def update(self,**kwargs):
        self.content.update(**kwargs)
        pass

#*** Main execution start here *************************************************
if __name__ == "__main__":
    print(str(sys.argv[0:]))
    print(snappyHexMeshDict("test"))



