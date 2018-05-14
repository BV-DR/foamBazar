import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join

"""
  Convenience class to simply write DecomposeParDict
"""

feature = '''
(
     {
         file "body.eMesh";
         level 1;
     }
)
'''

class SnappyHexMeshDict(WriteParameterFile):
    """
       SnappyHexMeshDict dictionary
    """
    def __init__(self , case, castellatedMesh=True, snap=True, addLayers=False, stlname="body", refinementLength=0.5, nSurfaceLayers=3, expansionRatio=1.3, finalLayerThickness=0.1, minThickness=0.02):
        WriteParameterFile.__init__(self,  name = join(case, "system" , "snappyHexMeshDict" )  )
      
        self["castellatedMesh"] = castellatedMesh
        self["snap"] = snap
        self["addLayers"] = addLayers

        
        geometry = DictProxy()
        body = DictProxy() 
        body["type"]  = "triSurfaceMesh"
        body["name"]  = stlname
        body["patchInfo"]  = { "type" : "wall" }
        geometry[stlname+".stl"] = body
        self["geometry"] = geometry
        
        if castellatedMesh:
            castel = DictProxy()
            castel["maxLocalCells"] = 1000000
            castel["maxGlobalCells"] = 10000000
            castel["minRefinementCells"] = 0
            castel["nCellsBetweenLevels"] = 4
        
            # toto = ' file '+stlname+'.eMesh level 1 '
            castel["features"] = [{'file': '"'+stlname+'.eMesh"', 'level': 1}]
            castel["refinementSurfaces"] = { stlname : { "level" : "(1 1)" } }
            
            castel["resolveFeatureAngle"] = 45
            
            castel["refinementRegions"] = { stlname : { "mode" : "distance", "levels" : "(({:.3f} 1))".format(refinementLength) } }
            
            castel["locationInMesh"] = "(1e-2 5.0 0.5)"
            castel["allowFreeStandingZoneFaces"] = True
            self["castellatedMeshControls"] = castel
        
        if snap:
            snapControl = DictProxy()
            snapControl["nSmoothPatch"] = 3
            snapControl["tolerance"] = 4.0
            snapControl["nSolveIter"] = 30
            snapControl["nRelaxIter"] = 5
            self["snapControls"] = snapControl
        
        if addLayers:
            addLayersControl = DictProxy()
            addLayersControl["relativeSizes"] = False
            addLayersControl["layers"] = { stlname :  { "nSurfaceLayers" : nSurfaceLayers } }
            addLayersControl["nMedialAxisIter"] = 10
            addLayersControl["expansionRatio"] = expansionRatio
            addLayersControl["finalLayerThickness"] = finalLayerThickness
            addLayersControl["minThickness"] = minThickness
            addLayersControl["nGrow"] = 0
            addLayersControl["featureAngle"] = 89
            addLayersControl["nRelaxIter"] = 5
            addLayersControl["nSmoothSurfaceNormals"] = 1
            addLayersControl["nSmoothNormals"] = 3
            addLayersControl["nSmoothThickness"] = 10
            addLayersControl["maxFaceThicknessRatio"] = 0.5
            addLayersControl["maxThicknessToMedialRatio"] = 0.3
            addLayersControl["minMedianAxisAngle"] = 90
            addLayersControl["nBufferCellsNoExtrude"] = 0
            addLayersControl["nLayerIter"] = 50
            self["addLayersControls"] = addLayersControl
            
            
        mesh = DictProxy()
        mesh["maxNonOrtho"] = 55
        mesh["maxBoundarySkewness"] = 20
        mesh["maxInternalSkewness"] = 4
        mesh["maxConcave"] = 80
        mesh["minVol"] = 1e-13
        mesh["minTetQuality"] = 1e-30
        mesh["minArea"] = -1
        mesh["minTwist"] = 0.05
        mesh["minDeterminant"] = 0.001
        mesh["minFaceWeight"] = 0.05
        mesh["minVolRatio"] = 0.01
        mesh["minTriangleTwist"] = -1
        mesh["nSmoothScale"] = 5
        mesh["errorReduction"] = 0.85
        self["meshQualityControls"] = mesh
        
        self["mergeTolerance"] = 1e-6
        

if __name__ == "__main__" : 
   print(SnappyHexMeshDict("test"))



