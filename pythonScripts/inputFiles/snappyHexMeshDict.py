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
    def __init__(self , case, castellatedMesh=True, snap=True, addLayers=False, stlname="body",
                        locationInMesh=[1e-2,5.0,0.5], nCellsBetweenLevels=4, edgeLvl=1, hullLvl=[1,1],
                        resolveFeatureAngle=45, refinementLength=None, allowFreeStandingZoneFaces=True,
                        snapTol=4., nSolveIter=30, relativeSizes=False, nSurfaceLayers=3, expansionRatio=1.3, finalLayerThickness=0.1,
                        minThickness=0.02, featureAngle=89, maxNonOrtho=55, minTwist=0.05, nSmoothScale=5,
                        errorReduction=0.85, shipPatches=None, noLayers=None):
        WriteParameterFile.__init__(self,  name = join(case, "system" , "snappyHexMeshDict" )  )
        
        stlname = stlname.split('.stl')[0] #remove .stl extension
        
        self["#inputMode"] = "overwrite"
      
        self["castellatedMesh"] = castellatedMesh
        self["snap"] = snap
        self["addLayers"] = addLayers

        geometry = DictProxy()
        body = DictProxy() 
        body["type"]  = "triSurfaceMesh"
        body["name"]  = stlname
        body["patchInfo"]  = { "type" : "wall" }
        geometry['"'+stlname+'.stl"'] = body
        self["geometry"] = geometry
        
        castel = DictProxy()
        castel["maxLocalCells"] = 1000000
        castel["maxGlobalCells"] = 10000000
        castel["minRefinementCells"] = 0
        castel["nCellsBetweenLevels"] = nCellsBetweenLevels
        
        castel["locationInMesh"] = "({} {} {})".format(*locationInMesh)
        castel["features"] = [{'file': '"'+stlname+'.eMesh"', 'level': edgeLvl}]
        castel["refinementSurfaces"] = { stlname : { "level" : "({} {})".format(*hullLvl) } }
        
        castel["resolveFeatureAngle"] = resolveFeatureAngle
        
        if refinementLength is not None:
            castel["refinementRegions"] = { stlname : { "mode" : "distance", "levels" : "(({:.3f} 1))".format(refinementLength) } }
        else:
            castel["refinementRegions"] = {}

        castel["allowFreeStandingZoneFaces"] = allowFreeStandingZoneFaces
        self["castellatedMeshControls"] = castel
        
        snapControl = DictProxy()
        snapControl["nSmoothPatch"] = 3
        snapControl["tolerance"] = snapTol
        snapControl["nSolveIter"] = nSolveIter
        snapControl["nRelaxIter"] = 5
        snapControl["nFeatureSnapIter"] = 10
        snapControl["implicitFeatureSnap"] = False
        snapControl["explicitFeatureSnap"] = True
        snapControl["multiRegionFeatureSnap"] = True
        self["snapControls"] = snapControl
            
        

        addLayersControl = DictProxy()
        addLayersControl["relativeSizes"] = relativeSizes
        
        patch = DictProxy()
        if shipPatches is not None:
            if  len(shipPatches)>1:
                if noLayers==None: noLayers = []
                elif len(noLayers)>0: print("    disable layers on patch(es):", noLayers)
                for sp in shipPatches:
                    if sp in noLayers: continue
                    patch["ship_"+str(sp)] = { "nSurfaceLayers" : nSurfaceLayers }
            else:
                patch["ship"] = { "nSurfaceLayers" : nSurfaceLayers }
        else:
            patch[stlname] = { "nSurfaceLayers" : nSurfaceLayers }
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
        self["addLayersControls"] = addLayersControl
            
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
        self["meshQualityControls"] = mesh

        self["mergeTolerance"] = 1e-6
        self["debug"] = 0

if __name__ == "__main__" : 
   print(SnappyHexMeshDict("test"))



