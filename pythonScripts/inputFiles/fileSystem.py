

from inputFiles import ControlDict, FvSchemes, FvSolution, DecomposeParDict, DynamicMeshDict,  TransportProperties, WaveProperties
from inputFiles.turbulenceProperties import TurbulenceProperties, RASProperties
from inputFiles import BoundaryPressure, BoundaryVelocity, BoundaryPointDisplacement, BoundaryAlpha,BlockMeshDict

default = {
        "controlDict"          : ("system/controlDict", ControlDict),
        "fvSchemes"            : ("system/fvSchemes", FvSchemes ),
        "fvSolution"           : ("system/fvSolution", FvSolution),
        "decomposeParDict"     : ("system/decomposeParDict", DecomposeParDict),
        "dynamicMeshDict"      : ("constant/dynamicMeshDict", DynamicMeshDict),
        "transportProperties"  : ("constant/transportProperties", TransportProperties),
        "waveProperties"       : ("constant/waveProperties", WaveProperties),
        "rasProperties"        : ("constant/RASProperties", RASProperties),
        "turbulenceProperties" : ("constant/turbulenceProperties", TurbulenceProperties),
        "boundaryPressure" : ("0/org/p_rgh", BoundaryPressure),
        "boundaryVelocity" : ("0/org/U", BoundaryVelocity),
        "boundaryPointDisplacement" : ("0/org/pointDisplacement", BoundaryPointDisplacement),
        "boundaryAlpha" : ("0/org/alpha.water", BoundaryAlpha),
        "blockMeshDict" : ("constant/polyMesh/blockMeshDict", BlockMeshDict),
        "snappyHexMeshDict" : ("constant/polyMesh/blockMeshDict", BlockMeshDict),
        "extrudeMeshDict" : ("constant/polyMesh/blockMeshDict", BlockMeshDict),
        "surfaceFeatureExtractDict" : ("constant/polyMesh/blockMeshDict", BlockMeshDict)        
              }

from copy import deepcopy
plusDefault = deepcopy(default)
plusDefault["blockMeshDict"] = ("system/blockMeshDict", BlockMeshDict)

def getFilePath( fileName , OFversion = 5) :
    
    return default[fileName][0]


def getFileClass( fileName) :
    return default[fileName][1]