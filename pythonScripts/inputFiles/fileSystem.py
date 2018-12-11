from copy import deepcopy

default = {
        "controlDict"               : "system/controlDict",
        "fvSchemes"                 : "system/fvSchemes",
        "fvSolution"                : "system/fvSolution",
        "decomposeParDict"          : "system/decomposeParDict",
        "snappyHexMeshDict"         : "system/snappyHexMeshDict",
        "extrudeMeshDict"           : "system/extrudeMeshDict",
        "surfaceFeatureExtractDict" : "system/surfaceFeatureExtractDict",
        "dynamicMeshDict"           : "constant/dynamicMeshDict",
        "transportProperties"       : "constant/transportProperties",
        "waveProperties"            : "constant/waveProperties",
        "RASProperties"             : "constant/RASProperties",
        "turbulenceProperties"      : "constant/turbulenceProperties", 
        "blockMeshDict"             : "constant/polyMesh/blockMeshDict",
        "boundaryPressure"          : "0/org/p_rgh", 
        "boundaryVelocity"          : "0/org/U", 
        "boundaryPointDisplacement" : "0/org/pointDisplacement", 
        "boundaryAlpha"             : "0/org/alpha.water",
        "boundaryOmega"             : "0/org/omega",
        "boundaryK"                 : "0/org/k",
        "BoundaryLevelSetDiff"      : "0/org/levelSetDiff",
        "BoundaryUdiff"             : "0/org/UDiff"
              }


#Modification of file path in plus version
plusVersion = deepcopy(default)
plusVersion["blockMeshDict"] = "system/blockMeshDict"


def getFilePath( fileName , OFversion = 5) :
    if "p" in str(OFversion).lower() :
        return plusVersion[fileName]
    else : 
        return default[fileName]
