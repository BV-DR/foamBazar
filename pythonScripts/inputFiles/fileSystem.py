from copy import deepcopy

default = {
        "controlDict"          : "system/controlDict",
        "fvSchemes"            : "system/fvSchemes",
        "fvSolution"           : "system/fvSolution",
        "decomposeParDict"     : "system/decomposeParDict",
        "dynamicMeshDict"      : "constant/dynamicMeshDict",
        "transportProperties"  : "constant/transportProperties",
        "waveProperties"       : "constant/waveProperties",
        "rasProperties"        : "constant/RASProperties",
        "turbulenceProperties" : "constant/turbulenceProperties", 
        "boundaryPressure" : "0/org/p_rgh", 
        "boundaryVelocity" : "0/org/U", 
        "boundaryPointDisplacement" : "0/org/pointDisplacement", 
        "boundaryAlpha" : "0/org/alpha.water", 
        "blockMeshDict" : "constant/polyMesh/blockMeshDict",
        "snappyHexMeshDict" : "constant/polyMesh/snappyHexMeshDict",
        "extrudeMeshDict" : "constant/polyMesh/extrudeMeshDict",
        "surfaceFeatureExtractDict" : "constant/polyMesh/surfaceFeatureExtractDict"
              }


#Modification of file path in plus version
plusVersion = deepcopy(default)
plusVersion["blockMeshDict"] = "system/blockMeshDict"
plusVersion["snappyHexMeshDict"] = "system/snappyHexMeshDict"


def getFilePath( fileName , OFversion = 5) :
    if "p" in OFversion.lower() :
        return plusVersion[fileName][0]
    else : 
        return default[fileName][0]
