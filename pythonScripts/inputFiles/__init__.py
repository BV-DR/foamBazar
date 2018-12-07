from inputFiles import readWriteFile
from inputFiles.fileSystem import getFilePath
ReadWriteFile = readWriteFile.ReadWriteFile
from inputFiles import gravity, refineMeshDict, transportProperties, controlDict
from inputFiles import decomposeParDict, dynamicMeshDict, turbulenceProperties
from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.dynamicMeshDict import DynamicMeshDict
from inputFiles.transportProperties import TransportProperties
from inputFiles.turbulenceProperties import TurbulenceProperties, RASProperties
from inputFiles.blockMeshDict import BlockMeshDict
from inputFiles import boundaryCondition


BoundaryPointDisplacement = boundaryCondition.BoundaryPointDisplacement
BoundaryPressure = boundaryCondition.BoundaryPressure
BoundaryVelocity = boundaryCondition.BoundaryVelocity
BoundaryAlpha = boundaryCondition.BoundaryAlpha
BoundaryK = boundaryCondition.BoundaryK
BoundaryOmega = boundaryCondition.BoundaryOmega

from inputFiles.waveProperties import RelaxZone, WaveCondition, WaveProperties, noWaves

from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.extrudeMeshDict import ExtrudeMeshDict

f = {
        "controlDict"          : ControlDict,
        "fvSchemes"            : FvSchemes,
        "fvSolution"           : FvSolution,
        "decomposeParDict"     : DecomposeParDict,
        "dynamicMeshDict"      : DynamicMeshDict,
        "transportProperties"  : TransportProperties,
        "waveProperties"       : WaveProperties,
        "rasProperties"        : RASProperties,
        "turbulenceProperties" : TurbulenceProperties,
        "boundaryPressure" : BoundaryPressure,
        "boundaryVelocity" : BoundaryVelocity,
        "boundaryPointDisplacement" : BoundaryPointDisplacement,
        "boundaryAlpha" : BoundaryAlpha,
        "blockMeshDict" : BlockMeshDict,
        "snappyHexMeshDict" : SnappyHexMeshDict,
        "extrudeMeshDict" : ExtrudeMeshDict,
        "surfaceFeatureExtractDict" : SurfaceFeatureExtractDict,
    }

def getFileClass(filename):
    return f[filename]
    
    