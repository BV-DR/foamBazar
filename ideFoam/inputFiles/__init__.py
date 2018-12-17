from ideFoam.inputFiles import readWriteFile
from ideFoam.inputFiles.fileSystem import getFilePath
ReadWriteFile = readWriteFile.ReadWriteFile
from ideFoam.inputFiles import gravity
from ideFoam.inputFiles.fvSchemes import FvSchemes
from ideFoam.inputFiles.fvSolution import FvSolution
from ideFoam.inputFiles.controlDict import ControlDict
from ideFoam.inputFiles.decomposeParDict import DecomposeParDict
from ideFoam.inputFiles.dynamicMeshDict import DynamicMeshDict
from ideFoam.inputFiles.transportProperties import TransportProperties
from ideFoam.inputFiles.turbulenceProperties import TurbulenceProperties, RASProperties
from ideFoam.inputFiles.blockMeshDict import BlockMeshDict
from ideFoam.inputFiles.setSelection import SetSelection
from ideFoam.inputFiles.refineMeshDict import RefineMeshDict
from ideFoam.inputFiles.waveProperties import RelaxZone, WaveCondition, WaveProperties, noWaves
from ideFoam.inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from ideFoam.inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from ideFoam.inputFiles.extrudeMeshDict import ExtrudeMeshDict
from ideFoam.inputFiles.sixDofDomainBody import SixDofDomainBody
from ideFoam.inputFiles.flexProperties import InitFlexDict, FlexFile
from ideFoam.inputFiles.waveProbes import createLinearWaveProbesList
from ideFoam.inputFiles import boundaryCondition
BoundaryPointDisplacement = boundaryCondition.BoundaryPointDisplacement
BoundaryPressure = boundaryCondition.BoundaryPressure
BoundaryVelocity = boundaryCondition.BoundaryVelocity
BoundaryAlpha = boundaryCondition.BoundaryAlpha
BoundaryK = boundaryCondition.BoundaryK
BoundaryOmega = boundaryCondition.BoundaryOmega

f = {
        "controlDict"               : ControlDict,
        "fvSchemes"                 : FvSchemes,
        "fvSolution"                : FvSolution,
        "decomposeParDict"          : DecomposeParDict,
        "dynamicMeshDict"           : DynamicMeshDict,
        "transportProperties"       : TransportProperties,
        "waveProperties"            : WaveProperties,
        "rasProperties"             : RASProperties,
        "turbulenceProperties"      : TurbulenceProperties,
        "boundaryPressure"          : BoundaryPressure,
        "boundaryVelocity"          : BoundaryVelocity,
        "boundaryPointDisplacement" : BoundaryPointDisplacement,
        "boundaryAlpha"             : BoundaryAlpha,
        "blockMeshDict"             : BlockMeshDict,
        "snappyHexMeshDict"         : SnappyHexMeshDict,
        "extrudeMeshDict"           : ExtrudeMeshDict,
        "surfaceFeatureExtractDict" : SurfaceFeatureExtractDict,
        "sixDofDomainBody"          : SixDofDomainBody,
    }

def getFileClass(filename):
    return f[filename]
    
    