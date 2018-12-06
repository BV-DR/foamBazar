from inputFiles import readWriteFile
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


from inputFiles.fileSystem import getFileClass, getFilePath