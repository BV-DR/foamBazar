from inputFiles import readWriteFile
ReadWriteFile = readWriteFile.ReadWriteFile

from inputFiles import fvSchemes, fvSolution, gravity, refineMeshDict, transportProperties, controlDict
from inputFiles import decomposeParDict, dynamicMeshDict, turbulenceProperties


FvSchemes = fvSchemes.FvSchemes
FvSolution = fvSolution.FvSolution
ControlDict = controlDict.ControlDict
DecomposeParDict = decomposeParDict.DecomposeParDict
DynamicMeshDict = dynamicMeshDict.DynamicMeshDict
TransportProperties = transportProperties.TransportProperties
TurbulenceProperties = turbulenceProperties.TurbulenceProperties
RASProperties = turbulenceProperties.RASProperties

from inputFiles import boundaryCondition

BoundaryPointDisplacement = boundaryCondition.BoundaryPointDisplacement
BoundaryPressure = boundaryCondition.BoundaryPressure
BoundaryVelocity = boundaryCondition.BoundaryVelocity
BoundaryAlpha = boundaryCondition.BoundaryAlpha
BoundaryK = boundaryCondition.BoundaryK
BoundaryOmega = boundaryCondition.BoundaryOmega

from inputFiles import waveProperties

RelaxZone = waveProperties.RelaxZone
WaveCondition = waveProperties.WaveCondition
WaveProperties = waveProperties.WaveProperties
noWaves = waveProperties.noWaves
