"""
   Prepare input for wave probes
"""

from ideFoam.inputFiles.compatOF import alpha
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import DictProxy


def setWaveProbes( waveProbesList , version = "foamStar" , writeProbesInterval = 0.1 ) :
    d =DictProxy()
    d["type"] = "surfaceElevation"
    d["fields"] = "({})".format(alpha[version])
    d["writePrecision"] = 6
    d["interpolationScheme"] = "cellPointFace"
    d["outputControl"] = "timeStep"
    d["outputInterval"] = 1

    if version != "foamStar" :
        d["file"] =   "surfaceElevation.dat"
        if writeProbesInterval is not None :
            d["surfaceSampleDeltaT"] =  writeProbesInterval
    
    d["sets"] =   [ "p_{0:05} {{start ({1:} {2:} {3:}); end ({1:} {2:} {4:});  type face; axis z; nPoints {5:}; }}".format( i ,  *p )  for i, p  in enumerate(waveProbesList)  ]
    
    return d


def createLinearWaveProbesList( xMin, xMax, nX, y, zMin, Zmax, nZ) :
    waveProbesList = []
    deltaX=(xMax-xMin)/(nX-1)
    for i in range(0,nX):
        waveProbesList.append([xMin+i*deltaX, y, zMin, Zmax, nZ])
    return waveProbesList

if __name__ == "__main__" :

    # works with tuples or with lists
    #waveProbesList = ( (10.,0.,-1.,+1 , 100) , (15.,0.,-1.,+1 , 100) )
    waveProbesList = [ [10.,0.,-1.,+1 , 100] , [15.,0.,-1.,+1 , 100] ]
    print(waveProbesList)
    
    waveProbesList=createLinearWaveProbesList( -100.0, 100.0, 201, 0.05, -3.0, 3.0, 100)
    print(waveProbesList)
    d = setWaveProbes ( waveProbesList, "foamStar" , writeProbesInterval = 0.01 )
    
    waveProbFile = WriteParameterFile("waveProb.inc")
    waveProbFile["functions"] = d
    waveProbFile.writeFile()


