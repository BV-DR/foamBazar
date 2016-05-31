"""
   Prepare input for wave probes
"""

from compatOF import surfaceElevation, alpha

def setWaveProbes( waveProbesList , version = "foamStar" , writeProbesInterval = 0.1 ) :

   d = { surfaceElevation[version] : { "type" : "surfaceElevation" ,
                                      "outputControl" : "timeStep" ,
                                      "outputInterval" : 10,
                                      "interpolationScheme" : "cellPointFace"
                                 }

           }
   if version != "foamStar" :
      d[surfaceElevation[version]]["file"] =   "surfaceElevation.dat"
      d[surfaceElevation[version]]["fields"] = alpha[version]
      if writeProbesInterval is not None :
         d[surfaceElevation[version]]["surfaceSampleDeltaT"] =  writeProbesInterval
   d[surfaceElevation[version]]["sets"] =   [ "p_{0:05} {{start ({1:} {2:} {3:}); end ({1:} {2:} {4:});  type face; axis z; nPoints {5:}; }}".format( i ,  *p )  for i, p  in enumerate(waveProbesList)  ]

   return d