"""
   Prepare input for wave probes
"""

from compatOF import surfaceElevation

def setWaveProbes( waveProbesList , version = "foamStar" ) :

   d = { surfaceElevation[version] : { "type" : "surfaceElevation" ,
                                      "outputControl" : "timeStep" ,
                                      "outputInterval" : 10
                                 }
           }

   d[surfaceElevation[version]]["sets"] =   [ "p_{0:05} {{start ({1:} {2:} {3:}); end ({1:} {2:} {4:});  type face; axis z; nPoints {5:}; }}".format( i ,  *p )  for i, p  in enumerate(waveProbesList)  ]
   
   return d