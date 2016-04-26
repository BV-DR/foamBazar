import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector,Field,Dimension
from os.path import join
from compatOF import alpha, p_rgh


"""
  Convenience class to simply write boudary condition for sea-keeping case
"""
namePatch = {
             "inlet" : "inlet" ,
             "outlet" : "outlet" ,
             "side" : "side" ,
             "inlet" : "inlet" ,
             "symmetryPlane" : "symmetryPlane" ,
             "bottom" : "bottom" ,
             "top" : "top" ,
             "structure" : "structure" ,
             "frontBack" : "frontBack" ,   #In case of  2D computation
            }

class BoundaryAlpha(WriteParameterFile) :
   """
      Alpha boundary
   """
   def __init__(self , case, symmetryPlane = "yes", namePatch = namePatch , version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" , alpha[version])  )
      self.header["class"] = "volScalarField"
      self["dimensions"] = Dimension(*[0,0,0,0,0,0,0])
      self["internalField"] = "uniform 0"
      self["boundaryField"] = { namePatch["inlet"]         : { "type" : "waveAlpha" , "value" : "uniform 0" } ,
                                namePatch["outlet"]        : { "type" : "waveAlpha" , "value" : "uniform 0" } ,
                                namePatch["bottom"]        : { "type" : "zeroGradient"  },
                                namePatch["top"]           : { "type" : "zeroGradient"  },
                                namePatch["structure"]     : { "type" : "zeroGradient"  },
                                "defaultFaces"     : { "type" : "empty"  },
                               }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }
      elif symmetryPlane == "no" :
         self["boundaryField"][namePatch["side"]] = { "type" : "waveAlpha" }



class BoundaryVelocity(WriteParameterFile) :
   """
      Velocity boundary
   """
   def __init__(self , case, speed, symmetryPlane = "yes", namePatch = namePatch, version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"U")  )
      self.header["class"] =  "volVectorField"
      self["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
      self["internalField"] = "uniform (0 0 0)"
      self["boundaryField"] = {
                                namePatch["inlet"]         : { "type" : "waveVelocity" , "value" : "uniform (0 0 0)" } ,
                                namePatch["outlet"]        : { "type" : "waveVelocity" , "value" : "uniform (0 0 0)" } ,
                                namePatch["bottom"]        : { "type" : "slip"  }         ,
                                namePatch["top"]           : { "type" : "pressureInletOutletVelocity" , "value" : "uniform (0 0 0)"  , "tangentialVelocity" : "uniform ({} 0. 0.)".format(speed) } ,
                                namePatch["structure"]     : { "type" : "slip"  }  ,
                                "defaultFaces"     : { "type" : "empty"  },
                              }
      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }
      elif symmetryPlane == "no" :
         self["boundaryField"][namePatch["side"]] = { "type" : "waveVelocity" , "value" : "uniform (0 0 0)" } ,
         



class BoundaryPressure(WriteParameterFile) :
   def __init__(self , case,  symmetryPlane = "yes", namePatch = namePatch , version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" , p_rgh[version])  )
      self.header["class"] = "volScalarField"
      self["dimensions"] = Dimension(*[1, -1, -2, 0, 0, 0, 0])
      self["internalField"] = "uniform 0"
      if version == "foamStar" :
         self["boundaryField"] = {
                                   namePatch["structure"]      : { "type" : "fixedFluxPressure" },
                                   namePatch["inlet"]          : { "type" : "fixedFluxPressure" },
                                   namePatch["outlet"]         : { "type" : "fixedFluxPressure" },
                                   namePatch["bottom"]         : { "type" : "fixedFluxPressure" },
                                   namePatch["top"]            : { "type" : "totalPressure" , "U" :  "U" , "phi" : "phi", "rho" : "rho",    "psi" : "none",   "gamma" :  0,   "p0"    :  "uniform 0",   "value" :  "uniform 0"  }  ,
                                   "defaultFaces"     : { "type" : "empty"  },
                                 }
      else :
         self["boundaryField"] = {
                                   namePatch["structure"]      : { "type" : "zeroGradient" },
                                   namePatch["inlet"]          : { "type" : "zeroGradient" },
                                   namePatch["outlet"]         : { "type" : "zeroGradient" },
                                   namePatch["bottom"]         : { "type" : "zeroGradient" },
                                   namePatch["top"]            : { "type" : "totalPressure" , "U" :  "U" , "phi" : "phi", "rho" : "rho",    "psi" : "none",   "gamma" :  0,   "p0"    :  "uniform 0",   "value" :  "uniform 0"  }  ,
                                   "defaultFaces"     : { "type" : "empty"  },
                                 }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }
      elif symmetryPlane == "no" :
         self["boundaryField"][namePatch["side"]] = { "type" : "fixedFluxPressure" , "value" : "uniform 0" } ,
      elif symmetryPlane == "2D" :
         self["boundaryField"][namePatch["frontBack"]] = { "type" : "empty" }
         
         
class BoundaryLevelSetdiff(WriteParameterFile) :
   """
      Velocity boundary
   """
   def __init__(self , case, speed, symmetryPlane = "yes", namePatch = namePatch, version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"U")  )
      self.header["class"] =  "volVectorField"
      self["dimensions"] = Dimension(*[ 0, 0, 0, 0, 0, 0, 0])
      self["internalField"] = "uniform 0"
      self["boundaryField"] = {
                                namePatch["inlet"]         : { "type" : "fixedValue" , "value" : "uniform 0" } ,
                                namePatch["outlet"]        : { "type" : "fixedValue" , "value" : "uniform 0" } ,
                                namePatch["bottom"]        : { "type" : "zeroGradient"  }         ,
                                namePatch["top"]           : { "type" : "inletOutlet" , "value" : "uniform 0"  , "inletValue" : "uniform ({} 0. 0.)".format(speed) } ,
                                namePatch["structure"]     : { "type" : "zeroGradient"  }  ,
                                "defaultFaces"     : { "type" : "empty"  },
                              }
      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }
      elif symmetryPlane == "no" :
         self["boundaryField"][namePatch["side"]] = { "type" : "fixedValue" , "value" : "uniform 0" } ,


class BoundaryUdiff(WriteParameterFile) :
   """
      Velocity boundary
   """
   def __init__(self , case, speed, symmetryPlane = "yes", namePatch = namePatch, version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"U")  )
      self.header["class"] =  "volVectorField"
      self["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
      self["internalField"] = "uniform (0 0 0)"
      self["boundaryField"] = {
                                namePatch["inlet"]         : { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
                                namePatch["outlet"]        : { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
                                namePatch["bottom"]        : { "type" : "slip"  }         ,
                                namePatch["top"]           : { "type" : "pressureInletOutletVelocity" , "value" : "uniform (0 0 0)"  , "tangentialVelocity" : "uniform ({} 0. 0.)".format(speed) } ,
                                namePatch["structure"]     : { "type" : "movingWallVelocity" , "value" : "uniform (0 0 0)"  }  ,
                                "defaultFaces"     : { "type" : "empty"  },
                              }
      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }
      elif symmetryPlane == "no" :
         self["boundaryField"][namePatch["side"]] = { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,

