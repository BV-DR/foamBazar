import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector,Field,Dimension
from os.path import join
from compatOF import alpha, p_rgh , waveAlpha, waveVelocity


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

class BoundaryOmega(WriteParameterFile):

   def __init__(self,  case, symmetryPlane = "yes", wallFunction = False, version = "foamStar", namePatch = namePatch , omega = 2. ) :

      WriteParameterFile.__init__(self,  name = join(case, "0" , "org", "omega")  )
      self.header["class"] = "volScalarField"
      self["dimensions"] = Dimension(*[0,0,-1,0,0,0,0])
      self["internalField"] = "uniform {}".format(omega)
      self["boundaryField"] = { namePatch["inlet"]         : { "type" : "fixedValue", "value" : omega } ,
                                namePatch["outlet"]        : { "type" : "fixedValue", "value" : omega } ,
                                namePatch["side"]        : { "type" : "fixedValue", "value" : omega } ,
                                namePatch["bottom"]        : { "type" : "fixedValue", "value" : omega  },
                                namePatch["top"]           : { "type" : "fixedValue", "value" : omega  },
                                namePatch["structure"]     : { "type" : "zeroGradient"   },
                                "defaultFaces"     : { "type" : "empty"  },
                               }

      if wallFunction :
         self["boundaryField"][ namePatch["structure"]] = { "type" : "omegaWallFunction"   }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }


class BoundaryK(WriteParameterFile):
   def __init__(self,  case, symmetryPlane = "yes", wallFunction = False, version = "foamStar", namePatch = namePatch , k = 0.00015 ) :

      WriteParameterFile.__init__(self,  name = join(case, "0" , "org", "k")  )
      self.header["class"] = "volScalarField"
      self["dimensions"] = Dimension(*[0,2, -2,0,0,0,0])
      self["internalField"] = "uniform {}".format(k)
      self["boundaryField"] = { namePatch["inlet"]         : { "type" : "fixedValue"  , "value" : k } ,
                                namePatch["outlet"]        : { "type" : "fixedValue"  , "value" : k } ,
                                namePatch["side"]        : { "type" : "fixedValue"  , "value" : k } ,
                                namePatch["bottom"]        : { "type" : "fixedValue"  , "value" : k },
                                namePatch["top"]           : { "type" : "fixedValue"  , "value" : k },
                                namePatch["structure"]     : { "type" : "zeroGradient"   },
                                "defaultFaces"     : { "type" : "empty"  },
                               }

      if wallFunction :
         self["boundaryField"][ namePatch["structure"]] = { "type" : "kqRWallFunction"   }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }


class BoundaryAlpha(WriteParameterFile) :
   """
      Alpha boundary
   """
   def __init__(self , case, symmetryPlane = "yes", namePatch = namePatch , version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" , "org", alpha[version])  )
      self.header["class"] = "volScalarField"
      self["dimensions"] = Dimension(*[0,0,0,0,0,0,0])
      self["internalField"] = "uniform 0"
      self["boundaryField"] = { namePatch["inlet"]         : { "type" : waveAlpha[version] , "value" : "uniform 0" } ,
                                namePatch["outlet"]        : { "type" : waveAlpha[version] , "value" : "uniform 0" } ,
                                namePatch["side"]          : { "type" : waveAlpha[version] , "value" : "uniform 0" } ,
                                namePatch["bottom"]        : { "type" : "zeroGradient"  },
                                namePatch["top"]           : { "type" : "zeroGradient"  },
                                namePatch["structure"]     : { "type" : "zeroGradient"  },

                                "defaultFaces"     : { "type" : "empty"  },
                               }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }



class BoundaryVelocity(WriteParameterFile) :
   """
      Velocity boundary
   """
   def __init__(self , case, speed, symmetryPlane = "yes", namePatch = namePatch, version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"org", "U")  )
      self.header["class"] =  "volVectorField"
      self["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
      self["internalField"] = "uniform (0 0 0)"
      self["boundaryField"] = {
                                namePatch["inlet"]         : { "type" : waveVelocity[version] , "value" : "uniform (0 0 0)" } ,
                                namePatch["outlet"]        : { "type" : waveVelocity[version] , "value" : "uniform (0 0 0)" } ,
                                namePatch["side"]          : { "type" : waveVelocity[version] , "value" : "uniform (0 0 0)" } ,
                                namePatch["bottom"]        : { "type" : "slip"  }         ,
                                namePatch["top"]           : { "type" : "pressureInletOutletVelocity" , "value" : "uniform (0 0 0)"  , "tangentialVelocity" : "uniform ({} 0. 0.)".format(speed) } ,
                                namePatch["structure"]     : { "type" : "movingWallVelocity"  }  ,
                                "defaultFaces"             : { "type" : "empty"  },
                              }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }


class BoundaryPressure(WriteParameterFile) :
   def __init__(self , case,  symmetryPlane = "yes", namePatch = namePatch , version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"org",  p_rgh[version])  )
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
                                   namePatch["side"]           : { "type" : "zeroGradient" },
                                   namePatch["outlet"]         : { "type" : "zeroGradient" },
                                   namePatch["bottom"]         : { "type" : "zeroGradient" },
                                   namePatch["top"]            : { "type" : "totalPressure" , "U" :  "U" , "phi" : "phi", "rho" : "rho",    "psi" : "none",   "gamma" :  0,   "p0"    :  "uniform 0",   "value" :  "uniform 0"  }  ,
                                   "defaultFaces"     : { "type" : "empty"  },
                                 }

      if symmetryPlane == "yes" :
         self["boundaryField"][namePatch["symmetryPlane"]] = { "type" : "symmetryPlane" }
      elif symmetryPlane == "2D" :
         self["boundaryField"][namePatch["frontBack"]] = { "type" : "empty" }


class BoundaryLevelSetDiff(WriteParameterFile) :
   """
      Velocity boundary
   """
   def __init__(self , case, symmetryPlane = "yes", namePatch = namePatch, version = "swenseFoam") :
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"org", "levelSetDiff")  )
      self.header["class"] =  "volScalarField"
      self["dimensions"] = Dimension(*[ 0, 0, 0, 0, 0, 0, 0])
      self["internalField"] = "uniform 0"
      self["boundaryField"] = {
                                namePatch["inlet"]         : { "type" : "zeroGradient" } ,
                                namePatch["outlet"]        : { "type" : "zeroGradient" } ,
                                namePatch["bottom"]        : { "type" : "zeroGradient" } ,
                                namePatch["top"]           : { "type" : "zeroGradient" } ,
                                namePatch["structure"]     : { "type" : "zeroGradient" }  ,
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
      WriteParameterFile.__init__(self,  name = join(case, "0" ,"org","UDiff")  )
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

class BoundaryUinc(BoundaryUdiff):
   def __init__(self , case, speed, **kwargs) :
      BoundaryUdiff.__init__(self, case, speed, **kwargs)
      self.name = join( case, "0" , "org", "UInc" )


def writeAllBoundaries(case, version , speed = 0.0,  symmetryPlane = "yes" , namePatch = namePatch) :

   a = BoundaryAlpha( case , symmetryPlane = symmetryPlane , version = version)
   a.writeFile()

   a = BoundaryVelocity( case , speed = speed, symmetryPlane = symmetryPlane , version = version)
   a.writeFile()

   a = BoundaryPressure( case , symmetryPlane = symmetryPlane , version = version)
   a.writeFile()

   if version == "swenseFoam" :
      a = BoundaryUdiff( case , speed = speed, symmetryPlane = symmetryPlane , version = version)
      a.writeFile()
      a = BoundaryUinc( case , speed = speed, symmetryPlane = symmetryPlane , version = version)
      a.writeFile()
      a = BoundaryLevelSetDiff( case , symmetryPlane = symmetryPlane , version = version)
      a.writeFile()

if __name__ == "__main__" :

   print  BoundaryOmega("test", wallFunction = True)

