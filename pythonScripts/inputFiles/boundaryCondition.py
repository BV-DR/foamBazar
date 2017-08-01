import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector, Field, Dimension, DictProxy
from os.path import join
from compatOF import alpha, p_rgh , waveAlpha, waveVelocity, pointDisp, namePatch

"""
  Convenience class to simply write boudary condition for sea-keeping case
"""

class BoundaryOmega(WriteParameterFile):

    def __init__(self,  case, symmetryPlane="yes", wallFunction=False, version="foamStar", namePatch=namePatch , omega = 2. ) :
    
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" , "org", "omega")  )
        self.header["class"] = "volScalarField"
        self["dimensions"] = Dimension(*[0,0,-1,0,0,0,0])
        self["internalField"] = "uniform {}".format(omega)
        bf = DictProxy()
        bf[patch["outlet"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["inlet"]] = { "type" : "fixedValue", "value" : omega }
        if version=="foamStar":
            bf[patch["side1"]] = { "type" : "fixedValue", "value" : omega }
            bf[patch["side2"]] = { "type" : "fixedValue", "value" : omega }
        else:
            bf[patch["side"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["bottom"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["top"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["structure"]] = { "type" : "zeroGradient"   }
        bf["defaultFaces"] = { "type" : "empty"  }
        self["boundaryField"] = bf
    
        if wallFunction :
            self["boundaryField"][ patch["structure"]] = { "type" : "omegaWallFunction"   }
    
        if symmetryPlane == "yes" :
            self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }


class BoundaryK(WriteParameterFile):
    def __init__(self,  case, symmetryPlane = "yes", wallFunction = False, version="foamStar", namePatch = namePatch , k = 0.00015 ) :
    
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" , "org", "k")  )
        self.header["class"] = "volScalarField"
        self["dimensions"] = Dimension(*[0,2, -2,0,0,0,0])
        self["internalField"] = "uniform {}".format(k)
        self["boundaryField"] = { patch["inlet"]         : { "type" : "fixedValue"  , "value" : k } ,
                                  patch["outlet"]        : { "type" : "fixedValue"  , "value" : k } ,
                                  patch["side"]        : { "type" : "fixedValue"  , "value" : k } ,
                                  patch["bottom"]        : { "type" : "fixedValue"  , "value" : k },
                                  patch["top"]           : { "type" : "fixedValue"  , "value" : k },
                                  patch["structure"]     : { "type" : "zeroGradient"   },
                                  "defaultFaces"     : { "type" : "empty"  },
                                }
    
        if wallFunction :
            self["boundaryField"][patch["structure"]] = { "type" : "kqRWallFunction"   }
    
        if symmetryPlane == "yes" :
            self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }

class BoundaryAlpha(WriteParameterFile) :
    """
        Alpha boundary
    """
    def __init__(self , case, symmetryPlane="yes", namePatch=namePatch, version="foamStar"):
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0", "org", alpha[version])  )
        self.header["class"] = "volScalarField"
        self["dimensions"] = Dimension(*[0,0,0,0,0,0,0])
        self["internalField"] = "uniform 0"
        bf = DictProxy()
        bf[patch["outlet"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
        bf[patch["inlet"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
        if version=="foamStar":
            bf[patch["side1"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
            bf[patch["side2"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
        else:
            bf[patch["side"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
        bf[patch["bottom"]] = { "type" : "zeroGradient"  }
        top=DictProxy()
        top["type"] = "inletOutlet"
        top["inletValue"] = "uniform 0"
        top["value"] = "uniform 0"
        bf[patch["top"]] = top
        bf[patch["structure"]] = { "type" : "zeroGradient"  }
        # bf["defaultFaces"] = { "type" : "empty"  }
        self["boundaryField"] = bf
    
        if symmetryPlane == "yes":
            if version=="foamStar":
                self["boundaryField"][patch["side1"]] = { "type" : "symmetryPlane" }
                self["boundaryField"][patch["side2"]] = { "type" : "symmetryPlane" }
            else:
                self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }

class BoundaryVelocity(WriteParameterFile) :
    """
        Velocity boundary
    """
    def __init__(self , case, speed, symmetryPlane = "yes", namePatch=namePatch, version = "foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org", "U")  )
        self.header["class"] =  "volVectorField"
        self["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
        self["internalField"] = "uniform ({} 0. 0.)".format(-speed)
        bf = DictProxy()
        bf[patch["outlet"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
        bf[patch["inlet"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
        if version=="foamStar":
            bf[patch["side1"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
            bf[patch["side2"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
        else:
            bf[patch["side"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
        bf[patch["bottom"]] = { "type" : "fixedValue", "value" : "uniform ({} 0. 0.)".format(-speed) }
        top=DictProxy()
        top["type"] = "pressureInletOutletVelocity"
        top["tangentialVelocity"] = "uniform ({} 0. 0.)".format(-speed)
        top["value"] = "uniform (0 0 0)"
        bf[patch["top"]] = top
        bf[patch["structure"]] = { "type" : "movingWallVelocity", "value" : "uniform (0 0 0)"  }
        # bf["defaultFaces"] = { "type" : "empty"  }
        self["boundaryField"] = bf
    
        if symmetryPlane == "yes" :
            if version=="foamStar":
                self["boundaryField"][patch["side1"]] = { "type" : "symmetryPlane" }
                self["boundaryField"][patch["side2"]] = { "type" : "symmetryPlane" }
            else:
                self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }

class BoundaryPressure(WriteParameterFile) :
    def __init__(self , case,  symmetryPlane="yes", namePatch=namePatch , version="foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org",  p_rgh[version])  )
        self.header["class"] = "volScalarField"
        self["dimensions"] = Dimension(*[1, -1, -2, 0, 0, 0, 0])
        self["internalField"] = "uniform 0"
        if version=="foamStar":
            bf = DictProxy()
            bf[patch["outlet"]] = { "type" : "fixedFluxPressure" }
            bf[patch["inlet"]] = { "type" : "fixedFluxPressure" }
            if symmetryPlane=="yes":
                bf[patch["side1"]] = { "type" : "symmetryPlane" }
                bf[patch["side2"]] = { "type" : "symmetryPlane" }
            bf[patch["bottom"]] = { "type" : "fixedFluxPressure" }
            top=DictProxy()
            top["type"] = "totalPressure"
            top["p0"] = "uniform 0"
            top["U"] = "U"
            top["phi"] = "phi"
            top["rho"] = "rho"
            top["psi"] = "none"
            top["gamma"] = 1
            top["value"] = "uniform 0"
            bf[patch["top"]] = top
            bf[patch["structure"]] = { "type" : "fixedFluxPressure"  }
            # bf["defaultFaces"] = { "type" : "empty"  }
            self["boundaryField"] = bf
        else :
            self["boundaryField"] = {
                                    patch["structure"]      : { "type" : "zeroGradient" },
                                    patch["inlet"]          : { "type" : "zeroGradient" },
                                    patch["side"]           : { "type" : "zeroGradient" },
                                    patch["outlet"]         : { "type" : "zeroGradient" },
                                    patch["bottom"]         : { "type" : "zeroGradient" },
                                    patch["top"]            : { "type" : "totalPressure" , "U" :  "U" , "phi" : "phi", "rho" : "rho",    "psi" : "none",   "gamma" :  0,   "p0"    :  "uniform 0",   "value" :  "uniform 0"  }  ,
                                    "defaultFaces"     : { "type" : "empty"  },
                                    }
    
            if symmetryPlane == "yes" :
                self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
            elif symmetryPlane == "2D" :
                self["boundaryField"][patch["frontBack"]] = { "type" : "empty" }

class BoundaryPointDisplacement(WriteParameterFile) :
    def __init__(self , case,  symmetryPlane="yes", namePatch=namePatch , version="foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org",  pointDisp[version])  )
        self.header["class"] = "pointVectorField"
        self["dimensions"] = Dimension(*[0, 1, 0, 0, 0, 0, 0])
        self["internalField"] = "uniform (0 0 0)"
        if version=="foamStar":
            bf = DictProxy()
            bf[patch["outlet"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["inlet"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            if symmetryPlane=="yes":
                bf[patch["side1"]] = { "type" : "symmetryPlane" }
                bf[patch["side2"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["bottom"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["top"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["structure"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf["InternalFaces"] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            self["boundaryField"] = bf

class BoundaryLevelSetDiff(WriteParameterFile) :
    """
        Velocity boundary
    """
    def __init__(self , case, symmetryPlane = "yes", namePatch = namePatch, version = "swenseFoam") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org", "levelSetDiff")  )
        self.header["class"] =  "volScalarField"
        self["dimensions"] = Dimension(*[ 0, 0, 0, 0, 0, 0, 0])
        self["internalField"] = "uniform 0"
        self["boundaryField"] = {
                                    patch["inlet"]         : { "type" : "zeroGradient" } ,
                                    patch["outlet"]        : { "type" : "zeroGradient" } ,
                                    patch["bottom"]        : { "type" : "zeroGradient" } ,
                                    patch["top"]           : { "type" : "zeroGradient" } ,
                                    patch["structure"]     : { "type" : "zeroGradient" }  ,
                                    "defaultFaces"     : { "type" : "empty"  },
                                }
        if symmetryPlane == "yes" :
            self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        elif symmetryPlane == "no" :
            self["boundaryField"][patch["side"]] = { "type" : "fixedValue" , "value" : "uniform 0" } ,


class BoundaryUdiff(WriteParameterFile) :
    """
        Velocity boundary
    """
    def __init__(self , case, speed, symmetryPlane = "yes", namePatch = namePatch, version = "foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org","UDiff")  )
        self.header["class"] =  "volVectorField"
        self["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
        self["internalField"] = "uniform (0 0 0)"
        self["boundaryField"] = {
                                    patch["inlet"]         : { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
                                    patch["outlet"]        : { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
                                    patch["bottom"]        : { "type" : "slip"  }         ,
                                    patch["top"]           : { "type" : "pressureInletOutletVelocity" , "value" : "uniform (0 0 0)"  , "tangentialVelocity" : "uniform ({} 0. 0.)".format(speed) } ,
                                    patch["structure"]     : { "type" : "movingWallVelocity" , "value" : "uniform (0 0 0)"  }  ,
                                    "defaultFaces"     : { "type" : "empty"  },
                                }
        if symmetryPlane == "yes" :
            self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        elif symmetryPlane == "no" :
            self["boundaryField"][patch["side"]] = { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,

class BoundaryUinc(BoundaryUdiff):
    def __init__(self , case, speed, **kwargs) :
        BoundaryUdiff.__init__(self, case, speed, **kwargs)
        self.name = join( case, "0" , "org", "UInc" )


def writeAllBoundaries(case, version, speed=0.0,  symmetryPlane="yes", namePatch=namePatch) :

    a = BoundaryAlpha( case, symmetryPlane=symmetryPlane, version=version)
    a.writeFile()
    
    a = BoundaryVelocity( case, speed=speed, symmetryPlane=symmetryPlane, version=version)
    a.writeFile()
    
    a = BoundaryPressure( case, symmetryPlane=symmetryPlane, version=version)
    a.writeFile()
    
    a = BoundaryPointDisplacement( case, symmetryPlane=symmetryPlane, version=version)
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

