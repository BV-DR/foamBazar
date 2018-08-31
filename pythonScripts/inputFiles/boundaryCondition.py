import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector, Field, Dimension, DictProxy
from os.path import join
from inputFiles.compatOF import alpha, p_rgh , waveAlpha, waveVelocity, pointDisp, namePatch

"""
  Convenience class to simply write boudary condition for sea-keeping case
"""

class BoundaryOmega(WriteParameterFile):

    def __init__(self,  case, symmetry=1, wallFunction=False, version="foamStar", namePatch=namePatch , omega = 2., case2D=False ) :
    
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
    
        if wallFunction: self["boundaryField"][ patch["structure"]] = { "type" : "omegaWallFunction"   }
    
        if symmetry==1: self["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }


class BoundaryK(WriteParameterFile):
    def __init__(self,  case, symmetry=1, wallFunction = False, version="foamStar", namePatch = namePatch , k = 0.00015, case2D=False ) :
    
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
    
        if wallFunction: self["boundaryField"][patch["structure"]] = { "type" : "kqRWallFunction"   }
    
        if symmetry==1: self["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }

class BoundaryAlpha(WriteParameterFile) :
    """
        Alpha boundary
    """
    def __init__(self , case, symmetry=1, namePatch=namePatch, case2D=False, relaxZone=False, struct='', version="foamStar"):
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0", "org", alpha[version])  )
        self.header["class"] = "volScalarField"
        self["dimensions"] = Dimension(*[0,0,0,0,0,0,0])
        self["internalField"] = "uniform 0"
        bf = DictProxy()
        if case2D:
            bf[patch["outlet"]] = { "type" : "empty" }
            bf[patch["inlet"]]  = { "type" : "empty" }
            if symmetry==1: bf[patch["side1"]] = { "type" : "symmetry" }
            elif symmetry==2: bf[patch["side1"]] = { "type" : "symmetryPlane" }
            else: bf[patch["side1"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
            if relaxZone: bf[patch["side2"]] = { "type" : waveAlpha[version], "value" : "uniform 0" }
            else: bf[patch["side2"]] = { "type" : "zeroGradient" }
        else:
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
        if len(struct)>0: bf[struct] = { "type" : "zeroGradient"  }
        else: bf[patch["structure"]] = { "type" : "zeroGradient"  }
        
        self["boundaryField"] = bf
    
class BoundaryVelocity(WriteParameterFile) :
    """
        Velocity boundary
    """
    def __init__(self , case, speed, symmetry=1, namePatch=namePatch, case2D=False, relaxZone=False, struct='', version = "foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org", "U")  )
        self.header["class"] =  "volVectorField"
        self["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
        self["internalField"] = "uniform ({} 0. 0.)".format(-speed)
        bf = DictProxy()
        if case2D:
            bf[patch["outlet"]] = { "type" : "empty" }
            bf[patch["inlet"]] = { "type" : "empty" }
            if symmetry==1: bf[patch["side1"]] = { "type" : "symmetry" }
            elif symmetry==2: bf[patch["side1"]] = { "type" : "symmetryPlane" }
            else: bf[patch["side1"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
            if relaxZone: bf[patch["side2"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
            else: bf[patch["side2"]] = { "type" : "fixedValue", "value" : "uniform (0. 0. 0.)" }
        else:
            bf[patch["outlet"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
            bf[patch["inlet"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
            if version=="foamStar":
                if symmetry==1: bf[patch["side1"]] = { "type" : "symmetry" }
                elif symmetry==2: bf[patch["side1"]] = { "type" : "symmetryPlane" }
                else: bf[patch["side1"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
                bf[patch["side2"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
            else:
                if symmetry==1: bf[patch["side"]] = { "type" : "symmetry" }
                elif symmetry==2: bf[patch["side"]] = { "type" : "symmetryPlane" }
                else: bf[patch["side"]] = { "type" : waveVelocity[version], "value" : "uniform (0 0 0)" }
                
        bf[patch["bottom"]] = { "type" : "fixedValue", "value" : "uniform ({} 0. 0.)".format(-speed) }
        top=DictProxy()
        top["type"] = "pressureInletOutletVelocity"
        if speed>0: top["tangentialVelocity"] = "uniform ({} 0. 0.)".format(-speed)
        top["value"] = "uniform (0 0 0)"
        bf[patch["top"]] = top
        if len(struct)>0: bf[struct] = { "type" : "movingWallVelocity", "value" : "uniform (0 0 0)"  }
        else: bf[patch["structure"]] = { "type" : "movingWallVelocity", "value" : "uniform (0 0 0)"  }
        self["boundaryField"] = bf


class BoundaryPressure(WriteParameterFile) :
    def __init__(self , case,  symmetry=1, namePatch=namePatch, case2D=False, struct='', version="foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org",  p_rgh[version])  )
        self.header["class"] = "volScalarField"
        self["dimensions"] = Dimension(*[1, -1, -2, 0, 0, 0, 0])
        self["internalField"] = "uniform 0"
        if version=="foamStar":
            bf = DictProxy()
            if case2D:
                bf[patch["outlet"]] = { "type" : "empty" }
                bf[patch["inlet"]]  = { "type" : "empty" }
                if symmetry==1: bf[patch["side1"]]  = { "type" : "symmetry" }
                elif symmetry==2: bf[patch["side1"]]  = { "type" : "symmetryPlane" }
                else: bf[patch["side1"]]  = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
                bf[patch["side2"]]  = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
            else:
                bf[patch["outlet"]] = { "type" : "fixedFluxPressure" }
                bf[patch["inlet"]]  = { "type" : "fixedFluxPressure" }
                if symmetry==1:
                    bf[patch["side1"]] = { "type" : "symmetry" }
                    bf[patch["side2"]] = { "type" : "symmetry" }
                elif symmetry==2:
                    bf[patch["side1"]] = { "type" : "symmetryPlane" }
                    bf[patch["side2"]] = { "type" : "symmetryPlane" }
                else:
                    bf[patch["side1"]] = { "type" : "fixedFluxPressure" }
                    bf[patch["side2"]] = { "type" : "fixedFluxPressure" }
            bf[patch["bottom"]] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
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
            if len(struct)>0: bf[struct] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
            else: bf[patch["structure"]] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
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
    
            if symmetry==1: self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetry" }
            elif symmetry==2: self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
            elif symmetry=="2D": self["boundaryField"][patch["frontBack"]] = { "type" : "empty" }

class BoundaryPointDisplacement(WriteParameterFile) :
    def __init__(self , case,  symmetry=1, namePatch=namePatch, case2D=False, version="foamStar") :
        
        patch = namePatch[version]
        WriteParameterFile.__init__(self,  name = join(case, "0" ,"org",  pointDisp[version])  )
        self.header["class"] = "pointVectorField"
        self["dimensions"] = Dimension(*[0, 1, 0, 0, 0, 0, 0])
        self["internalField"] = "uniform (0 0 0)"
        if version=="foamStar":
            bf = DictProxy()
            bf[patch["outlet"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["inlet"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            if case2D:
                bf[patch["side1"]] = { "type" : "empty" }
                bf[patch["side2"]] = { "type" : "empty" }
            elif symmetry==1:
                bf[patch["side1"]] = { "type" : "symmetry" }
                bf[patch["side2"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            elif symmetry==2:
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
    def __init__(self , case, symmetry=1, namePatch = namePatch, case2D=False, version = "swenseFoam") :
        
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
        if symmetry==1: self["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        else:
            self["boundaryField"][patch["side"]] = { "type" : "fixedValue" , "value" : "uniform 0" } ,


class BoundaryUdiff(WriteParameterFile) :
    """
        Velocity boundary
    """
    def __init__(self , case, speed, symmetry=1, namePatch = namePatch, case2D=False, version = "foamStar") :
        
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
        if symmetry==1: self["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: self["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        else: self["boundaryField"][patch["side"]] = { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,

class BoundaryUinc(BoundaryUdiff):
    def __init__(self , case, speed, **kwargs) :
        BoundaryUdiff.__init__(self, case, speed, **kwargs)
        self.name = join( case, "0" , "org", "UInc" )


def writeAllBoundaries(case, version, speed=0.0,  symmetry=1, case2D=False, relaxZone=False, struct='', namePatch=namePatch) :

    a = BoundaryAlpha( case, symmetry=symmetry, case2D=case2D, relaxZone=relaxZone, struct=struct, version=version)
    a.writeFile()
    
    a = BoundaryVelocity( case, speed=speed, symmetry=symmetry, case2D=case2D, relaxZone=relaxZone, struct=struct, version=version)
    a.writeFile()
    
    a = BoundaryPressure( case, symmetry=symmetry, case2D=case2D, struct=struct, version=version)
    a.writeFile()
    
    if not case2D:
        a = BoundaryPointDisplacement( case, symmetry=symmetry, case2D=case2D, version=version)
        a.writeFile()
    
    if version == "swenseFoam" :
        a = BoundaryUdiff( case , speed = speed, symmetry=symmetry , version = version)
        a.writeFile()
        a = BoundaryUinc( case , speed = speed, symmetry=symmetry , version = version)
        a.writeFile()
        a = BoundaryLevelSetDiff( case , symmetry=symmetry , version = version)
        a.writeFile()

if __name__ == "__main__" :

   print(BoundaryOmega("test", wallFunction = True))

