from PyFoam.Basics.DataStructures import Vector, Field, Dimension, DictProxy
from os.path import join
from ideFoam.inputFiles.compatOF import alpha, p_rgh , waveAlpha, waveVelocity, pointDisp, namePatch
from ideFoam.inputFiles import ReadWriteFile, getFilePath

"""
  Convenience class to simply write boudary condition for sea-keeping case
"""

class BoundaryOmega(ReadWriteFile):

    @classmethod
    def Build(cls,  case, symmetry=1, wallFunction=False, application="foamStar", namePatch=namePatch , omega = 2., case2D=False ) :
        
        patch = namePatch[application]
        res = cls( name = join(case, getFilePath("boundaryOmega") ), read = False )
        res.header["class"] = "volScalarField"
        res["dimensions"] = Dimension(*[0,0,-1,0,0,0,0])
        res["internalField"] = "uniform {}".format(omega)
        bf = DictProxy()
        bf[patch["outlet"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["inlet"]] = { "type" : "fixedValue", "value" : omega }
        if application=="foamStar":
            bf[patch["side1"]] = { "type" : "fixedValue", "value" : omega }
            bf[patch["side2"]] = { "type" : "fixedValue", "value" : omega }
        else:
            bf[patch["side"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["bottom"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["top"]] = { "type" : "fixedValue", "value" : omega }
        bf[patch["structure"]] = { "type" : "zeroGradient"   }
        bf["defaultFaces"] = { "type" : "empty"  }
        res["boundaryField"] = bf
    
        if wallFunction: res["boundaryField"][ patch["structure"]] = { "type" : "omegaWallFunction"   }
    
        if symmetry==1: res["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: res["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        return res


class BoundaryK(ReadWriteFile):
    
    @classmethod
    def Build(cls,  case, symmetry=1, wallFunction = False, application="foamStar", namePatch = namePatch , k = 0.00015, case2D=False ) :
        
        patch = namePatch[application]
        res = cls( name = join(case, getFilePath("boundaryK") ), read = False )
        res.header["class"] = "volScalarField"
        res["dimensions"] = Dimension(*[0,2, -2,0,0,0,0])
        res["internalField"] = "uniform {}".format(k)
        res["boundaryField"] = { patch["inlet"]         : { "type" : "fixedValue"  , "value" : k } ,
                                  patch["outlet"]        : { "type" : "fixedValue"  , "value" : k } ,
                                  patch["side"]        : { "type" : "fixedValue"  , "value" : k } ,
                                  patch["bottom"]        : { "type" : "fixedValue"  , "value" : k },
                                  patch["top"]           : { "type" : "fixedValue"  , "value" : k },
                                  patch["structure"]     : { "type" : "zeroGradient"   },
                                  "defaultFaces"     : { "type" : "empty"  },
                                }
    
        if wallFunction: res["boundaryField"][patch["structure"]] = { "type" : "kqRWallFunction"   }
    
        if symmetry==1: res["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: res["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        return res

class BoundaryAlpha(ReadWriteFile) :
    """
        Alpha boundary
    """
    @classmethod
    def Build(cls , case, symmetry=1, namePatch=namePatch, case2D=False, wave=True, relaxZone=False, struct='', application="foamStar"):
        
        patch = namePatch[application]
        if application=="foamStar": res = cls( name = join(case, getFilePath("boundaryAlpha") ), read = False )
        else: res = cls( name = join(case, "0" , "org", alpha[application]) , read = False )
        res.header["class"] = "volScalarField"
        res["dimensions"] = Dimension(*[0,0,0,0,0,0,0])
        res["internalField"] = "uniform 0"
        
        bf = DictProxy()
        if wave: wavePatch = { "type" : waveAlpha[application], "value" : "uniform 0" }
        else: wavePatch = { "type" : "zeroGradient" }
        
        if case2D:
            bf[patch["outlet"]] = { "type" : "empty" }
            bf[patch["inlet"]]  = { "type" : "empty" }
        else:
            bf[patch["outlet"]] = wavePatch
            bf[patch["inlet"]] = wavePatch
        
        if symmetry==1: bf[patch["side1"]] = { "type" : "symmetry" }
        elif symmetry==2: bf[patch["side1"]] = { "type" : "symmetryPlane" }
        else: bf[patch["side1"]] = wavePatch
        if relaxZone: bf[patch["side2"]] = wavePatch
        else: bf[patch["side2"]] = { "type" : "zeroGradient" }
        
        bf[patch["bottom"]] = { "type" : "zeroGradient"  }
        top=DictProxy()
        top["type"] = "inletOutlet"
        top["inletValue"] = "uniform 0"
        top["value"] = "uniform 0"
        bf[patch["top"]] = top
        if len(struct)>0: bf[struct] = { "type" : "zeroGradient"  }
        else: bf[patch["structure"]] = { "type" : "zeroGradient"  }
        
        res["boundaryField"] = bf
        return res
    
class BoundaryVelocity(ReadWriteFile) :
    """
        Velocity boundary
    """
    @classmethod
    def Build(cls , case, speed=0., symmetry=1, namePatch=namePatch, case2D=False, wave = True, relaxZone=False, struct='', application = "foamStar") :
        
        patch = namePatch[application]
        res = cls( name = join(case, getFilePath("boundaryVelocity") ), read = False )
        res.header["class"] =  "volVectorField"
        res["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
        res["internalField"] = "uniform ({} 0. 0.)".format(-speed)
        bf = DictProxy()
        
        if wave: wavePatch = { "type" : waveVelocity[application], "value" : "uniform (0 0 0)" }
        else: wavePatch = { "type" : "fixedValue", "value" : "uniform (0. 0. 0.)" }
        
        if case2D:
            bf[patch["outlet"]] = { "type" : "empty" }
            bf[patch["inlet"]] = { "type" : "empty" }
        else:
            bf[patch["outlet"]] = wavePatch
            bf[patch["inlet"]] = wavePatch
        
        if symmetry==1: bf[patch["side1"]] = { "type" : "symmetry" }
        elif symmetry==2: bf[patch["side1"]] = { "type" : "symmetryPlane" }
        else: bf[patch["side1"]] = wavePatch
        
        if relaxZone: bf[patch["side2"]] = wavePatch
        else: bf[patch["side2"]] = { "type" : "fixedValue", "value" : "uniform (0. 0. 0.)" }

        bf[patch["bottom"]] = { "type" : "fixedValue", "value" : "uniform ({} 0. 0.)".format(-speed) }
        
        top=DictProxy()
        top["type"] = "pressureInletOutletVelocity"
        if speed > 0 :
            top["tangentialVelocity"] = "uniform ({} 0. 0.)".format(-speed)
        top["value"] = "uniform (0 0 0)"
        bf[patch["top"]] = top
        
        
        if len(struct)>0: bf[struct] = { "type" : "movingWallVelocity", "value" : "uniform (0 0 0)"  }
        else: bf[patch["structure"]] = { "type" : "movingWallVelocity", "value" : "uniform (0 0 0)"  }
        res["boundaryField"] = bf
        return res


class BoundaryPressure(ReadWriteFile) :
    @classmethod
    def Build(cls , case,  symmetry=1, namePatch=namePatch, case2D=False, struct='', application="foamStar") :
        
        patch = namePatch[application]
        if application=="foamStar":
            res = cls( name = join(case, getFilePath("boundaryPressure") ), read = False )
        else:
            res = cls(name = join(case, "0" ,"org",  p_rgh[application]) , read = False )
        res.header["class"] = "volScalarField"
        res["dimensions"] = Dimension(*[1, -1, -2, 0, 0, 0, 0])
        res["internalField"] = "uniform 0"
        if application=="foamStar":
            bf = DictProxy()
            if case2D:
                bf[patch["outlet"]] = { "type" : "empty" }
                bf[patch["inlet"]]  = { "type" : "empty" }
            else:
                bf[patch["outlet"]] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
                bf[patch["inlet"]]  = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
                
            if symmetry==1: bf[patch["side1"]]  = { "type" : "symmetry" }
            elif symmetry==2: bf[patch["side1"]]  = { "type" : "symmetryPlane" }
            else: bf[patch["side1"]]  = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
            bf[patch["side2"]]  = { "type" : "fixedFluxPressure", "value" : "uniform 0" }

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
            res["boundaryField"] = bf
        else :
            res["boundaryField"] = {
                                    patch["structure"]      : { "type" : "zeroGradient" },
                                    patch["inlet"]          : { "type" : "zeroGradient" },
                                    patch["side"]           : { "type" : "zeroGradient" },
                                    patch["outlet"]         : { "type" : "zeroGradient" },
                                    patch["bottom"]         : { "type" : "zeroGradient" },
                                    patch["top"]            : { "type" : "totalPressure" , "U" :  "U" , "phi" : "phi", "rho" : "rho",    "psi" : "none",   "gamma" :  0,   "p0"    :  "uniform 0",   "value" :  "uniform 0"  }  ,
                                    "defaultFaces"     : { "type" : "empty"  },
                                    }
    
            if symmetry==1: res["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetry" }
            elif symmetry==2: res["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
            elif symmetry=="2D": res["boundaryField"][patch["frontBack"]] = { "type" : "empty" }
        return res

class BoundaryPointDisplacement(ReadWriteFile) :
    @classmethod
    def Build(cls , case,  symmetry=1, namePatch=namePatch, case2D=False, struct='', cpMorphing = False, application="foamStar") :
        
        patch = namePatch[application]
        if application=="foamStar":
            res = cls( name = join(case, getFilePath("boundaryPointDisplacement") ), read = False )
        else:
            res = cls(name = join(case, "0" ,"org",  pointDisp[application]) , read = False )
        res.header["class"] = "pointVectorField"
        res["dimensions"] = Dimension(*[0, 1, 0, 0, 0, 0, 0])
        res["internalField"] = "uniform (0 0 0)"
        if application=="foamStar":
            bf = DictProxy()
            if case2D:
                bf[patch["outlet"]] = { "type" : "empty" }
                bf[patch["inlet"]] = { "type" : "empty" }
            else:
                bf[patch["outlet"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
                bf[patch["inlet"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
                
            if symmetry==1: bf[patch["side1"]] = { "type" : "symmetry" }
            elif symmetry==2: bf[patch["side1"]] = { "type" : "symmetryPlane" }
            else: bf[patch["side1"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["side2"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["bottom"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            bf[patch["top"]] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            
            if len(struct)>0: bf[struct] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
            else: bf[patch["structure"]] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
            
            if cpMorphing :
                bf[patch["structure"]] = { "type" : "calculated" }
            else :
                if len(struct)>0: bf[struct] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
                else: bf[patch["structure"]] = { "type" : "fixedFluxPressure", "value" : "uniform 0" }
            
            bf["InternalFaces"] = { "type" : "fixedValue", "value" : "uniform (0 0 0)" }
            res["boundaryField"] = bf
        return res

class BoundaryLevelSetDiff(ReadWriteFile) :
    """
        Velocity boundary
    """
    @classmethod
    def Build(cls ,case, symmetry=1, namePatch = namePatch, case2D=False, application = "swenseFoam") :
        
        patch = namePatch[application]
        res = cls( name = join(case, getFilePath("BoundaryLevelSetDiff") ), read = False )
        res.header["class"] =  "volScalarField"
        res["dimensions"] = Dimension(*[ 0, 0, 0, 0, 0, 0, 0])
        res["internalField"] = "uniform 0"
        res["boundaryField"] = {
                                    patch["inlet"]         : { "type" : "zeroGradient" } ,
                                    patch["outlet"]        : { "type" : "zeroGradient" } ,
                                    patch["bottom"]        : { "type" : "zeroGradient" } ,
                                    patch["top"]           : { "type" : "zeroGradient" } ,
                                    patch["structure"]     : { "type" : "zeroGradient" }  ,
                                    "defaultFaces"     : { "type" : "empty"  },
                                }
        if symmetry==1: res["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: res["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        else:
            res["boundaryField"][patch["side"]] = { "type" : "fixedValue" , "value" : "uniform 0" } ,
        return res

class BoundaryUdiff(ReadWriteFile) :
    """
        Velocity boundary
    """
    @classmethod
    def Build(cls ,case, speed, symmetry=1, namePatch = namePatch, case2D=False, application = "foamStar") :
        patch = namePatch[application]
        res = cls( name = join(case, getFilePath("BoundaryUdiff") ), read = False )
        res.header["class"] =  "volVectorField"
        res["dimensions"] = Dimension(*[ 0, 1, -1, 0, 0, 0, 0])
        res["internalField"] = "uniform (0 0 0)"
        res["boundaryField"] = {
                                    patch["inlet"]         : { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
                                    patch["outlet"]        : { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
                                    patch["bottom"]        : { "type" : "slip"  }         ,
                                    patch["top"]           : { "type" : "pressureInletOutletVelocity" , "value" : "uniform (0 0 0)"  , "tangentialVelocity" : "uniform ({} 0. 0.)".format(speed) } ,
                                    patch["structure"]     : { "type" : "movingWallVelocity" , "value" : "uniform (0 0 0)"  }  ,
                                    "defaultFaces"     : { "type" : "empty"  },
                                }
        if symmetry==1: res["boundaryField"][patch["symmetry"]] = { "type" : "symmetry" }
        elif symmetry==2: res["boundaryField"][patch["symmetryPlane"]] = { "type" : "symmetryPlane" }
        else: res["boundaryField"][patch["side"]] = { "type" : "inletOutlet" , "value" : "uniform (0 0 0)" , "inletValue" : "uniform (0 0 0)"  } ,
        return res

class BoundaryUinc(BoundaryUdiff):
    @classmethod
    def Build(cls ,case, speed, **kwargs) :
        res = cls( case, speed, **kwargs)
        res.name = join( case, "0" , "org", "UInc" )
        return res


def writeAllBoundaries(case, application, speed=0.0,  symmetry=1, case2D=False, wave=True, relaxZone=False, struct='', namePatch=namePatch) :

    a = BoundaryAlpha.Build( case, symmetry=symmetry, case2D=case2D, wave=wave, relaxZone=relaxZone, struct=struct, application=application)
    a.writeFile()
    
    a = BoundaryVelocity.Build( case, speed=speed, symmetry=symmetry, case2D=case2D, wave=wave, relaxZone=relaxZone, struct=struct, application=application)
    a.writeFile()
    
    a = BoundaryPressure.Build( case, symmetry=symmetry, case2D=case2D, struct=struct, application=application)
    a.writeFile()
    
    # if not case2D:
        # a = BoundaryPointDisplacement( case, symmetry=symmetry, case2D=case2D, application=application)
        # a.writeFile()
    
    if application == "swenseFoam" :
        a = BoundaryUdiff.Build( case , speed = speed, symmetry=symmetry , application = application)
        a.writeFile()
        a = BoundaryUinc.Build( case , speed = speed, symmetry=symmetry , application = application)
        a.writeFile()
        a = BoundaryLevelSetDiff.Build( case , symmetry=symmetry , application = application)
        a.writeFile()

if __name__ == "__main__" :

   print(BoundaryOmega.Build("test", wallFunction = True))

