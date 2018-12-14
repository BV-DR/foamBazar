from inputFiles import ReadWriteFile, getFilePath
from os.path import join

"""
  Convenience class to simply write "TurbulenceProperties" and RASModel
"""

RASturbulenceModel = ["kOmegaSST" , ]  #navalFoam rhoKomegaSST

class TurbulenceProperties(ReadWriteFile) :
    """
    FvSchemes dictionnary
    """
    @classmethod
    def Build(cls , case, turbulenceModel = None ) :

        res = cls(  name = join(case, getFilePath("turbulenceProperties") ), read = False )

        if turbulenceModel == "laminar" :
            res[ "simulationType" ] =  "laminar"
        elif turbulenceModel in RASturbulenceModel:
            res[ "simulationType" ] =  "RASModel"
        else :
            print("Unknown turbulence model")
        
        return res


class RASProperties(ReadWriteFile) :
    """
        RASProperties dictionnary
    """
    @classmethod
    def Build(cls , case, turbulenceModel = "laminar" ) :
    
        res = cls(  name = join(case, getFilePath("RASProperties") ), read = False )
        
        res["RASModel"] = turbulenceModel
        
        if turbulenceModel == "laminar" :
            res["turbulence"] = False
            res["printCoeffs"] = False
        elif turbulenceModel in RASturbulenceModel :
            res["turbulence"] = True
            res["printCoeffs"] = True
        else :
            print("Unknown turbulence model")
            
        return res
         
def writeTurbulenceProperties( case , turbulenceModel ) :

    R = TurbulenceProperties.Build(case , turbulenceModel )
    R.writeFile()
    
    T = RASProperties.Build(case , turbulenceModel )
    T.writeFile()

if __name__ == "__main__" :
#    print  TurbulenceProperties("laminar" )
#    print  RASProperties("laminar" )
    
    print(TurbulenceProperties.Build(  ".", turbulenceModel = "kOmegaSST" ))
    print(RASProperties.Build("." , turbulenceModel = "kOmegaSST" ))


