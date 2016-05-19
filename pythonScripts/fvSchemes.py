import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from os.path import join
from compatOF import alpha, p_rgh

"""
  Convenience class to simply write "TransportProperties"
"""

class FvSchemes(WriteParameterFile) :
   """
      TransportProperties dictionnary
   """
   def __init__(self , case, version = "foamStar") :

      WriteParameterFile.__init__(self,  name = join(case, "system" , "fvSchemes" )  )
      self["ddtSchemes"] = {
                             "default"         : "CrankNicolson 0.9",
                             # Incident (explicit) schemes
                             "ddt(UInc)"       : "CrankNicolson 0.9",
                             "ddt(levelSetInc)" : "CrankNicolson 0.9",
                           }

      self["gradSchemes"] = { "default" : "Gauss linear" }
      if version == "swenseFoam" :
         self["gradSchemes"]["snGradCorr(pd)"] = "interfaceGauss linear"


      self["divSchemes"] = {
                                "div(rhoPhi,U)" : "Gauss linearUpwind GradU",
                                "div(phi,alpha)" : "Gauss vanLeer",
                                "div(phirb,alpha)" : "Gauss linear",
                                "div((muEff*dev(T(grad(U)))))" : "Gauss linear",

                                "div(phi,k)"  :   "Gauss upwind",
                                "div(phi,omega)" :  "Gauss upwind",

                                #FoamExtend
                                "div(phi,U)"  :  "Gauss linearUpwind Gauss linear",

                               #SwenseFoam
                               "div(phi,UDiff)"        : "Gauss linearUpwind Gauss linear",
                               "div(phi,levelSetDiff)" : "Gauss vanLeerDC",
                               "div(phi,UInc)"         :   "Gauss linear" ,
                               "div(phi,levelSetInc)"  :  "Gauss linear"  ,
                            }


      self["laplacianSchemes"] = {
                                    "default" : "Gauss linear limited 0.5"
                                 }
      if version == "swenseFoam" :
         self["laplacianSchemes"]["laplacian(rAU,pd)"] = "interfaceGauss linear interfaceLimited 0.5"

      self["interpolationSchemes"] =  {
                                         "default"        : "linear",
                                         "interpolate(y)" : "linear"
                                      }

      self["snGradSchemes"] =  {
                                  "default" : "limited 0.5"
                               }
      if version == "swenseFoam" :
         self["snGradSchemes"]["snGrad(pd)"] = "interfaceLimited 0.5"

      self["fluxRequired"] =  {
                                   "default" :  "no",
                                   p_rgh[version] : "",
                                   "pcorr" : "",
                                   alpha[version] : "",
                               }


if __name__ == "__main__" : 
   print  FvSchemes("test" , version = "swenseFoam")




