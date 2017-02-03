import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from os.path import join
from compatOF import alpha, p_rgh

"""
  Convenience class to simply write "fvSheme"
"""

class FvSchemes(WriteParameterFile) :
   """
      FvSchemes dictionnary
   """
   def __init__(self , case, version = "foamStar", prsJump = False,  orthogonalCorrection = False, blendCN = 0.9, steady = False ) :

      WriteParameterFile.__init__(self,  name = join(case, "system" , "fvSchemes" )  )


      #-------- ddtSchemes
      if not steady :
         self["ddtSchemes"] = {
                                 "default"         : "CrankNicolson {}".format(blendCN),
                              }
      else :
         self["ddtSchemes"] = {
                                 "default"         : "steadyState",
                                 "ddt(alpha)"      : "Euler",
                                 "ddt(U)"          : "Euler",
                              }
    

      #-------- gradSchemes
      self["gradSchemes"] = { "default" : "Gauss linear" }
      if prsJump :
         self["gradSchemes"]["snGradCorr(pd)"] = "interfaceGauss linear"


      #-------- divSchemes
      if steady  : bounded = "bounded"
      else : bounded = ""
      self["divSchemes"] = {
                               "div(rhoPhi,U)"  : "{} Gauss linearUpwind default".format(bounded),
                               "div(phi,U)"     : "Gauss linearUpwind default",  #FoamExtend
                               "div(phi,k)"     : "Gauss upwind",
                               "div(phi,omega)" : "Gauss upwind",
                            }

      if version == "foamStar" :
         self["divSchemes"]["div(phi,alpha)"]   = "Gauss vanLeer01"            #vanLeer01DC not in openCFD version
         self["divSchemes"]["div(phirb,alpha)"] = "Gauss interfaceCompression" #From Vuko's opinion "Gauss linear" should not be used
         self["divSchemes"]["div((muEff*dev(T(grad(U)))))"] = "Gauss linear"
      else :
         self["divSchemes"]["div(phi,alpha)"]   = "Gauss vanLeer01DC"
         self["divSchemes"]["div(phirb,alpha)"] = "Gauss vofCompression" #From Vuko's opinion "Gauss linear" should not be used
         if version == "swenseFoam" :
            self["divSchemes"]["div(phi,UDiff)"]        = "Gauss linearUpwind Gauss linear"
            self["divSchemes"]["div(phi,levelSetDiff)"] = "Gauss vanLeerDC"
            self["divSchemes"]["div(phi,UInc)"]         = "Gauss linear"
            self["divSchemes"]["div(phi,levelSetInc)"]  = "Gauss linear"


      #-------- laplacianSchemes
      if orthogonalCorrection :
         if type(orthogonalCorrection) == "float" :
             self["laplacianSchemes"] = { "default" : "Gauss linear limited {}".format(orthogonalCorrection) }
             if prsJump : 
                self["laplacianSchemes"]["laplacian(rAU,pd)"] = "interfaceGauss linear {}".format(orthogonalCorrection)
                self["laplacianSchemes"]["laplacian(rAU,p)"] = "interfaceGauss linear {}".format(orthogonalCorrection)
         elif orthogonalCorrection.lower() == "implicit" : 
             self["laplacianSchemes"] = { "default" : "Gauss linear corrected" }
         else : 
             print ("orthogonalCorrection" , orthogonalCorrection , "not reckognized")
      else :
         self["laplacianSchemes"] = { "default" : "Gauss linear uncorrected"}
         if prsJump : 
            self["laplacianSchemes"]["laplacian(rAU,pd)"] = "interfaceGauss linear interfaceUncorrected"
            self["laplacianSchemes"]["laplacian(rAU,p)"] = "interfaceGauss linear interfaceUncorrected"

      #-------- interpolationSchemes
      self["interpolationSchemes"] =  {
                                         "default"        : "linear",
                                         "interpolate(y)" : "linear"
                                      }
      #-------- snGradSchemes
      if orthogonalCorrection :
         if type(orthogonalCorrection) == "float" :
            self["snGradSchemes"] =  {  "default" : "limited {}".format(orthogonalCorrection)  }
            if prsJump : 
               self["snGradSchemes"]["snGrad(pd)"] = "interfaceLimited 0.5"
               
         elif orthogonalCorrection.lower() == "implicit" : 
            self["snGradSchemes"] =  {  "default" : "corrected"  }
      else :
         self["snGradSchemes"] =  {  "default" : "uncorrected"  }
         if prsJump : self["snGradSchemes"]["snGrad(pd)"] = "interfaceUncorrected"


      self["fluxRequired"] =  {
                                   "default" :  "no",
                                   p_rgh[version] : "",
                                   "pcorr" : "",
                                   alpha[version] : "",
                               }


if __name__ == "__main__" :
   print  FvSchemes("test" , version = "foamStar" , orthogonalCorrection = "implicit")
