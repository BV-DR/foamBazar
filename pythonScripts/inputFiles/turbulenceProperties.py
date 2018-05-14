from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from os.path import join
from compatOF import alpha, p_rgh

"""
  Convenience class to simply write "TurbulenceProperties" and RASModel
"""

RASturbulenceModel = ["kOmegaSST" , ]  #navalFoam rhoKomegaSST

class TurbulenceProperties(WriteParameterFile) :
   """
      FvSchemes dictionnary
   """
   def __init__(self , case, turbulenceModel = None ) :

      WriteParameterFile.__init__(self,  name = join(case, "constant" , "turbulenceProperties" )  )

      if turbulenceModel == "laminar" :
         self[ "simulationType" ] =  "laminar"
      elif turbulenceModel in RASturbulenceModel:
         self[ "simulationType" ] =  "RASModel"
      else :
         print("Unknown turbulence model")


class RASProperties(WriteParameterFile) :
   """
      RASProperties dictionnary
   """
   def __init__(self , case, turbulenceModel = "laminar" ) :
      WriteParameterFile.__init__(self,  name = join(case, "constant" , "RASProperties" )  )
      self["RASModel"] = turbulenceModel
      
      if turbulenceModel == "laminar" :
         self["turbulence"] = False
         self["printCoeffs"] = False
      elif turbulenceModel in RASturbulenceModel :
         self["turbulence"] = True
         self["printCoeffs"] = True
      else :
         print("Unknown turbulence model")
         
def writeTurbulenceProperties( case , turbulenceModel ) :

   R = TurbulenceProperties(case , turbulenceModel )
   R.writeFile()
   
   T = RASProperties(case , turbulenceModel )
   T.writeFile()

if __name__ == "__main__" :
#    print  TurbulenceProperties("laminar" )
#    print  RASProperties("laminar" )
   
   print(TurbulenceProperties(  ".", turbulenceModel = "kOmegaSST" ))
   print(RASProperties("." , turbulenceModel = "kOmegaSST" ))


