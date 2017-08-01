import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector
from os.path import join

"""
  Convenience class to simply write "g"
"""


class Gravity(WriteParameterFile) :
   """
      Gravity dictionnary
   """
   def __init__(self , case, g = 9.81,  version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "constant" , "g" )  )
      self.header["class"] = "uniformDimensionedVectorField"
      
      self["dimensions"] = Dimension(*[0,1,-2,0,0,0,0])
      self["value"] = Vector(*[0,0,-g])

if __name__ == "__main__" : 
   print  Gravity("test")



