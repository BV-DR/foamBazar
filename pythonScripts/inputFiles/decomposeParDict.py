import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector
from os.path import join

"""
  Convenience class to simply write DecomposeParDict
"""


class DecomposeParDict(WriteParameterFile) :
   """
      DecomposeParDict dictionnary
   """
   def __init__(self , case, nProcs = 1, method = "scotch", version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "system" , "decomposeParDict" )  )
      
      self["numberOfSubdomains"] = nProcs
      self["method"] = method

if __name__ == "__main__" : 
   print  DecomposeParDict("test")



