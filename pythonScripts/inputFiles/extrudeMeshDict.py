import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join

"""
  Convenience class to simply write DecomposeParDict
"""


class ExtrudeMeshDict(WriteParameterFile) :
   """
      ExtrudeMeshDict dictionnary
   """
   def __init__(self , case, sourceCase='"./"', exposedPatchName="front", nLayers=1, expansionRatio=1, thickness=0.1, mergeFaces=False, mergeTol=0):
        WriteParameterFile.__init__(self,  name = join(case, "system" , "extrudeMeshDict" )  )
      
        self["constructFrom"] = "patch"
        self["sourceCase"] = sourceCase
        self["sourcePatches"] = "(back)"
        
        self["exposedPatchName"] = exposedPatchName
        
        self["flipNormals"] = True
        
        self["extrudeModel"] = "linearNormal"
        self["nLayers"] = nLayers
        self["expansionRatio"] = expansionRatio
      
        linCoef = DictProxy()
        linCoef["thickness"]  = thickness
        self["linearNormalCoeffs"] = linCoef
        
        self["mergeFaces"] = mergeFaces
        self["mergeTol"] = mergeTol
         

if __name__ == "__main__" : 
   print(ExtrudeMeshDict("test"))



