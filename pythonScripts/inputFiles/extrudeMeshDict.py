import PyFoam
from inputFiles import ReadWriteFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join
from inputFiles.compatOF import namePatch

"""
  Convenience class to simply write DecomposeParDict
"""


class ExtrudeMeshDict(ReadWriteFile) :
   """
      ExtrudeMeshDict dictionnary
   """
   @classmethod
   def Build(cls , case, sourceCase='"./"', exposedPatchName="inlet", nLayers=1, expansionRatio=1, thickness=0.1, mergeFaces=False, mergeTol=0, version='foamStar'):
        res = cls ( name = join(case, "system" , "extrudeMeshDict" ), read = False )
      
        patch = namePatch[version]
        
        res["constructFrom"] = "patch"
        res["sourceCase"] = sourceCase
        res["sourcePatches"] = [patch["outlet"]]
        
        res["exposedPatchName"] = patch[exposedPatchName]
        
        res["flipNormals"] = True
        
        res["extrudeModel"] = "linearNormal"
        res["nLayers"] = nLayers
        res["expansionRatio"] = expansionRatio
      
        linCoef = DictProxy()
        linCoef["thickness"]  = thickness
        res["linearNormalCoeffs"] = linCoef
        
        res["mergeFaces"] = mergeFaces
        res["mergeTol"] = mergeTol
        return res
         

if __name__ == "__main__" : 
   print(ExtrudeMeshDict.Build("test"))



