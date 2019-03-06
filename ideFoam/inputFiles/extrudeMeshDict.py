from os.path import join
from ideFoam.inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import DictProxy
from ideFoam.inputFiles.compatOF import namePatch

"""Convenient class to simply write ExtrudeMeshDict

"""

class ExtrudeMeshDict(ReadWriteFile) :
    """ExtrudeMeshDict dictionnary
    
    """
    @classmethod
    def Build(cls , case, sourceCase='"./"', sourcePatch = "outlet", exposedPatchName="inlet", nLayers=1,
                          expansionRatio=1, thickness=0.1, mergeFaces=False, mergeTol=0, version='foamStar',
                          flipNormals = True):
            
        res = cls(  name = join(case, getFilePath("extrudeMeshDict") ), read = False )
      
        patch = namePatch[version]
        
        res["constructFrom"] = "patch"
        res["sourceCase"] = sourceCase

        res["sourcePatches"] = [patch[sourcePatch]]
        res["exposedPatchName"] = patch[exposedPatchName]

        if flipNormals :
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



