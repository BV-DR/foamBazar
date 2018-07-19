from inputFiles.ofDictionary import ofDictionary
from PyFoam.Basics.DataStructures import DictProxy

"""
  Convenience class to simply write DecomposeParDict
"""
class extrudeMeshDict(ofDictionary) :
    """
        ExtrudeMeshDict dictionnary
    """
    def __init__(self , root, fdir, fid="extrudeMeshDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)
        if self.exists: return

        self["constructFrom"] = "patch"
        self["sourceCase"] = "./"
        self["sourcePatches"] = "()"
        self["exposedPatchName"] = "NONE"
        self["flipNormals"] = True        
        self["extrudeModel"] = "linearNormal"
        self["nLayers"] = 1
        self["expansionRatio"] = 1.1
      
        linCoef = DictProxy()
        linCoef["thickness"]  = 0.1
        self["linearNormalCoeffs"] = linCoef
        self["mergeFaces"] = False
        self["mergeTol"] = 1e-8

        # update according to user input
        self.update(**kwargs)

if __name__ == "__main__" : 
   print(extrudeMeshDict("test"))



