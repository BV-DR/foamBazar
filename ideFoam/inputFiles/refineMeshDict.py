from ideFoam.inputFiles import ReadWriteFile
from PyFoam.Basics.DataStructures import DictProxy
from os.path import join

"""
  Convenience class to simply write DecomposeParDict
"""


class RefineMeshDict(ReadWriteFile) :
    """RefineMeshDict dictionary
    """
    
    @classmethod
    def Build(cls , case, orient=None, set="freesurface", refineUptoCellLevel=None, coordinateSystem="global", directions="normal", name=None, patch=None, useHexTopology=False, geometricCut=True, writeMesh=False):
        """Create RefineMeshDict from a few parameter
        """
        suffix = ''
        if orient is not None: suffix += '.'+orient
        if name is not None: suffix += '.'+name
        res = cls(join(case, "system" , "refineMeshDict"+suffix ), read = False  )
      
        res["set"] = set
        if refineUptoCellLevel is not None:
            res["refineUptoCellLevel"] = refineUptoCellLevel
        
        res["coordinateSystem"] = coordinateSystem
        
        globalCoef = DictProxy()
        globalCoef["tan1"]  = "(1 0 0)"
        globalCoef["tan2"]  = "(0 1 0)"
        res["globalCoeffs"] = globalCoef
        
        if orient is None:
            patchCoef = DictProxy()
            if patch=="outside":
                patchCoef["patch"] = patch
                patchCoef["tan1"]  = "(1 0 0)"
            else:
                patchCoef["patch"] = "patchName"
                patchCoef["tan1"]  = "(0 1 0)"
                patchCoef["tan2"]  = "(0 0 1)"
            res["patchLocalCoeffs"] = patchCoef
        
        res["directions"]  = '( {} )'.format(directions)
        
        res["useHexTopology"]  = useHexTopology
        res["geometricCut"]  = geometricCut
        res["writeMesh"]  = writeMesh
        return res

if __name__ == "__main__" : 
   print(RefineMeshDict.Build("test"))



