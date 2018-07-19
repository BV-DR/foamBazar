import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join

"""
  Convenience class to simply write DecomposeParDict
"""


class RefineMeshDict(WriteParameterFile) :
    """
       RefineMeshDict dictionary
    """
    def __init__(self , case, orient=None, set="freesurface", refineUptoCellLevel=None, coordinateSystem="global", directions="normal", patch=None, useHexTopology=False, geometricCut=True, writeMesh=False):
        if orient==None: suffix = ''
        else: suffix = '.'+orient
        
        WriteParameterFile.__init__(self,  name = join(case, "system" , "refineMeshDict"+suffix )  )
      
        self["set"] = set
        if refineUptoCellLevel is not None:
            self["refineUptoCellLevel"] = refineUptoCellLevel
        
        
        self["coordinateSystem"] = coordinateSystem
        
        globalCoef = DictProxy()
        globalCoef["tan1"]  = "(1 0 0)"
        globalCoef["tan2"]  = "(0 1 0)"
        self["globalCoeffs"] = globalCoef
        
        if orient is None:
            patchCoef = DictProxy()
            if patch=="outside":
                patchCoef["patch"] = patch
                patchCoef["tan1"]  = "(1 0 0)"
            else:
                patchCoef["patch"] = "patchName"
                patchCoef["tan1"]  = "(0 1 0)"
                patchCoef["tan2"]  = "(0 0 1)"
            self["patchLocalCoeffs"] = patchCoef
        
        self["directions"]  = '( {} )'.format(directions)
        
        self["useHexTopology"]  = useHexTopology
        self["geometricCut"]  = geometricCut
        self["writeMesh"]  = writeMesh
         

if __name__ == "__main__" : 
   print(RefineMeshDict("test"))



