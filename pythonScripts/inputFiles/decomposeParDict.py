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
        
        if method == "simple":
            self["simpleCoeffs"] = {"n" : "( 1 3 1 )",
                                    "delta" : 0.001 }
            self["hierarchicalCoeffs"] = {"n" : "( 3 2 1 )",
                                          "delta" : 0.001,
                                          "order" : "xzy" }
            self["manualCoeffs"] = {"dataFile" : '"cellDecomposition"'}
            
        elif method == "scotch":
            self["distributed"] = "no"
            self["roots"] = '()'

if __name__ == "__main__" : 
   print(DecomposeParDict("test"))



