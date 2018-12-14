

from inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import Dimension, Vector
from os.path import join

"""
  Convenience class to simply write DecomposeParDict
"""

class DecomposeParDict(ReadWriteFile) :
    """
        DecomposeParDict dictionnary
    """
    
    @classmethod
    def Build(cls, case, nProcs = 1, method = "scotch", application = "foamStar") :
        
        res = cls( name = join(case, getFilePath("decomposeParDict") ), read = False )
        
        res["numberOfSubdomains"] = nProcs
        res["method"] = method
        
        if method == "simple":
            res["simpleCoeffs"] = {"n" : "( 1 3 1 )",
                                    "delta" : 0.001 }
            res["hierarchicalCoeffs"] = {"n" : "( 3 2 1 )",
                                          "delta" : 0.001,
                                          "order" : "xzy" }
            res["manualCoeffs"] = {"dataFile" : '"cellDecomposition"'}
            
        elif method == "scotch":
            res["distributed"] = "no"
            res["roots"] = '()'
            
        return res

if __name__ == "__main__" : 
   print(DecomposeParDict.Build("test"))



