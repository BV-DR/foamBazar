from inputFiles import ReadWriteFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join

"""
  Convenience class to simply write SurfaceFeatureExtractDict
"""


class SurfaceFeatureExtractDict(ReadWriteFile) :
    """
       SurfaceFeatureExtractDict dictionary
    """
    
    @classmethod
    def Build(cls , case, stlname="body"):
        res = cls( name = join(case, "system" , "surfaceFeatureExtractDict"), read = False)
        
        stlname = stlname.split('.stl')[0] #remove .stl extension
        
        body = DictProxy()
        body["extractionMethod"]  = "extractFromSurface"
        body["extractFromSurfaceCoeffs"]  = { "includedAngle" : 150. }
        body["subsetFeatures"]  = { "nonManifoldEdges" : "yes",
                                    "openEdges" : "yes" }
        body["writeObj"]  = "yes"
        res[stlname+".stl"] = body
        return res
        
if __name__ == "__main__" : 
   print(SurfaceFeatureExtractDict.Build("test"))



