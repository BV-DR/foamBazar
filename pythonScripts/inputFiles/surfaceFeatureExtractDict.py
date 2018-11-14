#! /usr/bin/env upython3
import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join

"""
  Convenience class to simply write SurfaceFeatureExtractDict
"""


class SurfaceFeatureExtractDict(WriteParameterFile) :
    """
       SurfaceFeatureExtractDict dictionary
    """
    def __init__(self , case, stlname="body"):
        WriteParameterFile.__init__(self,  name = join(case, "system" , "surfaceFeatureExtractDict"))
        
        stlname = stlname.split('.stl')[0] #remove .stl extension
        
        body = DictProxy()
        body["extractionMethod"]  = "extractFromSurface"
        body["extractFromSurfaceCoeffs"]  = { "includedAngle" : 150. }
        body["subsetFeatures"]  = { "nonManifoldEdges" : "yes",
                                    "openEdges" : "yes" }
        body["writeObj"]  = "yes"
        self[stlname+".stl"] = body
        
if __name__ == "__main__" : 
   print(SurfaceFeatureExtractDict("test"))



