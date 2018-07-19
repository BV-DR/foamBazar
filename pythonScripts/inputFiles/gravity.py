from inputFiles.ofDictionary import ofDictionary
from PyFoam.Basics.DataStructures import Dimension, Vector

"""
  Convenience class to simply write "g"
"""
class gravity(ofDictionary) :
    """
    Gravity dictionnary
    """
    def __init__(self , root, fdir, fid="g", **kwargs):
        ofDictionary.__init__(self, root,fdir,fid,**kwargs)
        if self.exists: return
        self.header["class"] = "uniformDimensionedVectorField"
        self["dimensions"] = Dimension(*[0,1,-2,0,0,0,0])
        self["value"] = Vector(*[0,0,-9.81])

        # update according to user input
        self.update(**kwargs)

if __name__ == "__main__" : 
   print(gravity("test"))



