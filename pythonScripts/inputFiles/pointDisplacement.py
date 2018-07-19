from inputFiles.ofDictionary import ofDictionary, bcDictionary
from PyFoam.Basics.DataStructures import Dimension, Vector

"""
  Convenience class to simply write "g"
"""
class pointDisplacement(bcDictionary) :
    """
    Gravity dictionnary
    """
    def __init__(self , root, fdir, fid="pointDisplacement", **kwargs):
        bcDictionary.__init__(self , root,fdir,fid,**kwargs)
        if self.exists: return
        self.header['class']='pointVectorField'
        # update according to user input
        self.update(**kwargs)

if __name__ == "__main__" : 
   print(pointDisplacement("test"))



