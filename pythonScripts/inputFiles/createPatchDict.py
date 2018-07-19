from inputFiles.ofDictionary import ofDictionary
from PyFoam.Basics.DataStructures import DictProxy

"""
  Convenience class to simply write DecomposeParDict
"""
class createPatchDict(ofDictionary) :
    """
        ExtrudeMeshDict dictionnary
    """
    def __init__(self , root, fdir, fid="createPatchDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)
        if self.exists: return

        self["pointSync"] = False
        self["patches"] = "()"

        # update according to user input
        self.update(**kwargs)

if __name__ == "__main__" : 
   print(createPatchDict("test"))



