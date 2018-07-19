from inputFiles.ofDictionary import ofDictionary

"""
  Convenience class to simply write blockMeshDict
"""

class blockMeshDict(ofDictionary) :
    """
    blockMeshDict dictionary
    """
    def __init__(self , root, fdir, fid="blockMeshDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)
        if self.exists: return
        self["fastMerge"]  = "yes"
        self["convertToMeters"]  = 1
        self["xmin"]  = 0
        self["xmax"]  = 0
        self["ymin"]  = 0
        self["ymax"]  = 0
        self["zmin"]  = 0
        self["zmax"]  = 0
        self["vertices"] = []       
        self["blocks"]  = []
        self["edges"]   = []
        self["patches"] = []
        self["mergePatchPairs"]   = []

        # update according to user input
        self.update(**kwargs)

if __name__ == "__main__" : 
   print(blockMeshDict("test"))



