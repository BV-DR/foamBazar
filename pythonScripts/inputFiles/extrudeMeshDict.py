#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
"""
  Convenience class to simply write DecomposeParDict
"""
class extrudeMeshDict(ofDictionary) :
    """
        ExtrudeMeshDict dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('extrudeMeshDict'))
    def __init__(self , root, fdir, fid="extrudeMeshDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" : 
    fv = extrudeMeshDict(None,"test")
    print(fv)



