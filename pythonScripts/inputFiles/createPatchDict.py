#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings

"""
  Convenience class to simply write DecomposeParDict
"""
class createPatchDict(ofDictionary) :
    """
        ExtrudeMeshDict dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('createPatchDict'))
    def __init__(self , root, fdir, fid="createPatchDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" : 
    fv = createPatchDict(None,"test")
    print(fv)



