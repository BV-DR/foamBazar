#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings

"""
  Convenience class to simply write blockMeshDict
"""

class blockMeshDict(ofDictionary) :
    """
    blockMeshDict dictionary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('blockMeshDict'))
    def __init__(self , root, fdir, fid="blockMeshDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" : 
    fv = blockMeshDict(None,"test")
    print(fv)



