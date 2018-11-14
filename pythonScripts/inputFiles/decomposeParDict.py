#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings

"""
  Convenience class to simply write DecomposeParDict
"""

class decomposeParDict(ofDictionary) :
    """
        DecomposeParDict dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('decomposeParDict'))
    def __init__(self , root, fdir, fid="decomposeParDict",**kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" :
    fv = decomposeParDict(None,"test")
    fv.update(np=42, method='metis')
    print(fv)



