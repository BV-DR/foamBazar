#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings

"""
  Convenience class to simply write "g"
"""
class gravity(ofDictionary) :
    """
    Gravity dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('g'))
    def __init__(self , root, fdir, fid="g", **kwargs):
        ofDictionary.__init__(self, root,fdir,fid,**kwargs)
        self.header["class"] = "uniformDimensionedVectorField"

if __name__ == "__main__" :
    fv = gravity(None,"test")
    print(fv)



