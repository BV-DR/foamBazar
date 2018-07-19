from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings

"""
  Convenience class to simply write "dynamicMeshDict"
"""

class dynamicMeshDict(ofDictionary):
    """
        dynamicMeshDict dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('dynamicMeshDict'))
    def __init__(self , root, fdir, fid="dynamicMeshDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" : 
    fv = dynamicMeshDict(None,"test")
    fv.update(mesh='static')
    fv.update(mesh='multiBody')
    print(fv)
