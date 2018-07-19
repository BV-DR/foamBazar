from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
"""
    Convenience class to simply write "turbulenceProperties" and RASModel
"""

class turbulenceProperties(ofDictionary) :
    """
        turbulenceProperties dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('turbulenceProperties'))
    def __init__(self , root, fdir, fid="turbulenceProperties", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" :
    fv = turbulenceProperties(None,"test")
    fv.update(type='none')
    fv.update(type='laminar')
    fv.update(model='les')
    fv.update(model='ras')
    fv.update(model='kOmegaSST')
    print(fv)


    