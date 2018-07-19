from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
from PyFoam.Basics.DataStructures import DictProxy
"""
  Convenience class to simply write "transportProperties"
"""

class transportProperties(ofDictionary) :
    """
        TransportProperties dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('transportProperties'))
    def __init__(self , root, fdir, fid="transportProperties", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

if __name__ == "__main__" :
    fv = transportProperties(None,"test")
    fv.update(water='seaWater')
    print(fv)
    fv.update(water='basin')
    print(fv)



