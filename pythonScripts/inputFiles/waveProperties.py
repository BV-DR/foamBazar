from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
from PyFoam.Basics.DataStructures import Vector
from copy import deepcopy

# from Spectral.omega2wn import omega2omegae

"""

   Routines to create the waveProperties input file (+blendingZone setSet input)
   
   TODO : Extend to irregular sea-states

"""
def _nameCoeffs(key):
    if not key.endswith('Coeffs'): key += "Coeffs"
    return key
def _nameNoCoeffs(key):
    if key.endswith('Coeffs'): key = key[:-6]
    return key
class waveProperties( ofDictionary ) :
    """
    waveProperty foamStar file
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('waveProperties'))
    def __init__(self , root, fdir, fid="waveProperties", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid, **kwargs)
    def addPatch(self,name,**kwargs):
        def addme(key):
            key = _nameCoeffs(key)
            if key not in self: self[key]=deepcopy(self._db['addPatch'])
            self[key].update(**kwargs)
        if isinstance(name,list):
            for n in name: addme(n)
        elif isinstance(name,str):
            addme(name)
        else:
            raise(Exception('invalid data:',name,kwargs))
    def addZone(self,name,**kwargs):
        def addme(key):
            key = _nameCoeffs(key)
            if key not in self: self.addPatch(key)
            if 'relaxationZone' not in self[key]: self[key]['relaxationZone'] = deepcopy(self._db['addZone'])
            self[key]['relaxationZone'].update(**kwargs)
            key = _nameNoCoeffs(key)
            if key not in self['relaxationNames']: self['relaxationNames'].append(key)
            # FIXME: check origin, orientation and etc ...
        if isinstance(name,list):
            for n in name: addme(n)
        elif isinstance(name,str):
            addme(name)
        else:
            raise(Exception('invalid data:',name,kwargs))
        pass
    def rmPatch(self,name):
        def rmMe(key):
            if _nameNoCoeffs(key) in self['relaxationNames']: return
            key = _nameCoeffs(key)
            if key in self: del self[key]
        if isinstance(name,list):
            for n in name: rmMe(n)
        elif isinstance(name,str):
            rmMe(name)
        else:
            raise(Exception('cannot remove patch:',name))
    def rmZone(self,name):
        def rmMe(key):
            if _nameNoCoeffs(key) in self['relaxationNames']: self['relaxationNames'].remove(_nameNoCoeffs(key))
            key = _nameCoeffs(key)
            if key in self and 'relaxationZone' in self[key]: del self[key]['relaxationZone']
        if isinstance(name,list):
            for n in name: rmMe(n)
        elif isinstance(name,str):
            rmMe(name)
        else:
            raise(Exception('cannot remove Zone:',name))
    def rmZonePatch(self,name):
        self.rmZone(name)
        self.rmPatch(name)

if __name__ == "__main__" :
    fv = waveProperties(None,"test")
    fv.update(defaultWave='st5')
    fv.update(defaultWave={'height':5,'period':0.3,'depth':4,'rampTime':0.1})
    fv.update(init='noRamp')
    fv.addZone('inlet',zoneName='inletCellSet',origin=Vector(1,2,3),orientation=Vector(-1,0,0))
    fv.addZone('outlet',zoneName='outletCellSet',origin=Vector(-1,-2,-3),orientation=Vector(0,-1,0))
    fv.rmZonePatch('outlet')
    print(fv)

