from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
from PyFoam.Basics.DataStructures import TupleProxy
from copy import deepcopy
"""
  Convenience class to simply write "fvSheme"
"""

class fvSchemes(ofDictionary) :
    """
        fvSchemes dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('fvSchemes'))
    def __init__(self , root, fdir, fid="fvSchemes", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid, **kwargs)
    def _set_CN(self,keyname,key,data):
        if not isinstance(data,str): return
        if not data.startswith('CN'): return
        data=data[2:]
        val = TupleProxy()
        val.append('CrankNicolson')
        if "_ramp" in data:
            coeff,ramp=data.split("_ramp")
            val.append('ocCoeff');
            val.append({'type':'scale','scale':'linearRamp','duration':ramp,'value':coeff})
        elif "_blend_" in data:
            coeff,cellSet=data.split("_blend_")
            if key == 'default':
                val.append(coeff)
            else:
                val = TupleProxy(('BlendingCNEuler', coeff, cellSet))
        else:
            val.append(data)
        self.content[keyname][key] = val
        pass

if __name__ == "__main__" :
    fv = fvSchemes(None,'')
#    fv.update(ddt='CN0.9_ramp10')
#    fv.update(ddt='CN0.9_blend_EulerCell')
#    fv.update(ddt='CN0.9')
#    fv.update(ddt='Euler')
#    print(fv.__repr__())
#    fv.update(ddtScheme='CNEuler0.9')
#    fv.update(snGrad={'limitednn':1})
#    fv.update(flux='foamExtend')
#    fv.update(ddt='quickvvv')
#    fv.update(div='quick_tttt')
    print(fv)
    
    
    