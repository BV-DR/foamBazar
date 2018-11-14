#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
from PyFoam.Basics.DataStructures import Vector
from copy import deepcopy
import numpy as np

class controlDict( ofDictionary ):
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('controlDict'))
    def __init__(self , root, fdir, fid="controlDict", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid, **kwargs)
    def funcObj(self,obj,name,**kwargs):
        db = self._db['functions']
        out = None
        if name in self['functions']:
            # modify existing settings
            out = deepcopy(self['functions'][name])
        elif obj in db:
            out = deepcopy(db[obj]['_name_'])
            if '_exe_' in out: del out['_exe_']
            # add libs to database if not done already
            if 'libs' in out:
                for l in out['libs']:
                    if not l in self['libs']: self['libs'].append(l)
                del out['libs']
        if out is None:
            info = "unknown obj: "+obj+" \nKnown obj(s):"
            for k in db.keys(): info += "\n    " + k
            for k in db._regex: info += "\n    " + k[0]
            raise(Exception(info))
        # from here we have a valid object
        if obj.lower() in ['remove','delete','rm','del']:
            del self['functions'][name]
            return out
        if '_exe_' in db[obj]['_name_']:
            eval('self.'+db[obj]['_name_']['_exe_']+'(name,out,**kwargs)')
        else:
            out.update(**kwargs)
        self['functions'].update({name:out})
        # FIXME: check if patches exist, warning only if not exists
        return out
    def setWaveProbes(self,name,cfg,**kwargs):
        def defineProbeLocations(cfg,**kwargs):
            xList=[]
            yList=[]
            zList=[]
            # get user input
            if 'xList' in kwargs: xList=deepcopy(kwargs.pop('xList'))
            if 'yList' in kwargs: yList=deepcopy(kwargs.pop('yList'))
            if 'zList' in kwargs: zList=deepcopy(kwargs.pop('zList'))
            if 'x' in kwargs: xList=deepcopy(kwargs.pop('x'))
            if 'y' in kwargs: yList=deepcopy(kwargs.pop('y'))
            if 'z' in kwargs: zList=deepcopy(kwargs.pop('z'))
            # count for direction
            d=[0,0,0]
            if not isinstance(xList,(float,int)): d[0]=len(xList)
            if not isinstance(yList,(float,int)): d[1]=len(yList)
            if not isinstance(zList,(float,int)): d[2]=len(zList)
            if isinstance(xList,(float,int)): xList=[xList]
            if isinstance(yList,(float,int)): yList=[yList]
            if isinstance(zList,(float,int)): zList=[zList]        
            # get sample axis
            axisName='xyz'
            iVal = d.index(max(d))
            iFix = d.index(min(d))
            e=deepcopy(d)
            e[iVal]=-1e6
            iAxis=e.index(max(e))
            axis=axisName[iAxis]
            # evaluate
            size=max(d)
            p0=np.zeros((size,3),float)
            p0[:,iVal]=eval(axisName[iVal]+'List')
            p0[:,iFix]=eval(axisName[iFix]+'List')[0]
            p0[:,iAxis]=min(eval(axis+'List'))
            p1=deepcopy(p0)
            p1[:,iAxis]=max(eval(axis+'List'))
            return (p0,p1,axis,size,kwargs)
        def surfaceElevation(cfg,**kwargs):
            if not any(opt in kwargs for opt in ['n','prefix','x','y','z','xList','yList','zList']): return kwargs
            n=100
            prefix='p'
            # get user input
            if 'n' in kwargs: n=kwargs.pop('n')
            if 'prefix' in kwargs: prefix=kwargs.pop('prefix')
            if not prefix.endswith('_'): prefix += '_'
            p0,p1,axis,size,kwargs=defineProbeLocations(cfg,**kwargs)
            if 'samplingParams' in cfg:
                cfg['samplingParams'].update({'axis':axis,'nPoints':n})
            else:
                cfg['samplingParams'] = {'type':'uniform','axis':axis,'nPoints':n}
            val=np.concatenate((p0,p1),axis=1)
            strformat = prefix + "{0:0"+str(len(str(size)))+"} {{"
            strformat += "start ({1:} {2:} {3:}); end ({4:} {5:} {6:}); "
            strformat += "$samplingParams;}}"
            cfg["sets"] = [strformat.format(i+1, *p) for i, p in enumerate(val)]
            return kwargs
        def interfaceHeight(cfg,**kwargs):
            if not any(opt in kwargs for opt in ['x','y','z','xList','yList','zList']): return kwargs
            p0,p1,axis,size,kwargs=defineProbeLocations(cfg,**kwargs)
            cfg["locations"] = []
            for i,p in enumerate(p0):
                avg=0.5*(p0[i]+p1[i])
                cfg["locations"].append(Vector(avg[0],avg[1],avg[2]))
            return kwargs
        if cfg['type'] in ['surfaceElevation']: kwargs = surfaceElevation(cfg,**kwargs)
        if cfg['type'] in ['interfaceHeight']: kwargs = interfaceHeight(cfg,**kwargs)
        cfg.update(**kwargs)
    def setVBM(self,name,cfg,**kwargs):
        cfg.update(**kwargs)
        if 'motionData' in cfg and cfg['motionData'] is None: del cfg['motionData']
        pass

#*** Main execution start here ************************************************
if __name__ == "__main__" :
    fv = controlDict(None,"test")
    ## add force, wave, vbm, 
    fv.funcObj('force','moin',patches=['p1','"_p2.*"'],rhoInf=1025)
    fv.funcObj('waveProbe','moin1',xList=[0.1,1,5], y=0.5, z=[-1,1],n=51)
    fv.funcObj('vbm','cuts',donFileName='"moin.don"',patches='"ship.*"')
    fv.funcObj('vbm','cuts',updateme=1)
    fv.funcObj('interfaceHeight','moin2',xList=[0.1,1,5], y=0.5, z=[-1,1])
    fv.funcObj('remove','moin')
    print(fv)

    fv.funcObj('waveProbe','moin1',xList=[0.1,1,5], y=0.5, z=[-1,1],n=51)
    fv.funcObj('waveProbe','moin1',xList=[0.1,1,5], y=0.5, z=[-5,5],n=51)
    print(fv)
