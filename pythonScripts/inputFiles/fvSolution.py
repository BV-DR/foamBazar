#! /usr/bin/env upython3
from inputFiles.ofDictionary import ofDictionary, ofResourcePath, getDefaultSettings
from PyFoam.Basics.DataStructures import DictProxy
from copy import deepcopy
import foamBazar, pprint
"""
  Convenience class to simply write "fvSheme"
"""

class fvSolution(ofDictionary) :
    """
        fvSchemes dictionnary
    """
    _db,_alias,_defaults = getDefaultSettings(ofResourcePath('fvSolution'))
    def __init__(self , root, fdir, fid="fvSolution", **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)
    def _update_solvers(self,data):
        if not isinstance(data, dict):
            info = 'Syntax error: ' + pprint.pformat(data)
            info += '\n\n usage: solvers={"<fieldName>": "<solverName|solverDict>"}'
            info += '\n\n solverName: smoothSolver, CMULES, GAMG, PCG, ...'
            info += '\n\n SolverDict: ... as given in OpenFOAM using DictProxy() format ... '
            raise(Exception(info))
        def addme(key,val):
            solvers=self.content['solvers']
            if isinstance(val,str) and val in self._db['solvers']:
                val=self._db['solvers'][val]
            if isinstance(val,dict):
                tmp=DictProxy() # make sure the setting is a DictProxy-class
                tmp.update(val)
                solvers[key]=deepcopy(tmp)
                return
            solvers[key] = deepcopy(val)
        for f in data.keys(): addme(f,data[f])      # do std keys first
        if not isinstance(data,DictProxy): return
        for f in data._regex[::-1]: addme(f[0],f[2]) # then do regEx
    def _update_relaxationFactors(self,data):
        '''
        usage: relax={"<fieldName>": value}'
        '''
        if not isinstance(data,dict): return
        relax=self.content["relaxationFactors"]
        allkeys = [k for k in data.keys()]
        for key in allkeys:
            if key in relax['fields']: relax['fields'][key]=deepcopy(data.pop(key,None))
            if key in relax['equations']: relax['equations'][key]=deepcopy(data.pop(key,None))
        self.content["relaxationFactors"].update(data)

#*** Main execution start here ************************************************

if __name__ == "__main__" :
    fv = fvSolution(None,'')
    fv.update(solvers={'kOmega': 'PCG'})
    fv.update(pimple={'fsiTol': 1e-3})
    fv.update(relax={'U': 0.7, 'p_rgh':0.3, 'k':0.2, 'alpha':0.1})
    fv.update(relax='none')
    fv.update(relax={'fields': { 'pd':0.1}})
    print(fv)
    print(foamBazar._OFVERSION)


