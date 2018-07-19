from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile, ParsedParameterFile
from PyFoam.Basics.DataStructures import DictProxy, TupleProxy
from PyFoam.ThirdParty.six import PY3
from os.path import join
from copy import deepcopy
import pprint, timeit, inspect, os, inputFiles
import inputFiles.caseFolder as caseFolder

class _WriteParameterFile(WriteParameterFile):
    def __init__(self, path):
        WriteParameterFile.__init__(self,path)
    def readFile(self):
        '''
        Overload FileBasis.readFile()
        '''
        if not self.exists: return False
        super().readFile()
        return True
#        startTime = timeit.default_timer()
#        self.openFile()
#        txt=self.fh.read()
#        print("READ TXT:",timeit.default_timer()-startTime)
#        if PY3 and self.zipped:
#            txt=str(txt,"utf-8")
#        startTime = timeit.default_timer()
#        self.content=self.parse(txt)
#        print("PARSE TXT:",timeit.default_timer()-startTime)
#        self.closeFile()

class ofDictionary(_WriteParameterFile ):
    '''
    generic dictionary file
    '''
    def __init__(self , root, fdir, fid, **kwargs):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        if root is None: root = caseFolder.rootDir(None)
        _WriteParameterFile.__init__(self , join(root.root+fdir,fid))
        self._ref2Root = root
        self.listLengthUnparsed=100 # don't parse the list if it is too long
        if not self.readFile():
            if hasattr(self,'_defaults'): self.content.update(self._defaults)
        self.update(**kwargs)
    def __repr__(self):
        def contentStr(name, val,indent=0):
            out = ' '*indent + str(name) + " : "
            if not isinstance(val,dict):
                out += pprint.pformat(val) + '\n'
            else:
                out += "{\n"
                for key in val:
                    out += contentStr(key,val[key],indent=indent+4)
                if hasattr(val,'_regex'):
                    for f in val._regex:
                        out += contentStr(f[0],f[2],indent=indent+4)
                out += ' '*indent + '}\n'
            return out
        return contentStr(pprint.pformat(type(self)),self.content)
    def update(self,**kwargs):
        if hasattr(self,"_db"):
            allKeys = [key for key in kwargs.keys()]
            allFuncs = inspect.getmembers(self, predicate=inspect.ismethod)
            funcs = []
            for f in allFuncs:
                if f[0].startswith('_update_') and len(f[0])>8:
                    funcs.append(f[0][8:])
            for key in allKeys:
                keyname=key.lower()
                if hasattr(self,"_alias"):
                    if keyname in self._alias:
                        keyname=self._alias[keyname]
                    elif key in self._alias:
                        keyname=self._alias[key]
                if keyname in self._db:
                    if isinstance(kwargs[key],str) and kwargs[key] in self._db[keyname]:
                        data=kwargs.pop(key)
                        print(keyname,key,data)
                        if isinstance(self.content[keyname],dict):
                            self.content[keyname].update(self._db[keyname][data])
                        elif isinstance(self._db[keyname][data],(tuple,TupleProxy)):
                            self.content[keyname] = deepcopy(self._db[keyname][data][0])
                            for coeff in self._db[keyname][data][1]:
                                self.content[coeff] = deepcopy(self._db[keyname][data][1][coeff])
                        else:
                            self.content[keyname] = deepcopy(self._db[keyname][data])
                        if isinstance(self._db[keyname][data], dict):
                            for k in self._db[keyname][data]:
                                val = self._db[keyname][data][k]
                                if isinstance(val,str) and val.startswith('_set_'):
                                    eval("self."+val+"(keyname,k,data)")
                        if keyname in funcs: eval("self._update_"+keyname+"(data)")
                    elif isinstance(kwargs[key],(dict,str)):
                        if keyname in funcs:
                            eval("self._update_"+keyname+"(kwargs.pop(key))")
                        elif isinstance(kwargs[key],dict):
                            self.content[keyname].update(kwargs.pop(key))
                        else:
                            self.content[keyname] = kwargs.pop(key)
                    elif keyname in self.content: self.content[keyname] = kwargs.pop(key)
                else:
                    if keyname in self.content: self.content[keyname] = kwargs.pop(key)
        self.content.update(**kwargs)
        pass

class bcDictionary(ofDictionary):
    '''
    for boundaryCondition field e.g. U, p_rgh, alpha.water, ...etc
    '''
    def __repr__(self):
        info = pprint.pformat(type(self))
        for key in self.content.keys():
            info += '\n'+key+" : "
            if key in ['internalField']:
                if self.content[key].uniform:
                    info += pprint.pformat(self.content[key])
                else:
                    info += " nonuniform " + self.content[key].name
                continue
            if key in ['boundaryField']:
                tmp = deepcopy(self.content[key])
                for patch in tmp.keys():
                    if "value" in tmp[patch]:
                        if not tmp[patch]["value"].uniform:
                            tmp[patch]["value"] = 'Nonuniform ' + tmp[patch]["value"].name
                info += pprint.pformat(tmp)
                continue
            info += pprint.pformat(self.content[key])                    
        return info
    def __init__(self , root, fdir, fid, **kwargs):
        ofDictionary.__init__(self , root,fdir,fid,**kwargs)

def ofResourcePath(fname):
    return join(os.path.dirname(inputFiles.__file__), "resources", fname)
def getDefaultSettings(fname):
    db = DictProxy()
    alias = DictProxy()
    defaults = DictProxy()
    if not os.path.isfile(fname): return db,alias,defaults
    fid = ParsedParameterFile(fname,preserveComments=False,noHeader=True)
    db = deepcopy(fid['db'])
    alias = deepcopy(fid['alias'])
    defaults = deepcopy(fid['defaults'])
    return db,alias,defaults

#*** Main execution start here ************************************************      
if __name__ == "__main__" :
   print(ofDictionary("testDir",'','fname'))
