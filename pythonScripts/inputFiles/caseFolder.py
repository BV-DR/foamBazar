#!/usr/bin/env upython3
#########################################################################
# Filename: caseFolder.py                                               #
# Date:     2017-July-02                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:                                                              #
#########################################################################

import os, sys, pprint, inspect
from copy import deepcopy

# this one is from foamBazar!!!
import inputFiles

#from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
#from PyFoam.Basics.DataStructures import DictProxy
#from os.path import join
#from inputFiles.compatOF import alpha, p_rgh

"""
  Convenience class to simply write "fvSheme"
"""

_of_inputFiles = None

def _isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass

def _relPath(fname, root):
    '''
    return
    fdir  : relative path to root
    fname : path to filename relative to fdir
    '''
    if (fname.startswith("./")):
        fname=fname[2:]
    if not root.endswith('/'):
        root += '/'
    if (fname.startswith(root)):
        fname=fname.replace(root,"",1)
    if (fname.startswith("./")):
        fname=fname[2:]
    if len(fname.split("/")) == 1:
        fdir='./'
    else:
        fdir=fname.split("/")[0].lower()
    if _isNumber(fdir):
        fdir="time"
    return fdir,fname

def _which_class(fname):
    global _of_inputFiles
    if _of_inputFiles is None:
        _of_inputFiles = []
        for cls in inspect.getmembers(sys.modules["inputFiles"], inspect.isclass):
            _of_inputFiles.append(cls[0])
    fid=os.path.basename(fname)
    if fid.endswith(".gz"):
        fid=fid[:-3]
    if fid is '':
        raise(Exception('not of-inputFile: '+fname))
    cls=deepcopy(fid)
    # do some name aliasing here
    if cls in ['g']: cls='gravity'
    if cls in ['U']: cls='velocityField'
    if cls in ['p_rgh','p','pd']: cls='pressureField'
    if cls in ['alpha.water','alpha1']: cls='vofField'
    if cls in _of_inputFiles:
        return cls,fid
    elif _isNumber(cls):
        return "timeDir",fid
    else:
        # unknown class, default to ofDictionary 
        return "unknown",fid

class rootDir:
    '''
    root folder
    '''
    def __init__( self, root, opts='READ_IF_PRESENT'):
        if root is None:
            root='<None>'
            opts='NO_READ'
        assert not os.path.isfile(root), "cannot overwrite existed file: " + root
        if (not root.endswith('/')):
            root += "/"
        self.root = root
        self.opts = opts
        self.files = []
        self.unknown_files = []
        if opts in ['NO_READ']:
            # FIXME: we build an empty case ... nothing to be done here
            return
        
        ok = self.readExistingCaseFolder()
        
        if opts in ['MUST_READ']:
            if not ok:
                raise(Exception('Failed to read data from: '+root))
            return
        if opts in ['READ_IF_PRESENT']:
            if not ok:
                # default files go here ...
                self.addfile('./constant/RASProperties')
                self.addfile('./constant/dynamicMeshDict')
                self.addfile('./constant/g')
                self.addfile('./constant/transportProperties')
                self.addfile('./constant/turbulenceProperties')
                self.addfile('./system/controlDict')
                self.addfile('./system/fvSchemes')
                self.addfile('./system/fvSolution')
                self.addfile('./system/decomposeParDict')
            return
        # error if we reach this line
        raise(Exception('Unknown opts:',opts))
        pass
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\nroot     : '+self.root
        info += '\nopts     : '+self.opts
        if hasattr(self, 'time'):
            info += '\ntime     : '+pprint.pformat(type(self.time))
        if hasattr(self, 'constant'):
            info += '\nconstant : '+pprint.pformat(type(self.constant))
        if hasattr(self, 'system'):
            info += '\nsystem   : '+pprint.pformat(type(self.system))
        info += '\nfiles    : ' + pprint.pformat(self.files)
        return info

    def addfile(self,fname,**kwargs):
        '''
        addfile and create corresponding obj
        '''
        fdir,fname = _relPath(fname,self.root)
        if fdir in ['time', 'constant', 'system']:
            if hasattr(self, fdir):
                eval("self."+fdir+".addfile(fname,**kwargs)")
            else:
                self.__dict__[fdir] = eval('subDir(self,"'+fdir+'",fname=fname,**kwargs)')
            return
        # FIXME: unknown type ... simply add to fileName,
        # we should create a class for these unknown file also
        #if __debug__: print('DEBUG: add unknown file:',fdir,fname)
        self.unknown_files.append(fname)
        return
    def readExistingCaseFolder(self):
        ok = True
        if not os.path.isdir(self.root):
            ok=False
            return ok
        # we read 2 levels only ... i.e. ls -d */*
        for f in os.listdir(path=self.root):
            if os.path.isdir(self.root+f):
                for s in os.listdir(path=self.root+f):
                    fname=f+'/'+s
                    if fname.endswith(".gz"):
                        fname=fname[:-3]
                    cls,fid=_which_class(fname)
                    fdir,fname=_relPath(fname, self.root)
                    if cls is not 'unknown':
                        self.addfile(fname)
                    else:
                        if __debug__: print('DEBUG: ignore:',cls,fid,fname)
                        pass
                pass
            elif os.path.isfile(self.root+f):
                # FIXME
                pass
            else:
                raise(Exception('what is this file: ',f))

        return ok

class subDir:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n_fid : '+pprint.pformat(self._fid)
        info += '\nroot : '+pprint.pformat(self._ref2Root.root)
        for key in self.__dict__.keys(): info += '\n'+self._fid+'/'+key
        return info
    def __getitem__(self, key):
        return self.__dict__[key]
    def __init__(self , root, fid, fname=None, **kwargs):
        if isinstance(root,str): root = rootDir(root)
        self._fid=fid
        self._ref2Root = root
        if fname is not None:
            self.addfile(fname, **kwargs)
        elif (os.path.isdir(root.root + self._fid)):
            #if __debug__: print('DEBUG: read known file(s) in '+self._fid+'-folder')
            pass
        pass
    def addfile(self,fname,**kwargs):
        fdir,fname=_relPath(fname,self._ref2Root.root)
        cls,fid=_which_class(fname)
        absfdir=os.path.abspath(self._ref2Root.root)
        if fdir in ["./"]:
            fdir=self._fid
            fname=os.path.join(self._fid,fname)
        if (fdir == self._fid):
            if cls in ['unknown']: cls="ofDictionary"
            if fname not in self._ref2Root.files:
                if fdir in ['time']:
                    fdir=fname.split("/")[0]
                    self.__dict__[fname] = eval('inputFiles.'+cls+'(self._ref2Root,fdir="'+absfdir+'/'+fdir+'",fid="'+fid+'",**kwargs)')
                else:
                    self.__dict__[fid] = eval('inputFiles.'+cls+'(self._ref2Root,fdir="'+absfdir+'/'+self._fid+'",fid="'+fid+'",**kwargs)')
                self._ref2Root.files.append(fname)
                #print("DEBUG: addfile "+fname+",",type(self.__dict__[fid]))
            else:
                eval("self."+fid+".update(**kwargs)")
            return
        else:
            raise(Exception('Cannot add file to '+self._fid+'-folder: ',fname))

#*** Main execution start here ************************************************
if __name__ == "__main__" :
    #print(str(sys.argv[0:]))
    #print(_of_inputFiles)
    pass
