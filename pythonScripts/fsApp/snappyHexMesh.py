#!/usr/bin/env upython

#########################################################################
# Filename: snappyHexMesh.py                                            #
# Date:     2018-June-27                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:                                                              #
#########################################################################

import os, sys, math, tempfile, pathlib, subprocess, warnings, pprint, csv
import numpy as np

from copy import deepcopy
from io import StringIO
from subprocess import PIPE

import fsStl
from inputFiles import caseFolder


class snappyHexMesh(caseFolder.rootDir):
    def __repr__(self): # for printing to terminal
        info = super(snappyHexMesh,self).__repr__()
        return info
    def __init__( self, root, stl=None):
        caseFolder.rootDir.__init__(self,root)
        self.addfile('system/snappyHexMeshDict')
#        if (stl is not None):
#            def addSTL(obj):
#                if isinstance(obj,fsStl.stlObjClass):
#                    self.addfile(obj)
#                elif isinstance(obj,str):
#                    self.addfile(fsStl.buildFrom(obj))
#                self.system.snappyHexMeshDict.update
#                pass
#            if isinstance(stl,list):
#                for i in stl:
#                    addSTL(i)
#            else:
#                addSTL(stl)
        pass
#    def addfile(self,f):
#        if __debug__:
#            print('DEBUG: snappyHexMesh: addfile: ',f)
#        if isinstance(f,str):
#            super(snappyHexMesh,self).addfile(f)
#        elif isinstance(f,fsStl.stlObjClass):
#            super(snappyHexMesh,self).addfile(f.fout)
#            self.stlObjs.append(f)
#            # FIXME: we need to make sure that stl is added to snappyHexMeshDict
#        pass

#class moin(WriteParameterFile):
#    def __init__(self,root):
#        if isinstance(root,str):
#            root = caseFolder.rootDir(root,empty=True)
#        WriteParameterFile.__init__(self, name=join(root.root,"system","snappyHexMeshDict"))
#        self.ref2Root=root
#        pass

#*** Main execution start here *************************************************
if __name__ == "__main__":
    print(str(sys.argv[0:]))
#    hello = moin("test")
#    print(hello.ref2Root)

