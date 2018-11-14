#! /usr/bin/env upython3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 14:25:33 2018

@author: soseng
"""
from inputFiles.ofDictionary import ofDictionary
import inputFiles.caseFolder as caseFolder
from os.path import join
import pprint

from PyFoam.RunDictionary.ListFile import ListFile

class modeshapes:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\nname      : '+pprint.pformat(self.name)
        #info += '\n__dict__  : '+pprint.pformat(self.__dict__)
        return info    
    def __init__(self, root, fdir, fid="modeshapes", **kwargs):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        self._ref2Root = root        
        self.name = join(fdir,fid)
        pass

