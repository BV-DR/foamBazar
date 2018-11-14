#! /usr/bin/env upython3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  4 14:25:33 2018

@author: soseng
"""
from inputFiles.ofDictionary import ofDictionary
from os.path import join
import pprint

class processor:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n__dict__  : '+pprint.pformat(self.__dict__)
        return info    
    def __init__(self, fdir, fid="processor", **kwargs):
        self.name = join(fdir,fid)
        pass

