#!/usr/bin/env python3
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

class polyMesh:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\nname      : '+pprint.pformat(self.name)
        if self.boundary: info += '\nboundary  : '+pprint.pformat(self.boundary)
        if self.faces.nFaces: info += '\nfaces     : '+pprint.pformat(self.faces)
        if self.owner.nFaces: info += '\nowner     : '+pprint.pformat(self.owner)
        if self.neighbour.nInternalFaces: info += '\nneighbour : '+pprint.pformat(self.neighbour)
        if self.points.nPoints: info += '\npoints    : '+pprint.pformat(self.points)
        #info += '\n__dict__  : '+pprint.pformat(self.__dict__)
        return info    
    def __init__(self, root, fdir, fid="polyMesh", **kwargs):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        self._ref2Root = root
        self.name = join(fdir,fid)
        self.boundary = polyMeshBoundary(root,self.name)
        self.faces = polyMeshFaces(root,self.name)
        self.neighbour = polyMeshNeighbour(root,self.name)
        self.owner = polyMeshOwner(root,self.name)
        self.points = polyMeshPoints(root,self.name)
        pass

class polyMeshBoundary(ofDictionary):
    '''
    class for polyMesh/boundary
    note: pyFoam BoundaryDict does not do what we want ...
    '''
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.content)
        return info
    def __init__(self , root, fdir, fid="boundary", **kwargs):
        ofDictionary.__init__(self , root,fdir, fid, **kwargs)
        self.header['class']='polyBoundaryMesh'
        if self.exists: return
        self.update(**kwargs)
    def readFile(self):
        self.boundaryDict=True
        super().readFile()
    def update(self,**kwargs):
        #FIXME: raise(Exception('update polyMeshBoundary ... not implemenet'))
        pass

def _sizeOfListFile(fdir,fid):
    try:
        return ListFile(fdir,fid).getSize()
    except:
        return 0

class polyMeshFaces:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.__dict__)
        return info
    def __init__(self, root, fdir, fid="faces"):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        self._ref2Root = root
        self.name = join(fdir,fid)
        self.nFaces = _sizeOfListFile(fdir,fid)

class polyMeshNeighbour:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.__dict__)
        return info
    def __init__(self, root, fdir, fid="neighbour"):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        self._ref2Root = root        
        self.name = join(fdir,fid)
        self.nInternalFaces = _sizeOfListFile(fdir,fid)

class polyMeshOwner:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.__dict__)
        return info
    def __init__(self, root, fdir, fid="owner"):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        self._ref2Root = root
        self.name = join(fdir,fid)
        self.nFaces = _sizeOfListFile(fdir,fid)

class polyMeshPoints:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.__dict__)
        return info
    def __init__(self, root, fdir, fid="points"):
        if isinstance(root,str): root = caseFolder.rootDir(root)
        self._ref2Root = root
        self.name = join(fdir,fid)
        self.nPoints = _sizeOfListFile(fdir,fid)
