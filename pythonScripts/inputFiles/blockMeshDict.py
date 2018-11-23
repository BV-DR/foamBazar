import PyFoam
import numpy as np
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector, DictProxy
from os.path import join

"""
  Convenience class to simply write blockMeshDict
"""

vertices = '''
(
        ( $xmin $ymin $zmin )
        ( $xmax $ymin $zmin )
        ( $xmax $ymax $zmin )
        ( $xmin $ymax $zmin )
        
        ( $xmin $ymin $zmax )
        ( $xmax $ymin $zmax )
        ( $xmax $ymax $zmax )
        ( $xmin $ymax $zmax )
        
        ( $xmin $ymin $fsmin )
        ( $xmax $ymin $fsmin )
        ( $xmax $ymax $fsmin )
        ( $xmin $ymax $fsmin )
            
        ( $xmin $ymin $fsmax )
        ( $xmax $ymin $fsmax )
        ( $xmax $ymax $fsmax )
        ( $xmin $ymax $fsmax )
)
'''

default_patches = '''
(
        patch domainX0
        (
                (1 5 6 2)
        )
        patch domainX1
        (
                (0 4 7 3)
        )
        patch domainY0
        (
                (0 4 5 1)
        )
        patch domainY1
        (
                (3 7 6 2)
        )
        wall domainZ0
        (
                (0 1 2 3)
        )
        patch domainZ1
        (
                (4 5 6 7)
        )
)
'''

class BlockMeshDict(WriteParameterFile) :
   """
      blockMeshDict dictionary
   """
   def __init__(self , case, ndim = 3, waveMesh=False, xmin=-0.05, xmax=0.05, ymin=None, ymax=None, zmin=None, zmax=None,
                       fsmin=None, fsmax=None, Xcells=None, Ycells=12, Zcells=None, cellRatio=1, Zgrading=None, sym=False,
                       createPatch= True, patches=None, gridlvl=1, ofp=False):
        if ofp: WriteParameterFile.__init__(self,  name = join(case, "system" , "blockMeshDict" ) )
        else:   WriteParameterFile.__init__(self,  name = join(case, "constant" , "polyMesh", "blockMeshDict" ) )
        self["fastMerge"]  = "yes"
        self["convertToMeters"]  = 1
        
        self["xmin"]  = xmin
        self["xmax"]  = xmax
        if isinstance(ymin,float): self["ymin"]  = ymin
        if isinstance(ymax,float): self["ymax"]  = ymax
        if isinstance(zmin,float): self["zmin"]  = zmin
        if isinstance(zmax,float): self["zmax"]  = zmax
        
        if isinstance(fsmin,float): self["fsmin"]  = fsmin
        if isinstance(fsmax,float): self["fsmax"]  = fsmax
        
        if waveMesh:
            verticesWave = ''
            verticesWave += '($xmin $ymin $zmin)\n    ($xmax $ymin $zmin)\n    ($xmax $ymax $zmin)\n    ($xmin $ymax $zmin)\n'
            for z in zmax:
                verticesWave += '    ($xmin $ymin {0})\n    ($xmax $ymin {0})\n    ($xmax $ymax {0})\n    ($xmin $ymax {0})\n'. format(z)
            self["vertices"] = [verticesWave]
        else:
            self["vertices"] = vertices
        
        if waveMesh:
            if len(Zcells) != len(Zgrading):
                print('ERROR: number of "Zcells" must be equal to number of "Zgrading"')
                os_exit(1)
            blockstr = ''
            for i in range(len(Zcells)):
                blockstr += 'hex ({} {} {} {} {} {} {} {}) ( {:d} {:d} {:d} ) simpleGrading (1 1 {})\n    '.format(*range(4*i,4*i+8), Xcells, Ycells, Zcells[i],Zgrading[i])
            self["blocks"]  = [ blockstr ]
        else:
            ny = Ycells*(1+(not sym))
            rXY = (xmax - xmin)/(ymax - ymin)
            rYZ = (zmax - zmin)/(ymax - ymin)
            ny = int(round(ny*(np.sqrt(2)**(gridlvl-1))))
            nz = int(round(ny*rYZ*cellRatio))
            if ndim>2: nx= int(round(ny*rXY))
            else: nx = 1
            self["blocks"]  = '( hex (0 1 2 3 4 5 6 7) ( {:d} {:d} {:d} ) simpleGrading (1 1 1) )'.format(nx, ny,nz)
        
        self["edges"]   = '()'
       
        if createPatch:
            if patches is not None: self["patches"] = patches
            else: self["patches"] = default_patches
        else:
           faces = DictProxy()
           faces["type"] = 'patch'
           faces["faces"] = '()'
           self["boundary"] = ["defaultFaces", faces]
           
        self["mergePatchPairs"]   = '()'
        
if __name__ == "__main__" : 
   print(blockMeshDict("test"))



