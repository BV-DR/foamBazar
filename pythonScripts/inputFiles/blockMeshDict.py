import PyFoam
import numpy as np
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Dimension, Vector
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

patches = '''
(
        patch back
        (
                (0 4 7 3)
        )
        wall bottom
        (
                (0 1 2 3)
        )
        patch front
        (
                (1 5 6 2)
        )
        patch atmosphere
        (
                (4 5 6 7)
        )
        patch left
        (
                (0 4 5 1)
        )
        patch right
        (
                (3 7 6 2)
        )
)
'''

class BlockMeshDict(WriteParameterFile) :
   """
      blockMeshDict dictionary
   """
   def __init__(self , case, xmin=-0.05, xmax=0.05, ymin=None, ymax=None, zmin=None, zmax=None, fsmin=None, fsmax=None, sym=False, gridlvl=1 ):
        WriteParameterFile.__init__(self,  name = join(case, "constant" , "polyMesh", "blockMeshDict" ) )
        
        self["fastMerge"]  = "yes"
        self["convertToMeters"]  = 1
        
        self["xmin"]  = xmin
        self["xmax"]  = xmax
        self["ymin"]  = ymin
        self["ymax"]  = ymax
        self["zmin"]  = zmin
        self["zmax"]  = zmax
        
        self["fsmin"]  = fsmin
        self["fsmax"]  = fsmax
        
        self["vertices"] = vertices
        
        dy = ymax - ymin
        dz = zmax - zmin
        ratio = dz/dy
        
        if sym: ny = 12
        else:   ny = 24 
        ny = int(round(ny*(np.sqrt(2)**(gridlvl-1))))
        nz = int(round(ny*ratio))
        
        print(ny,nz)
        
        self["blocks"]  = '( hex (0 1 2 3 4 5 6 7) ( 1 {:d} {:d} ) simpleGrading (1 1 1) )'.format(ny,nz)
        self["edges"]   = '( )'
        self["patches"] = patches

        self["mergePatchPairs"]   = '( )'
        
if __name__ == "__main__" : 
   print(blockMeshDict("test"))



