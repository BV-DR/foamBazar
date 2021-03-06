from os.path import join
import numpy as np
from ideFoam.inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import DictProxy


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

class BlockMeshDict(ReadWriteFile) :
    """
      blockMeshDict dictionary
    """

    @classmethod
    def Build(cls , case, ndim = 3, waveMesh=False, xmin=-0.05, xmax=0.05, ymin=None, ymax=None, zmin=None, zmax=None,
                       fsmin=None, fsmax=None, Xcells=None, Ycells=12, Zcells=None, cellRatio=1, Zgrading=None, sym=False,
                       createPatch= True, patches=None, gridlvl=1, OFversion=5):
        """Build blockMeshDict file from a few parameters.

        Parameters
        ----------
        case : str
            Name of case

        ndim : int, default 3
            Number of dimentions of mesh (2 for 2D mesh or 3 for 3D mesh)
        sym : bool, default False
            Logical defining if symmetric mesh should be created

        waveMesh : bool, default False
            Logical defining if mesh should be defined according to wave Mesh good practices

        xmin : float, default -0.05
            Minimum position of X domain boundary
        xmax : float, default 0.05
            Maximum position of X domain boundary
        ymin : float
            Minimum position of Y domain boundary
        ymax : float
            Maximum position of Y domain boundary
        zmin : float
            Minimum position of Z domain boundary
        zmax : float
            Maximum position of Z domain boundary
        fsmin : float
            Minimum position of free surface zone
        fsmax : float
            Maximum position of free surface zone

        Xcells : int
            Number of cells in X direction
        Ycells : int, default 12
            Number of cells in Y direction
        Zcells : int
            Number of cells in Z direction
        Zgrading :
            TODO
        cellRatio : int, default 1
            If not waveMesh, Y/Z cell ratio

        createPatch : bool, default True
            Logical defining is patches should be created by blockMesh
        patches : list
            If createPatch is True, user can define a custom patch definition (refer to blockMesh help

        gridlvl : int, default 1
            TODO
        OFversion: int or str, default 5
            OpenFOAM version used

        """

        res = cls( name = join(case, getFilePath("blockMeshDict", OFversion) ), read = False )
        res["fastMerge"]  = "yes"
        res["convertToMeters"]  = 1

        res["xmin"]  = xmin
        res["xmax"]  = xmax
        if isinstance(ymin,float): res["ymin"]  = ymin
        if isinstance(ymax,float): res["ymax"]  = ymax
        if isinstance(zmin,float): res["zmin"]  = zmin
        if isinstance(zmax,float): res["zmax"]  = zmax

        if isinstance(fsmin,float): res["fsmin"]  = fsmin
        if isinstance(fsmax,float): res["fsmax"]  = fsmax

        if waveMesh:
            verticesWave = ''
            verticesWave += '($xmin $ymin $zmin)\n    ($xmax $ymin $zmin)\n    ($xmax $ymax $zmin)\n    ($xmin $ymax $zmin)\n'
            for z in zmax:
                verticesWave += '    ($xmin $ymin {0})\n    ($xmax $ymin {0})\n    ($xmax $ymax {0})\n    ($xmin $ymax {0})\n'. format(z)
            res["vertices"] = [verticesWave]
        else:
            res["vertices"] = vertices

        if waveMesh:
            if len(Zcells) != len(Zgrading):
                print('ERROR: number of "Zcells" must be equal to number of "Zgrading"')
                raise(Exception)
            blockstr = ''
            for i in range(len(Zcells)):
                blockstr += 'hex ({} {} {} {} {} {} {} {}) ( {:d} {:d} {:d} ) simpleGrading (1 1 {})\n    '.format(*range(4*i,4*i+8), Xcells, Ycells, Zcells[i],Zgrading[i])
            res["blocks"]  = [ blockstr ]
        else:
            ny = Ycells*(1+(not sym))
            rXY = (xmax - xmin)/(ymax - ymin)
            rYZ = (zmax - zmin)/(ymax - ymin)
            ny = int(round(ny*(np.sqrt(2)**(gridlvl-1))))
            nz = int(round(ny*rYZ*cellRatio))
            if ndim>2: nx= int(round(ny*rXY))
            else: nx = 1
            res["blocks"]  = '( hex (0 1 2 3 4 5 6 7) ( {:d} {:d} {:d} ) simpleGrading (1 1 1) )'.format(nx, ny,nz)

        res["edges"]   = '()'

        if createPatch:
            if patches is not None: res["patches"] = patches
            else: res["patches"] = default_patches
        else:
            faces = DictProxy()
            faces["type"] = 'patch'
            faces["faces"] = '()'
            res["boundary"] = ["defaultFaces", faces]

        res["mergePatchPairs"]   = '()'
        return res


if __name__ == "__main__" :
    print(BlockMeshDict.Build("test"))



