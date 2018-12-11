#!/usr/bin/env python

#########################################################################
# Filename: dropTestMesher.py                                           #
# Date:     2018-May-07                                                 #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################
#  This class can be used to generate a mesh for a drop test case i.e.  #
#  with no wave propagation. It work for both 2D and 3D configurations. #
#########################################################################

import re
import os
import math as mt
import numpy as np
import pandas as pd
from io import StringIO
from subprocess import call
from scipy import interpolate as interp
from fsTools import findBoundingBox, findSTLPatches

from ofMesher import OfMesher

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.extrudeMeshDict import ExtrudeMeshDict
from inputFiles.refineMeshDict import RefineMeshDict
from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.blockMeshDict import BlockMeshDict

class DropTestMesher( OfMesher ):
    """Class used to generate a CFD mesh for drop test cases (i.e. without wave propagation)
     This class functions should be used by a python script as presented in the following examples.

    Examples
    --------
    For 2D mesh generation (slamming sections)
    
    >>> import os
    >>> #
    >>> # Import meshing routine from foamBazar
    >>> # WARNING : foamBazar/pythonScripts must be added to PYTHONPATH !
    >>> from dropTestMesher import DropTestMesher
    >>> from fsTools import readSections, createSectionStl
    >>> #
    >>> # Input parameters
    >>> sections = [1,2,3,4,5,6,7,8,9,10]
    >>> sectionsFile = 'Slamming_sections.out'
    >>> genSTL = True
    >>> #
    >>> # Read sections (from Homer) and create corresponding STL files
    >>> if genSTL:
    >>>     sdict = readSections(sectionsFile,sections=sections)
    >>>     createSectionStl(sdict)
    >>> #
    >>> # Case name and additional optional parameters
    >>> for isect in sections:
    >>>     case = r'mesh/section_{}'.format(isect)
    >>>     myParams = {'OFversion' : 'P',
    >>>                 'nProcs' : 1,
    >>>                 'sectionsFile' : 'Slamming_sections.out',
    >>>                 'section' : isect,
    >>>                 'symmetry' : False,
    >>>                 'rot' : [0.0,0.0,0.0],
    >>>                 'zBounds' : [-6.,4.],
    >>>                 'nRefBoxes' : 3,
    >>>                 'zRefineBox' : [-2.0,-1.7,-1.5,3.0,2.5,2.0],
    >>>                 'zRefineBox' : [-2.0,-2.0,-2.0,3.5,3.5,3.5],
    >>>                 'nfsRefBoxes' : 3,
    >>>                 'fsRefineBox' : [-1.0,-0.7,-0.5,2.0,1.7,1.5],
    >>>                 'refineLength' : [0.15,0.375,0.675],
    >>>                 'cellRatio' : 2,
    >>>                 'layerLength' : 0.004,
    >>>                 'onLiger' : True
    >>>                 }
    >>> #
    >>> # Call routines for meshing here
    >>> drop = DropTestMesher.BuildFromAllParameters( case, **myParams )
    >>> fname = os.path.join(case,'log.input')
    >>> with open(fname,'w') as f: f.write(str(myParams))
    >>> drop.runInit()
    
    For 3D mesh generation
    
    >>> import os
    >>> #
    >>> # Import meshing routine from foamBazar
    >>> # WARNING : foamBazar/pythonScripts must be added to PYTHONPATH !
    >>> from dropTestMesher import DropTestMesher
    >>> #
    >>> # Case name and additional optional parameters
    >>> case = "mesh"
    >>> myParams = {'OFversion'     : 'P',
    >>>             'nProcs'        : 12,
    >>>             'ndim'          : 3,
    >>>             'stlFile'       : 'ship.stl',
    >>>             'symmetry'      : True,
    >>>             'domain'        : [-1.0,2.0,-9.0,9.0,-6.0,3.0],
    >>>             'nRefBoxes'     : 3,
    >>>             'xRefineBox'    : [-0.3,-0.2,-0.1,1.3,1.2,1.1],
    >>>             'zRefineBox'    : [-2.0,-1.7,-1.5,3.0,2.5,2.0],
    >>>             'nfsRefBoxes'   : 3,
    >>>             'fsRefineBox'   : [-1.0,-0.7,-0.5,2.0,1.7,1.5],
    >>>             'refineLength'  : [0.15,0.375,0.675],
    >>>             'layerLength'   : 0.004,
    >>>             'cellRatio'     : 2,
    >>>             'onLiger'       : True
    >>>            }
    >>> #
    >>> # Call routines for meshing here
    >>> drop = DropTestMesher.BuildFromAllParameters( case, **myParams )
    >>> drop.writeFiles()
    >>> fname = os.path.join(case,'log.input')
    >>> with open(fname,'w') as f: f.write(str(myParams))
    >>> drop.runInit()   
    
    """
    @classmethod
    def BuildFromAllParameters(cls,      case,
                                         nProcs           = 4,
                                         ndim             = 2,
                                         stlFile          = None,
                                         hullPatch        = None,
                                         section          = 1,
                                         gridLevel        = 1,
                                         symmetry         = False,
                                         trans            = [0.0,0.0,0.0],
                                         rot              = [0.0,0.0,0.0],
                                         domain           = [-0.5,0.5,-9.,9.,-2.,3.],
                                         cellRatio        = 1,
                                         fsBounds         = [-0.5,1.5],
                                         nRefBoxes        = 6,
                                         xRefineBox       = [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,1.6,1.5,1.4,1.3,1.2,1.1],
                                         yRefineBox       = [-6.5,-4.5,-3.0,-2.0,-1.5,-1.2,6.5,4.5,3.0,2.0,1.5,1.2],
                                         zRefineBox       = [-2.0,-2.0,-2.0,-1.5,-1.0,-0.5,3.0,2.5,2.0,1.8,1.6,1.4],
                                         nfsRefBoxes      = 5,
                                         fsRefineBox      = [-1.5,-1.0,-0.5,-0.3,-0.2,2.5,2.0,1.5,1.3,1.2],
                                         refineLength     = [0.1],
                                         layerLength      = 0.005,
                                         solver           = "snappyHexMesh",
                                         OFversion        = 5,
                                         onLiger          = False,
                                         clean            = False
                                         ):
        """Build mesh for CFD drop test mesk from a few parameters.
        
        Parameters
        ----------
        case : str
            Name of case to create
        nProcs : int, default 4
            Number of processors used to build the mesh
        ndim : int, default 2
            Number of dimentions of mesh (2 for 2D mesh or 3 for 3D mesh)
        stlFile : str
            STl files name to be used for snapping
        hullPatch : str
            Name of hull patch
        section : int
            Section number (for 2D slamming sections mesh)
        gridLevel : int, default 1
            TODO
        symmetry : bool, default False
            Logical defining if symmetric mesh (Y) should be created
        
        trans : list of float, dimension(3), default [0.0,0.0,0.0]
            Vector defining STl file translation in mesh
        rot : list of float, dimension(3), default [0.0,0.0,0.0]
            Vector defining STl file rotation in mesh (angles provided in degrees)

        domain : list of floats, dimension(6), default [-0.5,0.5,-9.,9.,-2.,3.]
            Coefficients defining overall domain size (relative to STL size) [Xmin, Xmax, Ymin, Ymax, Zmin, Zmax]
            Note: Xmin and Xmax are defined as absolute values for 2D meshes
        cellRatio : int, default 1
            Background mesh cell ratio
        fsBounds : list of floats, dimension(2), default [-0.5,1.5]
            Coefficients defining free surface area (relative to STL hight) [Zmin, Zmax]
        
        nRefBoxes : int, default 6
            Number of refinement boxes to create
        xRefineBox : list of floats, dimension(2*nRefBoxes), default [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,1.6,1.5,1.4,1.3,1.2,1.1]
            List of coefficients defining size of refinement boxes in X direction (relative to STL x-size)
            format: [Xmin_1,Xmin_2,Xmin_3,...,Xmax_1,Xmax_2,Xmax_3]
        yRefineBox : list of floats, dimension(2*nRefBoxes), default [-6.5,-4.5,-3.0,-2.0,-1.5,-1.2,6.5,4.5,3.0,2.0,1.5,1.2]
            List of coefficients defining size of refinement boxes in Y direction (relative to STL y-size)
            format: [Ymin_1,Ymin_2,Ymin_3,...,Ymax_1,Ymax_2,Ymax_3]
        zRefineBox : list of floats, dimension(2*nRefBoxes), default [-2.0,-2.0,-2.0,-1.5,-1.0,-0.5,3.0,2.5,2.0,1.8,1.6,1.4]
            List of coefficients defining size of refinement boxes in Z direction (relative to STL height)
            format: [Zmin_1,Zmin_2,Zmin_3,...,Zmax_1,Zmax_2,Zmax_3]
        nfsRefBoxes, int, default 5
            Number of free surface refinement boxes to create
        fsRefineBox : list of floats, dimension(2*nfsRefBoxes), default [-1.5,-1.0,-0.5,-0.3,-0.2,2.5,2.0,1.5,1.3,1.2]
            List of coefficients defining size of free surface refinement boxes (relative to STL height)
            example: [Zmin_1,Zmin_2,Zmin_3,...,Zmax_1,Zmax_2,Zmax_3]
        
        refineLength : list of float, default [0.1]
            Coeffiecients defining refinements around the body (relative to STL reference length defined as 'min(0.5*y-size,z-size)')
            The length of the list defines the number of proximity refinements.
        layerLength : float, default 0.005
            Coefficient defining length of boundary layer (relative to STL reference length defined as 'min(0.5*y-size,z-size)')
        
        solver : str, default 'snappyHexMesh'
            Solver to run
        OFversion : int or str, default 5
            OpenFOAM version
        onLiger : boolean, default False
            Logical defining if case is run on Liger cluster
        clean : boolean, default False
            Logical to force case overwrite
     
        """
        
        if hullPatch is None:
            if ndim==2: hullPatch='section_'+str(section)
            elif ndim==3: hullPatch = 'ship'
        
        #get STL size
        if (ndim==2) and (stlFile is None): stlFile = os.path.join('stl',hullPatch+'.stl')
        stlFile = os.path.join(os.getcwd(),stlFile)
        
        bBox = findBoundingBox(stlFile, verbose=True)
        Length = bBox[3]-bBox[0]
        Beam   = bBox[4]-bBox[1]
        Depth  = bBox[5]-bBox[2]
        print('Length = {}; Beam = {}; Depth = {}'.format(Length,Beam,Depth))
        
        print('Create system folder input files')
        #controlDict
        controlDict = ControlDict.Build(case                = case,
                                        version             = "snappyHexMesh",
                                        endTime             = 100,
                                        deltaT              = 1,
                                        writeControl        = "runTime",
                                        writeInterval       = 1,
                                        writeCompression    = "off",
                                        runTimeModifiable   = "true")
                                  
        #fvSchemes
        fvSchemes = FvSchemes.Build(case     = case,
                              simType  = "Euler" )

        #fvSolution
        fvSolution = FvSolution.Build(case = case )
        
        #decomposeParDict
        decomposeParDict = DecomposeParDict.Build(case   = case,
                                            nProcs = nProcs,
                                            method = "scotch")
        
        #extrudeMeshDict
        extrudeMeshDict = ExtrudeMeshDict.Build(case = case)
        
        #refineMeshDict
        refineMeshDict = RefineMeshDict.Build(case                = case,
                                        refineUptoCellLevel = 5)
                                        
        
        refineMeshDict1 = RefineMeshDict.Build(case           = case,
                                         orient         = 'z',
                                         set            = "c0",
                                         useHexTopology = True,
                                         geometricCut   = False)
        if ndim==2:
            orient = 'yz'
            directions = "tan2 normal"
        elif ndim==3:
            orient = 'xyz'
            directions = "tan1 tan2 normal"
         
        refineMeshDict2 = RefineMeshDict.Build(case         = case,
                                         orient         = orient,
                                         set            = "c0",
                                         directions     = directions,
                                         useHexTopology = True,
                                         geometricCut   = False)
                                         
        #read stl patches
        stlPatches = findSTLPatches(stlFile)
        
        #snappyHexMeshDict
        stlName = os.path.splitext(os.path.basename(stlFile))[0]
        referenceLength = min(Beam*0.5,Depth)
        snappyHexMeshDict = SnappyHexMeshDict.Build(case                = case,
                                              addLayers           = True,
                                              stlname             = stlName+'.stl',
                                              patchName           = hullPatch,
                                              stlPatches          = stlPatches,
                                              refinementLength    = [referenceLength*rf for rf in refineLength],
                                              nSurfaceLayers      = 3,
                                              expansionRatio      = 1.3,
                                              finalLayerThickness = referenceLength*layerLength,
                                              minThickness        = referenceLength*layerLength*0.1,
                                              OFversion           = OFversion)
                                              
        #surfaceFeatureExtractDict
        surfaceFeatureExtractDict = SurfaceFeatureExtractDict.Build(case = case,
                                                              stlname = stlName+'.stl')

        print('Create constant folder input files')
        #blockMeshDict
        if ndim==2:
            blockMeshDict = BlockMeshDict.Build( case      = case,
                                           ndim      = ndim,
                                                xmin      = domain[0],
                                                xmax      = domain[1],
                                                ymin      = domain[2]*Beam*0.5*(not symmetry),
                                                ymax      = domain[3]*Beam*0.5,
                                                zmin      = domain[4]*Depth,
                                                zmax      = domain[5]*Depth,
                                           fsmin     = fsBounds[0]*Depth,
                                           fsmax     = fsBounds[1]*Depth,
                                           sym       = symmetry,
                                           cellRatio = cellRatio,
                                           gridlvl   = gridLevel,
                                           OFversion = OFversion)
        elif ndim==3:
            blockMeshDict = BlockMeshDict.Build( case      = case,
                                           ndim      = ndim,
                                                xmin      = domain[0]*Length,
                                                xmax      = domain[1]*Length,
                                                ymin      = domain[2]*Beam*0.5*(not symmetry),
                                                ymax      = domain[3]*Beam*0.5,
                                                zmin      = domain[4]*Depth,
                                                zmax      = domain[5]*Depth,
                                           fsmin     = fsBounds[0]*Depth,
                                           fsmax     = fsBounds[1]*Depth,
                                           sym       = symmetry,
                                           cellRatio = cellRatio,
                                           gridlvl   = gridLevel,
                                           OFversion = OFversion)

        res =  cls( case, nProcs=nProcs,
                          controlDict=controlDict,
                          fvSchemes=fvSchemes,
                          fvSolution=fvSolution,
                          decomposeParDict=decomposeParDict,
                          extrudeMeshDict=extrudeMeshDict,
                          refineMeshDicts=[refineMeshDict,refineMeshDict1,refineMeshDict2],
                          snappyHexMeshDict=snappyHexMeshDict,
                          surfaceFeatureExtractDict=surfaceFeatureExtractDict,
                          blockMeshDict=blockMeshDict,
                          solver = solver,
                          isMesher = True,
                          clean = clean
                          )

        res.section = section
        res.stlFile = stlFile
        res.stlName = stlName
        res.ndim = ndim
        res.trans = trans
        res.rot = rot
        res.Length = Length
        res.Beam = Beam
        res.Depth = Depth
        res.hullPatch = hullPatch
        res.symmetry = symmetry
        res.nRefBoxes = nRefBoxes
        res.xRefineBox = xRefineBox
        res.yRefineBox = yRefineBox
        res.zRefineBox = zRefineBox
        res.nfsRefBoxes = nfsRefBoxes
        res.fsRefineBox = fsRefineBox
        res.OFversion = OFversion
        res.onLiger = onLiger
        
        res.writeFiles()
        return res
        
    def writeFiles(self) :
        OfMesher.writeFiles(self)
        
        #Write additional files
        


    def readSections(inputFile,sections=[]):
        
        sdict = {}
        
        with open(inputFile) as f:
            if len(sections)==0:
                sect = re.findall(r'(#Section \d+)',f.read())
                for s in sect:
                    sections += [int(i) for i in s.split() if i.isdigit()]
    
        with open(inputFile,'r') as f:
            lines = f.read()
        
        for isect in sections:
            string = r'(#Section '+str(isect)+r'\n)+([^#]*)'
            table = re.findall(string, lines)
            table = list(table[0])
            rows = StringIO(table[1])
            data = pd.read_csv(rows,header=None,sep=r'\s+',names=['y','z'])
            sdict[isect] = data[data.y>=0]
            sdict[isect] = sdict[isect].assign(x= sdict[isect]['y'].values * 0.0)
            sdict[isect] = sdict[isect].reindex(sorted(sdict[isect]),axis=1)
                
        return sdict

    def createSectionStl(sdict):
        if not os.path.exists('geo'): os.makedirs('geo')
        if not os.path.exists('stl'): os.makedirs('stl')
        
        for isect in sdict.keys():
            print('Section =',isect)
            name = 'section_'+str(isect)
            mysect = sdict[isect]
            
            fgeo=os.path.join(r'geo/',name+'.geo');
            fstl=os.path.join(r'stl/',name+'.stl');
        
            #Create .geo file
            thk = 10.0
            
            #Redistribute points
            s = np.sqrt((mysect.diff()**2).sum(axis=1)).cumsum()
            ds = 0.5 * np.sqrt((mysect.diff()*2).sum(axis=1)).mean()
            ss = np.linspace(s.min(), s.max(), mt.ceil((s.max()-s.min())/ds))
            
            #interpolate according to curvilinear absissa
            fx = interp.interp1d(s,mysect.x); x = fx(ss)
            fy = interp.interp1d(s,mysect.y); y = fy(ss)
            fz = interp.interp1d(s,mysect.z); z = fz(ss)
            
            #create symmetric
            interpsect = pd.DataFrame({ 'x': np.concatenate((x,x[1:])),
                                        'y': np.concatenate(( -1.0*y[:0:-1],y)),
                                        'z': np.concatenate(( z[:0:-1],z)) })
            ns = mt.ceil(16*len(ss))
    
            pSize = len(interpsect)
            
            lSize=1.0
            nx=1.0
            ny=2.0*mt.ceil((interpsect.max().y-interpsect.min().y)/(thk/nx))
    
            print('writing to file: ', fgeo)
            f = open(fgeo,'w')
            f.write('// Section {:d}\n'.format(isect))
            f.write('thk={:.6e}; xmin=-0.5*thk;\n'.format(thk))   
            f.write('p0=newp; pSize={:d};\n'.format(pSize))
            for i in range(pSize):
                f.write('Point(p0+{:d}) = {{xmin, {:.6e}, {:.6e}}};\n'.format(i,interpsect.loc[i].y,interpsect.loc[i].z))
            f.write('l0=newl; lSize={:.2f};\n'.format(lSize))
            f.write('Spline(l0+0) = {{p0:p0+{:d}}};\n'.format(pSize - 1))
            f.write('Line(l0+lSize) = {p0+pSize-1,p0};\n')
            f.write('Extrude {thk, 0, 0} {\n')
            f.write(' Line{l0:l0+lSize};\n')
            f.write(' Recombine;\n')
            f.write('}\n')
            f.write('Transfinite Line {{4, 5}} = {:.5f} Using Progression 1;\n'.format(nx))
            f.write('Transfinite Line {{2, 6}} = {:.5f} Using Progression 1;\n'.format(ny))
            f.write('Transfinite Line {{1, 3}} = {:.5f} Using Progression 1;'.format(ns))
            f.write('Transfinite Surface "*";\n')
            f.write('Recombine Surface "*";\n')
            f.close()
    
            #Create STL
            call('gmsh -2 -format stl -o "'+fstl+'" "'+fgeo+'"', shell=True)
            
            
    def writeAllinit(self):
        #Allinit
        print('Create run scripts')
        ainit = os.path.join(self.case,'Allinit')
        with open(ainit,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')
            
            #clearMeshInfo
            f.write('function clearMeshInfo()\n')
            f.write('{\n')
            f.write('    rm -fr background.msh 0/{ccx,ccy,ccz,*Level,polyMesh/cellMap}* constant/polyMesh/{sets,*Level*,*Zone*,refine*,surfaceIndex*}\n')
            f.write('}\n\n')
            
            #copySTL
            f.write('function copySTL()\n')
            f.write('{\n')
            fstlin = os.path.join(self.stlFile)
            fstlout = os.path.join('constant','triSurface',self.stlName+'.stl')
            fstltmp = os.path.join('constant','triSurface',self.stlName+'_tmp.stl')
            f.write('    cp {:s} {:s}\n'.format(fstlin,fstlout))
            if any(self.trans)!=0.0:
                f.write('    mv {:s} {:s}\n'.format(fstlout,fstltmp))
                f.write("    surfaceTransformPoints -translate '({:.6f} {:.6f} {:.6f})' {} {}\n".format(*self.trans,fstltmp,fstlout))
                f.write('    rm {:s}\n'.format(fstltmp))
            if any(self.rot)!=0.0:
                f.write('    mv {:s} {:s}\n'.format(fstlout,fstltmp))
                f.write("    surfaceTransformPoints -rollPitchYaw '({:.3f} {:.3f} {:.3f})' {} {}\n".format(*self.rot,fstltmp,fstlout))
                f.write('    rm {:s}\n'.format(fstltmp))
            f.write('    surfaceFeatureExtract\n')
            f.write('}\n\n')
            
            #refineBox
            f.write('function refineBox()\n')
            f.write('{\n')
            f.write('    echo "cellSet c0 new boxToCell $1" | setSet -latestTime\n')
            if self.ndim==2:   f.write('    eval BVrefineMesh -dict "system/refineMeshDict.yz"\n')
            elif self.ndim==3: f.write('    eval BVrefineMesh -dict "system/refineMeshDict.xyz"\n')
            f.write('}\n\n')
            
            #refineFS
            f.write('function refineFS()\n')
            f.write('{\n')
            f.write('    REFLEVEL=$1\n')
            f.write('    PARAM=$2\n')
            f.write('    shift; shift\n')
            f.write('    echo "cellSet c0 new boxToCell $PARAM" | setSet -latestTime\n')
            f.write('    eval BVrefineMesh -dict "system/refineMeshDict.z" -refineUpToLevel $REFLEVEL $@\n')
            f.write('}\n\n')
            
            #snap
            f.write('function snap()\n')
            f.write('{\n')
            if self.nProcs>1:
                f.write('    decomposePar -cellDist\n')
                if self.onLiger:
                    self.writeSbatch()
                    f.write('    sbatch run.sh\n')
                else:
                    f.write('    mpirun snappyHexMesh -parallel\n')
                    f.write('    reconstructParMesh -latestTime\n')
            else:
                f.write('    snappyHexMesh\n')
            f.write('}\n\n')
            
            # __main__
            f.write('(\n')
            f.write('    copySTL\n')
            f.write('    blockMesh\n')
            
            nx = int(len(self.xRefineBox)/2)
            ny = int(len(self.yRefineBox)/2)
            nz = int(len(self.zRefineBox)/2)
            nf = int(len(self.fsRefineBox)/2)
            if self.ndim==2: n = min(ny,nz)
            elif self.ndim==3: n = min(nx,ny,nz)
                
            if n < self.nRefBoxes:
                print('ERROR: Requested number of boxes is greater than refinement level provided!')
                os._exit(1)
            for i in range(self.nRefBoxes):
                if self.ndim==2:
                    f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(self.yRefineBox[i]*self.Beam*0.5,self.zRefineBox[i]*self.Depth,self.yRefineBox[i+ny]*self.Beam*0.5,self.zRefineBox[i+nz]*self.Depth))
                elif self.ndim==3:
                    f.write('    refineBox "({:.2f} {:.2f} {:.2f}) ({:.2f} {:.2f} {:.2f})"\n'.format(self.xRefineBox[i]*self.Length,self.yRefineBox[i]*self.Beam*0.5,self.zRefineBox[i]*self.Depth,self.xRefineBox[i+nx]*self.Length,self.yRefineBox[i+ny]*self.Beam*0.5,self.zRefineBox[i+nz]*self.Depth))
            
            if nf < self.nfsRefBoxes:
                print('ERROR: Requested number of boxes is greater than refinement level provided!')
                os._exit(1)
            for i in range(min(self.nfsRefBoxes,self.nRefBoxes)):
                if i < min(self.nfsRefBoxes,self.nRefBoxes)-1:
                    f.write('    refineFS {} "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(self.nRefBoxes,self.fsRefineBox[i]*self.Depth,self.fsRefineBox[i+nf]*self.Depth))
                elif i == min(self.nfsRefBoxes,self.nRefBoxes)-1:
                    f.write('    refineFS {} "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})" -noCellLevel\n'.format(self.nRefBoxes,self.fsRefineBox[i]*self.Depth,self.fsRefineBox[i+nf]*self.Depth))
            f.write('    snap\n')
            if self.ndim==2 and (not (self.onLiger and self.nProcs>1)): f.write('    extrudeMesh\n')
            f.write(') 2>&1 | tee log.mesh\n\n')
            if not (self.onLiger and self.nProcs>1): f.write('checkMesh -latestTime 2>&1 | tee log.checkMesh\n')
        os.chmod(ainit, 0o755)

    def writeAllclean(self):
        #Allclean
        aclean = os.path.join(self.case,'Allclean')
        with open(aclean,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('cd ${0%/*} || exit 1    # run from this directory\n\n')
            f.write('function clean_log()\n')
            f.write('{\n')
            f.write('    rm -fr log.*\n')
            f.write('}\n\n')
            f.write('function clean_mesh()\n')
            f.write('{\n')
            f.write('    rm -fr background.msh VTK\n')
            f.write('    rm -fr 0/{ccx,ccy,ccz,*Level,polyMesh/cellMap}*\n')
            f.write('    rm -fr constant/modeshapes/{cfd,mapper,modalInfo}_body\n')
            f.write('    rm -fr constant/extendedFeatureEdgeMesh/\n')
            f.write('    rm -fr constant/polyMesh/{sets,*Level*,*level*,*Index*,*History*}\n')
            f.write('}\n\n')
            f.write('function clean_parallel_mesh()\n')
            f.write('{\n')
            f.write('    rm -fr processor*\n')
            f.write('}\n\n')
            f.write('function clean_0()\n')
            f.write('{\n')
            f.write('    rm -fr 0/*.gz\n')
            f.write('}\n\n')
            f.write('#eval clean_log\n')
            f.write('eval clean_mesh\n')
            f.write('eval clean_parallel_mesh\n')
            f.write('eval clean_0\n')
        os.chmod(aclean, 0o755)
        
        
    def writeSbatch(self):
        #run.sh
        run = os.path.join(self.case,'run.sh')
        with open(run,'w') as f:
            f.write('#!/bin/bash -l\n')
            f.write('#SBATCH -J {}\n\n'.format(self.case))
            f.write('# 5 hour wall-clock\n')
            f.write('#SBATCH --account I1608251\n')
            f.write('#SBATCH -t 0-02:00:00\n')
            f.write('#SBATCH -n {:d}\n'.format(self.nProcs))
            f.write('#SBATCH -o log.run-%j\n\n')
            f.write('module load gcc/4.9.3 openmpi/1.8.4-gcc lapack/3.6.1/gcc/4.9.3\n')
            f.write('export FOAM_INST_DIR=/data/I1608251/OpenFOAM;\n')
            if   self.OFversion == 2 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-2.4.x/etc/bashrc;\n')
            elif self.OFversion == 3 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-3.0.x/etc/bashrc;\n')
            elif self.OFversion == 5 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-5.x/etc/bashrc;\n')
            elif 'p' in str(self.OFversion).lower(): f.write('source /data/I1608251/OpenFOAM/OpenFOAM-v1712/etc/bashrc;\n')
            f.write('export LC_ALL=C\n\n')
            f.write('mpirun {} -parallel\n'.format(self.solver))
