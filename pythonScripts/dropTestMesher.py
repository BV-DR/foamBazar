#!/usr/bin/env python

#########################################################################
# Filename: fsMesher2D.py                                               #
# Date:     2018-May-07                                                 #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import re
import os, sys, argparse
import shutil
import time, datetime
import math as mt
import numpy as np
import pandas as pd
from io import StringIO
from subprocess import call, Popen
from scipy import interpolate as interp
from fsTools import findBoundingBox, findSTLPatches

from ofCase import OfCase

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.extrudeMeshDict import ExtrudeMeshDict
from inputFiles.refineMeshDict import RefineMeshDict
from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.blockMeshDict import BlockMeshDict


class DropTestMesher( OfCase ):

    @classmethod
    def BuildFromAllParameters(cls,      case,
                                         nProcs           = 4,
                                         ndim             = 2,
                                         sectionsFile     = None,
                                         stlFile          = None,
                                         hullPatch        = None,
                                         section          = 1,
                                         gridLevel        = 1,
                                         symmetry         = False,
                                         trans            = [0.0,0.0,0.0],
                                         rot              = [0.0,0.0,0.0],
                                         xBounds          = [-0.5,0.5],
                                         yBounds          = [-9.,9.],
                                         zBounds          = [-2.,3.],
                                         fsBounds         = [-0.5,1.5],
                                         nRefBoxes        = 6,
                                         xRefineBox       = [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,1.6,1.5,1.4,1.3,1.2,1.1],
                                         yRefineBox       = [-6.5,-4.5,-3.0,-2.0,-1.5,-1.2,6.5,4.5,3.0,2.0,1.5,1.2],
                                         zRefineBox       = [-2.0,-2.0,-2.0,-1.5,-1.0,-0.5,3.0,2.5,2.0,1.8,1.6,1.4],
                                         nfsRefBoxes      = 5,
                                         fsRefineBox      = [-1.5,-1.0,-0.5,-0.3,-0.2,2.5,2.0,1.5,1.3,1.2],
                                         refineLength     = [0.1],
                                         layerLength      = 0.005,
                                         cellRatio        = 1,
                                         solver           = "snappyHexMesh",
                                         OFversion        = 3,
                                         onLiger          = False,
                                         overwrite        = False
                                         ):
        
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
        controlDict = ControlDict(case                = case,
                                  version             = "snappyHexMesh",
                                  endTime             = 100,
                                  deltaT              = 1,
                                  writeControl        = "runTime",
                                  writeInterval       = 1,
                                  writeCompression    = "off",
                                  runTimeModifiable   = "true")
                                  
    
        #fvSchemes
        fvSchemes = FvSchemes(case     = case,
                              simType  = "Euler" )

        #fvSolution
        fvSolution = FvSolution(case = case )
        
        #decomposeParDict
        decomposeParDict = DecomposeParDict(case   = case,
                                            nProcs = nProcs,
                                            method = "scotch")
        
        #extrudeMeshDict
        extrudeMeshDict = ExtrudeMeshDict(case = case)
        
        #refineMeshDict
        refineMeshDict = RefineMeshDict(case                = case,
                                        refineUptoCellLevel = 5)
                                        
        
        refineMeshDict1 = RefineMeshDict(case           = case,
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
         
        refineMeshDict2 = RefineMeshDict(case         = case,
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
        snappyHexMeshDict = SnappyHexMeshDict(case                = case,
                                              addLayers           = True,
                                              stlname             = stlName+'.stl',
                                              patchName           = hullPatch,
                                              stlPatches          = stlPatches,
                                              refinementLength    = [referenceLength*rf for rf in refineLength],
                                              nSurfaceLayers      = 3,
                                              expansionRatio      = 1.3,
                                              finalLayerThickness = referenceLength*layerLength,
                                              minThickness        = referenceLength*layerLength*0.1,
                                              ofp                 = OFversion=='P')
                                              
        #surfaceFeatureExtractDict
        surfaceFeatureExtractDict = SurfaceFeatureExtractDict(case = case,
                                                              stlname = stlName+'.stl')

        print('Create constant folder input files')
        #blockMeshDict
        if ndim==2:
            blockMeshDict = BlockMeshDict( case      = case,
                                           ndim      = ndim,
                                           xmin      = xBounds[0],
                                           xmax      = xBounds[1],
                                           ymin      = yBounds[0]*Beam*0.5*(not symmetry),
                                           ymax      = yBounds[1]*Beam*0.5,
                                           zmin      = zBounds[0]*Depth,
                                           zmax      = zBounds[1]*Depth,
                                           fsmin     = fsBounds[0]*Depth,
                                           fsmax     = fsBounds[1]*Depth,
                                           sym       = symmetry,
                                           cellRatio = cellRatio,
                                           gridlvl   = gridLevel,
                                           ofp       = OFversion=='P')
        elif ndim==3:
            blockMeshDict = BlockMeshDict( case      = case,
                                           ndim      = ndim,
                                           xmin      = xBounds[0]*Length,
                                           xmax      = xBounds[1]*Length,
                                           ymin      = yBounds[0]*Beam*0.5*(not symmetry),
                                           ymax      = yBounds[1]*Beam*0.5,
                                           zmin      = zBounds[0]*Depth,
                                           zmax      = zBounds[1]*Depth,
                                           fsmin     = fsBounds[0]*Depth,
                                           fsmax     = fsBounds[1]*Depth,
                                           sym       = symmetry,
                                           cellRatio = cellRatio,
                                           gridlvl   = gridLevel,
                                           ofp       = OFversion=='P')

        res =  cls( case, nProcs=nProcs,
                          controlDict=controlDict,
                          fvSchemes=fvSchemes,
                          fvSolution=fvSolution,
                          decomposeParDict=decomposeParDict,
                          extrudeMeshDict=extrudeMeshDict,
                          refineMeshDict=refineMeshDict,
                          refineMeshDict1=refineMeshDict1,
                          refineMeshDict2=refineMeshDict2,
                          snappyHexMeshDict=snappyHexMeshDict,
                          surfaceFeatureExtractDict=surfaceFeatureExtractDict,
                          blockMeshDict=blockMeshDict,
                          solver = solver,
                          isMesher = True,
                          overwrite = overwrite
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
        # if ndim==2:
            # res.section_name = 'section_'+str(section)
            # res.sdict = res.readSections(sectionsFile,sections=[section])
        # elif ndim==3:
            # res.section_name = 'ship'
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
        OfCase.writeFiles(self)


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
            elif self.OFversion == 'P' : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-v1712/etc/bashrc;\n')
            f.write('export LC_ALL=C\n\n')
            f.write('mpirun {} -parallel\n'.format(self.solver))
