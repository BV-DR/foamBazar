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
from fsTools import findBoundingBox
from droppy.Reader.readInput import readInput

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.extrudeMeshDict import ExtrudeMeshDict
from inputFiles.refineMeshDict import RefineMeshDict
from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.blockMeshDict import BlockMeshDict

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
        
def create2DMesh(param,sdict,runScript):
    for isect in sdict.keys():
        print('')
        print('Section =',isect)
        
        if param.symmetry: mysym = 'sym'
        else: mysym = 'r{:02d}'.format(int(param.rot))
        
        mysect = 'section_'+str(isect)
        
        for grid in param.gridLevel:
            if param.gridName is not None: mygrid = param.gridName 
            else: mygrid = 'grid_'+str(grid)
            mycase = os.path.join(mysym,mysect,mygrid)
            length = sdict[isect].y.max()
            height = sdict[isect].z.max()
            
            if not os.path.exists(mycase): os.makedirs(mycase)
            
            print('Create system folder input files')
            sysfolder = os.path.join(mycase,"system")
            if not os.path.exists(sysfolder): os.makedirs(sysfolder)
            
            #controlDict
            controlDict = ControlDict(case                = mycase,
                                      version             = "snappyHexMesh",
                                      endTime             = 100,
                                      deltaT              = 1,
                                      writeControl        = "runTime",
                                      writeInterval       = 1,
                                      writePrecision      = 6,
                                      writeCompression    = "off",
                                      runTimeModifiable   = "true")
            controlDict.writeFile()
    
            #fvSchemes
            fvSchemes = FvSchemes(case     = mycase,
                                simType  = "Euler" )
            fvSchemes.writeFile()
    
            #fvSolution
            fvSolution = FvSolution(case = mycase )
            fvSolution.writeFile()
        
            #decomposeParDict
            decomposeParDict = DecomposeParDict(case   = mycase,
                                                nProcs = 4,
                                                method = "simple")
            decomposeParDict.writeFile()
            
            #extrudeMeshDict
            extrudeMeshDict = ExtrudeMeshDict(case = mycase)
            extrudeMeshDict.writeFile()
            
            #refineMeshDict
            refineMeshDict = RefineMeshDict(case                = mycase,
                                            refineUptoCellLevel = 5)
            refineMeshDict.writeFile()
            
            refineMeshDictZ = RefineMeshDict(case           = mycase,
                                            orient         = 'z',
                                            set            = "c0",
                                            useHexTopology = True,
                                            geometricCut   = False)
            refineMeshDictZ.writeFile()
            
            refineMeshDictYZ = RefineMeshDict(case           = mycase,
                                            orient         = 'yz',
                                            set            = "c0",
                                            directions     = "tan2 normal",
                                            useHexTopology = True,
                                            geometricCut   = False)
            refineMeshDictYZ.writeFile()
            
            #snappyHexMeshDict
            snappyHexMeshDict = SnappyHexMeshDict(case                = mycase,
                                                  addLayers           = True,
                                                  stlname             = mysect,
                                                  refinementLength    = [length*rf for rf in param.refineLength],
                                                  nSurfaceLayers      = 3,
                                                  expansionRatio      = 1.3,
                                                  finalLayerThickness = length*param.layerLength,
                                                  minThickness        = length*param.layerLength*0.1,
                                                  ofp                 = param.OFplus)
            snappyHexMeshDict.writeFile()
        
            #surfaceFeatureExtractDict
            surfaceFeatureExtractDict = SurfaceFeatureExtractDict(case = mycase,
                                                                stlname = mysect)
            surfaceFeatureExtractDict.writeFile()
            
            print('Create constant folder input files')
            consfolder = os.path.join(mycase,"constant")
            polyfolder = os.path.join(mycase,"constant","polyMesh")
            trifolder = os.path.join(mycase,"constant","triSurface")
            if not os.path.exists(consfolder): os.makedirs(consfolder)
            if not os.path.exists(polyfolder): os.makedirs(polyfolder)
            if not os.path.exists(trifolder): os.makedirs(trifolder)
            
            #blockMeshDict
            print('WRITE BLOCKMESHDICT')
            blockMeshDict = BlockMeshDict( case      = mycase,
                                           ndim      = 2,
                                           xmin      = param.xBounds[0],
                                           xmax      = param.xBounds[1],
                                           ymin      = param.yBounds[0]*length*(not param.symmetry),
                                           ymax      = param.yBounds[1]*length,
                                           zmin      = param.zBounds[0]*height,
                                           zmax      = param.zBounds[1]*height,
                                           fsmin     = param.fsBounds[0]*height,
                                           fsmax     = param.fsBounds[1]*height,
                                           sym       = param.symmetry,
					   cellRatio = param.cellRatio,
                                           gridlvl   = grid,
                                           ofp       = param.OFplus)
            blockMeshDict.writeFile()
    
        
            print('Create 0 folder input files')
            zerofolder = os.path.join(mycase,"0")
            if not os.path.exists(zerofolder): os.makedirs(zerofolder)
        
            #Allrun
            print('Create run scripts')
            arun = os.path.join(mycase,'Allrun')
            with open(arun,'w') as f:
                f.write('#! /bin/bash\n')
                f.write('set -x\n\n')
                f.write('function clearMeshInfo()\n')
                f.write('{\n')
                f.write('    rm -fr background.msh 0/{ccx,ccy,ccz,*Level,polyMesh/cellMap}* constant/polyMesh/{sets,*Level*,*Zone*,refine*,surfaceIndex*}\n')
                f.write('}\n\n')
                f.write('function refineBox()\n')
                f.write('{\n')
                f.write('    echo "cellSet c0 new boxToCell $1" | setSet -latestTime\n')
                f.write('    eval BVrefineMesh -dict "system/refineMeshDict.yz"\n')
                f.write('}\n\n')
                f.write('function refineFS()\n')
                f.write('{\n')
                f.write('    REFLEVEL=$1\n')
                f.write('    PARAM=$2\n')
                f.write('    shift; shift\n')
                f.write('    echo "cellSet c0 new boxToCell $PARAM" | setSet -latestTime\n')
                f.write('    eval BVrefineMesh -dict "system/refineMeshDict.z" -refineUpToLevel $REFLEVEL $@\n')
                f.write('}\n\n')
                f.write('function snap()\n')
                f.write('{\n')
                f.write('    snappyHexMesh\n')
                f.write('}\n\n')
                f.write('(\n')
                f.write('    blockMesh\n')
                
                n = int(min([len(param.yRefineBox),len(param.zRefineBox)/2]))
                if n < param.nRefBoxes:
                    print('ERROR: Requested number of boxes is greater than refinement level provided!')
                    os._exit(1)
                for i in range(param.nRefBoxes):
                    f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(param.yRefineBox[i]*length,param.zRefineBox[i]*height,param.yRefineBox[i+n]*length,param.zRefineBox[i+n]*height))
                n = int(len(param.fsRefineBox)/2)
                if n < param.nfsRefBoxes:
                    print('ERROR: Requested number of boxes is greater than refinement level provided!')
                    os._exit(1)
                for i in range(min(param.nfsRefBoxes,param.nRefBoxes)):
                    if i < min(param.nfsRefBoxes,param.nRefBoxes)-1:
                        f.write('    refineFS {} "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(param.nRefBoxes,param.fsRefineBox[i]*height,param.fsRefineBox[i+n]*height))
                    elif i == min(param.nfsRefBoxes,param.nRefBoxes)-1:
                        f.write('    refineFS {} "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})" -noCellLevel\n'.format(param.nRefBoxes,param.fsRefineBox[i]*height,param.fsRefineBox[i+n]*height))
                f.write('    snap\n')
                f.write('    extrudeMesh\n')
                f.write(') 2>&1 | tee log.mesh\n\n')
                f.write('checkMesh -latestTime 2>&1 | tee log.checkMesh\n')
            os.chmod(arun, 0o755)
        
            #Allclean
            aclean = os.path.join(mycase,'Allclean')
            with open(aclean,'w') as f:
                f.write('#! /bin/bash\n')
                f.write('cd ${0%/*} || exit 1    # run from this directory\n\n')
                f.write('# Source tutorial clean functions\n')
                f.write('. $WM_PROJECT_DIR/bin/tools/CleanFunctions\n\n')
                f.write('cleanCase\n\n')
                f.write('rm -fr 0/alpha.water* 0/U* 0/p_rgh*\n')
            os.chmod(aclean, 0o755)
            
            print('Copy STL file, rotate and run surfaceFeatureExtract')
            #Copy STL in triSurface folder
            fstlin  = os.path.join(r'stl/',mysect+'.stl')
            fstlout = os.path.join(trifolder,mysect+'.stl')
            shutil.copyfile(fstlin,fstlout)
            
            #rotate section if needed
            if abs(param.rot)>1e-3:
                fstltmp = os.path.join(trifolder,mysect+'_tmp.stl')
                shutil.move(fstlout,fstltmp)
                p = Popen(["surfaceTransformPoints -yawPitchRoll '(0 0 {:.3f})' {} {}".format(param.rot,fstltmp,fstlout)], shell=True)
                p.wait()
                os.remove(fstltmp)
            
            #run surfaceFeatureExtract
            p = Popen(['surfaceFeatureExtract'], cwd=mycase)
            p.wait()
            
            #run Allrun script
            if runScript:
                p = Popen(['./Allclean'], cwd=mycase)
                p.wait()
                p = Popen(['./Allrun'], cwd=mycase)
                p.wait()
            
            #create file for Paraview
            open(os.path.join(mycase,'a.foam'), 'a').close()
            
def create3DMesh(param,runScript):
    if param.symmetry: mysym = 'sym'
    else: mysym = 'r{:02d}'.format(int(param.rot))
    
    mystl = param.stlFile.split('.stl')[0] #remove .stl extension
    
    for grid in param.gridLevel:
        if param.gridName is not None: mygrid = param.gridName
        else: mygrid = 'grid_'+str(grid)
        mycase = os.path.join(mysym,mygrid)
        
        #get STL size
        bBox = findBoundingBox(r'stl/'+mystl, verbose=True)
        Length = bBox[3]-bBox[0]
        Beam   = bBox[4]-bBox[1]
        Depth  = bBox[5]-bBox[2]
        print('Length = {}; Beam = {}; Depth = {}'.format(Length,Beam,Depth))
        
        if not os.path.exists(mycase): os.makedirs(mycase)
        
        print('Create system folder input files')
        sysfolder = os.path.join(mycase,"system")
        if not os.path.exists(sysfolder): os.makedirs(sysfolder)
        
        #controlDict
        controlDict = ControlDict(case                = mycase,
                                  version             = "snappyHexMesh",
                                  endTime             = 100,
                                  deltaT              = 1,
                                  writeControl        = "runTime",
                                  writeInterval       = 1,
                                  writePrecision      = 6,
                                  writeCompression    = "off",
                                  runTimeModifiable   = "true")
        controlDict.writeFile()
    
        #fvSchemes
        fvSchemes = FvSchemes(case     = mycase,
                              simType  = "Euler" )
        fvSchemes.writeFile()
    
        #fvSolution
        fvSolution = FvSolution(case = mycase )
        fvSolution.writeFile()
    
        #decomposeParDict
        decomposeParDict = DecomposeParDict(case   = mycase,
                                            nProcs = 4)
        decomposeParDict.writeFile()
        
        #refineMeshDict
        refineMeshDict = RefineMeshDict(case                = mycase,
                                        refineUptoCellLevel = 5)
        refineMeshDict.writeFile()
        
        refineMeshDictZ = RefineMeshDict(case           = mycase,
                                         orient         = 'z',
                                         set            = "c0",
                                         useHexTopology = True,
                                         geometricCut   = False)
        refineMeshDictZ.writeFile()
        
        refineMeshDictXYZ = RefineMeshDict(case           = mycase,
                                           orient         = 'xyz',
                                           set            = "c0",
                                           directions     = "tan1 tan2 normal",
                                           useHexTopology = True,
                                           geometricCut   = False)
        refineMeshDictXYZ.writeFile()
        
        #snappyHexMeshDict
        snappyHexMeshDict = SnappyHexMeshDict(case                = mycase,
                                              addLayers           = True,
                                              stlname             = mystl,
                                              refinementLength    = [Beam*0.5*rf for rf in param.refineLength],
                                              nSurfaceLayers      = 3,
                                              expansionRatio      = 1.3,
                                              finalLayerThickness = Beam*0.5*param.layerLength,
                                              minThickness        = Beam*0.5*param.layerLength*0.1,
                                              ofp                 = param.OFplus)
        snappyHexMeshDict.writeFile()
    
        #surfaceFeatureExtractDict
        surfaceFeatureExtractDict = SurfaceFeatureExtractDict(case = mycase,
                                                              stlname = mystl)
        surfaceFeatureExtractDict.writeFile()
        
        print('Create constant folder input files')
        consfolder = os.path.join(mycase,"constant")
        polyfolder = os.path.join(mycase,"constant","polyMesh")
        trifolder = os.path.join(mycase,"constant","triSurface")
        if not os.path.exists(consfolder): os.makedirs(consfolder)
        if not os.path.exists(polyfolder): os.makedirs(polyfolder)
        if not os.path.exists(trifolder): os.makedirs(trifolder)
        
        #blockMeshDict
        blockMeshDict = BlockMeshDict( case     = mycase,
                                       ndim     = 3,
                                       xmin     = param.xBounds[0]*Length,
                                       xmax     = param.xBounds[1]*Length,
                                       ymin     = param.yBounds[0]*Beam*0.5*(not param.symmetry),
                                       ymax     = param.yBounds[1]*Beam*0.5,
                                       zmin     = param.zBounds[0]*Depth,
                                       zmax     = param.zBounds[1]*Depth,
                                       fsmin    = param.fsBounds[0]*Depth,
                                       fsmax    = param.fsBounds[1]*Depth,
                                       sym      = param.symmetry,
                                       gridlvl  = grid,
                                       ofp      = param.OFplus)
        blockMeshDict.writeFile()
    
        print('Create 0 folder input files')
        zerofolder = os.path.join(mycase,"0")
        if not os.path.exists(zerofolder): os.makedirs(zerofolder)
    
        #Allrun
        print('Create run scripts')
        arun = os.path.join(mycase,'Allrun')
        with open(arun,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')
            f.write('function clearMeshInfo()\n')
            f.write('{\n')
            f.write('    rm -fr background.msh 0/{ccx,ccy,ccz,*Level,polyMesh/cellMap}* constant/polyMesh/{sets,*Level*,*Zone*,refine*,surfaceIndex*}\n')
            f.write('}\n\n')
            f.write('function refineBox()\n')
            f.write('{\n')
            f.write('    echo "cellSet c0 new boxToCell $1" | setSet -latestTime\n')
            f.write('    eval BVrefineMesh -dict "system/refineMeshDict.xyz"\n')
            f.write('}\n\n')
            f.write('function refineFS()\n')
            f.write('{\n')
            f.write('    REFLEVEL=$1\n')
            f.write('    PARAM=$2\n')
            f.write('    shift; shift\n')
            f.write('    echo "cellSet c0 new boxToCell $PARAM" | setSet -latestTime\n')
            f.write('    eval BVrefineMesh -dict "system/refineMeshDict.z" -refineUpToLevel $REFLEVEL $@\n')
            f.write('}\n\n')
            f.write('function snap()\n')
            f.write('{\n')
            f.write('    snappyHexMesh\n')
            f.write('}\n\n')
            f.write('(\n')
            f.write('    blockMesh\n')
            n = int(min([len(param.xRefineBox)/2, len(param.yRefineBox),len(param.zRefineBox)/2]))
            if n < param.nRefBoxes:
                print('ERROR: Requested number of boxes is greater than refinement level provided!')
                os._exit(1)
            for i in range(param.nRefBoxes):
                f.write('    refineBox "({:.2f} {:.2f} {:.2f}) ({:.2f} {:.2f} {:.2f})"\n'.format(param.xRefineBox[i]*Length,param.yRefineBox[i]*Beam*0.5,param.zRefineBox[i]*Depth,param.xRefineBox[i+n]*Length,param.yRefineBox[i+n]*Beam*0.5,param.zRefineBox[i+n]*Depth))
            n = int(len(param.fsRefineBox)/2)
            if n < param.nfsRefBoxes:
                print('ERROR: Requested number of boxes is greater than refinement level provided!')
                os._exit(1)
            for i in range(param.nfsRefBoxes):
                if i < param.nfsRefBoxes-1:
                    f.write('    refineFS {} "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(param.nRefBoxes,param.fsRefineBox[i]*Depth,param.fsRefineBox[i+n]*Depth))
                elif i == param.nfsRefBoxes-1:
                    f.write('    refineFS {} "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})" -noCellLevel\n'.format(param.nRefBoxes,param.fsRefineBox[i]*Depth,param.fsRefineBox[i+n]*Depth))
            f.write('    snap\n')
            f.write(') 2>&1 | tee log.mesh\n\n')
            f.write('checkMesh -latestTime 2>&1 | tee log.checkMesh\n')
        os.chmod(arun, 0o755)
    
        #Allclean
        aclean = os.path.join(mycase,'Allclean')
        with open(aclean,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('cd ${0%/*} || exit 1    # run from this directory\n\n')
            f.write('# Source tutorial clean functions\n')
            f.write('. $WM_PROJECT_DIR/bin/tools/CleanFunctions\n\n')
            f.write('cleanCase\n\n')
            f.write('rm -fr 0/alpha.water* 0/U* 0/p_rgh*\n')
        os.chmod(aclean, 0o755)
        
        print('Copy STL file, rotate and run surfaceFeatureExtract')
        #Copy STL in triSurface folder
        fstlin  = os.path.join(r'stl/'+mystl+'.stl')
        fstlout = os.path.join(trifolder,mystl+'.stl')
        shutil.copyfile(fstlin,fstlout)
        
        #rotate section if needed
        if abs(param.rot)>1e-3:
            fstltmp = os.path.join(trifolder,mystl+'_tmp.stl')
            shutil.move(fstlout,fstltmp)
            p = Popen(["surfaceTransformPoints -yawPitchRoll '(0 0 {:.3f})' {} {}".format(param.rot,fstltmp,fstlout)], shell=True)
            p.wait()
            os.remove(fstltmp)
        
        #run surfaceFeatureExtract
        p = Popen(['surfaceFeatureExtract'], cwd=mycase)
        p.wait()
        
        #run Allrun script
        if runScript:
            p = Popen(['./Allclean'], cwd=mycase)
            p.wait()
            p = Popen(['./Allrun'], cwd=mycase)
            p.wait()
        
        #create file for Paraview
        open(os.path.join(mycase,'a.foam'), 'a').close()
        
        
#*** Main execution start here *************************************************
    #Read input file
DEFAULT_PARAM = {'ndim'          : 2,
                 'sectionsFile'  : '',
                 'stlFile'       : '',
                 'sections'      : [],
                 'gridLevel'     : [1],
                 'rot'           : 0.0,
                 'symmetry'      : False,
                 'xBounds'       : [-0.5,0.5],
                 'yBounds'       : [-9.,9.],
                 'zBounds'       : [-2.,3.],
                 'fsBounds'      : [-0.5,1.5],
                 'nRefBoxes'     : 6,
                 'xRefineBox'    : [-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,1.6,1.5,1.4,1.3,1.2,1.1],
                 'yRefineBox'    : [-6.5,-4.5,-3.0,-2.0,-1.5,-1.2,6.5,4.5,3.0,2.0,1.5,1.2],
                 'zRefineBox'    : [-2.0,-2.0,-2.0,-1.5,-1.0,-0.5,3.0,2.5,2.0,1.8,1.6,1.4],
                 'nfsRefBoxes'   : 5,
                 'fsRefineBox'   : [-1.5,-1.0,-0.5,-0.3,-0.2,2.5,2.0,1.5,1.3,1.2],
                 'refineLength'  : [0.1],
		 'layerLength'   : 0.005,
                 'cellRatio'     : 1,
                 'OFplus'	 : False,
                 'gridName'      : None
                 }

DEFAULT_ARG = {'createStl' : True,
               'runScript' : True
              }

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

def fsMesher2D(userParam={}, userArgs={}):
    startTime = time.time()
    param = DEFAULT_PARAM
    param.update(userParam)
    param = Struct(**param)
    arg = DEFAULT_ARG
    arg.update(userArgs)
    arg =Struct(**arg)

    if abs(param.rot)>1e-3 and param.symmetry:
        print('ERROR: Mesh rotation cannot be applied with symmetry')
        os._exit(1)
    
    if param.ndim==2:
        sectDict = readSections(param.sectionsFile,sections=param.sections)
        if arg.createStl: createSectionStl(sectDict)
        create2DMesh(param,sectDict,arg.runScript)
    elif param.ndim==3:
        create3DMesh(param,arg.runScript)

    endTime = time.time()
    print('Mesh generation completed in '+str(datetime.timedelta(seconds=(endTime-startTime))))
