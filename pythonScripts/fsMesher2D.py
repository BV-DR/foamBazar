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
from Pluto.System.readInput import readInput

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.extrudeMeshDict import ExtrudeMeshDict
from inputFiles.refineMeshDict import RefineMeshDict
from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.blockMeshDict import BlockMeshDict

def cmdOptions(argv):
    global DEBUG
    # default global parameters
#    CMD_showLog = 'echo "\nlog-file: ./log.template"; tail -45 ./log.tempalte; echo "Please see log-file: ./log.template"'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='inputfile', help='read parameters from input file. Use this option to edit redefine parameters, e.g.: fsMesher2D.py -f myship.cfg ')
    parser.add_argument('-nostl', dest='createStl', action='store_false', help='disable STL files generation from Homer inputs')
    parser.add_argument('-norun', dest='runScript', action='store_false', help='disable run of Allrun script')
    args = parser.parse_args()
      
    if args.inputfile==None:
        call('echo "No option provided, please une -h for help"', shell=True)
        raise SystemExit('')
                
    return args


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
        
        #create symetric
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
        
def create2DMesh(sdict,gridlvl,symmetry,runScript):
    for isect in sdict.keys():
        print('')
        print('Section =',isect)
        
        mysect = 'section_'+str(isect)
        
        for grid in gridlvl:
        
            mygrid = 'grid_'+str(gridlvl)
            mycase = os.path.join(mysect,mygrid)
            length = sdict[isect].y.max()
            height = sdict[isect].z.max()
            
            if not os.path.exists(mysect): os.makedirs(mysect)
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
                                                refinementLength    = length/15.,   #mesh length/15
                                                nSurfaceLayers      = 3,
                                                expansionRatio      = 1.3,
                                                finalLayerThickness = length/200.,  #mesh length/200
                                                minThickness        = length/2000.) #mesh length/2000
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
            blockMeshDict = BlockMeshDict( case     = mycase,
                                        ymin     = -9.0*length*(not symmetry),
                                        ymax     =  9.0*length,
                                        zmin     = -2.0*height,
                                        zmax     =  3.0*height,
                                        fsmin    = -0.5*height,
                                        fsmax    =  1.5*height,
                                        sym     = symmetry,
                                        gridlvl = gridlvl )
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
                f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(-6.5*length,-2.0*height,6.5*length,3.0*height))
                f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(-4.5*length,-2.0*height,4.5*length,2.5*height))
                f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(-3.0*length,-2.0*height,3.0*length,2.0*height))
                f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(-2.0*length,-1.5*height,2.0*length,1.8*height))
                f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(-1.5*length,-1.0*height,1.5*length,1.6*height))
                f.write('    refineBox "(-1e6 {:.2f} {:.2f}) (1e6 {:.2f} {:.2f})"\n'.format(-1.1*length,-0.5*height,1.1*length,1.4*height))
                f.write('    refineFS 6 "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(-1.5*height,2.5*height))
                f.write('    refineFS 6 "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(-1.0*height,2.0*height))
                f.write('    refineFS 6 "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(-0.5*height,1.5*height))
                f.write('    refineFS 6 "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})"\n'.format(-0.3*height,1.3*height))
                f.write('    refineFS 6 "(-1e6 -1e6 {:.2f}) (1e6 1e6 {:.2f})" -noCellLevel\n'.format(-0.2*height,1.2*height))
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
            
            print('Copy STL file and run surfaceFeatureExtract')
            #Copy STL in triSurface folder
            fstlin  = os.path.join(r'stl/',mysect+'.stl')
            fstlout = os.path.join(trifolder,mysect+'.stl')
            shutil.copyfile(fstlin,fstlout)
            
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
DEFAULT_PARAM = {'sectionsFile'  : '',
                 'sections'      : [],
                 'gridLevel'     : [1],
                 'symmetry'      : False }
                 
iParam = ['gridLevel']
rParam = []
lParam = ['symmetry']
sParam = ['sectionsFile']
aParam = ['sections']

if __name__ == "__main__":
    startTime = time.time()    
    arg = cmdOptions(sys.argv)
    param = readInput(arg.inputfile,iParam=iParam,rParam=rParam,lParam=lParam,sParam=sParam,aParam=aParam,name='fsMesher2D',default=DEFAULT_PARAM)
    sectDict = readSections(param.sectionsFile,sections=param.sections)
    if arg.createStl: createSectionStl(sectDict)
    create2DMesh(sectDict,param.gridLevel,param.symmetry,arg.runScript)

    endTime = time.time()
    print('Mesh generation completed in '+str(datetime.timedelta(seconds=(endTime-startTime))))