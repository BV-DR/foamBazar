#!/usr/bin/env python

#########################################################################
# Filename: fsTemplate2D.py                                             #
# Date:     2018-May-14                                                 #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import re
import os, sys
import shutil
import time, datetime
import math as mt
import numpy as np
import pandas as pd
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from fsTools import *
from subprocess import call, Popen
from scipy import interpolate as interp

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.dynamicMeshDict import DynamicMeshDict
from inputFiles.waveProperties import WaveProperties, RelaxZone, WaveCondition
from inputFiles.boundaryCondition import writeAllBoundaries
from inputFiles.turbulenceProperties import writeTurbulenceProperties
from inputFiles.transportProperties import TransportProperties
from inputFiles.gravity import Gravity

   
def setBoundaries(param):
    mycase = os.path.join(param.case,'grid_'+str(param.gridLevel))
    
    boundfile = os.path.join(mycase,"constant","polyMesh","boundary")
    boundDict = ParsedParameterFile(boundfile,boundaryDict=True)
    
    nbound = int(len(boundDict)/2)
    for i in range(nbound):
        if boundDict[2*i] in ['domainX0','domainX1']:
            boundDict[2*i+1]['type'] = 'empty'
        elif param.symmetry and (boundDict[2*i] in ['domainY0']):
            boundDict[2*i+1]['type'] = 'symmetryPlane'
    
    boundDict.writeFile()

def copyMesh(param,overwrite=False):
    mycase = os.path.join(param.case,'grid_'+str(param.gridLevel))
    mymesh = os.path.join(param.meshDir,param.case,'grid_'+str(param.gridLevel))
    
    if os.path.exists(mycase):
        if overwrite:
            shutil.rmtree(mycase)
        else:
            valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
            res = input('Case "{}" already exists, do you want to overwrite ? (y/n) '.format(mycase)).lower()
            if valid.get(res,False):
                shutil.rmtree(mycase)
            else:
                print('Exiting')
                os._exit(1)
    
    os.makedirs(mycase)
    
    consfolder = os.path.join(mycase,"constant")

    if not os.path.exists(consfolder): os.makedirs(consfolder)

    if param.meshTime=='latestTime':
        timeFolders = getFoamTimeFolders(mymesh)
        meshTimeFolder = timeFolders[-1]
    elif param.meshTime=='constant':
        meshTimeFolder = 'constant'
    else:
        meshTimeFolder = param.meshTime
    
    print('Copy mesh from folder '+meshTimeFolder)
    shutil.copytree(mymesh+r'/'+meshTimeFolder+r'/polyMesh',mycase+r'/constant/polyMesh')
    shutil.copytree(mymesh+r'/constant/triSurface',mycase+r'/constant/triSurface')
    
    setBoundaries(param)
    
def create2DCase(param,runScript=True):
    mycase = os.path.join(param.case,'grid_'+str(param.gridLevel))
    
    print('Create system folder input files')
    sysfolder = os.path.join(mycase,"system")
    if not os.path.exists(sysfolder): os.makedirs(sysfolder)
    
    #controlDict
    if param.outputForces:
        if len(param.hullPatch)>0: forcesPatch = param.hullPatch
        else: forcesPatch = param.case
    else: forcesPatch = None
    controlDict = ControlDict( case                = mycase,
                               startFrom           = param.startTime,
                               endTime             = param.endTime,
                               deltaT              = param.timeStep,
                               writeInterval       = param.writeInterval,
                               purgeWrite          = param.purgeWrite,
                               writePrecision      = 15,
                               forcesPatch         = forcesPatch,
                               rhoWater            = 1025,
                               OFversion           = param.OFversion,
                               version             = "foamStar" )
    controlDict.writeFile()
    
    #fvSchemes
    fvSchemes = FvSchemes( case     = mycase,
                           simType  = param.scheme,
                           orthogonalCorrection = "implicit",
                           version  = "foamStar" )
    fvSchemes.writeFile()
    
    #fvSolution
    fvSolution = FvSolution( case    = mycase,
                             useEuler = param.scheme=='Euler',
                             version  = "foamStar" )
    fvSolution.writeFile()
        
    #decomposeParDict
    decomposeParDict = DecomposeParDict( case   = mycase,
                                         nProcs = param.nProcs )
    decomposeParDict.writeFile()
            
    print('Create constant folder input files')
    #waveProperties
    filename = os.path.join(mycase,'constant','waveProperties')
    waveCond  = WaveCondition( waveType   = param.wave )
    
    relaxZones = []
    if param.sideRelaxZone is not None:
        bBox = findCFDBoundingBox(mycase,False)
        if param.sideRelaxZone>0:
            relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0. , -1. , 0.], bound=param.sideRelaxZone)
        else:
            relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0. , -1. , 0.], bound=0.5*bBox[4])
        relaxZones += [relaxSide]

    waveProperties = WaveProperties( filename,
                                     initWaveCondition = waveCond,
                                     relaxZones        = relaxZones,
                                     version           = "foamStar" )
    waveProperties.writeFile()
    
    #cell.Set
    if param.sideRelaxZone is not None:
        filename = os.path.join(mycase,'cell.Set')
        waveProperties.writeBlendingZoneBatch(filename)

    #dynamicMeshDict
    shutil.copyfile(param.dispSignal,os.path.join(mycase,"dispSignal.dat"))
    dynamicMeshDict = DynamicMeshDict( case      = mycase,
                                       type      = 'solid',
                                       dispFile  = "dispSignal.dat",
                                       OFversion = param.OFversion,
                                       version   = "foamStar")
    dynamicMeshDict.writeFile()

    #g
    gravity = Gravity( case = mycase, g = param.gravity )
    gravity.writeFile()
    
    #turbulenceProperties, RASProperties
    writeTurbulenceProperties( mycase , "laminar" )
    
    #transportProperties
    transportProperties = TransportProperties( case = mycase,
                                               rhoWater = 1025,
                                               version = "foamStar")
    transportProperties.writeFile()
    
    print('Create 0 folder input files')
    zerofolder = os.path.join(mycase,"0")
    orgfolder = os.path.join(mycase,"0","org")
    if not os.path.exists(zerofolder): os.makedirs(zerofolder)
    if not os.path.exists(orgfolder): os.makedirs(orgfolder)
    
    
    #alpha water, p_rgh, U, pointDisplacement
    writeAllBoundaries( case  = mycase,
                        case2D = True,
                        symmetryPlane = param.symmetry,
                        struct = '"'+param.case+'|wetSurf"',
                        version = "foamStar" )
                        
    #Allrun
    print('Create run scripts')
    arun = os.path.join(mycase,'Allrun')
    with open(arun,'w') as f:
        f.write('#! /bin/bash\n')
        f.write('set -x\n\n')
        f.write('(\n')
        # f.write('    cp -r ../../snap/Grid1/constant/polyMesh/ constant/\n')
        # f.write('    cp constant/org/boundary constant/polyMesh/\n')
        if param.sideRelaxZone is not None:
            f.write('    setSet -batch cell.Set\n')
            f.write('    setsToZones -noFlipMap\n')
        f.write('    cp -rf 0/org/* 0/\n')
        if param.translateLength!=0.:
            f.write('    transformPoints -translate "( 0. 0. {})"\n'.format(param.translateLength))
        f.write('    decomposePar -cellDist\n')
        f.write('    mpirun -np {:d} initWaveField -parallel\n'.format(param.nProcs))
        f.write(') 2>&1 | tee log.init\n')
    os.chmod(arun, 0o755)
        
    #Allclean
    aclean = os.path.join(mycase,'Allclean')
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
        f.write('eval clean_log\n')
        f.write('eval clean_mesh\n')
        f.write('eval clean_parallel_mesh\n')
        f.write('eval clean_0\n')
    os.chmod(aclean, 0o755)
    
    #run Allrun script
    if runScript:
        p = Popen(['./Allclean'], cwd=mycase)
        p.wait()
        p = Popen(['./Allrun'], cwd=mycase)
        p.wait()
    
    #create file for Paraview
    open(os.path.join(mycase,'a.foam'), 'a').close()
    
    #run.sh
    run = os.path.join(mycase,'run.sh')
    with open(run,'w') as f:
        f.write('#!/bin/bash -l\n')
        f.write('#SBATCH -J {}\n\n'.format(mycase))
        f.write('# 5 hour wall-clock\n')
        f.write('#SBATCH --account I1608251\n')
        f.write('#SBATCH -t 3-00:00:00\n')
        f.write('#SBATCH -n {:d}\n'.format(param.nProcs))
        f.write('#SBATCH -o log.run-%j\n\n')
        f.write('module load gcc/4.9.3 openmpi/1.8.4-gcc lapack/3.6.1/gcc/4.9.3\n')
        f.write('export FOAM_INST_DIR=/data/I1608251/OpenFOAM;\n')
        if   param.OFversion==2: f.write('source /data/I1608251/OpenFOAM/OpenFOAM-2.4.x/etc/bashrc;\n')
        elif param.OFversion==3: f.write('source /data/I1608251/OpenFOAM/OpenFOAM-3.0.x/etc/bashrc;\n')
        elif param.OFversion==5: f.write('source /data/I1608251/OpenFOAM/OpenFOAM-5.x/etc/bashrc;\n')
        f.write('export LC_ALL=C\n\n')
        f.write('mpirun {} -parallel\n'.format(param.solver))
    
def tarCase(param):
    mycase = os.path.join(param.case,'grid_'+str(param.gridLevel))
    print('Creating archive {}.tar.gz'.format(mycase))
    process = 'tar czf {}.tar.gz {}'.format(os.path.join(param.case,param.case+'_g'+str(param.gridLevel)),mycase)
    subprocess.call(process, shell=True)
        
#*** Main execution start here *************************************************
    #Read input file
DEFAULT_PARAM = {'case'             : 'newCase',
                 'meshDir'          : 'mesh',
                 'meshTime'         : 'constant',
                 'gridLevel'        : [1],
                 'symmetry'         : False,
                 'outputForces'     : False,
                 'hullPatch'        : '',
                 'startTime'        : 'latestTime',
                 'endTime'          : 10,
                 'timeStep'         : 0.01,
                 'writeInterval'    : 1,
                 'purgeWrite'       : 0,
                 'scheme'           : 'Euler',
                 'nProcs'           : 4,
                 'nOuterCorrectors' : 5,
                 'wave'             : "noWaves",
                 'waveH'            : 0.0,
                 'waveT'            : 0.0,
                 'velocity'         : 0.0,
                 'depth'            : 100.,
                 'sideRelaxZone'    : None,
                 'dispSignal'       : None,
                 'solver'           : 'foamStar',
                 'OFversion'        : 3,
                 'translateLength'  : 0.0,
                 'gravity'          : 9.81
                }
                 
DEFAULT_ARGS = {'overwrite' : False,
                'runScript' : True,
                'tar'       : False
               }

class Struct:
    def __init__(self, **entries):
        self.__dict__.update(entries)

def fsTemplate2D(userParam={}, userArgs={}):
    startTime = time.time()

    param = DEFAULT_PARAM
    param.update(userParam)
    param = Struct(**param)
    arg = DEFAULT_ARGS
    arg.update(userArgs)
    arg = Struct(**arg)

    #mycase = os.path.join(param.case,'Grid_'+str(param.gridLevel))
    copyMesh(param,overwrite=arg.overwrite)
    create2DCase(param,runScript=arg.runScript)
    if arg.tar: tarCase(param)
    endTime = time.time()
    print('Mesh generation completed in '+str(datetime.timedelta(seconds=(endTime-startTime))))
