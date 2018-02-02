#!/usr/bin/env python

#########################################################################
# Filename: fsTemplate.py                                               #
# Date:     2017-May-02                                                 #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import os, subprocess, time, math
import numpy as np
import sys, argparse, configparser
import pprint
from fsTools import *
from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.waveProperties import WaveProperties, RelaxZone, WaveCondition
from inputFiles.waveProbes import createLinearWaveProbesList, setWaveProbes
from inputFiles.boundaryCondition import writeAllBoundaries
from inputFiles.turbulenceProperties import writeTurbulenceProperties
from inputFiles.transportProperties import TransportProperties
from inputFiles.sixDofDomainBody import SixDofDomainBody
from inputFiles.dynamicMeshDict import DynamicMeshDict
from inputFiles.flexProperties import writeFlexProperties
from inputFiles.gravity import Gravity
from inputFiles.decomposeParDict import DecomposeParDict

# abenhamou: 2017-may-16

DEBUG = False
DEFAULT_CFG_FILE = 'template.cfg'
DEFAULT_PARAMS = {
'caseDir' : 'newCase',
'meshDir' : 'mesh',
'meshTime' : 'latestTime',
'stlFile' : 'ship.stl',
'nProcs' : 4,
'solver' : 'foamStar',
'scheme' : 'Euler',
'startTime' : 'latestTime',
'endTime' : 1000,
'timeStep' : 0.01,
'writeInterval' : 1,
'purgeWrite' : 0,
'outputVBM' : False,
'outputWave' : False,
'outputMotions' : True,
'outputForces' : False,
'fsiTol' : 1e-6,
'hullPatch' : 'ship',
'donFile' : 'ship.don',
'wave' : 'noWaves',
'waveH' : 1.0,
'waveT' : 10.0,
'velocity' : 0., #in m/s
'depth' : 500,
'draft' : 0,
'waveSeaLvl' : 0,
'waveStartTime' : 0,
'waveRampTime' : 0,
'addDamping' : False,
'vtkOut' : True,
'localMotionPts' : [],
'mdFile' : None,
'modesToUse' : '',
'datFile' : None,
'dmigFile' : None,
'hmrUserOutput' : None,
'hmrScaling' : 0,
'EulerCellsDist' : 8,
'inletRelaxZone' : None,
'outletRelaxZone' : None,
'sideRelaxZone' : None,
'shipDamping' : None,
'waveProbes' : [],
'case2D' : False
}

#*** Data structure for user input data
class UserInput(object):
    def __init__(self, dir='', opts=DEFAULT_PARAMS, genCase=False, initCase=False, decomposeCase=False, disableBound=False, tar=False):
    
        if dir=='':
            self.caseDir =  opts['caseDir']
        else:
            self.caseDir =  dir
        self.meshDir        =  opts['meshDir']
        self.meshTime       =  opts['meshTime']
        self.stlFile        =  opts['stlFile']
        self.nProcs         =  opts['nProcs']
        self.solver         =  opts['solver']
        self.scheme         =  opts['scheme']
        self.startTime      =  opts['startTime']
        self.endTime        =  opts['endTime']
        self.timeStep       =  opts['timeStep']
        self.writeInterval  =  opts['writeInterval']
        self.purgeWrite     =  opts['purgeWrite']
        self.outputVBM      =  opts['outputVBM']
        self.outputWave     =  opts['outputWave']
        self.outputMotions  =  opts['outputMotions']
        self.outputForces   =  opts['outputForces']
        self.fsiTol         =  opts['fsiTol']
        self.hullPatch      =  opts['hullPatch']
        self.donFile        =  opts['donFile']
        self.donName        =  (self.donFile.rpartition('/')[-1]).partition('.don')[0]
        self.wave           =  opts['wave']
        self.waveH          =  opts['waveH']
        self.waveT          =  opts['waveT']
        self.velocity       =  opts['velocity']
        self.depth          =  opts['depth']
        self.draft          =  opts['draft']
        self.waveSeaLvl     =  opts['waveSeaLvl']
        self.waveStartTime  =  opts['waveStartTime']
        self.waveRampTime   =  opts['waveRampTime']
        self.addDamping     =  opts['addDamping']
        self.vtkOut         =  opts['vtkOut']
        self.localMotionPts =  opts['localMotionPts']
        self.waveProbes     =  opts['waveProbes']
        self.genCase        =  genCase
        self.initCase       =  initCase
        self.decomposeCase  =  decomposeCase
        self.disableBound   =  disableBound
        self.EulerCellsDist =  opts['EulerCellsDist']
        self.inletRelaxZone =  opts['inletRelaxZone']
        self.outletRelaxZone=  opts['outletRelaxZone']
        self.sideRelaxZone  =  opts['sideRelaxZone']
        self.mdFile         =  opts['mdFile']
        self.modesToUse     =  opts['modesToUse']
        self.datFile        =  opts['datFile']
        self.dmigFile       =  opts['dmigFile']
        self.hmrUserOutput  =  opts['hmrUserOutput']
        self.shipMass       =  0.
        self.hmrScaling     =  0.
        self.shipFreq       =  ''
        self.shipInertia    =  ''
        self.shipCOG        =  ''
        self.shipDamping    =  opts['shipDamping']
        self.case2D         =  opts['case2D']
        self.tar            =  tar
   
def cmdOptions(argv):
    global DEBUG
    # default global parameters
    CMD_showLog = 'echo "\nlog-file: ./log.template"; tail -45 ./log.tempalte; echo "Please see log-file: ./log.template"'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='case', help='if used without input file.  Use this option to define new case folder name, default parameters will be used for other purpose, e.g.: template.py -c newCase')
    parser.add_argument('-f', dest='inputfile', help='read parameters from input file. Use this option to edit redefine parameters, e.g.: template.py -f myship.cfg ')
    parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print all available input parameters (including theirs default values) to template.cfg file and exit. Use this option to generate a template for input files for -f')
    parser.add_argument('-i','-init', dest='initCase', nargs='?', default='void', help='enables automatic case initialisation (add <directory> if not used with -f or -c option)')
    parser.add_argument('-d', dest='decomposeCase', action='store_true', help='enables case decomposition for parallel run')
    parser.add_argument('-n', dest='disableBound', action='store_true', help='disable boundaries editing of mesh folder')
    parser.add_argument('-tar', dest='tar', action='store_true', help='activate compression of resulting folder')
    args = parser.parse_args()
    
    if args.showConfig:
        print "\nOutput default parameters to file:", DEFAULT_CFG_FILE
        subprocess.call('cat << EOF > '+DEFAULT_CFG_FILE + defaultParams_contents, shell=True)
        raise SystemExit('')
    
    if args.inputfile==None:
        if args.case==None:
            if args.initCase==None:
                subprocess.call('echo "Please provide case folder to initialize : -i <caseFolder>"', shell=True)
                raise SystemExit('')
            elif args.initCase=='void':
                subprocess.call('echo "No option provided, please une -h for help"', shell=True)
                raise SystemExit('')
            else:
                subprocess.call('echo "Initialising FoamStar case for case '+args.initCase+'"', shell=True)
                inputdata = UserInput(dir=str(args.initCase),initCase=True,decomposeCase=args.decomposeCase,disableBound=args.disableBound,tar=args.tar)
        else:
            if args.initCase=='void':
                inputdata = UserInput(dir=str(args.case),genCase=True,disableBound=args.disableBound,tar=args.tar)
            else:
                inputdata = UserInput(dir=str(args.case),genCase=True,initCase=True,decomposeCase=args.decomposeCase,disableBound=args.disableBound,tar=args.tar)
    else:
        # here we read all input parameters from files
        fid = str(args.inputfile)
        subprocess.call('echo "Input data are read from file '+fid+'"', shell=True)
        params = readInputParams(fid)
        if args.initCase=='void':
            inputdata = UserInput(opts=params,genCase=True,disableBound=args.disableBound,tar=args.tar)
        else:
            inputdata = UserInput(opts=params,genCase=True,initCase=True,decomposeCase=args.decomposeCase,disableBound=args.disableBound,tar=args.tar)
        pass
        
    return inputdata
    pass

def readInputParams(filename):
    params = dict(DEFAULT_PARAMS)
    config = configparser.ConfigParser()
    config.read(filename)
    name = 'template'
    
    print "Read input paramters from file:", filename
    
    # caseDir
    try:
        txt = str(config[name]['caseDir'])
        params['caseDir'] = txt
    except KeyError: pass
    
    # meshDir
    try:
        txt = str(config[name]['meshDir'])
        params['meshDir'] = txt
    except KeyError: pass
    
    # meshTime
    try:
        txt = str(config[name]['meshTime'])
        params['meshTime'] = txt
    except KeyError: pass
    
    # stl file
    try:
        txt = str(config[name]['stlFile'])
        params['stlFile'] = txt
    except KeyError: pass
    
    # nProcs
    try:
        txt = str(config[name]['nProcs'])
        params['nProcs'] = int(txt)
    except KeyError: pass
    
    # heading
    try:
        txt = str(config[name]['heading'])
        params['heading'] = float(txt) % float(360)
    except KeyError: pass

    # 'solver'
    try:
        txt = str(config[name]['solver'])
        params['solver'] = txt
    except KeyError: pass
    
    # 'scheme' : 'Euler', 'CrankNicolson 0.9'
    try:
        txt = str(config[name]['scheme'])
        params['scheme'] = txt
    except KeyError: pass
    
    # start time
    try:
        txt = str(config[name]['startTime'])
        params['startTime'] = txt
    except KeyError: pass
    
    # end time
    try:
        txt = str(config[name]['endTime'])
        params['endTime'] = float(txt)
    except KeyError: pass
    
    # time step
    try:
        txt = str(config[name]['timeStep'])
        params['timeStep'] = float(txt)
    except KeyError: pass
    
    # write interval
    try:
        txt = str(config[name]['writeInterval'])
        params['writeInterval'] = int(txt)
    except KeyError: pass
    
    # purge write
    try:
        txt = str(config[name]['purgeWrite'])
        params['purgeWrite'] = int(txt)
    except KeyError: pass
    
    # 'outputMotions' : output Motions or not?
    try:
        txt = str(config[name]['outputMotions'])
        params['outputMotions'] = getBool(txt)
    except KeyError: pass
    
    # 'outputForces' : output Forces or not?
    try:
        txt = str(config[name]['outputForces'])
        params['outputForces'] = getBool(txt)
    except KeyError: pass
    
    # 'outputVBM' : output VBM or not?
    try:
        txt = str(config[name]['outputVBM'])
        params['outputVBM'] = getBool(txt)
    except KeyError: pass
    
    # 'outputWave' : output Waves or not?
    try:
        txt = str(config[name]['outputWave'])
        params['outputWave'] = getBool(txt)
    except KeyError: pass
    
    # time step
    try:
        txt = str(config[name]['fsiTol'])
        params['fsiTol'] = float(txt)
    except KeyError: pass

    # hull patch name
    try:
        txt = str(config[name]['hullPatch'])
        params['hullPatch'] = txt
    except KeyError: pass
    
    # DON file name
    try:
        txt = str(config[name]['donFile'])
        params['donFile'] = txt
    except KeyError: pass

    # wave type
    try:
        txt = str(config[name]['wave'])
        params['wave'] = txt
    except KeyError: pass
    
    # wave height
    try:
        txt = str(config[name]['waveH'])
        params['waveH'] = float(txt)
    except KeyError: pass
    
    # wave period
    try:
        txt = str(config[name]['waveT'])
        params['waveT'] = float(txt)
    except KeyError: pass

    # velocity
    try:
        txt = str(config[name]['velocity'])
        params['velocity'] = float(txt)
    except KeyError: pass
    
    # depth
    try:
        txt = str(config[name]['depth'])
        params['depth'] = float(txt)
    except KeyError: pass
    
    # wave Start Time
    try:
        txt = str(config[name]['waveStartTime'])
        params['waveStartTime'] = float(txt)
    except KeyError: pass
    
    # wave Rmp Time
    try:
        txt = str(config[name]['waveRampTime'])
        params['waveRampTime'] = float(txt)
    except KeyError: pass

    # artificial damping for still water case
    try:
        txt = str(config[name]['addDamping'])
        params['addDamping'] = getBool(txt)
    except KeyError: pass
    
    # template for 2D case
    try:
        txt = str(config[name]['case2D'])
        params['case2D'] = getBool(txt)
    except KeyError: pass
    
    # output VTK files for modes shapes
    try:
        txt = str(config[name]['vtkOut'])
        params['vtkOut'] = getBool(txt)
    except KeyError: pass
    
    # list of local motion points
    try:
        txt = str(config[name]['localMotionPts'])
        if len(txt)>0:
            try:
                ptm = np.matrix(txt)
            except:
                sys.exit('ERROR: Unable to read list of local points. Please use the follwing format : \n \
                                 >>> localMotionPts = 1 2 3; 4 5 6; 7 8 9')
            if np.shape(ptm)[1]==3:
                params['localMotionPts'] = np.array(ptm)
            else:
                sys.exit('ERROR: Points for local motions must have 3 coordinates')
    except KeyError: pass
    
    # list of wave probes
    try:
        txt = str(config[name]['waveProbes'])
        if len(txt)>0:
            try:
                ptm = np.matrix(txt)
            except:
                sys.exit('ERROR: Unable to read list of waveProbes. Please use the follwing format : \n \
                                 >>> waveProbes = xMin1 xMax1 nX1 y1 zMin1 Zmax1 nZ1; xMin2 xMax2 nX2 y2 zMin2 Zmax2 nZ2')
            if np.shape(ptm)[1]==7:
                params['waveProbes'] = np.array(ptm)
            else:
                sys.exit('ERROR: Definition of waveProbes shoud include 7 parameters')
    except KeyError: pass
    
    # sea level
    try:
        txt = str(config[name]['waveSeaLvl'])
        params['waveSeaLvl'] = float(txt)
    except KeyError: pass
  
    # draft
    try:
        txt = str(config[name]['draft'])
        if not txt=='None':
            params['draft'] = math.fabs(float(txt))
    except KeyError: pass
    
    # EulerCellsDist
    try:
        txt = str(config[name]['EulerCellsDist'])
        params['EulerCellsDist'] = float(txt)
    except KeyError: pass
    
    # inletRelaxZone
    try:
        txt = str(config[name]['inletRelaxZone'])
        params['inletRelaxZone'] = float(txt)
    except KeyError: pass
    
    # outletRelaxZone
    try:
        txt = str(config[name]['outletRelaxZone'])
        params['outletRelaxZone'] = float(txt)
    except KeyError: pass
    
    # sideRelaxZone
    try:
        txt = str(config[name]['sideRelaxZone'])
        params['sideRelaxZone'] = float(txt)
    except KeyError: pass
    
    # mdFile
    try:
        txt = str(config[name]['mdFile'])
        params['mdFile'] = txt
    except KeyError: pass

    # modesToUse
    try:
        txt = str(config[name]['modesToUse'])
        params['modesToUse'] = txt
    except KeyError: pass
    
    # datFile
    try:
        txt = str(config[name]['datFile'])
        params['datFile'] = txt
    except KeyError: pass
    
    # dmigFile
    try:
        txt = str(config[name]['dmigFile'])
        params['dmigFile'] = txt
    except KeyError: pass
    
    # hmrUserOutput
    try:
        txt = str(config[name]['hmrUserOutput'])
        params['hmrUserOutput'] = txt
    except KeyError: pass
    
    # shipDamping
    try:
        txt = str(config[name]['shipDamping'])
        params['shipDamping'] = txt
    except KeyError: pass
    
    return params
    pass
    
def run_setSet():
    # FIXME: refineMesh is buggy when run in parallel
    runCommand(CMD_setSet + ' -batch .tmp_setSet ' + CMD_keepLog)

    #if NPROCS==1:
    #    runCommand(CMD_setSet + ' -batch .tmp_setSet ' + CMD_keepLog)
    #else:
    #    runCommand('mpirun -np '+str(NPROCS)+' '+ CMD_setSet + ' -parallel -batch .tmp_setSet ' + CMD_keepLog)
    pass

def readHomerData(data):
    print 'Reading data from Homer file : '+data.hmrUserOutput
    modes2use = [int(i) for i in data.modesToUse.split()]
    modes2use.sort()
    
    freq = []
    scaling = []
    inertia = [0.]*6
    
    f = open(data.hmrUserOutput,'r')
    itf = iter(f)
    for line in itf:
        if 'Location of center of gravity in global reference' in line:
            next(itf)
            pline = next(itf)
            cog = [float(i) for i in pline.split()]
            
        elif 'Mass =' in line:
            data.shipMass = float(line.split()[2])
            
        elif 'Roll Inertia =' in line:
            inertia[0] = float(line.split()[3])
            
        elif 'Pitch Inertia =' in line:
            inertia[1] = float(line.split()[3])
        
        elif 'Yaw Inertia =' in line:
            inertia[2] = float(line.split()[3])
            break
        elif len(modes2use)>0:
            for mode in modes2use:
                mstr = 'Mode '+str(mode)
                if mstr in line:
                    sline = line.split()
                    scaling.append(float(sline[3]))
                    freq.append(float(sline[6])**2)
    f.close()

    data.shipCOG    =  cog
    data.shipInertia = inertia
    if len(modes2use)>0:
        data.hmrScaling  = ('{:20.14e}').format(scaling[0])
        data.shipFreq    = ('{:20.14e} '*len(freq)).format(*freq)
    
def foamCase_template(data):
    print 'Create system folder input files'
    #controlDict
    if data.outputVBM:
        vbmPatch = [data.hullPatch,data.donName]
    else:
        vbmPatch = None
    if data.outputWave:
        # print 'WARNING: Wave probes cannot be included yet'
        wpList = []
        for wp in data.waveProbes:
            wpList += createLinearWaveProbesList(*wp)
    else:
        wpList = None
        
    if data.outputForces:
        forcesPatch = data.hullPatch
    else:
        forcesPatch = None

    controlDict = ControlDict( case                = data.caseDir,
                               startFrom           = data.startTime,
                               endTime             = data.endTime,
                               deltaT              = data.timeStep,
                               writeInterval       = data.writeInterval,
                               purgeWrite          = data.purgeWrite,
                               writePrecision      = 15,
                               outputMotions       = data.outputMotions,
                               outputLocalMotions  = len(data.localMotionPts)>0,
                               vbmPatch            = vbmPatch,
                               forcesPatch         = forcesPatch,
                               waveProbesList      = wpList,
                               version             = "foamStar" )
    controlDict.writeFile()

    #fvSchemes
    fvSchemes = FvSchemes( case     = data.caseDir,
                           simType  = data.scheme,
                           orthogonalCorrection = "implicit",
                           version  = "foamStar" )
    fvSchemes.writeFile()

    #fvSolution
    fvSolution = FvSolution( case    = data.caseDir,
                             fsiTol  = data.fsiTol,
                             useEuler = data.scheme=='Euler',
                             version  = "foamStar" )
    fvSolution.writeFile()
    
    #decomposeParDict
    decomposeParDict = DecomposeParDict( case=data.caseDir, nProcs=data.nProcs )
    decomposeParDict.writeFile()

    print 'Create constant folder input files'
    #waveProperties
    filename = data.caseDir+'/constant/waveProperties'
    waveCond  = WaveCondition( waveType   = data.wave,
                               height     = data.waveH,
                               period     = data.waveT,
                               U0         = -1*data.velocity,
                               depth      = data.depth,
                               refDirection = [-1,0,0],
                               startTime  = data.waveStartTime,
                               rampTime   = data.waveRampTime
                             )
    
    relaxZones = []
    bBox = findCFDBoundingBox(data.caseDir,False)
    if data.sideRelaxZone is not None:
        relaxSide   = RelaxZone( "side"  , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0. , -1. , 0.], bound=data.sideRelaxZone)
        relaxZones += [relaxSide]
    if data.outletRelaxZone is not None:
        relaxOutlet = RelaxZone( "outlet", relax=True, waveCondition=waveCond, origin=[bBox[0], 0., 0.], orientation = [  1. ,  0. , 0.], bound=data.outletRelaxZone)
        relaxZones += [relaxOutlet]
    if data.inletRelaxZone  is not None:
        relaxInlet  = RelaxZone( "inlet" , relax=True, waveCondition=waveCond, origin=[bBox[3], 0., 0.], orientation = [ -1. ,  0. , 0.], bound=data.inletRelaxZone)
        relaxZones += [relaxInlet]
    
    waveProperties = WaveProperties( filename, initWaveCondition=waveCond, relaxZones=relaxZones, version = "foamStar" )
    waveProperties.writeFile()
    
    #cell.Set
    filename = data.caseDir+'/cell.Set'
    waveProperties.writeBlendingZoneBatch(filename)
    if not data.case2D:
        bBox = findBoundingBox(data.caseDir+'/constant/triSurface/'+data.stlFile,False)
        bBox = [0.5*(bBox[0]+bBox[3]), 0.5*(bBox[1]+bBox[4])+math.fabs(bBox[1]-bBox[4]), 0.5*(bBox[2]+bBox[5])]
        outsidePoints = str(bBox[0])+" "+str(bBox[1])+" "+str(bBox[2])
        cSet = open(filename,'a')
        cSet.write('cellSet EulerCells new surfaceToCell "./constant/triSurface/'+data.stlFile+'" (('+outsidePoints+')) yes yes no '+str(data.EulerCellsDist)+' -1e6')
        cSet.close()
    
    #turbulenceProperties, RASProperties
    writeTurbulenceProperties( data.caseDir , "laminar" )
        
    #transportProperties
    transportProperties = TransportProperties( case = data.caseDir,
                                               nuWater = 1.14e-06,
                                               rhoWater = 1025,
                                               version = "foamStar"
                                             )
    transportProperties.writeFile()
    
    #dynamicMeshDict
    if data.case2D:
        dynamicMeshDict = DynamicMeshDict( case = data.caseDir,
                                           static = True,
                                           version = "foamStar"
                                         )
    else:
        if data.addDamping:
            bBox = findBoundingBox(data.caseDir+'/constant/triSurface/'+data.stlFile,False)
            LBP = round(0.95*abs(bBox[3]-bBox[0]),2)
            BC = round(abs(bBox[4]-bBox[1]),2)
        else:
            LBP = 0.
            BC = 0.
        dynamicMeshDict = DynamicMeshDict( case = data.caseDir,
                                        hullPatch = data.hullPatch,
                                        addDamping = data.addDamping,
                                        lpp = LBP,
                                        bc = BC,
                                        version = "foamStar"
                                        )
    dynamicMeshDict.writeFile()

    #g
    gravity = Gravity( case = data.caseDir )
    gravity.writeFile()
    
    print 'Create 0 folder input files'
    #alpha water, p_rgh, U, pointDisplacement
    writeAllBoundaries( case=data.caseDir,
                        speed=data.velocity,
                        case2D = data.case2D,
                        version = "foamStar" )
    
    #sixDofDomainBody, flexProperties
    if not data.case2D:
        sixDofDomainBody = SixDofDomainBody( case = data.caseDir,
                                            mass = data.shipMass,
                                            inertia = data.shipInertia,
                                            COG = data.shipCOG,
                                            nModes = data.modesToUse,
                                            donName = data.donName,
                                            version = "foamStar"
                                            )
        sixDofDomainBody.writeFile()
     
        if len(data.modesToUse)>0:
            print 'Create flexible properties input files'
            writeFlexProperties( case = data.caseDir,
                                 donName = data.donName,
                                 mdFile = data.mdFile,
                                 modes2use = data.modesToUse,
                                 datFile = data.datFile,
                                 dmigFile = data.dmigFile,
                                 draft = data.draft,
                                 scale = data.hmrScaling,
                                 vtkOut = data.vtkOut,
                                 hullPatch = data.hullPatch,
                                 localPts = data.localMotionPts,
                                 freq = data.shipFreq,
                                 damping = data.shipDamping,
                                 version = "foamStar"
                                )

    #run.sh
    filename = data.caseDir+'/run.sh'
    if not os.path.isfile('run.sh'):
        subprocess.call('cat << EOF > '+ filename + run_contents, shell=True)
        setValue(filename, 'CASE_NAME', data.caseDir)
        setValue(filename, 'NB_PROCS', data.nProcs)
        setValue(filename, 'SOLVER_NAME', data.solver)
        
def setBoundaries(data):
    rc = '\n'
    
    os.rename(data.caseDir+'/constant/polyMesh/boundary',data.caseDir+'/constant/polyMesh/boundary_old')
    obound = open(data.caseDir+'/constant/polyMesh/boundary_old','r')
    nbound = open(data.caseDir+'/constant/polyMesh/boundary','w')
    
    for line in obound:
        if (not data.case2D) and ('10' in line[:4]):
            nbound.write('7'+rc)
        elif data.case2D and '7' in line[:4]:
            nbound.write('6'+rc)
        elif 'defaultFaces' in line:
            for _ in xrange(5): obound.next()
        elif ('domainY0' in line) and data.case2D:
            nbound.write(line)
            line = obound.next()
            nbound.write(line)
            line = obound.next()
            nbound.write('        type            empty;'+rc)
            line = obound.next()
        elif 'domainY1' in line:
            nbound.write(line)
            line = obound.next()
            nbound.write(line)
            line = obound.next()
            if data.case2D:
                nbound.write('        type            empty;'+rc)
            else:
                nbound.write('        type            symmetryPlane;'+rc)
                nbound.write('        inGroups        1(symmetryPlane);'+rc)
        elif 'ship_hull' in line:
            while(not 'nFaces' in line): line = obound.next()
            nFaces = int(line.split()[-1][:-1])
            
            while(not 'startFace' in line): line = obound.next()
            startFace = int(line.split()[-1][:-1])
            
            while(not 'nFaces' in line): line = obound.next()
            nFaces += int(line.split()[-1][:-1])
            line = obound.next()
            
            while(not 'nFaces' in line): line = obound.next()
            nFaces += int(line.split()[-1][:-1])
            line = obound.next()
            line = obound.next()
            
            nbound.write('    ship'+rc)
            nbound.write('    {'+rc)
            nbound.write('        type            wall;'+rc)
            nbound.write('        inGroups        1(wall);'+rc)
            nbound.write('        nFaces          '+str(nFaces)+';'+rc)
            nbound.write('        startFace       '+str(startFace)+';'+rc)
            nbound.write('    }'+rc)
        else:
            nbound.write(line)
    
def copyMesh(data):

    subprocess.call('mkdir -p '+data.caseDir+'/{0/org,0/uniform,constant,system}', shell=True)

    if data.meshTime=='latestTime':
        timeFolders = getFoamTimeFolders(data.meshDir)
        meshTimeFolder = timeFolders[-1]
    else:
        meshTimeFolder = data.meshTime
    
    print 'Copy mesh from folder '+meshTimeFolder
    subprocess.call('cp -r '+data.meshDir+'/'+meshTimeFolder+'/polyMesh '+data.caseDir+'/constant', shell=True)
    
    if not data.disableBound: setBoundaries(data)
    
    if not data.case2D:
        print 'Copy '+data.stlFile+' file'
        subprocess.call('mkdir '+data.caseDir+'/constant/triSurface', shell=True)
        subprocess.call('cp -r '+data.meshDir+'/constant/triSurface/'+data.stlFile+' '+data.caseDir+'/constant/triSurface', shell=True)
        
        print 'Copy '+data.donName+'.don file'
        subprocess.call('cp '+data.donFile+' '+data.caseDir+'/.', shell=True)

def initCase(data):
    string=[]
    string.append('setSet -batch cell.Set')
    string.append('setsToZones -noFlipMap')
    string.append('cp -rf ./0/org/* ./0/')
    if len(data.modesToUse)>0 and not data.case2D: string.append('initFlx initFlexDict')
    string.append('initWaveField')
    rc = '\n'
    process = 'set -x'+rc+'('+rc 
    for s in string: process += s+rc
    process += ') 2>&1 | tee log.init'
    
    os.chdir(data.caseDir)
    subprocess.call(process, shell=True)
    if checkError('log.init'):
        print "ERROR : Something wrong happend in case initialization, check 'log.init' file"
        raise SystemExit('')
    os.chdir('..')
    
def decomposeCase(data):
    string=[]
    string.append('decomposePar -cellDist')
    rc = '\n'
    process = 'set -x'+rc+'('+rc 
    for s in string: process += s+rc
    process += ') 2>&1 | tee log.decompose'
    
    os.chdir(data.caseDir)
    subprocess.call(process, shell=True)
    os.chdir('..')

def tarCase(data):
    print 'Creating archive {}.tar.gz'.format(data.caseDir)
    process = 'tar czvf {}.tar.gz {}'.format(data.caseDir,data.caseDir)
    subprocess.call(process, shell=True)

#*** These are templates files *************************************************

fsHeader =r'''
/*--------------------------------*- C++ -*----------------------------------*\ 
| =========                 |                                                 | 
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | 
|  \\\\    /   O peration     | Version:  2.4.x                                 | 
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      | 
|    \\\\/     M anipulation  |                                                 | 
\*---------------------------------------------------------------------------*/ 
'''

run_contents='''
#!/bin/bash -l
#SBATCH -J CASE_NAME

# 5 hour wall-clock
#SBATCH --account I1608251
#SBATCH -t 3-00:00:00
#SBATCH -n NB_PROCS
#SBATCH -o log.run-%j

module load gcc/4.9.3 openmpi/1.8.4-gcc lapack/3.6.1/gcc/4.9.3
export FOAM_INST_DIR=/data/I1608251/OpenFOAM;
source /data/I1608251/OpenFOAM/OpenFOAM-2.4.x/etc/bashrc;
export LC_ALL=C

mpirun SOLVER_NAME -parallel

EOF
'''

defaultParams_contents='''
# this [template] header must exist to identify parameters for template.py
[template]

#Set name of case folder to be created
caseDir = newCase

#Set location of mesh folder, stl file name and time folder from which it should be retrieved ('0.08', '0.09' or 'latestTime')
meshDir = mesh
meshTime = latestTime

#Set name of stl file (should be located in "meshDir/constant/triSurface")
stlFile = ship.stl

#Set patch name given for hull in boundary file
hullPatch = ship

#Set OpenFOAM solver parameters 
solver = foamStar
nProcs = 24
scheme = Euler

#Set case control parameters
startTime = latestTime
endTime = 4000
timeStep = 0.05
writeInterval = 40
purgeWrite = 5
outputMotions = True
localMotionPts = 1 2 3; 4 5 6; 7 8 9
outputVBM = True
outputWave = False
waveProbes = xMin Xmax nX y zMin zMax nZ; xMin Xmax nX y zMin zMax nZ
fsiTol = 1e-6

#Set wave properties
wave = noWaves
depth = 500
draft = 0
velocity = 0.0
waveH = 1.0
waveT = 10.0
waveSeaLvl = 0
waveStartTime = 0
waveRampTime = 10
addDamping = False

#Set Relaxation zones (set 0 if no zone)
inletRelaxZone  = 400
outletRelaxZone = -250
sideRelaxZone   = 250

#Set Euler zones (set 0 if no zone)
EulerCellsDist  = 8

#Set structural data obtained with Homer
#all path should be written relatively to this input file location
datFile =  ../homer/ship.dat
donFile = ../homer/ship.don
dmigfile = ../homer/ship_dmig.pch
mdFile = ../homer/ship_md.pch
hmrUserOutput = ../homer/HmFEM.out
modesToUse = 7 8 9
shipDamping = 0.01 0.01 0.01

EOF
'''

#*** Main execution start here *************************************************
if __name__ == "__main__":
    startTime = time.time()    
    dat = cmdOptions(sys.argv)
    if dat.genCase:
        copyMesh(dat)
        if not dat.case2D: readHomerData(dat)
        foamCase_template(dat)
    if dat.initCase: initCase(dat)
    if dat.decomposeCase: decomposeCase(dat)
    if dat.tar: tarCase(dat)
    endTime = time.time()
    print 'Case generation completed in %d minutes' % ((endTime-startTime)/60)
    
