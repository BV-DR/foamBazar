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
'shipDamping' : None
}

#*** Data structure for user input data
class UserInput(object):
    def __init__(self, dir='', opts=DEFAULT_PARAMS, genCase=False, initCase=False, decomposeCase=False):
    
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
        self.outputForces  =  opts['outputForces']
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
        self.localMotionPts =  opts['localMotionPts']
        self.genCase        =  genCase
        self.initCase       =  initCase
        self.decomposeCase  =  decomposeCase
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
   
def cmdOptions(argv):
    global DEBUG
    # default global parameters
    CMD_showLog = 'echo "\nlog-file: ./log.template"; tail -45 ./log.tempalte; echo "Please see log-file: ./log.template"'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', dest='case', help='if used without input file.  Use this option to define new case folder name, default parameters will be used for other purpose, e.g.: template.py -c newCase')
    parser.add_argument('-f', dest='inputfile', help='read parameters from input file. Use this option to edit redefine parameters, e.g.: template.py -f myship.cfg ')
    parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print all available input parameters (including theirs default values) to template.cfg file and exit. Use this option to generate a template for input files for -f')
    parser.add_argument('-i','-init', dest='initCase', nargs='?', default='void', help='enables automatic case initialisation (add <directory> if not used with -f or -d option')
    parser.add_argument('-d', dest='decomposeCase', action='store_true', help='enables case decomposition for parallel run')
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
                inputdata = UserInput(dir=str(args.initCase),initCase=True,decomposeCase=args.decomposeCase)
        else:
            if args.initCase=='void':
                inputdata = UserInput(dir=str(args.case),genCase=True)
            else:
                inputdata = UserInput(dir=str(args.case),genCase=True,initCase=True,decomposeCase=args.decomposeCase)
    else:
        # here we read all input parameters from files
        fid = str(args.inputfile)
        subprocess.call('echo "Input data are read from file '+fid+'"', shell=True)
        params = readInputParams(fid)
        if args.initCase=='void':
            inputdata = UserInput(opts=params,genCase=True)
        else:
            inputdata = UserInput(opts=params,genCase=True,initCase=True,decomposeCase=args.decomposeCase)
        pass
        
    return inputdata
    pass

def getBool(string):
    if string in ['True','true','T','t','1']:
        return True
    elif string in ['False','false','F','f','0']:
        return False
    else:
        print 'Invalid boolean entry : '+str(string)
        raise SystemExit('')

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

def setValue(filename, variable, value):
    if type(value) is not str: value = str(value)
    if '"' in value: value = value.replace('"','\\"')
    if '/' in value: value = value.replace('/','\\/')
    # print 'sed -i "s/'+variable+'/'+value+'/g" '+filename
    subprocess.call('sed -i "s/'+variable+'/'+value+'/g" '+filename, shell=True)

# foamFile may be compressed i.e. ".gz"  
def foamFileExist(filename):
    found = os.path.isfile(filename)
    if not found:
        found = os.path.isfile(filename + '.gz')
    return found

def runCommand(cmd):
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        subprocess.call(CMD_showLog, shell=True)
        raise SystemExit('abort ...')
        pass
    except OSError:
        print cmd
        raise SystemExit('executable not found ... abort')
        pass
    pass
    
def run_setSet():
    # FIXME: refineMesh is buggy when run in parallel
    runCommand(CMD_setSet + ' -batch .tmp_setSet ' + CMD_keepLog)

    #if NPROCS==1:
    #    runCommand(CMD_setSet + ' -batch .tmp_setSet ' + CMD_keepLog)
    #else:
    #    runCommand('mpirun -np '+str(NPROCS)+' '+ CMD_setSet + ' -parallel -batch .tmp_setSet ' + CMD_keepLog)
    pass

def getFoamTimeFolders(dir,constant=False):
    found = []
    for val in os.listdir(dir):
        try:
            float(val)
            found.append(val)
        except ValueError:
            pass
    found.sort(key=float)
    timeFolders = [dir+'/constant/'] if constant else []
    if len(found):
        for val in found:
            timeFolders.append(val)
    return timeFolders
    pass
    
def findBoundingBox(stlFile, verbose=True):
    if verbose:
        print "Compute STL bounding box: " + stlFile
    p = subprocess.Popen("surfaceCheck "+stlFile+" | grep '^Bounding Box :' | sed \"s/.*: (//;s/[(,)]//g\" ", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    boundingBox,error = p.communicate()
    if error:
        print 'error: ', error
        raise SystemExit('abort ...')
    boundingBox = boundingBox.split(' ')
    boundingBox = [float(i) for i in boundingBox]
    if verbose:
        print "   ",boundingBox
    return boundingBox
    
def findCFDBoundingBox(case, verbose=True):
    if verbose:
        print "Compute CFD bounding box:"
    p = subprocess.Popen("fsBoundingBox -case "+case+" | grep 'Overall domain bounding box' | sed \"s/.*box (//;s/[(,)]//g\" ", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    boundingBox,error = p.communicate()
    if error:
        print 'error: ', error
        raise SystemExit('abort ...')
    boundingBox = boundingBox.split(' ')
    boundingBox = [float(i) for i in boundingBox]
    if verbose:
        print "   ", boundingBox
    return boundingBox 

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

    data.shipCOG    =  ('{:20.14e} '*3).format(*cog)
    data.shipInertia = ('{:20.14e} '*6) .format(*inertia)
    if len(modes2use)>0:
        data.hmrScaling  = ('{:20.14e}').format(scaling[0])
        data.shipFreq    = ('{:20.14e} '*len(freq)).format(*freq)

def checkError(file):
    error=False
    with open(file, 'r') as f:
        for line in f:
            if 'ERROR' in line: error=True
    
    return error
    
def foamCase_template(data):
    print 'Create system folder input files'
    #controlDict
    filename = data.caseDir+'/system/controlDict'
    if not os.path.isfile(filename):
        subprocess.call('cat << EOF > '+ filename + fsHeader + controlDict_contents, shell=True)
        setValue(filename, 'START_TIME', data.startTime)
        setValue(filename, 'END_TIME', data.endTime)
        setValue(filename, 'TIME_STEP', data.timeStep)
        setValue(filename, 'WRITE_INTERVAL', data.writeInterval)
        setValue(filename, 'PURGE_WRITE', data.purgeWrite)
        if data.outputMotions:
            setValue(filename, 'INCLUDE_MOTIONS', r'motionInfo { type motionInfo; }' )
        else:
            setValue(filename, 'INCLUDE_MOTIONS', '' )
        if data.outputVBM:
            setValue(filename, 'INCLUDE_VBM', r'#include "vbm.inc"')
            outfile = data.caseDir+'/system/vbm.inc'
            if not os.path.isfile(outfile):
                subprocess.call('cat << EOF > '+ outfile + vbm_contents, shell=True)
                setValue(outfile, 'HULL_PATCH', data.hullPatch)
                setValue(outfile, 'DON_FILENAME', data.donName)
        else:
            setValue(filename, 'INCLUDE_VBM', '')
        if data.outputWave:
            setValue(filename, 'INCLUDE_WAVE', r'#include "waveProbe.inc"')
            outfile = data.caseDir+'/system/waveProbe.inc'
            if not os.path.isfile(outfile):
                subprocess.call('cat << EOF > '+ outfile + waveProbe_contents, shell=True)
        else:
            setValue(filename, 'INCLUDE_WAVE', '')
        if data.outputForces:
            setValue(filename, 'INCLUDE_FORCES', r'#include "forces.inc"')
            outfile = data.caseDir+'/system/forces.inc'
            if not os.path.isfile(outfile):
                subprocess.call('cat << EOF > '+ outfile + forces_contents, shell=True)
                setValue(outfile, 'HULL_PATCH', data.hullPatch)
        else:
            setValue(filename, 'INCLUDE_FORCES', '')
        
        if len(data.localMotionPts)>0:
            setValue(filename, 'INCLUDE_LOCALMOTION', r'#include "localMotion.inc"')
            outfile = data.caseDir+'/system/localMotion.inc'
            if not os.path.isfile(outfile):
                subprocess.call('cat << EOF > '+ outfile + localMotion_contents, shell=True)
        else:
            setValue(filename, 'INCLUDE_LOCALMOTION', '')

    #fvSchemes
    filename = data.caseDir+'/system/fvSchemes'
    if not os.path.isfile('system/fvSchemes'):
        subprocess.call('cat << EOF > '+ filename + fvSchemes_contents, shell=True)
        setValue(filename, 'DDT_SCHEME', data.scheme)

    #fvSolution
    filename = data.caseDir+'/system/fvSolution'
    if not os.path.isfile('system/fvSolution'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + fvSolution_contents, shell=True)
        setValue(filename, 'FSI_TOL', data.fsiTol)
    
    #decomposeParDict
    filename = data.caseDir+'/system/decomposeParDict'
    if not os.path.isfile('system/decomposeParDict'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + decomposeParDict_contents, shell=True)
        setValue(filename, 'NB_PROCS', data.nProcs)

    print 'Create constant folder input files'
    #waveProperties
    filename = data.caseDir+'/constant/waveProperties'
    if not os.path.isfile('constant/waveProperties'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + waveProperties_contents, shell=True)
        setValue(filename, 'WAVE_TYPE', data.wave)
        setValue(filename, 'WAVE_START_TIME', data.waveStartTime)
        setValue(filename, 'WAVE_RAMP_TIME', data.waveRampTime)
        setValue(filename, 'WAVE_VELOCITY', data.velocity)
        setValue(filename, 'WAVE_DEPTH', data.depth)
        setValue(filename, 'WAVE_HEIGHT', data.waveH)
        setValue(filename, 'WAVE_PERIOD', data.waveT)
        setValue(filename, 'WAVE_SEALEVEL', data.waveSeaLvl)
        
        bBox = findCFDBoundingBox(data.caseDir,False)
        setValue(filename, 'BOUND_X0', bBox[0])
        setValue(filename, 'BOUND_X1', bBox[3])
        setValue(filename, 'BOUND_Y1', bBox[4])
        relaxDomain = ''
        if data.sideRelaxZone   is not None: relaxDomain += 'domainY1 '
        if data.outletRelaxZone is not None: relaxDomain += 'domainX0 '
        if data.inletRelaxZone  is not None: relaxDomain += 'domainX1 '
        setValue(filename, 'RELAX_DOMAIN', relaxDomain)
    
    #turbulenceProperties
    filename = data.caseDir+'/constant/turbulenceProperties'
    if not os.path.isfile('constant/turbulenceProperties'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + turbulenceProperties_contents, shell=True)
        
    #transportProperties
    filename = data.caseDir+'/constant/transportProperties'
    if not os.path.isfile('constant/transportProperties'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + transportProperties_contents, shell=True)
    
    #dynamicMeshDict
    filename = data.caseDir+'/constant/dynamicMeshDict'
    if not os.path.isfile('constant/dynamicMeshDict'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + dynamicMeshDict_contents, shell=True)
        setValue(filename, 'HULL_PATCH', data.hullPatch)
        if data.addDamping:
            bBox = findBoundingBox(data.caseDir+'/constant/triSurface/'+data.stlFile,False)
            LBP = round(0.95*abs(bBox[3]-bBox[0]),2)
            BC = round(abs(bBox[4]-bBox[1]),2)
            rc = '\\n'
            dmp = 'dampingCoeff'                        +rc+ \
                  '        {'                           +rc+ \
                  '        type dampingSinkageAndTrim;' +rc+ \
                  '        heaveDampingCoef 1000;'      +rc+ \
                  '        pitchDampingCoef 1000;'      +rc+ \
                  '        Lpp '+str(LBP)+';'           +rc+ \
                  '        Breadth '+str(BC)+';'        +rc+ \
                  '        }'
            setValue(filename, 'SW_DAMPING', dmp)
        else:
            setValue(filename, 'SW_DAMPING', '')

    #g
    filename = data.caseDir+'/constant/g'
    if not os.path.isfile('constant/g'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + g_contents, shell=True)
        
    #RASProperties
    filename = data.caseDir+'/constant/RASProperties'
    if not os.path.isfile('constant/RASProperties'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + RASProperties_contents, shell=True)
    
    #*.flex
    if len(data.modesToUse)>0:
        filename = data.caseDir+'/constant/'+data.donName+'.flex'
        if not os.path.isfile('constant/'+data.donName+'.flex'):
            subprocess.call('cat << EOF > '+ filename + flex_contents, shell=True)
            setValue(filename, 'SHIP_FREQ', data.shipFreq)
            setValue(filename, 'SHIP_DAMPING', data.shipDamping)
            dmig_name = ((data.dmigFile).rpartition('/')[-1]).partition('_dmig')[0]
            setValue(filename, 'DMIG_NAME', dmig_name)
            setValue(filename, 'HULL_PATCH', data.hullPatch) 
            if len(data.localMotionPts)>0:
                setValue(filename, 'LOCAL_MOTION', '#include "modeshapes/PTS_localMotion.flx"')
            else:
                setValue(filename, 'LOCAL_MOTION', '')

    print 'Create 0 folder input files'
    #alpha water
    filename = data.caseDir+'/0/org/alpha.water'
    if not os.path.isfile('0/org/alpha.water'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + alphawater_contents, shell=True)
        setValue(filename, 'HULL_PATCH', data.hullPatch)

    #p_rgh
    filename = data.caseDir+'/0/org/p_rgh'
    if not os.path.isfile('0/org/p_rgh'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + p_rgh_contents, shell=True)
        setValue(filename, 'HULL_PATCH', data.hullPatch)

    #U
    filename = data.caseDir+'/0/org/U'
    if not os.path.isfile('0/org/U'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + U_contents, shell=True)
        setValue(filename, 'HULL_PATCH', data.hullPatch)
        setValue(filename, 'WAVE_VELOCITY', data.velocity)
        
    #pointDisplacement
    filename = data.caseDir+'/0/org/pointDisplacement'
    if not os.path.isfile('0/org/pointDisplacement'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + pointDisplacement_contents, shell=True)
        setValue(filename, 'HULL_PATCH', data.hullPatch)
    
    #sixDofDomainBody
    filename = data.caseDir+'/0/uniform/sixDofDomainBody'
    if not os.path.isfile('0/uniform/sixDofDomainBody'):
        subprocess.call('cat << EOF > '+ filename + fsHeader + sixDofDomainBody_contents, shell=True)
        setValue(filename, 'SHIP_MASS', data.shipMass)
        setValue(filename, 'SHIP_INERTIA', data.shipInertia)
        setValue(filename, 'SHIP_COG', data.shipCOG)
        if len(data.modesToUse)>0:
            rc = '\\n'
            modal = 'modalData'                                        +rc+ \
                    '{'                                                +rc+ \
                    'readFromFile   "constant/'+data.donName+'.flex";' +rc+ \
                    '}'
            setValue(filename, 'MODAL_DATA', modal)
        else:
            setValue(filename, 'MODAL_DATA', '')
     
    print 'Create other input files'
    #initFlexDict
    if len(data.modesToUse)>0:
        filename = data.caseDir+'/initFlexDict'
        if not os.path.isfile('initFlexDict'):
            subprocess.call('cat << EOF > '+ filename + fsHeader + initFlexDict_contents, shell=True)
            setValue(filename, 'MD_FILE', data.mdFile)
            setValue(filename, 'MODES_TO_USE', data.modesToUse)
            setValue(filename, 'DAT_FILE', data.datFile)
            setValue(filename, 'DMIG_FILE', data.dmigFile)
            setValue(filename, 'HMR_SCALE', data.hmrScaling)
            setValue(filename, 'DRAFT', data.draft)
            setValue(filename, 'HULL_PATCH', data.hullPatch)
            if len(data.localMotionPts)>0:
                rc = '\\n'
                lms = '    pointList (localMotion);' + rc + \
                      '    localMotion' + rc + \
                      '    (' + rc
                for pt in data.localMotionPts: lms += '        ('+('{:8.6e} '*3).format(*pt.tolist())+')'+ rc
                lms +='     );' + rc
                setValue(filename, 'LOCAL_MOTION', lms)
            else:
                setValue(filename, 'LOCAL_MOTION', '') 
    
    #cell.Set
    filename = data.caseDir+'/cell.Set'
    if not os.path.isfile('cell.Set'):
        subprocess.call('cat << EOF > '+ filename + cellSet_contents, shell=True)
        setValue(filename, 'STL_FILE', data.stlFile)
        bBox = findBoundingBox(data.caseDir+'/constant/triSurface/'+data.stlFile,False)
        bBox = [0.5*(bBox[0]+bBox[3]), 0.5*(bBox[1]+bBox[4])+math.fabs(bBox[1]-bBox[4]), 0.5*(bBox[2]+bBox[5])]
        outsidePoints = str(bBox[0])+" "+str(bBox[1])+" "+str(bBox[2])
        setValue(filename, 'OUTSIDE_POINT', outsidePoints)
        setValue(filename, 'EULER_DIST', data.EulerCellsDist)
        if data.inletRelaxZone is not None:
            setValue(filename, 'INLET_X', data.inletRelaxZone)
        else:
            setValue(filename, 'cellSet inletZone', r'//cellSet inletZone')
        if data.outletRelaxZone is not None:
            setValue(filename, 'OUTLET_X', data.outletRelaxZone)
        else:
            setValue(filename, 'cellSet outletZone', r'//cellSet outletZone')
        if data.sideRelaxZone is not None:
            setValue(filename, 'SIDE_X', data.sideRelaxZone)
        else:
            setValue(filename, 'cellSet sideZone', r'//cellSet sideZone')
            
    #run.sh
    filename = data.caseDir+'/run.sh'
    if not os.path.isfile('run.sh'):
        subprocess.call('cat << EOF > '+ filename + run_contents, shell=True)
        setValue(filename, 'CASE_NAME', data.caseDir)
        setValue(filename, 'NB_PROCS', data.nProcs)
        setValue(filename, 'SOLVER_NAME', data.solver)
        
def copyMesh(data):

    subprocess.call('mkdir -p '+data.caseDir+'/{0/org,0/uniform,constant,system}', shell=True)

    if data.meshTime=='latestTime':
        timeFolders = getFoamTimeFolders(data.meshDir)
        meshTimeFolder = timeFolders[-1]
    else:
        meshTimeFolder = data.meshTime
        
    timeFolders = getFoamTimeFolders(data.meshDir)
    
    print 'Copy mesh from folder '+meshTimeFolder
    subprocess.call('cp -r '+data.meshDir+'/'+meshTimeFolder+'/polyMesh '+data.caseDir+'/constant', shell=True)

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
    if len(data.modesToUse)>0: string.append('initFlx initFlexDict')
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

controlDict_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by template on $(hostname) @ $(date)
*/

application     foamStar;

startFrom       START_TIME;

startTime       0;

stopAt          endTime;

endTime         END_TIME;

deltaT          TIME_STEP;

writeControl    timeStep;

writeInterval   WRITE_INTERVAL;

purgeWrite      PURGE_WRITE;

writeFormat     ascii;

writePrecision  15;

writeCompression compressed;

timeFormat      general;

timePrecision   6;

runTimeModifiable no;

adjustTimeStep  no;

maxCo           1;
maxAlphaCo      1;
maxDeltaT       1;

libs
(
    "libfoamStar.so"
);
  
// ************************************************************************* //

functions
{
    INCLUDE_MOTIONS
    INCLUDE_VBM
    INCLUDE_WAVE
    INCLUDE_FORCES
    INCLUDE_LOCALMOTION
}

EOF
'''

fvSchemes_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSchemes;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by template on $(hostname) @ $(date)
*/

ddtSchemes
{
    default      DDT_SCHEME;
    ddt(U)       Euler;
}

gradSchemes
{
    default     Gauss linear;
}

divSchemes
{
    div(rhoPhi,U)  Gauss linearUpwindV GradU;
    div(phi,alpha) Gauss vanLeer;
    div(phirb,alpha) Gauss interfaceCompression;
    div((muEff*dev(T(grad(U))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default         linear;
}

snGradSchemes
{
    default         corrected;
}

fluxRequired
{
    default         no;
    p_rgh;
    pcorr;
    alpha.water;

}

// ************************************************************************* //

EOF
'''

fvSolution_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      fvSolution;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by template on $(hostname) @ $(date)
*/

solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      3;
        nAlphaSubCycles 1;
        cAlpha          0.3;

        MULESCorr       yes;
        nLimiterIter    5;
        alphaApplyPrevCorr  no;

        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-8;
        relTol          0;
        minIter         2;
    }

    "pcorr.*"
    {
        solver          PCG;
        preconditioner
        {
            preconditioner  GAMG;
            tolerance       1e-7;
            relTol          0;
            smoother        DICGaussSeidel;
            nPreSweeps      0;
            nPostSweeps     2;
            nFinestSweeps   2;
            cacheAgglomeration false;
            nCellsInCoarsestLevel 10;
            agglomerator    faceAreaPair;
            mergeLevels     1;
        }

        tolerance       1e-08;
        relTol          0;
        maxIter         1000;
        minIter         2;

    };

    p_rgh
    {
        solver          GAMG;
        tolerance       1e-8;
        relTol          0;
        smoother        DIC;
        nPreSweeps      0;
        nPostSweeps     2;
        nFinestSweeps   2;
        cacheAgglomeration true;
        nCellsInCoarsestLevel 10;
        agglomerator    faceAreaPair;
        mergeLevels     1;
        minIter         2;

    };

    p_rghFinal
    {
        solver          PCG;
        preconditioner
        {
            preconditioner  GAMG;
            tolerance       1e-7;
            relTol          0;
            nVcycles        2;
            smoother        DICGaussSeidel;
            nPreSweeps      2;
            nPostSweeps     2;
            nFinestSweeps   2;
            cacheAgglomeration true;
            nCellsInCoarsestLevel 10;
            agglomerator    faceAreaPair;
            mergeLevels     1;
        }

        tolerance       1e-8;
        relTol          0;
        maxIter         1000;
        minIter         2;

    };
    
    "(U|k|epsilon)"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0;
        nSweeps         2;
        minIter         2;
    }

    "(U|k|epsilon)Final"
    {
        solver          smoothSolver;
        smoother        GaussSeidel;
        tolerance       1e-8;
        relTol          0;
        nSweeps         2;
        minIter         2;
    }

    "cellDisplacement.*"
    {
        solver          GAMG;
        tolerance       1e-8;
        relTol          0;
        smoother        GaussSeidel;
        cacheAgglomeration true;
        nCellsInCoarsestLevel 10;
        agglomerator    faceAreaPair;
        mergeLevels     1;

    }
}

PIMPLE
{
    EulerCells  EulerCells;
    momentumPredictor   yes;
    nOuterCorrectors    22;
    fsiTol              FSI_TOL;
    fsiMaxIter          23;
    nCorrectors         4;
    nNonOrthogonalCorrectors 0;
    correctPhi          no;
    moveMeshOuterCorrectors yes;
}

relaxationFactors
{
    equations
    {
        U       1;
        UFinal  1;
        p_rgh   1;
        p_rghFinal 1;
    }
}

// ************************************************************************* //

EOF
'''

decomposeParDict_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by template on $(hostname) @ $(date)
*/

numberOfSubdomains NB_PROCS;
method          scotch;

// ************************************************************************* //

EOF
'''

vbm_contents='''
vbm
{
    type internalLoads; outputControl timeStep; outputInterval 1;
    hullPatches         (HULL_PATCH);
    donFileName         "DON_FILENAME.don";
    wldFileName         "";
    log true;
    ySymmetry           yes;
}

// ************************************************************************* //

EOF
'''

forces_contents='''
forces
    {
        type forces;
        functionObjectLibs ( "libforces.so" );
        patches (HULL_PATCH);
        rhoInf  1000;
        rhoName rho;
        pName   p;
        UName   U;
        log     on;
        outputControl   timeStep;
        outputInterval  1;
        CofR (0 0 0);
    }

// ************************************************************************* //

EOF
'''

waveProbe_contents='''
waveProbe
{
    type surfaceElevation;
    fields              (alpha.water);     // default is alpha.water
    writePrecision      6;                 // default is 6
    interpolationScheme cellPointFace;     // default is cellPointFace
    outputControl timeStep; outputInterval 1;

    // using syntax from libsampling.so
    sets
    (
        wp1
        {
            start (315 1 -5.6);
            end (315 1 3.5);
            type face; axis z;
            nPoints 100;
        }

        wp2
        {
            start (340 1 -5.6);
            end (340 1 3.5);
            type face; axis z;
            nPoints 100;
        }
    );
}

// ************************************************************************* //

EOF
'''

localMotion_contents='''
localMotion
{
    type localMotion;
    motionData sixDofDomainBody;
    flxPTS localMotion;
}

// ************************************************************************* //

EOF
'''

waveProperties_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      environmentalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#inputMode overwrite;

mycase
{
    waveType    WAVE_TYPE; //streamFunction
    refTime     0;
    startTime   WAVE_START_TIME;
    rampTime    WAVE_RAMP_TIME;
    refDirection (-1 0 0);
    U0 (-WAVE_VELOCITY 0 0);
    depth WAVE_DEPTH;
    height WAVE_HEIGHT;
    period WAVE_PERIOD;
}

seaLevel    WAVE_SEALEVEL;

initCoeffs
{
    \$mycase
}

relaxationNames (RELAX_DOMAIN);

domainX1Coeffs
{
    \$mycase
    relaxationZone
    {
        zoneName         inletZone;
        origin           (BOUND_X1 0 0);
        orientation      (-1 0 0);
    }
}

domainX0Coeffs
{
    \$mycase
    relaxationZone
    {
        zoneName         outletZone;
        origin           (BOUND_X0 0 0);
        orientation      (1 0 0);
    }
}

domainY1Coeffs
{
    \$mycase
    relaxationZone
    {
        zoneName         sideZone;
        origin           (0 BOUND_Y1 0);
        orientation      (0 -1 0);
    }
}

// ************************************************************************* //

EOF
'''

turbulenceProperties_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType laminar; //RASModel;

// ************************************************************************* //

EOF
'''

transportProperties_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      transportProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

phases (water air);

"water|phase1"
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1.14e-06;
    rho             rho [ 1 -3 0 0 0 0 0 ] 1025;
}

"air|phase2"
{
    transportModel  Newtonian;
    nu              nu [ 0 2 -1 0 0 0 0 ] 1.48e-05;
    rho             rho [ 1 -3 0 0 0 0 0 ] 1;
}

sigma           sigma [ 1 0 -2 0 0 0 0 ] 0.0;

// ************************************************************************* //

EOF
'''

g_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       uniformDimensionedVectorField;
    location    "constant";
    object      g;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -2 0 0 0 0];
value           ( 0 0 -9.81 );

// ************************************************************************* //

EOF
'''

dynamicMeshDict_contents='''
FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    object          dynamicMeshDict;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//dynamicFvMesh   staticFvMesh;
dynamicFvMesh   sixDofDomainFvMesh;

motionSolverLibs ("libfvMotionSolvers.so");
solver              displacementLaplacian;
displacementLaplacianCoeffs
{
    diffusivity  invDistSqr (HULL_PATCH);
    //diffusivity  invDistSqr (InternalFaces);
}

sixDofDomainFvMeshCoeffs
{
    hullPatches     (HULL_PATCH);
    meshConfigType  4;

    solver          RKCK45;
    absTol          1e-8;
    relTol          0;

    rampTime        15;

    loads
    {
        gravity { type gravity; value (0 0 -9.81); }
        fluid { type fluidForce; patches (HULL_PATCH); dynRelax off; relaxCoeff 0.5; ySym yes;}
        SW_DAMPING
    }

    constraints
    {
        heavePitch { type cog6DOF; default "fixAll"; except (heave pitch); }
    }
}

// ************************************************************************* //

EOF
'''

RASProperties_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      RASProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

RASModel            laminar; //rhoKOmegaSST;

turbulence          off;

printCoeffs         on;

// ************************************************************************* //

EOF
'''

alphawater_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.water;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    "domainX1"
    {
        type            waveAlpha;
        value           uniform 0;
    }

    "domainX0"
    {
        type            waveAlpha;
        value           uniform 0;
    }

    "domainZ0"
    {
        type            zeroGradient;
    }

    "domainY0"
    {
        type            symmetryPlane;
    }

    "domainY1"
    {
        type            symmetryPlane;
    }

    "domainZ1"
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }

    "HULL_PATCH"
    {
        type            zeroGradient;
    }
}

// ************************************************************************* //

EOF
'''

p_rgh_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      p_rgh;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [1 -1 -2 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    "domainX1"
    {
        type            fixedFluxPressure;
    }

    "domainX0"
    {
        type            fixedFluxPressure;
    }

    "domainZ0"
    {
        type            fixedFluxPressure;
    }

    "domainY0"
    {
        type            symmetryPlane;
    }

    "domainY1"
    {
        type            symmetryPlane;
    }

    "domainZ1"
    {
        type            totalPressure;
        p0              uniform 0;
        U               U;
        phi             phi;
        rho             rho;
        psi             none;
        gamma           1;
        value           uniform 0;
    }

    "HULL_PATCH"
    {
        type            fixedFluxPressure;
    }
}

// ************************************************************************* //

EOF
'''

U_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       volVectorField;
    location    "0";
    object      U;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 -1 0 0 0 0];

minusFwdSpeed   (-WAVE_VELOCITY 0 0);

internalField   uniform \$minusFwdSpeed;

boundaryField
{
    "domainX1"
    {
        type            waveVelocity;
        value           uniform (0 0 0);
    }

    "domainX0"
    {
        type            waveVelocity;
        value           uniform (0 0 0);
    }

    "domainZ0"
    {
        type            fixedValue;
        value           uniform \$minusFwdSpeed;
    }

    "domainY0"
    {
        type            symmetryPlane;
    }

    "domainY1"
    {
        type            symmetryPlane;
    }

    "domainZ1"
    {
        type            pressureInletOutletVelocity;
        tangentialVelocity uniform \$minusFwdSpeed;
        value           uniform (0 0 0);
    }

    "HULL_PATCH"
    {
        type            movingWallVelocity;
        value           uniform (0 0 0);
    }
}

// ************************************************************************* //

EOF
'''

pointDisplacement_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       pointVectorField;
    location    "0.01";
    object      pointDisplacement;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 0 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{
    "domainX1"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "domainX0"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "domainZ0"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "domainY1"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "domainY0"
    {
        type            symmetryPlane;
    }

    "domainZ1"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "HULL_PATCH"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }

    "InternalFaces"
    {
        type            fixedValue;
        value           uniform (0 0 0);
    }
}

// ************************************************************************* //

EOF
'''

cellSet_contents='''
cellSet EulerCells new surfaceToCell "./constant/triSurface/STL_FILE" ((OUTSIDE_POINT)) yes yes no EULER_DIST -1e6
cellSet inletZone new boxToCell (INLET_X -1e6 -1e6) (1e6 1e6 1e6)
cellSet outletZone new boxToCell (-1e6 -1e6 -1e6) (OUTLET_X 1e6 1e6)
cellSet sideZone new boxToCell (-1e6 SIDE_X -1e6) (1e6 1e6 1e6)

EOF
'''

initFlexDict_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      initFlxDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// initFlx will read each block and create interpolated mode shapes for CFD mesh
//- src data is in vtuFile
//- interpolated data from patches is in "constant/modeshapes/CFD_<name>*"
//- interpolated data from point{Set,List} is in "constant/modeshapes/PTS_<name>*"

//- note: pointSet are fvMesh points, but pointList are user-defined pointField

FEM_STRUCTURALMESH_VTU
{
    mdFile "../MD_FILE"; selected (MODES_TO_USE);
    datFile "../DAT_FILE";
    dmigMfile "../DMIG_FILE";
    dmigKfile "../DMIG_FILE";
    pchCoordinate (0 0 -DRAFT 0 0 0);
    pchScaleMode  HMR_SCALE;
    pchLengthUnit 1;
    pchMassUnit   1;

    patches (HULL_PATCH); ySym (true);

LOCAL_MOTION

}
// ************************************************************************* //

EOF
'''

sixDofDomainBody_contents='''
FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    object          singleBody;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

mass SHIP_MASS;

momentOfInertia (SHIP_INERTIA);

cogInitial (SHIP_COG);

Xrel (0 0 0);

dotXrel (0 0 0);

omega (0 0 0);

EulerZYX
{
    rollPitchYaw (0 0 0);
}

MODAL_DATA

// ************************************************************************* //

EOF
'''

flex_contents='''
wn2 (SHIP_FREQ);
dampingRatio (SHIP_DAMPING);


#include "modeshapes/DMIG_NAME_dmig.prj"

// These folder should be named after hull patch name : CFD_{patch name}_fpt.flx
#include "modeshapes/CFD_HULL_PATCH_fpt.flx"
#include "modeshapes/CFD_HULL_PATCH_fps.flx"
#include "modeshapes/CFD_HULL_PATCH_mpt.flx"
LOCAL_MOTION

EOF
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

#Set location of mesh folder, stl file name and time folder from which it should be retrived ('0.08', '0.09' or 'latestTime')
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
outputVBM = True
outputWave = False
localMotionPts = 1 2 3; 4 5 6; 7 8 9
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

#Set Euler and Relaxation zones (set 0 if no zone)
EulerCellsDist  = 8
inletRelaxZone  = 400
outletRelaxZone = -250
sideRelaxZone   = 250

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
        readHomerData(dat)
        foamCase_template(dat)
    if dat.initCase: initCase(dat)
    if dat.decomposeCase: decomposeCase(dat)
    endTime = time.time()
    print 'Case generation completed in %d minutes' % ((endTime-startTime)/60)
    
