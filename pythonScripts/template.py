#!/usr/bin/env upython

#########################################################################
# Filename: template.py                                                 #
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

# abenhamou: 2017-may-02

DEBUG = False
DEFAULT_CFG_FILE = 'template.cfg'
DEFAULT_PARAMS = {
'caseDir' : 'newCase',
'meshDir' : '',
'meshTime' : 'latestTime',
'nProcs' : 4,
'scheme' : 'Euler',
'meshid' : None,
'startTime' : 'latestTime',
'endTime' : 1000,
'timeStep' : 0.01,
'writeInterval' : 1,
'purgeWrite' : 0,
'outputVBM' : False,
'outputWave' : False,
'outputMotions' : True,
'fsiTol' : 1e-6,
'hullPatch' : 'ship',
'donFile' : None,
'wave' : None,
'waveH' : 1.0,
'waveT' : 10.0,
'velocity' : 0., #in m/s
'depth' : 500,
'draft' : None,
'waveSeaLvl' : 0,
'waveStartTime' : 0,
'waveRampTime' : 0,
}

#*** Data structure for user input data
class UserInput(object):
    def __init__(self, opts=DEFAULT_PARAMS):
    
        self.caseDir        =  opts['caseDir']
        self.meshDir        =  opts['meshDir']
        self.meshTime       =  opts['meshTime']
        self.nProcs         =  opts['nProcs']
        self.scheme         =  opts['scheme']
        self.meshid         =  opts['meshid']
        self.startTime      =  opts['startTime']
        self.endTime        =  opts['endTime']
        self.timeStep       =  opts['timeStep']
        self.writeInterval  =  opts['writeInterval']
        self.purgeWrite     =  opts['purgeWrite']
        self.outputVBM      =  opts['outputVBM']
        self.outputWave     =  opts['outputWave']
        self.outputMotions  =  opts['outputMotions']
        self.fsiTol         =  opts['fsiTol']
        self.hullPatch      =  opts['hullPatch']
        self.donFile        =  opts['donFile']
        self.wave           =  opts['wave']
        self.waveH          =  opts['waveH']
        self.waveT          =  opts['waveT']
        self.velocity       =  opts['velocity']
        self.depth          =  opts['depth']
        self.draft          =  opts['draft']
        self.waveSeaLvl     =  opts['waveSeaLvl']
        self.waveStartTime  =  opts['waveStartTime']
        self.waveRampTime   =  opts['waveRampTime']
   
def cmdOptions(argv):
    global DEBUG
    # default global parameters
    CMD_showLog = 'echo "\nlog-file: ./log.template"; tail -45 ./log.tempalte; echo "Please see log-file: ./log.template"'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='inputfile', help='read parameters from input file. Use this option with a stl-file to generate a mesh with default parameters, e.g.: template.py -f myship.stl ')
    parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print all available input parameters (including theirs default values) to template.cfg file and exit. Use this option to generate a template for input files for -f')
    args = parser.parse_args()
    
    if args.showConfig:
        print "\nOutput default parameters to file:", DEFAULT_CFG_FILE
        subprocess.call('cat << EOF > '+DEFAULT_CFG_FILE + defaultParams_contents, shell=True)
        raise SystemExit('')

    if args.inputfile==None:
        subprocess.call('echo "No input file provided, default value taken into acount"', shell=True)
        inputdata = UserInput()
    else:
        # here we read all input parameters from files
        fid = str(args.inputfile)
        subprocess.call('echo "Input data are read from file '+fid+'"', shell=True)
        params = readInputParams(fid)
        inputdata = UserInput(opts=params)
        pass
        
    return inputdata
    pass

def readInputParams(filename):
    params = dict(DEFAULT_PARAMS)
    config = configparser.ConfigParser()
    config.read(filename)
    name = 'template'
    
    print "Read input paramters from file:", filename
    
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

    # 'scheme' : 'Euler', 'CrankNicolson 0.9'
    try:
        txt = str(config[name]['scheme'])
        params['scheme'] = txt
    except KeyError: pass

    # mesh id
    try:
        txt = str(config[name]['meshid'])
        params['meshid'] = float(txt)
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
        params['outputMotions'] = txt
    except KeyError: pass
    
    # 'outputVBM' : output VBM or not?
    try:
        txt = str(config[name]['outputVBM'])
        params['outputVBM'] = txt
    except KeyError: pass
    
    # 'outputWave' : output Waves or not?
    try:
        txt = str(config[name]['outputWave'])
        params['outputWave'] = txt
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

    return params
    pass

def setValue(filename, variable, value):
    subprocess.call('sed -i "s/'+variable+'/'+str(value)+'/g" '+filename, shell=True)

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
    
def findBoundingBox(stlFile, verbose=True):
    if verbose:
        print "Compute bounding box: " + stlFile
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

def caseAlreadyDecomposed():
    isDecomposed = False
    nProcs = 1
    for val in os.listdir('.'):
        if val.startswith('processor'):
            isDecomposed=True
            nProcs = max(nProcs, int(val[9:])+1)
        pass
    return isDecomposed, nProcs
    pass

def foamCase_template(data):
    global NPROCS
    subprocess.call('mkdir -p '+data.caseDir+'/{0,constant,system}', shell=True)
    
    #controlDict
    filename = data.caseDir+'/system/controlDict'
    if not os.path.isfile(filename):
        subprocess.call('cat << EOF > '+ filename + controlDict_contents, shell=True)
        setValue(filename, 'START_TIME', data.startTime)
        setValue(filename, 'END_TIME', data.endTime)
        setValue(filename, 'TIME_STEP', data.timeStep)
        setValue(filename, 'WRITE_INTERVAL', data.writeInterval)
        setValue(filename, 'PURGE_WRITE', data.purgeWrite)
        if not data.outputMotions: setValue(filename, 'motionInfo {', '//motionInfo {')
        if not data.outputVBM:
            setValue(filename, '#include "vbm.inc"', '//#include "vbm.inc"')
        else:
            outfile = data.caseDir+'/system/vbm.inc'
            if not os.path.isfile(outfile):
                subprocess.call('cat << EOF > '+ outfile + vbm_contents, shell=True)
                setValue(outfile, 'HULL_PATCH', data.hullPatch)
                setValue(outfile, 'DON_FILENAME', data.donFile)
        if not data.outputWave:
            setValue(filename, '#include "waveProbe.inc"', '//#include "waveProbe.inc"')
        else:
            outfile = data.caseDir+'/system/waveProbe.inc'
            if not os.path.isfile(outfile):
                subprocess.call('cat << EOF > '+ outfile + waveProbe_contents, shell=True)

    #fvSchemes
    filename = data.caseDir+'/system/fvSchemes'
    if not os.path.isfile('system/fvSchemes'):
        subprocess.call('cat << EOF > '+ filename + fvSchemes_contents, shell=True)
        setValue(filename, 'DDT_SCHEME', data.scheme)

    #fvSolution
    filename = data.caseDir+'/system/fvSolution'
    if not os.path.isfile('system/fvSolution'):
        subprocess.call('cat << EOF > '+ filename + fvSolution_contents, shell=True)
        setValue(filename, 'FSI_TOL', data.fsiTol)
    
    #decomposeParDict
    filename = data.caseDir+'/system/decomposeParDict'
    if not os.path.isfile('system/decomposeParDict'):
        subprocess.call('cat << EOF > '+ filename + decomposeParDict_contents, shell=True)
        setValue(filename, 'NB_PROCS', data.nProcs)

    #waveProperties
    filename = data.caseDir+'/constant/waveProperties'
    if not os.path.isfile('constant/waveProperties'):
        subprocess.call('cat << EOF > '+ filename + waveProperties_contents, shell=True)
        setValue(filename, 'WAVE_TYPE', data.wave)
        setValue(filename, 'WAVE_START_TIME', data.waveStartTime)
        setValue(filename, 'WAVE_RAMP_TIME', data.waveRampTime)
        setValue(filename, 'WAVE_VELOCITY', data.velocity)
        setValue(filename, 'WAVE_DEPTH', data.depth)
        setValue(filename, 'WAVE_HEIGHT', data.waveH)
        setValue(filename, 'WAVE_PERIOD', data.waveT)
        setValue(filename, 'WAVE_SEALEVEL', data.waveSeaLvl)
        #relaxation zones must be handeled manually
    
    #turbulenceProperties
    filename = data.caseDir+'/constant/turbulenceProperties'
    if not os.path.isfile('constant/turbulenceProperties'):
        subprocess.call('cat << EOF > '+ filename + turbulenceProperties_contents, shell=True)
    
    #dynamicMeshDict
    filename = data.caseDir+'/constant/dynamicMeshDict'
    if not os.path.isfile('constant/dynamicMeshDict'):
        subprocess.call('cat << EOF > '+ filename + dynamicMeshDict_contents, shell=True)
        
    #g
    filename = data.caseDir+'/constant/g'
    if not os.path.isfile('constant/g'):
        subprocess.call('cat << EOF > '+ filename + g_contents, shell=True)
        
    #RASProperties
    filename = data.caseDir+'/constant/RASProperties'
    if not os.path.isfile('constant/RASProperties'):
        subprocess.call('cat << EOF > '+ filename + RASProperties_contents, shell=True)
        
def copyMesh(data):

    if data.meshTime=='latestTime':
        timeFolders = getFoamTimeFolders(data.meshDir)
        print timeFolders
        meshTimeFolder = timeFolders[-1]
    else:
        meshTimeFolder = data.meshTime
        
    timeFolders = getFoamTimeFolders(data.meshDir)
    
    print 'Copy mesh from folder '+meshTimeFolder
    subprocess.call('cp -r '+data.meshDir+'/'+meshTimeFolder+'/polyMesh '+data.caseDir+'/constant/.', shell=True)


#*** These are templates files ************************************************* 

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

maxCo           0.75;
maxAlphaCo      0.5;
maxDeltaT       0.1;

libs
(
    "libfoamStar.so"
);
  
// ************************************************************************* //

functions
{
    motionInfo { type motionInfo; }
    #include "forces.inc"
    #include "vbm.inc"
    #include "waveProbe.inc"
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
    default     Gauss Linear;
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
    nCorrectors         4;
    nNonOrthogonalCorrectors 1;
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
}

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
    waveType    WAVE_TYPE //streamFunction
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
    $mycase
}

relaxationNames (domainX0);

domainX1Coeffs
{
    $mycase
    //relaxationZone
    //{
    //    zoneName         inletZone;
    //    origin           (737.328 0 0);
    //    orientation      (-1 0 0);
    //}
}

domainX0Coeffs
{
    $mycase
    relaxationZone
    {
        zoneName         outletZone;
        origin           (-867.091 0 0);
        orientation      (1 0 0);
    }
}

//domainY1Coeffs
//{
//    $mycase
//    relaxationZone
//    {
//        zoneName         sideZone;
//        origin           (0 566.266 0);
//        orientation      (0 -1 0);
//    }
//}

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
    diffusivity  invDistSqr (ship);
    //diffusivity  invDistSqr (InternalFaces);
}

sixDofDomainFvMeshCoeffs
{
    hullPatches     (ship);
    meshConfigType  4;

    solver          RKCK45;
    absTol          1e-8;
    relTol          0;

    rampTime        15;

    loads
    {
        gravity { type gravity; value (0 0 -9.81); }
        fluid { type fluidForce; patches (ship); dynRelax off; relaxCoeff 0.5; }
        //dampingCoeff
        //{
        //      type dampingSinkageAndTrim;
        //        heaveDampingCoef 1000;
        //        pitchDampingCoef 1000;
        //        Lpp 281;
        //        Breadth 26.84;
        // }


    }

    constraints
    {
        heavePitch { type cog6DOF; default "fixAll"; except (heave pitch); }
    }
}

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

EOF
'''

defaultParams_contents='''
# this [template] header must exist to identify parameters for template.py
[template]

# caseDir=
meshDir='mesh'
# meshTime : 'latestTime'
nProcs=24
scheme=Euler
# meshid
startTime=latestTime
endTime=4000
timeStep=0.05
writeInterval=40
purgeWrite=5
outputVBM=True
outputWave = False
outputMotions=True
# fsiTol = 1e-6
hullPatch=ship
donFile='IACS4400.don'
wave=None
waveH=1.0
waveT=10.0

#Forward velocity (in m/s)
velocity=0.

depth=500
# draft=0
# waveSeaLvl=0
waveStartTime=0
waveRampTime=10

EOF
'''

#*** Main execution start here *************************************************
if __name__ == "__main__":
    startTime = time.time()    
    dat = cmdOptions(sys.argv)
    foamCase_template(dat)
    copyMesh(dat)
    endTime = time.time()
    print 'Case generation completed in %d minutes' % ((endTime-startTime)/60)
    
