#!/usr/bin/env upython

#########################################################################
# Filename: fsMesher.py                                                 #
# Date:     2016-May-17                                                 #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:  Parametric mesh generation tool                             #
#########################################################################

import os, subprocess, time, math
import numpy as np
import sys, argparse, configparser
import pprint
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from fsTools import runCommand, foamFileExist, findBoundingBox, getFoamTimeFolders, findSTLPatches, translateStl, rotateStl, createBoxStl

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.blockMeshDict import BlockMeshDict
from inputFiles.refineMeshDict import RefineMeshDict
from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.compatOF import namePatch

# sose: 2016-june-03
# FIXME: There are still many conditions inwhich this script will fail
# FIXME: add a sanity check for input parameter(s)
# FIXME: add parallel support for refineMesh
# FIXME: adjust the number of box automatically
# FIXME: add user defined local refinement based on a stl-file
# FIXME: add user defined edge refinement
# FIXME: unify cell selection calls for a better performance
# FIXME: add selection of bad cells for further refinements
# FIXME: add selection of cells for coarsening
# FIXME: add selection of cells based on surface curvature
# FIXME: add restart options for debug mode
# FIXME: refineMesh fails to run in parallel (multiple times) perhaps problems with face-ordering of the processor patches

DEBUG = False
NPROCS = 4
DEFAULT_SHIP_STL = 'ship.stl'
DEFAULT_CFG_FILE = 'fsMesher.cfg'
DEFAULT_PARAMS = {
'heading' : 180,
'draft' : None,
'side' : 'port',
'domain' : [-3.0,2.5, -2.0,2.0, -1.5,0.5],
'LOA' : None,
'fs' : None,
'fsdZ' : None,
'fsCellRatio' : 4,
'refBoxType' : 'wave',
'refBoxData' : [3],
'refBoxGrad' : 3, 
'cellBuffer' : 4,
'refBow' : True,
'refStern' : True,
'refSurfExtra' : None,
'refFS' : True, 
'shipBL' : [3, 1.3, 0.7, 0.7],
'noLayers' : []
}

dataRefineProximity = {
'direction': 'xyz',
'stlFile': 'ship.stl',
'outsidePoints': [0.0, 0.0, 0.0],
'includeCutCells': True,
'includeInside': True,
'includeOutside': False,
'curvature': -1e6,
'distance': 0.0,
'boundingBox': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
}

# domain = [Xmin,Xmax, Ymin,Ymax, Zmin,Zmax]
# refBoxData = [n, Xmin,Xmax,Ymin,Ymax, ..., Xmin,Xmax,Ymin,Ymax, ... ]
# shipBL = [nLayers, layerGrowth, finalLayerThickness, minThicknessRatio]

#*** Data structure for user input data
class UserInput(object):
    def __init__(self, shipStl, opts=DEFAULT_PARAMS):
        # construct ship.stl and compute bounding box(es)
        if not isinstance(shipStl, list):
            shipStl = [shipStl]           
        filename = './constant/triSurface/'+DEFAULT_SHIP_STL
        # note: we always overwrite existing file, because we cannot check
        #       the validity of the existing file
        createShipStl(shipStl, filename, overwrite=True)
        shipPatches = findSTLPatches(filename)
        print("Found patches in stl: ", shipPatches)
        shipBB = findBoundingBox(filename)
        
        if opts['draft'] is not None:
            draft = math.fabs(float(opts['draft']))
            print("Set draft:",draft)
            move = -shipBB[2] - draft
            translateStl(filename, [0.0,0.0,move], filename)
            shipBB[2] += move
            shipBB[5] += move

        if opts['refSurfExtra'] is not None:
            nameOnly = os.path.basename(opts['refSurfExtra'])
            surfFile = "./constant/triSurface/"+nameOnly
            print("Create stl: "+surfFile)
            runCommand("cp -f "+opts['refSurfExtra']+" "+surfFile, CMD_showLog)
            if opts['draft'] is not None:
                translateStl(surfFile, [0.0,0.0,move], surfFile)
        
        if rotateStl(filename, opts['heading'], filename):
            shipBBRot = findBoundingBox(filename)
        else:
            shipBBRot = shipBB

        LOA = shipBB[3]-shipBB[0] if opts['LOA']==None else opts['LOA']
        fs = [0.02*LOA, 0.01*LOA] if opts['fs']==None else opts['fs']
        fsdZ = fs[1]/3 if opts['fsdZ']==None else opts['fsdZ']
        
        refBow=0
        if opts['refBow']==True:
            refBow=0.20*LOA

        refStern=0
        if opts['refStern']==True:
            refStern=0.20*LOA
        
        refFS = bool(opts['refFS'])
        
        domain = opts['domain']
        if opts['side']=='port':
            domain[2]=0
        elif opts['side']=='starboard':
            domain[3]=0
        domain = [ round(LOA*val, 3) for val in domain ]
        
        locationInMesh = [0.5*(domain[1]+shipBBRot[3]), 0.5*(domain[3]+domain[2]), 0.5*(domain[5]+domain[4])]

        # cut the cell ratio to the closest of 2^n
        fsCellRatio = opts['fsCellRatio']

        self.shipStl = list(shipStl)
        self.shipPatches = list(shipPatches)
        self.noLayers = opts['noLayers']
        self.shipBB = list(shipBB)
        self.shipBBRot = list(shipBBRot)
        self.draft = opts['draft']
        self.heading = opts['heading']
        self.side = opts['side']
        self.LOA = LOA
        self.domain = list(domain)
        self.fsZmin = fs[0]
        self.fsZmax = fs[1]
        self.fsdZ = fsdZ
        self.fsCellRatio = fsCellRatio
        self.refBoxType = opts['refBoxType']
        self.refBoxData = list(opts['refBoxData'])
        self.refBoxGrad = opts['refBoxGrad']
        self.refBoxZdata = None
        self.cellBuffer = opts['cellBuffer']
        self.refBow = refBow
        self.refStern = refStern
        self.refSurfExtra = opts['refSurfExtra']
        self.refFS = refFS
        self.cellWidth = None
        self.zCellSize = None
        self.zAllCut = None
        self.zAllCutNCells = None
        self.zAllCutRatio = None
        self.locationInMesh = list(locationInMesh)
        self.shipBL = list(opts['shipBL']) # [nLayers, layerGrowth, finalLayerThickness, minThicknessRatio]

        #    shipLow = 0             # Refinement level for snappyHexMesh
        #    shipHi = 0
        #    shipLayers = 3            # No. of layers for extrusion
        #    shipLayerGrowth = 1.3    # Growth ratio for layers
        #    shipFinalLayer = 0.7    # Final layer thickness
        #    shipMinThicknessRatio = 0.7    # ratio of the min. thickness to thinnest layer
   
def cmdOptions(argv):
    global DEBUG
    global EXEC_BLOCKMESH
    global EXEC_REFINEBOX
    global EXEC_REFINEPROXIMITY
    global EXEC_SNAP
    global EXEC_ADDLAYERS
    global CMD_keepLog
    global CMD_showLog
    global CMD_blockMesh
    global CMD_autoPatch
    global CMD_setSet
    global CMD_refineMesh
    global CMD_surfaceFeatureExtract
    global CMD_snappyHexMesh
    global CMD_decomposePar
    # default global parameters
    EXEC_BLOCKMESH = True
    EXEC_REFINEBOX = True
    EXEC_REFINEPROXIMITY = True
    EXEC_SNAP = True
    EXEC_ADDLAYERS = True
    CMD_keepLog = ' >> ./log.fsMesher 2>&1 '
    CMD_showLog = 'echo "\nlog-file: ./log.fsMesher"; tail -45 ./log.fsMesher; echo "Please see log-file: ./log.fsMesher"'
    CMD_blockMesh = 'blockMesh'
    CMD_autoPatch = 'autoPatch -overwrite 80'
    CMD_setSet = 'setSet -latestTime'
    CMD_refineMesh = 'refineMesh'
    CMD_surfaceFeatureExtract = 'surfaceFeatureExtract'
    CMD_snappyHexMesh = 'snappyHexMesh'
    CMD_decomposePar = 'decomposePar -force -latestTime'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--debug', action='store_true', help='run in debug mode')
    parser.add_argument('-f', dest='inputfile', help='read parameters from input file. Use this option with a stl-file to generate a mesh with default parameters, e.g.: fsMesher.py -f myship.stl ')
    parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print(all available input parameters (including theirs default values) to fsMesher.cfg file and exit. Use this option to generate a template for input files for -f')
    parser.add_argument('-n','--no-exec', dest='skip', action='store_true', help='do not exec. any foam tool')
    parser.add_argument('-r','--run', help='run only selected steps. Use comma to select multiple steps, e.g: "fyMesher.py -r bg,box,snap".  Available selections are \
                                            [bg]:blockMesh, \
                                            [box]:refineBox, \
                                            [prox]:refineProximity, \
                                            [snap]:snappyHexMesh, \
                                            [layer]:snappyHexMesh')
    parser.add_argument('-s','--skip', help='skip selected steps, same arguments as for -r. \
                                            Do not use -s and -r options simultaneously')
    parser.add_argument('-clz','--clear-zones', dest='clearZones', action='store_true', help='Clear empty {point,cell,face}Zones files and exit. With -n option no file-removal action will be made.')
    args = parser.parse_args()

    if args.showConfig:
        print("\nOutput default parameters to file:", DEFAULT_CFG_FILE)
        subprocess.call('cat << EOF > '+DEFAULT_CFG_FILE + defaultParams_contents, shell=True)
        raise SystemExit('')

    if args.clearZones:
        if args.skip==True:
            clearEmptyZonesFiles(dryrun=True)
        else:
            clearEmptyZonesFiles()
        raise SystemExit('')
    
    if args.inputfile==None:
        inputdata = UserInput(DEFAULT_SHIP_STL)
    else:
        fid = str(args.inputfile)
        if fid.endswith('.stl'):
            inputdata = UserInput(fid)
        else:
            # here we read all input parameters from files
            shipStl,params = readInputParams(fid)
            inputdata = UserInput(shipStl, opts=params)
            pass    
    
    if args.skip==True:
        EXEC_BLOCKMESH = False
        EXEC_REFINEBOX = False
        EXEC_REFINEPROXIMITY = False
        EXEC_SNAP = False
        EXEC_ADDLAYERS = False
    elif not args.skip==False:
        skip = str(args.skip).split(',')
        for val in skip:
            if val=="bg":
                EXEC_BLOCKMESH=False
            elif val=="box":
                EXEC_REFINEBOX=False
            elif val=="prox":
                EXEC_REFINEPROXIMITY=False
            elif val=="snap":
                EXEC_SNAP=False
            elif val=="layer":
                EXEC_ADDLAYERS=False
            else:
                print("\nUnknown arg.: --skip", val,"\n")
                parser.print_help()
                raise SystemExit()

    if not args.run==None:
        EXEC_BLOCKMESH = False
        EXEC_REFINEBOX = False
        EXEC_REFINEPROXIMITY = False
        EXEC_SNAP = False
        EXEC_ADDLAYERS = False
        run = str(args.run).split(',')
        for val in run:
            if val=="bg":
                EXEC_BLOCKMESH=True
            elif val=="box":
                EXEC_REFINEBOX=True
            elif val=="prox":
                EXEC_REFINEPROXIMITY=True
            elif val=="snap":
                EXEC_SNAP=True
            elif val=="layer":
                EXEC_ADDLAYERS=True
            else:
                print("\nUnknown arg.: --run", val,"\n")
                parser.print_help()
                raise SystemExit()

    if args.debug:
        DEBUG=True

    if not DEBUG:
        CMD_setSet += ' -noVTK'
        CMD_refineMesh += ' -overwrite'
        #CMD_snappyHexMesh += ' -overwrite' # snappyHexMesh will be run in parallel!!! FIXME
    else:
        print("Running in DEBUG mode ... ")

    subprocess.call('echo "Executing fsMesh.py: $(hostname) @ $(date)" '+CMD_keepLog, shell=True)
    return inputdata
    pass

def readInputParams(filename):
    params = dict(DEFAULT_PARAMS)
    config = configparser.ConfigParser()
    config.read(filename)
    name = 'fsMesher'
    
    print("Read input parameters from file:", filename)
    
    # read stlFile(s)
    done=False
    shipStl = []
    i = 0
    while not done:
        try:
            txt = str(config[name]['stlFile'+str(i)])
            shipStl.append(txt)
            i += 1
        except KeyError:
            done=True
            pass
    if i is not 0:
        print("    stlFile(s):", shipStl)

    # heading
    try:
        txt = str(config[name]['draft'])
        if not txt=='None':
            params['draft'] = math.fabs(float(txt))
    except KeyError: pass

    # heading
    try:
        txt = str(config[name]['heading'])
        params['heading'] = float(txt) % float(360)
    except KeyError: pass

    # side = 'port' | 'starboard' | 'both'
    try:
        txt = str(config[name]['side'])
        params['side'] = txt
    except KeyError: pass

    # domain = [X(min,max), Y(min,max), Z(min,max)] of LOA
    try:
        txt = str(config[name]['domain'])
        if not txt=='auto':
            params['domain'] = [float(val) for val in txt.split(',')]
    except KeyError: pass

    # LOA : automatic from stl
    try:
        txt = str(config[name]['LOA'])
        if not txt=='auto':
            params['LOA'] = float(txt)
    except KeyError: pass
    
    # fs : [fsZmin fsZmax] free surface zone
    try:
        txt = str(config[name]['fsZone'])
        if not txt=='auto':
            params['fs'] = [float(val) for val in txt.split(',')]
    except KeyError: pass

    # 'fsdZ' : cell height inside free surface zone 
    try:
        txt = str(config[name]['fsCellHeight'])
        if not txt=='auto':
            params['fsdZ'] = float(txt)
    except KeyError: pass

    # 'fsCellRatio' : cell ratio (length/height) ratio inside free surface zone 
    try:
        txt = str(config[name]['fsCellRatio'])
        if not txt=='auto':
            params['fsCellRatio'] = float(txt)
    except KeyError: pass
    
    # 'refBoxType' : 'wave', 'kelvin', 'both'
    try:
        txt = str(config[name]['refBoxType'])
        if not txt=='auto':
            params['refBoxType'] = txt
    except KeyError: pass
    
    # 'refBoxData' : [number, xmin,xmax,ymin,ymax, ..., ], 
    try:
        txt = str(config[name]['refBoxData'])
        if not txt=='auto':
            params['refBoxData'] = [float(val) for val in txt.split(",")]
    except KeyError: pass

    # 'refBoxGrad' : box size ratio (farfield/nearField)
    try:
        txt = str(config[name]['refBoxRatio'])
        if not txt=='auto':
            # params['refBoxGrad'] = [float(val) for val in txt.split(",")]
            params['refBoxGrad'] = float(txt)
    except KeyError: pass

    # 'cellBuffer' : buffer between nabo cells
    try:
        txt = str(config[name]['cellBuffer'])
        if not txt=='auto':
            params['cellBuffer'] = int(txt)
    except KeyError: pass

    # 'refBow' : scalar value, refine bow or not?
    try:
        txt = str(config[name]['refineBow'])
        if not txt=='auto':
            params['refBow'] = float(txt)
    except KeyError: pass

    # 'refStern' : scalar value, refine stern or not?
    try:
        txt = str(config[name]['refineStern'])
        if not txt=='auto':
            params['refStern'] = float(txt)
    except KeyError: pass

    # 'refSurfExtra' : stl-file for extra refinement
    try:
        txt = str(config[name]['refineSurf'])
        if txt is not None:
            params['refSurfExtra'] = txt
    except KeyError: pass

    # 'refFS' : bool, refine free surface in nearField ?
    try:
        txt = config[name].getboolean('refineFS')
        if (not txt=='auto') & (not txt==None):
            params['refFS'] = txt
    except KeyError: pass

    # 'shipBL' : boundary layer, [nLayers, layerGrowth, finalLayerThickness, minThicknessRatio] 
    try:
        txt = config[name]['layers']
        if not txt=='auto':
            params['shipBL'] = [float(val) for val in txt.split(',')]
    except KeyError: pass

    # 'shipBL' : boundary layer, [nLayers, layerGrowth, finalLayerThickness, minThicknessRatio] 
    try:
        txt = config[name]['disableLayers']
        if not txt==None:
            params['noLayers'] = [str(val).strip() for val in txt.split(',')]
    except KeyError: pass

    #
    global DEBUG
    global NPROCS
    global DEFAULT_SHIP_STL
    global EXEC_BLOCKMESH
    global EXEC_REFINEBOX
    global EXEC_REFINEPROXIMITY
    global EXEC_SNAP
    global EXEC_ADDLAYERS
    global CMD_keepLog
    global CMD_showLog
    global CMD_blockMesh
    global CMD_autoPatch
    global CMD_setSet
    global CMD_refineMesh
    global CMD_surfaceFeatureExtract
    global CMD_snappyHexMesh
    global CMD_decomposePar

    txt = config['fsMesher-control'].getboolean('DEBUG')
    if not txt=='None': DEBUG = txt
    txt = config['fsMesher-control'].getboolean('EXEC_BLOCKMESH')
    if not txt=='None': EXEC_BLOCKMESH = txt
    txt = config['fsMesher-control'].getboolean('EXEC_REFINEBOX')
    if not txt=='None': EXEC_REFINEBOX = txt
    txt = config['fsMesher-control'].getboolean('EXEC_REFINEPROXIMITY')
    if not txt=='None': EXEC_REFINEPROXIMITY = txt
    txt = config['fsMesher-control'].getboolean('EXEC_SNAP')
    if not txt=='None': EXEC_SNAP = txt
    txt = config['fsMesher-control'].getboolean('EXEC_ADDLAYERS')
    if not txt=='None': EXEC_ADDLAYERS = txt

    try:
        txt = config['fsMesher-control']['NPROCS']
        if not txt == None:
            NPROCS = int(txt)
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['DEFAULT_SHIP_STL']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        DEFAULT_SHIP_STL = txt
    except KeyError: pass

    try:
        txt = config['fsMesher-control']['CMD_keepLog']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_keepLog = txt
    except KeyError: pass

    try:
        txt = config['fsMesher-control']['CMD_showLog']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_showLog = txt
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_blockMesh']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_blockMesh = txt
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_autoPatch']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_autoPatch = txt
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_setSet']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_setSet = txt
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_refineMesh']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_refineMesh = txt
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_surfaceFeatureExtract']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_surfaceFeatureExtract = txt
    except KeyError: pass

    try:
        txt = config['fsMesher-control']['CMD_snappyHexMesh']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_snappyHexMesh = txt
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_decomposePar']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_decomposePar = txt
    except KeyError: pass    

    DEFAULT_SHIP_STL=os.path.basename(str(DEFAULT_SHIP_STL))
    return shipStl, params
    pass

def createShipStl(shipStl, outputStl, overwrite=False):
    print("Creating stl: "+outputStl)
    if os.path.isfile(outputStl) & (not overwrite):
        print("    file already exist ... reuse existing file")
        return

    subprocess.call("rm -fr "+outputStl, shell=True)
    if len(shipStl)==0:
        print("    stlFile not specified ...")
        raise SystemExit('abort ...')
    foo = 0
    while foo < len(shipStl):
        stlFile = shipStl[foo]
        foo += 1
        #if not stlFile[::-1][0:4][::-1]==".stl":
        if not stlFile.endswith('.stl'):
            stlFile = stlFile + '.stl'
        print("    append: " + stlFile)
        if not os.path.isfile(stlFile):
            print("\nFile not found:",stlFile,"\n")
            raise SystemExit('abort ...')
        subprocess.call("mkdir -p $(dirname "+outputStl+")", shell=True)
        runCommand("cat "+ stlFile + " >> "+outputStl, CMD_showLog)

# return a list of numbers btw. 0 .. 1
# dNd1_Ratio=dN/d1
# x[0]=0, x[N]=1
# x[1] = 1/sum(c^i)_(i=0..N-1) = (1-c)/(1-c^N) (for all c != 1)
# x[n] = x1*sum(c^i)_(i=0..n-1) = x1*(1-c^n)/(1-c) (for all c != 1, n=1..N)
# c = dNd1_Ratio^(1/(N-1))
def simpleGrading(N, dNd1_Ratio):
    x = [0,1]
    if N <= 1:
        return [0,1]
    c = math.pow(dNd1_Ratio, 1.0/(N-1))
    if math.fabs(c-1.0)<1e-6:
        x = np.linspace(0,1,N+1)
        return x.tolist()
    x[1]=(1-c)/(1-math.pow(c,float(N)))
    for n in range(2,N):
        x.append(x[1]*(1-math.pow(c,float(n)))/(1.0-c))
    x.append(1)
    return x

# given x1 compute N which produce the closest dNd1_Ratio 
def simpleGradingN(x1, dNd1_Ratio):
    x1pre,N=1.0,0
    while x1pre>x1:
        N += 1
        x = simpleGrading(N, float(dNd1_Ratio))
        if x[1]<=x1:
            if math.fabs(x1pre-x1) < math.fabs(x[1]-x1):
                N -= 1
            break
        else:
            x1pre=x[1]
    return N

def createBlockMeshDict(data):
    # note: blockMesh simeplGradient is defined as the ratio last/first cell
    # so we have,
    #   N : given number of cells
    #   r : given ration r = dxN/dx1
    #   c : cell expansion ratio c = r^(1/(N-1))
    #   dx1 : size of the first cell dx1 = 1 / sum(c^i)_(i=0..N-1)
    #   dxN : size of the last cell dxN = dx1*r
    #
    #   in a reverse problem we have:
    #   dx1 : given size of the first or last cell
    #   r : given ratio r = dxN/dx1
    #   N : the closest number of cell is computed from N = minimum( r - sum(r^(i/(N-1))_(i=0..N-1) )

    cellBuffer = data.cellBuffer  # cells between each level 
    XminDomain = data.domain[0]
    XmaxDomain = data.domain[1]
    YminDomain = data.domain[2]
    YmaxDomain = data.domain[3]
    ZminDomain = data.domain[4]
    ZmaxDomain = data.domain[5]
    ZminSurface = data.fsZmin
    ZmaxSurface = data.fsZmax

    shipBBxMin = data.shipBBRot[0]
    shipBByMin = data.shipBBRot[1]
    shipBBzMin = data.shipBBRot[2]
    shipBBxMax = data.shipBBRot[3]
    shipBByMax = data.shipBBRot[4]
    shipBBzMax = data.shipBBRot[5]
        
    # how many refinement box? minimum is 1
    refBoxData = list(data.refBoxData)
    nRefBox = int(refBoxData[0])
    
    # The inner cell size is defined by fsdZ
    dZ = [ 0 for i in range(nRefBox+2)]
    dZ[nRefBox+1] = round(data.fsdZ, 5)
    for i in range(nRefBox, -1, -1):
        dZ[i] = dZ[i+1]*2  

    # keep is in a reverse order: i.e. level [n+1, n, n-1, n-2,..., 0]
    dZ = dZ[::-1]
    data.zCellSize = list(dZ)

    # horizontal cells for the whole domain    
    cellWidth = dZ[-2]*data.fsCellRatio
    Xcells = int((XmaxDomain-XminDomain)/cellWidth)
    Ycells = int((YmaxDomain-YminDomain)/cellWidth)
    data.cellWidth = cellWidth

    # update domain data
    data.domain[0]= round(XmaxDomain - cellWidth*Xcells, 3)
    XminDomain = data.domain[0]
    if data.side=='port':
        data.domain[3] = round(YminDomain + cellWidth*Ycells, 3)
        YmaxDomain = data.domain[3]
    elif data.side=='starboard':
        data.domain[2] = round(YmaxDomain - cellWidth*Ycells, 3)
        YminDomain = data.domain[2]
    elif data.side=='both':
        if (Ycells % 2 != 0):
            Ycells += 1 # we need this to be even
        data.domain[2] = -round(cellWidth*Ycells/2, 3)
        data.domain[3] = round(cellWidth*Ycells/2, 3)
        YminDomain = data.domain[2]
        YmaxDomain = data.domain[3]
    else:
        print("\nUnknown parameters, side=",data.side)
        raise SystemExit('abort ...')
    
    # compute vertical thickness of each box
    # All boxes must be larger than the ship's bounding box
    fsCellTop = int(math.ceil(round(math.fabs(ZmaxSurface)/(data.fsdZ),1)))
    fsCellBottom = int(math.ceil(round(math.fabs(ZminSurface)/(data.fsdZ),1)))
    data.fsZmax = fsCellTop*data.fsdZ
    data.fsZmin = -fsCellBottom*data.fsdZ
    ZmaxSurface = data.fsZmax
    ZminSurface = data.fsZmin

    # z-cuts in lower block
    lowerCutNCells = cellBuffer
    lowerCut = ZminSurface - dZ[1]*lowerCutNCells;
    while (lowerCut>(shipBBzMin-cellBuffer*dZ[1])) | (lowerCutNCells < cellBuffer):
        lowerCut -= dZ[1]
        lowerCutNCells += 1
    lowerCutNCells = [lowerCutNCells]
    lowerCut = [lowerCut]   

    #for i in range(1, nRefBox):
    #    lowerCutNCells.append(cellBuffer)
    #    lowerCut.append(lowerCut[i-1] - cellBuffer*dZ[i+1])

    # z-cuts in upper block
    upperCutNCells = cellBuffer
    upperCut = ZmaxSurface + dZ[1]*upperCutNCells;
    while (upperCut<(shipBBzMax+cellBuffer*dZ[1])) | (upperCutNCells < cellBuffer):
        upperCut += dZ[1]
        upperCutNCells += 1
    upperCutNCells = [upperCutNCells]
    upperCut = [upperCut]

    #for i in range(1, nRefBox):
    #    upperCutNCells.append(cellBuffer)
    #    upperCut.append(upperCut[i-1] + cellBuffer*dZ[i+1])

    # compute grading at the top/bottom block(s)
    ZratioTop = cellWidth/dZ[-nRefBox-1]
    dx1 = dZ[-nRefBox-1]/math.fabs(ZmaxDomain-upperCut[-1])
    ZcellsTop = simpleGradingN(dx1, ZratioTop)
    ZratioBottom = cellWidth/dZ[-nRefBox-1]
    dx1 = dZ[-nRefBox-1]/math.fabs(ZminDomain-lowerCut[-1])
    ZcellsBottom = simpleGradingN(dx1, ZratioBottom)
    ZratioBottom = 1.0/ZratioBottom

    # collect all z-cut
    zAllCut = lowerCut[::-1] + [ZminSurface, 0, ZmaxSurface] + upperCut + [ZmaxDomain]
    zAllCutNCells = [ZcellsBottom] + lowerCutNCells[::-1] + [fsCellBottom, fsCellTop] + upperCutNCells + [ZcellsTop]
    zAllCutRatio = [ZratioBottom] + list(1 for i in range(len(lowerCut)*2+2)) + [ZratioTop]

    # compute vertical position of all grid points 
    zGrid = [ZminDomain]
    for i in range(len(zAllCut)):
        grad = simpleGrading(zAllCutNCells[i], zAllCutRatio[i])
        L = zAllCut[i] - zGrid[-1]
        zGrid = zGrid + [zGrid[-1]+L*grad[val] for val in range(1,len(grad)-1)]
        zGrid.append(zAllCut[i])
    zGridDelta = [zGrid[i+1]-zGrid[i] for i in range(0,len(zGrid)-1)]
    
    # compute vertical extension for all refBox
    dx = cellWidth/2.0
    refBoxZdata = [zGrid[0]]*nRefBox + [zGrid[-1]]*nRefBox
    for i in range(nRefBox):
        for j in range(len(zGridDelta)):
            if zGrid[j] >= ZminSurface:
                break
            if zGridDelta[j] < dx*1.25:
                refBoxZdata[i] = zGrid[j]
                break
        for j in range(len(zGridDelta)):
            if zGrid[-j-1] <= ZmaxSurface:
                break
            if zGridDelta[-j-1] < dx*1.25:
                refBoxZdata[-i-1] = zGrid[-j-1]
                break
        dx /= 2.0
    
    print('Domain bounding box:')
    print("   ", [XminDomain, YminDomain, ZminDomain, XmaxDomain, YmaxDomain, ZmaxDomain])

    polyfolder = os.path.join("constant","polyMesh")
    if not os.path.exists(polyfolder): os.makedirs(polyfolder)
    
    #Write blockMeshDict file
    blockMeshDict = BlockMeshDict( case     = '.',
                                   ndim     = 3,
                                   waveMesh = True,
                                   xmin     = XminDomain,
                                   xmax     = XmaxDomain,
                                   ymin     = YminDomain,
                                   ymax     = YmaxDomain,
                                   zmin     = ZminDomain,
                                   zmax     = zAllCut,
                                   Xcells   = Xcells,
                                   Ycells   = Ycells,
                                   Zcells   = zAllCutNCells,
                                   Zgrading = zAllCutRatio)
    blockMeshDict.writeFile()
    
    # update data
    data.zAllCut = zAllCut
    data.zAllCutNCells = zAllCutNCells
    data.zAllCutRatio = zAllCutRatio
    data.refBoxZdata = refBoxZdata

    # compute x,y data for refBox
    if (len(refBoxData) == 1):
        grad = simpleGrading(nRefBox+1, data.refBoxGrad)
        xGradMin = [math.fabs(val-1.0) for val in grad[::-1]]
        xGradMax = grad
        yGradMin = [math.fabs(val-1.0) for val in grad[::-1]]
        yGradMax = grad
        #
        xCutMin = [(shipBBxMin-XminDomain)*val for val in xGradMin]
        del xCutMin[0], xCutMin[-1]
        for i in range(len(xCutMin)):
            xCutMin[i] = cellWidth*int(round(xCutMin[i]/cellWidth))
        xCutMin = [XminDomain + val for val in xCutMin]
        #
        if (data.side=='port'):
            yCutMin = [-1e-3 for val in xCutMin]
        else:
            yCutMin = [(shipBByMin-YminDomain)*val for val in yGradMin]
            del yCutMin[0], yCutMin[-1]
            for i in range(len(yCutMin)):
                yCutMin[i] = cellWidth*int(round(yCutMin[i]/cellWidth))
            yCutMin = [YminDomain + val for val in yCutMin]
        #
        if data.refBoxType=='wave':
            xCutMax = [XmaxDomain+1e-3 for val in xCutMin]
        elif data.refBoxType=='kelvin':
            print("not yet implemented")
        else:
            xCutMax = [(XmaxDomain-shipBBxMax)*val for val in xGradMax]
            del xCutMax[0], xCutMax[-1]
            for i in range(len(xCutMax)):
                xCutMax[i] = cellWidth*int(round(xCutMax[i]/cellWidth))
            xCutMax = [shipBBxMax + val for val in xCutMax]
        #
        if (data.side=='starboard'):
            yCutMax = [1e-3 for val in xCutMin]
        else:
            yCutMax = [(YmaxDomain-shipBByMax)*val for val in yGradMax]
            del yCutMax[0], yCutMax[-1]
            for i in range(len(yCutMax)):
                yCutMax[i] = cellWidth*int(round(yCutMax[i]/cellWidth))
            yCutMax = [shipBByMax + val for val in yCutMax]
            yCutMax = yCutMax[::-1]
        #
        # update data.refBoxData
        for i in range(len(xCutMin)):
            data.refBoxData.append(xCutMin[i])
            data.refBoxData.append(xCutMax[i])
            data.refBoxData.append(yCutMin[i])
            data.refBoxData.append(yCutMax[i])
    #

def createBackGroundMesh(data):
    if EXEC_BLOCKMESH:
        print("blockMesh: create a base mesh ...")
        runCommand(CMD_blockMesh + CMD_keepLog, CMD_showLog)
        print("autoPatch: create domain boundaries ...")
        runCommand(CMD_autoPatch + CMD_keepLog, CMD_showLog)
        
        # Rename boundary patches
        boundfile = os.path.join("constant","polyMesh","boundary")
        boundDict = ParsedParameterFile(boundfile,boundaryDict=True)
        nbound = int(len(boundDict)/2)
        for i in range(nbound):
            if boundDict[2*i] in namePatch["foamStar"]:
                boundDict[2*i] = namePatch["foamStar"][boundDict[2*i]]
                if boundDict[2*i] == 'domainZ0':
                    boundDict[2*i+1]['type'] = 'wall'
                elif boundDict[2*i] == 'domainY0' and data.side=='port':
                    boundDict[2*i+1]['type'] = 'symmetryPlane'
                elif boundDict[2*i] == 'domainY1' and data.side=='starboard':
                    boundDict[2*i+1]['type'] = 'symmetryPlane'
        boundDict.writeFile()
        
#        if NPROCS>1:
#            print("decomposePar: nProcs =",NPROCS)
#            runCommand(CMD_decomposePar + CMD_keepLog, CMD_showLog)
    else:
        print("blockMesh: ... skip")

    # how many refinement box? minimum is 1
    refBoxData = list(data.refBoxData)
    nRefBox = int(refBoxData[0])
    del refBoxData[0]
    refBoxBB = []
    if len(data.refBoxData)>1:
        if nRefBox != (len(data.refBoxData)-1)/4.:
            raise SystemExit('Error: invalid data for refinement boxes, ', data.refBoxData)
        for i in range(nRefBox):
            refBoxBB.append([refBoxData[i*4], refBoxData[i*4+2], data.refBoxZdata[i], refBoxData[i*4+1], refBoxData[i*4+3], data.refBoxZdata[-i-1]])
    else:
        raise SystemExit('\nData for refinement box is missing.\nrefBoxData=[#n, #xmin,#xmax,#ymax,#ymax, #xmin,#xmax,#ymin,#ymax, ..., repeat n times]\nabort ...')

    shipBBxMin = data.shipBBRot[0]
    shipBByMin = data.shipBBRot[1]
    shipBBzMin = data.shipBBRot[2]
    shipBBxMax = data.shipBBRot[3]
    shipBByMax = data.shipBBRot[4]
    shipBBzMax = data.shipBBRot[5]

    # number of level(s) for uniform 'xy'-refinement
    nxy = int(math.log(math.fabs(float(data.fsCellRatio)), float(2))) - 1

    # refine free surface (using proximity method)
    distance = data.fsdZ*data.cellBuffer*2.0*(nxy-1.0 + nRefBox-1.0 + float(bool(data.refBow) | bool(data.refStern) | bool(data.refFS)))
    
    lastInnerBox = refBoxBB[-1]
    if lastInnerBox[0] > (shipBBxMin - distance - data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox))):
        diff = lastInnerBox[0] - (shipBBxMin - distance - data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox)))
        for i in range(len(refBoxBB)):
            refBoxBB[i][0] -= diff
    if lastInnerBox[1] > (shipBByMin - distance - data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox))):
        diff = lastInnerBox[1] - (shipBByMin - distance - data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox)))
        for i in range(len(refBoxBB)):
            refBoxBB[i][1] -= diff
    if lastInnerBox[2] > (shipBBzMin - distance - data.zCellSize[1]*data.cellBuffer):
        diff = lastInnerBox[2] - (shipBBzMin - distance - data.zCellSize[1]*data.cellBuffer)
        for i in range(len(refBoxBB)):
            refBoxBB[i][2] -= diff
    if lastInnerBox[3] < (shipBBxMax + distance + data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox))):
        diff = -lastInnerBox[3] + shipBBxMax + distance + data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox))
        for i in range(len(refBoxBB)):
            refBoxBB[i][3] += diff
    if lastInnerBox[4] < (shipBByMax + distance + data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox))):
        diff = -lastInnerBox[4] + shipBByMax + distance + data.cellWidth*data.cellBuffer/math.pow(2.0, float(nRefBox))
        for i in range(len(refBoxBB)):
            refBoxBB[i][4] += diff
    if lastInnerBox[5] < (shipBBzMax + distance + data.zCellSize[1]*data.cellBuffer):
        diff = -lastInnerBox[5] + shipBBzMax + distance + data.zCellSize[1]*data.cellBuffer
        for i in range(len(refBoxBB)):
            refBoxBB[i][5] += diff

    #print("debug: refBoxData: ", data.refBoxData)
    #for BB in refBoxBB:
    #    refineBox(BB, 'xy')
    #    pass    

    for BB in refBoxBB:
        tmp = BB[4]
        BB[4]=data.domain[3] # YmaxDomain
        refineBox(BB, 'x')
        BB[4] = tmp
        pass    

    for BB in refBoxBB:
        refineBox(BB, 'y')
        pass    
    
    # this is a point outside ship.stl
    outsidePoints = [0.5*(shipBBxMin+shipBBxMax), 0.5*(shipBByMin+shipBByMax)+math.fabs(shipBByMin-shipBByMax), 0.5*(shipBBzMin+shipBBzMax)]

    for i in range(nxy):
        selectProximity('new', DEFAULT_SHIP_STL, distance, outsidePoints=outsidePoints)
        refineProximity('xy')
        distance *= 0.5

    BB = [-1e6,-1e6,data.fsZmin,1e6,1e6,data.fsZmax]
    selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
    refineProximity('xy')
    BB = [-1e6,-1e6,-1e6,1e6,1e6,data.fsZmin]
    selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
    refineProximity('xyz')
    BB = [-1e6,-1e6,data.fsZmax,1e6,1e6,1e6]
    selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
    refineProximity('xyz')
    #
    XmaxDomain = data.domain[1]
    # align cutting locations
    refBow = data.refBow if bool(data.refBow) else 0.2*(shipBBxMax-shipBBxMin)
    refBow = shipBBxMax-refBow 
    refBow = XmaxDomain-math.ceil((XmaxDomain-refBow)/data.fsdZ)*data.fsdZ
    if bool(data.refBow):
        distance *= 0.5
        BB = [refBow,-1e6,-1e6,1e6,1e6,shipBBzMax-data.fsdZ*data.cellBuffer]
        selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('xyz')

    if data.refSurfExtra is not None:
        nameOnly = os.path.basename(data.refSurfExtra)
        BB = [refBow+0.5*data.fsdZ*data.cellBuffer,-1e6,-1e6,1e6,1e6,shipBBzMax-1.5*data.fsdZ*data.cellBuffer]
        selectProximity('new', nameOnly, 0.5*distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('xyz')

    # align cutting locations
    refStern = data.refStern if bool(data.refStern) else 0.2*(shipBBxMax-shipBBxMin)
    refStern = shipBBxMin+refStern
    refStern = XmaxDomain-math.floor((XmaxDomain-refStern)/data.fsdZ)*data.fsdZ        
    if bool(data.refStern):
        if not bool(data.refBow):
            distance *= 0.5
        BB = [-1e6,-1e6,-1e6,refStern,1e6,shipBBzMax-data.fsdZ*data.cellBuffer]
        selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('xyz')
    if data.refSurfExtra is not None:
        nameOnly = os.path.basename(data.refSurfExtra)
        BB = [-1e6,-1e6,-1e6,refStern-0.5*data.fsdZ*data.cellBuffer,1e6,shipBBzMax-1.5*data.fsdZ*data.cellBuffer]
        selectProximity('new', nameOnly, 0.5*distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('xyz')

    if bool(data.refFS) & (bool(data.refBow) | bool(data.refStern)):
        BB = [refStern,-1e6,data.fsZmin,refBow,1e6,data.fsZmax]
        selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('z')

def createSnappyMesh(data):

    #surfaceFeatureExtract (if not already done)
    sltname = DEFAULT_SHIP_STL.split('.stl')[0] 
    if not foamFileExist('./constant/triSurface/'+sltname+'.eMesh'):
        print("\nExtract surface features from file: ./constant/triSurface/"+DEFAULT_SHIP_STL)
        surfaceFeatureExtractDict = SurfaceFeatureExtractDict(case = '.',
                                                              stlname = DEFAULT_SHIP_STL)
        surfaceFeatureExtractDict.writeFile()
        runCommand(CMD_surfaceFeatureExtract + CMD_keepLog, CMD_showLog)
    else:
        print("Reuse existing file: "+'./constant/triSurface/'+sltname+'.eMesh')
    
    #snappyHexMesh
    if EXEC_SNAP:
        snappyHexMeshDict = SnappyHexMeshDict(case                       = '.',
                                              stlname                    = sltname,
                                              castellatedMesh            = True,
                                              snap                       = True,
                                              addLayers                  = False,
                                              locationInMesh             = data.locationInMesh,
                                              nCellsBetweenLevels        = 1,
                                              edgeLvl                    = 0,
                                              hullLvl                    = [0,0],
                                              resolveFeatureAngle        = 15,
                                              allowFreeStandingZoneFaces = False,
                                              snapTol                    = 0.75,
                                              nSolveIter                 = 100,
                                              maxNonOrtho                = 65,
                                              minTwist                   = 0.02,
                                              nSmoothScale               = 5,
                                              errorReduction             = 0.75)
        snappyHexMeshDict.writeFile()
    
        if NPROCS>1:
            isDecomposed,nProcs = caseAlreadyDecomposed()
            if not isDecomposed:
                nProcs = NPROCS
                print("decomposePar: nProcs =",nProcs)
                runCommand(CMD_decomposePar + CMD_keepLog, CMD_showLog)
            print("snappyHexMesh: snapping ... in parallel, nProcs =",nProcs)
            cmd = "mpirun -np "+str(nProcs)+" "+CMD_snappyHexMesh + " -parallel "+ CMD_keepLog
        else:
            print("snappyHexMesh: snapping ... ")
            cmd = CMD_snappyHexMesh + CMD_keepLog
        runCommand(cmd, CMD_showLog)
    else:
        print("snappyHexMesh: snap ... skip")
        pass
        
    if EXEC_ADDLAYERS:
        nLayers = data.shipBL[0]
        layerGrowth = data.shipBL[1]
        finalLayerThickness = data.shipBL[2]
        minThicknessRatio = data.shipBL[3]
        minThickness = float(minThicknessRatio)*finalLayerThickness/(pow(layerGrowth,float(nLayers-1)))
    
        snappyHexMeshDict = SnappyHexMeshDict(case                       = '.',
                                              stlname                    = sltname,
                                              castellatedMesh            = False,
                                              snap                       = False,
                                              addLayers                  = True,
                                              relativeSizes              = True,
                                              nSurfaceLayers             = nLayers,
                                              expansionRatio             = layerGrowth,
                                              finalLayerThickness        = finalLayerThickness,
                                              minThickness               = minThickness,
                                              featureAngle               = 60,
                                              shipPatches                = data.shipPatches,
                                              noLayers                   = data.noLayers,
                                              maxNonOrtho                = 65,
                                              minTwist                   = 0.02,
                                              nSmoothScale               = 5,
                                              errorReduction             = 0.75)
        snappyHexMeshDict.writeFile()
    
        if NPROCS>1:
            isDecomposed,nProcs = caseAlreadyDecomposed()
            if not isDecomposed:
                nProcs = NPROCS
                print("decomposePar: nProcs =",nProcs)
                runCommand(CMD_decomposePar + CMD_keepLog, CMD_showLog)
            print("snappyHexMesh: add layers ... in parallel, nProcs =",nProcs)
            cmd = "mpirun -np "+str(nProcs)+" "+CMD_snappyHexMesh + " -parallel "+ CMD_keepLog
        else:
            print("snappyHexMesh: add layers ... ")
            cmd = CMD_snappyHexMesh + CMD_keepLog
        runCommand(cmd, CMD_showLog)
    else:
        print("snappyHexMesh: add layers ... skip")
        pass
    pass

def run_setSet():
    # FIXME: refineMesh is buggy when run in parallel
    runCommand(CMD_setSet + ' -batch .tmp_setSet ' + CMD_keepLog, CMD_showLog)

    #if NPROCS==1:
    #    runCommand(CMD_setSet + ' -batch .tmp_setSet ' + CMD_keepLog, CMD_showLog)
    #else:
    #    runCommand('mpirun -np '+str(NPROCS)+' '+ CMD_setSet + ' -parallel -batch .tmp_setSet ' + CMD_keepLog, CMD_showLog)
    pass

def run_refineMesh():
    # FIXME: refineMesh is buggy when run in parallel
    runCommand(CMD_refineMesh + CMD_keepLog, CMD_showLog)
    
    #if NPROCS==1:
    #    runCommand(CMD_refineMesh + CMD_keepLog, CMD_showLog)
    #else:
    #    runCommand('mpirun -np '+str(NPROCS)+' '+CMD_refineMesh + ' -parallel ' + CMD_keepLog, CMD_showLog)
    pass

def selectBoxToCell(BB):
    xmin,ymin,zmin,xmax,ymax,zmax=BB[0],BB[1],BB[2],BB[3],BB[4],BB[5]
    BBtxt = '('+str(xmin)+' '+str(ymin)+' '+str(zmin)+') ('+str(xmax)+' '+str(ymax)+' '+str(zmax)+')'
    cmd = 'cellSet c0 new boxToCell '+BBtxt
    print("SelectBoxToCell:",cmd)
    if EXEC_REFINEBOX:
        runCommand('echo "'+cmd+'" > .tmp_setSet', CMD_showLog)
        run_setSet()
        runCommand('rm -f .tmp_setSet', CMD_showLog)

#surfaceToCell<surface> <outsidePoints> <cut> <inside> <outside> <near> <curvature>
def selectProximity(opts, stlFile, distance, BB=[0], outsidePoints=None):
    stlFile = './constant/triSurface/'+stlFile
    if outsidePoints==None:
        tmp = findBoundingBox(stlFile, False)
        tmp = [0.5*(tmp[0]+tmp[3]), 0.5*(tmp[1]+tmp[4])+math.fabs(tmp[1]-tmp[4]), 0.5*(tmp[2]+tmp[5])]
        outsidePoints = " (("+str(tmp[0])+" "+str(tmp[1])+" "+str(tmp[2])+"))"
    else:
        outsidePoints = " (("+str(outsidePoints[0])+" "+str(outsidePoints[1])+" "+str(outsidePoints[2])+"))"
    includeCutCells=' yes'
    includeInside=' yes'
    includeOutside=' no'
    curvature=' -1e6'
    distance = " " + str(distance)
    cmd = 'cellSet c0 '+opts+' surfaceToCell \\"'+ stlFile +'\\"'+ outsidePoints + includeCutCells + includeInside + includeOutside + distance + curvature
    if len(BB)==6:
        BBtxt = '('+str(BB[0])+' '+str(BB[1])+' '+str(BB[2])+') ('+str(BB[3])+' '+str(BB[4])+' '+str(BB[5])+')'
        cmd += "\n" + 'cellSet c0 subset boxToCell '+BBtxt
    print("selectProximity:",cmd)
    if EXEC_REFINEPROXIMITY:
        runCommand('echo "'+cmd+'" > .tmp_setSet', CMD_showLog)
        run_setSet()
        runCommand('rm -f .tmp_setSet', CMD_showLog)
    pass

def refineBox(BB, direction):
    selectBoxToCell(BB)
    if EXEC_REFINEBOX:
        print("refineBox:", direction)
        dirString = ''
        for c in direction:
            if c=='x': dirString += ' tan1'
            if c=='y': dirString += ' tan2'
            if c=='z': dirString += ' normal'

        refineMeshDict = RefineMeshDict(case           = '.',
                                        set            = 'c0',
                                        directions     = dirString,
                                        useHexTopology = True,
                                        geometricCut  = False)
        refineMeshDict.writeFile()
        run_refineMesh()
    else:
        print("refineBox:", direction," ... skip")
    pass

def refineProximity(direction):
    if EXEC_REFINEPROXIMITY:
        print("refineProximity:", direction)
        dirString = ''
        for c in direction:
            if c=='x': dirString += ' tan1'
            if c=='y': dirString += ' tan2'
            if c=='z': dirString += ' normal'

        refineMeshDict = RefineMeshDict(case           = '.',
                                        set            = 'c0',
                                        directions     = dirString,
                                        useHexTopology = True,
                                        geometricCut   = False)
        refineMeshDict.writeFile()
        run_refineMesh()
    else:
        print("refineProximity:", direction," ... skip")
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

def clearEmptyZonesFiles(dryrun=False):
    if not dryrun:
        print("Clear empty {point,cell,face}Zones ... ")
    else:
        print("Clear empty {point,cell,face}Zones ... dryrun")
    files = ['polyMesh/faceZones','polyMesh/cellZones','polyMesh/pointZones']
    folders = getFoamTimeFolders(constant=True)
    for pwd in folders:    
        for fid in files:
            filename = pwd + fid
            if os.path.isfile(filename+'.gz'):
                filename += '.gz'
                cmd='zcat '+filename
            else:
                cmd='cat '+filename
            # remove files if empty
            if os.path.isfile(filename):
                p = subprocess.Popen(cmd+" | sed -n '18,19{/^0$/q1}' ", stdout=subprocess.PIPE, shell=True)
                p.communicate()
                if p.returncode==1: # found
                    print("delete:",filename)
                    if not dryrun:
                        try:
                            os.remove(filename)
                        except OSError as e:
                            print(("Error: %s - %s." % (e.filename,e.strerror)))
    pass

def foamCase_template():
    global NPROCS
    if not os.path.exists('0'): os.makedirs('0')
    if not os.path.exists('constant'): os.makedirs('constant')
    if not os.path.exists('system'): os.makedirs('system')
    
    #controlDict
    controlDict = ControlDict(case              = '.',
                              version           = "foamStar",
                              endTime           = 1000,
                              deltaT            = 0.01,
                              writeControl      = "timeStep",
                              writeInterval     = 50,
                              writePrecision    = 15,
                              writeCompression  = "compressed",
                              runTimeModifiable = "true")
    controlDict.writeFile()
    
    #fvSchemes
    fvSchemes = FvSchemes(case        = '.',
                          simType     = "CrankNicolson",
                          limitedGrad = True,
                          orthogonalCorrection = "implicit")
    fvSchemes.writeFile()
    
    #fvSolution
    fvSolution = FvSolution(case = '.' )
    fvSolution.writeFile()

    #decomposeParDict
    decomposeParDict = DecomposeParDict(case   = '.',
                                        nProcs = NPROCS)
    decomposeParDict.writeFile()

#*** These are templates files ************************************************* 

defaultParams_contents='''
# this [fsMesher] header must exist to identify parameters for fsMesher.py
[fsMesher]

# Only stl in ASCII format. File(s) "stlFile#" will be merged into "constant/triSurface/ship.stl"
# Supports multiple stlFile(s): stlFile0, stlFile1, stlFile2, ..., etc.
# Patch names will be "ship_<stl solid name>"
stlFile0 = ./kcs-hull.stl
stlFile1 = ./kcs-deck.stl

# extra refinement in bow/stern areas with strong curvatures
# refineSurf = ./kcs-curves.stl

# draft (measured from keel)
# if defined, ship.stl will be moved into position
# if "None", no operation will be performed to ship.stl
draft = None

# 180 deg is head sea
heading = 180

# mesh both side or not? choices are "port", "starboard", "both"
side = port

# overall side of the domain (relative to LOA)
# 6 floating numbers: X(min,max) Y(min,max) Z(min,max)
# auto: -3.0, 2.5, -2.0, 2.0, -1.5, 0.5555
domain = auto 

# characteristic length
# auto: compute from "./constant/triSurface/ship.stl"
LOA = auto

# free surface zone
# 2 floating numbers: fsZmin, fsZmax
# auto: 0.02*LOA, 0.01*LOA
fsZone = auto

# cell height within the free surface zone
# auto: fsZmax/3
fsCellHeight = auto

# cell ratio length/height within the free surface zone
# auto: 4
fsCellRatio = auto

# Type of refinement box(es) in the domain
# available type: "wave", "kelvin", "both"
refBoxType = wave

# data for refinement boxes
# format: number of boxes, Xmin,Xmax,Ymin,Ymax, Xmin,Xmax,Ymin,Ymax, ...
# The first number is the number of boxes.
# Dimension of each box is given as 4 floating numbers.
# If only the number of boxes is given, each box will be automatically dimensioned
# according to "refBoxRatio"
# auto: 3
refBoxData = auto

# ratio smallest/largest refinement boxes
# (not use when "refBoxData" is valid)
# auto: 3
refBoxRatio = auto

# cell buffer between refinement level
# auto: 4
cellBuffer = auto

# Additional refinement at the bow region
# auto: 0.2*LOA
refineBow = auto

# Additional refinement at the stern region
# auto: 0.2*LOA
refineStern = auto

# Additional vert. refinement in the free surface zone near the ship
# Accept only "True", "False"
refineFS = True

# Dimension for boundary layers (relative to cell size)
# format: nLayers, layerGrowth, finalLayerThickness, minThicknessRatio
# auto: 3, 1.3, 0.7, 0.7
layers = auto

# by default we add layers to all ship patches "ship_<stl solid name>"
# use "disableLayers" option to disable layer cell on specific patch e.g. deck
# (comma separated list of <stl solid name>)
disableLayers = deck

# these are optional settings for fsMesher.py
[fsMesher-control]

# these switches can be called from command line
# settings called from command line overrule settings in this config-file
# note: NPROCS has no effect when "system/decomposeParDict" exists
DEBUG = True
NPROCS = 4
EXEC_BLOCKMESH = True
EXEC_REFINEBOX = True
EXEC_REFINEPROXIMITY = True
EXEC_SNAP = True
EXEC_ADDLAYERS = True

# main stl-filename in "./constant/triSurface/"
DEFAULT_SHIP_STL = ship.stl

# Self-explained settings ... 
# log-file for foam tools is defined in CMD_keepLog.
# output from python script is flushed to stdOut
CMD_keepLog = ' >> ./log.fsMesher 2>&1 '
CMD_showLog = 'echo \"\\nlog-file: ./log.fsMesher\"; tail -45 ./log.fsMesher; echo \"Please see log-file: ./log.fsMesher\"'
CMD_blockMesh = 'blockMesh'
CMD_autoPatch = 'autoPatch -overwrite 80'
CMD_setSet = 'setSet -latestTime'
CMD_refineMesh = 'refineMesh'
CMD_snappyHexMesh = 'snappyHexMesh'
CMD_decomposePar = 'decomposePar -force -latestTime'

EOF
'''

def computeProximityData(data):
    print("debug: computeProximityData()")

    # cell data
    cSize = data.fsdZ
    cBuffer = data.cellBuffer
    nlevel = int(math.log(math.fabs(float(data.fsCellRatio)), float(2)))

    shipBBxMin = data.shipBBRot[0]
    shipBByMin = data.shipBBRot[1]
    shipBBzMin = data.shipBBRot[2]
    shipBBxMax = data.shipBBRot[3]
    shipBByMax = data.shipBBRot[4]
    shipBBzMax = data.shipBBRot[5]
    outsidePoints = [0.5*(shipBBxMin+shipBBxMax), 0.5*(shipBByMin+shipBByMax)+math.fabs(shipBByMin-shipBByMax), 0.5*(shipBBzMin+shipBBzMax)]
    
#    proxData = []
    n = nlevel
    while n>=1:
        n -= 1
        rData = dict(dataRefineProximity)
        rData['direction'] = 'xy'
        rData['stlFile'] = str(DEFAULT_SHIP_STL)
        rData['outsidePoints'] = list(outsidePoints)

    print(cSize,cBuffer, nlevel)
    pprint.pprint(rData)
    pass

#*** Main execution start here *************************************************
if __name__ == "__main__":
    startTime = time.time()    
    dat = cmdOptions(sys.argv)
    foamCase_template()
    createBlockMeshDict(dat)
    createBackGroundMesh(dat)
    createSnappyMesh(dat)
    
    #print(dat.__dict__)
    #createBoxStl(dat.shipBB, 'boundingBox.stl') # create boundingBox.stl
    #rotateStl('boundingBox.stl', dat.heading, 'boundingBoxRotated.stl')

    endTime = time.time()
    print('Completed meshing in %d minutes' % ((endTime-startTime)/60))
    