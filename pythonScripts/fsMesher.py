#!/usr/bin/python -u

#########################################################################
# Filename: fsMesher.py                                                 #
# Date:     2016-May-17                                                 #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:  An adaptation of CRS Parametric mesh generation tool of     #
#           S. Tornros (Caterpillar Propulsion, tornros_simon@cat.com)  #
#########################################################################

import os, subprocess, time, math
import numpy as np
import sys, argparse, configparser

# sose: 2016-june-03
# FIXME: There are still many conditions inwhich this script will fail
# FIXME: add a sanity check for input parameter(s)
# FIXME: add parallel suuport for refineMesh
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
        print "Found patches in stl: ", shipPatches
        shipBB = findBoundingBox(filename)
        
        if opts['draft'] is not None:
            draft = math.fabs(float(opts['draft']))
            move = -shipBB[2] - draft
            translateStl(filename, [0.0,0.0,move], filename)
            shipBB[2] += move
            shipBB[5] += move

        if opts['refSurfExtra'] is not None:
            nameOnly = os.path.basename(opts['refSurfExtra'])
            surfFile = "./constant/triSurface/"+nameOnly
            print "Create stl: "+surfFile
            runCommand("cp -f "+opts['refSurfExtra']+" "+surfFile)
            if opts['draft'] is not None:
                translateStl(surfFile, [0.0,0.0,move], surfFile)
        
        if rotateStl(filename, opts['heading'], filename):
            shipBBRot = findBoundingBox(filename)
        else:
            shipBBRot = shipBB

        LOA = shipBB[3]-shipBB[0] if opts['LOA']==None else opts['LOA']
        fs = [0.02*LOA, 0.01*LOA] if opts['fs']==None else opts['fs']
        fsdZ = fs[1]/3 if opts['fsdZ']==None else opts['fsdZ']

        if opts['refBow']==None:
            refBow=0
        elif opts['refBow']==True:
            refBow=0.20*LOA

        if opts['refStern']==None:
            refStern=0
        elif opts['refStern']==True:
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
    parser.add_argument('-p','--print-config', dest='showConfig', action='store_true', help='print all available input parameters (including theirs default values) to fsMesher.cfg file and exit. Use this option to generate a template for input files for -f')
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
        print "\nOutput default parameters to file:", DEFAULT_CFG_FILE
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
                print "\nUnknown arg.: --skip", val,"\n"
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
                print "\nUnknown arg.: --run", val,"\n"
                parser.print_help()
                raise SystemExit()

    if args.debug:
        DEBUG=True

    if not DEBUG:
        CMD_setSet += ' -noVTK'
        CMD_refineMesh += ' -overwrite'
        #CMD_snappyHexMesh += ' -overwrite' # snappyHexMesh will be run in parallel!!! FIXME
    else:
        print "Running in DEBUG mode ... "

    subprocess.call('echo "Executing fsMesh.py: $(hostname) @ $(date)" '+CMD_keepLog, shell=True)
    return inputdata
    pass

def readInputParams(filename):
    params = dict(DEFAULT_PARAMS)
    config = configparser.ConfigParser()
    config.read(filename)
    name = 'fsMesher'
    
    print "Read input paramters from file:", filename
    
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
        print "    stlFile(s):", shipStl

    # heading
    try:
        txt = str(config[name]['draft'])
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
            params['refBoxGrad'] = [float(val) for val in txt.split(",")]
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

def rotateLogFile():
    if os.path.isfile(logFile):
        ext = os.path.splitext(logFile)[-1][1:]
        try:
            i = int(ext) + 1
            logFile = os.path.splitext(logFile)[0]+"."+str(i)
        except ValueError:
            i = 0
            logFile += "."+str(i)
        return rotateLogFile()
        pass
    else:
        print 'Log file:',logFile
        return logFile

def setValue(filename, variable, value):
    subprocess.call('sed -i "s/^'+variable+'[ \t].*;/'+variable+' '+str(value)+';/g" '+filename, shell=True)

# foamFile may be compressed i.e. ".gz"  
def foamFileExist(filename):
    found = os.path.isfile(filename)
    if not found:
        found = os.path.isfile(filename + '.gz')
    return found

# foamBlock exists or not?
def foamBlockExist(filename, blockName):
    # check one-line format: "name { ... }"
    cmd="sed -e '/^[ \\t]*[\"]*\<"+blockName+"\>[\"]*[ \\t]*{.\+}/"
    p = subprocess.Popen(cmd+"{q1}' "+filename, stdout=subprocess.PIPE, shell=True)
    p.communicate()
    if p.returncode==0: # not found
        # check multi-lines format: "name { \newline ... }"
        cmd="sed -e '/^[ \\t]*[\"]*\<"+blockName+"\>[\"]*[ \\t]*{/,/}/"        
        p = subprocess.Popen(cmd+"{q1}' "+filename, stdout=subprocess.PIPE, shell=True)
        p.communicate()
    if p.returncode==0: # still not found
        # check multi-lines format: "name \newline { .. }"
        cmd="sed -e '/^[ \\t]*[\"]*\<"+blockName+"\>[\"]*/,/}/"
        p = subprocess.Popen(cmd+"{q1}' "+filename, stdout=subprocess.PIPE, shell=True)
        p.communicate()
    if p.returncode==0: # still not found
        return False,cmd
    return True,cmd

def renameFoamBlock(filename, oldName, newName):
    cmd="sed -i -e '1h;2,$H;$!d;g' -e 's/\<"+oldName+"\>\([ \\n\\t]\+{\)/"+newName+"\\1/' "
    subprocess.call(cmd+filename, shell=True)
    pass

def modifyFoamBlock(filename, blockName, variable, val):
    OK,cmd = foamBlockExist(filename, blockName)
    if OK:
        cmd += "{s/\(\<"+variable+"\>[ \\t]*\)[^;]*/\\1"+str(val)+"/}' -i "
        subprocess.call(cmd+filename, shell=True)
    else: # not found
        print "modifyFoamBlock: Substitution failed"
        print "filename: ", filename
        print "blockName: ", blockName
        print "keyword: ", variable
        print "value: ", val
        raise SystemExit('abort ...')

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
    
def createShipStl(shipStl, outputStl, overwrite=False):
    print "Creating stl: "+outputStl
    if os.path.isfile(outputStl) & (not overwrite):
        print "    file already exist ... reuse existing file"
        return

    subprocess.call("rm -fr "+outputStl, shell=True)
    if len(shipStl)==0:
        print "    stlFile not specified ..."
        raise SystemExit('abort ...')
    foo = 0
    while foo < len(shipStl):
        stlFile = shipStl[foo]
        foo += 1
        #if not stlFile[::-1][0:4][::-1]==".stl":
        if not stlFile.endswith('.stl'):
            stlFile = stlFile + '.stl'            
        print "    append: " + stlFile
        if not os.path.isfile(stlFile):
            print "\nFile not found:",stlFile,"\n"
            raise SystemExit('abort ...')
        subprocess.call("mkdir -p $(dirname "+outputStl+")", shell=True)
        runCommand("cat "+ stlFile + " >> "+outputStl)

def findSTLPatches(stlFile):
    p = subprocess.Popen("grep '^[ \\t]*\<solid\>' "+stlFile+" | sed 's/solid//g' | tr '\n' ' ' | sed 's/^[ \t]*//;s/ \+/ /g;s/\s*$//g' "  , stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    patches,error = p.communicate()
    if error:
        print 'error: ', error
        raise SystemExit('abort ...')
    return patches.split(' ')
    
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

def translateStl(inputStl, val, outputStl):
    val = "("+str(val[0])+" "+str(val[1])+" "+str(val[2])+")"
    print "Translate stl by "+val+": " + inputStl
    subprocess.call("surfaceTransformPoints -translate '"+val+"' " + inputStl + " " + outputStl + " > /dev/null", shell=True)
    return True

def rotateStl(inputStl, heading, outputStl):
    yaw = heading-180.0
    if str(yaw) == '0.0': return False
    print "Rotate stl (0 0 " + str(yaw) + "): " + inputStl
    subprocess.call("surfaceTransformPoints -rollPitchYaw '(0 0 "+str(heading-180.0)+")' " + inputStl + " " + outputStl + " > /dev/null", shell=True)
    return True

def createBoxStl(BB,name):
    Xmin,Ymin,Zmin,Xmax,Ymax,Zmax=BB[0],BB[1],BB[2],BB[3],BB[4],BB[5]
    tol = (Xmax-Xmin)*1e-6
    filename = "./constant/triSurface/" + name
    print "Creating stl: " + filename
    print "   ",BB
    subprocess.call("surfaceTransformPoints -scale '("+str(Xmax-Xmin-tol)+" "+str(Ymax-Ymin-tol)+" "+str(Zmax-Zmin-tol)+")' fsMesher/fsMesher_box.stl "+filename+" > /dev/null", shell=True)
    subprocess.call("surfaceTransformPoints -translate '("+str(Xmin+0.5*tol)+" "+str(Ymin+0.5*tol)+" "+str(Zmin+0.5*tol)+")' "+filename+" "+filename+" > /dev/null", shell=True)

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

def createBlockMeshDict(data, fileName):
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
        print "\nUnknown parameters, side=",data.side
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
    refBoxZdata = [zGrid[0] for i in range(nRefBox), zGrid[-1] for i in range(nRefBox)]
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
    
    filename = 'constant/polyMesh/blockMeshDict'
    blockMeshDict_template(filename)
    setValue(filename, 'XminDomain', XminDomain)
    setValue(filename,'XmaxDomain', XmaxDomain)
    setValue(filename,'YminDomain', YminDomain)
    setValue(filename,'YmaxDomain', YmaxDomain)
    setValue(filename,'ZminDomain', ZminDomain)
    setValue(filename,'Xcells', Xcells)
    setValue(filename,'Ycells', Ycells)

    print 'Domain bounding box:'
    print "   ", [XminDomain, YminDomain, ZminDomain, XmaxDomain, YmaxDomain, ZmaxDomain]

    # add z-cuts points
    # get line number at ($XminDomain $YmaxDomain $ZminDomain)
    p = subprocess.Popen("grep -n '//VERTICES' "+filename, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    line,error = p.communicate()
    line = line.split(":")
    line = int(line[0])+1
    for val in zAllCut:
        subprocess.call('sed -i "'+str(line)+'i\    (\$XminDomain \$YminDomain '+str(val)+')" '+filename, shell=True)
        line += 1
        subprocess.call('sed -i "'+str(line)+'i\    (\$XmaxDomain \$YminDomain '+str(val)+')" '+filename, shell=True)
        line += 1
        subprocess.call('sed -i "'+str(line)+'i\    (\$XmaxDomain \$YmaxDomain '+str(val)+')" '+filename, shell=True)
        line += 1
        subprocess.call('sed -i "'+str(line)+'i\    (\$XminDomain \$YmaxDomain '+str(val)+')" '+filename, shell=True)
        line += 1

    # add hex block
    # get line number at "hex (...)"
    p = subprocess.Popen("grep -n '//HEX_BLOCK' "+filename, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    line,error = p.communicate()
    line = line.split(":")
    line = int(line[0])+1
    idx = 0
    for i,val in enumerate(zAllCutNCells):
        subprocess.call('sed -i "'+str(line)+'i\    hex ('\
        +str(idx)+' '+str(idx+1)+' '+str(idx+2)+' '+str(idx+3)+' '+str(idx+4)+' '+str(idx+5)+' '+str(idx+6)+' '+str(idx+7)+\
        ') (\$Xcells \$Ycells ' +str(val)+ ') simpleGrading (1 1 '+str(zAllCutRatio[i])+')" '+filename, shell=True)
        idx += 4
        line += 1

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
            print "not yet implemented"
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
        print "blockMesh: create a base mesh ..."
        runCommand(CMD_blockMesh + CMD_keepLog)
        print "autoPatch: create domain boundaries ..."
        runCommand(CMD_autoPatch + CMD_keepLog)
        renameFoamBlock('./constant/polyMesh/boundary', 'auto0', 'domainX0')
        renameFoamBlock('./constant/polyMesh/boundary', 'auto1', 'domainX1')
        renameFoamBlock('./constant/polyMesh/boundary', 'auto2', 'domainY0')
        renameFoamBlock('./constant/polyMesh/boundary', 'auto3', 'domainY1')
        renameFoamBlock('./constant/polyMesh/boundary', 'auto4', 'domainZ0')
        renameFoamBlock('./constant/polyMesh/boundary', 'auto5', 'domainZ1')
        modifyFoamBlock('./constant/polyMesh/boundary', 'domainZ0', 'type', 'wall')
        if data.side=='port':
            modifyFoamBlock('./constant/polyMesh/boundary', 'domainY0', 'type', 'symmetryPlane')
        elif data.side=='starboard':
            modifyFoamBlock('./constant/polyMesh/boundary', 'domainY1', 'type', 'symmetryPlane')
#        if NPROCS>1:
#            print "decomposePar: nProcs =",NPROCS
#            runCommand(CMD_decomposePar + CMD_keepLog)
    else:
        print "blockMesh: ... skip"

    # how many refinement box? minimum is 1
    refBoxData = list(data.refBoxData)
    nRefBox = int(refBoxData[0])
    del refBoxData[0]
    refBoxBB = []
    if len(data.refBoxData)>1:
        if (nRefBox != (len(data.refBoxData))/4):
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

    # refine free surface (using proximity method)
    distance = data.fsdZ*data.cellBuffer*2.0*(nRefBox-1.0 + float(bool(data.refBow) | bool(data.refStern) | bool(data.refFS)))    
    
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

    #print "debug: refBoxData: ", data.refBoxData
    for BB in refBoxBB:
        refineBox(BB, 'xy')
        pass    
    
    # this is a point outside ship.stl
    outsidePoints = [0.5*(shipBBxMin+shipBBxMax), 0.5*(shipBByMin+shipBByMax)+math.fabs(shipBByMin-shipBByMax), 0.5*(shipBBzMin+shipBBzMax)]

    BB = [-1e6,-1e6,data.fsZmin,1e6,1e6,data.fsZmax]
    selectProximity('new', DEFAULT_SHIP_STL, distance, outsidePoints=outsidePoints)
    refineProximity('xy')
    #
    distance *= 0.5
    selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
    refineProximity('xy')
    BB = [-1e6,-1e6,-1e6,1e6,1e6,data.fsZmin]
    selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
    refineProximity('xyz')
    BB = [-1e6,-1e6,data.fsZmax,1e6,1e6,1e6]
    selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
    refineProximity('xyz')
    #
    distance *= 0.5
    if bool(data.refBow):
        BB = [shipBBxMax-data.refBow,-1e6,-1e6,1e6,1e6,shipBBzMax-data.fsdZ*data.cellBuffer]
        selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('xyz')
        if data.refSurfExtra is not None:
            nameOnly = os.path.basename(data.refSurfExtra)
            BB = [shipBBxMax-data.refBow+0.5*data.fsdZ*data.cellBuffer,-1e6,-1e6,1e6,1e6,shipBBzMax-1.5*data.fsdZ*data.cellBuffer]
            selectProximity('new', nameOnly, 0.5*distance, BB=BB, outsidePoints=outsidePoints)
            refineProximity('xyz')
        
    if bool(data.refStern):
        BB = [-1e6,-1e6,-1e6,shipBBxMin+data.refStern,1e6,shipBBzMax-data.fsdZ*data.cellBuffer]
        selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('xyz')
        if data.refSurfExtra is not None:
            nameOnly = os.path.basename(data.refSurfExtra)
            BB = [-1e6,-1e6,-1e6,shipBBxMin+data.refStern-0.5*data.fsdZ*data.cellBuffer,1e6,shipBBzMax-1.5*data.fsdZ*data.cellBuffer]
            selectProximity('new', nameOnly, 0.5*distance, BB=BB, outsidePoints=outsidePoints)
            refineProximity('xyz')

    if bool(data.refFS):
        BB = [shipBBxMin+data.refStern,-1e6,data.fsZmin,shipBBxMax-data.refBow,1e6,data.fsZmax]
        selectProximity('new', DEFAULT_SHIP_STL, distance, BB=BB, outsidePoints=outsidePoints)
        refineProximity('z')

def createSnappyMesh(data):
    filename='./system/snappyHexMeshDict'
    snappyMesh_template(filename, data.shipPatches, noLayers=data.noLayers)
    setValue(filename,"locationInMeshX",data.locationInMesh[0])
    setValue(filename,"locationInMeshY",data.locationInMesh[1])
    setValue(filename,"locationInMeshZ",data.locationInMesh[2])
    if EXEC_SNAP:
        setValue(filename, "castellatedMesh", "true")
        setValue(filename, "snap", "true")
        setValue(filename, "addLayers", "false")
        if NPROCS>1:
            isDecomposed,nProcs = caseAlreadyDecomposed()
            if not isDecomposed:
                nProcs = NPROCS
                print "decomposePar: nProcs =",nProcs
                runCommand(CMD_decomposePar + CMD_keepLog)
            print "snappyHexMesh: snapping ... in parallel, nProcs =",nProcs
            cmd = "mpirun -np "+str(nProcs)+" "+CMD_snappyHexMesh + " -parallel "+ CMD_keepLog
        else:
            print "snappyHexMesh: snapping ... "
            cmd = CMD_snappyHexMesh + CMD_keepLog
        runCommand(cmd)
    else:
        print "snappyHexMesh: snap ... skip"
        pass
        
    if EXEC_ADDLAYERS:
        nLayers = data.shipBL[0]
        layerGrowth = data.shipBL[1]
        finalLayerThickness = data.shipBL[2]
        minThicknessRatio = data.shipBL[3]
        minThickness = float(minThicknessRatio)*finalLayerThickness/(pow(layerGrowth,float(nLayers-1)))
        setValue(filename, "castellatedMesh", "false")
        setValue(filename, "snap", "false")
        setValue(filename, "addLayers", "true")
        setValue(filename, "SHIP_BL_layers", nLayers)
        setValue(filename, "SHIP_BL_layerGrowth", layerGrowth)
        setValue(filename, "SHIP_BL_finalLayerThickness", finalLayerThickness)
        setValue(filename, "SHIP_BL_minThickness", minThickness)
        if NPROCS>1:
            isDecomposed,nProcs = caseAlreadyDecomposed()
            if not isDecomposed:
                nProcs = NPROCS
                print "decomposePar: nProcs =",nProcs
                runCommand(CMD_decomposePar + CMD_keepLog)
            print "snappyHexMesh: add layers ... in parallel, nProcs =",nProcs
            cmd = "mpirun -np "+str(nProcs)+" "+CMD_snappyHexMesh + " -parallel "+ CMD_keepLog
        else:
            print "snappyHexMesh: add layers ... "
            cmd = CMD_snappyHexMesh + CMD_keepLog
        runCommand(cmd)
    else:
        print "snappyHexMesh: add layers ... skip"
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

def run_refineMesh():
    # FIXME: refineMesh is buggy when run in parallel
    runCommand(CMD_refineMesh + CMD_keepLog)
    
    #if NPROCS==1:
    #    runCommand(CMD_refineMesh + CMD_keepLog)
    #else:
    #    runCommand('mpirun -np '+str(NPROCS)+' '+CMD_refineMesh + ' -parallel ' + CMD_keepLog)
    pass

def selectBoxToCell(BB):
    xmin,ymin,zmin,xmax,ymax,zmax=BB[0],BB[1],BB[2],BB[3],BB[4],BB[5]
    BBtxt = '('+str(xmin)+' '+str(ymin)+' '+str(zmin)+') ('+str(xmax)+' '+str(ymax)+' '+str(zmax)+')'
    cmd = 'cellSet c0 new boxToCell '+BBtxt
    print "SelectBoxToCell:",cmd
    if EXEC_REFINEBOX:
        runCommand('echo "'+cmd+'" > .tmp_setSet')
        run_setSet()
        runCommand('rm -f .tmp_setSet')

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
    print "selectProximity:",cmd
    if EXEC_REFINEPROXIMITY:
        runCommand('echo "'+cmd+'" > .tmp_setSet')
        run_setSet()
        runCommand('rm -f .tmp_setSet')
    pass

def refineBox(BB, direction):
    selectBoxToCell(BB)
    if EXEC_REFINEBOX:
        print "refineBox:", direction
        refineMesh_template('./system/refineMeshDict', direction)
        run_refineMesh()
    else:
        print "refineBox:", direction," ... skip"
    pass

def refineProximity(direction):
    if EXEC_REFINEPROXIMITY:
        print "refineProximity:", direction
        refineMesh_template('./system/refineMeshDict', direction)
        run_refineMesh()
    else:
        print "refineProximity:", direction," ... skip"
    pass

def getFoamTimeFolders(constant=False):
    found = []
    for val in os.listdir('.'):
        try:
            float(val)
            found.append(val)
        except ValueError:
            pass
    found.sort(key=float)
    timeFolders = ['./constant/'] if constant else []
    if len(found):
        for val in found:
            timeFolders.append('./'+val+'/')
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

def clearEmptyZonesFiles(dryrun=False):
    if not dryrun:
        print "Clear empty {point,cell,face}Zones ... "
    else:
        print "Clear empty {point,cell,face}Zones ... dryrun"
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
                    print "delete:",filename
                    if not dryrun:
                        try:
                            os.remove(filename)
                        except OSError, e:
                            print ("Error: %s - %s." % (e.filename,e.strerror))
    pass


def foamCase_template():
    global NPROCS
    subprocess.call('mkdir -p 0 constant system', shell=True)
    if not os.path.isfile('system/controlDict'):
        subprocess.call('cat << EOF > system/controlDict' + controlDict_contents, shell=True)
    if not os.path.isfile('system/fvSchemes'):
        subprocess.call('cat << EOF > system/fvSchemes' + fvSchemes_contents, shell=True)
    if not os.path.isfile('system/fvSolution'):
        subprocess.call('cat << EOF > system/fvSolution' + fvSolution_contents, shell=True)
    if not os.path.isfile('system/decomposeParDict'):
        subprocess.call('cat << EOF > system/decomposeParDict' + decomposeParDict_contents, shell=True)
        if NPROCS>1:
            setValue('system/decomposeParDict', 'numberOfSubdomains', NPROCS)
    else:
        p = subprocess.Popen("grep -E '^numberOfSubdomains' system/decomposeParDict | sed 's/^numberOfSubdomains[ \\t]*//;s/;.*//g' ", stdout=subprocess.PIPE, shell=True)
        txt,error = p.communicate()
        if error:
            print "error:",error
            raise SystemExit('abort ...')
        else:
            NPROCS = int(txt)
    pass

def blockMeshDict_template(filename):
    subprocess.call('mkdir -p $(dirname '+filename+")", shell=True)
    subprocess.call('cat << EOF > '+filename + blockMeshDict_contents, shell=True)
    pass

def refineMesh_template(filename, dirString):
    subprocess.call('mkdir -p $(dirname '+filename+")", shell=True)
    subprocess.call('cat << EOF > '+filename + refineMeshDict_contents, shell=True)
    directions = ''
    for c in dirString:
        if c=='x':
            directions += ' tan1'
        if c=='y':
            directions += ' tan2'
        if c=='z':
            directions += ' normal'
    subprocess.call('sed -i "s/^directions.*;/directions ('+directions+' );/" '+filename, shell=True)
    pass
    
def snappyMesh_template(filename, shipPatches, noLayers=None):
    # don't overwrite existing file
    if not os.path.isfile(filename):
        print "Creating:",filename
        subprocess.call('mkdir -p $(dirname '+filename+")", shell=True)
        subprocess.call('cat << EOF > '+filename + snappyHexMeshDict_contents, shell=True)
        renameFoamBlock(filename, "DEFAULT_SHIP_STL", '\\"'+DEFAULT_SHIP_STL+'\\"')
        # add ship_patches to snappyHexMeshDict
        # get line number at "//ADD_SHIP_PATCHES_HERE"
        p = subprocess.Popen("grep -n '//ADD_SHIP_PATCHES_HERE' "+filename, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        line,error = p.communicate()
        line = line.split(":")
        line = int(line[0])+1
        if (len(shipPatches)>1):
            if noLayers==None: noLayers = []
            if len(noLayers)>0: print "    disable layers on patch(es):",noLayers
            for val in shipPatches:
                if val in noLayers: continue
                subprocess.call('sed -i "'+str(line)+'i\        ship_'+str(val)+'  { nSurfaceLayers \$SHIP_BL_layers; }" '+filename, shell=True)
                line += 1
        else:
            subprocess.call('sed -i "'+str(line)+'i\        ship { nSurfaceLayers \$SHIP_BL_layers; }" '+filename, shell=True)
            line += 1

    # "snappyHexMeshDict" already exists, check whether or not the "geometry" is OK
    else:
        print "Reuse existing file:",filename
        OK,cmd = foamBlockExist(filename, DEFAULT_SHIP_STL)
        if not OK:
            print "\n    Warning: geomtry definition not found for:",DEFAULT_SHIP_STL
            print "\n    Please review:",filename
            print "\n"
            
    eMeshFile=os.path.splitext(DEFAULT_SHIP_STL)[0] if DEFAULT_SHIP_STL.endswith('.stl') else DEFAULT_SHIP_STL
    if not foamFileExist('./constant/triSurface/'+eMeshFile+'.eMesh'):
        print "\nExtract surface features from file: ./constant/triSurface/"+DEFAULT_SHIP_STL
        subprocess.call('cat << EOF > ./system/surfaceFeatureExtractDict' + surfaceFeatureExtractDict_contents, shell=True)
        subprocess.call('sed -i "s/^ship.stl/'+DEFAULT_SHIP_STL+'/" '+filename, shell=True)
        runCommand(CMD_surfaceFeatureExtract + CMD_keepLog)
    else:
        print "Reuse existing file: "+'./constant/triSurface/'+eMeshFile+'.eMesh'
    
    pass

#*** These are templates files ************************************************* 

blockMeshDict_contents = '''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by fsMesher on $(hostname) @ $(date)
*/

fastMerge yes;
convertToMeters 1;

XminDomain	moin;
XmaxDomain	moin;
YminDomain	moin;
YmaxDomain	moin;
ZminDomain  moin;
Xcells	moin;
Ycells	moin;

vertices
(
    (\$XminDomain \$YminDomain \$ZminDomain)
    (\$XmaxDomain \$YminDomain \$ZminDomain)
    (\$XmaxDomain \$YmaxDomain \$ZminDomain)
    (\$XminDomain \$YmaxDomain \$ZminDomain)
    //VERTICES
);

blocks
(
    //HEX_BLOCK
);

edges ();

boundary
(
    defaultFaces
    {
        type patch; faces ();
    }
);

mergePatchPairs ();

// ************************************************************************* //
EOF
'''

refineMeshDict_contents = '''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      refineMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by fsMesher on $(hostname) @ $(date)
*/
set             c0;
directions      ( tan1 tan2 normal );
coordinateSystem global;
globalCoeffs { tan1 ( 1 0 0 ); tan2 ( 0 1 0 ); }
patchLocalCoeffs { patch outside; tan1 ( 1 0 0 ); }
useHexTopology  yes;
geometricCut    no;
writeMesh       no;

// ************************************************************************* //
EOF
'''

surfaceFeatureExtractDict_contents = '''
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  3.0.x                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      surfaceFeatureExtractDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by fsMesher on $(hostname) @ $(date)
*/

ship.stl
{
    extractionMethod    extractFromSurface;
    writeObj            true;    
    extractFromSurfaceCoeffs
    {
        // Mark edges whose adjacent surface normals are at an angle less
        // than includedAngle as features
        // - 0  : selects no edges
        // - 180: selects all edges
        includedAngle   150;
    }
}

// ************************************************************************* //
EOF
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
*   created by fsMesher on $(hostname) @ $(date)
*/

application     foamStar;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         1000;

deltaT          0.01;

writeControl    timeStep;

writeInterval   50;

purgeWrite      0;

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
    //motionInfo { type motionInfo; }
    //#include "forces.inc"
    //#include "vbm.inc"
    //#include "waveProbe.inc"
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
*   created by fsMesher on $(hostname) @ $(date)
*/

ddtSchemes
{
    default      CrankNicolson 0.9;
}

gradSchemes
{
    default         cellLimited leastSquares 1;
    limitedGrad     cellLimited Gauss linear 1;
}

divSchemes
{
    div(rhoPhi,U)   Gauss linearUpwind grad(U);
    div(phi,alpha)  Gauss vanLeer;
    div(phirb,alpha) Gauss linear;
    div(phi,k)      Gauss linearUpwind limitedGrad;
    div(phi,omega)  Gauss linearUpwind limitedGrad;
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
    p_rgh; pcorr; alpha.water;
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
*   created by fsMesher on $(hostname) @ $(date)
*/

solvers
{
    "alpha.water.*"
    {
        nAlphaCorr      3;
        nAlphaSubCycles 1;
        cAlpha          0.5;
        icAlpha         0;

        MULESCorr       yes;
        nLimiterIter    10;
        alphaApplyPrevCorr  no;

        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-10;
        relTol          0;
        minIter         2;
    }

    "pcorr.*"
    {
        solver          GAMG;
        smoother        DIC;
        agglomerator    faceAreaPair;
        mergeLevels     1;
        nCellsInCoarsestLevel 10;
        cacheAgglomeration true;

        tolerance       1e-4;
        relTol          0;
    };

    p_rgh
    {
        solver          GAMG;
        smoother        DIC;
        agglomerator    faceAreaPair;
        mergeLevels     1;
        nCellsInCoarsestLevel 10;
        cacheAgglomeration true;

        tolerance       1e-8;
        relTol          1e-4;
    };

    p_rghFinal
    {
        \$p_rgh;
        relTol          0;
    }

    "(U|k|omega).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        nSweeps         1;

        tolerance       1e-8;
        relTol          0;
        minIter         1;
    };

    "cellDisplacement.*"
    {
        solver          GAMG;
        tolerance       1e-7;
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
    momentumPredictor   yes;
    nOuterCorrectors    4;
    nCorrectors         4;
    nNonOrthogonalCorrectors 1;
    correctPhi          no;
    moveMeshOuterCorrectors yes;
}

relaxationFactors
{
    fields
    {
        //U 0.7; p_rgh 0.3;
    }
    equations
    {
        //UEqn 0.7; prghEqn 0.3;
    }
}

cache
{
    grad(U);
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
*   created by fsMesher on $(hostname) @ $(date)
*/

numberOfSubdomains 4;
method          scotch;
distributed     no;
roots           ( );

EOF
'''

snappyHexMeshDict_contents='''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*
*   created by fsMesher on $(hostname) @ $(date)
*/

#inputMode overwrite;

// user-defined parameters
locationInMeshX moin;
locationInMeshY moin;
locationInMeshZ moin;

SHIP_BL_layers 3;
SHIP_BL_layerGrowth 1.3;
SHIP_BL_finalLayerThickness 0.7;
SHIP_BL_minThickness 0.289;

SHIP_EDGE_LVL 0;
SHIP_HULL_LVL_MIN 0;
SHIP_HULL_LVL_MAX 0;

castellatedMesh moin;
snap            moin;
addLayers       moin;

geometry
{
    DEFAULT_SHIP_STL { type triSurfaceMesh; name ship; patchInfo { type wall; } }
};

castellatedMeshControls
{
    maxLocalCells 10000000;
    maxGlobalCells 100000000;
    minRefinementCells 0;
    nCellsBetweenLevels 1;

    locationInMesh (\$locationInMeshX \$locationInMeshY \$locationInMeshZ); 
    
    features
    (
        { file "ship.eMesh"; level \$SHIP_EDGE_LVL; }
    );

    refinementSurfaces
    {
        ship { level (\$SHIP_HULL_LVL_MIN \$SHIP_HULL_LVL_MAX); }
    }

    resolveFeatureAngle 15;

    refinementRegions
    {
    }

    allowFreeStandingZoneFaces 	false;
}

snapControls
{
    nSmoothPatch 3;
    tolerance 0.75;
    nSolveIter 100;
    nRelaxIter 5;
    nFeatureSnapIter 10;
    implicitFeatureSnap false;
    explicitFeatureSnap true;
    multiRegionFeatureSnap true;
}

addLayersControls
{
    relativeSizes true;
    layers
    {
        //ADD_SHIP_PATCHES_HERE
    }

    expansionRatio \$SHIP_BL_layerGrowth;
    finalLayerThickness \$SHIP_BL_finalLayerThickness;
    minThickness \$SHIP_BL_minThickness;
    nGrow 0;

    featureAngle 60;
    nRelaxIter 5;
    nSmoothSurfaceNormals 1;
    nSmoothNormals 3;
    nSmoothThickness 10;
    maxFaceThicknessRatio 0.5;
    maxThicknessToMedialRatio 0.3;
    minMedianAxisAngle 90;
    nBufferCellsNoExtrude 0;
    nLayerIter 50;
    nRelaxedIter 20;
    nMedialAxisIter 10;
}
meshQualityControls

{
    maxNonOrtho 65;
    maxBoundarySkewness 20;
    maxInternalSkewness 4;
    maxConcave 80;
    minVol 1e-13;
    minTetQuality 1e-15;
    minArea -1;
    minTwist 0.02;
    minDeterminant 0.001;
    minFaceWeight 0.05;
    minVolRatio 0.01;
    minTriangleTwist -1;
    minVolCollapseRatio 0.1;
    //
    nSmoothScale 4;
    errorReduction 0.75;
    //
    relaxed
    {
        maxNonOrtho 75;
    }
}

mergeTolerance 1E-6;
debug 0;

EOF
'''

defaultParams_contents='''
# this [fsMesher] header must exist to identify parameters for fsMesher.py
[fsMesher]

# Only stl in ASCII format. File(s) "stlFile#" will be merged into "constant/triSurface/ship.stl"
# Supports multiple stlFile(s): stlFile0, stlFile1, stlFile2, ..., etc.
# Patch names will be "ship_<stl solid name>"
stlFile0 = ./kcs-hull.stl
stlFile1 = ./kcs-deck.stl

# extra refinement in bow/stern areas with strong curvatures
refineSurf = ./kcs-curves.stl

# draft (measured from keel)
# if defined, ship.stl will be moved into position
# if "None" or "not defined", no operation will be performed to ship.stl
draft = None

# 180 deg is headsea
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
DEBUG = False
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


#*** Main execution start here *************************************************
if __name__ == "__main__":
    startTime = time.time()    
    dat = cmdOptions(sys.argv)
    foamCase_template()
    createBlockMeshDict(dat, 'constant/polyMesh/blockMeshDict')
    createBackGroundMesh(dat)
    createSnappyMesh(dat)
    
    #print dat.__dict__
    #createBoxStl(dat.shipBB, 'boundingBox.stl') # create boundingBox.stl
    #rotateStl('boundingBox.stl', dat.heading, 'boundingBoxRotated.stl')

    endTime = time.time()
    print 'Completed meshing in %d minutes' % ((endTime-startTime)/60)
    
