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
DEFAULT_SHIP_STL = 'ship.stl'
DEFAULT_CFG_FILE = 'fsMesher.cfg'
DEFAULT_PARAMS = {
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
            print "Set draft:",draft
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
    global CMD_showLog
    # default global parameters
    EXEC_BLOCKMESH = True
    CMD_showLog = 'echo "\nlog-file: ./log.fsMesher"; tail -45 ./log.fsMesher; echo "Please see log-file: ./log.fsMesher"'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='inputfile', help='read parameters from input file. Use this option with a stl-file to generate a mesh with default parameters, e.g.: fsMesher.py -f myship.stl ')
    parser.add_argument('-s','--skip', help='skip selected steps, same arguments as for -r. \
                                            Do not use -s and -r options simultaneously')
    parser.add_argument('-clz','--clear-zones', dest='clearZones', action='store_true', help='Clear empty {point,cell,face}Zones files and exit. With -n option no file-removal action will be made.')
    args = parser.parse_args()

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


    subprocess.call('echo "Executing fsMesh.py: $(hostname) @ $(date)" '+CMD_keepLog, shell=True)
    return inputdata
    pass

def readInputParams(filename):
    params = dict(DEFAULT_PARAMS)
    config = configparser.ConfigParser()
    config.read(filename)
    name = 'template'
    
    print "Read input paramters from file:", filename
    
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
        if not txt=='auto':
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
    #
    global DEBUG
    global NPROCS  
    global EXEC_BLOCKMESH
    global CMD_showLog

    txt = config['fsMesher-control'].getboolean('EXEC_BLOCKMESH')
    if not txt=='None': EXEC_BLOCKMESH = txt

    try:
        txt = config['fsMesher-control']['NPROCS']
        if not txt == None:
            NPROCS = int(txt)
    except KeyError: pass    

    try:
        txt = config['fsMesher-control']['CMD_showLog']
        if (txt[0] == txt[-1]) and txt.startswith(("'", '"')):
             txt = txt[1:-1]
        CMD_showLog = txt
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

    #print "debug: refBoxData: ", data.refBoxData
    for BB in refBoxBB:
        refineBox(BB, 'xy')
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
    
    #controlDict
    filename = 'system/controlDict'
    if not os.path.isfile(filename):
        subprocess.call('cat << EOF > system/controlDict' + controlDict_contents, shell=True)
        setValue(filename, 'START_TIME', opts['startTime'])
        setValue(filename, 'END_TIME', opts['endTime'])
        setValue(filename, 'TIME_STEP', opts['timeStep'])
        setValue(filename, 'WRITE_INTERVAL', opts['writeInterval'])
        setValue(filename, 'PURGE_WRITE', opts['purgeWrite'])
        if not opts['outputMotions']: setValue(filename, 'motionInfo {', '//motionInfo {')
        if not opts['outputVBM']:
            setValue(filename, '#include "vbm.inc"', '//#include "vbm.inc"')
        else:
            if not os.path.isfile('system/vbm.inc'):
                subprocess.call('cat << EOF > system/vbm.inc' + vbm_contents, shell=True)
                setValue('system/vbm.inc', 'HULL_PATCH', opts['hullPatch'])
                setValue('system/vbm.inc', 'DON_FILENAME', opts['donFile'])
        if not opts['outputWave']:
            setValue(filename, '#include "waveProbe.inc"', '//#include "waveProbe.inc"')
        else:
            if not os.path.isfile('system/waveProbe.inc'):
                subprocess.call('cat << EOF > system/waveProbe.inc' + waveProbe_contents, shell=True)

    #fvSchemes
    filename = 'system/fvSchemes'
    if not os.path.isfile('system/fvSchemes'):
        subprocess.call('cat << EOF > system/fvSchemes' + fvSchemes_contents, shell=True)
        setValue(filename, 'DDT_SCHEME', opts['scheme'])

    #fvSolution
    filename = 'system/fvSolution'
    if not os.path.isfile('system/fvSolution'):
        subprocess.call('cat << EOF > system/fvSolution' + fvSolution_contents, shell=True)
        setValue(filename, 'FSI_TOL', opts['fsiTol'])
    
    #decomposeParDict
    filename = 'system/decomposeParDict'
    if not os.path.isfile('system/decomposeParDict'):
        subprocess.call('cat << EOF > system/decomposeParDict' + decomposeParDict_contents, shell=True)
        setValue(filename, 'NB_PROCS', opts['nProcs'])

    #waveProperties
    filename = 'constant/waveProperties'
    if not os.path.isfile('constant/waveProperties'):
        subprocess.call('cat << EOF > constant/waveProperties' + waveProperties_contents, shell=True)
        setValue(filename, 'WAVE_TYPE', opts['wave'])
        setValue(filename, 'WAVE_START_TIME', opts['waveStartTime'])
        setValue(filename, 'WAVE_RAMP_TIME', opts['waveRampTime'])
        setValue(filename, 'WAVE_VELOCITY', opts['velocity'])
        setValue(filename, 'WAVE_DEPTH', opts['depth'])
        setValue(filename, 'WAVE_HEIGHT', opts['waveH'])
        setValue(filename, 'WAVE_PERIOD', opts['waveU'])
        setValue(filename, 'WAVE_SEALEVEL', opts['waveSeaLvl'])
        #relaxation zones must be handeled manually
    
    #turbulenceProperties
    filename = 'constant/turbulenceProperties'
    if not os.path.isfile('constant/turbulenceProperties'):
        subprocess.call('cat << EOF > constant/turbulenceProperties' + turbulenceProperties_contents, shell=True)
    
    #dynamicMeshDict
    filename = 'constant/dynamicMeshDict'
    if not os.path.isfile('constant/dynamicMeshDict'):
        subprocess.call('cat << EOF > constant/dynamicMeshDict' + dynamicMeshDict_contents, shell=True)
        
    #g
    filename = 'constant/g'
    if not os.path.isfile('constant/g'):
        subprocess.call('cat << EOF > constant/g' + g_contents, shell=True)
        
    #RASProperties
    filename = 'constant/RASProperties'
    if not os.path.isfile('constant/RASProperties'):
        subprocess.call('cat << EOF > constant/RASProperties' + RASProperties_contents, shell=True)
        
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
# if "None", no operation will be performed to ship.stl
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

def computeProximityData(data):
    print "debug: computeProximityData()"

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
    
    proxData = []
    n = nlevel
    while n>=1:
        n -= 1
        rData = dict(dataRefineProximity)
        rData['direction'] = 'xy'
        rData['stlFile'] = str(DEFAULT_SHIP_STL)
        rData['outsidePoints'] = list(outsidePoints)

   
    print cSize,cBuffer, nlevel
    pprint.pprint(rData)
    pass

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
    
