#!/usr/bin/python -u

#########################################################################
# Filename: fsPlot.py                                                   #
# Date:     2016-July-01                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:  handle data from openfoam                                   #
#########################################################################

import inspect
import sys, os, glob, pandas, shlex, datetime
import numpy as np
import subprocess, warnings
from copy import deepcopy 
from subprocess import PIPE
from subprocess import Popen

# template for custom metadata to be added to pendas.dataframe
# this will go to: df.fsData
FSDATA = {
'label': None,  # label of this dataframe (user-defined)
'info' : None,  # info-string
'module': None, # module used during read operation, e.g.: loadInternalLoads, loadMotionInfo, .., etc
'args': None    # either a string or dictionary to keep the original input-arguments
}

def tryImport(name):
    # remove ".py"
    newname = os.path.splitext(name)[0]
    try:
        newmodule = __import__(newname)
    except:
        print "Failed to import",name,"... please check your installation"
        raise SystemExit('abort ...')
    return newmodule

# add Popen PIPE to list "p" 
def addPipe(p, cmd):
    try:
        if len(p)==0:
            p.append(Popen(cmd, stdout=PIPE, stderr=PIPE))
        else:
            prev=p[-1]
            p.append(Popen(cmd, stdin=prev.stdout, stdout=PIPE, stderr=PIPE))
            prev.stdout.close()
    except OSError as e:
        print ("cmd:",' '.join(cmd))
        print e
        raise SystemExit('cmd failed to execute ... abort')
        pass
    except:
        print "cmd:",' '.join(cmd)
        print "Unexpected error while constructing subprocess.Popen()"
        raise SystemExit('abort ...')
    
    
# run the command, wait for it to finish, exit if the return code indicates a failure
def runCommand(cmd):
    txt = shlex.split(cmd) if isinstance(cmd, basestring) else cmd
    mycmd=[]
    p=[]
    for cmdarg in txt:
        if (cmdarg=='|'): # open a new pipe
            if len(mycmd)==0:
                print "cmd invalid:",cmd
                raise SystemExit('abort ...')
            addPipe(p, mycmd)
            mycmd=[]
        else:
            mycmd.append(cmdarg)
    if len(mycmd)==0:
        print "cmd invalid:",cmd
        raise SystemExit('abort ...')
    addPipe(p, mycmd)
    txt,err = p[-1].communicate()
    if err:
        print "error:",err
        raise SystemExit('abort ...')
    return txt

def filesOnly(names):
    files=[]
    for f in names:
        if os.path.isfile(f):
            files.append(f)
    return files

def cmd2numpy(cmd, dtype=float):
    """
        exec. cmd and convert output txt to numpy array
    """
    txt=shlex.split(cmd)
    mycmd=[]
    p=[]
    for cmdarg in txt:
        if (cmdarg=='|'): # open a new pipe
            if len(mycmd)==0:
                print "cmd invalid:",cmd
                raise SystemExit('abort ...')
            addPipe(p, mycmd)
            mycmd=[]
        else:
            mycmd.append(cmdarg)
    if len(mycmd)==0:
        print "cmd invalid:",cmd
        raise SystemExit('abort ...')
    addPipe(p, mycmd)

    # using pipe is way faster than StringIO
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, append=1)
        val = np.loadtxt(p[-1].stdout, dtype=dtype)
    # check err(s)
    for proc in p:
        err = proc.stderr.readline()
        if (err!="" and err[0:16]!="tac: write error"):
            print "cmd:",cmd
            print err
            raise SystemExit("cmd failed to execute ... abort")
    return val

def usage():
    print '''
# template for loading openfoam data into pandas.dataframe
    
import sys
sys.path.append("/home/soseng/OpenFOAM/bv/foamBazar/pythonScripts/")
import fsData as fs
from matplotlib import pyplot as plt
if __name__ == "__main__":
    log = fs.loadLogData("-p res -w init,Ux,Uy,Uz", logfiles=['log.run','fsLog'])
    mot = fs.loadMotionInfo("motionInfo", root='./')
    vbm = fs.loadInternalLoads("vbm", root='./', fnames=['my','fz','acc'])    
    '''
    pass

def addslash(name):
    return name if name[-1]=="/" else name + "/"
    
def rmslash(name):
    return name if name[-1]!="/" else name[:-1]
    
# where is postProcessing folder
def postProcessingFolder(objName, root='./'):
    dname = str(objName)
    if os.path.isdir(dname): return addslash(dname)
    basename = os.path.basename(rmslash(objName))
    root = addslash(root)
    dname = root + basename
    if os.path.isdir(dname): return addslash(dname)
    dname = root + "postProcessing/" + basename
    if os.path.isdir(dname): return addslash(dname)
    print "Data directory not found:",objName
    raise SystemExit('abort ...')
    pass

def timeFolder(root='./', sort=True):
    namestr=[]
    nameval=[]
    for key in os.listdir(root):
        try:
            nameval.append(float(key))
            namestr.append(key)
        except:
            pass
    if (len(nameval) and sort==True):
        idx = np.argsort(nameval)
        namestr = [namestr[i] for i in idx]
    return namestr

def postProcessingDatFile(fname, objName=None, root='./'):
    if objName!=None:
        dataFolder = postProcessingFolder(objName, root=root)
        timeNames = timeFolder(root=dataFolder)
    else:
        dataFolder = addslash(root)
        timeNames = []
    if len(timeNames)==0: timeNames=['']    # at least check the current folder
    keyName = os.path.basename(rmslash(fname))
    keyName = os.path.splitext(keyName)[0]
    datFiles = []
    for subdir in timeNames:
        found = filesOnly(sorted(glob.glob(dataFolder + subdir + "/" + keyName + "*.dat")))
        for f in found: datFiles.append(f)
    return datFiles

# concat dataframe and optionally merge xAxis
# When overlap, either keep 'last', 'first', or 'False'
# list_of_data must be of type pandas.dataframe
def concat_and_merge(data, keep='last'):
    if len(data)==0: return data
    if len(data)==1: return data[0]
    # sort the data blockwise
    first = [i.first_valid_index() for i in data]
    idx = np.argsort(first)
    data = [ data[i] for i in idx]
    # convert and merge list to pandas
    data = pandas.concat(data);
    data = data[~data.index.duplicated(keep=keep)]
    return data

def setmetadata(data, label=None, info=None, module=None, args=None):
    data.fsData = deepcopy(FSDATA)
    data.fsData['label'] = label
    data.fsData['info'] = info if info!=None else 'last update: ' + datetime.date.today().strftime("%I:%M%p %B %d, %Y")
    data.fsData['module'] = module
    data.fsData['args'] = args if args!=None else {
        'lastUpdate':datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
    }

# load log data given a list of logfiles/folder
# data will be merged and return as pandas.dataframe
# cmd: is fsPlot.py command line arguments
# e.g.: loadLogData('-p res', logfiles=['log.run0','log.run1',''])
# When overlap, either keep 'last', 'first', or 'False'
def loadLogData(cmd, logfiles=[], root=None, keep='last'):
    """
        def loadLogData(cmd, logfiles=[], root=None, keep='last'):
        
        load data from log-file(s) given "cmd" to pass to fsPlot.py
    """
    fsPlot = tryImport('fsPlot.py')
    allfiles = deepcopy(logfiles) if root==None else [addslash(root) + log for log in logfiles]
    data = []
    for log in allfiles: data.append(fsPlot.loadData(cmd + ' ' + log, dtype='pandas'))
    if len(data):
        data = concat_and_merge(data, keep=keep)
        setmetadata(data, label='logdata', module='loadLogData')
        data.fsData['args'] = {
            'cmd':cmd,
            'logfiles':logfiles,
            'root':root,
            'keep':keep,
            'lastUpdate':datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
        }
    return data

# FIXME: save data for a faster reload, check time stamp on file(s) to reload
def loadMotionInfo(objName, root='./', fname='sixDofDomainBody.dat', keep='last'):
    """
        def loadMotionInfo(objName, root='./', fname='sixDofDomainBody.dat', keep='last'):
        
        load data from: ./postProcessing/<objName>/<time>/<fname>*.dat
    """
    dataFiles = postProcessingDatFile(fname, objName=objName, root=root)
    data = []
    for f in dataFiles:
        header = runCommand('sed -n \"/^# time /p\" ' + f).split()[2:]
        val = cmd2numpy('sed \"/#/d\" ' + f, dtype=float)
        if len(val):
            data.append(pandas.DataFrame(val[:,1:], index=val[:,0], columns=header))
            data[-1].index.name = "motion"
            data[-1].columns.name = "name"
    if len(data):
        data = concat_and_merge(data, keep=keep)
        setmetadata(data, label=objName, module='loadMotionInfo')
        data.fsData['args'] = {
            'objName':objName,
            'root':root,
            'fname':fname,
            'keep': keep,
            'lastUpdate':datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
        }
    return data

def loadInternalLoads(objName, root='./', fnames=None, keep='last'):
    """
        def loadInternalLoads(objName, root='./', fnames=None, keep='last'):

        load data from: ./postProcessing/<objName>/<time>/<fname>*.dat
    """
    def addInfo(data, label, fname):
        setmetadata(data, label=label, module='loadInternalLoads')
        data.fsData['args'] = {
            'objName':objName,
            'root':root,
            'fnames':[fname],
            'keep': keep,
            'lastUpdate':datetime.datetime.utcnow().strftime('%Y%m%d%H%M%S')
        }
    def headerxyz(f):
        header0 = runCommand('sed -n \"/^# t /p\" ' + f).split()[2:]
        headerx = runCommand('sed -n \"/^# x /p\" ' + f).split()[2:]
        headery = runCommand('sed -n \"/^# y /p\" ' + f).split()[2:]
        headerz = runCommand('sed -n \"/^# z /p\" ' + f).split()[2:]
        header0 = [str(val) for val in header0]
        headerx = [float(val) for val in headerx]
        headery = [float(val) for val in headery]
        headerz = [float(val) for val in headerz]
        return (header0,headerx,headery,headerz)
    def appenddata(cmd,header,data):
        val = cmd2numpy(cmd, dtype=float)
        if len(val):
            data.append(pandas.DataFrame(val[:,1:], index=val[:,0], columns=header))
        pass
        
    if fnames==None:
        fnames = ['my','fz']
    elif isinstance(fnames,basestring):
        if fnames.lower() in ['all']:
            fnames = ['fx','fy','fz','mx','my','mz','fCstr','mCstr','fFluid','mFluid','acc']
        else:
            print "loadInternalLoads: unknown options fnames=",fnames
            raise SystemExit('abort ...')
    allData = []
    for fname in fnames:
        if fname in ['fx','fy','fz','mx','my','mz']:
            dataFiles = postProcessingDatFile(fname, objName=objName, root=root)
            data = []
            for f in dataFiles:
                h0,hx,hy,hz = headerxyz(f)
                header = [h0, hx, hy, hz]
                appenddata('sed \"/#/d;s/[()]//g\" '+f,header,data)
            if len(data):
                data = concat_and_merge(data, keep=keep)
                addInfo(data, fname, fname)
                data.index.name = fname
                data.columns.names = ['n','x','y','z']
                allData.append(data)
        elif fname in ['fCstr','mCstr','fFluid','mFluid']:
            dataFiles = postProcessingDatFile(fname, objName=objName, root=root)
            data = {'x':[],'y':[],'z':[]}
            for f in dataFiles:
                h0,hx,hy,hz = headerxyz(f)
                header = [h0, hx, hy, hz]
                for i,cmpt in enumerate(['x','y','z']):
                    appenddata('sed \"/#/d;s/[()]//g\" '+f+' | awk \'{printf (\"%s\",$1); for (i=2+'+str(i)+';i<=NF;i+=3) printf (\" %s\",$i); printf (\"\\n\");}\'',header,data[cmpt])
            if len(data['x']):
                data['x'] = concat_and_merge(data['x'], keep=keep)
                addInfo(data['x'], fname + "x", fname)
                data['x'].index.name = fname + "(x)"
                data['x'].columns.names = ['n','x','y','z']
                allData.append(data['x'])
            if len(data['y']):
                data['y'] = concat_and_merge(data['y'], keep=keep)
                addInfo(data['y'], fname + "y", fname)
                data['y'].index.name = fname + "(y)"
                data['y'].columns.names = ['n','x','y','z']
                allData.append(data['y'])
            if len(data['z']):
                data['z'] = concat_and_merge(data['z'], keep=keep)
                addInfo(data['z'], fname + "z", fname)
                data['z'].index.name = fname + "(z)"
                data['z'].columns.names = ['n','x','y','z']
                allData.append(data['z'])
        elif fname in ['acc']:
            dataFiles = postProcessingDatFile(fname, objName=objName, root=root)
            data = []
            for f in dataFiles:
                header=[]
                for cmpt in runCommand('sed -n \"/^# t /p\" ' + f).split()[2:]:
                    cmpt = cmpt.split("(")[0]
                    if cmpt.lower() in ['linacc','angacc','omega']:
                        header.append(cmpt+'x')
                        header.append(cmpt+'y')
                        header.append(cmpt+'z')
                    elif cmpt in ['R^T']:
                        header.append('RTxx')
                        header.append('RTxy')
                        header.append('RTxz')
                        header.append('RTyx')
                        header.append('RTyy')
                        header.append('RTyz')
                        header.append('RTzx')
                        header.append('RTzy')
                        header.append('RTzz')
                    else:
                        print "acc: unknown format:",cmpt
                        raise SystemExit('abort ...')
                appenddata('sed \"/#/d;s/[()]//g\" '+f,header,data)
            if len(data):
                data = concat_and_merge(data, keep=keep)
                addInfo(data, fname, fname)
                data.index.name = fname
                data.columns.names = ['name']
                allData.append(data)
        else:
            print "unknown data:",fname
    if len(allData)==1:
        return allData[0]
    else:
        return allData

#*** Main execution start here *************************************************
if __name__ == "__main__":
    usage()




