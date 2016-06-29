#!/usr/bin/python -u

#########################################################################
# Filename: fsPlot.py                                                   #
# Date:     2016-June-22                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:  plot/monitor log-file data from fsLog.awk                   #
#########################################################################

import re, os, sys, time, math, glob, shlex, subprocess, argparse, configparser
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.animation as anim

from copy import deepcopy
from StringIO import StringIO
from subprocess import PIPE
from subprocess import Popen

DEBUG = False
VERBOSE = False

# with "cat, tail, head" we can manage partial data efficiently (important for large data sets)
TRIMINSECOND = [0, 0]   # 0: trim lines, 1: trim second, 2: trim at absolute time positions 
TRIMHEADTAIL = [0.0, 0.0]

LIMITSAMPLEINSECOND = 0
LIMITSAMPLEPOINTS = 0

# raw log-file is parsed through a log-analyzer (i.e. fsLog.awk)
CMD_fsLogAwk = "fsLog.awk"

# A database for keywords defined in argparser
# The "str" for each key must corresponds to the keywords defined in "DATA"
KEYWORD = {
'co' : "Courant",
'ph' : "Phase",
'res' : "Res",
'nIter' : "nIter",
'cont' : "contErr",
'f' : "fluidForce",
'm' : "fluidMoment",
't' : "timing"
}

# this is for printing usage info with "--help" option 
INFO_PLOTKEYS = '''
[co] Courant number (files: Courant_*)
[ph] Phase fraction (files: Phase_*)
[res] Residual (files: Res_*)
[nIter] No. Iterations (files: nIter_*)
[cont] Continuity error (files: contErr_*)
[t] timing data (files: timing_*)
[f] fluid forces (files: fluidForce_*)
[m] fluid forces (files: fluidMoment_*)
'''

# By default, we select only these quantities
# All are selected, if the list is empty
# Data files have names following this format: <keyword>_<varname>_<varname>_<varname>_...ect
KEYWORD_OPTS = {
'Courant' : ['max','interface'],
'Phase' : ['min','max'],
'Res': ['fsi','init', 'Ux', 'Uy', 'Uz', 'prgh', 'fsi'],
'nIter': ['Ux','Uy','Uz','PIMPLE','prgh'],
'contErr': ['cumu'],
'fluidForce': ['relax','x','y','z'],
'fluidMoment': ['relax','x','y','z'],
'timing' : ['curr','meshUpdate']
}

# 'plot' contains a list of keywords which represent a group of data to plot
CONTROL_OPTS = {
'plot' : ['co', 'res', 'cont', 'ph', 'f', 'm', 't', 'nIter'],
'col' : [-1],               # plot the last column (not use when "iter == True")
'showIter' : 0,             # show each iteration or not? 0: no, 1: yes, 2: plot continuously connected lines
'updateInterval' : 1e15,    # keep the windows open and update regularly (value in seconds)
'keyOpts' : deepcopy(KEYWORD_OPTS) # which varname to plot?
}

# file(s) are loaded in "def checkAvailable(data):"
DATA = {
'logdir' : "./fsLog/",  # default sub-folder to look for data
'Time' : "Time",        # file name containing time data
'Courant' : [],         # names of files: "Courant_<name>"
'Phase' : [],           # names of files: "Phase_<name>"
'Res' : [],             # names of files: "Res_{init,final}_<name>"
'nIter' : [],           # names of files: "nIter_<name>"
'contErr' : [],         # names of files: "contErr_<name>"
'fluidForce' : [],      # names of files: "fluidForce_<name>"
'fluidMoment' : [],      # names of files: "fluidMoment_<name>"
'timing' : [],          # names of files: "timing_<name>"
'opts' : deepcopy(CONTROL_OPTS) # what to plot, how to plot, ..., etc ...
}

LINEPROPERTYCOUNTER = {
'marker': 0,
'style': 0,
'width':0,
'color':0
}

LINEPROPERTY = {
'marker' : None,    # line marker
'style' : None,     # line style
'width' : None,     # line width
'color' : None      # line color
}

# this is a template for each data set
PLOTME = {
'fname' : "",           # filename
'title' : "",           # title (from data file)
'cmd' : "",             # command line to extract data
'curLine': None,        # keep track of lines read from file
'lastModified': 0,      # keep track of last modified time
'data' : [],            # data to plot
'xlim' : [],            # axis [x0,x1]
'ylim' : [],            # axis [y0,y1]
'xlabel' : "t [s]",
'ylabel' : "",
'line' : deepcopy(LINEPROPERTY)
}

# for a nicer output of "usage info" in argparser
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

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
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        print "Failed to exec.:",cmd
        raise SystemExit('abort ...')
        pass
    except OSError:
        print cmd
        raise SystemExit('executable not found ... abort')
        pass
    pass    

# check log-file and return "logdir"
def getlogdir(name):
    logdir=name
    if os.path.isfile(name):
        fname = os.path.basename(name)
        # where is data folder?
        logdir = os.path.dirname(name) + "./fsLog_" + fname + "/"
        if not os.path.isdir(logdir):
            logdir = "./fsLog_" + fname + "/"
            if not os.path.isdir(logdir):
                # data folder not found ... create from log-file using fsLog.awk
                runCommand([CMD_fsLogAwk, name])
        return logdir

    # name is neither a file or a folder
    # we check try to add prefix: fsLog_
    if not os.path.isdir(name):
        fname = os.path.basename(name)
        logdir = os.path.dirname(name) + "./fsLog_" + fname + "/"
        if not os.path.isdir(logdir):
            print "log-file/folder not found:",name
            raise SystemExit('abort ...')
        return logdir
      
    # name is a folder, just proceed
    
    # add trailing slash "/"
    if not logdir[-1]=="/":
        logdir = logdir + "/"
    return logdir

def filesOnly(names):
    files=[]
    for f in names:
        if os.path.isfile(f):
            files.append(f)
    return files

# check available data
def checkAvailable(data):
    logdir = str(data['logdir'])
    if not os.path.isdir(logdir):
        print "Log-data not found in folder:", logdir
        raise SystemExit('abort ...')
    for key in KEYWORD:
        keyName = KEYWORD[key]
        data[keyName] = filesOnly(glob.glob(logdir + keyName + "_*"))    
    return data

def cmdOptions(argv):
    parser = argparse.ArgumentParser(formatter_class=SmartFormatter)
    parser.add_argument('-d', '--debug', action='store_true', help='Run in DEBUG mode')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show more comprehensive info while running')
    parser.add_argument('logfile', metavar='logfile|folder',nargs='?', help='Read data from log-file/folder. Log-file, if given in its raw format, will be piped through "fsLog.awk". By default the data is read from: ./fsLog/')
    parser.add_argument('-p', '--plot', metavar='key', dest='plot', type=str, help='R|Quantities to plot (comma separated keywords). KEY can\nbe the exact name of the data file, or the follwing\npre-defined keywords. ' + "Default: " + re.sub(r'[ \'\[\]]', '', str(CONTROL_OPTS['plot']))  + INFO_PLOTKEYS)
    parser.add_argument('-w','--with', metavar='var', dest='withvar', type=str, help='Name of each selected variable (comma separated names). Avail. names are shown in file names after the first underscore, e.g.: "final,Ux,Uz" will select both files "Res_final_Ux" and "Res_final_Uz". Plot keyword can be specified using ":", e.g. -w res:fsi,init,Uz')
    parser.add_argument('-c','--column', metavar='n', dest='col', help='Column(s) to plot againt time (not use when option --iter is set). The first column is 1. The last column is 0 (default). Use comma to define multiple columns.')
    parser.add_argument('-u','--update', metavar='n', nargs='?', const=30, default=None, dest='updateInterval', action='store', type=float, help='Update each data set every N seconds, (default: 30). We are checking file time stamp. The minimum update interval is limited to 0.1 s (i.e. 10 fps)')
    parser.add_argument('-t','--trim', metavar='m,n', dest='trim', help='Trim data M lines HEAD and N lines tail, e.g.: -t 2,3 will trim HEAD 2 lines and TAIL 3 lines. Use "s" to trim in second, e.g. -t 0.2s,2s will trim HEAD 0.2 sec and TAIL 2 sec. Use "S" to trim at absolute time positions. Without "s" or "S" only int is allowed.')
    parser.add_argument('--limit', metavar='n', dest='limit', help='Limit the total number of data points to N. Use "s" to set the limit in second, e.g.: --limit 5s will plot only the last 5 sec. of data')
    parser.add_argument('-i','--iter', nargs='?', const='split', default=None, metavar='cont', dest='iter', action='store', help='Show values at each iteration. By default iterations between each time step are shown as separated line segments. With option "-i cont" all values are shown in one continuous line. ')
    parser.add_argument('-a','--axis', metavar='val', dest='axis', help='Set axis range xlim and/or ylim. The format is [x.y]:min,max e.g.: -a y:-1.5,2  or -a y:-1.5,2,x:0,14')
    parser.add_argument('-l','--line', metavar='key', dest='line', type=str, help='Set line propreties: style (key s), linewidth (key w), marker (key m). Any line style supported by matplotlib can be used here, e.g.: "-" solid, "--" dashed, ":" dotted, "-." dash-dotted, "x" marker, "o" marker, etc. For multiple lines use comma e.g.: -l s:-,--,w:0.5,1')
    parser.add_argument('-o','--output', metavar='prefix', dest='output', type=str, help='Output each data set to file')
    
    #FIXME: add option to output data to files/png/raw
    
    # get the template
    data = deepcopy(DATA)
    args = parser.parse_args()

    global DEBUG
    global VERBOSE
    global TRIMINSECOND
    global TRIMHEADTAIL
    global LIMITSAMPLEPOINTS
    global LIMITSAMPLEINSECOND
    global PLOTME
    global LINEPROPERTY
    
    if args.debug:
        DEBUG=True
        VERBOSE=True
        pass
    if args.verbose:
        VERBOSE=True
        pass
    if args.logfile!=None:
        logfile = str(args.logfile)
        data['logdir'] = getlogdir(logfile)
        pass
    if args.plot!=None:
        data['opts']['plot'] = [str(val) for val in args.plot.split(",")]
        pass
    if args.withvar!=None:
        plot = data['opts']['plot']
        var = [str(val) for val in args.withvar.split(",")]
        reset = {}
        key=None
        for i in var:
            j = i.split(":")
            if len(j)==2:
                if j[0] in KEYWORD:
                    key=KEYWORD[j[0]]
                    if key in reset:
                        del reset[key]
            if len(j)>2:
                print "Invalid option(s): ",i
                raise SystemExit('abort ...')
            i = j[-1]
            if key==None:
                for k in plot:
                    if k in KEYWORD:
                        k=KEYWORD[k]
                        if not k in reset:
                            data['opts']['keyOpts'][k]=[]
                            reset[k]=True
                        data['opts']['keyOpts'][k].append(i)
            else:
                if not key in reset:
                    data['opts']['keyOpts'][key]=[]
                    reset[key]=True
                data['opts']['keyOpts'][key].append(i)
                pass
        pass
    if args.col!=None:
        txt = args.col.split(",")
        data['opts']['col'] = []
        for i in txt:
            try:
                i = int(i)
            except ValueError:
                print "Warning: ignore invalid int value for --column option:",i
            if (i==0): i=-1
            data['opts']['col'].append(i)
        if not len(data['opts']['col']): data['opts']['col']=[-1]
        pass
    if args.iter!=None:
        data['opts']['showIter'] = 2 if (args.iter=='cont') else 1
        pass    
    if args.updateInterval!=None:
        data['opts']['updateInterval'] = args.updateInterval
    if args.trim!=None:
        txt = args.trim.split(",")
        allGood = True
        sec=list(TRIMINSECOND)
        tmp=list(TRIMHEADTAIL)
        if (len(txt) != 2): # we need 2 numbers HEAD and TAIL
            allGood=False
        else:
            i=0
            for val in txt:
                try:
                    if (val[-1]=='s'):  # trim in second
                        sec[i] = 1
                        val=val[:-1]
                    elif (val[-1]=='S'): # trim at absolute time in second
                        sec[i] = 2
                        val=val[:-1]
                    if (sec[i]): # float is allowed when trim in second
                        tmp[i]=float(val)
                    else: # int is required when trim lines
                        tmp[i]=int(val)
                    if (tmp[i]<0) : raise ValueError
                    # for 'S' the first number must be smaller than the second
                    if (i==1 and sec[0]==2 and sec[1]==2):
                        if (tmp[0]>=tmp[1]): raise ValueError
                    i += 1
                except ValueError:
                    allGood=False
                    break
        if (not allGood):
            print "Warning: ignore invalid option: --trim",args.trim
        else:
            TRIMINSECOND = list(sec)
            TRIMHEADTAIL = list(tmp)
        pass
    if args.limit!=None:
        val = args.limit
        sec=LIMITSAMPLEINSECOND
        tmp=LIMITSAMPLEPOINTS
        allGood = True
        try:
            if (val[-1]=='s'):  # trim in second
                sec = 1
                val=val[:-1]
            if (sec): # float is allowed when trim in second
                tmp=float(val)
            else: # int is required when trim lines
                tmp=int(val)
            if (tmp<=0): raise ValueError
        except ValueError:
            allGood=False
        if (not allGood):
            print "Warning: ignore invalid option: --limit",args.limit
        else:
            LIMITSAMPLEINSECOND = sec
            LIMITSAMPLEPOINTS = tmp
        pass
    if args.axis!=None:
        word={'x':'xlim','y':'ylim'}
        def invalidOption():
            print "Invalid option(s): --axis",args.axis
            raise SystemExit('abort ...')
        txt = args.axis.split(",")
        if not (len(txt)==2 or len(txt)==4): invalidOption()
        while len(txt):
            j=txt[0].split(":")
            if (len(j)!=2) or (not j[0] in word): invalidOption()
            if (j[1]=="") or (txt[1]==""): invalidOption()
            key=word[j[0]]
            try:
                minval = float(j[1])
                maxval = float(txt[1])
                if (minval>=maxval): raise 
            except:
                invalidOption()
            PLOTME[key]=[minval, maxval]
            del txt[0],txt[0]
        pass
    if args.line!=None:
        word={'s':'style', 'w':'width', 'm':'marker', 'c':'color'}
        wtype={'style':str, 'width':float, 'marker':str, 'color':str}
        def invalidOption(info=''):
            if info!='': print info
            print "Invalid option(s): --line",args.line
            raise SystemExit('abort ...')
        def addvalue(key, value, dtype=None):
            if key==None: invalidOption()
            LINEPROPERTY[key].append(value)
            pass       
        txt = args.line.split(",")
        key = None
        for val in txt:
            j=val.split(":")
            if len(j)>2: invalidOption()
            if len(j)==2:
                if (j[0]=='' and j[1]==''):
                    addvalue(key, ":")
                    continue
                if not j[0] in word: invalidOption("unknown keyword")
                key=word[j[0]]
                if LINEPROPERTY[key]==None: LINEPROPERTY[key]=[]
                addvalue(key, j[1])
            elif len(j)==1:
                addvalue(key,j[0])
            else:
                invalidOption("this should never happed")
        for key in LINEPROPERTY:
            if LINEPROPERTY[key]!=None:
                LINEPROPERTY[key]=map(wtype[key],LINEPROPERTY[key])
        pass

    if args.output!=None:
        print args
        print "output each data set to file: ... not yet implemented ... abort ..."
        raise SystemExit

    data = checkAvailable(data)    
    return data

# prepare cmd to read data
def prepareData(fname, opts):
    plotme = deepcopy(PLOTME)
    if opts['showIter']==1: # plot iter.(s) using disconnected lines
        # read all data including iter.        
        # awk 'NR==1{t=$1;next}{dt=$1-t;t=$1}' # this will compute deltaT
        # awk 'NR==1{t=$1;next} (NF>1) {dt=($1-t)/(NF-1);t=$1; for (i=2;i<=NF;i++) print t+0.5*dt*(i-2),$i}' fsLog/nIter_PIMPLE
        cmd = "awk 'FNR==1{t=$1;c=0.75;next} (NF>1) {dt=($1-t)/(NF-1);t=$1; for (i=2;i<NF;i++) {print t+c*dt*(i-2),$i}; print t+c*dt*(NF-2),\"nan\"}' "
        print "disconnected"
        pass
    elif opts['showIter']==2: # plot iter.(s) using one continuous line
        # read all data including iter.        
        cmd = "awk 'FNR==1{t=$1;c=0.75;next} (NF>1) {dt=($1-t)/(NF-1);t=$1; for (i=2;i<NF;i++) print t+c*dt*(i-2),$i}' "
        print "continuous"
        pass
    else:
        # read only specific col(s)
        # cmd = "awk '/#/{next;} (NF>=???) {print $1,$?,$?,...}' "
        cmd = 'print $1'
        for col in opts['col']:
            cmd += ',$NF' if col==-1 else ',$' + str(int(col))
        cmd = '{' + cmd + '}'
        maxCol = int(max(opts['col']))
        if (maxCol>0): cmd = "(NF>=" + str(maxCol) + ") " + cmd 
        cmd = "awk '/#/{next;} " + cmd + "'" 
        pass

    p = Popen(["head","-n1",fname], stdout=PIPE, stderr=PIPE)
    txt,err = p.communicate()
    if err:
        print "error:",err
        raise SystemExit('abort ...')
    txt = re.sub(r'[#_\n]', '', txt)
    txt = re.sub(r':.*', '', txt)
    plotme['title'] = txt
    plotme['fname']=fname
    plotme['cmd']=cmd

    # set line properties
    # the global has a list of all user-defined properties, and we just rotate them
    global LINEPROPERTYCOUNTER
    for key in LINEPROPERTY:
        if LINEPROPERTY[key]!=None:
            i = LINEPROPERTYCOUNTER[key]
            plotme['line'][key]=LINEPROPERTY[key][i]
            if len(LINEPROPERTY[key])>1:
                LINEPROPERTYCOUNTER[key] = (i+1) % len(LINEPROPERTY[key])

    return plotme

def numpyArrayFromTXTfile(cmd, dtype=float):
    if DEBUG:
        print "DEBUG:",cmd
        t0=time.time()
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

    for proc in p:
        err = proc.stderr.readline()
        if (err!="" and err[0:16]!="tac: write error"):
            print "cmd:",cmd
            print err
            raise SystemExit("cmd failed to execute ... abort")
   
    if DEBUG: print "DEBUG: timing",time.time()-t0
    return val

# how many lines of data for the given sec. of time
# the results is most likely NOT exact
# return 0 when failed
def countLinesInSecond(val, fname, fromTail=False):
    if val==0: return int(0)
    if fromTail:
        cmd="tac " + fname + " | awk 'FNR==1 {init=0;} /#/ {next;} {if (init==0) {t0=$1+1e-12; init=1;} } {if ((t0-$1+tol) >= "+str(val)+") {print NR; exit;}}' "
    else:
        cmd="awk 'FNR==1 {init=0;} /#/ {next;} {if (init==0) {t0=$1-1e-12; init=1;} } {if (($1-t0) >= "+str(val)+") {print NR; exit;}}' " + fname
    count = numpyArrayFromTXTfile(cmd, dtype=int)
    if not count: count = 0
    return int(count)

def countLinesToPosition(val, fname, fromTail=False):
    if val==0: return int(0)
    if fromTail:
        cmd="tac " + fname + " | awk 'function abs(x){return ((x<0.0)?-x:x)} /#/ {next;} {if ($1 > "+str(val)+"){ prev=$1; next; } } { if (abs($1-"+str(val)+")<abs(prev-"+str(val)+")) { print NR-1; exit; } else { print NR-2; exit; } }' "
    else:
        cmd="awk 'function abs(x){return ((x<0.0)?-x:x)} /#/ {next;} {if ($1 < "+str(val)+"){ prev=$1; next; } } { if (abs($1-"+str(val)+")<abs(prev-"+str(val)+")) { print NR; exit; } else { print NR-1; exit; } }' " + fname
    count = numpyArrayFromTXTfile(cmd, dtype=int)
    if not count: count = 0
    return int(count)

def countHeaderLines(fname):
    cmd="awk 'FNR==1 {count=0;} !NF || /^[ \t]*#/ {count+=1; next;} {print NR-1; exit;}' " + fname
    count = numpyArrayFromTXTfile(cmd, dtype=int)
    if DEBUG: print "DEBUG: header lines:",count
    return 0 if count.size==0 else count

def fileHasBeenModified(fname, lastModified):
    oldtime = lastModified
    lastModified = os.path.getmtime(fname);
    return lastModified>oldtime, lastModified
    
# At the first call full data will be loaded
# At subquence call(s), we read and append the latest data only
def readData(plotme, latestOnly=True, verbose=True):
    needUpdate,plotme['lastModified'] = fileHasBeenModified(plotme['fname'], plotme['lastModified'])
    needUpdate = needUpdate or (not latestOnly)
    if needUpdate==False: return False,0
    global TRIMHEADTAIL
    global TRIMINSECOND
    global LIMITSAMPLEPOINTS
    global LIMITSAMPLEINSECOND
    fname = plotme['fname']
    cmd = plotme['cmd']
    if not latestOnly or (plotme['curLine']==None):
        plotme['curLine']==None
        head = TRIMHEADTAIL[0]
        headInSec = TRIMINSECOND[0]
    else:
        head = plotme['curLine']
        headInSec = 0   # trim to current line
    tail = TRIMHEADTAIL[1]
    tailInSec = TRIMINSECOND[1]
    if verbose: print "Read data:",fname
    nHeaderLines = countHeaderLines(fname)
    if (head>0 or tail>0):
        if headInSec:
            if headInSec==1:
                info="trim HEAD for "+str(head)+" sec"
                head=countLinesInSecond(head, fname, fromTail=False)
            else:
                info="trim HEAD at "+str(head)+" sec"
                head=countLinesToPosition(head, fname)
            if verbose: print info,"("+str(head)+" lines)"
        else:
            head += nHeaderLines + 1
            if verbose: print "trim HEAD "+str(head)+" lines"
        #
        if tailInSec:
            if tailInSec==1:
                info="trim TAIL for "+str(tail)+" sec"
                tail=countLinesInSecond(tail, fname, fromTail=True)
            else:
                info="trim TAIL at "+str(tail)+" sec"
                tail=countLinesToPosition(tail, fname, fromTail=True)
            if verbose: print info,"("+str(tail)+" lines)"
        else:
            if verbose: print "trim TAIL "+str(tail)+" lines"

    cmd = "head -n-"+str(int(tail))+" "+fname+" | tail -n+"+str(int(head))+" | "+cmd

    foundNewData=True
    startIndex = 0
    if (plotme['curLine']==None):
        plotme['data'] = numpyArrayFromTXTfile(cmd)
        plotme['curLine'] = int(head) + len(plotme['data']) - nHeaderLines + 1
        if DEBUG: print "DEBUG: cururent Line",plotme['curLine']
    else:
        newdata = numpyArrayFromTXTfile(cmd)
        if newdata.size:
            startIndex=len(plotme['data'])
            if newdata.ndim==1: newdata=np.array([newdata])
            plotme['data'] = np.append(plotme['data'], newdata, axis=0)
            plotme['curLine'] += len(newdata)
            if DEBUG: print "DEBUG: cururent Line",plotme['curLine']
        else:
            # same old data, nothing else to do ...
            foundNewData=False
            return foundNewData,startIndex

    # data is loaded quick enough, so we simply cut the data afterward 
    limit = LIMITSAMPLEPOINTS
    limitInSec = LIMITSAMPLEINSECOND
    if (limit>0):
        alldata = plotme['data']
        oldSize = len(alldata)
        idx = int(limit)
        if (limitInSec):
            idx = np.where((alldata[-1,0]-alldata[:,0])<=limit)[0]
            idx = len(idx)
            if not idx:
                plotme['data']=[]
            else:
                plotme['data'] = alldata[-idx:,:]
            if verbose: print "limit data points to",limit,"sec"
            pass
        else:
            plotme['data'] = alldata[-idx:,:]
            if verbose: print "limit data points to",idx,"lines"
            pass
        startIndex = max([0, startIndex - oldSize + min([oldSize,idx])])

    if not len(plotme['data']):
        if verbose: print "cmd:",cmd
        if verbose: print "Warning: cmd returns no data ... skip"
        
    return foundNewData, startIndex

# load data into np.array for plot, return data
def createArray(data):
    plotOpts = data['opts']
    if not len(plotOpts['plot']): # nothing to do 
        return data
    data['plotme'] = []
    for key in plotOpts['plot']:
        if key in KEYWORD:
            key = KEYWORD[key]
            for fname in data[key]:
                ok = True
                for subkey in os.path.basename(fname).split("_")[1::]:
                    ok = ok & (subkey in plotOpts['keyOpts'][key])
                if ok:
                    data['plotme'].append(prepareData(fname, plotOpts))
            pass
        else:
            # try to read data directly from file
            if os.path.isfile(key):
                fname = key
            elif os.path.isfile(data['logdir'] + os.path.basename(key)):
                fname = data['logdir'] + os.path.basename(key)
            else:
                print "Data file not found:",fname
                raise SystemExit('abort ...')
            data['plotme'].append(prepareData(fname, plotOpts))

    return data

def setPlotAxes(ax, plotme, check=False):
    def setx(xlim):
        ax.set_xlim(xlim[0],xlim[1])
    def sety(ylim):
        ax.set_ylim(ylim[0],ylim[1])
    def compute(dat):
        minval = min(dat)
        maxval = max(dat)
        span = maxval-minval
        if (math.fabs(span)<1e-16): span=1e-14
        return minval,maxval,span

    hasX = False
    hasY = False
    if len(plotme['xlim'])==2:
        setx(plotme['xlim'])
        hasX = True
    if len(plotme['ylim'])==2:
        sety(plotme['ylim'])
        hasY = True
    xdata=plotme['data'][:,0]
    ydata=plotme['data'][:,1]
    expand = 0.01   # expand axis by 1%
    if not check:
        ax.set_xlabel(plotme['xlabel'])
        ax.set_ylabel(plotme['ylabel'])
        ax.set_title(plotme['title'])
        if not hasX:
            xmin,xmax,xspan = compute(xdata)
            xlim = [xmin-expand*xspan, xmax+expand*xspan]
            setx(xlim)
        if not hasY:
            ymin,ymax,yspan = compute(ydata)
            ylim = [ymin-expand*yspan, ymax + expand*yspan]
            sety(ylim)
        return

    # check for existing axis
    if (plotme['xlabel'] != ax.get_xlabel()):
        print "debug: xlabel has changed?"
    if (plotme['ylabel'] != ax.get_ylabel()):
        print "debug: ylabel has changed?"
        print "plotme['ylabel']:",plotme['ylabel']
        print "ax.get_ylabel():",ax.get_ylabel()
    oldXmin, oldXmax = ax.get_xlim()
    oldYmin, oldYmax = ax.get_ylim()
    if not hasX:
        xmin,xmax,xspan = compute(xdata)
        xmin = xmin-expand*xspan if (xmin<oldXmin) else oldXmin
        xmax = xmax+expand*xspan if (xmax>oldXmax) else oldXmax
        setx([xmin,xmax])
    if not hasY:        
        ymin,ymax,yspan = compute(ydata)
        ymin = ymin-expand*yspan if (ymin<oldYmin) else oldYmin
        ymax = ymax+expand*yspan if (ymax>oldYmax) else oldYmax
        sety([ymin,ymax])    
    pass

def showPlot(data):
    global DEBUG
    global VERBOSE
    createArray(data)

    # FIXME: can we do subplot? rather not manually, it makes the user-interface very complex
    nrows=1
    ncols=1
    plot_number=1
    
    # no more than 10 fps
    updateInterval = min([data['opts']['updateInterval'], 0.1])
    updateInterval *= 1e3 # value in millisec. (required by FuncAnimation)

    fig=plt.figure()
    ax=fig.add_subplot(nrows,ncols,plot_number)
    plt.ioff()  # don't show the plot until we call plt.show()
    ax.grid(True)
    ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%.2e'))
    lines = [None] * len(data['plotme'])
    nLines = [0]    # we must use list here to keep nLines known in "def update"
    def update(frame):
        if DEBUG: print "DEBUG: frame",frame
        for i,item in enumerate(data['plotme']):
            foundNewData,startIndex = readData(item, verbose=(DEBUG or VERBOSE))
            if (not foundNewData or item['data'].ndim==1): continue
            name = item['title'] + " last: " + "{:.2e}".format(item['data'][:,1][-1])
            setPlotAxes(ax, item, check=nLines[0])
            # startIndex seems to be unneccessary 
            #xdata = item['data'][startIndex:,0]
            #ydata = item['data'][startIndex:,1]
            xdata = item['data'][:,0]
            ydata = item['data'][:,1]
            if lines[i]==None:
                if DEBUG: print "DEBUG: create Line2D",i
                linemarker = item['line']['marker']
                linestyle = item['line']['style']
                linewidth = item['line']['width']
                linecolor = item['line']['color']
                lines[i], = plt.plot(xdata,ydata,label=name,marker=linemarker)
                nLines[0]+=1
                if linestyle!=None: lines[i].set_linestyle(linestyle)
                if linewidth!=None: lines[i].set_linewidth(linewidth)
                if linecolor!=None: lines[i].set_color(linecolor)
            else:
                lines[i].set_label(name)
                lines[i].set_xdata(item['data'][:,0])
                lines[i].set_ydata(item['data'][:,1])
            plt.legend(loc=0)
            plt.draw()
        pass

    a = anim.FuncAnimation(fig, update, repeat=False, interval=updateInterval)
    plt.show()

    return data

#*** Main execution start here *************************************************
if __name__ == "__main__":
    data = cmdOptions(sys.argv)
    showPlot(data)
    
    #    filename = "./fsLog/Courant_interface"
    #    data = np.loadtxt(filename, dtype=float)
    #    print data[:,0]
    #    plt.plot(data[:,0], data[:,1])
    #    plt.show()


