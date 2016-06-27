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
import matplotlib.pyplot as plt
import matplotlib.ticker as tkr
import matplotlib.animation as animation

import pprint
from StringIO import StringIO
from subprocess import PIPE
from subprocess import Popen

DEBUG = False

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
'col' : [-1],           # plot the last column (not use when "iter == True")
'iter' : False,         # show each iteration or not?
'keepAlive' : False,    # keep the windows open and update regularly
'keyOpts' : dict(KEYWORD_OPTS) # which varname to plot?
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
'opts' : dict(CONTROL_OPTS) # what to plot, how to plot, ..., etc ...
}

# this is a template for each data set
PLOTME = {
'fname' : "",   # filename
'title' : "",   # title (from data file)
'cmd' : "",     # command line to extract data
'data' : [],    # data to plot
'axis' : [],    # axis [x0,x1, y0,y1]
'xlabel' : "t [s]",
'ylabel' : ""
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
        print "cmd:",' '.join(cmd)
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
        logdir = os.path.dirname(name) + "/fsLog_" + fname + "/"
        if not os.path.isdir(logdir):
            logdir = "./fsLog_" + fname + "/"
            if not os.path.isdir(logdir):
                # data folder not found ... create from log-file using fsLog.awk
                runCommand([CMD_fsLogAwk, name])
    elif not os.path.isdir(name):
        print "log-file/folder not found:",name
        raise SystemExit('abort ...')
    
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
    parser.add_argument('logfile', nargs='?', help='Read data from file/folder. Log-file, if given in its raw format, will be piped through "fsLog.awk". By default the data is read from: ./fsLog/')
    parser.add_argument('-p', '--plot', metavar='KEY', dest='plot', type=str, help='R|Quantities to plot (comma separated keywords). KEY can\nbe the exact file name, or the follwing pre-defined\nkeywords. ' + "Default: " + re.sub(r'[ \'\[\]]', '', str(CONTROL_OPTS['plot']))  + INFO_PLOTKEYS)
    parser.add_argument('-v','--var', dest='var', type=str, help='R|Name of each selected variable (comma separated names).\nAvail. names are shown in file names after the first\nunderscore, e.g.: "final,Ux,Uz" will select both files\n"Res_final_Ux" and "Res_final_Uz". Plot keyword can be\nspecified using ":" as the separator. Defaults are:\n'
        'co:max,interface\n'
        'ph:min,max\n'
        'res:fsi,init,Ux,Uy,Uz,prgh,fsi\n'
        'nIter:Ux,Uy,Uz,PIMPLE,prgh\n'
        'cont:cumu\n'
        't:curr,meshUpdate\n'
        'f:relax,x,y,z\n'
        'm:relax,x,y,z')
    parser.add_argument('-c','--column', dest='col', help='Column(s) to plot againt time (not use when option --iter is set). The first column is 1. The last column is -1 or 0 (default). (FIXME: not yet implemented) Use negative int to count from the last column. Use comma to define multiple columns.')
    parser.add_argument('-i','--iter', dest='iter', action='store_true', help='Show values at each iteration')
    parser.add_argument('-k','--keep-alive', dest='keepAlive', action='store_true', help='Keep the plot windows open and update regularly')
    parser.add_argument('-t','--trim', metavar='M,N', dest='trim', help='Trim data M lines HEAD and N lines tail, e.g.: -t 2,3 will trim HEAD 2 lines and TAIL 3 lines. Use "s" to trim in second, e.g. -t 0.2s,2s will trim HEAD 0.2 sec and TAIL 2 sec. Use "S" to trim at absolute time positions. Without "s" or "S" only int is allowed.')
    parser.add_argument('-l','--limit', metavar='N', dest='limit', help='Limit the total number of data points to N. Use "s" to set the limit in second, e.g.: --limit 5s will plot only the last 5 sec. of data')

    #FIXME: add option to output data to files
    #FIXME: add option to overwrite ylim
    #FIXME: add option to define line thickness
    #FIXME: add option to define line marker (reduced marker)

    # get the template
    data = dict(DATA)
    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG=True
        pass
    if args.logfile!=None:
        logfile = str(args.logfile)
        data['logdir'] = getlogdir(logfile)
        pass
    if args.plot!=None:
        data['opts']['plot'] = [str(val) for val in args.plot.split(",")]
        pass
    if args.var!=None:
        plot = data['opts']['plot']
        var = [str(val) for val in args.var.split(",")]
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
    if args.iter:
        data['opts']['iter'] = True
        pass
    if args.keepAlive:
        data['opts']['keepAlive'] = True
    if args.trim!=None:
        global TRIMINSECOND
        global TRIMHEADTAIL
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
        global LIMITSAMPLEPOINTS
        global LIMITSAMPLEINSECOND
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

    data = checkAvailable(data)    
    return data

# prepare cmd to read data
def prepareData(fname, opts):
    plotme = dict(PLOTME)
    if opts['iter']:
        # read all data including iter.        
        # awk 'NR==1{t=$1;next}{dt=$1-t;t=$1}' # this will compute deltaT
        # awk 'NR==1{t=$1;next} (NF>1) {dt=($1-t)/(NF-1);t=$1; for (i=2;i<=NF;i++) print t+0.5*dt*(i-2),$i}' fsLog/nIter_PIMPLE
        cmd = "awk 'NR==1{t=$1;next} (NF>1) {dt=($1-t)/(NF-1);t=$1; for (i=2;i<=NF;i++) print t+0.5*dt*(i-2),$i}' "
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

    return plotme

def numpyArrayFromTXTfile(cmd, dtype=float):
    if DEBUG:
        t0=time.time()
        print "DEBUG:",cmd
    txt=shlex.split(cmd)
    mycmd=[]
    p=[]
    for val in txt:
        if (val=='|'): # open a new pipe
            if len(mycmd)==0:
                print "cmd invalid:",cmd
                raise SystemExit('abort ...')
            addPipe(p, mycmd)
            mycmd=[]
        else:
            mycmd.append(val)
            
    if len(mycmd)==0:
        print "cmd invalid:",cmd
        raise SystemExit('abort ...')
    addPipe(p, mycmd)

    # direct piping is way faster than StringIO
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, append=1)
        val = np.loadtxt(p[-1].stdout, dtype=dtype)

    for proc in p:
        err = proc.stderr.readline()
        if (err!="" and err[0:16]!="tac: write error"):
            print "cmd:",cmd
            print err
            raise SystemExit("cmd failed to execute ... abort")

    #txt,err = p[-1].communicate()
    #if err:
    #    print "cmd:",cmd
    #    print "error:",err
    #    raise SystemExit('abort ...')
    #if (len(txt)==0): return np.array([])
    #fakeStream = StringIO(txt)
    #val = np.loadtxt(fakeStream, dtype=dtype)
    
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

def readData(plotme):
    # FIXME: check timeStamp on files and update only when neccessary
    global TRIMHEADTAIL
    global TRIMINSECOND
    global LIMITSAMPLEPOINTS
    global LIMITSAMPLEINSECOND
    fname = plotme['fname']
    cmd = plotme['cmd']
    head = TRIMHEADTAIL[0]
    headInSec = TRIMINSECOND[0]
    tail = TRIMHEADTAIL[1]
    tailInSec = TRIMINSECOND[1]
    print "Read data:",fname
    if (head>0 or tail>0):
        if headInSec:
            if headInSec==1:
                info="trim HEAD for "+str(head)+" sec"
                head=countLinesInSecond(head, fname, fromTail=False)
            else:
                info="trim HEAD at "+str(head)+" sec"
                head=countLinesToPosition(head, fname)
            print info,"("+str(head)+" lines)"
        else:
            head += 2   # we always have one header line in data files
            print "trim HEAD "+str(head)+" lines"
        #
        if tailInSec:
            if tailInSec==1:
                info="trim TAIL for "+str(tail)+" sec"
                tail=countLinesInSecond(tail, fname, fromTail=True)
            else:
                info="trim TAIL at "+str(tail)+" sec"
                tail=countLinesToPosition(tail, fname, fromTail=True)
            print info,"("+str(tail)+" lines)"
        else:
            print "trim TAIL "+str(tail)+" lines"

    cmd = "head -n-"+str(int(tail))+" "+fname+" | tail -n+"+str(int(head))+" | "+cmd
    plotme['data'] = numpyArrayFromTXTfile(cmd)
    
    # data is loaded quick enough, so we simply cut the data afterward 
    limit = LIMITSAMPLEPOINTS
    limitInSec = LIMITSAMPLEINSECOND
    if (limit>0):
        if (limitInSec):
            idx = np.where((plotme['data'][-1,0]-plotme['data'][:,0])<=limit)[0]
            if not len(idx):
                plotme['data']=[]
            else:
                plotme['data'] = plotme['data'][-len(idx):,:]
            print "limit data points to",limit,"sec"
            pass
        else:
            plotme['data'] = plotme['data'][-limit:,:]
            print "limit data points to",limit,"lines"
            pass
    # we can't seem to plot a dot, report empty when we have 1 data  
    if len(plotme['data'])<=1: plotme['data'] = []
    if not len(plotme['data']):
        print "cmd:",cmd
        print "Warning: cmd returns no data ... skip"

    return 

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
            fname = data['logdir'] + key
            if not os.path.isfile(fname):
                print "Data file not found:",fname
                raise SystemExit('abort ...')
            print "Read data from file:",fname  
            data['plotme'].append(prepareData(fname, plotOpts))
    return data

def setPlotAxes(ax, plotme, check=False):
    if not check:
        ax.set_xlabel(plotme['xlabel'])
        ax.set_ylabel(plotme['ylabel'])
        ax.set_title(plotme['title'])
        xmin,xmax, = min(plotme['data'][:,0]), max(plotme['data'][:,0])
        ymin,ymax, = min(plotme['data'][:,1]), max(plotme['data'][:,1])
    else:
        if (plotme['xlabel'] != ax.get_xlabel()):
            print "debug: xlabel has changed?"
        if (plotme['ylabel'] != ax.get_ylabel()):
            print "debug: ylabel has changed?"
            print "plotme['ylabel']:",plotme['ylabel']
            print "ax.get_ylabel():",ax.get_ylabel()
        oldXmin, oldXmax = ax.get_xlim()
        oldYmin, oldYmax = ax.get_ylim()
        xmin, xmax, = min([oldXmin, min(plotme['data'][:,0])]), max([oldXmax, max(plotme['data'][:,0])])
        ymin, ymax, = min([oldYmin, max(plotme['data'][:,1])]), max([oldYmax, max(plotme['data'][:,1])])
        pass
    xspan, yspan = xmax-xmin, ymax-ymin
    if (math.fabs(xspan)<1e-16) : xspan=1e-14
    if (math.fabs(yspan)<1e-16) : yspan=1e-14
    xlim = [xmin-0.01*xspan, xmax + 0.01*xspan]
    ylim = [ymin-0.01*yspan, ymax + 0.01*yspan]
    ax.set_xlim(xlim[0],xlim[1])
    ax.set_ylim(ylim[0],ylim[1])
    pass

def showPlot(data):
    createArray(data)

    fig,ax=plt.subplots()

    plotOK = False    
    for item in data['plotme']:
        readData(item)
        if len(item['data'])==0: continue
        xdata = item['data'][:,0]
        ydata = item['data'][:,1]
        name = item['title'] + " last: " + "{:.2e}".format(ydata[-1])
        setPlotAxes(ax, item, plotOK)
        plt.plot(xdata, ydata, label=name, marker='')
        plotOK=True

    if plotOK:
        ax.yaxis.set_major_formatter(tkr.FormatStrFormatter('%.2e'))
        ax.grid()    
        plt.legend(loc=2)
        plt.show()
    else:
        print "Warning: nothing to plot ... skip"
        
    return data

#*** Main execution start here *************************************************
if __name__ == "__main__":
    data = cmdOptions(sys.argv)
    showPlot(data)
    #pprint.pprint(data)
    
    #    filename = "./fsLog/Courant_interface"
    #    data = np.loadtxt(filename, dtype=float)
    #    print data[:,0]
    #    plt.plot(data[:,0], data[:,1])
    #    plt.show()


