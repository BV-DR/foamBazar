#!/usr/bin/env python

#########################################################################
# Filename: fsClone.py                                               #
# Date:     2017-July-27                                                 #
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

# abenhamou: 2017-july-27

DEFAULT_PARAMS = {
'caseDir' : 'newCase'
}

#*** Data structure for user input data
class UserInput(object):
    def __init__(self, dir, SWcase=False, WVcase=False):
    
        self.dir        =  dir
        self.SWcase     =  SWcase
        self.WVcase     =  WVcase
   
def cmdOptions(argv):
    # default global parameters
    CMD_showLog = 'echo "\nlog-file: ./log.clone"; tail -45 ./log.clone; echo "Please see log-file: ./log.clone"'
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', dest='foldername', help='set folder to clone e.g.: fsClone.py -f MyCase ')
    parser.add_argument('-s', dest='SWcase', action='store_true', help='option to create CrankNicolson still-water inputs in cloned case')
    parser.add_argument('-w', dest='WVcase', action='store_true', help='option to create wave inputs in cloned case')
    args = parser.parse_args()
    
    if args.foldername==None:
        raise SystemExit('ERROR: no input folder defined ! Please use -f option to define folder to clone.')
    elif args.SWcase and args.WVcase:
        raise SystemExit('ERROR: -s and -w option cannot be used together, please pick only one !')
    elif not (args.SWcase or args.WVcase):
        raise SystemExit('ERROR: -s or -w option should be selected, please pick one !')
    else:
        inputdata = UserInput(dir=args.foldername, SWcase=args.SWcase, WVcase=args.WVcase)
        
    return inputdata

def cloneFolder(data,time):
    
    baseDir = data.dir.split('_')[0]
    if data.SWcase:
        data.newDir = baseDir+'_sw_CN'
    elif data.WVcase:
        data.newDir = baseDir+'_wave'
    else:
        raise SystemExit('ERROR: this error should not happen !')

    process = 'mkclone.sh '+data.dir+' '+data.newDir+' '+str(time)
    print process
    process += ' 2>&1 | tee log.clone'
    
    subprocess.call(process, shell=True)
    if checkError('log.clone'):
        raise SystemExit("ERROR : Something wrong happend in case cloning, check 'log.clone' file")
    else:
        os.remove('log.clone')

def getJobName(data):
    filename = data.newDir+'/run.sh'
    rf = open(filename)
    for line in rf:
        if '#SBATCH -J' in line:
            job = line.split()[-1]
            break
    rf.close()
    jsplit = job.split('_')
    bjob = jsplit[0]+'_'+jsplit[1]
    
    return job, bjob
    
def setInputs(data):
    rc = '\n'
    
    if data.SWcase:
        # fvSolution
        filename = data.newDir+'/system/fvSolution'
        setValue(filename, r'//EulerCells  EulerCells;', 'EulerCells  EulerCells;')
        # fvSchemes
        filename = data.newDir+'/system/fvSchemes'
        setValue(filename, 'default      Euler;', 'default      CrankNicolson 0.9;')
        # run.sh
        filename = data.newDir+'/run.sh'
        oname, bname = getJobName(data)
        setValue(filename, oname, bname+'_sw_CN')
    elif data.WVcase:
        # fvSolution
        filename = data.newDir+'/system/fvSolution'
        setValue(filename, r'//EulerCells  EulerCells;', 'EulerCells  EulerCells;')
        # fvSchemes
        filename = data.newDir+'/system/fvSchemes'
        setValue(filename, 'default      Euler;', 'default      CrankNicolson 0.9;')
        # dynamicMeshDict
        filename = data.newDir+'/constant/dynamicMeshDict'
        os.rename(filename,filename+'_old')
        odf = open(filename+'_old','r')
        ndf = open(filename,'w')
        for line in odf:
            if 'dampingCoeff' in line:
                ndf.write(r'//'+line)
                line = odf.next()
                while not '}' in line:
                    ndf.write(r'//'+line)
                    line = odf.next()
                ndf.write(r'//'+line)
            else:
                ndf.write(line)
        odf.close()
        ndf.close()
        # waveProperties
        filename = data.newDir+'/constant/waveProperties'
        os.rename(filename,filename+'_old')
        owp = open(filename+'_old','r')
        nwp = open(filename,'w')
        for line in owp:
            if 'waveType    noWaves;' in line:
                nwp.write(r'    waveType    streamFunction;'+rc)
            elif 'domainX0Coeffs' in line:
                nwp.write(line)
                line = owp.next()
                while not 'relaxationZone' in line:
                    nwp.write(line)
                    line = owp.next()
                nwp.write(r'    waveType    noWaves;'+rc)
                nwp.write(line)
            else:
                nwp.write(line)
        owp.close()
        nwp.close()
        # run.sh
        filename = data.newDir+'/run.sh'
        oname, bname = getJobName(data)
        setValue(filename, oname, bname+'_wave')
    
#*** Main execution start here *************************************************
if __name__ == "__main__":
    startTime = time.time()  
    
    dat = cmdOptions(sys.argv)
    timeFolders = getFoamTimeFolders(dat.dir,inprocs=True)
    cloneFolder(dat,timeFolders[-1])
    setInputs(dat)
    
    endTime = time.time()
    print 'Case cloning completed in %d minutes' % ((endTime-startTime)/60)
    print "Don't forget to edit controlDict file !!!" 
    
