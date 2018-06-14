#!/usr/bin/env python

#########################################################################
# Filename: fsTools.py                                                  #
# Date:     2017-July-27                                                #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import os, subprocess, time, math
import numpy as np
import sys, argparse, configparser
import pprint

# abenhamou: 2017-july-27

def setValue(filename, variable, value, pout=False):
    if type(variable) is not str: variable = str(variable)
    if '"' in variable: variable = variable.replace('"','\\"')
    if '/' in variable: variable = variable.replace('/','\\/')
    if type(value) is not str: value = str(value)
    if '"' in value: value = value.replace('"','\\"')
    if '/' in value: value = value.replace('/','\\/')
    if pout:  print('sed -i "s/'+variable+'/'+value+'/g" '+filename)
    subprocess.call('sed -i "s/'+variable+'/'+value+'/g" '+filename, shell=True)

# foamFile may be compressed i.e. ".gz"  
def foamFileExist(filename):
    found = os.path.isfile(filename)
    if not found:
        found = os.path.isfile(filename + '.gz')
    return found
    
def findBoundingBox(stlFile, verbose=True):
    stlFile = stlFile.split('.stl')[0] #remove .stl extension
    if verbose: print( "Compute STL bounding box: "+stlFile+'.stl')
    p = subprocess.Popen("surfaceCheck "+stlFile+'.stl'+" | grep '^Bounding Box :' | sed \"s/.*: (//;s/[(,)]//g\" ", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    boundingBox,error = p.communicate()
    if error:
        print( 'error: ', error)
        raise SystemExit('abort ...')
    boundingBox = boundingBox.decode('ascii').split(' ')
    boundingBox = [float(i) for i in boundingBox]
    if verbose: print( "   ",boundingBox)
    return boundingBox
    
def findCFDBoundingBox(case, verbose=True):
    if verbose: print( "Compute CFD bounding box:")
    p = subprocess.Popen("checkMesh -case "+case+" | grep 'Overall domain bounding box' | sed \"s/.*box (//;s/[(,)]//g\" ", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    boundingBox,error = p.communicate()
    if error:
        print( 'error: ', error)
        raise SystemExit('abort ...')
    boundingBox = boundingBox.decode('ascii').split(' ')
    boundingBox = [float(i) for i in boundingBox]
    if verbose: print( "   ", boundingBox)
    return boundingBox 

def runCommand(cmd,showlog):
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        subprocess.call(showlog, shell=True)
        raise SystemExit('abort ...')
        pass
    except OSError:
        print( cmd)
        raise SystemExit('executable not found ... abort')
        pass
    pass
    
def getFoamTimeFolders(dir,constant=False,inprocs=False):
    found = []
    if inprocs: dir += '/processor0/'
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
    
def getBool(string):
    if string in ['True','true','T','t','1']:
        return True
    elif string in ['False','false','F','f','0']:
        return False
    else:
        print( 'Invalid boolean entry : '+str(string))
        raise SystemExit('')
    
def checkError(file):
    error=False
    with open(file, 'r') as f:
        for line in f:
            if 'ERROR' in line: error=True
    
    return error
    
def findSTLPatches(stlFile):
    p = subprocess.Popen("grep '^[ \\t]*\<solid\>' "+stlFile+" | sed 's/solid//g' | tr '\n' ' ' | sed 's/^[ \t]*//;s/ \+/ /g;s/\s*$//g' "  , stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    patches,error = p.communicate()
    if error:
        print('error: ', error)
        raise SystemExit('abort ...')
    return patches.decode('ascii').split(' ')

def translateStl(inputStl, val, outputStl):
    val = "("+str(val[0])+" "+str(val[1])+" "+str(val[2])+")"
    print("Translate stl by "+val+": " + inputStl)
    subprocess.call("surfaceTransformPoints -translate '"+val+"' " + inputStl + " " + outputStl + " > /dev/null", shell=True)
    return True

def rotateStl(inputStl, heading, outputStl):
    yaw = heading-180.0
    if str(yaw) == '0.0': return False
    print("Rotate stl (0 0 " + str(yaw) + "): " + inputStl)
    subprocess.call("surfaceTransformPoints -rollPitchYaw '(0 0 "+str(heading-180.0)+")' " + inputStl + " " + outputStl + " > /dev/null", shell=True)
    return True

def createBoxStl(BB,name):
    Xmin,Ymin,Zmin,Xmax,Ymax,Zmax=BB[0],BB[1],BB[2],BB[3],BB[4],BB[5]
    tol = (Xmax-Xmin)*1e-6
    filename = "./constant/triSurface/" + name
    print("Creating stl: " + filename)
    print("   ",BB)
    subprocess.call("surfaceTransformPoints -scale '("+str(Xmax-Xmin-tol)+" "+str(Ymax-Ymin-tol)+" "+str(Zmax-Zmin-tol)+")' fsMesher/fsMesher_box.stl "+filename+" > /dev/null", shell=True)
    subprocess.call("surfaceTransformPoints -translate '("+str(Xmin+0.5*tol)+" "+str(Ymin+0.5*tol)+" "+str(Zmin+0.5*tol)+")' "+filename+" "+filename+" > /dev/null", shell=True)
