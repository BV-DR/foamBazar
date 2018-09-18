#!/usr/bin/env python

#########################################################################
# Filename: fsTools.py                                                  #
# Date:     2017-July-27                                                #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import os, subprocess, re
import sys, argparse, configparser
from io import StringIO
import numpy as np
import pandas as pd
import math as mt
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
        print( 'error: ', error.decode('ascii'))
        print('abort ...')
        os._exit(1)
    boundingBox = boundingBox.decode('ascii').split(' ')
    boundingBox = [float(i) for i in boundingBox]
    if verbose: print( "   ",boundingBox)
    return boundingBox
    
def findCFDBoundingBox(case, verbose=True):
    if verbose: print( "Compute CFD bounding box:")
    p = subprocess.Popen("checkMesh -case "+case+" | grep 'Overall domain bounding box' | sed \"s/.*box (//;s/[(,)]//g\" ", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    boundingBox,error = p.communicate()
    if error:
        print( 'error: ', error.decode('ascii'))
        print('abort ...')
        os._exit(1)
    boundingBox = boundingBox.decode('ascii').split(' ')
    boundingBox = [float(i) for i in boundingBox]
    if verbose: print( "   ", boundingBox)
    return boundingBox 

def runCommand(cmd,showlog):
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError:
        subprocess.call(showlog, shell=True)
        print('abort ...')
        os._exit()
        pass
    except OSError:
        print( cmd)
        print('executable not found ... abort')
        os._exit(1)
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
        os._exit(1)
    
def checkError(file):
    error=False
    with open(file, 'r') as f:
        for line in f:
            if 'ERROR' in line: error=True
    
    return error
    
def findSTLPatches(stlFile):
    res = []
    with open(stlFile) as f: lines = f.readlines()
    for item in lines:
         if item.lstrip().startswith('solid'):
            res.append(item.lstrip().split('solid')[-1].strip())
    return res

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
    
def readSections(inputFile,sections=[],sym=False):
    
    sdict = {}
    
    with open(inputFile) as f:
        if len(sections)==0:
            sect = re.findall(r'(#Section \d+)',f.read())
            for s in sect:
                sections += [int(i) for i in s.split() if i.isdigit()]

    with open(inputFile,'r') as f:
        lines = f.read()
    
    for isect in sections:
        string = r'(#Section '+str(isect)+r'\n)+([^#]*)'
        table = re.findall(string, lines)
        table = list(table[0])
        rows = StringIO(table[1])
        data = pd.read_csv(rows,header=None,sep=r'\s+',names=['y','z'])
        
        sdict[isect] = data[data.y>=0]
        sdict[isect] = sdict[isect].assign(x= sdict[isect]['y'].values * 0.0)
        sdict[isect] = sdict[isect].reindex(sorted(sdict[isect]),axis=1)
        
        if sym:
            rsect = rsec = ((sdict[isect]*[1,-1,1]).iloc[1:,:]).iloc[::-1]
            sdict[isect] = rsect.append(sdict[isect],ignore_index=True)
            
    return sdict

def createSectionStl(sdict,redistribute=False):
    if not os.path.exists('geo'): os.makedirs('geo')
    if not os.path.exists('stl'): os.makedirs('stl')
    
    for isect in sdict.keys():
        print('Section =',isect)
        name = 'section_'+str(isect)
        mysect = sdict[isect]
        
        fgeo=os.path.join(r'geo/',name+'.geo');
        fstl=os.path.join(r'stl/',name+'.stl');
    
        #Create .geo file
        thk = 10.0
        
        if redistribute:
            #Redistribute points
            s = np.sqrt((mysect.diff()**2).sum(axis=1)).cumsum()
            ds = 0.5 * np.sqrt((mysect.diff()*2).sum(axis=1)).mean()
            ss = np.linspace(s.min(), s.max(), mt.ceil((s.max()-s.min())/ds))
            
            #interpolate according to curvilinear abscissa
            fx = interp.interp1d(s,mysect.x); x = fx(ss)
            fy = interp.interp1d(s,mysect.y); y = fy(ss)
            fz = interp.interp1d(s,mysect.z); z = fz(ss)
            
            ns = mt.ceil(16*len(ss))
        else:
            x = mysect.x
            y = mysect.y
            z = mysect.z
            
            ns = mt.ceil(16*len(z))

        #create symmetric
        interpsect = pd.DataFrame({ 'x': np.concatenate((x,x[1:])),
                                    'y': np.concatenate(( -1.0*y[:0:-1],y)),
                                    'z': np.concatenate(( z[:0:-1],z)) })
        

        pSize = len(interpsect)
        
        lSize=1.0
        nx=1.0
        ny=2.0*mt.ceil((interpsect.max().y-interpsect.min().y)/(thk/nx))

        print('writing to file: ', fgeo)
        f = open(fgeo,'w')
        f.write('// Section {:d}\n'.format(isect))
        f.write('thk={:.6e}; xmin=-0.5*thk;\n'.format(thk))   
        f.write('p0=newp; pSize={:d};\n'.format(pSize))
        for i in range(pSize):
            f.write('Point(p0+{:d}) = {{xmin, {:.6e}, {:.6e}}};\n'.format(i,interpsect.loc[i].y,interpsect.loc[i].z))
        f.write('l0=newl; lSize={:.2f};\n'.format(lSize))
        f.write('Spline(l0+0) = {{p0:p0+{:d}}};\n'.format(pSize - 1))
        f.write('Line(l0+lSize) = {p0+pSize-1,p0};\n')
        f.write('Extrude {thk, 0, 0} {\n')
        f.write(' Line{l0:l0+lSize};\n')
        f.write(' Recombine;\n')
        f.write('}\n')
        f.write('Transfinite Line {{4, 5}} = {:.5f} Using Progression 1;\n'.format(nx))
        f.write('Transfinite Line {{2, 6}} = {:.5f} Using Progression 1;\n'.format(ny))
        f.write('Transfinite Line {{1, 3}} = {:.5f} Using Progression 1;'.format(ns))
        f.write('Transfinite Surface "*";\n')
        f.write('Recombine Surface "*";\n')
        f.close()

        #Create STL
        subprocess.call('gmsh -2 -format stl -o "'+fstl+'" "'+fgeo+'"', shell=True)
