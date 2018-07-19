#!/usr/bin/env upython

#########################################################################
# Filename: fsPlot.py                                                   #
# Date:     2016-June-22                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:  plot/monitor log-file data from fsLog.awk                   #
#########################################################################

import os, sys, math, tempfile, pathlib, gmshScript, subprocess
import warnings, pprint
import numpy as np
import csv

from copy import deepcopy
from io import StringIO
from subprocess import PIPE

class _stl_files:
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.__dict__)
        return info
    def __init__( self, fname):
        pass
    
class _gmsh_GEO_FROM_YZ_DATA:
    '''
    dat-file contains offset points defined in a table of (yz) or (xyz)
    Given (yz) we create a stl on the extruded surface. Possible app(s) are
    e.g.: 2D forced roll, 2D water entry/exit
    Note: we don't check input file here ... we simply assume that the
    given file is valid.
    '''
    GMSH_CMD='gmsh -2 -format stl -o #FOUT #FGEO'
    GEO_TEMPLATE = r'''// fname: #FNAME
Geometry.OldNewReg = 0;
thk=#THK; xmin=-0.5*thk; p0=newp; pSize=#PSIZE;
//Point(p0+#PID) = {xmin, #Y, #Z};
l0=newl; lSize=1; Spline(l0+0) = {p0:p0+pSize-1};
Line(l0+lSize) = {p0+pSize-1,p0};
Extrude {thk, 0, 0} { Line{l0:l0+lSize}; Recombine; }
Transfinite Line {4, 5} = #NX Using Progression 1;
Transfinite Line {2, 6} = #NY Using Progression 1;
Transfinite Line {1, 3} = #NS Using Progression 1;
Transfinite Surface "*"; Recombine Surface "*";
Line Loop(1) = {6, 3}; Plane Surface(3) = {1};
Line Loop(2) = {2, 1}; Plane Surface(4) = {2};
Physical Surface("correct_surface_normal") = {1,2,3,-4};
EOF
'''
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\n'+pprint.pformat(self.__dict__)
        return info
    def __init__( self, fname):
        assert os.path.isfile(fname), "file not found: " + fname
        numCols=0
        with open(fname) as f:
            reader = csv.reader(f, delimiter=' ', skipinitialspace=True)
            firstRow = next(reader)
            numCols = len(firstRow)
        if (numCols==2):
            data=np.loadtxt(fname)
            data=self.redistribute_data_points(data)
            ymax,zmax=data.max(axis=0)
            ymin,zmin=data.min(axis=0)
            pSize=len(data)
            thk=math.ceil(ymax)
            nx=1
            ny=2*math.ceil((ymax-ymin)/(thk/nx))
            ns=math.ceil(10*pSize)
            geo=deepcopy(_gmsh_GEO_FROM_YZ_DATA.GEO_TEMPLATE)
            geo=geo.replace("#FNAME",fname)
            geo=geo.replace("#THK",str(thk))
            geo=geo.replace("#PSIZE",str(pSize))
            geo=geo.replace("#NX",str(nx))
            geo=geo.replace("#NY",str(ny))
            geo=geo.replace("#NS",str(ns))
        else:
            raise(Exception('unknown data: '+fname))
        self.fname=fname
        self.data=data
        self.geo=geo
    def redistribute_data_points(self,origData=None):
        if (origData is None):
            data=deepcopy(self.data)
        else:
            data=deepcopy(origData)
        ymin,zmin=data.min(axis=0)
        ymax,zmax=data.max(axis=0)
        if (math.fabs(ymin)<1e-14):
            sym=deepcopy(data[::-1][:-1])
            sym[:,0] *= -1;
            data=np.concatenate((sym,data), axis=0)
        return data

    def write(self,fname):
        ffile=os.path.basename(fname)
        fdir=os.path.dirname(fname)
        geo=deepcopy(self.geo)
        s='//Point(p0+#PID) = {xmin, #Y, #Z};'
        pointdata = '';
        for row,x in enumerate(self.data):
            tmp=deepcopy(s[2:])
            tmp=tmp.replace('#PID',str(row))
            tmp=tmp.replace('#Y',str(x[0]))
            tmp=tmp.replace('#Z',str(x[1]))
            pointdata += tmp + "\n"
        geo=geo.replace(s,pointdata[:-1])
        geofile = fname + ".geo";
        subprocess.call('mkdir -p "' + fdir + '"', shell=True)
        subprocess.call('cat << EOF > "' + geofile + '"\n' + geo, shell=True)
        fd1,tmp=tempfile.mkstemp(suffix='.msh')
        cmd=deepcopy(_gmsh_GEO_FROM_YZ_DATA.GMSH_CMD)
        cmd=cmd.replace("format stl","format msh")
        cmd=cmd.replace("#FOUT",tmp)
        cmd=cmd.replace("#FGEO",geofile)
        subprocess.call(cmd, shell=True)
        fd2,tmp2=tempfile.mkstemp(suffix='.geo')
        subprocess.call('cat << EOF > "'+tmp2+'"\n Merge "'+tmp+'";', shell=True)
        cmd=deepcopy(_gmsh_GEO_FROM_YZ_DATA.GMSH_CMD)
        cmd=cmd.replace("#FOUT",fname)
        cmd=cmd.replace("#FGEO",tmp2)
        subprocess.call(cmd, shell=True)
        os.close(fd1)
        os.close(fd2)
        os.remove(tmp)
        os.remove(tmp2)
        pass
        
class stlObjClass:
    '''
    class for common stl-operations in OpenFOAM
    '''
    def __init__(self, builder, fout='./constant/triSurface/stlObj.stl'):
        self.builder=deepcopy(builder)
        self.fout=fout
    def __repr__(self): # for printing to terminal
        info = pprint.pformat(type(self))
        info += '\nbuilder: ' + pprint.pformat(type(self.builder))
        info += '\nfout   : ' + self.fout
        return info
    def write(self,fname=None):
        if (fname is None):
            fname=self.fout
        self.builder.write(fname);
        pass
    
def buildFrom(fname):
    '''
    create stl from file(s): filetype are detected based on its extension
    '''
    # guess from extention
    if fname.lower().endswith('.stl'):
        return stlObjClass(_stl_files(fname))
    elif fname.lower().endswith('.dat'):
        return stlObjClass(_gmsh_GEO_FROM_YZ_DATA(fname))
    else:
        raise(Exception('Unknown file-format: '+fname))
        pass

#*** Main execution start here *************************************************
if __name__ == "__main__":
    print(str(sys.argv[0:]))
    

