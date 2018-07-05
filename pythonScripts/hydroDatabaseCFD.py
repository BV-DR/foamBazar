#!/usr/bin/env python

#########################################################################
#                                                                       #
#   This script aims at launching CFD computations                      #
#   on the whole hydrodynamics database                                 #
#                                                                       #
#########################################################################

import sys
import subprocess
from Sourcesv2.envSettings import env
#Set environment
env["hstar_path"] = r'/data/I1608251/HydroStar/bin'
env["database"] =    r'/data/I1608251/Ship_Database_june2018/database_argos'
env["output_path"] = r'/home/qitjofivv/hydroDatabase'
#env["hmr_path"] = r'D:\_code_drsvn\Homer2\fortran\

from Sourcesv2.createHST import createHST, getTrim
from Sourcesv2.reg_code import code_to_reg
#from Sourcesv2.launchHMR import runFreq , runTime, createFreq, runSpec

import os
from os.path import join
from Pluto.Mesh.MeshHstar import Mesh




createSTL = False
createMeshes = True
createfoamCases = True
createLaunchScript = True
launchCFD = False


speed = 5.

lShips = [
             ['C01', 'Full'],
#              ['C03', 'Full'],
#             ['C05', 'Full'],
#             ['C06', 'Full'],
#             ['C07', 'Full'],
#             ['C09', 'Full'],
#              ['C10', 'Full'],
#             ['C11', 'Full'],
#             ['C12', 'Full'],
#             ['C13', 'Full'],
#             ['C14', 'Full'],
#             ['C15', 'Full'],
#             ['C17', 'Full'],
#             ['C19', 'Full'],
#             ['C26', 'Full'],
#             ['C29', 'Full'],
         ]


config = { "name" : "1" ,
           "eqss" : 0.01,
           "NB_HSTEP" : 8,
           "dt" : 0.1 ,
        
   "nslam" : 15,
           "SlammingAngle" : 45,
           "GWMVELOCITY" : 2.0,
           "SlammingPart" : 0.25,
           "IDSS_DURATION" : 20.,
           "weather" : "rec34"
         }


for code , loading in lShips :
    reg = code_to_reg(code)

    #runFreq( code , loading, speed  , AmgParam = {"MPAR" : "50 20\n"}, wadim = [ 0.3, 2.4 , 0.3 ] , nbproc = 2 )
    #runSpec( code, loading, speed, config, option = "time" )
    #runSpec( code, loading, speed, config, option = "spec" )
    if (createSTL == True):
      print (reg)
      
      # load intel mpi libraries for HydroStar 
      proc=subprocess.Popen(['module load intel/2015.3.187 intelmpi/5.0.3.048 && env -0'], stdout = subprocess.PIPE, shell = True)
      a = proc.stdout.read()
      for b in a.split(b"\x00")[:-1]:
         (key,  value) = b.decode().split('=',1)
         os.environ[key] = value


      createHST(reg, loading, database=env["database"] , withTrim = False )

      folder = join(env["output_path"], code, loading, r'Input_Files')
  
      OFfolder =  join(env["output_path"], code, loading, r'foamMesh')
      if not os.path.exists(OFfolder):
          os.makedirs(OFfolder)
  
      mesh=Mesh(join(folder,loading+".hst"))
      mesh.makeConformal(0)
      mesh.sym()
      mesh.write(join(folder,"mesh_conformal.hst"))
      mesh.writeSTL(join(OFfolder,r"ship_without_trim.stl"))
    
    # unload intelmpi libraries (incompatible with OpenFOAM mpi library)
      print ("Unload intelmpi libraries")
      proc=subprocess.Popen(['module unload intelmpi/5.0.3.048 && env -0'], stdout = subprocess.PIPE, shell = True)
      a = proc.stdout.read()
      for b in a.split(b"\x00")[:-1]:
         (key,  value) = b.decode().split('=',1)
         os.environ[key] = value

    # load of5.x
      print ("load OpenFOAM environment")
      proc=subprocess.Popen(['module load gcc/4.9.3 openmpi/1.8.4-gcc lapack/3.6.1/gcc/4.9.3 && export FOAM_INST_DIR=/data/I1608251/OpenFOAM && source /data/I1608251/OpenFOAM/OpenFOAM-5.x/etc/bashrc && export LC_ALL=C && env -0'], stdout = subprocess.PIPE, shell = True)
      a = proc.stdout.read()
      for b in a.split(b"\x00")[:-1]:
         (key,  value) = b.decode().split('=',1)
         os.environ[key] = value
 
    # rotate stl to exact trim
      trim=getTrim(reg, loading, database=env["database"] )
      print ("Apllied trim angle: ", trim)
    
      OFRotateCommand = 'surfaceTransformPoints -yawPitchRoll \'(0 ' + str(trim) + ' 0)\' ' + join(OFfolder,r"ship_without_trim.stl") + ' ' + join(OFfolder,r"ship.stl")
      print (OFRotateCommand)
  
      subprocess.call(['surfaceTransformPoints -yawPitchRoll \'(0 ' + str(trim) + ' 0)\' ' + join(OFfolder,r"ship_without_trim.stl") + ' ' + join(OFfolder,r"ship.stl")], shell = True)

 #subprocess.call('transformPoints

    if (createMeshes == True):
       pass

    # launch fsMesher
    
    
    if (createfoamCases == True):
       pass
    # launch fsCall
