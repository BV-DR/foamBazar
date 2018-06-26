#!/usr/bin/env python

#########################################################################
#                                                                       #
#   This script aims at launching CFD computations                      #
#   on the whole hydrodynamics database                                 #
#########################################################################

import sys
sys.path.append(r"C:\drsvn\drtools\hyel_db_git")
from Sourcesv2.envSettings import env
#Set environment
env["hstar_path"] = r'C:\Program Files\Bureau Veritas\Hstar4Experts\v8.00\bin'
env["database"] =    r'C:\drsvn\dr-tools\Ship_Database\database_argos'
env["output_path"] = r'D:\Etudes\hydroDatabaseCFD'
#env["hmr_path"] = r'D:\_code_drsvn\Homer2\fortran\

from Sourcesv2.createHST import createHST, createHST4OpenFOAM
from Sourcesv2.reg_code import code_to_reg
from Sourcesv2.launchHMR import runFreq , runTime, createFreq, runSpec

import os
from os.path import join
from Pluto.Mesh.MeshHstar import Mesh





createMeshes = True
createfoamCases = True
createLaunchScript = True
launchCFD = False


speed = 5.

lShips = [
             ['C01', 'Full'],
#             ['C03', 'Full'],
#             ['C05', 'Full'],
#             ['C06', 'Full'],
#             ['C07', 'Full'],
#             ['C09', 'Full'],
#             ['C10', 'Full'],
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
    if (createMeshes == True):
      createHST4OpenFOAM(reg, loading, database=r'C:\drsvn\dr-tools\Ship_Database\database_argos')
      folder = join(env["output_path"], code, loading, r'Input_Files')
  
      OFfolder =  join(env["output_path"], code, loading, r'foamMesh')
      if not os.path.exists(OFfolder):
          os.makedirs(OFfolder)
  
      #mesh=Mesh(r"Full.hst")
      mesh=Mesh(join(folder,loading+".hst"))
      print (mesh)
      mesh.makeConformal(0)
      mesh.writeSTL(join(OFfolder,r"ship.stl"))
    
    # next steps
    # symmetrize stl
    # rotate stl to exact trim
    # launch fsMesher
    
    
    if (createfoamCases == True):
       pass
    # launch fsCall
