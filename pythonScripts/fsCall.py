#!/usr/bin/env python

#########################################################################
# Filename: fsCall.py                                                   #
# Date:     2018-June-15                                                #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################
#                                                                       #
#   This file is a sample script to call meshing or template routines   #
#   from foamBazar. Parameters and arguments are defined as dict that   #
#       can be edited within this script before calling routines.       #      
#                                                                       #
#########################################################################
import os, sys
import numpy as np

#Import meshing or template routine from foamBazar
# WARNING : foamBazar/pythonScripts must be added to PYTHONPATH !
from fsTemplate2D import fsTemplate2D


#This whole section can be scripted by user

#Sample of parameters dict
myParams = {'nProcs' : 12,
            'OFversion' : 5,
            'case' : 'section_10',
            'gridLevel' : 1,
            'symmetry' : True,
            'meshDir' : '../Sections/sym',
            'endTime' : 4,
            'timeStep' : 0.01,
            'writeInterval' : 10,
            'outputForces' : True,
            'dispSignal' : '../dispSignals/case2.dat',
            'sideRelaxZone' : -1,
            'translateLength' : 0.5,
            'gravity' : 0.0
            }

#Sample of arguments dict
myArgs = {'overwrite' = True }


#Calling routine for meshing or template here
fsTemplate2D(userParam = myParams, userArgs = myArgs)
