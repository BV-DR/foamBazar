#!/usr/bin/env python

#########################################################################
# Filename: fsRead.py                                                   #
# Date:     2017-July-06                                                #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import os
import numpy as np
import pandas as pd
import waveTimeFrequencyAnalysis.Reader as rd

dicoPost = {
             "motions" : 'motions',
             "sixDofDomainBody" : 'motionInfo',
             "surfaceElevation" : 'waveProbes',
             "fx" : 'vbm',
             "fy" : 'vbm',
             "fz" : 'vbm',
             "mx" : 'vbm',
             "my" : 'vbm',
             "mz" : 'vbm',
             "PTS_localMotion_pos" : 'localMotion',
             "PTS_localMotion_vel" : 'localMotion',
             "PTS_localMotion_acc" : 'localMotion',
             "forces" : 'forces',
             "fFluid" : 'vbm',
             "mFluid" : 'vbm',
             "fCstr" : 'vbm',
             "mCstr" : 'vbm',
             "acc" : 'vbm'
             }

#-----------------------------------------------------------------------------#
#                                   fsRead                                    #
#-----------------------------------------------------------------------------#
# This function can be used to read any foamStar output from postProcessing   #
# folder and store results in a DataFrame. The following options can be used: #
# - split: set True to parse results when multiple time steps are used in CFD #
#           case. With this option, the function returns a list of time steps #
#            and a list of DataFrame (one for each time step.)                #
# - csv : set True to create a CSV file corresponding to each DataFrame.      #
#-----------------------------------------------------------------------------#

def fsRead(case,res,split=False,csv=False):
    #read list of time directories and sort in ascending order
    tlist = np.array(os.listdir(os.path.join(case,'postProcessing', dicoPost[res])))
    tlist = tlist[np.argsort(tlist.astype(np.float))]
    ntlist = len(tlist)
    
    #read of data and remove overlapping parts
    data = [{} for i in range(ntlist)]
    if split: tdt = np.empty(0)
    for t in xrange(ntlist):
        data[t] = rd.dfRead(os.path.join(case,'postProcessing',dicoPost[res],str(tlist[t]),res+'.dat'),reader="openFoamReader")
        if t>0: data[t] = data[t][data[t].index>data[t-1].index[-1]]
        if split:
            dt = round(data[t].index[-1] - data[t].index[-2],10)
            tdt = np.append(tdt,dt)
    
    #concatenate and store data in pandas DataFrame
    if split:
        unik = np.unique(tdt)[::-1]
        dataT = []
        for idt in unik:
            dataTmp = pd.concat([data[t] for t in np.where(tdt==idt)[0]])
            if csv: dataTmp.to_csv(os.path.join(case,res+'_'+str(idt)+'.csv'),sep=';')
            dataT.append(dataTmp)
        return unik, dataT
    else:
        dataC = pd.concat([data[i] for i in xrange(ntlist)])
        if csv: dataC.to_csv(os.path.join(case,res+'.csv'),sep=';')
        return dataC
