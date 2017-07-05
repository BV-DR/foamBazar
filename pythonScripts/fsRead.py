import os
import numpy as np
import pandas as pd
from waveTimeFrequencyAnalysis.TimeDomain import TimeSignals as ts

dicoPost = {
             "motions" : 'motions',
             "sixDofDomainBody" : 'motionInfo',
             "surfaceElevation.dat" : 'waveProbe',
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

def fsRead(case,res):
    tlist = os.listdir(os.path.join(case,'postProcessing', dicoPost[res]))
    #sort tlist in ascending order
    ntlist = len(tlist)
    
    data = [{} for i in range(ntlist)]
    tdt = np.empty(0)
    for t in xrange(ntlist):
        data[t] = ts.read(os.path.join(case,'postProcessing',dicoPost[res],str(tlist[t]),res+'.dat'),reader="openFoamReader")
        dt = round(data[t].index[-1]- data[t].index[-2],10)
        tdt = np.append(tdt,dt)
    print tlist
    print tdt
    for idt in np.unique(tdt):
        dataC = pd.concat([data[i] for i in np.where(tdt==idt)[0]])
        dataC.set_index(dataC.index.values.round(10),inplace=True)
        #remove index duplicates (or erase new indexes ??)
        dataC = dataC.groupby(dataC.index).first()
        dataC.to_csv(os.path.join(case,res+'_C.csv'),sep=';')
        #add data to tuple or any other kind of object
    
    return dataC
    
    
