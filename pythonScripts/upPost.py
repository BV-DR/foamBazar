#!/usr/bin/env upython

#########################################################################
# Filename: upPost.py                                                 #
# Date:     2017-May-02                                                 #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################

import os, subprocess, time, math
import numpy as np
import sys
import pandas as pd
from TimeDomain import upCross
from Reader import openFoam

#Input parameters


#Code
dir_path = os.getcwd()
dir_file = os.listdir(dir_path)

inp_list = []
for file in dir_file:
    if file.endswith('.dat'):
        inp_list += [file[:-4]]

print 'Choose time serie to post-process:'
for ts in inp_list: print '- '+ts

while True:
    input = raw_input('>>>')
    if input in inp_list:
        break
    elif input in ['exit','break','cancel']:
        sys.exit()
    else:
        print 'Input is invalid, please select time serie from list.'

if input in ['fx','fy','fz','mx','my','mz']:
    data = openFoam.foamStarMotion(os.path.join(dir_path,input+'.dat'))
else:
    data = openFoam.openFoamReader(os.path.join(dir_path,input+'.dat'))

resu = pd.DataFrame()
for idx in data.columns.values:
    mini = min(data.loc[:,idx])
    maxi = max(data.loc[:,idx])
    if mini!=maxi:
        res = upCross.upCrossMinMax(data.loc[:,idx])
        tmp = pd.DataFrame( data={ input+'_'+idx+'_min': res.Minimum, input+'_'+idx+'_max': res.Maximum} )
        resu = pd.concat([resu,tmp],axis=1)

resu = resu.transpose()
resu.to_csv(input+'_upcross.csv',sep=';')
