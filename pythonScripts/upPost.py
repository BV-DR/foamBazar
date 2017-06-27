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
        print input
        break
    elif input in ['exit','break','cancel']:
        sys.exit()
    else:
        print 'Input is invalid, please select time serie from list.'

# data = pd.read_table(os.path.join(dir_path,input+'.dat'),sep=' ',comment='#',header=None)
# data = openFoam.openFoamReader(os.path.join(dir_path,input+'.dat'))
data = openFoam.foamStarMotion(os.path.join(dir_path,input+'.dat'))

# print data[['3']]

upCross.upCrossMinMax(data.loc[:,'3'])

# for col in xrange(data.shape[1]):
    # res = upCross.upCrossMinMax(data[[col+1]])
# print res.values




