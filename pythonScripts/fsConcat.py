#!/usr/bin/env python

#########################################################################
# Filename: fsConcat.py                                                 #
# Date:     2018-September-07                                           #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################
# This function can be used to concatenate foamStar results files from  #
#                different time folder into a single file               #
#########################################################################

import os, sys, argparse, shutil
import numpy as np
import pandas as pd
import droppy.Reader.openFoam as of

def cmdOptions(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--folder', dest='folder', default='postProcessing', help='Path to postProcessing folder')
    parser.add_argument('-r','--reverse', dest='reverse', action='store_true', help='is used, the new runs will overlap previous runs')
    args = parser.parse_args()
    return args
        
def concatData(args):
    flist = np.array(os.listdir(args.folder))
    for f in flist:
        print(f)
        tlist = np.array(os.listdir(os.path.join(args.folder,f)))
        newtlist = []
        for t in tlist:
            try:
                newtlist.append(t)
            except:
                continue
        newtlist = np.array(newtlist)
        newtlist = newtlist[np.argsort(newtlist.astype(np.float))]
        
        orgdir = os.path.join(args.folder,f,'ORG')
        if not os.path.exists(orgdir): os.makedirs(orgdir)
        
        for i, t in enumerate(newtlist):
            f0 = of.OpenFoamReadForce(os.path.join(args.folder,f,t,'forces.dat'))
            if i>0:
                if args.reverse: fc = f0.combine_first(fc)
                else: fc = fc.combine_first(f0)
            else:
                fc = f0  
            shutil.move(os.path.join(args.folder,f,t), os.path.join(orgdir,t))
       
        newdir = os.path.join(args.folder,f,'0')
        if not os.path.exists(newdir): os.makedirs(newdir)
        of.OpenFoamWriteForce(fc,os.path.join(newdir,'forces.dat'))

#*** Main execution start here *************************************************
if __name__ == "__main__":
    args = cmdOptions(sys.argv[1:])
    concatData(args)