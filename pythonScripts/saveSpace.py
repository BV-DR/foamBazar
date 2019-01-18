#!/usr/bin/python
import shutil
import os
from os.path import join
from seaKeepingVtkPicture import is_number


def cleanFolder( caseDir , tmin , tmax, verbose = False ) :
    for f in os.listdir(caseDir) :
        if os.path.isdir( join(caseDir, f) ) and is_number(f) :
            if tmin < float(f) < tmax :
                r = join( caseDir , f )
                if verbose : print ("cleaning" + r)
                os.system("rm -rf {}".format(r))
                #shutil.rmtree( r ) too slow !


def cleanCase( caseDir , tmin, tmax, verbose , parallel = "auto" ) :
    """Remove data between tmin and tmax
    """
    #Automatically choose ( parallel vs (serial/notRun))
    if parallel == "auto" :
        if os.path.exists( os.path.join(caseDir , "processor0") ) : parallel = True
        else : parallel = False

    if parallel :
        for p in os.listdir( caseDir ) :
            if "processor" in p : cleanFolder( os.path.abspath( join(caseDir,p )) , tmin , tmax, verbose = verbose )
    else :
        cleanFolder( os.path.abspath( caseDir ) , tmin , tmax , verbose = verbose )


if __name__ == "__main__" :

    import argparse
    parser = argparse.ArgumentParser(description='Remove time step data')
    parser.add_argument( '-tmin'  ,  help='Tmin ' , type = float )
    parser.add_argument( '-tmax'  ,  help='Tmax ' , type = float )
    parser.add_argument( '-case' ,  help='Quantity to plot' , type = str,  default = ".")
    parser.add_argument( '-verbose', '-v', help='verbose ouptut' , action="store_true")

    args = parser.parse_args()
    cleanCase(  caseDir = args.case , tmin =  args.tmin , tmax = args.tmax , verbose = args.verbose  )
