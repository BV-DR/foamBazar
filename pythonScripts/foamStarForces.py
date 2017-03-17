#!/usr/bin/python
"""
   Plot force output for foamStar post-processing output files
"""

import matplotlib.pyplot as plt

try :
   from TimeDomain import TimeSignals as ts
except ImportError as  e:
   print "Can not load TimeSignals, please add 'waveTimeFrequencyAnalysis' (from repo https://github.com/cmonroy/waveTimeFrequencyAnalysis.git) to your PYTHON_PATH\n"
   raise Exception(e)


if __name__ == "__main__" :

   """
      Example of use :

      FoamStarForces forces.dat                          #->plot all the forces components
      FoamStarForces forces.dat -indexName Fx Fy         #->plot the selected components, based on labels
      FoamStarForces forces.dat -index     0  1          #->plot the selected components, based on indexes

   """
   import argparse
   parser = argparse.ArgumentParser(description='foamStar forces plot')
   parser.add_argument( "forceFile" )
   parser.add_argument('-indexName',  nargs='+', type = str , help='Index to plot' )
   parser.add_argument('-index',  nargs='+', type = int , help='Index to plot' )
   args = parser.parse_args()

   a = ts.read( args.forceFile , reader = "openFoamReader" , field = "total") 

   if args.indexName :
      a[args.indexName].plot()
   elif args.index :
	  a.iloc[:,args.index].plot()
   else:
      a.plot()

   plt.show();
