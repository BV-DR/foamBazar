#!/usr/bin/python
"""
   Plot Fourier output from foamStar
"""

try :
   from TimeDomain import TimeSignals as ts
except ImportError as  e:
   print "Can not load TimeSignals, please add 'shared-code/Python' (from repo http://svn-drsvn.eua.bvcorp.corp/dr/shared-code) to your PYTHON_PATH\n"
   raise Exception(e)


if __name__ == "__main__" :

   """
      Example of use :

      FoamStarFourier forces.dat -index 1 -period 12.3    
		#->plot all the harmonics of Fy (index 1) which has a period of 12.3 s
     
      FoamStarFourier forces.dat -index 0 -period 12.3 -harmo 1  
		# ->plot the 1st harmonics (harmo 1) of Fx (index 0) which has a period of 12.3 s

   """
   import argparse
   parser = argparse.ArgumentParser(description='FoamStar Fourier plot')
   parser.add_argument( "forceFile" )
   parser.add_argument('-index',  nargs='+', type = int , help='Index to plot' )
   parser.add_argument('-period',  nargs='+', type = float , help='Index to plot' )
   parser.add_argument('-harmo',  nargs='+', type = int , help='Index to plot' )
   args = parser.parse_args()

   a = ts.read( args.forceFile , reader = "openFoamReader" , field = "total") 

   slFFT= ts.slidingFFT( a, args.period, signalIndex = args.index[0] )

   if args.harmo :
      slFFT.plotTS( args.harmo )
   else:
      slFFT.plotTS( "all" )
