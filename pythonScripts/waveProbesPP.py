"""

  Post-process foamStar wave probes  data file. ("surfaceElevation.dat")

"""
from  matplotlib import pyplot as plt
from  matplotlib import animation
import numpy as np
import pandas as pd
import os


def makeAnimation(filename , movieName = None, nShaddow = 20, xRatio= 1.0, version = "foamStar" , rate = 1) :
   """
      Make a 2D animation of the wave field.
   """
   print "Making animation file : " , movieName
   nShaddow = max(1,nShaddow)
   
   data = pd.read_csv(filename, comment = "#" , header=None , delim_whitespace=True, dtype = float).as_matrix()

   if version == "foamStar" :
      with open(filename , "r") as f :
         f.readline() ; f.readline() ; xVal = [float(x) for x in f.readline().split()[2:]]
      #Remove x = 0.0
      if 0.0 in xVal :
         i = xVal.index(0.0)
         xVal.pop(i)
         data = np.delete(data, i+1 , axis = 1)
      xVal = np.array(xVal)
      xVal /= xRatio

   else :
      xVal = np.linspace(0,1 , len(data[1,:]) -1 )

   fig, ax = plt.subplots()
   ls = []
   for i in range(nShaddow) :
      if i == 0 :
         color = "black"
      else :
         color = "blue"
      ltemp,  = ax.plot([], [], lw=1 , alpha = 1-i*1./nShaddow , color = color)
      ls.append(ltemp)
   ax.grid(True)
   ax.set_xlim( min(xVal) , max(xVal) )
   ax.set_ylim( -1.5 , +1.5 )
   ax.set_xlabel("x")
   ax.set_ylabel("Elevation(m)")

   def run(itime):
     ax.set_title("{}s".format(data[itime*rate,0]) )
     for s in range(nShaddow):
        if itime > s :
           ls[s].set_data( xVal , data[ rate*(itime - s),1:] )
     return ls

   ani = animation.FuncAnimation( fig , run, range(len(data)), blit=True, interval=30, repeat=False)

   if movieName is None :
      plt.show()
   else :
      mywriter = animation.FFMpegWriter(fps = 25 , codec="libx264")
      ani.save( movieName +'.mp4',writer=mywriter)


def getTimeSignal(filename, version = "foamStar", label = None) :
   if label is None :
      with open(filename , "r") as f :
         if version == "foamStar" :
            f.readline() ; f.readline() ; label = f.readline().split()[2:]
         else :
            label = [ "{:0002}".format(i+1) for i in range(len(f.readline().split())-1) ]

   if not os.path.exists(filename) :
      print filename, 'does not exist'
   return  pd.read_csv( filename , header = 0 , names = label , index_col = 0 , engine = "c"  , delim_whitespace = True , comment = "#")


if __name__ == "__main__" :

   print "Run"
   #makeAnimation( r"\\10.67.24.192\bigr\openFoam\2Dwave\Run_revI\g1_foamExtend_Speed_0.0\surfaceElevation.dat" ,  rate = 5 , version = "foamExtend" , nShaddow = 0)

   test = getTimeSignal ( r"\\10.67.24.192\bigr\openFoam\2Dwave\Run_rev1\g1_foamStar_Speed_0.0\postProcessing\waveProbes\0\surfaceElevation.dat" )

   try :
      from pyplotTools import dfSlider
      test.columns = map(float , test.columns)
      print test
      dfSlider(test)
   except :
      pass
