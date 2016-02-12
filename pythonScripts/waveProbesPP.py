"""

  Post-process foamStar wave probes  data file. ("surfaceElevation.dat")

"""


from  matplotlib import pyplot as plt
from  matplotlib import animation
import numpy as np


def makeAnimation(filename , movieName = None, nShaddow = 20, xRatio= 1.0) :
   """
      Make a 2D animation of the wave field.
   """
   print "Making animation file : " , movieName
   nShaddow = max(1,nShaddow)
   data = np.loadtxt(filename)
   with open(filename , "r") as f :
      f.readline() ; f.readline() ; xVal = [float(x) for x in f.readline().split()[2:]]

   #Remove x = 0.0
   if 0.0 in xVal :
      i = xVal.index(0.0)
      xVal.pop(i)
      data = np.delete(data, i+1 , axis = 1)
   xVal = np.array(xVal)
   xVal /= xRatio

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
     ax.set_title("{}s".format(data[itime,0]) )
     for s in range(nShaddow):
        if itime > s :
           ls[s].set_data( xVal , data[itime-s,1:] )
     return ls

   ani = animation.FuncAnimation( fig , run, range(len(data)), blit=True, interval=30, repeat=False)

   if movieName is None :
      plt.show()
   else :
      mywriter = animation.FFMpegWriter(fps = 25 , codec="libx264")
      ani.save( movieName +'.mp4',writer=mywriter)


def getTimeSignal(filename) :
   from TimeDomain import TimeSignals as ts
   with open(filename , "r") as f :
      f.readline() ; f.readline() ; label = f.readline().split()[2:]
   d = np.loadtxt(filename)
   data = d[:,1:]
   data.shape = -1 , len(label)
   xAxis = d[:,0]
   return ts.TimeSignals( data , xAxis = xAxis , columnTitles = label )