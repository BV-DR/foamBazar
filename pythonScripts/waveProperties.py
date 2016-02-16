import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector
import numpy as np
import os


"""

   Routines to create the waveProperties input file (+blendingZone setSet input)
   
   TODO : Extend to irregular sea-states

"""


class RelaxZone(object) :
   """
      Relaxation zone definition
   """
   def __init__(self, name , relax , waveCondition , origin = None , orientation = None, length  = None) :
      self.name = name
      self.waveCondition = waveCondition
      self.origin = origin
      self.orientation = [int(i) for i in orientation]  #Assume (1,0,0) ...
      self.length = length
      self.relax = relax  #If false no relaxation (jsut used for boundary condition)

   def pyFoamDict(self) :
      d = {}
      for key, val in self.waveCondition.pyFoamDict().items() :
         d[key] = val

      if self.relax :
         d["relaxationZone"] = { "zoneName"    : self.name + "Zone" ,
                                 "origin"      : Vector( *self.origin ) ,
                                 "orientation" : Vector( *self.orientation ) ,
                                 "zoneName"    :     self.name + "Zone" ,
                               }
      return d


   def zoneBatch(self) :
      """
      Return bat to define the relaxation zone in the mesh
      TODO : More general/elegant way. (right now, works only for x or y direction)
      """
      bigVal = 1e6
      c1 = [ -bigVal, -bigVal , -bigVal ] ;  c2 = [ +bigVal, +bigVal , +bigVal ]
      if -1 in self.orientation :
         i = self.orientation.index(-1)
         c1[ i ] = self.origin[ i ] - self.length
      elif +1 in self.orientation :
         i = self.orientation.index(+1)
         c2[ i ] = self.origin[ i ] + self.length
      return "cellSet {}Zone      new boxToCell ({} {} {}) ({} {} {});".format(self.name ,  *(c1+c2) )
      

class WaveCondition(object) :
   """
      Wave condition for foamStar run
   """
   def __init__( self , waveType = "stokes5th" , height = 1.0 , period = 6.0 , U0 = 0.0 , depth = 60. , rampTime = 0.0 , startTime = 0. , refDirection = [1,0,0]):
      self.waveType = waveType
      self.height  = height
      self.period  = period
      self.U0      = U0
      self.depth   = depth
      self.rampTime  = rampTime   # Ramp time at wave start-up
      self.startTime = startTime
      self.refDirection = refDirection


   def pyFoamDict(self) :
      d =  { "waveType" : self.waveType ,
             "height"   : self.height ,
             "period"   : self.period ,
             "U0"       : Vector( self.U0 ,0.,0. ),
             "depth"    : self.depth ,
             "rampTime" : self.rampTime ,
             "refTime" : self.startTime ,
             "startTime" : self.startTime ,
             "waveType" : self.waveType ,
             "iterateDistance" : "yes" ,
             "refDirection" :  Vector( *self.refDirection)
           }
      return d


class WaveProperties( WriteParameterFile ) :
   """
   waveProperty foamStar file
   """

   def __init__(self , name , initWaveCondition , waveCondition , relaxZones = [] ) :
      WriteParameterFile.__init__(self , name)
      self.relaxZones = relaxZones
      self ["#inputMode"]=  "overwrite";
      self["relaxationNames"] =  [ relax.name for relax in relaxZones if relax.relax]
      self["initCoeffs"] =  initWaveCondition.pyFoamDict()
      #self["initCoeffs"]["relaxZones"] = [ relax.name for relax in relaxZones if relax.relax]
      for relax  in relaxZones :
         self[relax.name + "Coeffs"] = relax.pyFoamDict()
         

   def writeBlendingZoneBatch(self , filename = None ) :
      if filename is None  :  filename = os.path.realpath( os.path.join( os.path.dirname(self.name) , ".." , "blendingZone.batch") )
      with  open( filename , "w" ) as f :
         for relax in self.relaxZones :
            f.write( relax.zoneBatch() + "\n" )
         f.write("quit")
      return

if __name__ == "__main__" :

   """
     Example of use
   """

   waveCond    = WaveCondition( height = 1.0 , period = 6.0  , refDirection = [1.,0.,0.], U0 = 0.0 , depth = 60 )
   relaxInlet  = RelaxZone( name = "inlet"   , relax = True , waveCondition = waveCond , origin = [0 , 0, 0] , orientation = [1 , 0, 0 ] , length = 50. )
   relaxOutlet = RelaxZone( name = "outlet"  , relax = True , waveCondition = waveCond , origin = [600, 0, 0] , orientation = [-1 , 0, 0 ], length = 50 )
   
   waveProperties = WaveProperties("waveProperties" , waveCond  ,  waveCond, relaxZones = (relaxInlet , relaxOutlet))

   print waveProperties
   print relaxInlet.zoneBatch()
   print relaxOutlet.zoneBatch()

   #to write the file
   waveProperties.writeFile()




