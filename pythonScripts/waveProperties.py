import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector
import numpy as np
import os
from compatOF import waveTypeDict
from math import pi 
from Spectral.omega2wn import omega2omegae

"""

   Routines to create the waveProperties input file (+blendingZone setSet input)
   
   TODO : Extend to irregular sea-states

"""


class RelaxZone(object) :
   """
      Relaxation zone definition
   """
   def __init__(self, name , relax , waveCondition , origin = None , orientation = None, length  = None , patchNames = None) :
      self.name = name
      self.waveCondition = waveCondition
      self.origin = origin
      self.orientation = [int(i) for i in orientation]  #Assume (1,0,0) ...
      self.length = length
      self.relax = relax  #If false no relaxation (jsut used for boundary condition)
      self.patchNames = patchNames

   def pyFoamDict(self, version = "foamStar") :
      d = {}
      if version == "foamStar" :
         for key, val in self.waveCondition.pyFoamDict(version = version).items() :
            d[key] = val
      else :
         d["waveTheoryName"] = waveTypeDict[self.waveCondition.waveType][version]

      if self.relax :
         if self.patchNames is None :
            d["relaxationZone"] = { "zoneName"    : self.name + "Zone" ,
                                    "relaxationScheme" : "spatial" ,
                                    "origin"      : Vector( *self.origin ) ,
                                    "orientation" : Vector( *self.orientation ) ,
                                    "relaxationShape"  : "rectangular"
                                  }
         else :
            d["relaxationZone"] = {
                                    "relaxationScheme" : "farfield"       ,
                                    "zoneName"    : self.name + "Zone"    ,
                                    "farfieldDistance"      : self.length      ,
                                    "blendingDistance"      : self.length*0.95 ,
                                    "farfieldPatchNames" : self.patchNames

                                  }
      return d

   def zoneBatch(self) :
      """
      Return bat to define the relaxation zone in the mesh
      TODO : More general/elegant way. (right now, works only for x or y direction)
      """
      if self.patchNames is None :
         bigVal = 1e6
         c1 = [ -bigVal, -bigVal , -bigVal ] ;  c2 = [ +bigVal, +bigVal , +bigVal ]
         if -1 in self.orientation :
            i = self.orientation.index(-1)
            c1[ i ] = self.origin[ i ] - self.length
         elif +1 in self.orientation :
            i = self.orientation.index(+1)
            c2[ i ] = self.origin[ i ] + self.length
         return "cellSet {}Zone      new boxToCell ({} {} {}) ({} {} {});".format(self.name ,  *(c1+c2) )
      else :
         return ""


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


   def pyFoamDict(self , version = "foamStar") :


      d =  { "waveType" : waveTypeDict[self.waveType ][version] ,
             "height"   : self.height ,
             "period"   : self.period ,
             "depth"    : self.depth ,
             "rampTime" : self.rampTime ,
             "refTime" : self.startTime ,
             "startTime" : self.startTime ,
             "iterateDistance" : "yes" ,
           }

      if self.waveType == "streamFunction" :
         d["order"] = 25

      if version == "foamExtend" or  version == "swenseFoam":
         d["wind"] = Vector( 0. ,0.,0. )
         d["currentType"] = "constantCurrent"
         d["U0"] = Vector( 0, 0. ,0. )
         d["setEulerianCurrent"] =  True
         d["EulerianCurrent"] =  self.U0
         d["waveDirection"] =  Vector( *self.refDirection)
         d["phi"] =  0.
         #In wave2foam, with "EulerianCurrent", input is encounter period
         d["period"] = 2*pi / omega2omegae(2*pi / self.period , v = self.U0, beta = 180.)

      elif version == "foamStar" :
         d["refDirection"] =  Vector( *self.refDirection)
         d["U0"] = Vector( self.U0 ,0.,0. )
         d["EulerianCurrent"] =  0.0
      else : 
         raise(Exception("Version not known {}".format(version)))

      return d


class WaveProperties( WriteParameterFile ) :
   """
   waveProperty foamStar file
   """

   def __init__(self , name , initWaveCondition , relaxZones = [] , version = "foamStar") :
      WriteParameterFile.__init__(self , name)
      self.relaxZones = relaxZones
      self ["#inputMode"]=  "overwrite";

      if version == "foamStar" :
         self["relaxationNames"] =  [ relax.name for relax in relaxZones if relax.relax]
         self["initCoeffs"] =  initWaveCondition.pyFoamDict(version = version)
         for relax  in relaxZones :
            self[relax.name + "Coeffs"] = relax.pyFoamDict(version = version)

      elif version == "foamExtend" :
         self["relaxationNames"] =  [ relax.name for relax in relaxZones if relax.relax]
         self ["seaLevel"]=  0.0;
         self["initWaveTheory"] = waveTypeDict[initWaveCondition.waveType][version]
         #Write the different wave condition
         uniqueWave = set([ initWaveCondition ] +  [i.waveCondition for i in relaxZones ])
         self["waveTheories"] = [waveTypeDict[i.waveType][version] for i in uniqueWave]
         for w in uniqueWave :
            self[ waveTypeDict[w.waveType][version] + "Coeffs"] =  w.pyFoamDict(version = version)
         for relaxZone in relaxZones :
           self[relaxZone.name + "Coeffs"] = relaxZone.pyFoamDict(version = version)

      elif version == "swenseFoam" :
         self["relaxationName"] =  "farfield"
         self["limitMagUInc"] = 3
         self["seaLevel"]=  0.0
         self["farfieldCoeffs"] = relaxZones[0].waveCondition.pyFoamDict(version = version)
         self["farfieldCoeffs"] = dict(self["farfieldCoeffs"] , **relaxZones[0].pyFoamDict(version = version) )


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

   noWaveCond  = WaveCondition( waveType = "noWaves" , U0 = 0.0)
   waveCond    = WaveCondition( height = 1.0 , period = 6.0  , refDirection = [1.,0.,0.], U0 = 0.0 , depth = 60 )
   relaxInlet  = RelaxZone( name = "inlet"   , relax = True , waveCondition = waveCond , origin = [0 , 0, 0] , orientation = [1 , 0, 0 ] , length = 50. , patchNames = ["inlet", "outlet"] )
   #relaxOutlet = RelaxZone( name = "outlet"  , relax = True , waveCondition = waveCond , origin = [600, 0, 0] , orientation = [-1 , 0, 0 ], length = 50 )
   relaxOutlet = RelaxZone( name = "outlet"  , relax = True , waveCondition = noWaveCond , origin = [600, 0, 0] , orientation = [-1 , 0, 0 ], length = 50 )

   waveProperties = WaveProperties("waveProperties" , waveCond  ,  relaxZones = (relaxInlet , relaxOutlet) , version = "swenseFoam" )

   print waveProperties



   #print relaxInlet.zoneBatch()
   #print relaxOutlet.zoneBatch()

   #to write the file
   waveProperties.writeFile()




