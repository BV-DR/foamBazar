from os.path import join
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector
from waveProbes import setWaveProbes
from compatOF import application

class ControlDict( WriteParameterFile ) :

   def __init__(self , case ,  startTime =0. , endTime = 60. , deltaT =  0.1 , writeInterval = 0.1 , writeProbesInterval = None,  waveProbesList = None, adjustTimeStep = None, forcesPatch = None , version = "foamStar"):

      WriteParameterFile.__init__(self , join(case , "system" , "controlDict" ))
      self ["startTime"] =  startTime
      self ["deltaT"]    =  deltaT
      self ["writeInterval"]    =  writeInterval
      self ["endTime"]    =  endTime

      if not adjustTimeStep is None :
         self["adjustTimeStep"] = "yes"
         self["maxCo"]          = adjustTimeStep[0]
         self["maxAlphaCo"]     = adjustTimeStep[1]
         self["maxDeltaT"]      = adjustTimeStep[2]
      else :
         self["adjustTimeStep"] = "no"
         self["maxCo"]          = 0.  #Unused
         self["maxAlphaCo"]     = 0.  #Unused
         self["maxDeltaT"]      = 0.  #Unused

      #Constant parameter
      self ["application"]=  application[version]
      self ["startFrom"]         =  "latestTime"
      self ["stopAt"]            =  "endTime"
      self ["writeControl"]      =  "adjustableRunTime"
      self ["purgeWrite"]        =  0
      self ["writePrecision"]    =  7
      self ["writeCompression"]  =  "compressed"
      self ["timeFormat"]        =  "general"
      self ["timePrecision"]     =  6
      self ["runTimeModifiable"] =  "yes"
      
      if version == "foamStar" :
         self ["libs"] =  ['"libCRSwaves.so"' , '"libforces.so"' , ]
      else :
         self ["libs"] =  [ '"libforces.so"' , ]

      #Construct waveProbe dict from waveProbes list  [ ( x,y,z_min,nb_point ) , ... ]
      if waveProbesList is not None  :
         self["functions"] = setWaveProbes( waveProbesList , version = version , writeProbesInterval = writeProbesInterval )

      #Get forces on patch list
      if forcesPatch is not None :
         self["functions"]["forcesTotal"] = {
                                               "type" :  "forces" ,
                                               #"functionObjectLibs" : ,
                                               "patches" :  forcesPatch ,
                                               "rhoInf"  :  1000 ,
                                               "rhoName" :  "rho" ,
                                               "pName"   :   "p" ,
                                               "UName"   :  "U"   ,
                                               "log"     :   True  ,
                                               "outputControl" : "timeStep" ,
                                               "outputInterval" : 1 ,
                                               "CofR" : Vector( *[0, 0, 0] )
                                           }


if __name__ == "__main__" :

   test = ControlDict( "test" , waveProbesList = ( (10.,0.,-1.,+1 , 100) , (15.,0.,-1.,+1 , 100) ) , forcesPatch = ["test"] )

   print test
