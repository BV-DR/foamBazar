from os.path import join
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector, DictProxy
from waveProbes import setWaveProbes
from compatOF import application

class ControlDict( WriteParameterFile ) :

    def __init__(self , case ,  startFrom="latestTime", startTime=0 , endTime=60. , deltaT=None , autoDeltaTCo=None , writeInterval=0.1 , purgeWrite=0, writePrecision=7, runTimeModifiable="no", writeProbesInterval=None,  waveProbesList=None, adjustTimeStep=None, outputMotions=False, vbmPatch=None, forcesPatch=None, outputLocalMotions=False, version="foamStar"):

        WriteParameterFile.__init__(self , join(case , "system" , "controlDict" ))
        self["application"]       = application[version]
        self["startFrom"]         = startFrom
        self["startTime"]         = startTime
        self["stopAt"]            = "endTime"
        self["endTime"]           = endTime
        self["deltaT"]            = deltaT
        self["writeControl"]      = "timeStep"
        self["writeInterval"]     = writeInterval
        self["purgeWrite"]        = purgeWrite
        self["writeFormat"]       = "ascii"
        self["writePrecision"]    = writePrecision
        self["writeCompression"]  = "compressed"
        self["timeFormat"]        = "general"
        self["timePrecision"]     = 6
        self["runTimeModifiable"] = runTimeModifiable
      
      
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
      
        if version == "foamStar" :
            self ["libs"] =  ['"libfoamStar.so"' , ]
        else :
            self ["libs"] =  [ '"libforces.so"' , ]
        
        # Set functions
        fDict = DictProxy()

        # Motions
        if outputMotions:
            motionDict = DictProxy()
            motionDict["type"] = "motionInfo"
            fDict["motionInfo"] = motionDict
            
        # Internal loads
        if vbmPatch is not None :
            vbmDict = DictProxy()
            vbmDict["type"]               = "internalLoads"
            vbmDict["outputControl"]      = "timeStep"
            vbmDict["outputInterval"]     = 1
            vbmDict["hullPatches"]        = vbmPatch[0]
            vbmDict["donFileName"]        = '"'+vbmPatch[1]+'.don"'
            vbmDict["wldFileName"]        = '""'
            vbmDict["log"]                = "true"
            vbmDict["ySymmetry"]          = "yes"
            fDict["vbm"] = vbmDict
        
        #Get forces on patch list
        if forcesPatch is not None :
            forcesDict = DictProxy()
            forcesDict["type"]               = "forces"
            forcesDict["functionObjectLibs"] = '"libforces.so"'
            forcesDict["patches"]            = forcesPatch
            forcesDict["rhoInf"]             = 1000
            forcesDict["rhoName"]            = "rho"
            forcesDict["pName"]              = "p"
            forcesDict["UName"]              = "U"
            forcesDict["log"]                = True
            forcesDict["outputControl"]      = "timeStep"
            forcesDict["outputInterval"]     = 1
            forcesDict["CofR"]               = "( 0 0 0 )"
            fDict["forces"] = forcesDict
            
        # localMotion
        if outputLocalMotions:
            localDict = DictProxy()
            localDict["type"]       = "localMotion"
            localDict["motionData"] = "sixDofDomainBody"
            localDict["flxPTS"]     = "localMotion"
            fDict["localMotion"] = localDict
        
        #Construct waveProbe dict from waveProbes list  [ ( x,y,z_min,nb_point ) , ... ]
        # if waveProbesList is not None  :
            # self["functions"] = setWaveProbes( waveProbesList , version = version , writeProbesInterval = writeProbesInterval )
            
        self["functions"] = fDict


if __name__ == "__main__" :

   test = ControlDict( "test" , waveProbesList = ( (10.,0.,-1.,+1 , 100) , (15.,0.,-1.,+1 , 100) ) , forcesPatch = ["test"] )

   print test
