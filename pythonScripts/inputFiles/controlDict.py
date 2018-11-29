from os.path import join
from inputFiles import ReadWriteFile
from PyFoam.Basics.DataStructures import DictProxy
from inputFiles.waveProbes import setWaveProbes
from inputFiles.compatOF import application, surfaceElevation

class ControlDict( ReadWriteFile ) :

    @classmethod
    def Build(cls, case, deltaT, endTime, startFrom="latestTime", startTime=0., autoDeltaTCo=None, writeControl="timeStep", 
              writeInterval=1, purgeWrite=0, writePrecision=7, writeCompression="compressed", runTimeModifiable="no",
              writeProbesInterval=None,  waveProbesList=None, adjustTimeStep=None, outputInterval=1,
              outputMotions=False, vbmPatch=None, forcesPatch=None, pressuresPatch=None, outputLocalMotions=False, rhoWater = 1000,
              OFversion=5, version="foamStar"):

        res = cls(join(case , "system" , "controlDict" ), read = False)
        res["application"]       = application[version]
        res["startFrom"]         = startFrom
        res["startTime"]         = startTime
        res["stopAt"]            = "endTime"
        res["endTime"]           = endTime
        res["deltaT"]            = deltaT
        res["writeControl"]      = writeControl
        res["writeInterval"]     = writeInterval
        res["purgeWrite"]        = purgeWrite
        res["writeFormat"]       = "ascii"
        res["writePrecision"]    = writePrecision
        res["writeCompression"]  = writeCompression
        res["timeFormat"]        = "general"
        res["timePrecision"]     = 6
        res["runTimeModifiable"] = runTimeModifiable
        if OFversion==5: res["fileHandler"] = "collated"
      
      
        if adjustTimeStep is not None :
            res["adjustTimeStep"] = "yes"
            res["maxCo"]          = adjustTimeStep[0]
            res["maxAlphaCo"]     = adjustTimeStep[1]
            res["maxDeltaT"]      = adjustTimeStep[2]
        else :
            res["adjustTimeStep"] = "no"
            res["maxCo"]          = 0.5
            res["maxAlphaCo"]     = 0.5
            res["maxDeltaT"]      = 1.
      
        if version == "foamStar" :
            res ["libs"] =  ['"libfoamStar.so"' , '"libBVtabulated6DoFMotion.so"', ]
        elif version == "snappyHexMesh" :
            pass
        else :
            res ["libs"] =  [ '"libforces.so"' , ]
        
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
            if OFversion==5:
                vbmDict["writeControl"]      = "timeStep"
                vbmDict["writeInterval"]     = 1
            else:
                vbmDict["outputControl"]      = "timeStep"
                vbmDict["outputInterval"]     = 1
            vbmDict["hullPatches"]        = '({})'.format(vbmPatch[0])
            vbmDict["donFileName"]        = '"{}.don"'.format(vbmPatch[1])
            vbmDict["wldFileName"]        = '""'
            vbmDict["log"]                = "true"
            vbmDict["ySymmetry"]          = "yes"
            fDict["vbm"] = vbmDict
        
        #Get forces on patch list
        if forcesPatch is not None :
            for fpatch in forcesPatch:
                forcesDict = DictProxy()
                forcesDict["type"]               = "forces"
                forcesDict["functionObjectLibs"] = ['"libforces.so"']
                forcesDict["patches"]            = [fpatch]
                if OFversion==5:
                    forcesDict["rhoInf"]         = rhoWater
                    forcesDict["rhoName"]        = "rho"
                    forcesDict["pName"]          = "p"
                    forcesDict["UName"]          = "U"
                    forcesDict["writeControl"]   = "timeStep"
                    forcesDict["writeInterval"]  = outputInterval
                else:
                    forcesDict["rhoInf"]         = rhoWater
                    forcesDict["rhoName"]        = "rho"
                    forcesDict["pName"]          = "p"
                    forcesDict["UName"]          = "U"
                    forcesDict["outputControl"]  = "timeStep"
                    forcesDict["outputInterval"] = outputInterval
                forcesDict["CofR"]               = "( 0 0 0 )"
                forcesDict["log"]                = True
                fDict["forces_"+fpatch] = forcesDict
            
        #Get forces on patch list
        if pressuresPatch is not None :
            pressuresDict = DictProxy()
            if OFversion==5:
                pressuresDict["type"]               = "surfaceFieldValue"
                pressuresDict["libs"]               = ['"libfieldFunctionObjects.so"']
                pressuresDict["regionType"]         = "patch"
                pressuresDict["name"]               = pressuresPatch[0]
                pressuresDict["surfaceFormat"]      = "foam"
                pressuresDict["writeControl"]       = "timeStep"
                pressuresDict["writeInterval"]      = outputInterval
                pressuresDict["writeFields"]        = True
            else:
                pressuresDict["type"]               = "faceSource"
                pressuresDict["functionObjectLibs"] = ['"libfieldFunctionObjects.so"']
                pressuresDict["valueOutput"]        = "true"
                pressuresDict["source"]             = "patch"
                pressuresDict["sourceName"]         = pressuresPatch[0]
                pressuresDict["surfaceFormat"]      = "foamFile"
                pressuresDict["outputControl"]      = "timeStep"
                pressuresDict["outputInterval"]     = outputInterval            

            pressuresDict["operation"]          = "none"
            pressuresDict["log"]                = True
            pressuresDict["fields"]             = ['p','p_rgh']
            fDict["pressures"] = pressuresDict
            
        # localMotion
        if outputLocalMotions:
            localDict = DictProxy()
            localDict["type"]       = "localMotion"
            localDict["motionData"] = "sixDofDomainBody"
            localDict["flxPTS"]     = "localMotion"
            fDict["localMotion"] = localDict
        
        #Construct waveProbe dict from waveProbes list  [ ( x,y,z_min,nb_point ) , ... ]
        if waveProbesList is not None:
            fDict[surfaceElevation[version]] = setWaveProbes( waveProbesList , version = version , writeProbesInterval = writeProbesInterval )
         
        if len(fDict)>0:
            res["functions"] = fDict
        return res


if __name__ == "__main__" :
   print(ControlDict.Build( "test" , deltaT = 0.10, endTime = 60., waveProbesList = ( (10.,0.,-1.,+1 , 100) , (15.,0.,-1.,+1 , 100) ) , forcesPatch = ["test"] ))
