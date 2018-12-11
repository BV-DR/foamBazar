import PyFoam
from inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import Vector, DictProxy
import numpy as np
import os
from inputFiles.compatOF import waveTypeDict, foamStarPatch
from math import pi 
# from Spectral.omega2wn import omega2omegae

"""

   Routines to create the waveProperties input file (+blendingZone setSet input)
   
   TODO : Extend to irregular sea-states

"""


class RelaxZone(object) :
    """
    Relaxation zone definition
    """
    def __init__(self, name, relax, waveCondition, origin=None, orientation=None, length=None, bound=None, patchNames=None) :
        """
        name : name of the zone
        waveCondition : WaveCondition object (define height, period...)
        origin : Relaxation zone starting point
        orientation : Orientation of the relaxation zone
        length : Length of the relaxation zone
        bound : End of the relaxation zone (either length or bound should be provided) 
        patchNmes : Name of the patch assiociated to relaxation zone (Old openFoam syntax. Does not work with foamStar 5).
        """
        self.name = name
        self.waveCondition = waveCondition
        self.origin = origin
        self.orientation = [int(i) for i in orientation]  #Assume (1,0,0) ...
        self.length = length
        self.bound = bound
        self.relax = relax  #If false no relaxation (just used for boundary condition)
        self.patchNames = patchNames
    
    def pyFoamDict(self, version = "foamStar"):
        d = DictProxy()
        if version == "foamStar" :
            # for key, val in self.waveCondition.pyFoamDict(version = version).items():
                # d[key] = val
            d["$mycase"] = ""
        else :
            d["waveTheoryName"] = waveTypeDict[self.waveCondition.waveType][version]
        
        if self.relax :
            r = DictProxy()
            if self.patchNames is None :
                r["zoneName"]         = self.name + "Zone"
                # r["relaxationScheme"] = "spatial"
                r["origin"]           = Vector( *self.origin )
                r["orientation"]      = Vector( *self.orientation )
                # r["relaxationShape"]  = "rectangular"
            else :
                r["relaxationScheme"]   = "farfield"
                r["zoneName"]           = self.name + "Zone"
                r["farfieldDistance"]   = self.length
                r["blendingDistance"]   = self.length*0.95
                r["farfieldPatchNames"] = self.patchNames
            d["relaxationZone"] = r
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
                if self.length is not None: c1[ i ] = self.origin[ i ] - self.length
                elif self.bound is not None: c1[ i ] = self.bound
                else: print('ERROR : Please provide at least length of bound for relaxZone')
            elif +1 in self.orientation :
                i = self.orientation.index(+1)
                if self.length is not None: c2[ i ] = self.origin[ i ] + self.length
                elif self.bound is not None: c2[ i ] = self.bound
                else: print('ERROR : Please provide at least length of bound for relaxZone')
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
        d = DictProxy()
        d["waveType"]  = waveTypeDict[self.waveType][version]
        if self.waveType != "noWaves":
            d["height"]    = self.height
            d["period"]    = self.period
            d["depth"]     = self.depth
            d["refTime"]   = self.startTime
            d["startTime"] = self.startTime
        d["rampTime"]  = self.rampTime
        
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
            # d["period"] = 2*pi / omega2omegae(2*pi / self.period , v = self.U0, beta = 180.)
    
        elif version == "foamStar" :
            if self.waveType != "noWaves": d["refDirection"] =  Vector( *self.refDirection)
            d["U0"] = Vector( self.U0 ,0.,0. )
            d["EulerianCurrent"] =  0.0
        else : 
            raise(Exception("Version not known {}".format(version)))
    
        return d


class WaveProperties( ReadWriteFile ) :
    """
    waveProperty foamStar file
    """
    
    @classmethod
    def Build(cls, case, initWaveCondition, relaxZones=[], seaLevel=0., version="foamStar") :
    
        res = cls( name = join(case, getFilePath("waveProperties") ), read = False )
        
        res.case = case
        res.relaxZones = relaxZones
        res ["#inputMode"]=  "overwrite";
    
        if version == "foamStar" :
            res["mycase"] = initWaveCondition.pyFoamDict(version = version)
            res["seaLevel"]=  0.0
            res["initCoeffs"] = {"$mycase" : ''}
            if len(relaxZones)>0:
                res["relaxationNames"] =  [ relax.name for relax in relaxZones if relax.relax]
                for relax  in relaxZones :
                    if version=="foamStar": 
                        relaxName = '"({}|{})Coeffs"'.format(foamStarPatch[relax.name],relax.name)
                    else:
                        relaxName = "{}Coeff".foamat(relax.name)
                    res[relaxName] = relax.pyFoamDict(version = version)
    
        elif version == "foamExtend" :
            res["relaxationNames"] =  [ relax.name for relax in relaxZones if relax.relax]
            res ["seaLevel"]=  0.0;
            res["initWaveTheory"] = waveTypeDict[initWaveCondition.waveType][version]
            #Write the different wave condition
            uniqueWave = set([ initWaveCondition ] +  [i.waveCondition for i in relaxZones ])
            res["waveTheories"] = [waveTypeDict[i.waveType][version] for i in uniqueWave]
            for w in uniqueWave:
                res[ waveTypeDict[w.waveType][version] + "Coeffs"] =  w.pyFoamDict(version = version)
            for relaxZone in relaxZones:
                res[relaxZone.name + "Coeffs"] = relaxZone.pyFoamDict(version = version)
    
        elif version == "swenseFoam" :
            res["relaxationName"] =  "farfield"
            res["limitMagUInc"] = 3
            res["seaLevel"]=  0.0
            res["farfieldCoeffs"] = relaxZones[0].waveCondition.pyFoamDict(version = version)
            res["farfieldCoeffs"] = dict(res["farfieldCoeffs"] , **relaxZones[0].pyFoamDict(version = version) )
            
        return res


    def writeBlendingZoneBatch(self , filename = None ) :
        if filename is None  :  filename = os.path.join( self.case , "blendingZone.batch")
        with open( filename , "w" ) as f :
            for relax in self.relaxZones :
                f.write( relax.zoneBatch() + "\n" )
            # f.write("quit")
        return
    
    def writeFile(self, *args, **kwargs):
        ReadWriteFile.writeFile(self, *args, **kwargs)
        try : 
            self.writeBlendingZoneBatch()
        except AttributeError : 
            print("ReParsing relaxation zone not possible yet")
        
        

noWaves = WaveCondition( waveType = "noWaves" , U0 = 0.0)        
        
if __name__ == "__main__" :
    """Example of use
    """

    waveCond    = WaveCondition( height = 1.0 , period = 6.0  , refDirection = [1.,0.,0.], U0 = 0.0 , depth = 60 )
    relaxInlet  = RelaxZone( name = "inlet"   , relax = True , waveCondition = waveCond , origin = [0 , 0, 0] , orientation = [1 , 0, 0 ] , length = 50. , patchNames = ["inlet", "outlet"] )
    #relaxOutlet = RelaxZone( name = "outlet"  , relax = True , waveCondition = waveCond , origin = [600, 0, 0] , orientation = [-1 , 0, 0 ], length = 50 )
    relaxOutlet = RelaxZone( name = "outlet"  , relax = True , waveCondition = noWaves , origin = [600, 0, 0] , orientation = [-1 , 0, 0 ], length = 50 )
    
    waveProperties = WaveProperties.Build("test" , waveCond  ,  relaxZones = (relaxInlet , relaxOutlet) , version = "foamStar" )
    
    print(waveProperties)
    
    #to write the file
    #waveProperties.writeFile()




