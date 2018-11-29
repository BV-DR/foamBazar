import os
import numpy as np
from PyFoam.Basics.DataStructures import DictProxy, Vector, Field, TupleProxy
from os.path import join
from inputFiles import ReadWriteFile

"""
  Convenience class to simply write "DynamicMeshDict"
"""



template_free = """
FoamFile
{
    version     5.0;
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}

dynamicFvMesh       multiBodyFvMesh;

multiBodyFvMeshCoeffs
{
    meshMotion TOFILL;
    
    solver          RKCK45;
    absTol          1e-8;
    relTol          0;

    rampTime        -1;
    releaseTime     -1;

    ship
    {
	type rigidBody;
    
    innerDistance       TOFILL;  // For cpMorphing only
    outerDistance       TOFILL;  // For cpMorphing only
    
    motionPatches TOFILL;
    
	mass                TOFILL;
	momentOfInertia     (-1 -1 -1 -1 -1 -1);
	CoRInitial          (-1 -1 -1); 

    	loads
    	{
        	gravity { type gravity; value (0 0 -9.81); }
        	fluid { type fluidForce; patches TOFILL; dynRelax no; relaxCoeff 0.5; }
    	}
    	constraints
    	{
        	heave_pitch_only
        	{
		    	type        sixDOF;
		    	default     fixAll;         // fixAll | freeAll
		    	except      (heave pitch);  // overwrite default
        	}
    	}
    }

}
"""



class DynamicMeshDict( ReadWriteFile ):
    
   
    @classmethod
    def Build_static( cls, case ):
        res = cls( join(case, "constant", "dynamicMeshDict" ), read = False )
        res.header["class"] = "dictionary"
        res["dynamicFvMesh"] = "staticFvMesh"
        return res
     
    @classmethod
    def Build_free( cls, case, mass, cog, inertia, rampTime, releaseTime, hullPatch="(ship)", meshMotion = "cpMorphing", innerDistance=None, outerDistance=None, OFversion=5, version="foamStar"):
        res = cls( join(case, "constant", "dynamicMeshDict" ), read = False )
        res.header["class"] = "dictionary"
        
        import tempfile
        with tempfile.NamedTemporaryFile("w", delete = False) as f :
            f.write( template_free)
        res.content = ReadWriteFile(  f.name  ).content
        os.unlink(f.name)
         
        res["multiBodyFvMeshCoeffs"]["meshMotion"] = meshMotion
        res.setMechanics(mass, cog, inertia,hullPatch)
        
        if meshMotion == "cpMorphing":
            res["multiBodyFvMeshCoeffs"]["ship"]["innerDistance"] = innerDistance
            res["multiBodyFvMeshCoeffs"]["ship"]["outerDistance"] = outerDistance
        else : 
            raise(Exception(meshMotion + " Not implented yet (in scripts)"))
            
        res["multiBodyFvMeshCoeffs"]["rampTime"] = rampTime
        res["multiBodyFvMeshCoeffs"]["releaseTime"] = releaseTime
        
        return res
            
         
    def setMechanics(self, mass, cog , inertia,hullPatch):
        self["multiBodyFvMeshCoeffs"]["ship"]["mass"] = mass
        self["multiBodyFvMeshCoeffs"]["ship"]["CoRInitial"] = Vector(*cog)
        self["multiBodyFvMeshCoeffs"]["ship"]["momentOfInertia"] = tuple(inertia).__str__().replace(",", " ")
        self["multiBodyFvMeshCoeffs"]["ship"]["motionPatches"] = hullPatch
        self["multiBodyFvMeshCoeffs"]["ship"]["loads"]["fluid"]["patches"] = hullPatch
       
        
    def parseBuild(self):
        """Return argument necessary to build the object
        """
        res = { "mass" : self["multiBodyFvMeshCoeffs"]["ship"]["mass"],
                "cog" : self["multiBodyFvMeshCoeffs"]["ship"]["CoRInitial"] ,
                "inertia" : self["multiBodyFvMeshCoeffs"]["ship"]["momentOfInertia"] ,
                "hullPatch" : self["multiBodyFvMeshCoeffs"]["ship"]["motionPatches"]
              }
        return res
        
    @classmethod
    def Build_imposed(cls, case,  dispFile='', COG=[0.,0.,0.], OFversion=5, version="foamStar") :
        res = cls( join(case, "constant", "dynamicMeshDict" ), read = False )
        if OFversion==5:
            res["dynamicFvMesh"] = "dynamicMotionSolverFvMesh"
            res["motionSolver"] = "solidBody"
            res["solidBodyMotionFunction"] = "tabulated6DoFMotion"
            tab = DictProxy()
            tab["timeDataFileName"] = '"'+dispFile+'"'
            tab["CofG"] = '( {:.6f} {:.6f} {:.6f} )'.format(*COG)
            res["tabulated6DoFMotionCoeffs"] = tab
        else:
            res["dynamicFvMesh"] = "solidBodyMotionFvMesh"
            sdc = DictProxy()
            sdc["solidBodyMotionFunction"] = "BVtabulated6DoFMotion"
            tab = DictProxy()
            tab["timeDataFileName"] = '"'+dispFile+'"'
            tab["CofG"] = '( {:.6f} {:.6f} {:.6f} )'.format(*COG)
            sdc["BVtabulated6DoFMotionCoeffs"] = tab
            res["solidBodyMotionFvMeshCoeffs"] = sdc
        return res


if __name__ == "__main__" :
    #print(DynamicMeshDict.Build_free("test", inertia = [1,1,1,1,1,1], cog = [1,1,1], mass = 123, releaseTime = 0.0, rampTime = 0.0, hullPatch='ship',  version = "foamExtend"))
    
    a = ReadWriteFile(r"/mnt/d/Etudes/Model_test_ECN/VBM_MODEL_TEST/OFRun/test/constant/dynamicMeshDict", read = True)
    b = DynamicMeshDict(r"/mnt/d/Etudes/Model_test_ECN/VBM_MODEL_TEST/OFRun/test/constant/dynamicMeshDict", read = True)
