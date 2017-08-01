import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import DictProxy
from os.path import join

"""
  Convenience class to simply write "DynamicMeshDict"
"""


class DynamicMeshDict(WriteParameterFile):
    """
        DynamicMeshDict dictionnary
    """
    def __init__(self , case, static=False, hullPatch=None, addDamping=False, lpp=0, bc=0, version="foamStar") :
        WriteParameterFile.__init__(self,  name = join(case, "constant", "dynamicMeshDict" ) )
        self.header["class"] = "dictionary"
    
        if static:
            self["dynamicFvMesh"] = "staticFvMesh"
        else:
            self["dynamicFvMesh"] = "sixDofDomainFvMesh"
        self["motionSolverLibs"] = ['"libfvMotionSolvers.so"', ]
        self["solver"] = "displacementLaplacian"
        self["displacementLaplacianCoeffs"] = { "diffusivity" : "invDistSqr ({})".format(hullPatch) }
        sdc = DictProxy()
        sdc["hullPatches"] = "({})".format(hullPatch)
        sdc["meshConfigType"] = 4
        sdc["solver"] = "RKCK45"
        sdc["absTol"] = 1e-8
        sdc["relTol"] = 0
        sdc["rampTime"] = 15
        ldc = DictProxy()
        ldc["gravity"] = { "type" : "gravity", "value" : "(0 0 -9.81)" }
        ff = DictProxy()
        ff["type"] = "fluidForce"
        ff["patches"] = "({})".format(hullPatch)
        ff["dynRelax"] = "off"
        ff["relaxCoeff"] = 0.5
        ff["ySym"] = "yes"
        ldc["fluid"] = ff
        if addDamping:
            dmp = DictProxy()
            dmp["type"] = "dampingSinkageAndTrim"
            dmp["heaveDampingCoef"] = 1000
            dmp["pitchDampingCoef"] = 1000
            dmp["Lpp"] = lpp
            dmp["Breadth"] = bc
            ldc["dampingCoeff"] = dmp
        sdc["loads"] = ldc
        hp = DictProxy()
        hp["type"] = "cog6DOF"
        hp["default"] = '"fixAll"'
        hp["except"] = "(heave pitch)"
        sdc["constraints"] = { "heavePitch" : hp }
        self["sixDofDomainFvMeshCoeffs"] = sdc
        
if __name__ == "__main__" : 
    print  DynamicMeshDict("test", hullPatch='ship', addDamping=True, lpp=100., bc=20., version = "foamExtend")
