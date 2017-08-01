import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector, SymmTensor
from os.path import join

"""
  Convenience class to simply write "SixDofDomainBody"
"""


class SixDofDomainBody(WriteParameterFile):
    """
        SixDofDomainBody dictionnary
    """
    def __init__(self , case, mass, inertia, COG, nModes=0, donName=None, version="foamStar") :
        WriteParameterFile.__init__(self,  name = join(case, "0", "uniform", "sixDofDomainBody" ) )
        self.header["class"] = "dictionary"
        self.header["object"] = "singleBody"
    
        self["mass"] = mass
        self["momentOfInertia"] = SymmTensor( *inertia )
        self["cogInitial"] = Vector( *COG )
        self["Xrel"] = "(0 0 0)"
        self["dotXrel"] = "(0 0 0)"
        self["omega"] = "(0 0 0)"
        self["EulerZYX"] = { "rollPitchYaw" : "(0 0 0)" }
        
        if nModes>0: self["modalData"] = { "readFromFile" : "constant/{}.flex".format(donName) }


if __name__ == "__main__" : 
    print  SixDofDomainBody("test", 1000., 2000., '1 2 3', version = "foamExtend")



