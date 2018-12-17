from os.path import join
from ideFoam.inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import Vector, SymmTensor


"""
  Convenience class to simply write "SixDofDomainBody"
"""


class SixDofDomainBody(ReadWriteFile):
    """SixDofDomainBody dictionnary
    """
    
    @classmethod
    def Build(cls , case, mass, inertia, COG, nModes=0, donName=None) :
        res = cls(  name = join(case, getFilePath("sixDofDomainBody") ), read = False )
        
        res.header["class"] = "dictionary"
        res.header["object"] = "singleBody"
    
        res["mass"] = mass
        res["momentOfInertia"] = SymmTensor( *inertia )
        res["cogInitial"] = Vector( *COG )
        res["Xrel"] = "(0 0 0)"
        res["dotXrel"] = "(0 0 0)"
        res["omega"] = "(0 0 0)"
        res["EulerZYX"] = { "rollPitchYaw" : "(0 0 0)" }
        
        if nModes>0: res["modalData"] = { "readFromFile" : "constant/{}.flex".format(donName) }

        return res

if __name__ == "__main__" : 
    print(SixDofDomainBody.Build("test", 1000., 2000., '1 2 3'))



