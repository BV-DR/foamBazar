from ideFoam.inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import  DictProxy
from os.path import join
from ideFoam.inputFiles.compatOF import water, air

"""
  Convenience class to simply write "TransportProperties"
"""


class TransportProperties(ReadWriteFile) :
    """
        TransportProperties dictionnary
    """
    @classmethod
    def Build(cls , case, rhoWater = 1000 , nuWater = 1e-6, rhoAir = 1. , nuAir = 1.48e-05, sigma = 0.0 ,  application = "foamStar") :
        
        res = cls(  name = join(case, getFilePath("transportProperties") ), read = False )
        
        res.header["class"] = "dictionary"

        if application == "foamStar" : res["phases"] = ["water" , "air"]
        
        dw = DictProxy()
        dw["transportModel"] =  "Newtonian"
        dw["nu"]             =  "nu [0 2 -1 0 0 0 0] {}".format(nuWater)
        dw["rho"]            =  "rho [1 -3 0 0 0 0 0] {}".format(rhoWater)
        res['"'+water[application]+'"'] = dw

        da = DictProxy() 
        da["transportModel"] =  "Newtonian",
        da["nu"]             =  "nu [0 2 -1 0 0 0 0] {}".format(nuAir)
        da["rho"]            =  "rho [1 -3 0 0 0 0 0] {}".format(rhoAir)
        res['"'+air[application]+'"'] = da

        res[r"sigma"] = "sigma [1 0 -2 0 0 0 0] {}".format(sigma)
        
        return res

if __name__ == "__main__" : 
    print(TransportProperties.Build("test" , application = "foamExtend"))



