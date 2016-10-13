import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector,Field,Dimension
from os.path import join
from compatOF import water, air

"""
  Convenience class to simply write "TransportProperties"
"""


class TransportProperties(WriteParameterFile) :
   """
      TransportProperties dictionnary
   """
   def __init__(self , case, rhoWater = 1000 , nuWater = 1e-6, rhoAir = 1. , nuAir = 1.48e-05, sigma = 0.0 ,  version = "foamStar") :
      WriteParameterFile.__init__(self,  name = join(case, "constant" , "transportProperties" )  )
      self.header["class"] = "dictionary"

      if version == "foamStar" : self["phases"] = ["water" , "air"]

      self['"'+water[version]+'"'] = {
                           "transportModel" :   "Newtonian",
                           "nu"             :   "nu [0 2 -1 0 0 0 0] {}".format(nuWater),
                           "rho"             :  "rho [1 -3 0 0 0 0 0] {}".format(rhoWater)
                         }

      self['"'+air[version]+'"'] = {
                                  "transportModel" : "Newtonian",
                                  "nu"             : "nu [0 2 -1 0 0 0 0] {}".format(nuAir),
                                  "rho"            : "rho [1 -3 0 0 0 0 0] {}".format(rhoAir)
                               }

      self[r"sigma"] = "sigma [1 0 -2 0 0 0 0] {}".format(sigma)

if __name__ == "__main__" : 
   print  TransportProperties("test" , version = "foamExtend")



