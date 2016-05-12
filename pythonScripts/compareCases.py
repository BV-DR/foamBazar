#!/usr/bin/python
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import DictProxy
import os
from os.path  import join

def compareFile(file1 ,file2, name1 = "1" , name2 = "2" , valName = "") :
   f1 = ParsedParameterFile(file1).getValueDict()
   f2 = ParsedParameterFile(file2).getValueDict()
   comparePyFoam(f1 , f2 , name1 = name1, name2 = name2 , valName = valName)

def compareCase( dir1 , dir2 , version = "s"  ) :

   fileList = [
               r"0/org/pd" ,
               r"0/org/alpha1" ,
               r"0/org/U" ,
               r"constant/waveProperties" ,
               r"constant/polyMesh/blockMeshDict" ,
               r"system/fvSchemes" ,
               r"system/fvSolution" ,
               r"system/controlDict" ,
               ]

   print "Comparing {} with {}".format(dir1 , dir2)
   for f in fileList :
      print "!------- {} -------!".format(os.path.basename(f))
      compareFile( join(dir1,f) , join(dir2,f) , name1 = "1" , name2 = "2", valName = os.path.basename(f) )
      print "!------------------!\n"


def comparePyFoam(d1 , d2 , name1 = "1" , name2 = "2", valName = ""):

   if type(d1) == type(d2) == DictProxy or type(d1) == type(d2) == dict:
      s1 = set(d1.keys())
      s2 = set(d2.keys())
      for k in s1.intersection(s2) :
         comparePyFoam( d1[k] , d2[k] , name1 , name2,  valName + "/" + k )
      for k in s1.difference(s2) :
         print "{} {} not in {} and in {}".format(valName, k,name1, name2)
      for k in s2.difference(s1) :
         print "{} {} not in {} and in {}".format(valName, k,name1, name2)

   elif type(d1)== type(d2) == list :
      if len(d1) == len(d2) :
         for i in range(len(d1)) :
            comparePyFoam(d1[i] , d2[i] , name1 = "1" , name2 = "2", valName = valName + "/" + str(i+1))
      else :
         print "list {} is different {} != {}".format( valName, d1 , d2 )


   else :
      if d1.__repr__() != d2.__repr__() :
         print "value {} is different {} != {}".format( valName, d1 , d2 )


if __name__ == "__main__" :

   import argparse
   parser = argparse.ArgumentParser(description='FoamStar log parser')
   parser.add_argument( "file1" )
   parser.add_argument( "file2" )
   args = parser.parse_args()
   
   compareCase(args.file1 ,args.file2)

   #file1 = r"\\10.67.24.192\bigr\openFoam\2Dwave\RunFine\CrankNicolson_0.9_Speed_swenseFoam0.0"
   #file2 = r"\\10.67.24.192\bigr\openFoam\2Dwave\Run\swenseFoam_g1_Speed_0.0"

   



