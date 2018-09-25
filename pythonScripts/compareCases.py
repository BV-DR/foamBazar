#!/usr/bin/python
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import DictProxy
import os
from os.path  import join

import copy

def compareFile(file1 ,file2, name1 = "1" , name2 = "2" , valName = "") :
   f1 = ParsedParameterFile(file1)
   f2 = ParsedParameterFile(file2)

   comparePyFoam(f1.content , f2.content , name1 = name1, name2 = name2 , valName = valName)

def compareCase( dir1 , dir2 , version = "s"  ) :

   fileList = [
                r"0/org/pd" ,
                r"0/org/alpha1" ,
                r"0/org/alpha.water" ,
                r"0/org/U" ,
                r"0/org/UDiff" ,
                r"0/org/UInc" ,
                r"0/org/LevelSetDiff" ,
                r"0/org/pd" ,
                r"0/org/p_rgh" ,
                
                r"0/pd" ,
                r"0/alpha1" ,
                r"0/alpha.water" ,
                r"0/U" ,
                r"0/UDiff" ,
                r"0/UInc" ,
                r"0/LevelSetDiff" ,
                r"0/pd" ,
                r"0/p_rgh" ,

                r"constant/RASProperties" ,
                r"constant/waveProperties" ,
                r"constant/polyMesh/blockMeshDict" ,
                r"constant/turbulenceProperties" ,
                r"constant/g" ,
		r"constant/wavesProperties", 
		r"constant/dynamicMeshDict",
                r"constant/transportProperties" ,
                r"system/fvSchemes" ,
                r"system/decomposeParDict" ,
                r"system/fvSolution" ,
                r"system/controlDict" ,
               ]

   #print "Comparing {} with {}".format(dir1 , dir2)
   for f in fileList :
      f1 =  join(dir1,f)
      f2 =  join(dir2,f)
      if os.path.exists(f1) and os.path.exists(f2) :
         print("!------- {:^40} -------!".format(os.path.basename(f)))
         compareFile( f1 , f2 , name1 = "1" , name2 = "2", valName = os.path.basename(f) )
         print("!--------------------------------------------------------!\n")
      elif not os.path.exists(f1) and os.path.exists(f2) :
         print("{} not in {}".format( f, dir1 ))
      elif not os.path.exists(f2) and os.path.exists(f1) :
         print("{} not in {}".format( f, dir2 ))



def dictProxyAllKeys(dictProxy):
   """
      Return all dict proxy available keys (including key1 and key2 when "key1|key2")
   """
   allKeys = [f for f in dictProxy.keys()]
   for a in dictProxy._regex :
      allKeys.append(a[0][1:-1])    # Remove the double quote, so that dictProxy[ regexPatern ] works
   return allKeys



def comparePyFoam(d1 , d2 , name1 = "1" , name2 = "2", valName = ""):
   #print "Comp" , valName

   if type(d1) == type(d2) == DictProxy or type(d1) == type(d2) == dict:
      s1 = set(d1.keys())
      s2 = set(d2.keys())
      for k in s1.intersection(s2) :
         comparePyFoam( d1[k] , d2[k] , name1 , name2,  valName + "/" + str(k) )
      for k in s1.difference(s2) :
         print("{:48s} {:25s} {:20s} != {:20s}".format(valName, str(k), str(d1[k]) , "-"))
      for k in s2.difference(s1) :
         print("{:48s} {:25s} {:20s} != {:20s}".format(valName, str(k), "-"  , str(d2[k])))

      #Regex stuff :
      s1 = set([l[0] for l in d1._regex])
      s2 = set([l[0] for l in d2._regex])
      for k in s1.intersection(s2) :
         iReg1 = [l[0] for l in d1._regex].index(k)
         iReg2 = [l[0] for l in d2._regex].index(k)
         comparePyFoam( d1._regex[iReg1][2], d2._regex[iReg2][2] , name1 , name2,  valName + "/" + str(k) )
      for k in s1.difference(s2) :
         iReg1 = [l[0] for l in d1._regex].index(k)
         print("{:48s} {:25s} {:20s} != {:20s}".format(valName, str(k), str(d1._regex[iReg1][2]) , "-"))
      for k in s2.difference(s1) :
         iReg2 = [l[0] for l in d2._regex].index(k)
         print("{:48s} {:25s} {:20s} != {:20s}".format(valName, str(k), "-"  , str(d2._regex[iReg2][2])))


   elif type(d1)== type(d2) == list :
      if len(d1) == len(d2) :
         for i in range(len(d1)) :
            comparePyFoam(d1[i] , d2[i] , name1 = "1" , name2 = "2", valName = valName + "/" + str(i+1))
      else :
         print("list {:s20} is different {:} != {:}".format( valName, d1 , d2 ))

   else :
      if d1.__repr__() != d2.__repr__() :
         try :
            if float(d1) != float(d2) :
               print("{:74s} {:20s} != {:20s}".format( valName, str(d1) , str(d2) ))
         except :
            print("{:74s} {:20s} != {:20s}".format( valName, str(d1) , str(d2) ))


if __name__ == "__main__" :

   if True :
      import argparse
      parser = argparse.ArgumentParser(description='FoamStar log parser')
      parser.add_argument( "file1" )
      parser.add_argument( "file2" )
      args = parser.parse_args()
      compareCase(args.file1 ,args.file2)
   
   else:
      file1 = r"\\10.67.24.192\bigr\openFoam\2Dwave\Run_rev1\g1_foamExtend_Speed_0.0"
      file2 = r"\\10.67.24.192\bigr\openFoam\2Dwave\Run\g1_foamExtendVuko_Speed_0.0"

      compareCase(file1 ,file2)
