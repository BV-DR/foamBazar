#!/usr/bin/env upython3

#########################################################################
# Filename: foamStar.py                                                 #
# Date:     2018-June-27                                                #
# Version:  1.                                                          #
# Author:   Sopheak Seng                                                #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    sopheak.seng@bureauveritas.com                              #
# Project:                                                              #
#########################################################################

import os, sys
import inputFiles.caseFolder as caseFolder


class foamStar(caseFolder.rootDir):
    def __init__(self, case):
        if not isinstance(case,str): raise(Exception('invalid case dir.:',case))
        if not case: raise(Exception('invalid case dir.:',case))
        caseFolder.rootDir.__init__(self,case,opts='READ_IF_PRESENT')
        self.addfile('./constant/dynamicMeshDict')
        self.addfile('./constant/g')
        self.addfile('./constant/transportProperties')
        self.addfile('./constant/turbulenceProperties')
        self.addfile('./system/decomposeParDict')
        self.addfile('./system/fvSchemes')
        self.addfile('./system/controlDict', application='foamStar')
        self.addfile('./system/fvSolution')
        pass

#*** Main execution start here ************************************************
if __name__ == "__main__":
    #print(str(sys.argv[0:]))
    of=foamStar("test")
    #print(of.system.fvSolution)
    print(of.system.fvSchemes)
    #print(of.system.controlDict)
    #print(of,'\n')
    #print(of.constant.dynamicMeshDict)

