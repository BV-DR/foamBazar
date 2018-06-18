
import os
import shutil
import subprocess
from inputFiles.turbulenceProperties import TurbulenceProperties, RASProperties
from inputFiles.gravity import Gravity
from fsTools import getFoamTimeFolders

class OfCase(object):
    """ Base class for openFoam case
    Can be sub-classed to deal with more specific situation (example : DropTestCase)
    """

    def __init__(self, case, nProcs=1, controlDict=None,
                                       fvSchemes=None,
                                       fvSolution=None,
                                       waveProperties=None,
                                       dynamicMeshDict=None,
                                       transportProperties=None,
                                       decomposeParDict=None,
                                       turbulenceModel = "laminar",
                                       symmetry = False,
                                       gravity = 9.81,
                                       OFversion = 5,
                                       version = "foamStar"
                                       ) :


        self.case = case  # path to case

        self.clean()
        self.nProcs = nProcs

        self.controlDict = controlDict
        self.fvSchemes = fvSchemes
        self.fvSolution = fvSolution
        self.waveProperties = waveProperties
        self.dynamicMeshDict = dynamicMeshDict
        self.transportProperties = transportProperties
        self.decomposeParDict = decomposeParDict

        self.turbulenceModel = turbulenceModel
        self.turbulenceProperties = TurbulenceProperties(case , turbulenceModel )
        self.RASProperties = RASProperties(case , turbulenceModel )

        self.symmetry = symmetry
        self.gravity = gravity

        self.OFversion = OFversion
        self.version = "foamStar"


        print('Create file tree')
        self.sysfolder_ = os.path.join(self.case,"system")
        self.zerofolder_ = os.path.join(self.case,"0")
        self.orgfolder_ = os.path.join(self.case,"0","org")
        self.constfolder_ = os.path.join(self.case,"constant")
        if not os.path.exists(self.sysfolder_):
            os.makedirs(self.sysfolder_)
        if not os.path.exists( os.path.join(self.zerofolder_, "org")):
            os.makedirs( os.path.join(self.zerofolder_, "org"))
        if not os.path.exists( self.zerofolder_):
            os.makedirs( self.zerofolder_)
        if not os.path.exists( self.constfolder_):
            os.makedirs( self.constfolder_)

        self.createParaviewFile()


    def clean(self,overwrite = False) :
        if os.path.exists(self.case):
            if overwrite:
                shutil.rmtree(self.case)
            else:
                valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
                res = input('Case "{}" already exists, do you want to overwrite ? (y/n) '.format(self.case)).lower()
                if valid.get(res,False):
                    shutil.rmtree(self.case)
                else:
                    print('Exiting')
                    os._exit(1)


    @classmethod
    def Clone( cls, case, source ):
        """Read file from existing case folder
        """
        #TODO
        return cls( case )


    def createParaviewFile(self):
        #create file for Paraview
        open(os.path.join(self.case,'view.foam'), 'a').close()

    @classmethod
    def BuildFromAllParameters(cls, case):
        return OfCase(case)


    def writeFiles(self) :

        gravity = Gravity( case = self.case, g = self.gravity )
        gravity.writeFile()

        self.controlDict.writeFile()
        self.fvSchemes.writeFile()
        self.fvSolution.writeFile()
        self.turbulenceProperties.writeFile()
        self.RASProperties.writeFile()
        self.decomposeParDict.writeFile()
        self.transportProperties.writeFile()
        self.dynamicMeshDict.writeFile()

        if self.waveProperties is not None :
            self.waveProperties.writeFile()





    def runInit(self) :
        #run Allrun script
        if not os.path.exists( "Allclean" ):
            self.writeAllClean()

        if not os.path.exists( "AllInit" ):
            self.writeAllInit()

        p = subprocess.Popen(['./Allclean'], cwd=self.case)
        p.wait()

        p = subprocess.Popen(['./AllInit'], cwd=self.case)
        p.wait()

    def run(self) :
        #run Allrun script
        if not os.path.exists( "run" ):
            self.writeRun()
        p = subprocess.Popen(['./run'], cwd=self.case)
        p.wait()

    def writeAllInit(self):
        """To be implemented in subclass"""
        raise(NotImplementedError)

    def writeRun(self):
        """To be implemented in subclass"""
        raise(NotImplementedError)

    def Allclean(self):
        """To be implemented in subclass"""
        raise(NotImplementedError)

    def copyMesh(self, meshDir, meshTime, overwrite=False):
        if meshTime=='latestTime':
            timeFolders = getFoamTimeFolders(meshDir)
            meshTimeFolder = timeFolders[-1]
        elif meshTime=='constant':
            meshTimeFolder = 'constant'
        else:
            meshTimeFolder = meshTime

        print('Copy mesh from folder ' + meshTimeFolder)

        shutil.copytree( os.path.join( meshDir , meshTimeFolder ,'polyMesh') , os.path.join( self.case , "constant/polyMesh"))
        shutil.copytree( os.path.join( meshDir , "constant" ,'triSurface') , os.path.join( self.case  , "constant/triSurface"))

