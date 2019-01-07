
import os
from os.path import join, abspath, exists
import shutil
import subprocess

from ideFoam.inputFiles import getFileClass, getFilePath
from pythonScripts.fsTools import getFoamTimeFolders

class OfCase(object):
    """ Base class for openFoam case
    Can be sub-classed to deal with more specific situation (example : DropTestCase, SeakeepingCase)
    """

    handledFiles =  [
                    "controlDict",
                    "decomposeParDict",
                    "fvSchemes",
                    "fvSolution"
                    ]

    def __init__(self, case, nProcs=1, controlDict      = None,
                                       fvSchemes        = None,
                                       fvSolution       = None,
                                       decomposeParDict = None,
                                       symmetry         = False,
                                       OFversion        = 5,
                                       application      = 'foamStar',   #Solver name (foamStar, navalFoam)
                                       executable       = None,   # => default to application name
                                       isMesher         = False,
                                       clean            = False,  #True => Remove case folder and go on. #False: ask user interactively
                                       meshFolder       = None
                                       ) :


        self.case = abspath(case)  # path to case
        self.clean(clean)
        self.nProcs = nProcs

        #system
        self.controlDict = controlDict
        self.decomposeParDict = decomposeParDict
        self.fvSchemes = fvSchemes
        self.fvSolution = fvSolution

        self.meshFolder = meshFolder

        self.symmetry = symmetry
        self.OFversion = OFversion
        self.application = application
        self.isMesher = isMesher

        self.sysfolder_ = join(self.case,"system")
        self.zerofolder_ = join(self.case,"0")
        self.orgfolder_ = join(self.case,"0","org")
        self.constfolder_ = join(self.case,"constant")

        self.executable = executable
        if self.executable :
            self.executable = self.application


        self.writeFolders()

    def writeFolders(self) :
        print('Create file tree')
        if not exists(self.sysfolder_): os.makedirs(self.sysfolder_)
        if not exists( join(self.zerofolder_, "org")): os.makedirs( join(self.zerofolder_, "org"))
        if not exists( self.zerofolder_): os.makedirs( self.zerofolder_)
        if not exists( self.constfolder_): os.makedirs( self.constfolder_)

        if self.isMesher:
            self.polyfolder_ = join(self.case,"constant","polyMesh")
            self.trifolder_ = join(self.case,"constant","triSurface")
            if not exists( self.polyfolder_): os.makedirs( self.polyfolder_)
            if not exists( self.trifolder_): os.makedirs( self.trifolder_)

        self.createParaviewFile()


    def clean(self,clean = False) :
        if exists(self.case):
            if clean:
                print('Overwriting case "{}"'.format(self.case))
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
    def Read( cls, case, source, application = "foamStar" ):
        """Read file from existing folder
        """
        source = abspath(source)
        case = abspath(case)
        fileDict = {}
        for f in cls.handledFiles:
            fname =  join(source, getFilePath(f)  )
            if exists(fname) :
                tmpobj = getFileClass()( join(source, getFilePath(f) ) , read = True)
                tmpobj.case = case
                tmpobj.name = tmpobj.name.replace(source, case)
                fileDict[ f ] = tmpobj
            else:
                print (fname, "does not exists")
        return cls( case, **fileDict , application = application, meshFolder = join(source, "constant") )


    def createParaviewFile(self):
        #create file for Paraview
        open(join(self.case,'view.foam'), 'a').close()

    @classmethod
    def BuildFromAllParameters(cls, case, isMesher=False):
        return OfCase(case, isMesher)


    def writeFiles(self) :
        """Write all input file in case folder
        """
        #write all handeled files
        for f in self.handledFiles:
            file_ = getattr(self, f)
            if file_ is not None:
                print ("Writting :", file_.name)
                file_.writeFile()

        #write Allinit and Allclean scripts
        if not exists( "Allclean" ): self.writeAllclean()
        if not exists( "Allinit" ): self.writeAllinit()

    def runInit(self) :
        #run Allclean and Allinit script
        p = subprocess.Popen(['./Allclean'], cwd=self.case)
        p.wait()

        p = subprocess.Popen(['./Allinit'], cwd=self.case)
        p.wait()

    def run(self) :
        #run Allrun script
        if not exists( "run" ):
            self.writeRun()
        p = subprocess.Popen(['./run'], cwd=self.case)
        p.wait()

    def writeAllinit(self):
        """To be implemented in subclass"""
        raise(NotImplementedError)

    def writeAllclean(self):
        """To be implemented in subclass"""
        raise(NotImplementedError)

    def writeRun(self):
        """To be implemented in subclass"""
        raise(NotImplementedError)

    def copyMesh(self, meshDir, meshTime, overwrite=False):
        meshTime = str(meshTime)
        if meshTime=='latestTime':
            timeFolders = getFoamTimeFolders(meshDir)
            meshTimeFolder = timeFolders[-1]
        elif meshTime=='constant':
            meshTimeFolder = 'constant'
        else:
            meshTimeFolder = meshTime

        print('Copy mesh from folder ' + meshTimeFolder)

        shutil.copytree( join( meshDir , meshTimeFolder ,'polyMesh') , join( self.case , "constant/polyMesh"))
        shutil.copytree( join( meshDir , "constant" ,'triSurface') , join( self.case  , "constant/triSurface"))

    def writeSbatch(self):
        #run.sh
       	print('Write run.sh batch file')
        run = os.path.join(self.case,'run.sh')
        with open(run,'w') as f:
            f.write('#!/bin/bash -l\n')
            f.write('#SBATCH -J {}\n\n'.format(self.case))
            f.write('# 5 hour wall-clock\n')
            f.write('#SBATCH --account I1608251\n')
            if self.isMesher:
                f.write('#SBATCH -t 0-02:00:00\n')
            else:
                f.write('#SBATCH -t 3-00:00:00\n')
            f.write('#SBATCH -n {:d}\n'.format(self.nProcs))
            f.write('#SBATCH -o log.run-%j\n\n')
            f.write('module load gcc/4.9.3 openmpi/1.8.4-gcc lapack/3.6.1/gcc/4.9.3\n')
            f.write('export FOAM_INST_DIR=/data/I1608251/OpenFOAM;\n')
            if   self.OFversion == 2 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-2.4.x/etc/bashrc;\n')
            elif self.OFversion == 3 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-3.0.x/etc/bashrc;\n')
            elif self.OFversion == 5 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-5.x/etc/bashrc;\n')
            f.write('export LC_ALL=C\n\n')
            f.write('mpirun {} -parallel\n'.format(self.executable))

    def runSbatch(self):
        self.writeSbatch()
        subprocess.call(['sbatch', 'run.sh'], cwd=self.case)
