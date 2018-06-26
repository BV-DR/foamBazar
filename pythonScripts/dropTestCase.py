import os
import shutil


from ofCase import OfCase

from fsTools import findCFDBoundingBox

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.dynamicMeshDict import DynamicMeshDict
from inputFiles.waveProperties import WaveProperties, RelaxZone, WaveCondition
from inputFiles.boundaryCondition import writeAllBoundaries
from inputFiles.transportProperties import TransportProperties
from inputFiles.gravity import Gravity

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile


class DropTestCase( OfCase ) :

    @classmethod
    def BuildFromAllParameters(cls,      case,
                                         section_name ,
                                         meshDir          = "mesh",
                                         meshTime         = "constant",
                                         symmetry         = False,
                                         outputForces     = False,
                                         outputPressures  = False,
                                         outputInterval   = 1,
                                         hullPatch        = "",
                                         startTime        = "latestTime",
                                         endTime          = 10,
                                         timeStep         = 0.01,
                                         writeInterval    = 1,
                                         purgeWrite       = 0,
                                         scheme           = "Euler",
                                         nProcs           = 4,
                                         nOuterCorrectors = 5,
                                         wave             = "noWaves",
                                         waveH            = 0.0,
                                         waveT            = 0.0,
                                         velocity         = 0.0,
                                         depth            = 100.,
                                         sideRelaxZone    = None,
                                         dispSignal       = None,
                                         solver           = "foamStar",
                                         OFversion        = 3,
                                         translateLength  = 0.0,
                                         gravity          = 0.0,
                                         turbulenceModel  = "laminar",
                                         ):

        print('Create system folder input files')
        #controlDict
        if outputForces:
            if len(hullPatch) > 0: forcesPatch = hullPatch
            else: forcesPatch = section_name
        else:
            forcesPatch = None
            
        if outputPressures:
            if len(hullPatch) > 0: pressuresPatch = hullPatch
            else: pressuresPatch = section_name
        else:
            pressuresPatch = None

        controlDict = ControlDict( case                = case,
                                   startFrom           = startTime,
                                   endTime             = endTime,
                                   deltaT              = timeStep,
                                   writeInterval       = writeInterval,
                                   purgeWrite          = purgeWrite,
                                   writePrecision      = 15,
                                   forcesPatch         = forcesPatch,
                                   pressuresPatch      = pressuresPatch,
                                   outputInterval      = outputInterval,
                                   rhoWater            = 1025,
                                   OFversion           = OFversion,
                                   version             = solver )

        #fvSchemes
        fvSchemes = FvSchemes( case     = case,
                               simType  = scheme,
                               orthogonalCorrection = "implicit",
                               version  = solver )

        #fvSolution
        fvSolution = FvSolution( case    = case,
                                 useEuler = scheme == 'Euler',
                                 version  = solver )

        #decomposeParDict
        decomposeParDict = DecomposeParDict( case   = case,
                                             nProcs = nProcs )


        #dynamicMeshDict
        dynamicMeshDict = DynamicMeshDict( case      = case,
                                           type      = 'solid',
                                           dispFile  = "dispSignal.dat",
                                           OFversion = OFversion,
                                           version   = solver )



        #transportProperties
        transportProperties = TransportProperties( case = case,
                                                   rhoWater = 1025,
                                                   version = solver )



        res =  cls( case, nProcs=nProcs,
                          controlDict=controlDict,
                          fvSchemes=fvSchemes,
                          fvSolution=fvSolution,
                          #waveProperties=waveProperties,
                          dynamicMeshDict=dynamicMeshDict,
                          transportProperties=transportProperties,
                          decomposeParDict=decomposeParDict,
                          turbulenceModel = turbulenceModel ,
                          gravity = gravity,
                          solver = solver
                        )

        res.dispSignal = dispSignal
        res.sideRelaxZone = sideRelaxZone
        res.translateLength = translateLength
        res.section_name = section_name
        res.symmetry = symmetry

        res.copyMesh( meshDir, meshTime )

        #Write controlDict to be able to run checkMesh
        controlDict.writeFile()


        #waveProperties
        filename = os.path.join(case,'constant','waveProperties')
        waveCond  = WaveCondition( waveType   = wave )

        relaxZones = []
        if sideRelaxZone is not None:
            bBox = findCFDBoundingBox(case,False)
            if sideRelaxZone > 0:
                relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0. , -1. , 0.], bound=sideRelaxZone)
            else:
                relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0. , -1. , 0.], bound=0.5*bBox[4])
            relaxZones += [relaxSide]

        res.waveProperties = WaveProperties( filename,
                                         initWaveCondition = waveCond,
                                         relaxZones        = relaxZones,
                                         version           = solver )

        #cell.Set
        if sideRelaxZone is not None:
            filename = os.path.join(case,'cell.Set')
            res.waveProperties.writeBlendingZoneBatch(filename)

        res.writeFiles()
        return res


    def writeFiles(self) :
        OfCase.writeFiles(self)
        shutil.copyfile(self.dispSignal, os.path.join(self.case,"dispSignal.dat"))

        self.setBoundaries()

        #alpha water, p_rgh, U, pointDisplacement
        writeAllBoundaries( case  = self.case,
                            case2D = True,
                            symmetryPlane = self.symmetry,
                            struct = '"' + self.section_name + '|wetSurf"',
                            version = self.solver )


    def writeAllinit(self):
        #Allinit
        print('Create init script')
        ainit = os.path.join(self.case,'Allinit')
        with open(ainit,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')
            f.write('(\n')
            # f.write('    cp -r ../../snap/Grid1/constant/polyMesh/ constant/\n')
            # f.write('    cp constant/org/boundary constant/polyMesh/\n')
            if self.sideRelaxZone is not None:
                f.write('    setSet -batch cell.Set\n')
                f.write('    setsToZones -noFlipMap\n')
            f.write('    cp -rf 0/org/* 0/\n')
            if self.translateLength != 0. :
                f.write('    transformPoints -translate "( 0. 0. {})"\n'.format(self.translateLength))
            f.write('    decomposePar -cellDist\n')
            f.write('    mpirun -np {:d} initWaveField -parallel\n'.format(self.nProcs))
            f.write(') 2>&1 | tee log.init\n')
        os.chmod(ainit, 0o755)

    def writeAllclean(self):
        #Allclean
        aclean = os.path.join(self.case,'Allclean')
        with open(aclean,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('cd ${0%/*} || exit 1    # run from this directory\n\n')
            f.write('function clean_log()\n')
            f.write('{\n')
            f.write('    rm -fr log.*\n')
            f.write('}\n\n')
            f.write('function clean_mesh()\n')
            f.write('{\n')
            f.write('    rm -fr background.msh VTK\n')
            f.write('    rm -fr 0/{ccx,ccy,ccz,*Level,polyMesh/cellMap}*\n')
            f.write('    rm -fr constant/modeshapes/{cfd,mapper,modalInfo}_body\n')
            f.write('    rm -fr constant/extendedFeatureEdgeMesh/\n')
            f.write('    rm -fr constant/polyMesh/{sets,*Level*,*level*,*Index*,*History*}\n')
            f.write('}\n\n')
            f.write('function clean_parallel_mesh()\n')
            f.write('{\n')
            f.write('    rm -fr processor*\n')
            f.write('}\n\n')
            f.write('function clean_0()\n')
            f.write('{\n')
            f.write('    rm -fr 0/*.gz\n')
            f.write('}\n\n')
            f.write('eval clean_log\n')
            f.write('eval clean_mesh\n')
            f.write('eval clean_parallel_mesh\n')
            f.write('eval clean_0\n')
        os.chmod(aclean, 0o755)


    def setBoundaries(self):
        boundfile = os.path.join(self.constfolder_,"polyMesh","boundary")
        boundDict = ParsedParameterFile(boundfile,boundaryDict=True)
        nbound = int(len(boundDict)/2)
        for i in range(nbound):
            if boundDict[2*i] in ['domainX0','domainX1']:
                boundDict[2*i+1]['type'] = 'empty'
            elif self.symmetry and (boundDict[2*i] in ['domainY0']):
                boundDict[2*i+1]['type'] = 'symmetryPlane'
        boundDict.writeFile()



    def writeSbatch(self):
        #run.sh
        run = os.path.join(self.case,'run.sh')
        with open(run,'w') as f:
            f.write('#!/bin/bash -l\n')
            f.write('#SBATCH -J {}\n\n'.format(self.case))
            f.write('# 5 hour wall-clock\n')
            f.write('#SBATCH --account I1608251\n')
            f.write('#SBATCH -t 3-00:00:00\n')
            f.write('#SBATCH -n {:d}\n'.format(self.nProcs))
            f.write('#SBATCH -o log.run-%j\n\n')
            f.write('module load gcc/4.9.3 openmpi/1.8.4-gcc lapack/3.6.1/gcc/4.9.3\n')
            f.write('export FOAM_INST_DIR=/data/I1608251/OpenFOAM;\n')
            if   self.OFversion == 2 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-2.4.x/etc/bashrc;\n')
            elif self.OFversion == 3 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-3.0.x/etc/bashrc;\n')
            elif self.OFversion == 5 : f.write('source /data/I1608251/OpenFOAM/OpenFOAM-5.x/etc/bashrc;\n')
            f.write('export LC_ALL=C\n\n')
            f.write('mpirun {} -parallel\n'.format(self.solver))
