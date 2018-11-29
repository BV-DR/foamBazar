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

    def __init__(self, case, *kwargs) :
        OfCase.__init__(self, case, extrudeMeshDict=None,
                                    refineMeshDict=None,
                                    refineMeshDict1=None,
                                    refineMeshDict2=None,
                                    snappyHexMeshDict=None,
                                    surfaceFeatureExtractDict=None,
                                    blockMeshDict=None,
                                    **kwargs  #Argument of OfCase
                                    ):

        #Mesh stuff
        self.extrudeMeshDict = extrudeMeshDict
        self.refineMeshDict = refineMeshDict
        self.refineMeshDict1 = refineMeshDict1
        self.refineMeshDict2 = refineMeshDict2
        self.snappyHexMeshDict = snappyHexMeshDict
        self.surfaceFeatureExtractDict = surfaceFeatureExtractDict
        self.blockMeshDict = blockMeshDict



    @classmethod
    def BuildFromAllParameters(cls,      case,
                                         meshDir          = "mesh",
                                         meshTime         = "constant",
                                         ndim             = 2,
                                         symmetry         = 0,          # 0 = None ; 1 = symmetry ; 2 = symmetryPlane
                                         outputForces     = False,
                                         forcesPatch      = None,
                                         outputPressures  = False,
                                         pressuresPatch   = None,
                                         outputInterval   = 1,
                                         hullPatch        = "ship",
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
                                         frontRelaxZone    = None,
                                         backRelaxZone    = None,
                                         sideRelaxZone    = None,
                                         dispSignal       = None,
                                         solver           = "foamStar",
                                         OFversion        = 3,
                                         translate        = [0.0,0.0,0.0],
                                         rotate           = [0.0,0.0,0.0],
                                         COG              = [0.0,0.0,0.0],
                                         gravity          = 0.0,
                                         turbulenceModel  = "laminar",
                                         ):

        #controlDict
        if outputForces and (forcesPatch is None): forcesPatch = [hullPatch]
        if outputPressures and (pressuresPatch is None): pressuresPatch = [hullPatch]

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
                                           COG       = COG,
                                           OFversion = OFversion,
                                           version   = solver )



        #transportProperties
        transportProperties = TransportProperties( case = case,
                                                   rhoWater = 1025,
                                                   version = solver )



        res =  cls( case, nProcs=nProcs,
                          OFversion=OFversion,
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
        res.translate = translate
        res.rotate = rotate
        res.COG = COG
        res.hullPatch = hullPatch
        res.symmetry = symmetry
        res.ndim = ndim

        res.copyMesh( meshDir, meshTime )

        #Write controlDict to be able to run checkMesh
        controlDict.writeFile()

        #waveProperties
        filename = os.path.join(case,'constant','waveProperties')
        waveCond  = WaveCondition( waveType   = wave )

        relaxZones = []
        if any(x is not None for x in [frontRelaxZone,backRelaxZone,sideRelaxZone]):
            bBox = findCFDBoundingBox(case,False)
            
        if frontRelaxZone is not None:
            if frontRelaxZone > 0: relaxFront = RelaxZone( "front" , relax=True, waveCondition=waveCond, origin=[bBox[0], 0., 0.], orientation = [ -1., 0., 0.], bound=frontRelaxZone)
            else: relaxFront = RelaxZone( "front" , relax=True, waveCondition=waveCond, origin=[bBox[0], 0., 0.], orientation = [ -1., 0., 0.], bound=0.1*bBox[0])
            relaxZones += [relaxFront]
            
        if backRelaxZone is not None:
            if backRelaxZone > 0: relaxBack = RelaxZone( "back" , relax=True, waveCondition=waveCond, origin=[bBox[3], 0., 0.], orientation = [ 1., 0., 0.], bound=backRelaxZone)
            else: relaxBack = RelaxZone( "back" , relax=True, waveCondition=waveCond, origin=[bBox[3], 0., 0.], orientation = [ 1., 0., 0.], bound=0.1*bBox[3])
            relaxZones += [relaxBack]
        
        if sideRelaxZone is not None:
            if sideRelaxZone > 0: relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0., -1., 0.], bound=sideRelaxZone)
            else: relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[4], 0.], orientation = [  0., -1., 0.], bound=0.5*bBox[4])
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

        if self.extrudeMeshDict is not None:            self.extrudeMeshDict.writeFile()
        if self.refineMeshDict is not None:             self.refineMeshDict.writeFile()
        if self.refineMeshDict1 is not None:            self.refineMeshDict1.writeFile()
        if self.refineMeshDict2 is not None:            self.refineMeshDict2.writeFile()
        if self.snappyHexMeshDict is not None:          self.snappyHexMeshDict.writeFile()
        if self.surfaceFeatureExtractDict is not None:  self.surfaceFeatureExtractDict.writeFile()
        if self.blockMeshDict is not None:              self.blockMeshDict.writeFile()

        shutil.copyfile(self.dispSignal, os.path.join(self.case,"dispSignal.dat"))

        self.setBoundaries()

        #alpha water, p_rgh, U, pointDisplacement
        writeAllBoundaries( case  = self.case,
                            case2D = (self.ndim==2),
                            symmetry = self.symmetry,
                            wave = False,
                            struct = '"' + self.hullPatch + '.*"',
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
            if any(self.rotate)!=0.: f.write('    BVtransformPoints -EulerZYX "( {:.3f} {:.3f} {:.3f})"  -CoR "( {:.3f} {:.3f} {:.3f})" \n'.format(*self.rotate,*self.COG))
            if any(self.translate)!=0.: f.write('    BVtransformPoints -translate "( {:.3f} {:.3f} {:.3f})"\n'.format(*self.translate))
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
                if self.ndim>2: boundDict[2*i+1]['type'] = 'patch'
                else: boundDict[2*i+1]['type'] = 'empty'
            elif self.symmetry>0 and (boundDict[2*i] in ['domainY0']):
                if self.symmetry==1: boundDict[2*i+1]['type'] = 'symmetry'
                elif self.symmetry==2: boundDict[2*i+1]['type'] = 'symmetryPlane'
        boundDict.writeFile()

    def writeSbatch(self):
        #run.sh
       	print('Write run.sh batch file')
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
