import os
from ofRun import OfRun

from inputFiles import ControlDict, FvSchemes, FvSolution, DecomposeParDict
from inputFiles import DynamicMeshDict, WaveProperties, RelaxZone, WaveCondition, TransportProperties
from inputFiles import BoundaryPressure, BoundaryVelocity, BoundaryPointDisplacement, BoundaryAlpha

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from meshTools import getBounds


class DropTestCase( OfRun ) :

    additionalFiles = ["boundaryPressure" , "boundaryVelocity" ,"boundaryPointDisplacement" , "boundaryAlpha", "sixDofDomainBody"]
    handledFiles = OfRun.handledFiles + additionalFiles

    def __init__(self, *args, **kwargs):
        """Same arguments as OfRun + some boundary input files
        """

        for fattr in self.additionalFiles:
            if fattr in kwargs.keys():
                setattr(self, fattr,  kwargs.pop(fattr))

        OfRun.__init__(self, *args, **kwargs)

    @classmethod
    def BuildFromParams(cls, case,
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
                             waveType         = "noWaves",
                             waveH            = 0.0,
                             waveT            = 0.0,
                             velocity         = 0.0,
                             depth            = 100.,
                             inletRelaxZone   = None,
                             outletRelaxZone  = None,
                             sideRelaxZone    = None,
                             dispSignal       = None,
                             rhoWater         = 1025.,
                             application      = "foamStar",
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

        controlDict = ControlDict.Build(case                = case,
                                        startFrom           = startTime,
                                        endTime             = endTime,
                                        deltaT              = timeStep,
                                        writeInterval       = writeInterval,
                                        purgeWrite          = purgeWrite,
                                        writePrecision      = 15,
                                        forcesPatch         = forcesPatch,
                                        pressuresPatch      = pressuresPatch,
                                        outputInterval      = outputInterval,
                                        rhoWater            = rhoWater,
                                        OFversion           = OFversion,
                                        application         = application )

        #fvSchemes
        fvSchemes = FvSchemes.Build(case                    = case,
                                    simType                 = scheme,
                                    orthogonalCorrection    = "implicit",
                                    application             = application )

        #fvSolution
        fvSolution = FvSolution.Build(case          = case,
                                      useEuler      = (scheme=='Euler'),
                                      application   = application )

        #decomposeParDict
        decomposeParDict = DecomposeParDict.Build(case   = case,
                                                  nProcs = nProcs )
                                                  
        #waveProperties
        waveCond  = WaveCondition( waveType   = waveType )

        relaxZones = []
        bBox = getBounds(os.path.join(meshDir,str(meshTime),'polyMesh','points.gz'))
            
        if inletRelaxZone is not None:
            if inletRelaxZone > 0: relaxFront = RelaxZone( "inlet" , relax=True, waveCondition=waveCond, origin=[bBox[0][0], 0., 0.], orientation = [ -1., 0., 0.], bound=inletRelaxZone)
            else: relaxFront = RelaxZone( "inlet" , relax=True, waveCondition=waveCond, origin=[bBox[0][0], 0., 0.], orientation = [ -1., 0., 0.], bound=0.1*bBox[0][0])
            relaxZones += [relaxFront]
            
        if outletRelaxZone is not None:
            if outletRelaxZone > 0: relaxBack = RelaxZone( "outlet" , relax=True, waveCondition=waveCond, origin=[bBox[0][1], 0., 0.], orientation = [ 1., 0., 0.], bound=outletRelaxZone)
            else: relaxBack = RelaxZone( "outlet" , relax=True, waveCondition=waveCond, origin=[bBox[0][1], 0., 0.], orientation = [ 1., 0., 0.], bound=0.1*bBox[0][1])
            relaxZones += [relaxBack]
        
        if sideRelaxZone is not None:
            if sideRelaxZone > 0: relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[1][1], 0.], orientation = [  0., -1., 0.], bound=sideRelaxZone)
            else: relaxSide = RelaxZone( "side" , relax=True, waveCondition=waveCond, origin=[0., bBox[1][1], 0.], orientation = [  0., -1., 0.], bound=0.5*bBox[1][1])
            relaxZones += [relaxSide]
            
        waveProperties = WaveProperties.Build(case              = case,
                                              initWaveCondition = waveCond,
                                              relaxZones        = relaxZones,
                                              application       = application )

        #dynamicMeshDict
        dynamicMeshDict = DynamicMeshDict.Build_imposed(case        = case,
                                                        dispFile    = "dispSignal.dat",
                                                        cog         = COG,
                                                        OFversion   = OFversion,
                                                        application = application )



        #transportProperties
        transportProperties = TransportProperties.Build(case        = case,
                                                        rhoWater    = rhoWater,
                                                        application = application )
                                                        
        #alpha water, p_rgh, U, pointDisplacement
        boundaryAlpha = BoundaryAlpha.Build(case        = case,
                                            symmetry    = symmetry,
                                            case2D      = (ndim==2),
                                            wave        = False,
                                            struct = '"' + hullPatch + '.*"',
                                            application = application)
    
        boundaryVelocity = BoundaryVelocity.Build(case          = case,
                                                  symmetry      = symmetry,
                                                  case2D      = (ndim==2),
                                                  wave        = False,
                                                  struct = '"' + hullPatch + '.*"',
                                                  application   = application)
    
        boundaryPressure = BoundaryPressure.Build(case          = case,
                                                  symmetry      = symmetry,
                                                  case2D      = (ndim==2),
                                                  struct = '"' + hullPatch + '.*"',
                                                  application   = application)
                            
        boundaryPointDisplacement = BoundaryPointDisplacement.Build(case = case,
                                                                    symmetry = symmetry,
                                                                    case2D      = (ndim==2),
                                                                    application   = application)


        res =  cls( case, nProcs                    = nProcs,
                          OFversion                 = OFversion,
                          controlDict               = controlDict,
                          fvSchemes                 = fvSchemes,
                          fvSolution                = fvSolution,
                          dynamicMeshDict           = dynamicMeshDict,
                          transportProperties       = transportProperties,
                          decomposeParDict          = decomposeParDict,
                          turbulenceModel           = turbulenceModel,
                          waveProperties            = waveProperties,
                          boundaryAlpha             = boundaryAlpha,
                          boundaryVelocity          = boundaryVelocity,
                          boundaryPressure          = boundaryPressure,
                          boundaryPointDisplacement = boundaryPointDisplacement,
                          gravity                   = gravity,
                          application               = application
                        )

        res.dispSignal = dispSignal
        res.nRelaxZones = len(relaxZones)
        res.translate = translate
        res.rotate = rotate
        res.COG = COG
        res.hullPatch = hullPatch
        res.symmetry = symmetry
        res.ndim = ndim

        res.copyMesh(meshDir, meshTime)
        res.setBoundaries()

        res.writeFiles()
        return res

    def writeAllinit(self):
        #Allinit
        print('Create init script')
        ainit = os.path.join(self.case,'Allinit')
        with open(ainit,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')
            
            #copy disp file
            if self.dispSignal is not None:
                f.write('function copyDisp()\n')
                f.write('{\n')
                fin  = os.path.join(os.getcwd(),self.dispSignal)
                fout = os.path.join(os.getcwd(),self.case,"dispSignal.dat")
                f.write('    cp {:s} {:s}\n'.format(fin,fout))
                f.write('}\n\n')
            
            #set Relax blending zone
            if self.nRelaxZones>0:
                f.write('function setRelax()\n')
                f.write('{\n')
                f.write('    setSet -batch {}\n'.format(os.path.join('system','setSet.relax')))
                f.write('    setsToZones -noFlipMap\n')
                f.write('}\n\n')
            
            f.write('(\n')
            if self.dispSignal is not None: f.write('    copyDisp\n')
            if self.nRelaxZones>0: f.write('    setRelax\n')
            f.write('    cp -rf 0/org/* 0/\n')
            if any(self.rotate)!=0.: f.write('    BVtransformPoints -EulerZYX "( {:.3f} {:.3f} {:.3f})"  -CoR "( {:.3f} {:.3f} {:.3f})" \n'.format(*self.rotate,*self.COG))
            if any(self.translate)!=0.: f.write('    BVtransformPoints -translate "( {:.3f} {:.3f} {:.3f})"\n'.format(*self.translate))
            if self.nProcs>1:
                f.write('    decomposePar -cellDist\n')
                f.write('    mpirun -np {} initWaveField -parallel\n'.format(self.nProcs))
            else:
                f.write('    initWaveField\n')
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