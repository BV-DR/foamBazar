import os
from ideFoam.ofRun import OfRun
from ideFoam.inputFiles import FvSchemes, FvSolution, ControlDict, DecomposeParDict, TransportProperties, DynamicMeshDict
from ideFoam.inputFiles import BoundaryPressure, BoundaryVelocity, BoundaryAlpha
from ideFoam.inputFiles import RelaxZone, noWaves, WaveCondition, WaveProperties
from ideFoam.inputFiles import createLinearWaveProbesList
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

#from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from pythonScripts.meshTools import getBounds

class Wave2DCase(OfRun):
    """Class used to generate a CFD seakeeping case.
    It should be called by a python script as presented in the following example

    Example
    -------

    For hydro-elastic case

    >>> import os
    >>> #
    >>> # Import case routine from foamBazar
    >>> # WARNING : foamBazar must be added to PYTHONPATH !
    >>> from ideFoam.seakeepingCase import SeakeepingCase
    >>> #
    >>> # Case name and additional optional parameters
    >>> case = "run"
    >>> myParams = {'nProcs'            : 12,
    >>>             'meshDir'           : 'mesh',
    >>>             'meshTime'          : 13,
    >>>             'endTime'           : 1.0,
    >>>             'timeStep'          : 0.075,
    >>>             'writeInterval'     : 40,
    >>>             'purgeWrite'        : 5,
    >>>             'outputWave'        : False,
    >>>             'fsiTol'            : 5e-5,
    >>>             'wave'              : None,
    >>>             'depth'             : 443,
    >>>             'draft'             : 11.75,
    >>>             'speed'             : 2.5722,
    >>>             'waveType'          : 'noWaves',
    >>>             'waveH'             : 11.1,
    >>>             'waveT'             : 12.75,
    >>>             'waveStartTime'     : 0,
    >>>             'waveRampTime'      : 12,
    >>>             'addDamping'        : True,
    >>>             'EulerCellsDist'    : 4,
    >>>             'inletRelaxZone'    : 255.,
    >>>             'outletRelaxZone'   : 51.,
    >>>             'sideRelaxZone'     : 150.,
    >>>             'datFile'           : 'homer/4400_Bulk.dat',
    >>>             'donFile'           : 'homer/4400.don',
    >>>             'dmigFile'          : 'homer/4400_dmig.pch',
    >>>             'mdFile'            : 'homer/4400_md.pch',
    >>>             'hmrUserOutput'     : 'homer/HmFEM.out',
    >>>             'modesToUse'        : [7,8,9],
    >>>             'shipDamping'       : [0.0,0.0,0.0],
    >>>             'clean'             : True
    >>>             }
    >>> #
    >>> # Call routines for meshing here
    >>> run = SeakeepingCase.BuildFromParams( case, **myParams )
    >>> fname = os.path.join(case,'log.input')
    >>> with open(fname,'w') as f: f.write(str(myParams))
    >>> run.runInit()

    """

    additionalFiles = ["boundaryPressure" , "boundaryVelocity" , "boundaryAlpha"]
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
                             meshDir            = 'mesh',
                             meshTime           = 'constant',
                             startTime          = 'latestTime',
                             endTime            = 1000,
                             timeStep           = 0.01,
                             adjustTimeStep     = None,
                             writeInterval      = 1,
                             purgeWrite         = 0,
                             writeFormat        = "ascii",
                             symmetry           = 1,
                             turbulenceModel    = "laminar",

                             waveProbes         = [],

                             outputInterval     = 1,

                             inletRelaxZone     = None,
                             outletRelaxZone    = None,
                             outletRelaxTarget  = 'still',
                             sideRelaxZone      = None,

                             wave               = None,
                             waveType           = None,
                             waveH              = None,
                             waveT              = None,
                             speed              = 0.,
                             depth              = 500,
                             draft              = 0,
                             waveStartTime      = 0,
                             waveRampTime       = 0,
                             addDamping         = False,
                             vtkOut             = True,

                             ddtScheme          = 'Euler',
                             rhoWater           = 1025.,
                             nProcs             = 1,
                             OFversion          = 5,
                             application        = 'foamStar',
                             clean              = False,
                             ):
        """Build mesh for CFD seakeeping case from a few parameters.

        Parameters
        ----------
        case : str
            Name of case to create
        stlFile : str, default 'ship.stl'
            Path to ship STL file (used for mesh snapping)
        meshDir : str, default 'mesh'
            Path to mesh directory
        meshTime : int or str, default'constant'
            Name of folder containing mesh in meshDir
        startTime : str, default 'latestTime'
            String defining simulation startTime parameter
        endTime : float, default 1000
            End simulation time
        timeStep : float, default 0.01
            Simulation time step
        adjustTimeStep : list of float, default None
            Definition of adjustable time step : None for constant time step, [maxCo, maxAlphaCo, maxDeltaT] otherwise
        writeInterval : int, default 1
            Define interval at wich simualtion results are written
        purgeWrite : int, default 0
            Define number of results time steps after wich results are erased
        writeFormat : str, default 'ascii'
            TODO
        symmetry : int, default 1
            Define symmetry type : 0 = None, 1 = symmetry, 2 = symmetryPlane

        outputWave : bool, default False
            Logical defining if wave probes are output
        waveProbes : TODO, default []
            TODO
        outputMotions : bool, default True
            Logical defining if global motions are output
        localMotionPts : TODO, default []
            Define list of points for local motion output
        outputForces : bool, default False
            Logical defining if forces are output on forcesPatch
        forcesPatch : list of str, default None
            List of patch where forces are output
        outputPressures : bool, default False
            Logical defining if pressures are output on pressuresPatch
        pressuresPatch : list of str, default None
            List of patch where pressures are output
        outputInterval : int, default 1
            Define interval at wich simualtion outputs are written

        inletRelaxZone : float, default None
            Length of inlet relaxation zone. None for no relaxation zone.
        outletRelaxZone : float, default None
            Length of outlet relaxation zone. None for no relaxation zone.
        outletRelaxTarget : str, default 'still'
            Target for outlet relaxation zone. Either 'still' (wave damping zone) or 'incident'.
        sideRelaxZone : float, default None
            Length of side relaxation zone. None for no relaxation zone.

        wave : ideFoam.waveProperties.WaveCondition, default None
            Wave condition defined by ideFoam.waveProperties.WaveCondition object. Otherwise, define wave condition with the following arguments (waveType, waveH, waveT, ... )
        waveType : str, default None
            Wave type used for simulation.
            options available: "noWave", ""
        waveH : float, default None
            Wave height
        waveT : float, default  None
            Wave period
        speed : float, default 0.
            Ship speed (in m/s)
        depth : float, default 500
            Depth
        draft : float, default 0.
            Draft (measured from keel)
            if defined, ship.stl will be moved into position
            if "None", no operation will be performed to ship.stl
        waveStartTime : float, default 0.
            Wave start time
        waveRampTime : float, default 0.
            Duration of ramp applied to wave
        addDamping : bool, default False
            Logical defining if additional artificial damping should be added. This is typically used to speed up ship balancing.
        vtkOut : bool, default True
            TODO
        ddtScheme, str, default 'Euler'
            TODO
        rhoWater : float, default 1025.
            Water density
        nProcs : int, default 1
            Number of processors
        OFversion : int or str, default 5
            OpenFOAM version number
        application : str, default 'foamStar'
            Application to run
        nProcs : int, default 4
            Number of processors used to build the mesh
        OFversion : int or str, default 5
            OpenFOAM version
        onLiger : boolean, default False
            Logical defining if case is run on Liger cluster
        clean : boolean, default False
            Logical to force case overwrite

        """

        #Read data from Homer or from inputs


        wpList = None
        if len(waveProbes) > 0:
            print('WARNING: Wave probes cannot be included yet')
            wpList = waveProbes
            #wpList = []
            #for wp in waveProbes:
            #    wpList += createLinearWaveProbesList(*wp)

        symmetry = "2D"

        controlDict = ControlDict.Build(case                = case,
                                        startFrom           = startTime,
                                        endTime             = endTime,
                                        deltaT              = timeStep,
                                        adjustTimeStep      = adjustTimeStep,
                                        writeInterval       = writeInterval,
                                        purgeWrite          = purgeWrite,
                                        writePrecision      = 15,
                                        writeFormat         = writeFormat,
                                        outputInterval      = outputInterval,
                                        rhoWater            = rhoWater,
                                        waveProbesList      = wpList,
                                        application         = application )

        #fvSchemes
        fvSchemes = FvSchemes.Build(case                    = case,
                                    simType                 = ddtScheme,
                                    limitedGrad             = False,
                                    orthogonalCorrection    = "implicit",
                                    application             = application )

        #fvSolution
        fvSolution = FvSolution.Build(case          = case,
                                      application   = application )

        #decomposeParDict
        decomposeParDict = DecomposeParDict.Build(case    = case,
                                                  nProcs  = nProcs )

        #waveProperties
        if wave is not None:
            waveCond = wave
        elif all((waveType,waveH,waveT)):
            waveCond  = WaveCondition(waveType   = waveType,
                                      height     = waveH,
                                      period     = waveT,
                                      U0         = -1.*speed,
                                      depth      = depth,
                                      refDirection = [-1,0,0],
                                      startTime  = waveStartTime,
                                      rampTime   = waveRampTime
                                      )
        else:
            raise(Exception('No or incomplete wave properties provided. Please provide "wave" or all the following parameters ("waveType","waveH","waveT").'))

        relaxZones = []
        bBox = getBounds(os.path.join(meshDir,str(meshTime),'polyMesh'))
        if sideRelaxZone is not None:
            relaxSide   = RelaxZone( "side"  , relax=True, waveCondition=waveCond, origin=[0., bBox[1][1], 0.], orientation = [  0. , -1. , 0.], length=sideRelaxZone)
            relaxZones += [relaxSide]
        if outletRelaxZone is not None:
            if outletRelaxTarget == "still":
                relaxOutlet = RelaxZone( "outlet", relax=True, waveCondition=noWaves, origin=[bBox[0][0], 0., 0.], orientation = [  1. ,  0. , 0.], length=outletRelaxZone)
            elif outletRelaxTarget == "incident":
                relaxOutlet = RelaxZone( "outlet", relax=True, waveCondition=waveCond, origin=[bBox[0][0], 0., 0.], orientation = [  1. ,  0. , 0.], length=outletRelaxZone)
            else:
                raise(ValueError('"still" of "incident" expected for outletRelaxTarget'))
            relaxZones += [relaxOutlet]
        if inletRelaxZone is not None:
            relaxInlet  = RelaxZone( "inlet" , relax=True, waveCondition=waveCond, origin=[bBox[0][1], 0., 0.], orientation = [ -1. ,  0. , 0.], length=inletRelaxZone)
            relaxZones += [relaxInlet]

        waveProperties = WaveProperties.Build(case              = case,
                                              initWaveCondition = waveCond,
                                              relaxZones        = relaxZones,
                                              application       = application )

        #transportProperties
        transportProperties = TransportProperties.Build(case        = case,
                                                        nuWater     = 1.14e-06,
                                                        rhoWater    = rhoWater,
                                                        application = application )


        #alpha water, p_rgh, U, pointDisplacement
        boundaryAlpha = BoundaryAlpha.Build(case        = case,
                                            symmetry    = symmetry,
                                            application = application)

        boundaryVelocity = BoundaryVelocity.Build(case          = case,
                                                  speed         = speed,
                                                  symmetry      = symmetry,
                                                  application   = application)

        boundaryPressure = BoundaryPressure.Build(case          = case,
                                                  symmetry      = symmetry,
                                                  application   = application)

        dynamicMeshDict = DynamicMeshDict.Build_static(case = case)


        res =  cls( case, nProcs                    = nProcs,
                          OFversion                 = OFversion,
                          controlDict               = controlDict,
                          turbulenceModel           = turbulenceModel,
                          fvSchemes                 = fvSchemes,
                          fvSolution                = fvSolution,
                          transportProperties       = transportProperties,
                          decomposeParDict          = decomposeParDict,
                          dynamicMeshDict           = dynamicMeshDict,
                          waveProperties            = waveProperties,
                          boundaryAlpha             = boundaryAlpha,
                          boundaryVelocity          = boundaryVelocity,
                          boundaryPressure          = boundaryPressure,
                          clean                     = clean,
                          application               = application )


        res.symmetry = symmetry
        res.meshDir = meshDir
        res.nRelaxZones = len(relaxZones)
        res.copyMesh(meshDir,meshTime)
        res.setBoundaries()

        res.writeFiles()
        return res

    def __str__(self):
        s = OfRun.__str__(self)
        return "SeakeepingCase : " + s

    def setBoundaries(self):
        boundfile = os.path.join(self.constfolder_,"polyMesh","boundary")
        boundDict = ParsedParameterFile(boundfile,boundaryDict=True)
        nbound = int(len(boundDict)/2)
        #read and modify entries
        for i in range(nbound):
            if boundDict[2*i] in ['domainY0', 'domainY1']:
                boundDict[2*i+1]['type'] = 'empty'
                if "inGroups" in boundDict[2*i+1].keys() :
                    del boundDict[2*i+1]["inGroups"]

        boundDict.writeFile()

    def writeAllinit(self, batchName="Allinit"):
        """Write bash script 'Allinit' to initialize case
        """
        print('Create init script')
        ainit = os.path.join(self.case,'Allinit')
        with open(ainit,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')

            #set Relax blending zone
            if self.nRelaxZones > 0:
                f.write('function setRelax()\n')
                f.write('{\n')
                f.write('    setSet -batch {}\n'.format(os.path.join('system','setSet.relax')))
                f.write('    setsToZones -noFlipMap\n')
                f.write('}\n\n')

            # __main__
            f.write('(\n')
            if self.nRelaxZones > 0: f.write('    setRelax\n')
            f.write('    cp -rf ./0/org/* ./0/\n')
            if self.nProcs > 1:
                f.write('    decomposePar -cellDist\n')
                f.write('    mpirun -np {} initWaveField -parallel\n'.format(self.nProcs))
            else:
                f.write('    initWaveField\n')
            f.write(') 2>&1 | tee log.init\n\n')
        os.chmod(ainit, 0o755)

    def writeAllclean(self):
        """Write bash script 'Allclean' to clean case folder(s)
        """
        aclean = os.path.join(self.case,'Allclean')
        with open(aclean,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('cd ${0%/*} || exit 1    # run from this directory\n\n')
            f.write('function clean_log()\n')
            f.write('{\n')
            f.write('    rm -fr log.*\n')
            f.write('}\n\n')
            f.write('function clean_parallel()\n')
            f.write('{\n')
            f.write('    rm -fr processor*\n')
            f.write('}\n\n')
            f.write('function clean_0()\n')
            f.write('{\n')
            f.write('    rm -fr 0/*.gz\n')
            f.write('}\n\n')
            f.write('#eval clean_log\n')
            f.write('eval clean_parallel\n')
            f.write('eval clean_0\n')
        os.chmod(aclean, 0o755)

    def writeRun(self):
        """To be implemented in subclass"""
        with open(os.path.join(self.case, "Allrun"), "w") as f:
            if self.nProcs > 1:
                f.write("mpirun -n {} {} -parallel 2>&1 | tee foamStar.log".format(self.nProcs, self.executable))
            else:
                f.write("{} 2>&1 | tee foamStar.log".format(self.executable))

    def runInit(self):
        import subprocess
        print("Run initialization")
        os.chdir(self.case)
        subprocess.call(["/bin/bash", "Allinit"], shell=False)

    def run(self):
        import subprocess
        print("Run initialization")
        os.chdir(self.case)
        subprocess.call(["/bin/bash", "Allrun"], shell=False)
