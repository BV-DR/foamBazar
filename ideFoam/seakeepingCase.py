import os
from ideFoam.ofRun import OfRun
from ideFoam.inputFiles import FvSchemes, FvSolution, ControlDict, DecomposeParDict, TransportProperties, SetSelection
from ideFoam.inputFiles import BoundaryPressure, BoundaryVelocity, BoundaryPointDisplacement, BoundaryAlpha
from ideFoam.inputFiles import RelaxZone, noWaves, WaveCondition, WaveProperties, DynamicMeshDict
from ideFoam.inputFiles import SixDofDomainBody, InitFlexDict, FlexFile
from ideFoam.inputFiles import createLinearWaveProbesList

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from pythonScripts.meshTools import getBounds
from pythonScripts.fsTools import findBoundingBox

class SeakeepingCase(OfRun):
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

    additionalFiles = ["boundaryPressure" , "boundaryVelocity" ,"boundaryPointDisplacement" , "boundaryAlpha", "sixDofDomainBody", "initFlexDict", "flexFile"]
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

                             donFile            = None,
                             outputWave         = False,
                             waveProbes         = [],
                             localMotionPts     = [],
                             outputForces       = False,
                             forcesPatch        = None,
                             outputPressures    = False,
                             pressuresPatch     = None,
                             outputInterval     = 1,
                             hullPatch          = 'ship',

                             mdFile             = None,
                             modesToUse         = [],
                             datFile            = None,
                             dmigFile           = None,
                             hmrUserOutput      = None,
                             hmrScaling         = 0,
                             shipDamping        = None,

                             EulerCellsDist     = None,

                             inletRelaxZone     = None,
                             outletRelaxZone    = None,
                             outletRelaxTarget  = 'still',
                             sideRelaxZone      = None,

                             mass               = None,
                             inertia            = None,
                             COG                = None,

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

                             meshMotion         = 'cpMorphing',
                             innerDistance      = None,
                             outerDistance      = None,

                             ddtScheme          = 'Euler',
                             fsiTol             = 1e-6,
                             rhoWater           = 1025.,
                             nProcs             = 1,
                             OFversion          = 5,
                             application        = 'foamStar',
                             clean              = False,
                             stlFile            = None,    # Require only for EulerCell and addDamping
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

        hullPatch : str, default 'ship'
            Set patch name given for hull in boundary file
        donFile : str, default 'ship.don'
            Path to *.don file generated by Homer
        mdFile : str, default None
            Path to *_md.pch file generated by Nastran (via Homer)
        datFile : str, default None
            Path to *.dat file generated by Homer
        dmigFile : str, default None
            Path to *_dmig.pch file generated by Nastran (via Homer)
        hmrUserOutput : str, default None
            Path to 'HmFEM.out' file generated by Homer
        modesToUse : list of int, default []
            List of structural modes to use
        hmrScaling : list of float, deault None
            List of scaling factors used by Homer for each structural mode
        shipDamping : list of float, deault None
            List of additional damping coefficients for each structural mode

        EulerCellsDist : flaot, default None
            Distance from ship used for Euler cells blending
        inletRelaxZone : float, default None
            Length of inlet relaxation zone. None for no relaxation zone.
        outletRelaxZone : float, default None
            Length of outlet relaxation zone. None for no relaxation zone.
        outletRelaxTarget : str, default 'still'
            Target for outlet relaxation zone. Either 'still' (wave damping zone) or 'incident'.
        sideRelaxZone : float, default None
            Length of side relaxation zone. None for no relaxation zone.

        mass : float, default None
            If Homer is not used, mass of ship.
        inertia : list of floats, dimension(3), default None
            If Homer is not used, ship inertias.
        COG : list of floats, dimension(3), default None
            If Homer is not used, ship center of grativty

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

        meshMotion : TODO, default 'cpMorphing'
            TODO
        innerDistance : TODO
            TODO
        outerDistance : TODO
            TODO

        ddtScheme, str, default 'Euler'
            TODO
        fsiTol : TODO, default 1e-6
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
        nModesToUse = len(modesToUse)
        if hmrUserOutput is not None:
            print('Reading mass properties from Homer file : "{}"'.format(hmrUserOutput))
            modesToUse.sort()

            freq = []
            scaling = []
            inertia = [0.]*6

            with open(hmrUserOutput,'r') as f:
                itf = iter(f)
                for line in itf:
                    if 'Location of center of gravity in global reference' in line:
                        next(itf)
                        pline = next(itf)
                        COG = [float(i) for i in pline.split()]
                    elif 'Mass =' in line:
                        mass = float(line.split()[2])
                    elif 'Roll Inertia =' in line:
                        inertia[0] = float(line.split()[3])
                    elif 'Pitch Inertia =' in line:
                        inertia[1] = float(line.split()[3])
                    elif 'Yaw Inertia =' in line:
                        inertia[2] = float(line.split()[3])
                        break
                    elif nModesToUse > 0:
                        for mode in modesToUse:
                            mstr = 'Mode '+str(mode)
                            if mstr in line:
                                sline = line.split()
                                scaling.append(float(sline[3]))
                                freq.append(float(sline[6])**2)

            if nModesToUse > 0:
                hmrScaling  = ('{:20.14e}').format(scaling[0])
                shipFreq    = ('{:20.14e} '*len(freq)).format(*freq)
        elif mass is not None and inertia is not None and COG is not None:
            print('Reading mass properties  from user inputs')
        else:
            raise(Exception('No or incomplete mass properties provided. Please provide "hmrUserOutput" (for Homer) or all the following parameters ("mass","inertia","COG").'))

        #controlDict
        vbmPatch = None
        if donFile is not None:
            donName = os.path.splitext(os.path.basename(donFile))[0]
            vbmPatch = [hullPatch,donName]
        else:
            donName = None

        wpList = None
        if len(waveProbes) > 0:
            print('WARNING: Wave probes cannot be included yet')
            wpList = waveProbes

        if outputForces and (forcesPatch is None): forcesPatch = [hullPatch]
        if outputPressures and (pressuresPatch is None): pressuresPatch = [hullPatch]

        controlDict = ControlDict.Build(case                = case,
                                        startFrom           = startTime,
                                        endTime             = endTime,
                                        deltaT              = timeStep,
                                        adjustTimeStep      = adjustTimeStep,
                                        writeInterval       = writeInterval,
                                        purgeWrite          = purgeWrite,
                                        writePrecision      = 15,
                                        writeFormat         = writeFormat,
                                        outputLocalMotions  = len(localMotionPts)>0,
                                        vbmPatch            = vbmPatch,
                                        forcesPatch         = forcesPatch,
                                        pressuresPatch      = pressuresPatch,
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
                                      fsiTol        = fsiTol,
                                      useEulerCells = (EulerCellsDist is not None),
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

        #Euler blending cell selection
        setSelections = []
        if EulerCellsDist is not None:
            outsidePoints = [0.5*(bBox[0][0]+bBox[0][1]), 0.5*bBox[1][1], 0.5*(bBox[2][0]+bBox[2][1])]
            setSelections.append(SetSelection(case          = case,
                                             name          = 'euler',
                                             selType       = 'proximity',
                                             selName       = 'EulerCells',
                                             opts          = 'new',
                                             stlFile       = stlFile,
                                             outsidePoints = outsidePoints,
                                             distance      = EulerCellsDist ))


        #transportProperties
        transportProperties = TransportProperties.Build(case        = case,
                                                        nuWater     = 1.14e-06,
                                                        rhoWater    = rhoWater,
                                                        application = application )

        #dynamicMeshDict
        if addDamping:
            bBox = findBoundingBox(os.path.join(meshDir,'constant','triSurface',stlFile),False)
            LBP = round(0.95*abs(bBox[3]-bBox[0]),2)
            BC = round(abs(bBox[4]-bBox[1]),2)
        else:
            LBP = 0.
            BC = 0.

        # dynamicMeshDict
        # for use with Homer
        if hmrUserOutput is not None:
            dynamicMeshDict = DynamicMeshDict.Build_elastic(case        = case,
                                                            hullPatch   = hullPatch,
                                                            addDamping  = addDamping,
                                                            lpp         = LBP,
                                                            bc          = BC,
                                                            application = application )
        # for use without Homer
        elif mass is not None and inertia is not None and COG is not None:
            dynamicMeshDict = DynamicMeshDict.Build_free(case           = case,
                                                         mass           = mass,
                                                         cog            = COG,
                                                         inertia        = inertia,
                                                         rampTime       = 0.1,
                                                         releaseTime    = 0.0,
                                                         hullPatch      = "(ship)",
                                                         meshMotion     = meshMotion,
                                                         innerDistance  = innerDistance,
                                                         outerDistance  = outerDistance)

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

        boundaryPointDisplacement = BoundaryPointDisplacement.Build(case = case,
                                                                    symmetry = symmetry,
                                                                    cpMorphing = (meshMotion.lower() == "cpmorphing") )

        #sixDofDomainBody, flexProperties
        if nModesToUse > 0 :
            sixDofDomainBody = SixDofDomainBody.Build(case    = case,
                                                      mass    = mass,
                                                      inertia = inertia,
                                                      COG     = COG,
                                                      nModes  = nModesToUse,
                                                      donName = donName )
        else:
            sixDofDomainBody=None


        if nModesToUse > 0:
            #initFlexDict
            initFlexDict = InitFlexDict.Build(case      = case,
                                              mdFile    = mdFile,
                                              modes2use = modesToUse,
                                              datFile   = datFile,
                                              dmigFile  = dmigFile,
                                              draft     = draft,
                                              scale     = hmrScaling,
                                              vtkOut    = vtkOut,
                                              hullPatch = hullPatch,
                                              localPts  = localMotionPts )

            #flexFile
            flexFile = FlexFile.Build(case      = case,
                                      donName   = donName,
                                      freq      = shipFreq,
                                      damping   = shipDamping,
                                      dmigFile  = dmigFile,
                                      hullPatch = hullPatch,
                                      localPts  = localMotionPts )
        else:
            initFlexDict = None
            flexFile = None

        res =  cls( case, nProcs                    = nProcs,
                          OFversion                 = OFversion,
                          controlDict               = controlDict,
                          turbulenceModel           = turbulenceModel,
                          fvSchemes                 = fvSchemes,
                          fvSolution                = fvSolution,
                          dynamicMeshDict           = dynamicMeshDict,
                          transportProperties       = transportProperties,
                          decomposeParDict          = decomposeParDict,
                          waveProperties            = waveProperties,
                          sixDofDomainBody          = sixDofDomainBody,
                          initFlexDict              = initFlexDict,
                          flexFile                  = flexFile,
                          boundaryAlpha             = boundaryAlpha,
                          boundaryVelocity          = boundaryVelocity,
                          boundaryPressure          = boundaryPressure,
                          boundaryPointDisplacement = boundaryPointDisplacement,
                          setSelections             = setSelections,
                          clean                     = clean,
                          application               = application )

        res.meshMotion = meshMotion
        res.hullPatch = hullPatch
        res.symmetry = symmetry
        res.meshDir = meshDir
        res.stlFile = stlFile
        res.donFile = donFile
        res.EulerCellsDist = EulerCellsDist
        res.nRelaxZones = len(relaxZones)
        res.nModesToUse = nModesToUse

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
        idMerge = []
        nFacesMerge = 0
        #read and modify entries
        for i in range(nbound):
            if self.symmetry > 0 and (boundDict[2*i] in ['domainY0']):
                if self.symmetry == 1 : boundDict[2*i+1]['type'] = 'symmetry'
                elif self.symmetry == 2 : boundDict[2*i+1]['type'] = 'symmetryPlane'
            elif self.hullPatch in boundDict[2*i]:
                idMerge.append(2*i)
                nFacesMerge += boundDict[2*i+1]['nFaces']

        #merge all hull patches
        boundDict[idMerge[0]] = self.hullPatch
        boundDict[idMerge[0]+1]['nFaces'] = nFacesMerge
        for id in idMerge[:0:-1]:
            del boundDict[id+1]
            del boundDict[id]

        boundDict.writeFile()

    def writeAllinit(self, batchName="Allinit"):
        """Write bash script 'Allinit' to initialize case
        """
        print('Create init script')
        ainit = os.path.join(self.case,'Allinit')
        with open(ainit,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')

            #copy DON
            if self.donFile is not None:
                f.write('function copyDON()\n')
                f.write('{\n')
                fin  = self.donFile
                fout = self.case
                f.write('    cp {:s} {:s}\n'.format(fin,fout))
                f.write('}\n\n')

            #set Euler blending zone
            if self.EulerCellsDist is not None:
                #copy STL
                f.write('function copySTL()\n')
                f.write('{\n')
                fin  = os.path.join(os.getcwd(),self.meshDir,'constant','triSurface',self.stlFile)
                fout = os.path.join(os.getcwd(),self.case,'constant','triSurface')
                f.write('    mkdir constant/triSurface \n'.format(fin,fout))
                f.write('    cp -r {:s} {:s}\n'.format(fin,fout))
                f.write('}\n\n')

                f.write('function setEuler()\n')
                f.write('{\n')
                f.write('    copySTL\n')
                f.write('    setSet -batch {}\n'.format(os.path.join('system','setSet.euler')))
                f.write('    setsToZones -noFlipMap\n')
                f.write('}\n\n')

            #set Relax blending zone
            if self.nRelaxZones > 0:
                f.write('function setRelax()\n')
                f.write('{\n')
                f.write('    setSet -batch {}\n'.format(os.path.join('system','setSet.relax')))
                f.write('    setsToZones -noFlipMap\n')
                f.write('}\n\n')

            # __main__
            f.write('(\n')
            if self.donFile is not None: f.write('    copyDON\n')
            if self.EulerCellsDist is not None: f.write('    setEuler\n')
            if self.nRelaxZones > 0: f.write('    setRelax\n')
            f.write('    cp -rf ./0/org/* ./0/\n')
            if self.nModesToUse > 0: f.write('    initFlx initFlexDict\n')
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

    def runInit(self):
        import subprocess
        print("Run initialization")
        os.chdir(self.case)
        subprocess.call(["/bin/bash", "Allinit"], shell=False)

    def run(self):
        import subprocess
        print("Run")
        os.chdir(self.case)
        subprocess.call(["/bin/bash", "Allrun"], shell=False)
