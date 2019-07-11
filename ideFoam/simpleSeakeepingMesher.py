#!/usr/bin/env python

#########################################################################
#  This class can be used to generate a mesh for a wave test case i.e.  #
#                  with ship and 3D wave propagation.                   #
#                       Uniform backgroung mesh                         #
#########################################################################

import os
import shutil
import math as mt
from pythonScripts.fsTools import findBoundingBox, findSTLPatches, translateStl, rotateStl, simpleGrading, simpleGradingN

from ideFoam.ofMesher import OfMesher
from ideFoam.inputFiles import ControlDict, FvSchemes, FvSolution, DecomposeParDict,ExtrudeMeshDict
from ideFoam.inputFiles import BlockMeshDict, SetSelection, RefineMeshDict, SurfaceFeatureExtractDict, SnappyHexMeshDict

from scipy.optimize import fsolve

class simpleSeakeepingMesher( OfMesher ):
    """Class used to generate a CFD mesh for seakeeping cases.
    It should be called by a python script as presented in the following example

    Example
    -------
    >>> import os
    >>> #
    >>> # Import meshing routine from foamBazar
    >>> # WARNING : foamBazar must be added to PYTHONPATH !
    >>> from ideFoam.seakeepingMesher import SeakeepingMesher
    >>> #
    >>> # Case name and additional optional parameters
    >>> case = "mesh"      
    >>> waveLength=150
    >>> Lpp=355
    >>> d=waveLength/Lpp
    >>> myParams = {'nProcs'        : 8,
    >>>             'stlFiles'      : ['DTC-fullScale.stl'],
    >>>             'draft'         : 14.5,
    >>>             'domain'        : [(-2*d),1.0+d,0,1.5*d,-0.6,0.4],
    >>>             'fsZone'        : [-1.5*H,H],
    >>>             'waveLength'    : waveLength,
    >>>             'waveHeight'    : H,  
    >>>             'nxPerWaveLength': 100,
    >>>             'nzPerWaveHeight': 10,
    >>>             'nzBelowFS'      : 40,
    >>>             'nzAboveFS'      : 20,
    >>>             'onLiger'       : False
    >>>             }
    >>> #
    >>> # Call routines for meshing here
    >>> mesh = simpleSeakeepingMesher.BuildFromParams( case, **myParams )
    >>> fname = os.path.join(case,'log.input')
    >>> with open(fname,'w') as f: f.write(str(myParams))
    >>> mesh.runInit()

    """
    @classmethod
    def BuildFromParams(cls,    case,
                                nProcs           = 4,
                                stlFiles         = None,
                                stlName          = 'ship',
                                waveLength       = 300,
                                waveHeight       = 10,
                                LOA              = None,                                
                                domain           = [-3.0,2.5, -2.0,2.0, -1.5,0.5],
                                fsZone           = [-15., 10.],
                                side             = 'port',
                                draft            = None,
                                heading          = 180,
                                nxPerWaveLength  = 100,
                                nzPerWaveHeight  = 10,
                                nzBelowFS        = 40,
                                nzAboveFS        = 30,
                                                             
                                refBow           = False,
                                refBowLength     = None,
                                refStern         = False,
                                refSternLength   = None,

                                refBox4Snap      = 'False',

                                shipBL           = [3, 1.3, 0.7, 0.7],
                                noLayers         = [],
                                application      = "snappyHexMesh",
                                OFversion        = 5,
                                onLiger          = False,
                                clean            = "i"
                                ):
        """Build mesh for CFD seakeeping mesh form a few parameters.

        Parameters
        ----------
        case : str
            Name of case to create
        nProcs : int, default 4
            Number of processors used to build the mesh
        stlFiles : list of strings
            List of STl files names to be used for snapping
        stlName : str, default 'ship'
            Name used for final merged STL patch
        domain : list of floats, default [-3.0,2.5,-2.0,2.0,-1.5,0.5]
            Coeficients defining overall domain size (relative to LOA) [Xmin, Xmax, Ymin, Ymax, Zmin, Zmax]
        side : str, default 'port'
            Side to mesh ('port', 'starboard' or 'both')
        LOA : float, default value computed from STL
            Characteristic length (computed from STL file by default)
        draft : float
            Draft (measured from keel)
                if defined, ship.stl will be moved into position
                if "None", no operation will be performed to ship.stl
        heading : float, default 180.
            Heading of ship (180 is head sea)

        fsZone : list of floats, default [-1.5*H, H]
            Free surface zone [fsZmin, fsZmax]


        refBow : boolean, default False
            Logical defining if bow area should be refined
        refBowLength : float
            Coefficient defining lengh of selection for bow refinement (relatively to LOA)
        refStern : boolean, default False
            Logical defining if stern area should be refined
        refSternLength : boolean, default False
            Coefficient defining lengh of selection for stern refinement (relatively to LOA)


        refBox4Snap: boolean
            Refinement in ship region to ensure better snapping

        shipBL : list of float, default [3, 1.3, 0.7, 0.7]
             Dimension for boundary layers (relative to cell size)
             format: [nLayers, layerGrowth, finalLayerThickness, minThicknessRatio]
        noLayers: list of str
            Disable layers creation on selected patches (e.g. deck)

        application : str, default 'snappyHexMesh'
            application to run
        OFversion : int or str, default 5
            OpenFOAM version
        onLiger : boolean, default False
            Logical defining if case is run on Liger cluster
        clean : boolean, default False
            Logical to force case overwrite

        """

        ###FORMER ROUTINE : UserInput
        # create ship.stl and compute bounding box(es)
        if not isinstance(stlFiles, list): stlFiles = [stlFiles]
        stlFiles = [os.path.abspath( f ) for f in stlFiles]

        cwd_save = os.getcwd()

        #Work in the case directory, to get temporary file there.
        simpleSeakeepingMesher.makeCaseFolder(case, clean)

        filename = os.path.join(case, stlName+'_tmp.stl')
        if os.path.isfile(filename): os.remove(filename)

        with open(filename, 'w') as outfile:
            for fstl in stlFiles:
                if not fstl.endswith('.stl'): fstl += '.stl'
                if not os.path.isfile(fstl): raise SystemExit('File not found : {}'.format(fstl))
                with open(fstl) as infile:
                    for line in infile:
                        outfile.write(line)

        shipPatches = findSTLPatches(filename)
        print("Found patches in stl: ", shipPatches)
        shipBB = findBoundingBox(filename)

        if draft is not None:
            draft = abs(float(draft))
            print("Set draft:",draft)
            move = -shipBB[2] - draft
            translateStl(filename, [0.0,0.0,move], filename)
            shipBB[2] += move
            shipBB[5] += move

#        if refSurfExtra is not None:
#            nameOnly = os.path.basename(refSurfExtra)
#            surfFile = "./constant/triSurface/"+nameOnly
#            print("Create stl: "+surfFile)
#            shutil.copyfile(refSurfExtra,surfFile)
#            if draft is not None:
#                translateStl(surfFile, [0.0,0.0,move], surfFile)

        if rotateStl(filename, heading, filename):
            shipBBRot = findBoundingBox(filename)
        else:
            shipBBRot = shipBB

        LOA = shipBB[3] - shipBB[0] if LOA is None else LOA
        fsZone = [-1.5*waveHeight, waveHeight] if fsZone is None else fsZone

        #set default refinement values
        if (refBow and refBowLength is None): refBowLength = 0.20 * LOA
        if (refStern and refSternLength is None): refSternLength = 0.20 * LOA

        if side == 'port': domain[2] = 0
        elif side == 'starboard': domain[3] = 0
        domain = [ round(LOA*val, 3) for val in domain ]

        locationInMesh = [0.5*(domain[1]+shipBBRot[3]), 0.5*(domain[3]+domain[2]), 0.5*(domain[5]+domain[4])]

        print('Create system folder input files')
        #controlDict
        controlDict = ControlDict.Build(case              = case,
                                        application       = "snappyHexMesh",
                                        endTime           = 100,
                                        deltaT            = 1,
                                        writeControl      = "timeStep",
                                        writeInterval     = 1,
                                        writeCompression  = "compressed",
                                        runTimeModifiable = "true")

        #fvSchemes
        fvSchemes = FvSchemes.Build(case        = case,
                                    simType     = "CrankNicolson",
                                    limitedGrad = True,
                                    orthogonalCorrection = "implicit")

        #fvSolution
        fvSolution = FvSolution.Build(case = case )

        #decomposeParDict
        decomposeParDict = DecomposeParDict.Build(case   = case,
                                                  nProcs = nProcs)



        # horizontal cells for the whole domain
        longitudinalCellWidth = waveLength/nxPerWaveLength
        print('longitudinalCellWidth', longitudinalCellWidth)
        Xcells = int((domain[1]-domain[0])/longitudinalCellWidth)
        Ycells = int((domain[3]-domain[2])/longitudinalCellWidth)

        # update domain data
        domain[0] = round(domain[1] - longitudinalCellWidth*Xcells, 3)
        if side == 'port':
            domain[3] = round(domain[2] + longitudinalCellWidth*Ycells, 3)
        elif side == 'starboard':
            domain[2] = round(domain[3] - longitudinalCellWidth*Ycells, 3)
        elif side == 'both':
            if (Ycells % 2 != 0):
                Ycells += 1 # we need this to be even
            domain[2] = -round(longitudinalCellWidth*Ycells/2, 3)
            domain[3] = round(longitudinalCellWidth*Ycells/2, 3)
        else:
            print("\nUnknown parameters, side=",side)
            raise SystemExit('abort ...')

 
    
       # vertical cells for the whole domain
        verticalCellWidthFS=waveHeight/nzPerWaveHeight
        nzFS=int((fsZone[1]-fsZone[0])/verticalCellWidthFS)
        
        u0=verticalCellWidthFS
        
        LBelow = fsZone[0]-domain[4]
        LAbove = domain[5]-fsZone[1]
            
        funcBelow = lambda q : LBelow - u0* ((1.0 - q**(nzBelowFS+1))/(1.0 - q)) 
        funcAbove = lambda q : LAbove - u0* ((1.0 - q**(nzAboveFS+1))/(1.0 - q)) 

        q_initial_guess = 1.1
        q_below = fsolve(funcBelow, q_initial_guess)
        grading_below=q_below[0]**nzBelowFS


        q_initial_guess = 1.1
        q_above = fsolve(funcAbove, q_initial_guess)
        grading_above=q_above[0]**nzAboveFS
                                
                                

        print('Domain bounding box:')
        print("   ", [domain[0], domain[2], domain[4], domain[1], domain[3], domain[5]])

        ptype = {}
        ptype['X0'], ptype['X1'] = "patch", "patch"
        if side=='port': ptype['Y0'], ptype['Y1'] = "symmetryPlane", "patch"
        elif side=='starboard': ptype['Y0'], ptype['Y1'] = "patch", "symmetryPlane"
        else: ptype['Y0'], ptype['Y1'] = "patch", "patch"
        ptype['Z0'], ptype['Z1'] = "wall", "patch"

        nz = 3
        patches = []
        for ps in ['X0','X1','Y0','Y1','Z0','Z1']:
            patches.append("{} domain{}".format(ptype[ps],ps))
            vert = []
            for i in range(nz):
                if ps=='X0': vert.append([4*i+0, 4*(i+1)+0, 4*(i+1)+3, 4*i+3])
                if ps=='X1': vert.append([4*i+1, 4*(i+1)+1, 4*(i+1)+2, 4*i+2])
                if ps=='Y0': vert.append([4*i+0, 4*(i+1)+0, 4*(i+1)+1, 4*i+1])
                if ps=='Y1': vert.append([4*i+3, 4*(i+1)+3, 4*(i+1)+2, 4*i+2])
            if ps=='Z0': vert.append([0, 1, 2, 3])
            if ps=='Z1': vert.append([4*nz+0, 4*nz+1, 4*nz+2, 4*nz+3])
            patches.append(vert)


        #Write blockMeshDict file
        blockMeshDict = BlockMeshDict.Build(case        = case,
                                            ndim        = 3,
                                            waveMesh    = True,
                                            xmin        = domain[0],
                                            xmax        = domain[1],
                                            ymin        = domain[2],
                                            ymax        = domain[3],
                                            zmin        = domain[4],
                                            zmax        = [fsZone[0],fsZone[1],domain[5]],
                                            Xcells      = Xcells,
                                            Ycells      = Ycells,
                                            Zcells      = [nzBelowFS, nzFS ,nzAboveFS],
                                            Zgrading    = [1/grading_below, 1, grading_above],
                                            createPatch = True,
                                            patches     = patches,
                                            OFversion   = OFversion)


        refineMeshDicts = []
        setSelections = []


        if (refBox4Snap == True):

             # this is a point outside ship.stl
             outsidePoints = [0.5*(shipBBRot[0]+shipBBRot[3]), 0.5*(shipBBRot[1]+shipBBRot[4])+abs(shipBBRot[1]-shipBBRot[4]), 0.5*(shipBBRot[2]+shipBBRot[5])]
            
            
             distance = 5 * longitudinalCellWidth
             #BB = [refBowLength,-1e6,-1e6,1e6,1e6,shipBBRot[5]+0.1*shipBBRot[5]]
             setSelections.append(SetSelection(case          = case,
                                               selType       = 'proximity',
#                                              BB            = BB,
                                               stlFile       = stlName,
                                               opts          = 'new',
                                               distance      = distance,
                                               outsidePoints = outsidePoints,
                                               name          = 'xyz'))
             
             refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                         set            = 'c0',
                                                         name           = 'xyz',
                                                         directions     = 'tan1 tan2 normal',
                                                         useHexTopology = True,
                                                         geometricCut   = False)     )       
            
            

        distance *= 0.5 
        # align cutting locations
        if not refBow: refBowLength = 0.2*(shipBBRot[3]-shipBBRot[0])
        refBowLength = shipBBRot[3]-refBowLength
        refBowLength = domain[1]-mt.ceil((domain[1]-refBowLength)/longitudinalCellWidth)*longitudinalCellWidth
        if refBow:
            distance *= 0.5
            BB = [refBowLength,shipBBRot[2],shipBBRot[4],shipBBRot[1],shipBBRot[3],shipBBRot[5]]
            setSelections.append(SetSelection(case          = case,
                                              selType       = 'proximity',
                                              BB            = BB,
                                              stlFile       = stlName,
                                              opts          = 'new',
                                              distance      = distance,
                                              outsidePoints = outsidePoints,
                                              name          = 'xyz3'))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'xyz3',
                                                        directions     = 'tan1 tan2 normal',
                                                        useHexTopology = True,
                                                        geometricCut   = False))


        # align cutting locations
        if not refStern: refSternLength = 0.2*(shipBBRot[3]-shipBBRot[0])
        refSternLength = shipBBRot[0]+refSternLength
        refSternLength = domain[1]-mt.floor((domain[1]-refSternLength)/longitudinalCellWidth)*longitudinalCellWidth
        if refStern:
            if not refBow: distance *= 0.5
            BB = [0.,shipBBRot[2],shipBBRot[4],refSternLength,shipBBRot[3],shipBBRot[5]]
            setSelections.append(SetSelection(case          = case,
                                              selType       = 'proximity',
                                              BB            = BB,
                                              stlFile       = stlName,
                                              opts          = 'new',
                                              distance      = distance,
                                              outsidePoints = outsidePoints,
                                              name          = 'xyz5'))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'xyz5',
                                                        directions     = 'tan1 tan2 normal',
                                                        useHexTopology = True,
                                                        geometricCut   = False))





        ###FORMER ROUTINE : createSnappyMesh
        #surfaceFeatureExtract
        surfaceFeatureExtractDict = SurfaceFeatureExtractDict.Build(case = case,
                                                                    stlname = stlName)

        #snappyHexMesh
        minThickness = float(shipBL[3])*shipBL[2]/(pow(shipBL[1],float(shipBL[0]-1)))
        snappyHexMeshDict = SnappyHexMeshDict.Build(case                       = case,
                                                    stlname                    = stlName,
                                                    castellatedMesh            = True,
                                                    snap                       = True,
                                                    addLayers                  = True,
                                                    relativeSizes              = True,
                                                    locationInMesh             = locationInMesh,
                                                    nCellsBetweenLevels        = 1,
                                                    edgeLvl                    = 0,
                                                    hullLvl                    = [0,0],
                                                    resolveFeatureAngle        = 15,
                                                    allowFreeStandingZoneFaces = False,
                                                    snapTol                    = 0.75,
                                                    nSolveIter                 = 100,
                                                    nSurfaceLayers             = shipBL[0],
                                                    expansionRatio             = shipBL[1],
                                                    finalLayerThickness        = shipBL[2],
                                                    minThickness               = minThickness,
                                                    featureAngle               = 60,
                                                    stlPatches                 = shipPatches,
                                                    noLayers                   = noLayers,
                                                    maxNonOrtho                = 65,
                                                    minTwist                   = 0.02,
                                                    nSmoothScale               = 5,
                                                    errorReduction             = 0.75,
                                                    OFversion                  = OFversion)

        res =  cls( case, nProcs=nProcs,
                          controlDict=controlDict,
                          fvSchemes=fvSchemes,
                          fvSolution=fvSolution,
                          decomposeParDict=decomposeParDict,
                          refineMeshDicts=refineMeshDicts,
                          snappyHexMeshDict=snappyHexMeshDict,
                          surfaceFeatureExtractDict=surfaceFeatureExtractDict,
                          blockMeshDict=blockMeshDict,
                          setSelections=setSelections,
                          application = application,
                          isMesher = True,
                          clean = False
                          )

        res.stlName = stlName
        res.refBox4Snap = refBox4Snap
#        res.nxy = nxy
#        res.refBoxBB = refBoxBB
        res.refFS = refFS
        res.refBow = refBow
        res.refStern = refStern
        res.refSurfExtra = refSurfExtra
        res.OFversion = OFversion
        res.onLiger = onLiger

        res.writeFiles()

        os.chdir(cwd_save)

        return res


    def get2DCase(self, case2D=None, step = 10):
        """Create a 2D mesh (extrusion of the y=0 slice)
        """
        if case2D is None :
            case2D = os.path.join(self.case, "wave2D")
        #Copy mesh in
        return ExtrudeMeshWave2D.BuildWave2D(case = case2D, sourceCase= self.case, step = step)

    def writeAllinit(self):
        """Write bash script 'Allinit' to create mesh
        """
        print('Create run scripts')
        ainit = os.path.join(self.case,'Allinit')
        with open(ainit,'w') as f:
            f.write('#! /bin/bash\n')
            f.write('set -x\n\n')

            #moveSTL
            f.write('function moveSTL()\n')
            f.write('{\n')
            fstlin = os.path.join(self.case,self.stlName+'_tmp.stl')
            fstlout = os.path.join(self.case,'constant','triSurface',self.stlName+'.stl')
            f.write('    mv {:s} {:s}\n'.format(fstlin,fstlout))
            f.write('}\n\n')

#            #decomp
            if self.nProcs>1:
                f.write('function decomp()\n')
                f.write('{\n') 
                f.write('    decomposePar -force -latestTime\n')
                f.write('}\n\n')


#            #refineBox
            f.write('function refineBox()\n')
            f.write('{\n') 
            if (self.refBox4Snap==True):            

                if self.nProcs>1:
                    f.write('    mpirun -np {} setSet -parallel -latestTime -batch "system/setSet.xyz"\n'.format(self.nProcs))
                    f.write('    mpirun -np {} refineMesh -parallel -dict "system/refineMeshDict.xyz"\n'.format(self.nProcs))
                else:
                    f.write('setSet -latestTime -batch "system/setSet.xyz"\n')
                    f.write('refineMesh -dict "system/refineMeshDict.xyz"\n')               
            if self.refBow:
                f.write('setSet -latestTime -batch "system/setSet.xyz3"\n')               
                if self.nProcs>1:
                    f.write('    mpirun -np {} refineMesh -parallel -dict "system/refineMeshDict.xyz3"\n'.format(self.nProcs))
                else:
                    f.write('refineMesh -dict "system/refineMeshDict.xyz3"\n')                   

            if self.refStern:
                f.write('setSet -latestTime -batch "system/setSet.xyz5"\n')
                if self.nProcs>1:
                    f.write('    mpirun -np {} refineMesh -parallel -dict "system/refineMeshDict.xyz5"\n'.format(self.nProcs))
                else:
                    f.write('refineMesh -dict "system/refineMeshDict.xyz5"\n')  
            f.write('}\n\n')

            #snap
            f.write('function snap()\n')
            f.write('{\n')
            if self.nProcs>1:
                #f.write('    decomposePar -force -latestTime\n')
                if self.onLiger:
                    self.writeSbatch()
                    f.write('    sbatch run.sh\n')
                else:
                    f.write('    mpirun -np {} snappyHexMesh -parallel\n'.format(self.nProcs))
                    f.write('    reconstructParMesh -latestTime\n')
            else:
                f.write('    snappyHexMesh\n')
            f.write('}\n\n')

            # __main__
            f.write('(\n')
            f.write('    moveSTL\n')
            f.write('    blockMesh\n')
            if self.nProcs>1:
                f.write('    decomp\n')                
            if (self.refBox4Snap==True or self.refBow==True or self.refStern==True):
                f.write('    refineBox\n')
            f.write('    surfaceFeatureExtract\n')
            f.write('    snap\n')
            f.write(') 2>&1 | tee log.mesh\n\n')
            if not (self.onLiger and self.nProcs>1): f.write('checkMesh -latestTime 2>&1 | tee log.checkMesh\n')
        os.chmod(ainit, 0o755)

    def writeAllclean(self):
        """Write bash script 'Allclean' to clean mesh folder(s)
        """
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
            f.write('#eval clean_log\n')
            f.write('eval clean_mesh\n')
            f.write('eval clean_parallel_mesh\n')
            f.write('eval clean_0\n')
        os.chmod(aclean, 0o755)




class ExtrudeMeshWave2D(OfMesher):

    def __init__(self, case, step=None, **kwargs):
        OfMesher.__init__(self, case, **kwargs)
        self.sourceCase = self.extrudeMeshDict["sourceCase"]
        self.step = step

    @classmethod
    def BuildWave2D(cls, case, sourceCase, step):
        """ Build 2D mesh from seakeeping mesh
        """
        extrudeMeshDict = ExtrudeMeshDict.Build( case = case, sourceCase = sourceCase, sourcePatch = "side1",  exposedPatchName = "side2")
        controlDict = ControlDict.Build(case = case, deltaT=1, endTime = 100)
        return cls( case, extrudeMeshDict = extrudeMeshDict, controlDict = controlDict, step = step )

    def writeAllinit(self):

        #Copy required step in step "1000"
        str_ = "#! /bin/bash\n"
        if self.step is not None :
            str_ += "cp -r {} {}\n".format( os.path.join(self.sourceCase, str(self.step) ), os.path.join(self.sourceCase, "1000" ) )
        str_ += "extrudeMesh"
        with open( os.path.join(self.case, "Allinit"), "w" ) as fil:
            fil.write(str_)






