#!/usr/bin/env python

#########################################################################
# Filename: seakeepingMesher.py                                         #
# Date:     2018-Dec-06                                                 #
# Version:  1.                                                          #
# Author:   Alexis Benhamou                                             #
# Org.:     Bureau Veritas, (HO, France)                                #
# Email:    alexis.benhamou@bureauveritas.com                           #
#########################################################################
#  This class can be used to generate a mesh for a wave test case i.e.  #
#                  with ship and 3D wave propagation.                   #
#########################################################################

import os
import shutil
import math as mt
import numpy as np
from fsTools import findBoundingBox, findSTLPatches, translateStl, rotateStl, simpleGrading, simpleGradingN

from ofMesher import OfMesher

from inputFiles.fvSchemes import FvSchemes
from inputFiles.fvSolution import FvSolution
from inputFiles.controlDict import ControlDict
from inputFiles.decomposeParDict import DecomposeParDict
from inputFiles.refineMeshDict import RefineMeshDict
from inputFiles.snappyHexMeshDict import SnappyHexMeshDict
from inputFiles.surfaceFeatureExtractDict import SurfaceFeatureExtractDict
from inputFiles.blockMeshDict import BlockMeshDict
from inputFiles.setSelection import SetSelection

class SeakeepingMesher( OfMesher ):
    """Class used to generate a CFD mesh for seakeeping cases.
    It should be called by a python script as presented in the following example
    
    Example
    -------
    >>> import os
    >>> #
    >>> # Import meshing routine from foamBazar
    >>> # WARNING : foamBazar/pythonScripts must be added to PYTHONPATH !
    >>> from seakeepingMesher import SeakeepingMesher
    >>> #
    >>> # Case name and additional optional parameters
    >>> case = "mesh"
    >>> myParams = {'nProcs'        : 12,
    >>>             'stlFiles'      : ['hull.stl','deck.stl','stern.stl'],
    >>>             'draft'         : 11.75,
    >>>             'domain'        : [-2.78,3.,0,2.,-0.855,0.45],
    >>>             'fsZone'        : [8,7],
    >>>             'fsCellHeight'  : 1,
    >>>             'refFS'         : False,
    >>>             'onLiger'       : True
    >>>             }
    >>> #
    >>> # Call routines for meshing here
    >>> mesh = SeakeepingMesher.BuildFromParams( case, **myParams )
    >>> fname = os.path.join(case,'log.input')
    >>> with open(fname,'w') as f: f.write(str(myParams))
    >>> mesh.runInit()
    
    """
    @classmethod
    def BuildFromParams(cls,    case,
                                nProcs           = 4,
                                stlFiles         = None,
                                stlName          = 'ship',
                                refCellBuffer       = 4,
                                domain           = [-3.0,2.5, -2.0,2.0, -1.5,0.5], 
                                side             = 'port',
                                LOA              = None,
                                draft            = None,
                                heading          = 180,
                                fsZone           = None,
                                fsCellHeight     = None,
                                fsCellRatio      = 4,
                                refBow           = False,
                                refBowLength     = None,
                                refStern         = False,
                                refSternLength   = None,
                                refFS            = True,
                                refBoxType       = 'wave',
                                refBoxData       = [3],
                                refBoxRatio       = 3,
                                refSurfExtra     = None,
                                shipBL           = [3, 1.3, 0.7, 0.7],
                                noLayers         = [],
                                application           = "snappyHexMesh",
                                OFversion        = 5,
                                onLiger          = False,
                                clean            = False
                                ):
        """Build mesh for CFD seakeeping mesh form a few parameters.
        
        Parameters
        ----------
        case : str
            Name of case to create
        nProcs : int, default 4
            Number of processors used to build the mesh
        stlFiles : list of stings
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
            
        fsZone : list of floats, default [0.02*LOA, 0.01*LOA]
            Free surface zone [fsZmin, fsZmax]
        fsCellHeight : float, default fsZmax/3
            Cell height within the free surface zone
        fsCellRatio : int, default 4
            Cell ratio length/height within the free surface zone
            
        refBow : boolean, default False
            Logical defining if bow area should be refined
        refBowLength : flaot
            Coefficient defining lengh of selection for bow refinement (relatively to LOA)
        refStern : boolean, default False
            Logical defining if stern area should be refined
        refSternLength : boolean, default False
            Coefficient defining lengh of selection for stern refinement (relatively to LOA)
        refFS: boolean, default True
            Additional vertical refinement in the free surface zone near the ship
            
        refBoxType : str, default 'wave'
            Type of refinement box(es) in the domain ('wave', 'kelvin' or 'both')
        refBoxData : list of int, default [3]
            Data for defining refinement boxes
            The first number is the number of boxes.
            Dimension of each box is given as 4 floating numbers.
            If only the number of boxes is given, each box will be automatically dimensioned according to "refBoxRatio"
            format: [number of boxes]
                    [number of boxes, Xmin, Xmax, Ymin, Ymax, Xmin,Xmax,Ymin,Ymax, ...]
        refBoxRatio, int, default 3
            Ratio smallest/largest refinement boxes
        refCellBuffer : int, default 4
            Number of cells between each refinement level
        refSurfExtra : str
            STL file of regions where extra refinement in bow/stern areas with strong curvatures will be applied.
            
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
        filename = stlName+'_tmp.stl'
        if os.path.isfile(filename): os.remove(filename)
        if len(stlFiles)==0: raise SystemExit('No STL file specified.')
        
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
            draft = mt.fabs(float(draft))
            print("Set draft:",draft)
            move = -shipBB[2] - draft
            translateStl(filename, [0.0,0.0,move], filename)
            shipBB[2] += move
            shipBB[5] += move

        if refSurfExtra is not None:
            nameOnly = os.path.basename(refSurfExtra)
            surfFile = "./constant/triSurface/"+nameOnly
            print("Create stl: "+surfFile)
            shutil.copyfile(refSurfExtra,surfFile)
            if draft is not None:
                translateStl(surfFile, [0.0,0.0,move], surfFile)
        
        if rotateStl(filename, heading, filename):
            shipBBRot = findBoundingBox(filename)
        else:
            shipBBRot = shipBB

        LOA = shipBB[3]-shipBB[0] if LOA==None else LOA
        fsZone = [0.02*LOA, 0.01*LOA] if fsZone==None else fsZone
        fsCellHeight = fsZone[1]/3 if fsCellHeight==None else fsCellHeight
        
        #set default refinement values
        if (refBow and refBowLength==None): refBowLength=0.20*LOA
        if (refStern and refSternLength==None): refSternLength=0.20*LOA
        
        if side=='port': domain[2]=0
        elif side=='starboard': domain[3]=0
        domain = [ round(LOA*val, 3) for val in domain ]
        
        locationInMesh = [0.5*(domain[1]+shipBBRot[3]), 0.5*(domain[3]+domain[2]), 0.5*(domain[5]+domain[4])]

        ###FORMER ROUTINE : foamCase_template
        print('Create system folder input files')
        #controlDict
        controlDict = ControlDict.Build(case              = case,
                                        application           = "snappyHexMesh",
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
        
        ###FORMER ROUTINE : createBlockMeshDict
        nRefBox = int(refBoxData[0])    # how many refinement box? minimum is 1
    
        # The inner cell size is defined by fsCellHeight
        dZ = [ 0 for i in range(nRefBox+2)]
        dZ[nRefBox+1] = round(fsCellHeight, 5)
        for i in range(nRefBox, -1, -1):
            dZ[i] = dZ[i+1]*2  

        # keep is in a reverse order: i.e. level [n+1, n, n-1, n-2,..., 0]
        dZ = dZ[::-1]
        zCellSize = list(dZ)

        # horizontal cells for the whole domain    
        cellWidth = dZ[-2]*fsCellRatio
        Xcells = int((domain[1]-domain[0])/cellWidth)
        Ycells = int((domain[3]-domain[2])/cellWidth)

        # update domain data
        domain[0]= round(domain[1] - cellWidth*Xcells, 3)
        if side=='port':
            domain[3] = round(domain[2] + cellWidth*Ycells, 3)
        elif side=='starboard':
            domain[2] = round(domain[3] - cellWidth*Ycells, 3)
        elif side=='both':
            if (Ycells % 2 != 0):
                Ycells += 1 # we need this to be even
            domain[2] = -round(cellWidth*Ycells/2, 3)
            domain[3] = round(cellWidth*Ycells/2, 3)
        else:
            print("\nUnknown parameters, side=",side)
            raise SystemExit('abort ...')
    
        # compute vertical thickness of each box
        # All boxes must be larger than the ship's bounding box
        fsCellTop = int(mt.ceil(round(mt.fabs(fsZone[1])/(fsCellHeight),1)))
        fsCellBottom = int(mt.ceil(round(mt.fabs(fsZone[0])/(fsCellHeight),1)))
        fsZmax = fsCellTop*fsCellHeight
        fsZmin = -fsCellBottom*fsCellHeight
        fsZone[1] = fsZmax
        fsZone[0] = fsZmin

        # z-cuts in lower block
        lowerCutNCells = refCellBuffer
        lowerCut = fsZone[0] - dZ[1]*lowerCutNCells;
        while (lowerCut>(shipBBRot[2]-refCellBuffer*dZ[1])) | (lowerCutNCells < refCellBuffer):
            lowerCut -= dZ[1]
            lowerCutNCells += 1
        lowerCutNCells = [lowerCutNCells]
        lowerCut = [lowerCut]   

        #for i in range(1, nRefBox):
        #    lowerCutNCells.append(refCellBuffer)
        #    lowerCut.append(lowerCut[i-1] - refCellBuffer*dZ[i+1])
    
        # z-cuts in upper block
        upperCutNCells = refCellBuffer
        upperCut = fsZone[1] + dZ[1]*upperCutNCells;
        while (upperCut<(shipBBRot[5]+refCellBuffer*dZ[1])) | (upperCutNCells < refCellBuffer):
            upperCut += dZ[1]
            upperCutNCells += 1
        upperCutNCells = [upperCutNCells]
        upperCut = [upperCut]

        #for i in range(1, nRefBox):
        #    upperCutNCells.append(refCellBuffer)
        #    upperCut.append(upperCut[i-1] + refCellBuffer*dZ[i+1])

        # compute grading at the top/bottom block(s)
        ZratioTop = cellWidth/dZ[-nRefBox-1]
        dx1 = dZ[-nRefBox-1]/mt.fabs(domain[5]-upperCut[-1])
        ZcellsTop = simpleGradingN(dx1, ZratioTop)
        ZratioBottom = cellWidth/dZ[-nRefBox-1]
        dx1 = dZ[-nRefBox-1]/mt.fabs(domain[4]-lowerCut[-1])
        ZcellsBottom = simpleGradingN(dx1, ZratioBottom)
        ZratioBottom = 1.0/ZratioBottom

        # collect all z-cut
        zAllCut = lowerCut[::-1] + [fsZone[0], 0, fsZone[1]] + upperCut + [domain[5]]
        zAllCutNCells = [ZcellsBottom] + lowerCutNCells[::-1] + [fsCellBottom, fsCellTop] + upperCutNCells + [ZcellsTop]
        zAllCutRatio = [ZratioBottom] + list(1 for i in range(len(lowerCut)*2+2)) + [ZratioTop]

        # compute vertical position of all grid points 
        zGrid = [domain[4]]
        for i in range(len(zAllCut)):
            grad = simpleGrading(zAllCutNCells[i], zAllCutRatio[i])
            L = zAllCut[i] - zGrid[-1]
            zGrid = zGrid + [zGrid[-1]+L*grad[val] for val in range(1,len(grad)-1)]
            zGrid.append(zAllCut[i])
        zGridDelta = [zGrid[i+1]-zGrid[i] for i in range(0,len(zGrid)-1)]
    
        # compute vertical extension for all refBox
        dx = cellWidth/2.0
        refBoxZdata = [zGrid[0]]*nRefBox + [zGrid[-1]]*nRefBox
        for i in range(nRefBox):
            for j in range(len(zGridDelta)):
                if zGrid[j] >= fsZone[0]:
                    break
                if zGridDelta[j] < dx*1.25:
                    refBoxZdata[i] = zGrid[j]
                    break
            for j in range(len(zGridDelta)):
                if zGrid[-j-1] <= fsZone[1]:
                    break
                if zGridDelta[-j-1] < dx*1.25:
                    refBoxZdata[-i-1] = zGrid[-j-1]
                    break
            dx /= 2.0
    
        print('Domain bounding box:')
        print("   ", [domain[0], domain[2], domain[4], domain[1], domain[3], domain[5]])
        
        ptype = {}
        ptype['X0'], ptype['X1'] = "patch", "patch"
        if side=='port': ptype['Y0'], ptype['Y1'] = "symmetryPlane", "patch"
        elif side=='starboard': ptype['Y0'], ptype['Y1'] = "patch", "symmetryPlane"
        else: ptype['Y0'], ptype['Y1'] = "patch", "patch"
        ptype['Z0'], ptype['Z1'] = "wall", "patch"
        
        nz = len(zAllCut)
        patches = []
        for ps in ['X0','X1','Y0','Y1','Z0','Z1']:
            patches.append("{} domain{}".format(ptype[ps],ps))
            vert = []
            for i in range(nz):
                if ps=='X0': vert.append([4*i+1, 4*(i+1)+1, 4*(i+1)+2, 4*i+2])
                if ps=='X1': vert.append([4*i+0, 4*(i+1)+0, 4*(i+1)+3, 4*i+3])
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
                                            zmax        = zAllCut,
                                            Xcells      = Xcells,
                                            Ycells      = Ycells,
                                            Zcells      = zAllCutNCells,
                                            Zgrading    = zAllCutRatio,
                                            createPatch = True,
                                            patches     = patches,
                                            OFversion   = OFversion)

        # compute x,y data for refBox
        if (len(refBoxData) == 1):
            grad = simpleGrading(nRefBox+1, refBoxRatio)
            xGradMin = [mt.fabs(val-1.0) for val in grad[::-1]]
            xGradMax = grad
            yGradMin = [mt.fabs(val-1.0) for val in grad[::-1]]
            yGradMax = grad
            #
            xCutMin = [(shipBBRot[0]-domain[0])*val for val in xGradMin]
            del xCutMin[0], xCutMin[-1]
            for i in range(len(xCutMin)):
                xCutMin[i] = cellWidth*int(round(xCutMin[i]/cellWidth))
            xCutMin = [domain[0] + val for val in xCutMin]
            #
            if side=='port':
                yCutMin = [-1e-3 for val in xCutMin]
            else:
                yCutMin = [(shipBBRot[1]-domain[2])*val for val in yGradMin]
                del yCutMin[0], yCutMin[-1]
                for i in range(len(yCutMin)):
                    yCutMin[i] = cellWidth*int(round(yCutMin[i]/cellWidth))
                yCutMin = [domain[2] + val for val in yCutMin]
            #
            if refBoxType=='wave':
                xCutMax = [domain[1]+1e-3 for val in xCutMin]
            elif refBoxType=='kelvin':
                print("not yet implemented")
            else:
                xCutMax = [(domain[1]-shipBBRot[3])*val for val in xGradMax]
                del xCutMax[0], xCutMax[-1]
                for i in range(len(xCutMax)):
                    xCutMax[i] = cellWidth*int(round(xCutMax[i]/cellWidth))
                xCutMax = [shipBBRot[3] + val for val in xCutMax]
            #
            if side=='starboard':
                yCutMax = [1e-3 for val in xCutMin]
            else:
                yCutMax = [(domain[3]-shipBBRot[4])*val for val in yGradMax]
                del yCutMax[0], yCutMax[-1]
                for i in range(len(yCutMax)):
                    yCutMax[i] = cellWidth*int(round(yCutMax[i]/cellWidth))
                yCutMax = [shipBBRot[4] + val for val in yCutMax]
                yCutMax = yCutMax[::-1]
            #
            # update refBoxData
            for i in range(len(xCutMin)):
                refBoxData.append(xCutMin[i])
                refBoxData.append(xCutMax[i])
                refBoxData.append(yCutMin[i])
                refBoxData.append(yCutMax[i])
                
                
        ###FORMER ROUTINE : createBackGroundMesh
        # how many refinement boxes ? minimum is 1 
        refBoxBB = []
        if len(refBoxData)>1:
            if nRefBox != (len(refBoxData)-1)/4.:
                raise SystemExit('Error: invalid data for refinement boxes, ', refBoxData)
            for i in range(nRefBox):
                refBoxBB.append([refBoxData[i*4+1], refBoxData[i*4+3], refBoxZdata[i], refBoxData[i*4+2], refBoxData[i*4+4], refBoxZdata[-i-1]])
        else:
            raise SystemExit('\nData for refinement box is missing.\nrefBoxData=[#n, #xmin,#xmax,#ymax,#ymax, #xmin,#xmax,#ymin,#ymax, ..., repeat n times]\nabort ...')

        # number of level(s) for uniform 'xy'-refinement
        nxy = int(mt.log(mt.fabs(float(fsCellRatio)), float(2))) - 1

        # refine free surface (using proximity method)
        distance = fsCellHeight*refCellBuffer*2.0*(nxy-1.0 + nRefBox-1.0 + float(refBow | refStern | refFS))
    
        lastInnerBox = refBoxBB[-1]
        if lastInnerBox[0] > (shipBBRot[0] - distance - cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox))):
            diff = lastInnerBox[0] - (shipBBRot[0] - distance - cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox)))
            for i in range(len(refBoxBB)):
                refBoxBB[i][0] -= diff
        if lastInnerBox[1] > (shipBBRot[1] - distance - cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox))):
            diff = lastInnerBox[1] - (shipBBRot[1] - distance - cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox)))
            for i in range(len(refBoxBB)):
                refBoxBB[i][1] -= diff
        if lastInnerBox[2] > (shipBBRot[2] - distance - zCellSize[1]*refCellBuffer):
            diff = lastInnerBox[2] - (shipBBRot[2] - distance - zCellSize[1]*refCellBuffer)
            for i in range(len(refBoxBB)):
                refBoxBB[i][2] -= diff
        if lastInnerBox[3] < (shipBBRot[3] + distance + cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox))):
            diff = -lastInnerBox[3] + shipBBRot[3] + distance + cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox))
            for i in range(len(refBoxBB)):
                refBoxBB[i][3] += diff
        if lastInnerBox[4] < (shipBBRot[4] + distance + cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox))):
            diff = -lastInnerBox[4] + shipBBRot[4] + distance + cellWidth*refCellBuffer/mt.pow(2.0, float(nRefBox))
            for i in range(len(refBoxBB)):
                refBoxBB[i][4] += diff
        if lastInnerBox[5] < (shipBBRot[5] + distance + zCellSize[1]*refCellBuffer):
            diff = -lastInnerBox[5] + shipBBRot[5] + distance + zCellSize[1]*refCellBuffer
            for i in range(len(refBoxBB)):
                refBoxBB[i][5] += diff
        
        refineMeshDicts = []
        setSelections = []
        
        for i, BB in enumerate(refBoxBB):
            tmp = BB[4]
            BB[4]=domain[3] # YmaxDomain
            setSelections.append(SetSelection(case    = case,
                                              selType = 'box',
                                              BB      = BB,
                                              name    = 'x_'+str(i)))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'x_'+str(i),
                                                        directions     = 'tan1',
                                                        useHexTopology = True,
                                                        geometricCut   = False))
            BB[4] = tmp
            pass
    
        for i, BB in enumerate(refBoxBB):
            setSelections.append(SetSelection(case    = case,
                                              selType = 'box',
                                              BB      = BB,
                                              name    = 'y_'+str(i)))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'y_'+str(i),
                                                        directions     = 'tan2',
                                                        useHexTopology = True,
                                                        geometricCut   = False))
            pass
            
        # this is a point outside ship.stl
        outsidePoints = [0.5*(shipBBRot[0]+shipBBRot[3]), 0.5*(shipBBRot[1]+shipBBRot[4])+mt.fabs(shipBBRot[1]-shipBBRot[4]), 0.5*(shipBBRot[2]+shipBBRot[5])]

        for i in range(nxy):
            setSelections.append(SetSelection(case          = case,
                                              selType       = 'proximity',
                                              stlFile       = stlName,
                                              opts          = 'new',
                                              distance      = distance,
                                              outsidePoints = outsidePoints,
                                              name          = 'xy_'+str(i)))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'xy_'+str(i),
                                                        directions     = 'tan1 tan2',
                                                        useHexTopology = True,
                                                        geometricCut   = False))
            distance *= 0.5

        BB = [-1e6,-1e6,fsZmin,1e6,1e6,fsZmax]
        setSelections.append(SetSelection(case          = case,
                                          selType       = 'proximity',
                                          BB            = BB,
                                          stlFile       = stlName,
                                          opts          = 'new',
                                          distance      = distance,
                                          outsidePoints = outsidePoints,
                                          name          = 'xy'))
        refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                    set            = 'c0',
                                                    name           = 'xy',
                                                    directions     = 'tan1 tan2',
                                                    useHexTopology = True,
                                                    geometricCut   = False))
        
        BB = [-1e6,-1e6,-1e6,1e6,1e6,fsZmin]
        setSelections.append(SetSelection(case          = case,
                                          selType       = 'proximity',
                                          BB            = BB,
                                          stlFile       = stlName,
                                          opts          = 'new',
                                          distance      = distance,
                                          outsidePoints = outsidePoints,
                                          name          = 'xyz1'))
        refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                    set            = 'c0',
                                                    name           = 'xyz1',
                                                    directions     = 'tan1 tan2 normal',
                                                    useHexTopology = True,
                                                    geometricCut   = False))
                                              
        BB = [-1e6,-1e6,fsZmax,1e6,1e6,1e6]
        setSelections.append(SetSelection(case          = case,
                                          selType       = 'proximity',
                                          BB            = BB,
                                          stlFile       = stlName,
                                          opts          = 'new',
                                          distance      = distance,
                                          outsidePoints = outsidePoints,
                                          name          = 'xyz2'))
        refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                    set            = 'c0',
                                                    name           = 'xyz2',
                                                    directions     = 'tan1 tan2 normal',
                                                    useHexTopology = True,
                                                    geometricCut   = False))

        # align cutting locations
        if not refBow: refBowLength = 0.2*(shipBBRot[3]-shipBBRot[0])
        refBowLength = shipBBRot[3]-refBowLength 
        refBowLength = domain[1]-mt.ceil((domain[1]-refBowLength)/fsCellHeight)*fsCellHeight
        if refBow:
            distance *= 0.5
            BB = [refBowLength,-1e6,-1e6,1e6,1e6,shipBBRot[5]-fsCellHeight*refCellBuffer]
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
    
        if refSurfExtra is not None:
            nameOnly = os.path.basename(refSurfExtra)
            BB = [refBow+0.5*fsCellHeight*refCellBuffer,-1e6,-1e6,1e6,1e6,shipBBRot[5]-1.5*fsCellHeight*refCellBuffer]
            setSelections.append(SetSelection(case          = case,
                                              selType       = 'proximity',
                                              BB            = BB,
                                              stlFile       = nameOnly,
                                              opts          = 'new',
                                              distance      = 0.5*distance,
                                              outsidePoints = outsidePoints,
                                              name          = 'xyz4'))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'xyz4',
                                                        directions     = 'tan1 tan2 normal',
                                                        useHexTopology = True,
                                                        geometricCut   = False))

        # align cutting locations
        if not refStern: refSternLength = 0.2*(shipBBRot[3]-shipBBRot[0])
        refSternLength = shipBBRot[0]+refSternLength
        refSternLength = domain[1]-mt.floor((domain[1]-refSternLength)/fsCellHeight)*fsCellHeight
        if refStern:
            if not refBow: distance *= 0.5
            BB = [-1e6,-1e6,-1e6,refSternLength,1e6,shipBBRot[5]-fsCellHeight*refCellBuffer]
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

        if refSurfExtra is not None:
            nameOnly = os.path.basename(refSurfExtra)
            BB = [-1e6,-1e6,-1e6,refSternLength-0.5*fsCellHeight*refCellBuffer,1e6,shipBBRot[5]-1.5*fsCellHeight*refCellBuffer]
            setSelections.append(SetSelection(case          = case,
                                              selType       = 'proximity',
                                              BB            = BB,
                                              stlFile       = nameOnly,
                                              opts          = 'new',
                                              distance      = 0.5*distance,
                                              outsidePoints = outsidePoints,
                                              name          = 'xyz6'))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'xyz6',
                                                        directions     = 'tan1 tan2 normal',
                                                        useHexTopology = True,
                                                        geometricCut   = False))
    
        if refFS & (refBow | refStern):
            BB = [refSternLength,-1e6,fsZmin,refBowLength,1e6,fsZmax]
            setSelections.append(SetSelection(case          = case,
                                              selType       = 'proximity',
                                              BB            = BB,
                                              stlFile       = stlName,
                                              opts          = 'new',
                                              distance      = distance,
                                              outsidePoints = outsidePoints,
                                              name          = 'z'))
            refineMeshDicts.append(RefineMeshDict.Build(case           = case,
                                                        set            = 'c0',
                                                        name           = 'z',
                                                        directions     = 'normal',
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
                          clean = clean
                          )
                          
        res.stlName = stlName
        res.nxy = nxy
        res.refBoxBB = refBoxBB
        res.refFS = refFS
        res.refBow = refBow
        res.refStern = refStern
        res.refSurfExtra = refSurfExtra
        res.OFversion = OFversion
        res.onLiger = onLiger
        
        res.writeFiles()
        return res
        
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
            fstlin = os.path.join(os.getcwd(),self.stlName+'_tmp.stl')
            fstlout = os.path.join(self.case,'constant','triSurface',self.stlName+'.stl')
            f.write('    mv {:s} {:s}\n'.format(fstlin,fstlout))
            f.write('}\n\n')

            #refineBox
            f.write('function refineBox()\n')
            f.write('{\n')
            for dir in ['x','y']:
                for i in range(len(self.refBoxBB)):
                    f.write('setSet -latestTime -batch "system/setSet.{}_{}"\n'.format(dir,i))
                    f.write('refineMesh -dict "system/refineMeshDict.{}_{}"\n'.format(dir,i))
            for i in range(self.nxy):
                f.write('setSet -latestTime -batch "system/setSet.xy_{}"\n'.format(i))
                f.write('refineMesh -dict "system/refineMeshDict.xy_{}"\n'.format(i))
            f.write('setSet -latestTime -batch "system/setSet.xy"\n')
            f.write('refineMesh -dict "system/refineMeshDict.xy"\n')
            f.write('setSet -latestTime -batch "system/setSet.xyz1"\n')
            f.write('refineMesh -dict "system/refineMeshDict.xyz1"\n')
            f.write('setSet -latestTime -batch "system/setSet.xyz2"\n')
            f.write('refineMesh -dict "system/refineMeshDict.xyz2"\n')
            if self.refBow:
                f.write('setSet -latestTime -batch "system/setSet.xyz3"\n')
                f.write('refineMesh -dict "system/refineMeshDict.xyz3"\n')
            if self.refSurfExtra is not None:
                f.write('setSet -latestTime -batch "system/setSet.xyz4"\n')
                f.write('refineMesh -dict "system/refineMeshDict.xyz4"\n')
            if self.refStern:
                f.write('setSet -latestTime -batch "system/setSet.xyz5"\n')
                f.write('refineMesh -dict "system/refineMeshDict.xyz5"\n')
            if self.refSurfExtra is not None:
                f.write('setSet -latestTime -batch "system/setSet.xyz6"\n')
                f.write('refineMesh -dict "system/refineMeshDict.xyz6"\n')
            if self.refFS & (self.refBow | self.refStern):
                f.write('setSet -latestTime -batch "system/setSet.z"\n')
                f.write('refineMesh -dict "system/refineMeshDict.z"\n')
            f.write('}\n\n')
            
            #snap
            f.write('function snap()\n')
            f.write('{\n')
            if self.nProcs>1:
                f.write('    decomposePar -force -latestTime\n')
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
