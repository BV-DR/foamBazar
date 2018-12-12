from ofCase import OfCase

class OfMesher( OfCase ):

    
    additionalFiles = ["blockMeshDict", "snappyHexMeshDict", "extrudeMeshDict", "surfaceFeatureExtractDict"]
    handledFiles = OfCase.handledFiles + additionalFiles
        
    def __init__(self,   *args, extrudeMeshDict=None,
                                snappyHexMeshDict=None,
                                surfaceFeatureExtractDict=None,
                                blockMeshDict=None,
                                refineMeshDicts=[],
                                setSelections=[],
                                **kwargs):
        """
        Same arguments as OfCase + some boundary input files
        """

        OfCase.__init__(self, *args, **kwargs)
        
        self.blockMeshDict = blockMeshDict
        self.snappyHexMeshDict = snappyHexMeshDict
        self.surfaceFeatureExtractDict = surfaceFeatureExtractDict
        self.extrudeMeshDict = extrudeMeshDict
        
        self.refineMeshDicts = refineMeshDicts
        self.setSelections = setSelections
        
    def writeFiles(self):
        OfCase.writeFiles(self)
        if self.blockMeshDict : self.blockMeshDict.writeFile()
        if self.snappyHexMeshDict : self.snappyHexMeshDict.writeFile()
        if self.extrudeMeshDict : self.extrudeMeshDict.writeFile()
        if self.surfaceFeatureExtractDict : self.surfaceFeatureExtractDict.writeFile()
        
        
        for i in self.refineMeshDicts : 
            i.writeFile()
            
        for i in self.setSelections : 
            i.writeFile()
        
        
        
        
        
