from ofCase import OfCase

class OfMesher( OfCase ):

    
    additionalFiles = ["blockMeshDict", "snappyHexMeshDict", "extrudeMeshDict", "surfaceFeatureExtractDict"]
    handledFiles = OfCase.handledFiles + additionalFiles
        
    def __init__(self,   *args, extrudeMeshDict=None,
                                snappyHexMeshDict=None,
                                surfaceFeatureExtractDict=None,
                                blockMeshDict=None,
                                refineMeshDicts=None,
                                setSelections=None,
                                **kwargs):
        """
        Same arguments as OfCase + some boundary input files
        """

        OfCase.__init__(self, *args, **kwargs)
        
        self.blockMeshDict = extrudeMeshDict
        self.snappyHexMeshDict = snappyHexMeshDict
        self.surfaceFeatureExtractDict = surfaceFeatureExtractDict
        self.extrudeMeshDict = extrudeMeshDict
        
        
        self.refineMeshDicts = refineMeshDicts
        self.setSelections = setSelections
        
        
        
