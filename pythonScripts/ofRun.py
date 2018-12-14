from ofCase import OfCase
from inputFiles.gravity import Gravity
from inputFiles.turbulenceProperties import TurbulenceProperties, RASProperties

class OfRun( OfCase ):

    additionalFiles =   [
                        "dynamicMeshDict",
                        "transportProperties",
                        "waveProperties"
                        ]
    handledFiles = OfCase.handledFiles + additionalFiles
    
    #Always written
    # "turbulenceProperties"
    # "gravity"
    # "rasProperties"
    
    def __init__(self, *args, dynamicMeshDict       = None,
                              transportProperties   = None,
                              waveProperties        = None,
                              setSelections         = [],
                              turbulenceModel       = "laminar",
                              gravity               = 9.81,
                              **kwargs):
        """
        Same arguments as OfCase + some boundary input files
        """
        
        OfCase.__init__(self, *args, **kwargs)
        
        self.dynamicMeshDict = dynamicMeshDict
        self.transportProperties = transportProperties
        self.waveProperties = waveProperties
        self.turbulenceModel = turbulenceModel
        
        self.setSelections = setSelections
        self.gravity = gravity
        
    def writeFiles(self):
        OfCase.writeFiles(self)
        
        gravity = Gravity( case = self.case, g = self.gravity )
        gravity.writeFile()
        
        self.turbulenceProperties = TurbulenceProperties.Build(self.case, self.turbulenceModel )
        self.turbulenceProperties.writeFile()
        
        self.RASProperties = RASProperties.Build(self.case , self.turbulenceModel )
        self.RASProperties.writeFile()
        
        for i in self.setSelections: 
            i.writeFile()
        