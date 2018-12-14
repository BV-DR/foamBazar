from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import DictProxy
import os


class ReadWriteFile(ParsedParameterFile):
    """A specialization that is used to only write to the file"""
    def __init__(self,
                 name,
                 backup=False,
                 className="dictionary",
                 objectName=None,
                 createZipped=False,
                 read = True,  #If True read, if false create empty file
                 **kwargs):
        
        ParsedParameterFile.__init__(self,
                                     name,
                                     backup=backup,
                                     dontRead=not read,
                                     createZipped=createZipped,
                                     **kwargs)
        if not read : 
            if objectName==None:
                objectName=os.path.basename(name)
    
            self.content=DictProxy()
            self.header={"version":"2.0",
                         "format":"ascii",
                         "class":className,
                         "object":objectName}