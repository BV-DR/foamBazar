import os, re
from io import StringIO
import pandas as pd
import gzip


def readPoints(polyMeshDir):
    """Read points from a polymesh
    """
    pointsFile = os.path.join(polyMeshDir,'points')
    pointsFileGZ = os.path.join(polyMeshDir,'points.gz')
    #points in compressed format (*.gz)
    if os.path.isfile(pointsFileGZ):
        with gzip.open(pointsFileGZ, "rb") as f : 
            data = f.read().decode()
    #points in ascii format
    elif os.path.isfile(pointsFile):
        with open(pointsFile,'r') as f:
            data = f.read()
    else:
        raise(FileNotFoundError)
    pattern = r"\(\n(.*)\)\n\)"
    s = StringIO( re.search(pattern, data, re.DOTALL).group(1).replace("(", "").replace(")", ""))
    return pd.read_csv(s, delim_whitespace = True, header = None, names = ["x", "y", "z"] )


def getBounds(pointsFile):
    
    points = readPoints(pointsFile)
    
    return ( (float(points.x.min()) , float(points.x.max())),
             (float(points.y.min()) , float(points.y.max())),
             (float(points.z.min()) , float(points.z.max())),
           )
