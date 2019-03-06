import os, re
from io import StringIO
import pandas as pd
import gzip
import numpy as np


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
        raise(FileNotFoundError( 1, "{} or {} does not exists".format(pointsFile,pointsFileGZ) ))
    pattern = r"\(\n(.*)\)\n\)"
    s = StringIO( re.search(pattern, data, re.DOTALL).group(1).replace("(", "").replace(")", ""))
    return pd.read_csv(s, delim_whitespace = True, header = None, names = ["x", "y", "z"] )


def readPointsBin( polyMeshDir ):
    """Read points from a polymesh, binary format
    """

    pointsFile = os.path.join(polyMeshDir,'points')
    with open(pointsFile , "rb") as f:
        data = f.read( 2000 )

    header = data.split(b"(")[0]
    nbPoints = int(header.split(b"\n")[-2].decode())

    with open(pointsFile , "rb") as f:
        data = f.read(len(header)+1)
        next_ = f.read()

    xyz = np.frombuffer( next_[ :nbPoints*8*3 ]).reshape(-1,3)

    return pd.DataFrame( data = xyz, columns = ["x","y","z"] )



def getBounds(polyMeshDir):

    #TODO replace try/except, by reading the "format" field in the file header
    try :
        points = readPoints(polyMeshDir)
    except :
        points = readPointsBin(polyMeshDir)

    return ( (float(points.x.min()) , float(points.x.max())),
             (float(points.y.min()) , float(points.y.max())),
             (float(points.z.min()) , float(points.z.max())),
           )
