import re
from io import StringIO
import pandas as pd
import gzip


def readPoints(pointsFile):
    """Read points from a polymesh
    """
    with gzip.open(pointsFile, "rb") as f : 
        data = f.read().decode()
    pattern = r"\(\n(.*)\)\n\)"
    s = StringIO( re.search(pattern, data, re.DOTALL).group(1).replace("(", "").replace(")", ""))
    return pd.read_csv(s, delim_whitespace = True, header = None, names = ["x", "y", "z"] )


def getBounds(pointsFile):
    
    points = readPoints(pointsFile)
    
    return ( (float(points.x.min()) , float(points.x.max())),
             (float(points.y.min()) , float(points.y.max())),
             (float(points.z.min()) , float(points.z.max())),
           )
