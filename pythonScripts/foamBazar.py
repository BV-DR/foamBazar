import PyFoam
import os
import numpy as np
import fsStl
import fsApp
#import fsPlot

#__all__ = ['os',]

def foamVersion():
    '''
    return OpenFOAM version, ("name","version")
    e.g.
        ('foamStar', '5.x')
        ('OpenFOAM+', 'v1712')
        ('foam', '4.0')
    '''
    if not "WM_PROJECT_VERSION" in os.environ: return (None,None)
    name=os.environ["WM_PROJECT"]
    vers=os.environ["WM_PROJECT_VERSION"]
    from shutil import which
    if which('foamStar') is not None: name='foamStar'
    if vers.startswith('v'): name += '+'
    return (name,vers)

_OFVERSION=foamVersion()



