import PyFoam
from PyFoam.RunDictionary.ParsedParameterFile import WriteParameterFile
from PyFoam.Basics.DataStructures import Vector, SymmTensor, DictProxy
from os.path import join
from compatOF import alpha, p_rgh , waveAlpha, waveVelocity, pointDisp, namePatch

"""
  Convenience class to simply write boudary condition for sea-keeping case
"""

class InitFlexDict(WriteParameterFile) :
    """
        InitFlexDict file
    """
    def __init__(self , case, mdFile=None, modes2use=None, datFile=None, dmigFile=None, draft=0., scale=1., vtkOut=True, hullPatch=None, localPts=None, version="foamStar"):
        
        WriteParameterFile.__init__(self,  name = join(case, "initFlexDict")  )
        self.header["class"] = "dictionary"
        
        fem = DictProxy()
        fem["mdFile"] = '"../{}"; selected ({})'.format(mdFile,modes2use)
        fem["datFile"] = '"../{}"'.format(datFile)
        fem["dmigMfile"] = '"../{}"'.format(dmigFile)
        fem["dmigKfile"] = '"../{}"'.format(dmigFile)
        fem["pchCoordinate"] = SymmTensor(*[0,0,-draft,0,0,0])
        fem["pchScaleMode"] = scale
        fem["pchLengthUnit"] = 1
        fem["pchMassUnit"] = 1
        if not vtkOut: fem["outputToVTK"] = "no"
        fem["patches"] = "({})".format(hullPatch)
        fem["ySym"] = "(true)"
        if localPts is not None:
            if len(localPts)>0:
                fem["pointList"] = "(localMotion)"
                fem["localMotion"] = [ Vector(*pt) for pt in localPts ]
        
        self["FEM_STRUCTURALMESH_VTU"] = fem
        
class FlexFile(WriteParameterFile) :
    """
        *.flex file
    """
    def __init__(self , case, donName, freq={}, damping=[], dmigFile=None, hullPatch=None, localPts=None, version="foamStar"):
        
        WriteParameterFile.__init__(self,  name = join(case, "constant", donName+".flex"), noHeader=True )
        
        self["wn2"] = "({})".format(freq)
        self["dampingRatio"] = "({})".format(damping)
        dmig_name = ((dmigFile).rpartition('/')[-1]).partition('_dmig')[0]
        self['#include "modeshapes/{}_dmig.prj"'.format(dmig_name)] = ""
        self['#include "modeshapes/CFD_{}_fpt.flx"'.format(hullPatch)] = ""
        self['#include "modeshapes/CFD_{}_fps.flx"'.format(hullPatch)] = ""
        self['#include "modeshapes/CFD_{}_mpt.flx"'.format(hullPatch)] = ""
        if localPts is not None:
            if len(localPts)>0:
                self['#include "modeshapes/PTS_localMotion.flx"'] =  ""

def writeFlexProperties(case, donName, version, freq=[], damping=[], mdFile=None, modes2use=None, datFile=None, dmigFile=None, draft=0., scale=1., vtkOut=True, hullPatch=None, localPts=None) :

    a = InitFlexDict( case, mdFile, modes2use, datFile, dmigFile, draft, scale, vtkOut, hullPatch, localPts, version="foamStar")
    a.writeFile()
    
    a = FlexFile( case, donName, freq, damping, dmigFile, hullPatch, localPts, version="foamStar")
    a.writeFile()

if __name__ == "__main__" :

   print(BoundaryOmega("test", wallFunction = True))

