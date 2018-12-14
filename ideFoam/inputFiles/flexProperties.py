from os.path import join
from ideFoam.inputFiles import ReadWriteFile, getFilePath
from PyFoam.Basics.DataStructures import Vector, SymmTensor, DictProxy

"""
  Convenience class to simply write flexible properties for sea-keeping case
"""

class InitFlexDict(ReadWriteFile) :
    """InitFlexDict file
    """
    
    @classmethod
    def Build(cls , case, mdFile=None, modes2use=None, datFile=None, dmigFile=None, draft=0., scale=1., vtkOut=True, hullPatch=None, localPts=None):
        
        res = cls( name = join(case, getFilePath("initFlexDict")), read = False )

        res.header["class"] = "dictionary"
        
        fem = DictProxy()
        fem["mdFile"] = '"../{}"; selected ( '.format(mdFile)+len(modes2use)*'{} '.format(*modes2use)+')'
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
        
        res["FEM_STRUCTURALMESH_VTU"] = fem
        
        return res
        
class FlexFile(ReadWriteFile) :
    """*.flex file
    """
    
    @classmethod
    def Build(cls , case, donName, freq={}, damping=[], dmigFile=None, hullPatch=None, localPts=None):
        
        res = cls( name = join(case, "constant", donName+".flex"), read = False )
        
        res["wn2"] = "({})".format(freq)
        res["dampingRatio"] = "({})".format(damping)
        dmig_name = ((dmigFile).rpartition('/')[-1]).partition('_dmig')[0]
        res['#include "modeshapes/{}_dmig.prj"'.format(dmig_name)] = ""
        res['#include "modeshapes/CFD_{}_fpt.flx"'.format(hullPatch)] = ""
        res['#include "modeshapes/CFD_{}_fps.flx"'.format(hullPatch)] = ""
        res['#include "modeshapes/CFD_{}_mpt.flx"'.format(hullPatch)] = ""
        if localPts is not None:
            if len(localPts)>0:
                res['#include "modeshapes/PTS_localMotion.flx"'] =  ""
        
        return res

# def writeFlexProperties(case, donName, freq=[], damping=[], mdFile=None, modes2use=None, datFile=None, dmigFile=None, draft=0., scale=1., vtkOut=True, hullPatch=None, localPts=None) :

    # a = InitFlexDict.Build( case, mdFile, modes2use, datFile, dmigFile, draft, scale, vtkOut, hullPatch, localPts)
    # a.writeFile()
    
    # a = FlexFile.Build( case, donName, freq, damping, dmigFile, hullPatch, localPts)
    # a.writeFile()

if __name__ == "__main__" :

   print(InitFlexDict("test"))
   print(FlexFile("test",donName='test.don'))

