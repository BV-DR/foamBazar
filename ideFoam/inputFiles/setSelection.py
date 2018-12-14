import math as mt
from os.path import join
from pythonScripts.fsTools import findBoundingBox

"""
  Convenience class to simply write set selection input file
"""

class SetSelection() :
    """SetSelection dictionary
    """
    def __init__(self , case, selType='box', BB=None, stlFile=None, opts=None, distance=None, outsidePoints=None, name='', selName='c0'):
        #set selection parameters
        if selType=='box':
            self.cmd = 'cellSet {} new boxToCell ({} {} {}) ({} {} {})'.format(selName,*BB)
        elif selType=='proximity':
            stlFile = './constant/triSurface/'+stlFile
            if not stlFile.endswith('.stl'): stlFile += '.stl'
            if outsidePoints==None:
                tmp = findBoundingBox(stlFile, False)
                tmp = [0.5*(tmp[0]+tmp[3]), 0.5*(tmp[1]+tmp[4])+mt.fabs(tmp[1]-tmp[4]), 0.5*(tmp[2]+tmp[5])]
                outsidePoints = tmp
            includeCutCells=' yes'
            includeInside=' yes'
            includeOutside=' no'
            curvature=' -1e6'
            self.cmd = 'cellSet {} {} surfaceToCell "{}"  (({} {} {})) {} {} {} {} {}'.format(selName,opts,stlFile,*outsidePoints,includeCutCells,includeInside,includeOutside,distance,curvature)
            if BB is not None:
                self.cmd += '\n cellSet {} subset boxToCell ({} {} {}) ({} {} {})'.format(selName,*BB)
        
        #set command file name
        if case==None: self.path = '.'
        else: self.path = join(case,'system')
        self.name = 'setSet'
        if name is not None: self.name += '.'+name

    def writeFile(self):
        fname = join(self.path,self.name)
        with open(fname,'w') as f:
            f.write(self.cmd)

if __name__ == "__main__" : 
   print(SetSelection("test"))



