
import math as mt
from os.path import join

"""
  Convenience class to simply write set selection input file
"""


class SetSelection() :
    """
       RefineMeshDict dictionary
    """
    def __init__(self , case, selType='box', BB=None, stlFile=None, opts=None, distance=None, outsidePoints=None, name=''):
        #set selection parameters
        if selType=='box':
            BBtxt = '({} {} {}) ({} {} {})'.format(*BB)
            self.cmd = 'cellSet c0 new boxToCell '+BBtxt
        elif selType=='proximity':
            stlFile = './constant/triSurface/'+stlFile
            if not stlFile.endswith('.stl'): stlFile += '.stl'
            if outsidePoints==None:
                tmp = findBoundingBox(stlFile, False)
                tmp = [0.5*(tmp[0]+tmp[3]), 0.5*(tmp[1]+tmp[4])+mt.fabs(tmp[1]-tmp[4]), 0.5*(tmp[2]+tmp[5])]
                outsidePoints = " (({} {} {}))".format(*tmp)
            else:
                outsidePoints = " (({} {} {}))".format(*outsidePoints)
            includeCutCells=' yes'
            includeInside=' yes'
            includeOutside=' no'
            curvature=' -1e6'
            distance = " " + str(distance)
            self.cmd = 'cellSet c0 '+opts+' surfaceToCell "'+ stlFile +'"'+ outsidePoints + includeCutCells + includeInside + includeOutside + distance + curvature
            if BB is not None:
                BBtxt = '({} {} {}) ({} {} {})'.format(*BB)
                self.cmd += '\n' + 'cellSet c0 subset boxToCell '+BBtxt
        
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



