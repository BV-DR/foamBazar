import numpy as np
from .misc import *

'''
A point on its own has a position (pos) and a characteristic length (lc), and no idx. The idx will be generated when the point is added to geo-file.   
'''
class Point(object):
    _FORMAT_FLOAT = '%.6e'
    _FLOAT_TO_STR = '{:.6e}'
    _ID_NAME = '_POINT_ID'
    _DB_NAME = '_EXISTING_POINTS'
    def __init__(self, x, y=None, z=None, lc=None):
        p=[0,0,0]
        assert isinstance(x, (np.ndarray,list,int,long,float))
        if not isinstance(x, (np.ndarray,list)):
            assert isinstance(y, (int,long,float))
            assert isinstance(z, (int,long,float))
            p[0]=float(x)
            p[1]=float(y)
            p[2]=float(z)
        else:
            p[0]=float(x[0])
            p[1]=float(x[1])
            p[2]=float(x[2])            
        if lc is not None: assert isinstance(lc, (int,long,float))
        self.pos = np.array(p)
        self.lc = lc
        return

    # for printing to terminal
    def __repr__(self):
        return "p("+str(self.pos[0])+','+str(self.pos[1])+','+str(self.pos[2])+')'

    def code(self, geom):
        '''
        Return the code for use in the geo-file
        '''
        # we do not allow the same point to be added twice
        # self.exist(...) should return a (new) idx if not found 
        found,idx = exist(geom,self)
        if found: return ''
        fmt = ', '.join(len(self.pos) * [Point._FORMAT_FLOAT])
        if self.lc is not None:
            return '\n'.join([('Point(%d) = {' + fmt + ', ' + Point._FORMAT_FLOAT + '};') % ((idx,) + tuple(self.pos) + (self.lc,))])
        else:
            return '\n'.join([('Point(%d) = {' + fmt + '};') % ((idx,) + tuple(self.pos))])
        pass 
        
    def key(self, master=False):
        keystr = remove_bracket(np.array2string(self.pos, separator=',', prefix=')', formatter={'float_kind':lambda x: Point._FORMAT_FLOAT % x}))
        return keystr

    # this is an alternative constructor which can be called directly as "Point.fromkey(keystr)"
    @classmethod
    def fromkey(cls, keystr):
        pos=[float(i) for i in keystr.split(',')]
        return Point(pos)
    @classmethod
    def masterDBKeys(cls, geom):
        return getDB(geom,cls).keys()
        
