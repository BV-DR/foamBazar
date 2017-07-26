import copy
from .point import Point 
from .misc import *

'''
Line is defined using two point(s).
'''
class Line(object):
    _ID_NAME = '_LINE_ID'
    _DB_NAME = '_EXISTING_LINES'
    def __init__(self, geom, p0, p1):
        def check(p):
            if geom is None: return p
            if isinstance(p, Point):
                found,pid = exist(geom,p)
                if found: return pid
            else:
                if geom.get(Point,p) is not None: return p
            return None
        assert isinstance(p0, (Point, int, long))
        assert isinstance(p1, (Point, int, long))
        self.pid = [check(p0), check(p1)]
        if self.pid[0] is None: raise RuntimeError("Line: Point p0 does not exist in geo-file")
        if self.pid[1] is None: raise RuntimeError("Line: Point p1 does not exist in geo-file")
        if self.pid[0] == self.pid[1]: raise RuntimeError("Line: Cannot construct lines of zero length")
        return

    # for printing to terminal
    def __repr__(self):
        return "l("+remove_bracket(str(self.dataFromKey(self.key())))+")"

    def code(self, geom):
        '''
        Return the code for use in the geo-file
        '''
        # we do not allow the same line to be added twice
        # self.exist(...) should return a (new) idx if not found 
        found,idx = exist(geom,self)
        if found: return ''
        return '\n'.join([('Line(%d) = {%d,%d};') % (idx,self.pid[0], self.pid[1])])

    # NOTE: for uniqueness the sorted idx is used as "key" in the database
    def key(self, master=False):
        keystr=remove_bracket(str(sorted(map(abs,self.pid)) + self.pid))
        if master:
            return remove_bracket(str(sorted(map(abs,self.pid))))
        return keystr

    # this is an alternative constructor which can be called directly as "Line.fromkey(keystr)"
    @classmethod
    def fromkey(cls, keystr):
        pid=cls.dataFromKey(keystr)
        return Line(None, pid[0], pid[1])
    @classmethod
    def masterDBKeys(cls, geom):
        subkeys=copy.deepcopy(getDB(geom,cls).keys())
        for i in range(0,len(subkeys)):
            tmp=subkeys[i].split(',')
            subkeys[i]=",".join(tmp[:len(tmp)/2])
        return subkeys
    @staticmethod
    def dataFromKey(keystr):
        return [int(i) for i in keystr.split(',')][2:]
        
