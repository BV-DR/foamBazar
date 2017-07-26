import copy
from .misc import *
from .point import Point
from .line import Line

'''
Surface(s) are defined using a list of connected line(s). A line loop is created
(normally) with negative label if necessary to reverse the end points
'''
class Surface(object):
    _ID_NAME = '_SURF_ID'
    _DB_NAME = '_EXISTING_SURFS'
    def __init__(self, geom, lineList):
        def check(l):
            if geom is None: return l
            if isinstance(l, Line):
                found,pid = exist(geom,l)
                if found: return pid
            else:
                if geom.get(Line,l) is not None: return l
            raise RuntimeError("Surface: line not found: " + str(l))
            return None
        for i in lineList: assert isinstance(i, (Line, int, long))
        lid = unique_and_keep_order([check(i) for i in lineList])
        if len(lid) < 3: return RuntimeError("Surface: need at least 3 lines")        
        # check all line id(s) and reverse it if needed due to the point connectivity
        self.lid = [lid.pop(0)]
        if geom is None: # cannot check if the geo-data is not provided
            for i in lid: self.lid.append(i)
            return
        l = geom.get(Line, self.lid[0])
        p0 = l.pid[0] if self.lid[0]>0 else l.pid[1]
        plast = l.pid[1] if self.lid[0]>0 else l.pid[0]
        while len(lid):
            failed = True
            for i,val in enumerate(lid):
                l = geom.get(Line, val)
                if l.pid[0]==plast:
                    plast=l.pid[1]
                    self.lid.append(lid.pop(i))
                    failed = False
                    break
                elif l.pid[1]==plast:
                    plast=l.pid[0]
                    self.lid.append(-lid.pop(i))
                    failed = False
                    break
            if failed: raise RuntimeError("Surface: Lines are not connected")
        if not (plast == p0): raise RuntimeError("Surface: lines are not closed")
        return

    # for printing to terminal
    def __repr__(self):
        return "s("+remove_bracket(str(self.dataFromKey(self.key())))+")"

    # NOTE: for uniqueness the sorted idx is used as "key" in the database
    def key(self, master=False):
        keystr=remove_bracket(str(sorted(map(abs,self.lid)) + self.lid))
        if master:
            return remove_bracket(str(sorted(map(abs,self.lid))))
        return keystr

    # this is an alternative constructor which can be called directly as "Surface.fromkey(keystr)"
    @classmethod
    def fromkey(cls, keystr):
        return Surface(None, cls.dataFromKey(keystr))
    @classmethod
    def masterDBKeys(cls, geom):
        subkeys=copy.deepcopy(getDB(geom,cls).keys())
        for i in range(0,len(subkeys)):
            tmp=subkeys[i].split(',')
            subkeys[i]=",".join(tmp[:len(tmp)/2])
        return subkeys
    @staticmethod
    def dataFromKey(keystr):
        lid=[int(i) for i in keystr.split(',')]
        return lid[len(lid)/2:]

