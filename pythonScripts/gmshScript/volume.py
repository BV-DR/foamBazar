from .misc import *
from .point import Point
from .line import Line
from .surface import Surface

'''
Volume(s) are defined using a list of connected surface(s). A surface loop is created
(normally) to define the orientation of the mesh generated in the volume 
'''
class Volume(object):
    _ID_NAME = '_VOL_ID'
    _DB_NAME = '_EXISTING_VOLS'
    def __init__(self, geom, surfList):
        def check(s):
            if geom is None: return s
            if isinstance(s, Surface):
                found,idx = exist(geom,s)
                if found: return idx
            else:
                if geom.get(Surface, s) is not None: return s
            raise RuntimeError("Volume: Surface not found: " + str(s))
            return None
        for s in surfList: assert isinstance(s, (Surface, int, long))
        sid = unique_and_keep_order([check(s) for s in surfList])
        if len(sid) < 3: return RuntimeError("Volume: need at least 3 surfaces")
        # check all line id(s) and reverse it if needed due to the point connectivity
        self.sid = sid;        
        return

    # for printing to terminal
    def __repr__(self):
        return "v("+remove_bracket(str(self.dataFromKey(self.key())))+")"

    # NOTE: for uniqueness the sorted idx is used as "key" in the database
    def key(self, master=False):
        keystr=remove_bracket(str(sorted(map(abs,self.sid)) + self.sid))
        if master:
            return remove_bracket(str(sorted(map(abs,self.sid))))
        return keystr

    # this is an alternative constructor which can be called directly as "Volume.fromkey(keystr)"
    @classmethod
    def fromkey(cls, keystr):
        return Volume(None, cls.dataFromKey(keystr))
    @classmethod
    def masterDBKeys(cls, geom):
        return getDB(geom,cls).keys()
    @staticmethod
    def dataFromKey(keystr):
        sid=[int(i) for i in keystr.split(',')]
        return sid[len(sid)/2:]

