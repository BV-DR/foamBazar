import numpy as np
import copy
from .misc import *
from .point import Point
from .line import Line
from .surface import Surface
from .volume import Volume

"""
    class to handle gmsh geo-file(s)
"""

class extdb(dict):
    '''
    Extrude database, this is for conveniently accessing dict-keys by calling as attribute
    '''
    def __getattr__(self, attr):
        return self[attr]
        
class geo(object):
    def __init__(self):
        '''
        GMSH requires users to provide unique ID(s) for point(s), line(s), etc.
        and we need to keep track of these ID(s) manually
        '''
        self.__dict__[Point._ID_NAME] = 0
        self.__dict__[Line._ID_NAME] = 0
        self.__dict__[Surface._ID_NAME] = 0
        self.__dict__[Volume._ID_NAME] = 0
        self.__dict__[Point._DB_NAME] = dict()
        self.__dict__[Line._DB_NAME] = dict()
        self.__dict__[Surface._DB_NAME] = dict()
        self.__dict__[Volume._DB_NAME] = dict()
        self._EXTRUDE_ID = 0
        self._PHYS_IDS = []     # array of physical group id(s)
        self._CODE = [
        '/* This script was generated using fsMesher.gmshScript */',
        'Geometry.OldNewReg=0;'
        ]
        return

    # for printing to terminal
    def __repr__(self):
        quickinfo = "geo(p:"+str(len(getDB(self,Point)))
        if self.hasDB(Line):
            quickinfo += ",l:" + str(len(getDB(self,Line)))
        if self.hasDB(Surface):
            quickinfo += ",s:" + str(len(getDB(self,Surface)))
        return quickinfo + ")"

    def printDB(self):
        if not self.hasDB(Point):
            print('no data')
            return
        self._print_db(getDB(self,Point), prefix='p')
        print('next p:', getIDX(self,Point) + 1)
        self._print_db(getDB(self,Line), prefix='l')
        print('next l:', getIDX(self,Line) + 1)
        self._print_db(getDB(self,Surface), prefix='s')
        print('next s:', getIDX(self,Surface) + 1)
        self._print_db(getDB(self,Volume), prefix='v')
        print('next v:', getIDX(self,Volume) + 1)
        print
        self.printScript()
        return

    def _print_db(self, db, prefix=''):
        idx = sorted(db, key=db.get)
        for i in idx:
            print(prefix + str(db[i]), ':', i)
        return

    def printScript(self):
        tmp = self._CODE
        for i in tmp:
            print(i)
        return

    def add(self, obj):
        '''
        Add a geometrical object to the code ... the actual code is generated in
        obj.code(self) where the arg. self is needed for a proper check of id(s) 
        '''
        obj_code = obj.code(self)
        if obj_code:
            self._CODE.append(obj_code)
            self._db_insert(obj)
        return
        
    def addPoint(self, x, y, z, lc=None):
        p = Point(x,y,z,lc)
        self.add(p)
        return p
    def addLine(self, p0, p1):
        l = Line(self,p0,p1)
        self.add(l)
        return l

    def extrude(self, obj, dx, dy, dz, layers=1, opts=None):
        '''
        Extrude "point|line|surface" along translation axis
        '''
        # we need the object in a list format
        objList = obj if isinstance(obj, list) else [obj]
        if len(objList) == 0 or objList[0] is None: return
        assert isinstance(dx, (int,long,float))
        assert isinstance(dy, (int,long,float))
        assert isinstance(dz, (int,long,float))
        assert isinstance(layers, (str,int,list,np.ndarray))

        #The layers are defined using two arrays i.e. Layers {{nElem[]},{nCut[]}}
        #The first array nElem[]={1,1,1,(n elements),1,1,1} defines the number of element created between each cut.
        #The second array nCut[]={0.1,0.2,(n cuts),...,1} defines the cut location (normalized) where the last cut must be at 100% i.e. 1
        layers_str='1'
        if isinstance(layers, (int, long)):
            layers_str=str(layers)
        elif isinstance(layers, str):
            # user(s) need to provide a valid format here
            # e.g: '#n' or '{n,n,n,n}, {float,float,float,1}'
            layers_str=layers
        elif isinstance(layers, (np.ndarray,list)):
            layerList = copy.deepcopy(layers)
            layerList.sort()
            maxVal = max(layerList) # for normalization
            # assume each cut has 1 element, and use only cut locations to control the extrude
            nElem_str = ','.join(str(1) for i in layerList)
            cut_str = ','.join(Point._FLOAT_TO_STR.format(float(i)/maxVal) for i in layerList)
            layers_str = '{' + nElem_str + '},{' + cut_str + '}'

        #
        # Scan the object list and determine the type
        # All element must be of the same type i.e. either Point|Line|Surface
        objtype = objList[0].__class__
        for i in objList:
            if not isinstance(i, objtype):
                raise RuntimeError("extrude: all extruded obj must be of the same type")
        #
        if isinstance(objList[0], Point):
            return self._extrude_points(objList, [dx,dy,dz], layers_str, opts=opts)
        elif isinstance(objList[0], Line):
            return self._extrude_lines(objList, [dx,dy,dz], layers_str, opts=opts)
        elif isinstance(objList[0], Surface):
            return self._extrude_surfaces(objList, [dx,dy,dz], layers_str, opts=opts)
        else:
            raise RuntimeError('The object to be extruded must be of type Point|Line|Surface')        
        return
    
    def hasDB(self,obj):
        return bool(getDB(self,obj))

    def incIDX(self,obj,n):
        self.__dict__[obj._ID_NAME] += n
        return

    def get(self, obj, idx):
        db=getDB(self,obj)
        allIdx=db.values()
        if not abs(idx) in allIdx: return None
        return obj.fromkey(db.keys()[allIdx.index(abs(idx))])

    def _create_idx_str(self, objList):
        idx = []
        for obj in objList:
            if not obj.key() in getDB(self,obj):
                raise RuntimeError('id not found: ' + str(obj))
            idx.append(getDB(self,obj)[obj.key()])
        return ','.join(str(i) for i in idx)

    def _db_insert(self, obj):
        found,idx = exist(self,obj)
        self.incIDX(obj,1)  # gmsh always keeps incrementing the id by 1 !!!
        if not found:
            getDB(self,obj)[obj.key()] = getIDX(self,obj)
            return True     # insert successful
        else:
            return False    # no need to insert, the obj already exists
        
    def _extrude_points(self, pointList, axis, layers, opts=None):
        '''
        line[] = Extrude{dx, dy, dz} { Point{#ID}; Layers{{1,..(nElem)..,1},{0.1,..(nCut)..,1}}; };
        For each point extruded, 1 new point and 1 new line are created
        '''        
        out = extdb({
        'newPoints': [],
        'newLines': []
        })
        ok_to_extrude=False
        for i in pointList:
            newpoint = Point(np.asarray(axis) + i.pos)
            if self._db_insert(newpoint): ok_to_extrude=True
            newline = Line(self,i,newpoint)
            if self._db_insert(newline): ok_to_extrude=True
            out['newPoints'].append(newpoint)
            out['newLines'].append(newline)
        if ok_to_extrude:
            idx_str = self._create_idx_str(pointList)
            axis_str = ','.join(Point._FLOAT_TO_STR.format(i) for i in axis)
            self._EXTRUDE_ID += 1
            self._CODE.append(
            'ex%d[] = Extrude {%s} { Point{%s}; Layers{%s}; };' %
            (self._EXTRUDE_ID, axis_str, idx_str, layers)
            )
        return out
        
    def _extrude_lines(self, lineList, axis, layers, opts=None):
        '''
        surface[] = Extrude{dx, dy, dz} { Line{#ID}; Layers{{1,..(nElem)..,1},{0.1,..(nCut)..,1}}; };
        For each line extruded, 2 new points, 3 new lines and 1 surface are created
        '''
        out = extdb({
        'newPoints': [],
        'newLines': [],
        'newSurfaces': []
        })
        axis_as_nparray = np.asarray(axis)
        ok_to_extrude=False
        for i in lineList:
            # 2 old point(s),
            oldPoint0 = self.get(Point, i.pid[0])
            oldPoint1 = self.get(Point, i.pid[1])
            # 2 new points
            newpoint0 = Point(axis_as_nparray + oldPoint0.pos)
            newpoint1 = Point(axis_as_nparray + oldPoint1.pos)
            # create 3 new lines
            if self._db_insert(newpoint0): ok_to_extrude=True
            if self._db_insert(newpoint1): ok_to_extrude=True            
            newline1 = Line(self,newpoint0,newpoint1)
            if self._db_insert(newline1): ok_to_extrude=True
            #
            self.incIDX(Point,2) # stupid gmsh
            newline2 = Line(self,oldPoint0,newpoint0)
            if self._db_insert(newline2): ok_to_extrude=True
            #
            self.incIDX(Point,2) # stupid gmsh
            newline3 = Line(self,oldPoint1,newpoint1)
            if self._db_insert(newline3): ok_to_extrude=True
            # create 1 new surface
            newsurf = Surface(self,[i,newline3,newline1,newline2])
            if self._db_insert(newsurf): ok_to_extrude=True
            out['newPoints'].append(newpoint0)
            out['newPoints'].append(newpoint1)
            out['newLines'].append(newline1)
            out['newLines'].append(newline2)
            out['newLines'].append(newline3)
            out['newSurfaces'].append(newsurf)
        if ok_to_extrude:
            idx_str = self._create_idx_str(lineList)
            axis_str = ','.join(Point._FLOAT_TO_STR.format(i) for i in axis)        
            opts_str = opts if opts is not None else 'Recombine;'
            self._EXTRUDE_ID += 1
            self._CODE.append(
            'ex%d[] = Extrude {%s} { Line{%s}; Layers{%s}; %s};' %
            (self._EXTRUDE_ID, axis_str, idx_str, layers, opts_str)
            )
        return out

    def _extrude_surfaces(self, surfList, axis, layers, opts=None):
        '''
        volume[] = Extrude{dx, dy, dz} { Surface{#ID}; Layers{{1,..(nElem)..,1},{0.1,..(nCut)..,1}}; };
        If the surface has n lines, we will create
            n new points,
            2*n new lines,
            n+1 new surfaces,
            and 1 volume
        '''
        out = extdb({
        'newPoints': [],
        'newLines': [],
        'newSurfaces': [],
        'newVolumes': [],
        })
        axis_as_nparray = np.asarray(axis)
        ok_to_extrude=False
        newp=out['newPoints']
        newl=out['newLines']
        news=out['newSurfaces']
        newv=out['newVolumes']
        for s in surfList:
            # extract ordered surface points
            sp=[]
            for i in s.lid:
                l=self.get(Line, i)
                if (i<0):
                    sp.append(self.get(Point,l.pid[1]))
                else:
                    sp.append(self.get(Point, l.pid[0]))
            n = len(sp) # the total number of point(s) on this surface
            # create line(s) parallel to old lines
            # treat 1st line (stupid gmsh), 2 newp, 1 newl
            newp.append(Point(axis_as_nparray + sp[0].pos))
            if self._db_insert(newp[-1]): ok_to_extrude=True
            newp.append(Point(axis_as_nparray + sp[1].pos))
            if self._db_insert(newp[-1]): ok_to_extrude=True
            newl.append(Line(self,newp[-2],newp[-1]))
            if self._db_insert(newl[-1]): ok_to_extrude=True
            # treat internal line(s), 1 newp, 1 newl for each internal line
            for i in sp[2:]:
                newp.append(Point(axis_as_nparray + i.pos))
                self.incIDX(Point,3)    # stupid gmsh
                if self._db_insert(newp[-1]): ok_to_extrude=True
                newl.append(Line(self,newp[-2],newp[-1]))
                if self._db_insert(newl[-1]): ok_to_extrude=True
            #
            # Important Note to myself:
            # Do not change self.incIDX(Point,???) before this line
            #
            # treat last line, no newp, 1 newl
            self.incIDX(Point,18)    # stupid gmsh
            newl.append(Line(self,newp[-1],newp[-n]))
            if self._db_insert(newl[-1]): ok_to_extrude=True
            # create lines in the extruded direction, n newl
            # the first two lines are treated differently (stupid gmsh)
            self.incIDX(Line,1) # stupid gmsh
            newl.append(Line(self, sp[0], newp[-n]))
            if self._db_insert(newl[-1]): ok_to_extrude=True
            newl.append(Line(self, sp[1], newp[-n+1]))
            if self._db_insert(newl[-1]): ok_to_extrude=True
            for i in range(2,n):
                self.incIDX(Point,6)    # stupid gmsh
                self.incIDX(Line,2) # stupid gmsh
                newl.append(Line(self, sp[i], newp[-n+i]))
                if self._db_insert(newl[-1]): ok_to_extrude=True
            #
            # Important Note to myself:
            # Do not change self.incIDX(Line,???) before this line
            #            
            # create n+1 new surfaces
            self.incIDX(Line,3) # stupid gmsh
            self.incIDX(Surface,1)  # stupid gmsh
            for i in range(0,n-1):
                news.append(Surface(self,[s.lid[i],newl[-2*n+i],newl[-n+i],newl[-n+i+1]]))
                if self._db_insert(news[-1]): ok_to_extrude=True                
            news.append(Surface(self,[s.lid[-1],newl[-n-1],newl[-n],newl[-1]]))
            if self._db_insert(news[-1]): ok_to_extrude=True
            lList=[] # last surface
            for i in range(0,n): lList.append(newl[-2*n+i])
            news.append(Surface(self,lList))
            if self._db_insert(news[-1]): ok_to_extrude=True                    
            # create 1 volume
            newv.append(Volume(self, [s,news[-1]] + news[-n-1:-1]))
            if self._db_insert(newv[-1]): ok_to_extrude=True                                   
        if ok_to_extrude:
            idx_str = self._create_idx_str(surfList)
            axis_str = ','.join(Point._FLOAT_TO_STR.format(i) for i in axis)        
            opts_str = opts if opts is not None else 'Recombine;'
            self._EXTRUDE_ID += 1
            self._CODE.append(
            'ex%d[] = Extrude {%s} { Surface{%s}; Layers{%s}; %s};' %
            (self._EXTRUDE_ID, axis_str, idx_str, layers, opts_str)
            )
        return out
        
        
                
        
