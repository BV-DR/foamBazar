import numpy as np
import math
from scipy.optimize import fsolve


def getIDX(geom, obj):
    return getattr(geom,obj._ID_NAME)
    
def getDB(geom, obj):
    return getattr(geom,obj._DB_NAME)

def unique_and_keep_order(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]    

def remove_bracket(input_str):
    return input_str.replace('[','').replace(']','').replace(' ','')
    
def exist(geom, obj):
    db=getDB(geom,obj)
    masterDBKeys=obj.masterDBKeys(geom)
    masterKey=obj.key(master=True)
    # we need the option to match only a substring of keys()
    if masterKey in masterDBKeys: return True,db[db.keys()[masterDBKeys.index(masterKey)]]
    return False, getIDX(geom,obj) + 1

def stretch(method, n=None, r=None, opts=None):
    '''
    Return stretched coordinate i.e. a list of numbers between 0 ... 1
    The supported stretched method(s) are as follows:
    
        'ratio'     :   input n=(int, default:11), r=(float, default:2)
                        Given the ratio between the first and the last element
                        compute the position(s) for n element(s)
        
        'rate'      :   input n=(int, default:11), r=(float, default:1.1)
                        r is the growth/shrink rate
                        E.g.: for r=1.1 the next cell is 1.1 larger

    Both of the above method support the following options

        opts={'firstCell': (float)}
                        The first cell is fixed. Adjust other cells based on
                        (approximated value of) either "n" or "r". Either "n" or
                        "r" must be set to "None". A RuntimeError will be raised
                        if both "n" and "r" are defined.

        opts={'lastCell': (float)}
                        Same as firstCell, but instead the last cell is fixed.
                        
        Note: Avoid both 'firstCell' and 'lastCell' can be used simultaneously
    '''
    assert isinstance(method,str)
    if n is not None: assert isinstance(n, (int,long))
    if r is not None: assert isinstance(r, (int,long,float))
    if opts is not None: assert isinstance(opts, dict)
    if n<=0: n=None
    if r<=0: r=None
    if method.lower() not in {'ratio', 'rate'}:
        raise RuntimeError('stretch: unkown stretching method')
    if method.lower() in {'ratio'}: return _stretch_ratio(n,r,opts)
    if method.lower() in {'rate'}: return _stretch_rate(n,r,opts)   
    return

def _stretch_ratio(n, r, opts):
    print "DEBUG: Ratio-based stretching "
    return
    
def _stretch_rate(n, r, opts):
    #print "DEBUG: Rate-based stretching "
    def stretch_func(h, n, r, inverse=False):
        if inverse:
            out=[1.0/h]
            for i in range(1,n): out.append(out[-1]+pow(float(r),i)/h)
            return out
        else:
            out=[h]
            for i in range(1,n): out.append(out[-1]+h*pow(float(r),i))
            return out
        return
    def check_limit(h):
        if (h<=1e-10) or (h>1):
            raise RuntimeError('_stretch_rate: firstCell must be a positive float less than 1')
        else:
            return h
    if opts is None:
        # n,r are fixed. Compute the first cell
        if n is None: n=11
        if r is None: r=1.1
        oneOver_h=1.0;
        for i in range(1,n): oneOver_h += pow(r,i)
        return stretch_func(oneOver_h,n,r,inverse=True)
    
    #from this point on either 'n' or 'r' must be defined
    if ((n is not None) and (r is not None)) or ((n is None) and (r is None)):
        raise RuntimeError('_stretch_rate: when fixing the (first|last)Cell pls. define either "n" or "r" (not both)')

    firstCell=None
    lastCell=None
    
    # firstCell is fix.
    if ('firstCell' in opts) and ('lastCell' not in opts):
        h=check_limit(float(opts['firstCell']))
        #
        # if h==1, obviously the only solution is to have n=1, r=1
        if math.fabs(h-1.0) < 1e-10:
            return stretch_func(1.0, 1, 1.0)
        #
        # define(firstCell, n), free(r)
        if r is None:
            #print "DEBUG: n="+str(n), "r=free", "firstCell: ", h
            # if n==1, the only solution is to have h=1 too, so we need to
            # recompute "n" to preserve h
            if n<=1:
                return _stretch_rate(2, None, opts)
            # check if r=1 is good enough
            if math.fabs(1.0/n - h) < 1e-10:
                return stretch_func(h, n, 1.0)
            else:
                # we need to solve a nonlinear equation, h*r^n - r - h + 1 = 0
                # first we estimate the starting point
                def f(x):
                   return float(h) * x ** float(n) - x + 1.0 - h
                def x0():
                   x=math.exp(math.log(1.0/n/h)/(n - 1))
                   if x > 1.0:
                      return x + 0.5*(x-1.0)
                   else:
                      return 0.5 * x
                r=fsolve(f, x0())[0]
                #print "DEBUG: r: " + str(r)
                return stretch_func(h, n, r)
            pass
        #
        # define(firstCell, r), free(n)
        # Note: we "may" need to change r to preserve "h", that's why we need to
        # call _stretch_rate(...) with free r !!!
        if n is None:
            #print "DEBUG: n=free, r="+str(r), "firstCell: ", h
            #
            # check if r==1.0
            if math.fabs(r-1.0) < 1e-10:
                n=max(2,int(math.floor(1.0/h))) # the minimum n is 2
                return _stretch_rate(n,None,opts)
            #
            # if r <= 1-h, we don't have a solution, so we are free to decide n,r
            # here we take n closest to 1/h but smaller than 10
            if math.fabs(r+h-1.0) < 1e-10:
                n=min(10,max(2,int(math.floor(1.0/h)))) # min. 2, max. 10
                return _stretch_rate(n,None,opts)
            #
            # we try the close-form expression for n
            # FIXME: why r=0.95, h=0.02 failed?
            
            n=max(2,int(round(math.log((r+h-1.0)/h)/math.log(r))))
            return _stretch_rate(n,None,opts)
        pass
        
    # lastCell is fix.
    if ('lastCell' in opts) and ('firstCell' not in opts):
        # we compute it as it was the firstCell and then reverse it !!!
        if r is None:
            out=_stretch_rate(n,None,opts={'firstCell': float(opts['lastCell'])})
        else:
            out=_stretch_rate(n,2.0-r,opts={'firstCell': float(opts['lastCell'])})
        return np.asarray((1.0-np.array([0.0] + out))[::-1][1:])
    
    # firstCell and lastCell are defined simultaneously,
    # for now, we don't support this
    
    return 
   
    
    
