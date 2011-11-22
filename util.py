import numpy, os, random, math
from operator import mul
from scipy import ndimage
import bsddb3

def get_db(fname, new=False, fsync=False):
    env = bsddb3.db.DBEnv()
    if fsync:
        env.set_flags(bsddb3.db.DB_TXN_NOSYNC, 0)
    else:
        env.set_flags(bsddb3.db.DB_TXN_NOSYNC, 1)
    env.set_flags(bsddb3.db.DB_NOLOCKING, 1)
    env.set_flags(bsddb3.db.DB_AUTO_COMMIT, 0)
    env.open('.', bsddb3.db.DB_PRIVATE | bsddb3.db.DB_CREATE | bsddb3.db.DB_INIT_MPOOL)
    db = bsddb3.db.DB(env)
    if new and fname is not None and os.path.isfile(fname):
        os.unlink(fname)
    db.open(fname,
            dbtype=bsddb3.db.DB_HASH,
            flags=bsddb3.db.DB_CREATE,
            mode=0666)
    return db



def subarray(arr, coords):
    """
    given coordinates of interest, returns the smallest subarray that contains
    all points in coords

    coords: a list containing points that lie on the bounding box
    return subarray, f:subarray coord -> arr coord
    """
    if len(coords) == 0:
        return numpy.array([])
    ndim = len(coords[0])
    mins = numpy.zeros(ndim)
    maxs = numpy.zeros(ndim)
    for dim in xrange(ndim):
        mins[dim] = max(min([x[dim] for x in coords]), 0)
        maxs[dim] = min(max([x[dim] for x in coords]), arr.shape[dim])
    # add 1 to compensate for xrange
    shape = map(int, maxs-mins+1)

    def gen_indices(shape, minv, dim):
        outertimes = reduce(mul, shape[:dim], 1)
        innertimes = reduce(mul, shape[dim+1:], 1)
        return [i + minv for i in xrange(shape[dim]) for j in xrange(innertimes)] * outertimes

    nparr = numpy.array(arr)
    indices = [gen_indices(shape, minv, dimidx)
               for (dimidx, dimsize), minv in zip(enumerate(shape), mins)]
    newarr = nparr[indices].reshape(shape)
    newarr = newarr#.tolist()
    
    def addmins(subcoord, mins):
        return tuple(map(int,map(lambda x: x[0]+x[1], zip(subcoord, mins))))

    def submins(subcoord, mins):
        return tuple(map(int,map(lambda x: x[0]-x[1], zip(subcoord, mins))))

    #submins = lambda p: tuple(map(int,map(lambda x: x[0]-x[1], zip(p[0], p[1]))))
            
    return (newarr, 
            lambda subcoord: submins(subcoord, mins),
            lambda subcoord: addmins(subcoord, mins))
    


def ncells(arr):
    # TODO: actual checksum
    if len(arr) == 0:
        return 0
    if type(arr[0]) == list:
        return sum(map(lambda item: ncells(item), arr))
    return len(arr)


def print_matrix(varname, m):
    s = ';'.join([', '.join([str(c) for c in row]) for row in m])
    s = '%s = [%s]' % (varname, s)
    #print s
    return s

def write_fits(arr, fname):
    import pyfits
    nparr = numpy.array(arr)
    hdu = pyfits.PrimaryHDU(nparr)
    outname = os.path.basename(os.path.abspath('./%s' % fname))
    if os.path.exists(outname): os.remove(outname)
    hdu.writeto(outname)    

def zipf(n, l = 1.5):
    hns = [1.0 / math.pow(i + 1, l) for i in xrange(n)]
    c = sum(hns)
    probs = [c / math.pow(i+1, l) for i in xrange(n)]
    sumprobs = sum(probs)
    probs = [prob / sumprobs for prob in probs]
    ret = [sum(probs[:i+1]) for i in xrange(n)]
    return ret


def log_prov(log, prov=[]):
    # summarize using a single value (total number of coordinates)
    log.info("Num prov coords\t%d" , len(prov))
    return

    
    for path, vals in prov.items():
        total = sum([len(coords) for vals in prov.values() for arrid, coords in vals ])
        log.info("Num prov coords\t%d", total)
    return
    
    
    
    for path, vals in prov.items():
        if len(path) == 0: continue
        if len(path) == 1:
            log.info(path[0])
        else:
            log.info("%s...%s", path[0], path[-1])
        for idx, val in enumerate(vals):
            if val:
                (arrid, coords) = val
                log.info("     arg(%d)\t%s\t%d",idx, arrid, len(coords))
            else:
                log.info("     arg(%d)\t%s\t0", idx, arrid)
        log.info("\t")


def enc(coord, dims):
    """
    linearize the coordinate based on the dimension sizes
    """
    return sum(map(lambda x: x[0] * x[1], zip(coord, dims)))

def dec(v, dims):
    """
    extract coordinates from integer value
    """
    insize = dims
    coord = []
    for i in xrange(len(insize)-1, -1, -1):
        if i > 0:
            mod = v % insize[i-1]
        else:
            mod = v
        coord.append(mod / insize[i])
        v -= mod
    coord.reverse()
    return coord


        
     
def gen_data(xsize, ysize, nstars=3, starradius=10, brightness=2000):
    # 1) lots of stars, big
    # 2) lots of tiny stars
    # 3) few stars, big
    # 4) few stars, tiny

    footprint = ndimage.generate_binary_structure(2,1)

    ret = numpy.zeros((xsize, ysize))
    for star in xrange(nstars):
        xcenter = random.randint(0, xsize-1)
        ycenter = random.randint(0, ysize-1)
        for x in xrange(xcenter-1, xcenter+2):
            for y in xrange(ycenter-1, ycenter+2):
                if x >= 0 and y >= 0 and x < xsize and y < ysize:
                    ret[x,y] = brightness / 3
        ret[xcenter, ycenter] = brightness
    for i in xrange(starradius):
        ret = ndimage.grey_dilation(ret, footprint=footprint)

    # add some cosmic rays (single points)
    for i in xrange(30):
        xcenter = random.randint(0, xsize-1)
        ycenter = random.randint(0, ysize-1)
        ret[xcenter, ycenter] = brightness

    return ret








if __name__ == '__main__':

    arr = numpy.arange(25).reshape((5,5))
    subarr, convert = subarray(arr, [(1,1), (4,4)])
    print subarr[1,2], arr[convert((1,2))]
    
    exit()
    
    
    dims = [7, 7, 7]


    mults = [reduce(lambda x,y: x*y, dims[i+1:], 1)
             for i in xrange(len(dims))]
    mults[-1] = 1
    dims = mults

    print dims

    for x in xrange(7):
        for y in xrange(7):
            for z in xrange(7):            
                v = enc((x,y,z), dims)
                xx, yy, zz = tuple(dec(v, dims))
                if xx != x or yy != y or zz != z:
                    print "UHOH", x, y, z, v, (xx, yy, zz)
    


