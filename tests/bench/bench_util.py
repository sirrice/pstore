import sys, os, random
sys.path.append("../")
sys.path.append("../../")
import struct
import numpy as np
from op import *
from runtime import *
import psycopg2 as pg

CACHEDIR = "./_cache"


def cluster(arr):
    threshold = 0.5
    ids = np.zeros(arr.shape)
    nrow, ncol = arr.shape        
    fixups = {}
    nextid = 1
    clusters = {}

    for x in xrange(nrow):
        for y in xrange(ncol):
            if arr[x,y] > threshold:
                if y > 0 and ids[x,y-1]:
                    oid = ids[x,y-1]
                elif y > 0 and x > 0 and ids[x-1,y-1]:
                    oid = ids[x-1,y-1]
                elif x > 0 and ids[x-1,y]:
                    oid = ids[x-1,y]
                elif x > 0 and y < len(arr[0])-1 and ids[x-1,y+1]:
                    oid = ids[x-1,y+1]
                else:
                    oid = nextid
                    nextid += 1

                ids[x,y] = oid

                if x > 0 and y < len(arr[0])-1:# and oldid and oid != oldid:
                    oldid = ids[x-1,y+1]                        
                    if oldid and oid > oldid:
                        fixups[oid] = oldid
                        self.fix(fixups, oid)

                if oid not in clusters:
                    clusters[oid] = []
                clusters[oid].append((x,y))
    return clusters

class BenchOp(Op):
    class Wrapper(object):
        def __init__(self, shape):
            self.nargs = 1
            self.shape = shape
            self._arr = None

        def get_input_shapes(self, run_id):
            return [self.shape]
        
        def get_input_shape(self, run_id, arridx):
            return self.shape

        def get_inputs(self, run_id):
            return [self._arr]

    
    def __init__(self, arr, strat, noutput, fanin, oclustsize, density):
        super(BenchOp, self).__init__()
        self.wrapper = BenchOp.Wrapper(arr.shape)
        self.workflow = None
        self.noutput = noutput
        self.fanin = fanin
        self.oclustsize = oclustsize
        self.density = density
        self.strat = strat

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id, 0)

    
    def run(self, inputs, run_id):

        arr = inputs[0]
        
        fprefix = pstore_name(arr.shape, self.strat, self.noutput, self.fanin,
                              self.oclustsize, self.density)
        pstore = self.pstore(run_id)

        tmparr = np.zeros(arr.shape)
        ingen = inter_to_inlist(arr, self.noutput, self.fanin, self.oclustsize, self.density)
        outgen = inter_to_outlist(arr, self.noutput, self.fanin, self.oclustsize, self.density)

        if pstore.strat.mode == Mode.PTR:
            for incoords, outcoords in zip(ingen, outgen):
                pstore.write(outcoords, incoords)
        if pstore.strat.mode == Mode.PT_MAPFUNC:
            for incoords, outcoords, in zip(ingen, outgen):
                pstore.write(outcoords, '')

        pstore.close()
        return fprefix, pstore

    def fmap_obj(self, obj, run_id, arridx):
        return obj[0]

    def bmap_obj(self, obj, run_id, arridx):
        return obj[0]

    def fmap(self, coord, run_id, arridx):
        return [(0,0)]

    def bmap(self, coord, run_id, arridx):
        return [(0,0)]

    def supported_model(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]
        
    

class BoxTestOp(Op):
    class Wrapper(object):
        def __init__(self, shape, noutput):
            self.nargs = 1
            self.shape = shape
            self.noutput = noutput
            self._arr = None

        def get_input_shape(self, run_id, arridx):
            return self.shape

        def get_inputs(self, run_id):
            if self._arr == None:
                self._arr = create_array(self.shape, self.noutput)
            return [ self._arr ] 

    def __init__(self, runtime, noutput, fanin, density, shape):
        super(BoxTestOp, self).__init__()
        self.runtime = runtime
        self.fanin = fanin
        self.density = density
        self.noutput = noutput
        self.strat = Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True)
        self.shape = shape
        self.wrapper = BoxTestOp.Wrapper(shape, noutput)
        self.workflow = None


    def run(self, inputs, run_id):
        random.seed(0)
        arr = inputs[0]
        arrsize = float(reduce(mul, arr.shape))
        totsize = float(reduce(mul, self.shape))
        tosleep = float(self.runtime) * arrsize / totsize

        time.sleep(tosleep)

        pstore = self.pstore(run_id)        
        start = time.time()
        cost = 0.0
        for outcoords, incoords in box_prov(arr, self.runtime, self.fanin, self.density):
            pstore.write(outcoords, incoords)
        rpstore = pstore.close()
        cost += time.time() - start

        #print "sleeping %f\t%d\t%d\tcost %f" % ( tosleep, arrsize, reduce(mul,self.shape), cost)
        return pstore.fprefix, pstore
    
        

def box_prov(arr, runtime, fanin, density):
    shape = arr.shape
    min_insize = [int(math.ceil(math.pow(fanin, 0.5)))] * 2
    exp_insize = [int(min_insize[0] / math.pow(density, 0.5))] * 2
    insize = map(min, zip(map(max, zip(min_insize, exp_insize)), shape))
    for outcoord in np.argwhere(arr):
        offsets = map(lambda p: max(0,p[0]-p[1]) ,
                      zip(map(sum, zip(outcoord, insize)), shape))
        mins = map(lambda p: p[1] - p[0], zip(offsets, outcoord))
        maxs = map(lambda p: p[1] - p[0], zip(offsets, map(sum, zip(outcoord, insize))))
        incoords = sample(mins, maxs, fanin)
        yield [tuple(outcoord)], map(tuple,incoords)


def sample(mins, maxs, n, rand=random):
    #random.seed(0)
    dims = map(lambda p: p[1] - p[0], zip(mins, maxs))
    maxv = dims[0] * dims[1]
    if n > maxv:
        return [(i,j) for i in range(mins[0], maxs[0]) for j in range(mins[1], maxs[1])]
    if n == 0:
        return []
    if n == 1:
        return [mins]
    if n == 2:
        return [mins, maxs]

    if True or (maxv >= n ** 2 / 10 and maxv >= n * 10):
        encs = set()
        while len(encs) < n-2:
            encs.add(random.randrange(1, maxv-1))
        encs = list(encs)
    else:
        encs = range(1, maxv-1)
        rand.shuffle(encs)
        encs = encs[:n-2]
    coords = []
    encs.extend([0,maxv-1])
    coords = [(int(math.floor(enc / dims[1])) + mins[0], enc % dims[1] + mins[1]) for enc in encs]
    return coords


def create_array(shape, noutput):
    random.seed(0)
    outsize = [1,1]
    arr = np.zeros(shape, int)
    noutputgenerated = 0
    for i in xrange(noutput):
        outcoord = [random.randint(0,(dim-minus)) for dim, minus in zip(shape, outsize)]
        arr[outcoord[0],outcoord[1]] = 1
    return arr
    



def arrcache(postfix=""):
    def dec(fn):
        def w(arr, noutput, fanin, oclustsize, iareasize):
            fname = "%s/%dx%d_%d_%d_%d_%f.%s.npy" % (CACHEDIR, arr.shape[0], arr.shape[1], noutput, fanin, oclustsize, iareasize, postfix)

            if not os.path.exists(fname):
                ret = fn(arr, noutput, fanin, oclustsize, iareasize)
                ret.tofile(fname)
            else:
                ret = numpy.fromfile(fname).reshape(arr.shape)
            return ret
        return w
    return dec


def gencache(postfix=""):
    def dec(fn):
        def w(arr, noutput, fanin, oclustsize, iareasize):

            fname = "%s/%dx%d_%d_%d_%d_%f.%s.cache" % (CACHEDIR, arr.shape[0], arr.shape[1], noutput, fanin, oclustsize, iareasize, postfix)
            
            if not os.path.exists(fname):
                f = file(fname, 'wb')
                allcoords = []
                for coords in fn(arr, noutput, fanin, oclustsize, iareasize):
                    allcoords.append(coords)

                f.write(struct.pack("I", len(allcoords)))
                for coords in allcoords:
                    buf = []
                    map(buf.extend, coords)
                    f.write(struct.pack("%dI" % (2*len(coords)+1), len(coords), *buf))
                f.close()
            else:
                f = file(fname, 'rb')
                totaln, = struct.unpack("I", f.read(4))
                for x in xrange(totaln):
                    n, = struct.unpack("I", f.read(4))
                    coords = []
                    for i in xrange(n):
                        coords.append(struct.unpack("II", f.read(4*2)))
                    yield coords
                f.close()
        return w
    return dec
            
            
            
@arrcache("interarr")
def gen_interarr(blankarr, noutput, fanin, oclustsize, iareasize):
    """
    Generates an array with ground truth of the centroids for each output cluster
    output array and input array are deterministically derived from this.
    """
    random.seed(0)

    shape = blankarr.shape
    arr = np.zeros(shape)

    # calculate oclust size
    shape = map(lambda x: x-1, shape)
    outedgesize = int(math.ceil(math.pow(oclustsize, 0.5)))
    xblocks, yblocks = (shape[0] - outedgesize) / float(outedgesize), \
                       (shape[1] - outedgesize) / float(outedgesize)

    if xblocks * yblocks * outedgesize * outedgesize < noutput:
        print xblocks , yblocks , outedgesize ,noutput        
        raise RuntimeError

    outsize = [int(math.ceil(math.pow(oclustsize, 0.5)))] * 2        
    curxblock, curyblock = 0,0

    noutputgenerated = 0
    while noutputgenerated < noutput:
        coord = (curxblock * outedgesize, curyblock*outedgesize)
        arr[coord] = 1

        curxblock += 1
        if curxblock >= xblocks: curxblock, curyblock = 0, curyblock + 1
        noutputgenerated += oclustsize
    
    return arr

            
@gencache("inlist")
def inter_to_inlist(arr, noutput, fanin, oclustsize, iareasize):

    r = random.Random()
    r.seed(0)
    intercoords = map(tuple, np.argwhere(arr))

    min_insize = [int(math.ceil(math.pow(fanin, 0.5)))] * 2
    exp_insize = [int(min_insize[0] / iareasize), int(min_insize[1] / iareasize)]
    insize = map(min, zip(map(max, zip(min_insize, exp_insize)), arr.shape))
    insize = map(max, zip(insize, [1]*len(insize)))

    for intercoord in intercoords:
        intercoord = tuple(intercoord)
        ur = map(sum, zip(insize, intercoord))
        overflow = [max(0, v-bound) for v, bound in zip(ur, arr.shape)]
        ur = map(min, zip(ur, arr.shape))
        ll = [v - w for v,w in zip(intercoord, overflow)]
        incoords = sample(ll, ur, fanin, rand=r)
        yield incoords

@arrcache("inarr")
def inter_to_inarr(arr, noutput, fanin, oclustsize, iareasize):
    ret = np.zeros(arr.shape)
    for incoords in inter_to_inlist(arr, noutput, fanin, oclustsize, iareasize):
        for incoord in incoords:
            ret[tuple(incoord)] = 1
    return ret

@gencache("outlist")
def inter_to_outlist(arr, noutput, fanin, oclustsize, iareasize):
    r = random.Random()
    r.seed(0)
    intercoords = map(tuple, np.argwhere(arr))
    outsize = [int(math.ceil(math.pow(oclustsize, 0.5)))] * 2

    for intercoord in intercoords:
        coords = sample(intercoord, map(sum, zip(outsize, intercoord)), oclustsize, rand=r)
        yield coords


@arrcache("outarr")
def inter_to_outarr(arr, noutput, fanin, oclustsize, iareasize):
    ret = np.zeros(arr.shape)
    for coord in inter_to_outlist(arr, noutput, fanin, oclustsize, iareasize):
        ret[tuple(coord)] = 1
    return ret
    
        



def gen_forward_qs(arr, nqs, qsize, perc_matches, noutput, fanin, oclustsize, iareasize):
    incoordsall = set()
    for coords in inter_to_inlist(arr, noutput, fanin, oclustsize, iareasize):
        incoordsall.update(coords)

    nfound = int(math.ceil(perc_matches * qsize))
    notfounds = []
    j = 0
    for x in xrange(1000):
        if len(notfounds) >= qsize - nfound: break
        for y in xrange(1000):
            if len(notfounds) >= qsize - nfound: break
            if (x,y) not in incoordsall:
                notfounds.append((x,y))

                    
    
    i = 0
    incoordsall = list(incoordsall)[:nqs]
    random.shuffle(incoordsall)
    qs = []
    for incoord in incoordsall:
        qs.append(incoord)
        if len(qs) >= nfound:
            qs.extend(notfounds)
            yield qs
            qs = []
        i += 1
        if i >= nqs: break
    if len(qs):
        yield qs



def gen_backward_qs(arr, nqs, qsize, perc_matches, noutput, fanin, oclustsize, iareasize):

    outall = set()
    for coords in inter_to_outlist(arr, noutput, fanin, oclustsize, iareasize):
        outall.update(coords)


    nfound = int(math.ceil(perc_matches * qsize))
    notfounds = []
    j = 0
    for x in xrange(1000):
        if len(notfounds) >= qsize - nfound: break
        for y in xrange(1000):
            if len(notfounds) >= qsize - nfound: break
            if (x,y) not in outall:
                notfounds.append((x,y))

    i = 0
    outall = list(outall)[:nqs]
    random.shuffle(outall)
    qs = []
    for outcoord in outall:
        qs.append(outcoord)
        if len(qs) >= nfound:
            qs.extend(notfounds)
            yield qs
            qs = []
        i += 1
        if i >= nqs: break
    if len(qs):
        yield qs



# def viz(xs, ys, strats = [STRAT_BOX, STRAT_F, STRAT_PSET, STRAT_PGRID, STRAT_BULK]):
#     import matplotlib.pyplot as plt
#     #plt.ion()
#     colors = ['b', 'r', 'o', 'y', 'g', 'p']
#     symbols = []
#     for color in colors:
#         symbols.extend([color+'-', color+'--', color+'-.'])
    
#     plt.cla()        
#     plt.clf()
#     for idx in xrange(len(ys)):
#         plt.plot(xs, ys[idx], symbols[idx], label=str(strats[idx]))
#     plt.legend()
#     plt.show()
#     plt.draw()


def pstore_name(shape, strat, noutput, fanin, oclustsize, density):
    fprefix = './pstore/%dx%d_%d_%d_%d_%f_%d_%s' % (shape[0], shape[1],
                                                    noutput, fanin, oclustsize,
                                                    density, 0,
                                                    str(strat).strip())
    return fprefix



def get_db():
    conn = pg.connect("dbname=microbench")
    setup(conn)
    return conn

def setup(conn):
    c = conn.cursor()
    try:
        c.execute("""create table runs (id serial unique primary key,
                     tstamp timestamp default now(),
                     arrsize int,
                     shape text,
                     real bool,
                     notes text)""")
        c.execute("""create table stats (id serial,
                     rid int references runs(id),
                     strat text,
                     fanin int,
                     fanout int,
                     density float,
                     noutput int,
                     nqs int,
                     sercost float,
                     writecost float,
                     overhead float,
                     opcost float,
                     disk int,
                     fcost float,
                     fnres int,
                     bcost float,
                     bnres int
                     )""")
        conn.commit()
    except Exception, e:
        print e
        conn.rollback()
    finally:
        c.close()
    


def add_run(conn, arrsize, shape, real=True,notes=''):
    c = conn.cursor()

    c.execute("insert into runs values (default, default, %s, %s, %s, %s) returning id",
              (arrsize, 'x'.join(map(str, shape)), real,notes))
    rowid = c.fetchone()[0]
    conn.commit()
    c.close()
    return rowid

def add_datapoint(conn, runid, params, wcosts, rcosts):
    # add to sqlite database
    q = "insert into stats values(default, %s)" % (','.join(["%s"] * 16))

    args = [runid]
    args.append(str(params[0]).strip())
    args.extend(params[1:])
    args.extend(wcosts)
    args.extend(rcosts)
    c = conn.cursor()
    c.execute(q, tuple(args))
    conn.commit()
    c.close()
    
    

    
    
def str_to_strat(s):
    arr = [x.strip() for x in s[2:-2].split(',')]
    mode = int(arr[0])
    outcoords, payload = tuple([int(x.strip()) for x in arr[1].split('/')])
    spec = Spec(outcoords, payload)
    backward = bool(arr[-1])
    return Strat.single(mode, spec, backward)
    


if __name__ == '__main__':
    import sys
    sys.path.append("../")
    sys.path.append("../../")
    from runtime import *

    results = {}
    strats = set()

    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:
        fname = './benchmark_model'
    
    f = file(fname, "r")
    for l in f:
        arr = l.strip().split('\t')
        strat = str_to_strat(arr[0])

        params = tuple(map(float,arr[1:4]))
        if params[1] == 9: continue
        #(serializecost, writecost, runtime, disk_size)
        wcosts = tuple(map(float, arr[4:]))

        row = results.get(params, {})
        strats.add(strat)
        row[strat] = wcosts
        results[params] = row
    f.close()

    strats = sorted(strats)

    
    def f(wlabels, extract):
        if len(wlabels) != 1: raise RuntimeError
        
        print "====%s====" % (wlabels[0])
        
        wheader = ['fanin', 'oclustsize', 'density']
        wheader.extend([str(strat) for strat in strats])
        #wheader.extend(['%s_%s' % (strat, l) for strat in strats for l in wlabels])
        wheader = '\t'.join(wheader)
        print wheader

        for params, row in results.items():
            pstr =  '\t'.join(map(str, params))
            dat = []
            for strat in strats:
                costs = row.get(strat, None)
                if costs == None:
                    dat.append(0)
                else:
                    dat.append(extract(costs))


            dstr = '\t'.join(map(lambda x: '%f'%x, dat))
            print '\t'.join([pstr, dstr])

    


    f(['runtime'], lambda costs: costs[5])
    f(['disk'], lambda costs: costs[7])
    f(['fq'], lambda costs: costs[-4])
    f(['bq'], lambda costs: costs[-2])



    
