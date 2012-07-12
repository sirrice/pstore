import pdb
import random
import sys
sys.path.extend(['../', '../../'])

from models import *
from provstore import *
from op import *
from runtime import *
from arraystore import ArrayStore
from random import randint


import numpy as np
from scipy.optimize import curve_fit
from cStringIO import StringIO
import timeit




class BenchOp(Op):
    class Wrapper(object):
        def __init__(self, arr):
            shape = arr.shape
            self.nargs = 1
            self.shape = shape
            self._arr = arr

        def get_input_shapes(self, run_id):
            return [self.shape]

        def get_input_shape(self, run_id, arridx):
            return self.shape

        def get_inputs(self, run_id):
            return [self._arr]


    def __init__(self, arr, fanin, fanout):
        super(BenchOp, self).__init__()
        shape = arr.shape
        self.wrapper = BenchOp.Wrapper(arr)
        self.workflow = None
        self.shape = shape
        self.fanin = fanout
        self.fanout = fanout

    def run(self, inputs, run_id):
        pass

    def output_shape(self, run_id):
        return self.shape


    def gen_provenance(self, coord, prov_size):
        """
        @coord centroid of output
        @param prov_size either fanin or fanout value
        @return a list of coordinates that represent simulated
        provenance of input coordinate        
        """
        radius = max(int(math.ceil(math.sqrt(prov_size))), 1)
        x,y = tuple(coord)
        ix = min(self.shape[0], max(0, x - radius))
        iy = min(self.shape[1], max(0, y - radius))
        ret = []
        
        for i in xrange(prov_size):
            ix += 1
            if ix > min(self.shape[0], x+radius):
                ix = min(self.shape[0], max(0, x - radius))
                iy += 1
            if iy >= min(self.shape[1], y+radius):
                break

            ret.append((ix, iy))

        return ret


    def fmap_obj(self, obj, run_id, arridx):
        coords = obj[0]
        return self.gen_provenance(coords[0], self.fanin)

    
    def bmap_obj(self, obj, run_id, arridx):
        coords = obj[0]
        return self.gen_provenance(coords[0], self.fanin)







def setup_table(db):
    cur = db.cursor()
    create = """create table stats(id integer primary key autoincrement, strat varchar(128),
                fanin int, fanout int, noutput int, payload_size int,
                wcost float, incache float, outcache float, flush float, serout float, serin float,
                mergecost float, bdbcost float, disk int, idx int)"""
    try:
        cur.execute(create)
        db.commit()
    except Exception, e:
        print e
        db.rollback()

    create = """create table qcosts(sid int references stats(id), qsize int, backward bool,
                cost float, nres int,
                parsecost float, keycost float, datacost float, extractcost float, nhits int)"""
    try:
        cur.execute(create)
        db.commit()
    except  Exception, e:
        print e
        db.rollback()
        
    cur.close()



def get_prov(config, arr_shape):
    strat, fanin, fanout, noutput, payload_size = config    
    side = int(math.ceil(math.pow(float(fanin), 0.5)))
    oside = int(math.ceil(fanout ** 0.5))
    prov = []
    incoords = []
    outcoords = []
    fanin_n = 0
    maxx, maxy = arr_shape[0], arr_shape[1]

    # centroid of output coords
    optx, opty = randint(0, maxx - oside), randint(0, maxy - oside)
    n = 0
    nyielded = 0

    while True:

        outcoords = []
        incoords = []
        optx, opty = randint(0, maxx - oside), randint(0, maxy-oside) 


        while len(outcoords) < fanout and n < noutput:
            outcoord = (randint(0,oside) + optx,
                         randint(0,side) + opty ) 
            outcoords.append(outcoord)

            n += 1


        if len(outcoords) == 0:
            break

        # generate input coords either deterministically or
        # randomly

        if 'random' != 'random':
            ptx, pty = randint(0, maxx-side), randint(0, maxy-side)
            for i in xrange(fanin):
                incoord = (randint(0,side) + ptx,
                           randint(0,side) + pty) 
                incoords.append(incoord)
                fanin_n += 1
        else:
            radius = int(math.ceil(math.sqrt(float(fanin)) / 2.))
            ix, iy = optx - radius, opty - radius
            ix = max(0, min(maxx - side, ix))
            iy = max(0, min(maxy - side, iy))

            for i in xrange(fanin):
                incoords.append((ix, iy))
                fanin_n += 1

                ix += 1
                if ix >= opty + radius or ix >= maxx:
                    ix = max(0, min(maxx - side, ix))
                    iy += 1
                if iy >= maxy:
                    break

        yield (outcoords, incoords)

    return


        

def run_exp(db, configs, qsizes, arr_shape):
    

    arr = np.zeros(arr_shape)

    runid = 1
    setup_table(db)
    cur = db.cursor()    

    for config in configs:
        random.seed(0)
        print 'running %s\t%d\t%d\t%d\t%d' % tuple(config)
        strat, fanin, fanout, noutput, payload_size = config        
        
        op = BenchOp(arr, fanin, fanout)
        sid, pstore = run_config(cur, config, op, runid, arr_shape)
        run_queries(cur, qsizes, sid, pstore, arr_shape)
        db.commit()


    cur.close()


def run_config(cur, config, op, runid, arr_shape, niter=1):
    strat, fanin, fanout, noutput, payload_size = config
    Runtime.instance().set_strategy(op, strat)
    runid = 1

    # data loading costs
    costs = []
    for iterid in xrange(niter):

        pstore = op.pstore(runid)
        runid += 1

        
        # write all of the provenance
        start = time.time()
        for outcoords, incoords in get_prov(config, arr_shape):
            if pstore.uses_mode(Mode.PTR):
                pstore.write(outcoords, incoords)
            else:
                pstore.write(outcoords, 's' * payload_size)
        pstore.close()
        end = time.time()


        incache = pstore.get_stat('incache', 0)
        outcache = pstore.get_stat('outcache', 0)
        flush = pstore.get_stat('flush', 0)
        serin = pstore.get_stat('serin', 0)
        serout = pstore.get_stat('serout', 0)
        mergecost = pstore.get_stat('mergecost', 0)
        bdbcost = pstore.get_stat('bdbcost', 0)
        wcost = end - start #pstore.get_stat('write', 0)
        disk = pstore.disk()
        idx = pstore.indexsize()
        runcosts = (wcost, incache, outcache, flush, serout, serin, mergecost, bdbcost, disk, idx)
        costs.append(runcosts)

    stds = map(np.std, zip(*costs))
    costs = map(np.mean, zip(*costs))
    

    print 'provwrite cost\t%.4f\t%.4f' % (costs[0], stds[0])

    stratname = str(strat)
    if Mode.PT_MAPFUNC in strat.modes():
        stratname = '%s_%d' % (stratname, payload_size)
    
    params = [stratname, fanin, fanout, noutput, payload_size]
    params.extend(costs)
    sql = """insert into stats(strat, fanin, fanout, noutput, payload_size,
             wcost, incache, outcache, flush, serout, serin,
             mergecost, bdbcost, disk, idx ) values(%s)""" % (','.join(["?"]*(len(params)))) 
    cur.execute(sql, tuple(params))
    sid = cur.lastrowid

    return sid, pstore

def run_queries(cur, qsizes, sid, pstore, arr_shape):

    # query the provenance store
    # vary query size
    qcosts = []
    for qsize in qsizes:

        qcosts.append([])

        for backward in (True, False):

            stats = run_query(qsize, backward, pstore, arr_shape)

            params = [sid, qsize, backward]
            params.extend(stats)


            sql = """insert into qcosts values (%s)""" % (','.join(['?']*len(params)))
            cur.execute(sql, tuple(params))

def run_query(qsize, backward, pstore, arr_shape, niter=50):
    all_stats = []

    for i in xrange(niter):

        pstore.clear_stats()
        
        qcoords = [(random.randint(0, arr_shape[0]-1),
                    random.randint(0, arr_shape[1]-1))
                    for i in xrange(qsize)]
        q = Scan(qcoords)

        
        nres = 0
        start = time.time()

        for res in pstore.join(q, 0, backward):
            nres += 1

        cost = time.time() - start
	if cost > 15 and len(all_stats) > 5:
	    break
	    


        stats = [cost, nres]
        stats.extend([pstore.get_stat('parsecost', 0),
                      pstore.get_stat('keycost', 0),
                      pstore.get_stat('datacost', 0),
                      pstore.get_stat('extractcost', 0),
                      pstore.get_stat('nhits', 0) ])
        all_stats.append(stats)

    stds= map(np.std, zip(*all_stats))        
    all_stats = map(np.mean, zip(*all_stats))

    print backward, '\t', qsize, '\t', ('%.6f\t' * 4) % (all_stats[0], all_stats[1], stds[0], stds[1])
    
    return all_stats

if __name__ == '__main__':
    import sqlite3

    


    strats = [
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.GRID), True),            
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.GRID), True),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.BOX), True),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),



	# used strategies in the experiments
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),

        
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),            

        Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.BINARY), True),
        Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True),


        # forward optimized
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),        
        


    ]

    noutputs = (5000, 50000) #(1000, 10000, 100000)
    fanins = [1, 10, 50, 100]
    fanouts = [100, 50, 1]
    qsizes = [1, 100, 1000]
    payload_sizes = [50]#0, 10, 50, 100]
    arr_shape = (1000, 1000)


    def gen_configs(strats, noutputs, fanins, fanouts, payload_sizes):

        for noutput in noutputs:
            for strat in strats:

                # payload can ignore fanout values!
                if Mode.PT_MAPFUNC in strat.modes():
                    fanin = 10 

                    for fanout in fanouts:
                        for payload_size in payload_sizes:
                            yield (strat, fanin, fanout, noutput, payload_size)

                else:
                    for fanout in fanouts:
                        if fanout > noutput:
                            continue
                                                
                        for fanin in fanins:
                            yield (strat, fanin, fanout, noutput, 0)


    if len(sys.argv) <= 1:
        print "python pstore_micro.py [db path]"
        exit()
    dbname = sys.argv[1]
    db = sqlite3.connect(dbname)

    run_exp(db,
            gen_configs(strats, noutputs, fanins, fanouts, payload_sizes),
            qsizes,
            arr_shape)

    db.close()
