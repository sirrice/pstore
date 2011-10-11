import sys, random

from bench_util import *
sys.path.insert(0, "../protos/")
sys.path.insert(0, "../")
from models import *
from provstore import *
from op import *
from runtime import *

import numpy as np


def run_exp(strat, runtime, noutput, fanin, density, shape, qsize):
    arr, pstore, wcosts = setup_box(strat, runtime, noutput, fanin, density, shape)
    qcosts = run_box(runtime, arr, pstore, qsize)
    return wcosts, qcosts

def run_model(strat, runtime, noutput, fanin, density, shape, qsize):
    disk = disk_model(strat, fanin, 1, density, 100)
    overhead = write_model(strat, fanin, 1, density, 100)
    fcost = forward_model(strat, fanin, 1, density, 100, runtime, qsize, 1.0,
                          reduce(mul, shape))
    bcost = backward_model(strat, fanin, 1, density, 100, runtime, qsize, 1.0,
                           reduce(mul, shape))
    wcosts = (0,0,overhead,runtime,disk)
    rcosts = (fcost, 0, bcost, 0)
    return wcosts, rcosts
        

def setup_box(strat, runtime, noutput, fanin, density, shape):
    Runtime.instance().clear()    
    boxop = BoxTestOp(runtime, noutput, fanin, density, shape)
    Runtime.instance().set_strategy(boxop, strat)
    arr = boxop.wrapper.get_inputs(0)[0]

    start = time.time()
    pstore = boxop.run([arr], 0)
    opcost = (time.time() - start) - pstore.get_stat('write', 0)
    writecosts = (pstore.get_stat('_serialize'),
                  0,
                  pstore.get_stat('write'),
                  opcost,
                  pstore.disk())
    
    return arr, pstore, writecosts

def run_box(runtime, arr, rpstore, qsize, perc_missing = 0.0):
    rpstore.op.runtime = runtime
    
    shape = arr.shape
    coords = set(map(tuple,np.argwhere(arr)))
    # notfounds = []
    # nfound = int((1.0-perc_missing) * qsize)
    # 
    # for x in xrange(shape[0]):
    #     for y in xrange(shape[1]):
    #         if (x,y) not in coords:
    #             notfounds.append((x,y))
    #         if len(notfounds) > qsize-nfound: break
    #     if len(notfounds) > qsize-nfound: break
    #q.extend(notfounds)

    q = list(coords)[:qsize]

    cost = 0.0
    fcosts = []
    for i in xrange(5):
        start = time.time()
        nres = 0
        for coord in rpstore.join(q, 0, False):
            nres += 1
        if i <= 1: continue
        cost = time.time() - start
        fcosts.append((cost, nres))
    fcosts = tuple(map(numpy.mean, zip(*fcosts)))        

    bcosts = []
    for i in xrange(5):
        start = time.time()
        nres = 0
        for coord in rpstore.join(q, 0, True):
            nres += 1
        if i <= 1: continue
        cost = time.time() - start
        bcosts.append((cost, nres))

    bcosts = tuple(map(numpy.mean, zip(*bcosts)))    
        
        
    return [fcosts[0], fcosts[1], bcosts[0], bcosts[1]]

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 2:
        print """python pstore_box.py  [run type]  [run_id]
        
args:
    runtype: model|exp
             default=exp
             
    runid:   >0 to specify runid
             -1 to not store
             nothing to create new runid"""
        exit()

    
    qsizes = [1,10,50]
    densities = [1, 0.9, 0.5, 0.1, 0.05]
    runtimes = [0,1,5,10]
    fanin = 1000
    noutput = 100
    oclustsize = 1
    shape = (1000, 1000)
    ys = []
    labels = []
    strats = [
        Strat.single(Mode.QUERY, Spec(Spec.NONE, Spec.NONE), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),        
        ]

    emptyarr = np.empty(shape)
    arrid = ArrayStore.instance().add(emptyarr)
    inputsize = ArrayStore.instance().size(arrid)
    del emptyarr

    conn = get_db()

    if sys.argv[1].strip() == 'model':
        experiment = run_model
        notes = 'box_model'
    else:
        experiment = run_exp
        notes = 'box'

    if len(sys.argv) > 2:
        runid = int(sys.argv[2])
    else:
        runid = add_run(conn, inputsize, shape, notes=notes)

    print "runid: ", runid
    print "experiment: ", sys.argv[1]


    # qsizes = [50]
    # densities = [0.05]
    # runtimes = [10]
    # strats = [Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True)]
    #           #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True)    ]
    for qsize in qsizes:
        for density in densities:
            labels.append("%f\t%d" % (density, qsize))
            ys.append([])        
            for strat in strats:
                for runtime in runtimes:
                    params = (strat, fanin, oclustsize, density, noutput, qsize)

                    wcosts, rcosts = experiment(strat, runtime, noutput, fanin, density, shape, qsize)

                    if runid > 0:
                        add_datapoint(conn, runid, params, wcosts, rcosts)
                    
                    print "%f\t%f\t%d\t%s\t%f\t%d\t%f\t%d" % (density, runtime, qsize, strat,
                                                              rcosts[0], rcosts[1], rcosts[2], rcosts[3])
                    #rpstore.close()
                    #del arr
                    #del rpstore



    conn.close()

