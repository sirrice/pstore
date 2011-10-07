import sys, random

from bench_util import *
sys.path.insert(0, "../protos/")
sys.path.insert(0, "../")
from models import *
from provstore import *
from op import *
from runtime import *

import numpy as np




def setup_box(strat, runtime, noutput, fanin, density, shape):
    boxop = BoxTestOp(runtime, noutput, fanin, density, shape) 
    Runtime.instance().set_strategy(boxop, strat)
    arr = boxop.wrapper.get_inputs(0)[0]

    start = time.time()
    fprefix, wpstore, rpstore = boxop.run([arr], 0)
    opcost = (time.time() - start) - wpstore.runtime
    writecosts = (wpstore.serializecost,
                  wpstore.writecost, 
                  wpstore.runtime,
                  opcost,
                  #runtime * reduce(mul, shape) / 100.0,
                  rpstore.disk_size())
    
    return arr, rpstore, writecosts

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

    res = PQDenseResult(arr.shape)
    cost = 0.0
    fcosts = []
    for i in xrange(5):
        res.clear()
        start = time.time()
        rpstore.forward_set(0, q, res, reduce(mul, shape))
        if i <= 1: continue
        cost = time.time() - start
        fcosts.append((cost, len(res)))
    fcosts = tuple(map(numpy.mean, zip(*fcosts)))        

    bcosts = []
    for i in xrange(5):
        res.clear()
        start = time.time()
        rpstore.backward_set(0, q, res, reduce(mul, shape))
        if i <= 1: continue
        cost = time.time() - start
        bcosts.append((cost, len(res)))
    bcosts = tuple(map(numpy.mean, zip(*bcosts)))    
        
        
    return [fcosts[0], fcosts[1], bcosts[0], bcosts[1]]

if __name__ == '__main__':
    qsizes = [1,10,50]
    densities = [1, 0.9, 0.5, 0.1, 0.05]
    runtimes = [0,1,5,10]
    fanin = 1000
    noutput = 100
    oclustsize = 1
    shape = (1000, 1000)
    ys = []
    labels = []
    strats = [STRAT_BOX, STRAT_PSET]

    emptyarr = np.empty(shape)
    arrid = ArrayStore.instance().add(emptyarr)
    inputsize = ArrayStore.instance().size(arrid)
    del emptyarr

    conn = get_db()
    runid = add_run(conn, inputsize, shape, notes= 'box_model')

    # qsizes = [10]
    # densities = [0.05]
    # runtimes = [1]
    for qsize in qsizes:
        for density in densities:
            labels.append("%f\t%d" % (density, qsize))
            ys.append([])        
            for strat in strats:
                for runtime in runtimes:
                    disk = disk_model(strat, fanin, 1, density, 100)
                    overhead = write_model(strat, fanin, 1, density, 100)
                                         
                    fcost = forward_model(strat, fanin, 1, density, 100, runtime, qsize, 1.0,
                                          reduce(mul, shape))
                    bcost = backward_model(strat, fanin, 1, density, 100, runtime, qsize, 1.0,
                                           reduce(mul, shape))
                    wcosts = (0,0,overhead,runtime,disk)
                    rcosts = (fcost, 0, bcost, 0)

                    params = (strat, fanin, oclustsize, density, noutput, qsize)
                    #arr, rpstore, wcosts = setup_box(strat, runtime, noutput, fanin, density, shape)
                    #rcosts = run_box(runtime, arr, rpstore, qsize, 0.0)

                    add_datapoint(conn, runid, params, wcosts, rcosts)
                    
                    print "%f\t%f\t%d\t%s\t%f\t%f" % (density, runtime, qsize, strat, rcosts[0], rcosts[2])
                    #rpstore.close()
                    #del arr
                    #del rpstore



    conn.close()

