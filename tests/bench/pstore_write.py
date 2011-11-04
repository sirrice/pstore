import sys, random
sys.path.insert(0, "../")
sys.path.insert(0, "../../")

from bench_util import *
from models import *
from provstore import *
from op import *
from runtime import *
from arraystore import ArrayStore

import numpy as np


def run_model(arr, strat, noutput, fanin, oclustsize, density):
    runtime = 0.0038750171661377
    disk = disk_model(strat, fanin, oclustsize, density, noutput)
    cost = write_model(strat, fanin, oclustsize, density, noutput, runtime)
    fcost = forward_model(strat, fanin, oclustsize, density, noutput, runtime,
                          10, 1.0, inputarea=reduce(mul, arr.shape))
    bcost = backward_model(strat, fanin, oclustsize, density, noutput, runtime,
                           10, 1.0, inputarea=reduce(mul, arr.shape))
    return (0,0,cost,runtime,disk), (fcost,0,bcost,0)

    
def run_exp(arr, strat, noutput, fanin, oclustsize, density):
    rpstore, writecosts = run_op(arr, strat, noutput, fanin, oclustsize, density)
    fcosts = run_qs(rpstore, arr, strat, noutput, fanin, oclustsize, density, backward=False)
    bcosts = run_qs(rpstore, arr, strat, noutput, fanin, oclustsize, density, backward=True)
    rpstore.close()
    del rpstore
    return (writecosts, [fcosts[0], fcosts[1], bcosts[0], bcosts[1]])

def run_op(arr, strat, noutput, fanin, oclustsize, density):
    writecosts = []
    pstore = None
    for i in xrange(3):
        if pstore: pstore.close()
        Runtime.instance().clear()    
        
        op = BenchOp(arr, strat, noutput, fanin, oclustsize, density)
        Runtime.instance().set_strategy(op, strat)
        start = time.time()
        op.wrapper._arr = arr
        pstore, gencost = op.run([arr], i)
        opcost = (time.time() - start) - pstore.get_stat('write', 0)
        

        if i == 0: continue
        updatecost = pstore.get_stat('update_stats', 0)
        bdbcost = pstore.get_stat('bdb', 0)
        serin = pstore.get_stat('serin', 0)
        serout = pstore.get_stat('serout', 0)
        ser = pstore.get_stat('_serialize', 0)
        wcost = pstore.get_stat('write', 0)
        runcosts = (ser, wcost, gencost, updatecost, bdbcost, serin, serout)
        # runcosts = (ser,
        #             pstore.get_stat('write', 0) / pstore.nsampled / (fanin + oclustsize), 
        #             pstore.get_stat('write', -1) + pstore.get_stat('close', 0),
        #             opcost,
        #             pstore.disk())
        writecosts.append(runcosts)
    writecosts = map(numpy.mean, zip(*writecosts))
    return pstore, writecosts

def run_qs(rpstore, arr, strat, noutput, fanin, oclustsize, density, backward=True):
    costs = []
    try:
        for i in xrange(3):
            for coords in gen_forward_qs(arr, 50, 10, 1.0, noutput, fanin, oclustsize, density):
                rpstore.clear_stats()
                start = time.time()
                nres = 0
                for coord in rpstore.join(coords, 0, backward=backward):
                    nres += 1
                if i > 0:
                    #costs.append((time.time()-start, rpstore.get_stat('_parse')))
                    costs.append((time.time()-start, nres))
        costs = tuple(map(numpy.mean, zip(*costs)))
        return costs
    except:
        return (-1, -1)

    

def printit(params, wcosts, rcosts):
    rcosts = []
    pstr = '\t'.join(map(str, params))
    wstr = '\t'.join(map(lambda x: '%f'%x, tuple(wcosts)))
    rstr = '\t'.join(map(lambda x: '%f'%x, tuple(rcosts)))
    print '\t'.join([pstr, wstr, rstr])




if __name__ == "__main__":
    import sys

    conn = get_db()


    random.seed(0)
    shape = (1000,1000)
    emptyarr = np.empty(shape)
    arrid = ArrayStore.instance().add(emptyarr)
    inputsize = ArrayStore.instance().size(arrid)

    if len(sys.argv) < 2:
        print """python pstore_write.py  [run type]  [run_id]
        
args:
    runtype: model|exp
             default=exp
             
    runid:   >0 to specify runid
             -1 to not store
             nothing to create new runid"""
        exit()

    if sys.argv[1].strip() == 'model':
        experiment = run_model
        notes = 'ptr_model'
    else:
        experiment = run_exp
        notes = 'ptr'
        
    if len(sys.argv) > 2:
        runid = int(sys.argv[2])
    else:
        runid = add_run(conn, inputsize, shape, notes=notes)

    print "runid: ", runid
    print "experiment", sys.argv[1]

    strat, noutput, fanin, oclustsize, density = None, 100, 5,5, 1.0
    fanins = [1,25,64,100,169]
    oclustsizes = [1,25,64,100]
    densities = [0.01, 0.1, 0.5, 0.9]
    noutputs = [1000]


    strats = [
        Strat.single(Mode.FULL_MAPFUNC, Spec(Spec.NONE, Spec.NONE), True),
        Strat.single(Mode.QUERY, Spec(Spec.NONE, Spec.NONE), True),
        
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
        
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),        
        
        ]

    #noutputs = [10, 100, 1000, 5000]
    #fanins = [1, 10]
    #oclustsizes = [1, 50, 100]
    #strats = [Strat.stat()]
    fanins = [1, 10]
    densities = []
    for noutput in noutputs:
        for strat in strats:
            for fanin in fanins:
                for oclustsize in oclustsizes:
                    for density in [1.0]:
                        arr = gen_interarr(emptyarr, noutput, fanin, oclustsize, density)
                        wcosts, rcosts = experiment(arr, strat, noutput, fanin, oclustsize, density)
                        params = (strat, fanin, oclustsize, density, noutput, 10)
                        if runid > 0:
                            add_datapoint(conn, runid, params, wcosts, rcosts)
                        printit(params, wcosts, rcosts)
                        del arr

            fanin, oclustsize = 100, 100
            for density in densities:
                arr = gen_interarr(emptyarr, noutput, fanin, oclustsize, density)
                wcosts, rcosts = experiment(arr, strat, noutput, fanin, oclustsize, density)
                params = (strat, fanin, oclustsize, density, noutput, 10)
                if runid > 0:                
                    add_datapoint(conn, runid, params, wcosts, rcosts)
                printit(params, wcosts, rcosts)
                del arr



    conn.close()
