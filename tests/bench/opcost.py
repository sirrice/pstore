"""
compute query costs for simple LSST-like workflows
"""

import sys, random

from bench_util import *
sys.path.insert(0, "../protos/")
sys.path.insert(0, "../")
from models import *
from provstore import *
from op import *
from lsst import *
from runtime import *
from arraystore import ArrayStore
from util import *
import numpy as np
import pyfits


if __name__ == '__main__':
    random.seed(0)
    shape = 2000
    img = gen_data(shape, shape, nstars=30, starradius=4, brightness=2000)
    hdulist = pyfits.open('../data/imsim_85408556_R23_S11_C04_E001.fits')    
    img = hdulist[0].data
    img = img[0:100,0:100]
    kernel = np.array([[0,-1,0], [-1,4,-1], [0,-1,0]])

    
    stats = Stats.instance(fname='./opcost_stats.db')
    stats.add_exec(shape, shape, -1, '', '', 'opcost benchmark')
    
    
    w = Workflow()

    # mean = Mean()
    # diff = Mean()
    # w.register(mean, 1)
    # w.register(diff, 1)
    # w.connect(mean, diff, 0)


    mean = Mean()
    rmbg = Diff()
    cr = CRDetect()
    rmcr = RemoveCRs()
    #op = Cluster()
    w.register(mean, 1)
    w.register(rmbg, 2)
    w.register(cr, 2)
    w.register(rmcr, 2)
    w.connect(mean, rmbg,1)
    w.connect(rmbg, cr, 0)
    w.connect(cr, rmcr, 1)
    
    kernel = np.array([[0,-1,0], [-1,4,-1], [0,-1,0]])
    #w.add_static_input(op, 0, img)
    #w.add_static_input(op, 1, numpy.array([[1]]))
    #w.add_static_input(op, 1, kernel)
    #w.add_static_input(op, 1, numpy.array([[10000]]))
    #Runtime.instance().set_strategy(op, STRAT_SQ)



    def run_setup(stratid, strats):
        Stats.instance().add_exec(img.shape[0], img.shape[1],
                                  -1, '_'.join(map(lambda p:str(p[1]).strip(), strats)),
                                  './',
                                  "benchmark-%d" % stratid)
    
        
        for op, strat in strats:
            Runtime.instance().set_strategy(op, strat)

        w.add_static_input(mean, 0, img)
        w.add_static_input(rmbg, 0, img)
        #w.add_static_input(op, 0, img)
        w.add_static_input(cr, 1, kernel)
        w.add_static_input(rmcr, 0, img)
        w.run()
        runid = w._runid-1

        write_fits(cr.wrapper.get_output(runid).astype(int), './generated1.fits')
        write_fits(rmcr.wrapper.get_output(runid), './generated2.fits')


        random.seed(0)
        qs = []
        bpath = fpath = [(rmcr,0), (cr,0)]
        bpath = fpath = [(cr,0)]

        for nq in [1, 10, 100, 500]:
            qcoords = [(random.randint(0, img.shape[0]-1), random.randint(0, img.shape[1]-1))
                       for i in xrange(nq)]
            qs.append([qcoords, runid, bpath, 'backward'])
            qs.append([qcoords, runid, fpath, 'forward'])
            

        costs = []
        for maxres in [1, -1]:
            for inputs, runid, path, direction in qs:
                nres = 0        
                start = time.time()
                if direction == 'forward':
                    res = w.forward_path(inputs, runid, path)
                else:
                    res = w.backward_path(inputs, runid, path)
                for coord in res:
                    nres += 1
                    if maxres != -1 and nres >= maxres: break

                qcost = time.time()-start
                costs.append(qcost)
                print "nres", nres, '\t',qcost

                path_ops = [x[0] for x in path]
                pqid = Stats.instance().add_pq(runid, path_ops, direction,
                                               [len(inputs)], qcost, maxres)
                for op, arridx in path:
                    Stats.instance().add_iq(pqid, op, arridx,
                                            Runtime.instance().get_strategy(op, runid),
                                            direction, len(inputs), 0, 0)

        Stats.instance().finish_exec()                
        return costs
        

    allstrats = []
    #allstrats.append([(mean, STRAT_F), (diff, STRAT_F) ])
    allstrats.append([(mean, STRAT_F), (rmbg, STRAT_F), (cr, STRAT_DIFF), (rmcr, STRAT_DIFF)])
    allstrats.append([(mean, STRAT_F), (rmbg, STRAT_F), (cr, STRAT_FDIFF), (rmcr, STRAT_FDIFF)])


    for sidx, s in enumerate(allstrats):
        costs = run_setup(sidx, s)
        print "strat-%d\t%s" % (sidx, '\t'.join(map(str,costs)))

    Stats.instance().close()
