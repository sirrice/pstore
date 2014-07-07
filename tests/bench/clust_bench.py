import sys, random
from bench_util import *
sys.path.insert(0, "../protos/")
sys.path.insert(0, "../")
sys.path.insert(0, "../../")
from models import *
from provstore import *
from op import *
from lsst import *
from common import *
from strat import *
from runtime import *


import numpy as np


def runbench(op, strat, arrs):
    # setup Runtime to do the right thing
    runtime = Runtime.instance()
    runtime.reset()


#    w = Workflow()
#    w.register(op, 2)

    runtime.set_strategy(op, strat)
    #thresh = np.array([[100,100]])
    wpstore = op.pstore(0)
    wop = op.wrapper
    op.run(arrs, 0)
    #wop.run(0)
    return wpstore




import pyfits, numpy
hdulist1 = pyfits.open('../../data/imsim_85408556_R23_S11_C04_E000.fits')
scidata1 = hdulist1[0].data
smallshape = (1000,1000)

for opclass in [CRDetect]:
    for size in [10, 50, 100, 500, 1000]:
        arr = scidata1[100:100+size,100:100+size]
        kernel = np.array([[0,-1,0], [-1,4,-1], [0,-1,0]])
        #arr = arr - numpy.mean(arr)
        for strat in [Mode.FULL_MAPFUNC, Mode.NOOP]:
            op = bench_op(opclass)(strat, [arr, kernel])
            start = time.time()
            wpstore =  runbench(op, Mode.NOOP, (arr, kernel))
            cost = time.time() - start
            print "%s\t%s\t%d\t%f\t%d\t%f\t%s" % (opclass.__name__, strat, size, cost,
                                                  wpstore.nptrs, wpstore.outsize(),
                                                  str(wpstore.fanins()))

    
