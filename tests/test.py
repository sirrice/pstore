import sys
sys.path.append("../")

from runtime import *
from op import *
from provstore import *
from strat import *
from lsst import *


        


if __name__ == '__main__':
    import pyfits
    
    op = Convolve()
    op = Cluster()
    #op = CRDetect()
    w = Workflow()
    w.register(op, 2)

    hdulist = pyfits.open('../data/imsim_85408556_R23_S11_C04_E001.fits')    
    scidata = hdulist[0].data
    size = 50
    scidata = scidata[:size,:size]

    arr = numpy.empty((10, 10))
    kernel = np.array([[0,-1,0], [-1,4,-1], [0,-1,0]])

    strats = [
        # Strat.single(Mode.PT_MAPFUNC|Mode.FULL_MAPFUNC,
        #              Spec(Spec.COORD_ONE, Spec.BINARY), False),
        # Strat.single(Mode.PT_MAPFUNC|Mode.FULL_MAPFUNC,
        #              Spec(Spec.COORD_MANY, Spec.BINARY), False),

        # Strat.single(Mode.QUERY, Spec(Spec.NONE, Spec.NONE), True),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), False),                
        # Strat.single(Mode.PT_MAPFUNC|Mode.FULL_MAPFUNC,
        #              Spec(Spec.COORD_ONE, Spec.BINARY), True),
        #Strat.single(Mode.PT_MAPFUNC|Mode.FULL_MAPFUNC,
        #             Spec(Spec.COORD_MANY, Spec.BINARY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),                
        ]

    for strat in strats:
        w.add_static_input(op, 0, scidata)
        w.add_static_input(op, 1, numpy.array([[1650]]))
        #w.add_static_input(op, 1, kernel)
        
        Runtime.instance().set_strategy(op, strat)
        w.run()

        pstore = op.pstore(w._runid-1)

        coords = np.argwhere(op.wrapper.get_output(1))[:20]
        q = (coords, w._runid-1, [(op, 0)])
        n = 0
        start = time.time()
        for coord in w.forward_path(*q):
            n += 1
        print n, (time.time()-start), pstore.stats.get('write',0), pstore.disk(), strat

