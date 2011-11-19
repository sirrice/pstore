import sys
sys.path.append("../")
from lsst import *
from util import log_prov
from runtime import *

if __name__ == '__main__':

    mlog = logging.getLogger('main')
    logging.basicConfig()
    mlog.setLevel(logging.INFO)

    Stats.instance('_output/execsmall.db')

    runmode = -1
            

    def run_img(img, threshold):
        kernel = np.array([[1,1,1], [1,1,1], [1,1,1]])
        
        
        w.add_static_input(mean, 0, img)
        w.add_static_input(diff, 0, img)
        w.add_static_input(conv, 1, kernel)
        w.add_static_input(clus, 1, threshold)
        w.run()
        
        #Runtime.instance().set_exec_mode(Runtime.NORMAL)
        print
        runid = w._runid-1


        queries = []
        queries.append([[(0,3)], runid, [(clus,0)]])
        queries.append([[(0,3)], runid, [(mean,0), (diff, 0), (conv, 0)]])
        queries.append([[(0,3)], runid, [(mean,0), (diff,0), (conv,0), (clus,0)]])
        queries.append([[(3,3)], runid, [(mean,0), (diff,0), (conv,0), (clus,0)]])
        for q in queries:
            print len( w.forward_path(*q)[0])
        output = clus.wrapper.get_output(runid)
        indices = map(tuple,np.argwhere(output))
        print len(w.backward_path([ indices[0] ], runid, [(clus,0)])[0])


        
    
    mean = Mean()
    diff = Diff()
    conv = Convolve()
    clus = Cluster()

    w = Workflow()
    w.register(mean, 1)
    w.register(diff, 2)  # need to pass in original
    w.register(conv, 2)  # need to pass in kernel
    w.register(clus, 2)

    w.connect(mean, diff, 1)  
    w.connect(diff, conv, 0)  
    w.connect(conv, clus, 0)

    w.default_strategy()
    strat = Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False)
    print strat
    Runtime.instance().set_strategy(clus, strat)

    if len(sys.argv) <= 1:
        print "python exec_lsst_small.py [small | big]"
        exit()
    elif sys.argv[1] == 'small':

        # small image!
        img = numpy.array([[0,0,1,1],
                           [0,0,0,1],
                           [1,0,0,0],
                           [1,1,0,0]], dtype=np.float)
    elif sys.argv[1] == 'big':
        # big image!
        img = numpy.zeros((500,500), dtype=float)
        for i in xrange(50):
            for j in xrange(50):
                img[i,j] = 1
    else:
        print "python exec_lsst_small.py [small | big]"
        exit()


    Stats.instance().add_exec(10, 10,
                              runmode, 'grid', '_output', "lsst_grid")
    
    run_img(img, numpy.array([[0.5, 0.5]]))

    Stats.instance().finish_exec()    
