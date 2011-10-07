import sys
sys.path.append("../")
from lsst import *
from util import log_prov
from runtime import *

if __name__ == '__main__':

    mlog = logging.getLogger('main')
    logging.basicConfig()
    mlog.setLevel(logging.INFO)


    runtype = sys.argv[1]
    runtypes = ['box', 'bulk', 'set', 'grid', 'test', 'query', 'sql', 'run',
                'noop', 'super', 'sbox', 'squery']
    if runtype not in runtypes:
        print "%s is not a recognized runtype" % runtype
        exit()

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
        queries.append([[(0,3)], runid, [(mean,0), (diff,0), (conv,0), (clus,0)]])
        queries.append([[(3,3)], runid, [(mean,0), (diff,0), (conv,0), (clus,0)]])
        for q in queries:
            print len( w.forward_path(*q))

        output = clus.wrapper.get_output(runid)
        indices = np.argwhere(output)
        print len(w.backward_path([ indices[0] ], runid, [(clus,0)]))


        
    
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
    #strat = Strategy.from_str(runtype)
    # Runtime.instance().set_strategy(clus, strat)

    img = numpy.array([[0,0,1,1],
                       [0,0,0,1],
                       [1,0,0,0],
                       [1,1,0,0]], dtype=np.float)
    run_img(img, numpy.array([[0.5, 0.5]]))
