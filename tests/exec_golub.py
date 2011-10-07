from golub import *
import numpy as np
import random
from scipy import ndimage
from util import write_fits, log_prov
from guppy import hpy
from nlp import *
from stats import *

def set_all_strategy(strat):
    def set_strategy(wop):
        if wop.op.implements_mapfunctions():
            Runtime.instance().set_strategy(wop.op, STRAT_F)
        else:
            Runtime.instance().set_strategy(wop.op, strat)
    w.visit(set_strategy)


if __name__ == "__main__":
    import logging, os, sys


    if len(sys.argv) < 2:
        print "exec_lsst_full.py [opt|box|list|set|grid|test] [log output dir] "
        exit()

    runtype, logdir = sys.argv[1], sys.argv[2]
    runtypes = ['opt', 'box', 'list', 'set', 'grid', 'test']
    if runtype not in runtypes:
        exit()
    if os.path.exists(logdir) and not os.path.isdir(logdir):
        print "%s exists and is not a directory"  % logdir
        exit()

    if not os.path.exists(logdir):
        os.mkdir(logdir)


    logging.basicConfig(filename="%s/out_%s" % (logdir, runtype), filemode = 'w')
    mlog = logging.getLogger('main')
    mlog.setLevel(logging.INFO)


    random.seed(0)
    # mem profile object
    hp = hpy()


    def print_stats(w):
        stats = Stats.instance()
        for tup in stats.runtime_stats:
            mlog.info( '%s\t'*8 % tuple(tup))
        mlog.info('\n')
        
        for tup in stats.pq_stats:
            mlog.info( '%s\t' * 5 % tuple(tup))

        mlog.info( "avg prov q %f\t%f",
                   numpy.mean([tup[-1] for tup in stats.pq_stats]),
                   numpy.std([tup[-1] for tup in stats.pq_stats]))
    



    clean1 = Clean()
    clean2 = Clean()
    filt = Filter()
    normtrain = NormalizeTrain()
    normtest = NormalizeTest()
    corr = Correlate()
    ttest = TTest()
    stats = GenStats()
    mask_corr = PredictorMask()
    mask_ttest = PredictorMask()
    predict_corr = Predict()
    predict_ttest = Predict()

    w = Workflow()
    w.register(clean1, 1)
    w.register(clean2, 1)
    w.register(filt, 1)
    w.register(normtrain, 1)
    w.register(normtest, 2)
    w.register(corr, 2)
    w.register(ttest, 2)
    w.register(stats, 2)
    w.register(mask_corr, 2)
    w.register(mask_ttest, 2)
    w.register(predict_corr, 3)
    w.register(predict_ttest, 3)
    
    w.connect(clean1, filt, 0)
    w.connect(clean1, normtrain, 0)
    #w.connect(normtrain, corr, 0)
    w.connect(clean1, corr, 0)
    w.connect(normtrain, ttest, 0)
    w.connect(normtrain, stats, 0)
    w.connect(clean2, normtest, 0)
    w.connect(stats, normtest, 1)
    w.connect(filt, mask_corr, 0)
    w.connect(corr, mask_corr, 1)
    w.connect(filt, mask_ttest, 0)
    w.connect(ttest, mask_ttest, 1)
    w.connect(normtest, predict_corr, 0)
    w.connect(corr, predict_corr, 1)
    w.connect(mask_corr, predict_corr, 2)
    w.connect(normtest, predict_ttest, 0)
    w.connect(ttest, predict_ttest, 1)
    w.connect(mask_ttest, predict_ttest, 2)


    w.default_strategy()


    strategies = {'list' : STRAT_PLIST, 'set' : STRAT_PSET,
                  'test' : STRAT_PTEST, 'grid' : STRAT_PGRID,
                  'box' : STRAT_BOX}
    if runtype in strategies:
        strat = strategies[runtype]
        set_all_strategy(strat)

        if runtype =='test':
            Stats.instance().save('%s/stats.dat' % logdir)
    else:
        set_all_strategy(STRAT_F)
        if os.path.exists('%s/stats.dat' % logdir):
            stats = Stats.instance()
            stats.load('%s/stats.dat' % logdir)
            strategies = run_nlp(stats, w, 5, 3000)

            for op, strats in strategies.items():
                for s in strats:
                    print op, s
                    Runtime.instance().set_strategy(op, s)

            stats.clear()
            exit()
        elif os.path.exists('nlpresult.dat'):
            assignmentmatrix = parse_nlpresult()
            strategies = assign_strategies(assignmentmatrix,
                                           w.get_optimizable_ops(),
                                           w.get_matstrats())

            for op, strats in strategies.items():
                for s in strats:
                    Runtime.instance().set_strategy(op, s)
        else:
            raise RuntimeError, "I don't know how to optimize!"
        

        

    
    names, training, test = load_gene_data('data_aml_all.txt')
    npatients = numpy.array([27])



    handler = logging.FileHandler('%s/provres_%s.log' % (logdir, runtype), 'w')
    handler.setLevel(logging.DEBUG)
    mlog.addHandler(handler)

    for niter in xrange(1):
        w.add_static_input(clean1, 0, training)
        w.add_static_input(clean2, 0, test)
        w.add_static_input(corr, 1, npatients)
        w.add_static_input(ttest, 1, npatients)
        w.add_static_input(stats, 1, npatients)
        w.run()
        corr_res = predict_corr.wrapper.get_output(w._runid-1)
        for rowid in xrange(1):#len(corr_res)):
            log_prov(mlog, w.backward([(rowid,0)], predict_corr, w._runid-1, 5 ))


    mlog.removeHandler(handler)
    handler.close()

    handler = logging.FileHandler('%s/perf_%s.log' % (logdir, runtype), 'w')
    handler.setLevel(logging.DEBUG)
    mlog.addHandler(handler)

    print_stats(w)


    mlog.removeHandler(handler)
    handler.close()


    #for crow in corr_res:
    #    print '%d\t%f' % (crow[0], crow[1])

    # mask = mask_corr.wrapper.get_output(w._runid-1)

    # tmp = []
    # for name, b, weight, row in zip(names, mask, corr.wrapper.get_output(w._runid-1),
    #                                 normtrain.wrapper.get_output(w._runid-1)):
    #     if b:
    #         tmp.append(( name, weight, row))
    # tmp.sort(key=lambda v: v[1])
    # for x in tmp: print '%s\t%f' % (x[0].ljust(25), x[1])



