if __name__ == "__main__":
    import logging, os, sys
    sys.path.append("../")
    import pyfits, numpy
    from runtime import *

    if len(sys.argv) < 3:
        print "exec_lsst_full.py [opt|box|list|set|grid|test] [log output dir] [run mode: 0, 1] [shape size]"
        exit()
    print sys.argv
    runtype, logdir, runmode = sys.argv[1], sys.argv[2], int(sys.argv[3])

    if len(sys.argv) > 3:
        smallshape = (int(sys.argv[4]), int(sys.argv[4]))
    else:
        smallshape = (2000,2000)

    if os.path.exists(logdir) and not os.path.isdir(logdir):
        print "%s exists and is not a directory"  % logdir
        exit()

    if not os.path.exists(logdir):
        os.mkdir(logdir)


    logging.basicConfig(filename="%s/out_%s" % (logdir, runtype), filemode = 'w')
    mlog = logging.getLogger('main')
    mlog.setLevel(logging.INFO)



    from lsst import *
    import numpy as np
    import random
    from scipy import ndimage
    from util import write_fits, log_prov
    #from guppy import hpy
    from nlp import *
    from modelpredict import ModelPredictor


    random.seed(0)
    # mem profile object
    #hp = hpy()


    hdulist1 = pyfits.open('../data/imsim_85408556_R23_S11_C04_E000.fits')
    hdulist2 = pyfits.open('../data/imsim_85408556_R23_S11_C04_E001.fits')    
    scidata1 = hdulist1[0].data
    scidata2 = hdulist1[0].data
    print hdulist1.info()

    # go through all subsets of the large fits image
    if runmode == 1:
        bigshape = scidata1.shape
        smallshape = [min(x,y) for x,y in zip(bigshape, smallshape)]
        factor = tuple(numpy.floor(numpy.array(bigshape) / numpy.array(smallshape)).astype(int).tolist())
    
        print factor
    else:
        bigshape = smallshape
        factor = [1,1]
    iters, maxiters = 0, 1
    darkmask1 = biasmask1 = satmask1 = flatmask1 = None
    darkmask2 = biasmask2 = satmask2 = flatmask2 = None


    
    Stats.instance('%s/lsstfull.db' % logdir)
    

    


    def run_img(exp1, exp2, threshold):
        kernel = np.array([[0,-1,0], [-1,4,-1], [0,-1,0]])

        
        #Runtime.instance().set_exec_mode(Runtime.WARMUP)
        w.add_static_input(linearity1, 0, exp1)
        w.add_static_input(badpixrm1, 1, satmask1)
        w.add_static_input(biascorrect1, 1, biasmask1)
        w.add_static_input(darkcorrect1, 1, darkmask1)
        w.add_static_input(flatcorrect1, 1, flatmask1)
        #w.add_static_input(cr1, 0, img)
        w.add_static_input(cr1, 1, kernel)
        #w.add_static_input(rmcr1, 0, img)
        w.add_static_input(linearity2, 0, exp2)
        w.add_static_input(badpixrm2, 1, satmask2)
        w.add_static_input(biascorrect2, 1, biasmask2)
        w.add_static_input(darkcorrect2, 1, darkmask2)
        w.add_static_input(flatcorrect2, 1, flatmask2)
        #w.add_static_input(cr1, 0, img)
        w.add_static_input(cr2, 1, kernel)

        w.add_static_input(clust, 1, threshold)


        start = time.time()
        w.run()



    def run_prov_workload(maxoutput=-1):
        print "========= running provenance queries ========"
        queries = gen_prov_workload(w._runid - 1)
        totalcost = 0.0
        for inputs, runid, path, direction in queries:
            start = time.time()

            if direction == 'forward':
                res = w.forward_path(inputs, runid, path)
            else:
                res = w.backward_path(inputs, runid, path)

            n = 0
            for coord in res:
                n+=1
                if maxoutput != -1 and n >= maxoutput: break
            qcost = time.time() - start
            totalcost += qcost

            path_ops = [x[0] for x in path]
            pqid = Stats.instance().add_pq(w._runid-1, path_ops, direction,
                                           [len(inputs)], qcost, maxoutput)
            print "results", n, qcost
        return totalcost
            

    def gen_prov_workload(runid):
        # initially all 4 queries.
        # random bursts of 1a, 1b queries
        # regular 3 and 4 queries
        # annually, s2 queries
        random.seed(0)
        queries = []
        queries.extend(scenario2(runid))
        queries.extend(scenario1a(runid))
        queries.extend(scenario1b(runid))
        #queries.extend(scenario3(runid))
        queries.extend(scenario4(runid))
        return queries

        if runid < 10:
            queries.extend(scenario4(runid))
        else:
            if random.random() > 0.01:
                queries.extend(scenario1a(runid))
            if random.random() > 0.01:
                queries.extend(scenario1b(runid))
            if runid % 4:
                if random.random() > 0.5:
                    queries.extend(scenario3(runid))
                if random.random() > 0.5:
                    queries.extend(scenario4(runid))


    def scenario1a(runid):

        # # scenario 1a: too many masked pixels from Cosmic Ray
        # # backward from masked pixels to beginning
        # # high probability for past couple days
        # # overscan to sources
        mlog.info("Scenario 1a")
        mask = w.wrapper(rmcr).get_output(runid)
        coords = np.argwhere(mask)

        path1 = [(rmcr,1), (maskexp,0), (cr1, 0), (rmbg1, 0), (flatcorrect1, 0),
                (darkcorrect1, 0), (biascorrect1, 0), (rmoverscan1, 0), (badpixrm1, 0),
                (addmaskattr1, 0), (linearity1, 0)]
        path2 = [(rmcr,1), (maskexp,1), (cr2, 0), (rmbg2, 0), (flatcorrect2, 0),
                 (darkcorrect2, 0), (biascorrect2, 0), (rmoverscan2, 0), (badpixrm2, 0),
                 (addmaskattr2, 0), (linearity2, 0)]

        queries = []
        random.shuffle(coords)
        coords = coords[:10]
        queries.append([coords, runid, path1, 'backward'])
        queries.append([coords, runid, path2, 'backward'])
        return queries

        
        if len(coords):
            for i in xrange(1):
                coord = tuple(coords[random.randint(0, len(coords)-1)].tolist())
                queries.append([ [coord], runid, path1, 'backward' ])
                queries.append([ [coord], runid, path2, 'backward' ])
        return queries

    def scenario1b(runid):
        # # scenario 1b: too many stars (suspect false positives)
        # # backward from Sources to beginning
        # # higher probability on recent days
        mlog.info("Scenario 1b")
        output = clust.wrapper.get_output(runid)
        indices = np.argwhere(output).tolist()
        random.shuffle(indices)


        path1 = [(clust,0), (rmcr,1), (maskexp,0), (cr1, 0), (rmbg1, 0), (flatcorrect1, 0),
                (darkcorrect1, 0), (biascorrect1, 0), (rmoverscan1, 0), (badpixrm1, 0),
                (addmaskattr1, 0), (linearity1, 0)]
        path2 = [(clust,0), (rmcr,1), (maskexp,1), (cr2, 0), (rmbg2, 0), (flatcorrect2, 0),
                 (darkcorrect2, 0), (biascorrect2, 0), (rmoverscan2, 0), (badpixrm2, 0),
                 (addmaskattr2, 0), (linearity2, 0)]

        queries = []
        # random distribution of num stars
        #for star in indices[:1]:
        queries.append([ indices[:10], runid, path1, 'backward' ] )
        queries.append([ indices[:10], runid, path2, 'backward' ] )
        return queries
            
    def scenario2(runid):
        # scenario 2: portion of exposure is faulty
        # forward from ISR using subset of exposure to end
        # high probability for last week, drops off for older

        path = [(linearity1,0), (addmaskattr1,0), (badpixrm1,0), (rmoverscan1,0), (biascorrect1,0)]
        path = [(linearity1,0), (addmaskattr1,0), (badpixrm1,0), (overscanimg1,0),
                (meanoverscan1,0), (rmoverscan1, 1)]        

        
        mlog.info("Scenario 2")
        region = [(x,y) for x in xrange(0, 10) for y in xrange(0, 10) ]
        #log_prov(mlog, w.forward_path(region, runid, path))
        queries = []
        queries.append([ region, runid, path, 'forward' ])
        return queries
        


    def scenario3(runid):
        # # scenario 3: fast pass to ensure no false positives
        # # backward from a random subset of sources
        # # high prob for current night, drops zipfian further back
        mlog.info("Scenario 3")
        path1 = [(clust,0), (rmcr,1), (maskexp,0), (cr1, 0), (rmbg1, 0), (flatcorrect1, 0)]
        path2 = [(clust,0), (rmcr,1), (maskexp,1), (cr2, 0), (rmbg2, 0), (flatcorrect2, 0)]
        
        output = clust.wrapper.get_output(runid)
        indices = np.argwhere(output).tolist()
        random.shuffle(indices)
        queries = []
        for star in indices[:5]:
            queries.append([ [ star ], runid, path1, 'backward'])
            queries.append([ [ star ], runid, path2, 'backward'])
        return queries
        
    def scenario4(runid):
        # # scenario 4: a star just disappeared.  what happened?
        # # backward from pixel that should contain a star but doesn't
        # # low probability overall, higher for current night.
        mlog.info("Scenario 4")
        path1 = [(rmcr,1)]#, (maskexp,0), (cr1, 0), (rmbg1, 0), (flatcorrect1, 0),
        #(darkcorrect1, 0), (biascorrect1, 0), (rmoverscan1, 0), (badpixrm1, 0),
         #       (addmaskattr1, 0), (linearity1, 0)]
        path2 = [(rmcr,1)]#, (maskexp,1), (cr2, 0), (rmbg2, 0), (flatcorrect2, 0),
#                 (darkcorrect2, 0), (biascorrect2, 0), (rmoverscan2, 0), (badpixrm2, 0),
 #                (addmaskattr2, 0), (linearity2, 0)]

        pixels = [ (x,y) for x in xrange(25, 50) for y in xrange(30, 50)  ]
        queries = []
        queries.append([pixels, runid, path1, 'backward'])
        queries.append([pixels, runid, path2, 'backward'])
        return queries


    def execute(exp1, exp2, runprovqs=True):
        global runmode
        global darkmask1, biasmask1, satmask1, flatmask1
        global darkmask2, biasmask2, satmask2, flatmask2        
        if runmode == 0:
            exp1 = exp2 = numpy.array([[0,0,1,1],
                                       [0,0,0,1],
                                       [1,0,0,0],
                                       [1,1,0,0]], dtype=np.float)
            thresh = np.array([[0.5, 0.5]])
        elif runmode == 1:
            thresh = np.array([[100,100]])
        elif runmode == 2:
            print "runmode 2 not implemented"
        elif runmode >= 3:
            if runmode == 3: # lots of large stars
                exp1 = exp2 = gen_data(smallshape[0], smallshape[1], nstars=10, starradius=4)
            elif runmode == 4: # small number of large stars
                exp1 = exp2 = gen_data(smallshape[0], smallshape[1], nstars=2, starradius=4)
            elif runmode == 5: # lots of tiny stars
                exp1 = exp2 = gen_data(smallshape[0], smallshape[1], nstars=10, starradius=1)
            elif runmode == 6: # small number of tiny stars
                exp1 = exp2 = gen_data(smallshape[0], smallshape[1], nstars=2, starradius=1)
            elif runmode == 7: # lots of gigantic stars                
                exp1 = exp2 = gen_data(smallshape[0], smallshape[1], nstars=10, starradius=10)
            thresh = np.array([[100,100]])


        satmask1 = biasmask1 = darkmask1 = np.zeros(exp1.shape, bool)
        flatmask1 = np.ones(exp1.shape, bool)
        satmask2 = biasmask2 = darkmask2 = np.zeros(exp2.shape, bool)
        flatmask2 = np.ones(exp2.shape, bool)
        run_img(exp1, exp2, thresh)
        if runprovqs:
            run_prov_workload()

        del exp1, exp2
        del satmask1, biasmask1, darkmask1, flatmask1
        del satmask2, biasmask2, darkmask2, flatmask2        
        


    def run_iterations(logdir, runtype, scidata1, scidata2, runprovqs=False):
        global mlog, w

        handler = logging.FileHandler('%s/provres_%s.log' % (logdir, runtype), 'w')
        handler.setLevel(logging.DEBUG)
        mlog.addHandler(handler)

        iters = 0
        for blockx in xrange(factor[0]):
            for blocky in xrange(factor[1]):
                x = smallshape[0] * blockx
                y = smallshape[1] * blocky
                exp1 = scidata1[x:x+smallshape[0],y:y+smallshape[1]]
                exp2 = scidata2[x:x+smallshape[0],y:y+smallshape[1]]
                execute(exp1, exp2, runprovqs)
                iters += 1
                if iters >= maxiters: break
            if iters >= maxiters: break


        mlog.removeHandler(handler)
        handler.close()

        handler = logging.FileHandler('%s/perf_%s.log' % (logdir, runtype), 'w')
        handler.setLevel(logging.DEBUG)
        mlog.addHandler(handler)


        mlog.removeHandler(handler)
        handler.close()


    def set_all_strategy(strat, w):

        def set_strategy(wop):
            if Mode.FULL_MAPFUNC in wop.op.supported_modes():
                Runtime.instance().set_strategy(wop.op, Strat.full())
            else:
                Runtime.instance().set_strategy(wop.op, strat)
        w.visit(set_strategy)


    def onetoone(x):
        return x
    linearity1 = OneToOne(onetoone)
    addmaskattr1 = OneToOne(onetoone)
    badpixrm1 = BadPixRemoval()
    overscanimg1 = Subsample(box=numpy.array([[0,1],[0,1]]))
    meanoverscan1 = MeanSingleVal()
    rmoverscan1 = DiffSingleVal()
    biascorrect1 = Diff()
    darkcorrect1 = Diff()
    flatcorrect1 = Diff()
    bgmean1 = Mean()
    rmbg1 = Diff()
    cr1 = CRDetect()
    

    linearity2 = OneToOne(onetoone)
    addmaskattr2 = OneToOne(onetoone)
    badpixrm2 = BadPixRemoval()
    overscanimg2 = Subsample(box=numpy.array([[0,1],[0,1]]))
    meanoverscan2 = MeanSingleVal()
    rmoverscan2 = DiffSingleVal()
    biascorrect2 = Diff()
    darkcorrect2 = Diff()
    flatcorrect2 = Diff()
    bgmean2 = Mean()
    rmbg2 = Diff()
    cr2 = CRDetect()

    sumexp = Add()
    maskexp = LogicalAnd()
    rmcr = RemoveCRs()
    clust = Cluster()
    

    w = Workflow()
    # 1st exposure
    w.register(linearity1, 1)
    w.register(addmaskattr1, 1)
    w.register(badpixrm1, 2)
    w.register(overscanimg1, 1)
    w.register(meanoverscan1, 1)
    w.register(rmoverscan1, 2)  # array and meanoverscan
    w.register(biascorrect1, 2) 
    w.register(darkcorrect1, 2)
    w.register(flatcorrect1, 2)
    w.register(bgmean1, 1)
    w.register(rmbg1, 2)
    w.register(cr1, 2)
    # 2nd exposure
    w.register(linearity2, 1)
    w.register(addmaskattr2, 1)
    w.register(badpixrm2, 2)
    w.register(overscanimg2, 1)
    w.register(meanoverscan2, 1)
    w.register(rmoverscan2, 2)  # array and meanoverscan
    w.register(biascorrect2, 2) 
    w.register(darkcorrect2, 2)
    w.register(flatcorrect2, 2)
    w.register(bgmean2, 1)
    w.register(rmbg2, 2)
    w.register(cr2, 2)
    # merge the exposures
    w.register(sumexp, 2)
    w.register(maskexp, 2)    
    w.register(rmcr, 2)
    w.register(clust, 2)  # array and threshold

    w.connect(linearity1, addmaskattr1, 0)
    w.connect(addmaskattr1, badpixrm1, 0)    
    w.connect(badpixrm1, overscanimg1, 0)
    w.connect(overscanimg1, meanoverscan1, 0)
    w.connect(badpixrm1, rmoverscan1, 0)    
    w.connect(meanoverscan1, rmoverscan1, 1)
    w.connect(rmoverscan1, biascorrect1, 0)
    w.connect(biascorrect1, darkcorrect1, 0)
    w.connect(darkcorrect1, flatcorrect1, 0)    
    w.connect(flatcorrect1, bgmean1, 0)
    w.connect(flatcorrect1, rmbg1, 0)
    w.connect(bgmean1, rmbg1, 1)
    w.connect(rmbg1, cr1, 0)

    w.connect(linearity2, addmaskattr2, 0)
    w.connect(addmaskattr2, badpixrm2, 0)    
    w.connect(badpixrm2, overscanimg2, 0)
    w.connect(overscanimg2, meanoverscan2, 0)
    w.connect(badpixrm2, rmoverscan2, 0)    
    w.connect(meanoverscan2, rmoverscan2, 1)
    w.connect(rmoverscan2, biascorrect2, 0)
    w.connect(biascorrect2, darkcorrect2, 0)
    w.connect(darkcorrect2, flatcorrect2, 0)    
    w.connect(flatcorrect2, bgmean2, 0)
    w.connect(flatcorrect2, rmbg2, 0)
    w.connect(bgmean2, rmbg2, 1)
    w.connect(rmbg2, cr2, 0)

    w.connect(rmbg1, sumexp, 0)
    w.connect(rmbg2, sumexp, 1)
    w.connect(cr1, maskexp, 0)
    w.connect(cr2, maskexp, 1)
    w.connect(sumexp, rmcr, 0)
    w.connect(maskexp, rmcr, 1)    
    w.connect(rmcr, clust, 0)

    w.default_strategy()

    # save = False
    # if save:
    #     run_iterations(logdir, runtype, scidata1, scidata2, False)
    #     w.save('./w.state')
    # else:
    #     w.load('./w.state')
    #     run_prov_workload()

    # exit()


    # noop
    # set_all_strategy(Strat.noop(), w)
        
    # Stats.instance().add_exec(smallshape[0], smallshape[1],
    #                           runmode, 'noop', logdir, "lsst_noop")
    # run_iterations(logdir, runtype, scidata1, scidata2, False)
    # Stats.instance().finish_exec()    

    # # stats
    # set_all_strategy(Strat.stat(), w)
        
    # Stats.instance().add_exec(smallshape[0], smallshape[1],
    #                           runmode, 'stats', logdir, "lsst_stats")
    # run_iterations(logdir, runtype, scidata1, scidata2, False)
    # Stats.instance().finish_exec()

    # # query everything
    # set_all_strategy(Strat.query(), w)
    # for op in [cr1, cr2, clust, rmcr]:
    #     Runtime.instance().set_strategy(op, Strat.query())
        
    # Stats.instance().add_exec(smallshape[0], smallshape[1],
    #                           runmode, 'query', logdir, "lsst_queryall")
    # run_iterations(logdir, runtype, scidata1, scidata2, False)

    # for maxcount in [-1]:#, 1, 10, 100]:
    #     print "query results: %d\t%f" % (maxcount, run_prov_workload(maxcount))
    # Stats.instance().finish_exec()        

    # # PSET everything
    # set_all_strategy(Strat.query(), w)
    # for op in [cr1, cr2, clust, rmcr]:
    #     Runtime.instance().set_strategy(op, Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY)))
        
    # Stats.instance().add_exec(smallshape[0], smallshape[1],
    #                           runmode, 'one_key', logdir, "lsst_one_key")
    # run_iterations(logdir, runtype, scidata1, scidata2, False)
    # for maxcount in [-1]:#, 1, 10, 100]:
    #     print "query results: %d\t%f" % (maxcount, run_prov_workload(maxcount))
    # Stats.instance().finish_exec()


    # many_key everything
    set_all_strategy(Strat.query(), w)
    for op in [cr1, cr2, clust, rmcr]:
        Runtime.instance().set_strategy(op, Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY)))
        
    Stats.instance().add_exec(smallshape[0], smallshape[1],
                              runmode, 'many_key', logdir, "lsst_many_key")
    run_iterations(logdir, runtype, scidata1, scidata2, False)
    for maxcount in [-1]:#, 1, 10, 100]:
        print "query results: %d\t%f" % (maxcount, run_prov_workload(maxcount))
    Stats.instance().finish_exec()



    # optimal backwards
    set_all_strategy(Strat.query(), w)
    FDIFF = Strat([Bucket([ Desc(Mode.FULL_MAPFUNC, Spec.default(), True),
                            Desc(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True) ])])

    Runtime.instance().set_strategy(cr1, FDIFF)
    Runtime.instance().set_strategy(cr2, FDIFF)
    Runtime.instance().set_strategy(clust, Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY)))
    Runtime.instance().set_strategy(rmcr, FDIFF)

    Stats.instance().add_exec(smallshape[0], smallshape[1],
                              runmode, 'opt_manual', logdir, "lsst_opt_manual")
    run_iterations(logdir, runtype, scidata1, scidata2, False)
    for maxcount in [-1]:#, 1, 10, 100]:
        print "query results: %d\t%f" % (maxcount, run_prov_workload(maxcount))

    write_fits(rmbg1.wrapper.get_output(w._runid-1), './generated.fits')
    write_fits(rmcr.wrapper.get_output(w._runid-1), './generated2.fits')    
    Stats.instance().finish_exec()
    exit()
    
    


    if runtype == 'opt':
        set_all_strategy(Strat.noop(),w)
        run_iterations(logdir, runtype, scidata1, scidata2, False)
        run_id = w._runid-1

        
        noop_eids = Stats.instance().get_matching_noops(runmode, smallshape)
        print "noop eids:", noop_eids
        queries = gen_prov_workload(run_id)
            
        for q in queries:
            q[0] = len(q[0])
        mp = ModelPredictor(noop_eids, w, queries)
        
        disksizes = [1, 10, 100]
        runconstraint = 100

        for disksize in disksizes:
            set_all_strategy(Strat.single(Mode.FULL_MAPFUNC, Spec.default()), w)

            # check if we can load a strategy from the current directory
            # assumes TEST + optimize was run already
            strategies = run_nlp(Stats.instance(), w, mp, disksize, runconstraint)

            Stats.instance().add_exec(smallshape[0], smallshape[1],
                                      runmode, runtype, logdir, "lsst_full",
                                      disksize, runconstraint, 0)
            Stats.instance().add_opt_mappings(strategies)

            set_all_strategy(Strat.query(), w)
            for op, strats in strategies.items():
                for s in strats:
                    if s != Strat.full():
                        print op, s
                    Runtime.instance().set_strategy(op, s)


            run_iterations(logdir, runtype, scidata1, scidata2, True)
            Stats.instance().finish_exec()

    else:
        Stats.instance().add_exec(smallshape[0], smallshape[1],
                                  runmode, runtype, logdir,
                                  "lsst_full")
        raise RuntimeError, "not implemented"
        strat = Strategy.from_str(runtype)
        set_all_strategy(strat, w)


        # cache settings
        # 
        # cache all
        # no caching
        # no caching + top K
        
        #w.load('./w.state')
        #run_prov_workload()
        

        # if strat == STRAT_NOOP or strat == STRAT_F or strat == STRAT_STATS:
        #     run_iterations(logdir, runtype, scidata1, scidata2, False)            
        # else:
        #     run_iterations(logdir, runtype, scidata1, scidata2, True)

    Stats.instance().finish_exec()
    Stats.instance().close()

    
