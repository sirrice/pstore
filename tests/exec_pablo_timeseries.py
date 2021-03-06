import sys
sys.path.append('../')
import numpy as np
import math
from pablo import *
from common import *
from strat import *
from modelpredict import *



model_type = "D"


target_feature = "RELAPSE"
clinical = ["HISTO", "MET"]
pathway = ["YU_CMYC"]
pathway_dep_subclass = [
       "TRANSLATION_FACTORS",                             # c1
       "BRCA_BRCA1_POS",                                  # c1
       "FLUMAZENILPATHWAY",                               # c1
       "AGUIRRE_PANCREAS_CHR8",                           # c1
       "MENSSEN_MYC_UP",                                  # c1
       "PENG_RAPAMYCIN_DN",                               # c1
       "DSRNA",                                           # c2
       "HINATA_NFKB_DN",                                  # c2
       "P53_DN.v2",                                       # c2
       "IL22BPPATHWAY",                                   # c2
       "ABBUD_LIF_GH3_DN",                                # c2
       "NGUYEN_KERATO",                                   # c2
       "POMEROY_DESMOPLASIC_VS_CLASSIC_MD_UP",            # c3
       "HSA04340_HEDGEHOG_SIGNALING_PATHWAY",             # c3
       "HINATA_NFKB",                                     # c3
       "RADAEVA_IFNA_DN",                                 # c3
       "RORIE_ES_PNET_UP",                                # c3
       "SARCOMAS_LIPOSARCOMA",                            # c3
       "GPCRS_CLASS_C_METABOTROPIC_GLUTAMATE_PHEROMONE",  # c4
       "DFOSB_BRAIN_2WKS_UP",                             # c4
       "PARKINPATHWAY",                                   # c4
       "HSA05020_PARKINSONS_DISEASE",                     # c4
       "MITRPATHWAY",                                     # c4
       "HSA04080_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION", # c4
       "PHOTO_DN.v1",                                     # c5
       "CRX_DN.v1",                                       # c5
       "CtIP_DN.v1",                                      # c5
       "NRL_DN.v1",                                       # c5
       "WNT_UP.v1",                                       # c5
       "CHREBPPATHWAY",                                   # c5
       "NEUROTRANSMITTERSPATHWAY",                        # c6
       "BCAT_UP.v1_UP",                                   # c6
       "ST_WNT_BETA_CATENIN_PATHWAY",                     # c6
       "LEF1_UP.v1_UP",                                   # c6
       "TGFBETA_LATE_UP",                                 # c6
       "UVB_NHEK3_C6"]                                    # c6
subclass = "C6_CLUSTER"
pathway_dep_outcome = [ "mTOR_UP.v1",               # c1
               "HOGERKORP_CD44_UP",        # c2
               "COLLER_MYC",               # c3
              "HISTIDINE_METABOLISM",      # c4
              "GLI1_UP.v1_DN",             # c5
             "RIBAVIRIN_RSV"]              # c6

genom = ["Amp_8q24.21",      # C-MYC  (score: )
         "Amp_2p24.3",       # N-MYC
         "Del_6q",           # monosomy 6
         "Del_16q",
         "Del_16q23.3",
         "Amp_7q21.3",
         "Amp_3q26.32"]


features = []
for feats in (target_feature, clinical, pathway, subclass,
              pathway_dep_outcome, genom):
    if isinstance(feats, list):
        features.extend(feats)
    else:
        features.append(feats)
ftypes = []
ftypes += ['target.feature'] 
ftypes += ['clinical'] * len(clinical)
ftypes += ['pathway'] * len(pathway)
ftypes += ['subclass']
ftypes += ['pathway.dep.outcome'] * len(pathway_dep_outcome)
ftypes += ['genom'] * len(genom)

sfeatures = [subclass] + pathway_dep_subclass
stypes = ['subclass'] + ['pathway.dep.subclass'] * len(pathway_dep_subclass)


features = np.array([features])
ftypes = np.array([ftypes])


def create_workflow():

    gn = GetNames()
    ss = Subsample(((0, ds.shape[0]), (4, ds.shape[1])))
    tr = Transpose()
    en = ExtractNonsubclass()
    cm = CreateModel()
    pr = Predict()
    cum = CumOdds()
    prob = OneToOne( lambda lo: 1 / (1 + math.exp(-lo)) )
    klass = GetClasses()
    val = Validate()




    w = Workflow()
    w.register(gn, 1) # names
    w.register(ss, 1)
    w.register(tr, 1) # m2
    w.register(en, 3) # md
    w.register(cm, 4) # ev_list
    w.register(pr, 5)
    w.register(cum, 1)
    w.register(prob, 1)
    w.register(klass, 1)
    w.register(val, 2)


    w.connect(ss, tr, 0)
    w.connect(tr, en, 0)
    w.connect(gn, en, 1)
    w.connect(en, cm, 0)
    w.connect(tr, cm, 1)
    w.connect(gn, cm, 2)
    w.connect(cm, pr, 0)
    w.connect(en, pr, 1)
    w.connect(tr, pr, 2)
    w.connect(gn, pr, 4)
    w.connect(pr, cum, 0)
    w.connect(cum, prob, 0)
    w.connect(prob, klass, 0)
    w.connect(klass, val, 0)
    w.connect(tr, val, 1)

    def run_model(ds, runmode, runtype, disk=0, runcost=0, eids = [-1]):
        runtype = '%s_m' % runtype
        eids = Stats.instance().get_matching_noops(runmode, ds.shape) or [-1]
        Stats.instance().add_exec(ds.shape[0], ds.shape[1],
                                  runmode, runtype, './_output',
                                  "notes", disk, runcost, eids[0])
        eid = Stats.instance().eid

        mp = ModelPredictor(eids, w)
        for op, s in Runtime.instance().cur_strats.items():
            disk = mp.get_disk(op,s) * 1048576.0
            overhead = mp.get_provcost(op, s)
            opcost = mp.get_opcost(op, s)

            Stats.instance().add_modelrun(eid, op, s, disk, overhead, eids)

        qs = map(list, get_qs())
        for q in qs:
            coords, runid, path, direction = q
            mp = ModelPredictor(eids, w)
            qcost = 0.0
            for op, s in Runtime.instance().cur_strats.items():
                qcost += mp.get_pqcost(op,s)
            nres = 0
            path_ops = [x[0] for x in path]
            Stats.instance().add_pq(runid, path_ops, direction, [len(coords)], qcost, nres) 

        Stats.instance().finish_exec()   





    def run(ds, runmode, runtype, disk=0, runcost=0, eids=[-1]):
        if runtype == 'stats':
            w.default_strategy(userstrat=Strat.single(Mode.STAT, Spec.default()))

        
        random.seed(0)
        Stats.instance().add_exec(ds.shape[0], ds.shape[1],
                                  runmode, runtype, './_output',
                                  "notes", disk, runcost, eids[0])
        
        w.add_static_input(gn, 0, ds)
        w.add_static_input(ss, 0, ds)
        w.add_static_input(en, 2, features)
        w.add_static_input(cm, 3, ftypes)
        w.add_static_input(pr, 3, ftypes)
        
        w.run()

        Stats.instance().finish_exec()


    def get_fqs(runid):
        qs = []
        # bad hospital query
        path = [ (ss, 0), (tr, 0), (en, 0), (cm, 0)]
        qs.append( [ [(42, random.randint(1, 10))], runid, path, 'forward' ] )

        path = [ (ss, 0), (tr, 0), (en, 0), (cm, 0), (pr, 0), (cum, 0), (prob, 0), (klass, 0), (val, 0) ]
        #qs.append( [ [(j, i) for i in xrange(10, 15) for j in (5, 42)], runid, path, 'forward' ] )
        minv = random.randint(5, 32)
        qs.append( [ [(j, i) for i in xrange(10, 15) for j in xrange(minv, minv + 10)], runid, path, 'forward' ] )
        return qs

    def get_bqs(runid):
        qs = []
        # high likelihood query
        path = [ (cm,0), (en, 0), (tr, 0), (ss, 0) ]
        qs.append([ [(0, random.randint(0, 50))], runid, path, 'backward' ])

        # Relapse query
        path = [ (klass, 0), (prob, 0), (cum, 0), (pr, 0), (cm,0), (en, 0), (tr, 0), (ss, 0) ]
        qs.append([ [(0, random.randint(0, 50))], runid, path, 'backward' ])        
        return qs

    def get_qs(iteridx=3):
        random.seed(iteridx)
        runid = w._runid-1
        qs = []

        for i in xrange(5):
            if i < iteridx:
                qs.extend(get_fqs(runid))
            else:
                qs.extend(get_bqs(runid))

        # === Forward Queries ===


        # path = [ (ss, 0), (tr, 0), (cm, 1), (pr, 0), (cum, 0), (prob, 0), (klass, 0) ]#, (val, 0) ]
        # for i in xrange(5,10):
        #     qs.append( [ [(42, i)], runid, path, 'forward' ] )

        # path = [ (gn, 0), (en, 1), (cm, 0), (pr, 0), (cum, 0), (prob, 0), (klass, 0) ] #, (val, 0) ]
        # coords = [ (i, 0)  for i in xrange(20)]
        # qs.append( [ coords, runid, path, 'forward' ] )

        # path = [ (gn, 0) , (cm, 2), (pr, 0), (cum, 0), (prob, 0), (klass, 0) ]#, (val, 0) ]
        # coords = [ (i, 0)  for i in xrange(20)]
        # qs.append( [ coords, runid, path, 'forward' ] )


        # path = [ (pr, 2), (cum, 0), (prob, 0), (klass, 0) ]#, (val, 0) ]
        # coords = [ (i, 0)  for i in xrange(20)]
        # qs.append( [ coords, runid, path, 'forward' ] )


        #
        # backward queries
        #
        

        # path = [(klass, 0), (prob, 0), (cum, 0), (pr, 0), (cm,0) ]
        # for i in xrange(4):
        #     qs.append([ [(0, i)], runid, path, 'backward' ])

        # path = [(klass, 0), (prob, 0), (cum, 0), (pr, 1), (en, 0), (tr, 0) ]
        # for i in xrange(4):
        #     qs.append([ [(0, i)], runid, path, 'backward' ])

            
            

        # path = [(pr, 0), (cm,1), (tr, 0), (ss, 0) ]
        # for i in xrange(10):
        #     qs.append([ [(random.randint(0,50), random.randint(0,50))  ], runid, path, 'backward' ])

        # path = [(pr, 0), (cm,0), (en, 0), (tr, 0), (ss, 0) ]
        # qs.append([ [(random.randint(0,50), random.randint(0,50)) for i in xrange(10) ],
        #             runid, path, 'backward' ])

        # path = [(val, 0), (klass, 0), (prob, 0), (cum, 0), (pr, 0), (cm,0), (en, 0), (tr, 0), (ss, 0) ]
        # qs.append([ [(0, 0)], runid, path, 'backward' ])        
            
        return qs

    def get_strats():

        def set_all(s):
            w.default_strategy(userstrat=s)

        def noop():
            set_all(Strat.noop())
            return 'noop'

        def stat():
            set_all(Strat.single(Mode.STAT, Spec.default()))
            return 'stats'

        def query_all():
            set_all(Strat.single(Mode.QUERY, Spec.default()))
            return 'q_all'

        def query_opt():
            for op in w.ops.keys():
                if Mode.FULL_MAPFUNC in op.supported_modes():
                    Runtime.instance().set_strategy(op, Strat.single(Mode.FULL_MAPFUNC, Spec.default()))
                else:
                    Runtime.instance().set_strategy(op, Strat.single(Mode.QUERY, Spec.default()))
            return 'q_opt'

        def set_ptr_wrapper_opt(s):
            for op in w.ops.keys():
                matches = True
                for mode in s.modes():
                    if mode not in op.supported_modes():
                        matches = False

                if Mode.FULL_MAPFUNC in op.supported_modes():
                    Runtime.instance().set_strategy(op, Strat.single(Mode.FULL_MAPFUNC, Spec.default()))
                elif matches:
                    Runtime.instance().set_strategy(op, s)
                else:
                    Runtime.instance().set_strategy(op, Strat.single(Mode.QUERY, Spec.default()))


        def set_ptr_wrapper(s):
            for op in w.ops.keys():
                matches = True
                for mode in s.modes():
                    if mode not in op.supported_modes():
                        matches = False

                if matches:
                    Runtime.instance().set_strategy(op, s)
                else:
                    Runtime.instance().set_strategy(op, Strat.single(Mode.QUERY, Spec.default()))


        def pt1():
            strat = Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True)
            set_ptr_wrapper_opt(strat)
            return '2_ONE_KEY_B'

        def pt2():
            strat = Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.BINARY), True)
            set_ptr_wrapper_opt(strat)
            return '2_MANY_MANY_B'

        def pt3():
            buckets = [Bucket([Desc(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True)]),
                       Bucket([Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False)]) ]
            s = Strat(buckets)
            set_ptr_wrapper_opt(s)
            return '2_F_B'

        def pt4():
            buckets = [Bucket([Desc(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True)]),
                       Bucket([Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False)]) ]
            s = Strat(buckets)
            set_ptr_wrapper_opt(s)
            return '2_F_B'


        def ptr1():
            strat = Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True)
            set_ptr_wrapper_opt(strat)
            return '3_ONE_KEY_B'

        def ptr2():
            strat = Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True)
            set_ptr_wrapper_opt(strat)
            return '3_MANY_MANY_B'

        def ptr25():
            strat = Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.GRID), True)
            set_ptr_wrapper_opt(strat)
            return '3_MANY_GRID_B'

        def ptr3():
            strat = Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False)
            set_ptr_wrapper_opt(strat)
            return '3_ONE_KEY_F'

        def ptr4():
            strat = Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False)
            set_ptr_wrapper_opt(strat)
            return '3_MANY_MANY_F'

        def ptr5():
            buckets = [Bucket([Desc(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True)]),
                       Bucket([Desc(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False)])]
            s = Strat(buckets)
            set_ptr_wrapper_opt(s)
            return '3_F_B'

        def ptr6():
            buckets = [Bucket([Desc(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True)]),
                       Bucket([Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False)])]
            s = Strat(buckets)
            set_ptr_wrapper_opt(s)
            return '3_F_B_2'

        def custom():
            mkf = Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False)
            pt = Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True)
            Runtime.instance().set_strategy(cm , mkf)
            Runtime.instance().set_strategy(val , mkf)
            Runtime.instance().set_strategy(pr , pt)
            return "CUSTOM"

        def opt(ds, eids, runmode, disk, runcost):
            mp = ModelPredictor(eids, w, disk, runcost)
            strategies, torm = run_nlp(Stats.instance(), w, mp, disk,runcost)
            for op in sorted(strategies.keys()):
                Runtime.instance().set_strategy(op, strategies[op][0])

            # remove old strategies from cache
            for op, runid in torm:
                b = Runtime.instance().delete_pstore(op, runid)

            return 'opt_%.1f_%d' % (disk, runcost)

        return [opt, query_opt, noop]

    
    return w, run, run_model, get_strats, get_qs

def run_qs(w, qs, bmodel):
    ret = []
    for coords, runid, path, direction in qs:
        start = time.time()        
        if direction == 'forward':
            res,optcost,qsizes = w.forward_path(coords, runid, path)
        elif direction == 'backward':
            res,optcost,qsizes = w.backward_path(coords, runid, path)

        nres = 0
        for coord in res:
            #print coord
            nres += 1
        #nres = len(res)

        end = time.time()
        qcost = end - start - optcost
        
        path_ops = [x[0] for x in path]
        pqid = Stats.instance().add_pq(runid, path_ops, direction, [len(coords)], qcost, nres)
        for (op, arridx), (opcost, strat, insize, outsize) in qsizes.items():
            Stats.instance().add_iq(pqid, op, arridx, strat, direction, insize, outsize, opcost)

        ret.append(  (nres, qcost) )
    return ret





if __name__ == '__main__':
    fname = '../data/MD_train_set.txt'
    runmode = 1
    bmodel = bdynamic = False    
    dbname = '_output/pablostats.db'

    for arg in sys.argv:
        if arg == 'model':
            bmodel = True
        elif arg == 'dynamic':
            bdynamic = True
        elif '.db' in arg:
            dbname = arg
        else:
            try:
                runmode = int(arg)
                if runmode > 1:
                    fname = '../data/MD_train_set_%dx.txt' % runmode
            except:
                pass

    with file(fname, 'r') as f:
        ds = np.array([l.strip().split('\t') for l in f])[1:,:]
    print dbname, bmodel, bdynamic, runmode
    Stats.instance(dbname)
    Stats.instance().timeseries = True
    w, run, run_model, get_strats, get_qs = create_workflow()


    def run_opt(ds, runmode, runtype, set_strat, get_qs, bmodel=False):
        basesize = ds.shape[0] * ds.shape[1] * 8 / 1048576.0 * 2
        disksizes = [1, 10,20,50,100]
        #disksizes = [100]
        runcost = 10000
        
        w.boptimize = bdynamic
        
        for disk in disksizes:
            Runtime.instance().restore_pstores() # this resets the experiment
            qmixes = [0,0,1,3,5,5,5,5,5,5]
            #qmixes = [5,5,5,5]

            # gotta get some stats
            run(ds, runmode, 'stats', disk, runcost, [-1])
            eids = Stats.instance().get_matching_noops(runmode, ds.shape, disk, runcost)
            print 'eids', eids
            
            
            for iteridx in xrange(len(qmixes)):
                qmix = qmixes[iteridx]
                runtype = set_strat(ds, eids, runmode, disk, runcost)
                print runtype, disk
                if bmodel:
                    run_model(ds, runmode, runtype, disk, runcost, eids)
                else:
                    run(ds, runmode, runtype, disk, runcost, eids)
                    totaldisk = Runtime.instance().get_total_disk()
                    Stats.instance().add_exec_stats(totaldisk)
                    
                    mp = ModelPredictor(eids, w, disk, runcost)
                    w.mp = mp
                    
                    for x in run_qs(w, get_qs(qmix), bmodel):
                        print x
                    print
        w.boptimize = False
        w.mp = None

    def run_normal(ds, runmode, runtype, set_strat, get_qs, bmodel=False):
        runtype = set_strat()
        print runtype

        eids = Stats.instance().get_matching_noops(runmode, ds.shape)
        mp = ModelPredictor(eids, w)
        w.boptimize = bdynamic
        w.mp = mp

        # print "disk\t", sum( [mp.get_disk(op, strat) for op, strat in Runtime.instance().cur_strats.items()] )
        # print "cost\t", sum( [mp.get_pqcost(op, strat) for op, strat in Runtime.instance().cur_strats.items()] )
        # #return
        if bmodel:
            run_model(ds, runmode, runtype)
        else:
            run(ds, runmode, runtype)
            if runtype in ('noop', 'stats'):
                return
            for x in run_qs(w, get_qs(), bmodel):
                print x
            print


    runtype = None
    for set_strat in get_strats():
        try:
            runtype = set_strat()
            run_all = run_normal
        except:
            run_all = run_opt

        run_all(ds, runmode, runtype, set_strat, get_qs, bmodel)
