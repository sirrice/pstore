import numpy as np
from models import *
from stats import *
from operator import mul
from util import zipf


class ModelPredictor(object):

    def __init__(self, eids, workflow, diskconstraint=None, runconstraint=None, l = 1.9):
        self.workflow = workflow
        self.counts = {}
        self.cache = {}
        self.fprobs = {}
        self.bprobs = {}
        self.qprobs = {}
        self.fqsizes = {}
        self.bqsizes = {}
        self.opqcount = {}
        self.eids = eids

        self.qfanouts = {}  # forward qs.  (op, arridx) -> query fanout.  default = 100000
        self.qfanins = {}   # backward qs. (op, arridx) -> query fanout.  default = 100000

        qeids = set()
        for eid in self.eids:
            qeids.update(Stats.instance().get_similar_eids(eid, diskconstraint, runconstraint))
        self.qeids = qeids
        self._cache_stats()
        if len(qeids) > 0:
            self._cache_fanouts(qeids, diskconstraint, runconstraint)
            self._calc_probs(qeids, diskconstraint, runconstraint)



        keys = set(self.fqsizes.keys())
        keys.update(self.bqsizes.keys())
        for key in keys:
            print "qfanouts\t", key[0], key[1], '\t', \
                '%s\t%s\t%.3f\t%.3f' % (str(self.fqsizes.get(key,-1)).ljust(10),
                                        str(self.bqsizes.get(key,-1)).ljust(10),
                                        self.fprobs.get(key, -1),
                                        self.bprobs.get(key, -1))
        print 

        cumprobs = zipf(workflow._runid, l)
        probs = []
        prob = 0.0
        for p in cumprobs:
            probs.append(p-prob)
            prob = p
        self.probs = dict([(workflow._runid-i, p) for i, p in enumerate(probs)])

        self.debug = True

    def _wma(self, qsizes, alpha=0.1, default=1000000):
        if len(qsizes) == 0:
            return default
        if qsizes[0] is None:
            return self._wma(qsizes[1:], alpha, default)
        if len(qsizes) == 1:
            return qsizes[0]
        return (1 - alpha) * qsizes[0] + alpha * self._wma(qsizes[1:], alpha, default)


    def _calc_probs(self, qeids, diskc, runc):
        qeids = sorted(list(qeids), reverse=True)
        qeids.remove(self.workflow._runid-1)
        keys = []
        allfstats = []
        allbstats = []
         
        def f(w):
            op = w.op
            for arridx in xrange(w.nargs):
                key = (op, arridx)
                fstats, bstats = Stats.instance().get_iq_stat(qeids, op.oid, arridx)
                # if 'Extract' in str(op) and 0 == arridx:
                #     print "EXTRACT\t", map(lambda s: s[1], fstats)
                #     print "EXTRACT\t", map(lambda s: s[1], bstats)
                allfstats.append(fstats)
                allbstats.append(bstats)
                keys.append(key)

        self.workflow.visit(f)

        allfprobs, allbprobs, allqprobs, allcounts = [], [], [], []
        for tstep in xrange(len(allfstats[0])):
            fcounts = [stats[tstep][1] for stats in allfstats]
            bcounts = [stats[tstep][1] for stats in allbstats]

            fsum = sum(fcounts)
            bsum = sum(bcounts)
            total = float(fsum + bsum)
                
            fprobs = [fc+bc > 0 and fc / float(fc+bc) or 0 for fc, bc in zip(fcounts, bcounts)]
            bprobs = [fc+bc > 0 and bc / float(fc+bc) or 0 for fc, bc in zip(fcounts, bcounts)]
            counts = [fc+bc for fc, bc in zip(fcounts, bcounts)]
            qprobs = [total > 0 and (fc+bc) / total or 0 for fc, bc in zip(fcounts, bcounts)]
            allfprobs.append(fprobs)
            allbprobs.append(bprobs)
            allqprobs.append(qprobs)
            allcounts.append(counts)            
            
        for keyidx in xrange(len(keys)):
            key = keys[keyidx]
            fprobs = [row[keyidx] for row in allfprobs]
            bprobs = [row[keyidx] for row in allbprobs]
            qprobs = [row[keyidx] for row in allqprobs]
            counts = [row[keyidx] for row in allcounts]
            fprob = self._wma(fprobs, default=0)
            bprob = self._wma(bprobs, default=0)
            qprob = self._wma(qprobs, default=0)
            count = self._wma(counts, default=1)
            self.fprobs[key] = fprob
            self.bprobs[key] = bprob
            self.qprobs[key] = qprob
            self.opqcount[key] = count
            # if 'Extract' in str(key[0]) and 0 == key[1]:
            #     print "EXTRACT\t",fprobs
            #     print "EXTRACT\t",bprobs
            #     print "EXTRACT\t",qprobs
            #     print "EXTRACT\t", counts
                
                
        # # probability of forward query
        # self.fprobs = dict([(key, val[0] / float(sum(val))) for key, val in self.counts.items() if sum(val) > 0])
        # self.bprobs = dict([(key, val[1] / float(sum(val))) for key, val in self.counts.items() if sum(val) > 0])


        # # probability of querying an op
        # total = sum(map(sum, self.counts.values()))
        # if total > 0:
        #     qprobs = {}
        #     for (op,arridx), count in self.counts.items():
        #         if op not in qprobs: qprobs[op] = 0
        #         qprobs[op] += sum(count)
        #     self.opqcount = dict(qprobs)
        #     self.qprobs = dict([(op, v / float(total)) for op, v in qprobs.items()])
        # else:
        #     self.opqcount = {}
        #     self.qprobs = {}
        
            
    def _cache_fanouts(self, qeids, diskc, runc):
            
        def f(w):
            op = w.op
            for arridx in xrange(w.nargs):
                key = (op, arridx)
                # try to get from database
                fstats, bstats = Stats.instance().get_iq_stat(qeids, op.oid, arridx)
                fscales, bscales = [stat[0] for stat in fstats], [stat[0] for stat in bstats]
                fscale = self._wma(fscales)
                bscale = self._wma(bscales)
                self.qfanouts[key] = fscale
                self.qfanins[key] = bscale

                fcounts, bcounts =  [stat[1] for stat in fstats], [stat[1] for stat in bstats]
                fcount = self._wma(fcounts, default=0)
                bcount = self._wma(bcounts, default=0)
                self.counts[key] = [fcount, bcount]

                fsizes, bsizes =  [stat[2] for stat in fstats], [stat[2] for stat in bstats]
                fsize = self._wma(fsizes, default=100000)
                bsize = self._wma(bsizes, default=100000)
                self.fqsizes[key] = fsize
                self.bqsizes[key] = bsize
#                print '_cache\t',op, arridx, fscale, fsize, fcount, bscale, bsize, bcount
                
        self.workflow.visit(f)

        
            

    def _cache_stats(self):
        def f(w):
            op = w.op
            for arridx in xrange(w.nargs):
                if (op,arridx) not in self.cache:
                    stats = Stats.instance().get_op_stats(self.eids, op.oid, arridx)

                    if not stats or len(stats) == 0:
                        stats = [1.0] * 8
                        stats.append(0.00001)
                    if stats[0] == 0:
                        print stats
                        print op
                        print arridx
                        exit()
                    
                    self.cache[(op,arridx)] = stats
        self.workflow.visit(f)
        
    
    def get_input_shape(self, op, arridx):
        return self.cache[(op, arridx)][7]

    def get_output_shape(self, op, arridx):
        return self.cache[(op, arridx)][6]

    def get_provcost(self, op, strat, runid = None):
        opcost, prov, disk, qcosts = self.est_cost(op,strat, runid=runid)
        return prov#opcost + prov

    def get_opcost(self, op, strat, runid = None):
        opcost, prov, disk, qcosts = self.est_cost(op,strat, runid=runid)
        return opcost

    def get_disk(self, op, strat, runid = None):
        opcost, prov, disk, qcosts = self.est_cost(op,strat, runid=runid)
        return disk / 1048576.0

    def get_pqcost(self, op, strat, runid = None):
        # if 'PTMAP_ONE_KEY' in str(strat) and 'CreateM' in str(op) and len(self.qeids) > 1:
        #     import pdb
        #     pdb.set_trace()
        opcost, prov, disk, qcosts = self.est_cost(op,strat,runid=runid)
#        print "\t", op, '\t', '%.8f' % qcosts, '\t', strat
        return qcosts
        

    def est_cost(self, op, strat, runid = None):
        opcosts = []
        provs = []
        disks = []
        qcosts = []

        # weights = []
        # for arridx in xrange(op.wrapper.nargs):
        #     key = (op, arridx)
        #     if self.opqcount.get(op, 1.0) == 0:
        #         weights.append(0.0)
        #     else:
        #         weight = sum(self.counts.get(key,[0])) / float(self.opqcount.get(op,1.0))
        #         weights.append(weight)

        # if sum(weights) < 1.0:
        #     weights = [1.0 / op.wrapper.nargs] * op.wrapper.nargs

        for arridx in xrange(op.wrapper.nargs):
            prov, disk, fcost, bcost, opcost = self.est_arr_cost(op, strat, runid, arridx)            
            key = (op, arridx)
            fprob = self.fprobs.get(key,0)
            bprob = self.bprobs.get(key,0)
            qprob = self.qprobs.get(key, 0)
            qcost = qprob * ((fcost * fprob) + (bcost * bprob))

            modes = strat.modes()
            for mode in modes:
                if mode != Mode.QUERY and mode not in op.supported_modes():
                    qcost = 1000.0

            opcosts.append(opcost)
            provs.append(prov)
            disks.append(disk)
            qcosts.append(qcost)

        opcosts = np.mean(opcosts)
        provs = sum(provs)
        disks = sum(disks)
        qcosts = sum(qcosts)

        if runid == None:
            runid = self.workflow._runid
        prob = self.probs.get(runid, 0.0)
        # if 'KEY' in str(strat) and self.debug:
        #     print 'RUN', op, strat, prob, opcosts, provs, runid, self.workflow._runid
        #     import pdb
        #     pdb.set_trace()
        #     self.debug = False

        return opcosts, provs, disks, qcosts

        if runid == self.workflow._runid:
            return prob * opcosts, prob * provs, disks, prob * qcosts
        return 0, 0, disks, prob * qcosts

            
    def est_arr_cost(self, op, strat, runid, arridx, fqsize=None, bqsize=None):
        key = (op, arridx)
        stats = self.cache[key]
        fanin, area, density, oclustsize, nptrs, noutcells, outputsize, inputsize, opcost = stats
        #print "stats",  op, strat, arridx, weight, stats

        if fqsize is None:
            fqsize = self.fqsizes.get((op, arridx), 1000000)
        if bqsize is None:
            bqsize = self.bqsizes.get((op, arridx), 1000000)
        boxperc = area / inputsize
        prov = write_model(strat, fanin, oclustsize, density, noutcells, opcost) 
        disk = disk_model(strat, fanin, oclustsize, density, noutcells)
        fcost = forward_model(strat, fanin, oclustsize, density, noutcells, opcost, fqsize, 1.0, inputarea=inputsize)
        bcost = backward_model(strat, fanin, oclustsize, density, noutcells, opcost, bqsize, 1.0, inputarea=inputsize)
        # idx,key,parse,extract = forward_model_desc(list(strat.descs())[0], fanin, oclustsize,
        #                                            density, noutcells, opcost, fqsize, 1.0, inputarea=inputsize)

        return prov, disk, fcost, bcost, opcost



    def _proc_fpath(self, ncoords, run_id, path):
        """
        count the number of queries and query size for each operator in the path
        """
        op, arridx = tuple(path[0])
        wop = op.wrapper
        mininputsize = ncoords
        #print "fq", op, arridx, ncoords

        key = (op,arridx)
        if key not in self.counts:
            self.counts[key] = [0,0] # f, b
            self.fqsizes[key] = []
        self.counts[key][0] += 1
        self.fqsizes[key].append(ncoords)        

        path = path[1:]
        if len(path) == 0:  return

        found = False
        for child, idx in wop.children():
            if child == path[0][0] and idx == path[0][1]:
                found = True
        if not found:
            for child, idx in wop.children():
                print child, idx
            raise RuntimeError, ("could not find %s\t%d" % ( path[0][0], path[0][1] ))

        fanout = self.qfanouts.get((op, arridx), 100000)
        mininputsize = ncoords * fanout
        mininputsize = min(mininputsize, self.get_output_shape(op, arridx))
        
        self._proc_fpath(mininputsize, run_id, path)

    def _proc_bpath(self, ncoords, run_id, path):
        op, arridx = path[0]
        wop = op.wrapper
        #print "bq", op, arridx, ncoords

        if op not in self.bqsizes:
            self.bqsizes[(op, arridx)] = []
        self.bqsizes[(op, arridx)].append(ncoords)


        if (op, arridx) not in self.counts:
            self.counts[(op, arridx)] = [0,0]
        self.counts[(op,arridx)][1] += 1


        # multiply by fanin
        fanin = self.qfanins.get((op, arridx), 100000)
        ncoords *= fanin
        ncoords = min(ncoords, self.get_input_shape(op, arridx))


        path = path[1:]
        if len(path) == 0: return 

        # verify the next in the path
        found = False
        for par, idx in wop.parents():
            if par == path[0][0]:
                found = True
        if not found:  raise RuntimeError
        self._proc_bpath(ncoords, run_id, path)



if __name__ == '__main__':

    if len(sys.argv) < 3:
        print "exec_lsst_full.py [opt|box|list|set|grid|test] [log output dir] [run mode: 0, 1]"
        exit()

    runtype, logdir, runmode = sys.argv[1], sys.argv[2], int(sys.argv[3])
    runtypes = ['box', 'bulk', 'set', 'grid', 'test', 'query', 'sql', 'run', 'noop', 'opt']
    if runtype not in runtypes:
        print "%s is not a recognized runtype" % runtype
        exit()
    if os.path.exists(logdir) and not os.path.isdir(logdir):
        print "%s exists and is not a directory"  % logdir
        exit()

    if not os.path.exists(logdir):
        os.mkdir(logdir)


    smallshape = (50,50)    
    Stats.instance('%s/stats.db' % logdir)
    print runmode
    noop_eids = Stats.instance().get_matching_noops(runmode, smallshape)
    print noop_eids

    stats = Stats.instance().get_op_stats(noop_eids, 11, 0)
    print stats
    


    

    for wmodel, dmodel, bmodel, fmodel in zip(write_models, disk_models, back_models, forw_models):
        for eid, opid, op, arridx, fanin, area, density, nptrs, oclustsize, noutcells, opcost in stats:

            wcost = wmodel(fanin, oclustsize, density, noutcells)
            dcost = dmodel(fanin, oclustsize, density, noutcells)
            bcost = fcost = 0,0

            tmp = Stats.instance().get_array_sizes(eid, opid, arridx, 0)
            if tmp != None:
                incount, outcount = tuple(tmp)
                qsize = 10
                boxperc = area / incount
                bcost = bmodel(fanin, oclustsize, density, noutcells, opcost, boxperc, qsize)

            tmp = Stats.instance().get_array_sizes(eid, opid, arridx, 1)
            if tmp != None:
                incount, outcount = tuple(tmp)
                qsize = 10
                boxperc = area / incount
                fcost = fmodel(fanin, oclustsize, density, noutcells, opcost, qsize, incount)
            
            print ('%f\t' * 4) % ( wcost, dcost, bcost, fcost)
    
