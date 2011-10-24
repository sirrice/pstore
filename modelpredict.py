import numpy as np
from models import *
from stats import *


class ModelPredictor(object):

    def __init__(self, eids, workflow, queries):
        self.workflow = workflow
        self.queries = queries
        self.counts = {}
        self.cache = {}
        self.fprobs = {}
        self.qprobs = {}
        self.fqsizes = {}
        self.bqsizes = {}
        self.eids = eids

        self._cache_stats()
        self._proc_queries()

        self.fqsizes = dict([(k,np.mean(v)) for k,v in self.fqsizes.items()])
        self.bqsizes = dict([(k,np.mean(v)) for k,v in self.bqsizes.items()])



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

    def get_provcost(self, op, strat):
        opcost, prov, disk, qcosts = self.est_cost(op,strat)
        return prov#opcost + prov

    def get_opcost(self, op, strat):
        opcost, prov, disk, qcosts = self.est_cost(op,strat)
        return opcost

    def get_disk(self, op, strat):
        opcost, prov, disk, qcosts = self.est_cost(op,strat)
        return disk / 1048576.0

    def get_pqcost(self, op, strat):
        opcost, prov, disk, qcosts = self.est_cost(op,strat)
        return qcosts
        

    def est_cost(self, op, strat):
        opcosts = []
        provs = []
        disks = []
        qcosts = []

        weights = []
        for arridx in xrange(op.wrapper.nargs):
            key = (op, arridx)
            weight = sum(self.counts.get(key,[0])) / float(self.opqcount.get(op,1.0))
            weights.append(weight)
        if sum(weights) < 1.0:
            weights = [1.0 / op.wrapper.nargs] * op.wrapper.nargs
        
        for arridx, weight in zip(xrange(op.wrapper.nargs), weights):
            key = (op, arridx)
            fprob = self.fprobs.get(key,0)

            stats = self.cache[key]

            #print "stats",  op, strat, arridx, weight, stats
            
            fanin, area, density, oclustsize, nptrs, noutcells, outputsize, inputsize, opcost = stats

            fqsize = self.fqsizes.get(op, 1)
            bqsize = self.bqsizes.get(op, 1)
            boxperc = area / inputsize

            prov = write_model(strat, fanin, oclustsize, density, noutcells) 
            disk = disk_model(strat, fanin, oclustsize, density, noutcells)
            fcost = forward_model(strat, fanin, oclustsize, density, noutcells, opcost, fqsize, 1.0, inputarea=inputsize)
            bcost = backward_model(strat, fanin, oclustsize, density, noutcells, opcost, bqsize, 1.0, inputarea=inputsize)
            #print fcost, bcost, fprob, op, strat
            qcost = (fcost * fprob) + (bcost * (1.0 - fprob))
            #print qcost, op, strat


            modes = strat.modes()
            for mode in modes:
                if mode != Mode.QUERY and mode not in op.supported_modes():
                    qcost = 10000.0

            opcosts.append(opcost)
            provs.append(prov)
            disks.append(disk)
            qcosts.append(weight*qcost)
        return sum(opcosts), sum(provs), sum(disks), sum(qcosts)
        return np.mean(opcosts), np.mean(provs), np.mean(disks), np.mean(qcosts)


    def _proc_queries(self):
        for query in self.queries:
            inputs, run_id, path, direction = query
            if direction == "forward":
                self._proc_fpath(inputs, run_id, path)
            else:
                self._proc_bpath(inputs, run_id, path)

        # probability of forward query
        self.fprobs = dict([(key, val[0] / float(sum(val))) for key, val in self.counts.items()])

        # probability of querying an op
        total = sum(map(sum, self.counts.values()))
        qprobs = {}
        for (op,arridx), count in self.counts.items():
            if op not in qprobs: qprobs[op] = 0
            qprobs[op] += sum(count)
        self.opqcount = dict(qprobs)
        self.qprobs = dict([(op, v / float(total)) for op, v in qprobs.items()])

        

    def _proc_fpath(self, ncoords, run_id, path):
        """
        count the number of queries and query size for each operator in the path
        """
        op, arridx = tuple(path[0])
        wop = op.wrapper
        mininputsize = ncoords


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
        
        self._proc_fpath(mininputsize, run_id, path)

    def _proc_bpath(self, ncoords, run_id, path):
        op, arridx = path[0]
        wop = op.wrapper

        if op not in self.bqsizes:
            self.bqsizes[op] = []
        self.bqsizes[op].append(ncoords)
        

        if (op, arridx) not in self.counts:
            self.counts[(op, arridx)] = [0,0]
        self.counts[(op,arridx)][1] += 1

        # multiply by fanin
        fanin = self.cache.get((op,arridx), (1,1))[0]
        ncoords *= fanin

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
    
