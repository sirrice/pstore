import numpy as np
from models import *
from stats import *
from operator import mul
from util import zipf


class ModelPredictor(object):

    def __init__(self, eids, workflow, queries, l = 1.9):
        self.workflow = workflow
        self.queries = queries
        self.counts = {}
        self.cache = {}
        self.fprobs = {}
        self.bprobs = {}
        self.qprobs = {}
        self.fqsizes = {}
        self.bqsizes = {}
        self.eids = eids

        self._cache_stats()
        self._proc_queries()
        self._cache_qsizes()


        cumprobs = zipf(workflow._runid, l)
        probs = []
        prob = 0.0
        for p in cumprobs:
            probs.append(p-prob)
            prob = p
        self.probs = dict([(workflow._runid-i, p) for i, p in enumerate(probs)])

        self.debug = True

    def _cache_qsizes(self):
        try:
            qeids = Stats.instance().get_similar_eids(self.eids[0])
        except Exception, e:
            print e
            self.fqsizes = dict([(k,np.mean(v)) for k,v in self.fqsizes.items()])
            self.bqsizes = dict([(k,np.mean(v)) for k,v in self.bqsizes.items()])
            return
            
        def f(w):
            op = w.op
            for arridx in xrange(w.nargs):
                key = (op, arridx)
                # try to get from database
                fqsize, bqsize = Stats.instance().get_iq_stat(qeids, op.oid, arridx)

                if fqsize is not None:
                    self.fqsizes[key] = fqsize
                elif key in self.fqsizes:
                    self.fqsizes[key] = np.mean(self.fqsizes[key])
                else:
                    self.fqsizes[key] = 1

                if bqsize is not None:
                    self.bqsizes[op] = bqsize
                elif op in self.bqsizes:
                    self.bqsizes[op] = np.mean(self.bqsizes[op])
                else:
                    self.bqsizes[op] = 1
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
        opcost, prov, disk, qcosts = self.est_cost(op,strat,runid=runid)
#        print "\t", op, '\t', '%.8f' % qcosts, '\t', strat
        return qcosts
        

    def est_cost(self, op, strat, runid = None):
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
            prov, disk, fcost, bcost, opcost = self.est_arr_cost(op, strat, runid, arridx)            
            key = (op, arridx)
            fprob = self.fprobs.get(key,0)
            bprob = self.bprobs.get(key,0)
            qprob = self.qprobs.get(op, 0)
            qcost = qprob * ((fcost * fprob) + (bcost * bprob))

            modes = strat.modes()
            for mode in modes:
                if mode != Mode.QUERY and mode not in op.supported_modes():
                    qcost = 10.0

            opcosts.append(opcost)
            provs.append(prov)
            disks.append(disk)
            qcosts.append(weight*qcost)

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

        if runid == self.workflow._runid:
            return prob * opcosts, prob * provs, disks, prob * qcosts
        return 0, 0, disks, prob * qcosts

            
    def est_arr_cost(self, op, strat, runid, arridx, fqsize=None, bqsize=None):
        key = (op, arridx)
        stats = self.cache[key]
        fanin, area, density, oclustsize, nptrs, noutcells, outputsize, inputsize, opcost = stats
        #print "stats",  op, strat, arridx, weight, stats

        if fqsize is None:
            fqsize = self.fqsizes.get(op, 1)
        if bqsize is None:
            bqsize = self.bqsizes.get(op, 1)
        boxperc = area / inputsize
        prov = write_model(strat, fanin, oclustsize, density, noutcells, opcost) 
        disk = disk_model(strat, fanin, oclustsize, density, noutcells)
        fcost = forward_model(strat, fanin, oclustsize, density, noutcells, opcost, fqsize, 1.0, inputarea=inputsize)
        bcost = backward_model(strat, fanin, oclustsize, density, noutcells, opcost, bqsize, 1.0, inputarea=inputsize)
        # idx,key,parse,extract = forward_model_desc(list(strat.descs())[0], fanin, oclustsize,
        #                                            density, noutcells, opcost, fqsize, 1.0, inputarea=inputsize)

        return prov, disk, fcost, bcost, opcost

    def _proc_queries(self):
        for query in self.queries:
            inputs, run_id, path, direction = query
            if direction == "forward":
                self._proc_fpath(inputs, run_id, path)
            else:
                self._proc_bpath(inputs, run_id, path)

        # probability of forward query
        self.fprobs = dict([(key, val[0] / float(sum(val))) for key, val in self.counts.items()])
        self.bprobs = dict([(key, val[1] / float(sum(val))) for key, val in self.counts.items()])


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

        fanout = self.cache.get((op, arridx), (1,1))[3]
        mininputsize = ncoords * fanout
        mininputsize = min(mininputsize, self.get_output_shape(op, arridx))
        
        self._proc_fpath(mininputsize, run_id, path)

    def _proc_bpath(self, ncoords, run_id, path):
        op, arridx = path[0]
        wop = op.wrapper
        #print "bq", op, arridx, ncoords

        if op not in self.bqsizes:
            self.bqsizes[op] = []
        self.bqsizes[op].append(ncoords)


        if (op, arridx) not in self.counts:
            self.counts[(op, arridx)] = [0,0]
        self.counts[(op,arridx)][1] += 1


        # multiply by fanin
        fanin = self.cache.get((op,arridx), (1,1))[0]
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
    
