import sqlite3, time, logging, sys, numpy, os, operator, random
from operator import mul
from provstore import *
from arraystore import *
from runtime import *
from util import *
from nlp import *
from queryresult import *
from stats import Stats

wlog = logging.getLogger('workflow')
logging.basicConfig()
wlog.setLevel(logging.ERROR)


class Wrapper(object):
    """
    Manages metadata and execution on a per-operator basis
      configuration files and external data are modeled as wrappers without parents
      that return the data
    """
    def __init__(self, op, nargs):
        self.op = op
        op.wrapper = self
        self.nargs = nargs

        self.slots = [None] * nargs
        self.inputs = {}
        self.outputs = {}

    def serialize(self):
        arr = [self.nargs, self.inputs, self.outputs]
        return pickle.dumps(arr)

    def parse(self, s):
        arr = pickle.loads(s)
        self.nargs, self.inputs, self.outputs = tuple(arr)
        

    def log_stats(self, stats):
        self.stats.append(stats)

    def clean_inputs(self):
        self.slots = [None] * self.nargs

    def add_input(self, idx, arrid):
        self.slots[idx] = arrid

    def get_input(self, run_id, arridx):
        if run_id not in self.inputs:
            raise Exception("run #%d does not exist" % run_id)
        if arridx >= len(self.inputs[run_id]):
            raise Exception("..")
        return ArrayStore.instance().get(self.inputs[run_id][arridx])

    def get_inputs(self, run_id):
        if run_id not in self.inputs:
            raise Exception("run #%d does not exist" % run_id)
        return map(lambda x: ArrayStore.instance().get(x), self.inputs[run_id])

    def get_input_shapes(self, run_id):
        if run_id not in self.inputs:
            raise Exception, "run #%d does not exist" % run_id
        ret = []
        for idx in xrange(self.nargs):
            ret.append(ArrayStore.instance().shape(self.inputs[run_id][idx]))
        return ret
        

    def get_input_shape(self, run_id, arridx):
        if run_id not in self.inputs:
            raise Exception, "run #%d does not exist" % run_id
        return ArrayStore.instance().shape(self.inputs[run_id][arridx])


    def get_output(self, run_id):
        if run_id not in self.outputs:
            raise Exception, "run #%d does not exist" % run_id
        return ArrayStore.instance().get(self.outputs[run_id])

    def get_output_shape(self, run_id):
        if run_id not in self.outputs:
            raise Exception, "run #%d does not exist" % run_id
        return ArrayStore.instance().shape(self.outputs[run_id])

    def parents(self):
        return self.workflow.parents(self.op)

    def children(self):
        return self.workflow.children(self.op)

    def ready(self):
        return None not in self.slots

    

    def run(self, run_id):
        # queue input into the write input slot,
        # execute oprator if all inputs are available
        if not self.ready(): return

        newinputs = map(ArrayStore.instance().get, self.slots)
        inputshapes = [arr.shape for arr in newinputs]
        self.inputs[run_id] = list(self.slots)
        self.slots = [None] * self.nargs

        wlog.info('%s.run(%d)', self.op, run_id)
        pstore = self.op.pstore(run_id)


        start = time.time()        
        output, stats = self.op.run(newinputs, run_id)
        runtime = time.time() - start

        newinputs = None

        # calculate runtime and provenance overheads
        pstore.close()

        if pstore.get_stat('write') > 0:
            for f in ('outcache', 'incache', 'serin', 'serout', 'mergcost', 'bdbcost',
                      'keycost', 'idxcost', 'write', 'flush'):
                wlog.debug( "%s\t%f", f, pstore.get_stat(f) )

        # store outputs
        outputid = ArrayStore.instance().add(output)
        outputshape = output.shape
        del output
        self.outputs[run_id] = outputid
        outputdisk = ArrayStore.instance().size(outputid)

        Stats.instance().add_wrun(run_id, self.op, runtime, inputshapes, outputshape, outputdisk, pstore)
        return outputid


    def __str__(self):
        return str(self.op)


class Workflow(object):
    """
    Manage DAG of wrappers, and metadata between wrappers
    """
    _runid = 1
    
    def __init__(self):
        self.ops = {}         # op -> wrapper
        self.id2op = {}       # oid -> op
        self.par2child = {}   # [par][child] -> list(argidx)
        self.child2par = {}   # [child] -> list(par)

        # flag for if we want to optimize query execution
        # set to True and create a model prediction object to optimize
        self.boptimize = False
        self.mp = None
        self.beststrats = None


    def serialize(self):
        arr = []
        arr.append(self._runid)
        arr.append(map(lambda op: op.oid, self.ops.keys()))
        arr.append(map(lambda op: op.wrapper.serialize(), self.ops.keys()))
        arr.append(Runtime.instance().serialize())
        return pickle.dumps(arr)

    def parse(self, s):
        arr = pickle.loads(s)
        self._runid = arr[0]
        for oid, serstr in zip(arr[1], arr[2]):
            self.id2op[oid].wrapper.parse(serstr)
        Runtime.instance().parse(arr[3], self)
        
    def save(self, fname):
        f = file(fname, 'w')
        s = self.serialize()
        f.write(s)
        f.close()

    def load(self, fname):
        f = file(fname, 'r')
        self.parse(f.read())
        f.close()
        
            


            

    def register(self, op, nargs):
        wrapper = Wrapper(op, nargs)
        self.ops[op] = wrapper
        self.id2op[op.oid] = op
        op.workflow = self
        wrapper.workflow = self

    def nargs(self, op):
        return self.ops[op].nargs

    def wrapper(self, op):
        return self.ops.get(op, None)

    def clear_stats(self):
        pass


    def connect(self, par, child, argidx):
        "Add parent-child relationship between two operators in the workflow"
        parents = self.child2par.get(child, {})
        argidxs = parents.get(par,[])
        if argidx not in argidxs: argidxs.append(argidx)
        parents[par] = argidxs
        self.child2par[child] = parents

        children = self.par2child.get(par, {})
        argidxs = children.get(child, [])
        if argidx not in argidxs: argidxs.append(argidx)
        children[child] = argidxs
        self.par2child[par] = children


    def children(self, op):
        "returns list of child,arg pairs"
        ret = []
        children = self.par2child.get(op, {})
        for child, args in children.items():
            for arg in args:
                ret.append((child, arg))
                #yield (child, arg)
        return ret

    def parents(self, op):
        "return list of parent operators, input argument idx pairs"
        ret = []
        parents = self.child2par.get(op, {})
        for par, args in parents.items():
            for arg in args:
                ret.append((par, arg))
                #yield (par, arg)
        return ret


    def roots(self):
        "return list of (root, [indexes that need inputs]) pairs"
        ret = []
        for op, opw in self.ops.items():
            inputs = [0] * opw.nargs
            needinputs = [] # indexes that need inputs
            for par, arg in self.parents(op):
                inputs[arg] = 1
            for idx, v in enumerate(inputs):
                if not v: needinputs.append(idx)
            if needinputs:
                ret.append((op, needinputs))
        return ret

    def nodes(self):
        "return list of wrappers of workflow's operators"
        nodes = set()
        def f(w):
            nodes.add(w)
        self.visit(f)
        return list(nodes)

    def add_static_input(self, op, arg, data):
        """
        Add an input to an open slot of a root operator
        open: not connected to another operator's output
        """
        # check if this is allowed
        roots = self.roots()
        root, args = None, None
        for pair in roots:
            if op == pair[0]:
                root, args = pair
                break

        if not root:
            raise Exception("%s is not a root node and cannot accept static arguments" % op)

        if arg not in args:
            raise Exception("%s arg #%d is not a  static argument" % (op, arg))

        arrid = ArrayStore.instance().add(data)
        self.ops[op].add_input(arg, arrid)
        
    
    def prepare(self):
        """
        prepares the workflow for a new execution
         - goes through the operators and sets their inputs slots to empty (or default value)
        """
        for wop in self.ops.values():
            wop.clean_inputs()

    def default_strategy(self, userstrat = None):
        # decide what provenance strategy to use for the nodes
        # rules based : if fmap, then use FuncPStore, else PtrPStore (which one?)
        runtime = Runtime.instance()
        def set_strategy(w):
            if userstrat != None:
                strat = userstrat
            elif Mode.FULL_MAPFUNC in w.op.supported_modes():
                strat = Strat.single(Mode.FULL_MAPFUNC, Spec(Spec.NONE, Spec.NONE))
            else:
                strat = Strat.single(Mode.QUERY, Spec(Spec.NONE, Spec.NONE))

            runtime.set_strategy(w.op, strat)
        self.visit(set_strategy)


    def run(self):
        """
        recursively execute workflow from roots
        roots: [(op, inputs)*]
        """
        run_id = self._runid
        self._runid += 1

        wlog.info("workflow #%d" , run_id)
       

        jobq = []
        for root, static_args in self.roots():
            w = self.ops[root]
            for static_arg in static_args:
                if w.slots[static_arg] is None:
                    raise Exception("%s not all inputs provided" % root)
            jobq.append(w)

        
        wlog.info("running jobq")

        while len(jobq) > 0:
            job = jobq.pop(0)

            if job.ready():

                outputid = job.run(run_id)

                # check if output is used by any operators
                wlog.info("piping to children %s", self.children(job.op))
                for childop, argidx in self.children(job.op):
                    childjob = self.ops[childop]
                    childjob.add_input(argidx, outputid)
                    if childjob not in jobq:
                        jobq.append(childjob)
            else:
                jobq.append(job)
                #wlog.info("%s is not ready\t%s", job.op,
                #          map(lambda x: x is not None and 1 or 0, job.slots))


        wlog.info("workflow #%d done!", run_id)


    def pick_forward_strat(self, qsize, op, arridx, run_id):
        fanin, area, density, oclustsize, nptrs, noutcells, outputsize, inputsize, opcost = self.mp.cache[(op,arridx)]
        wlog.debug( 'PrepFQ\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%f', str(op), fanin, oclustsize, nptrs, noutcells, inputsize, outputsize, density, opcost  )
        beststrat = None
        bestcost = None
        for strat in Runtime.instance().available_strats(op, run_id):
            _,_,fcost,_,opcost = self.mp.est_arr_cost(op, strat, run_id, arridx, fqsize=qsize)
            wlog.debug( 'PrepFQ\t%s\t%d\t%s\t%d\t%f', op, arridx, strat, qsize, fcost)
            if beststrat == None or fcost < bestcost:
                beststrat = strat
                bestcost = fcost
        return beststrat


    def forward_path(self, incoords, run_id, path):
        """
        Compile query into a query plan.
        Return: query plan
        """
        optcost = 0.0
        qsizes = {}
        child = Scan(incoords)
        for op, arridx in path:
           
            wop = op.wrapper
            shape = wop.get_output_shape(run_id)

            if self.boptimize:
                optstart = time.time()
                strat = self.pick_forward_strat(len(child), op, arridx, run_id)
                optcost += time.time() - optstart
                pstore = Runtime.instance().get_query_pstore(op, run_id, strat)
            else:
                strat = Runtime.instance().get_strategy(op, run_id)
                pstore = Runtime.instance().get_pstore(op, run_id)
            if pstore is None: raise RuntimeError

            insize = len(child)
            start = time.time()
            child = pstore.join(child, arridx, backward=False)
            child = DedupQuery(child, shape)
            end = time.time()
            outsize = len(child)
            wlog.debug( 'Fpathdedup\t%f\t%s\t%d\t%s', end-start, op, len(child), str(pstore.strat) )
            qsizes[(op, arridx)] = (end-start, strat, insize, outsize)

        end = time.time()
        return DedupQuery(child, shape), optcost, qsizes
        return NBDedupQuery(q)

    def pick_backward_strat(self, qsize, op, arridx, run_id):
        beststrat = None
        bestcost = None
        for strat in Runtime.instance().available_strats(op, run_id):
            _,_,_,bcost,opcost = self.mp.est_arr_cost(op, strat, run_id, arridx, bqsize=qsize)
            wlog.debug( 'PrepBQ\t%s\t%d\t%s\t%d\t%f', op, arridx, strat, qsize, bcost)
            if beststrat == None or bcost < bestcost:
                beststrat = strat
                bestcost = bcost
        return beststrat
        


    def backward_path(self, outcoords, run_id, path):
        """
        return: query plan
        """
        wlog.debug( 'New Query\t%s', '\t'.join(map(lambda s: s.strip(), map(str, map(lambda p: p[0], path)))))
        optcost = 0.0
        qsizes = {}
        child = Scan(outcoords)
        for op, arridx in path:
            wop = op.wrapper
            shape = wop.get_input_shape(run_id, arridx)

            if self.boptimize:
                optstart = time.time()
                strat = self.pick_backward_strat(len(child), op, arridx, run_id)
                optcost += time.time() - optstart
                pstore = Runtime.instance().get_query_pstore(op, run_id, strat)
            else:
                strat = Runtime.instance().get_strategy(op, run_id)
                pstore = Runtime.instance().get_pstore(op, run_id)
            if pstore is None: raise RuntimeError

            insize = len(child)
            start = time.time()
            child = pstore.join(child, arridx, backward=True)
            child = DedupQuery(child, shape)
            end = time.time()
            outsize = len(child)
            wlog.debug( 'Bpathdedup\t%f\t%s\t%d\t%s', end-start, op, len(child), str(pstore.strat) )
            qsizes[(op, arridx)] = (end-start, strat, insize, outsize)
            
        return DedupQuery(child, shape), optcost, qsizes

    def rewrite(self, root):
        if isinstance(root, Scan): return root
        
        child = self.rewrite(root.child)
        if hasattr(child, "pstore") and child.pstore.op.alltoall(child.arridx):
            return child
        return root
        
 

    def split_forward_pq(self, op, nsteps, f=lambda x:x):
        children = self.children(op)
        if nsteps <= 0 or not children:
            return
        f(op)
        for child, arg in children:
            self.split_forward_pq(child, nsteps-1)


    def split_backward_pq(self, op, nsteps, f=lambda x:x):
        parents = self.parents(op)
        if nsteps <= 0 or not parents:
            return
        f(op)
        for par, arg in parents:
            self.split_backward_pq(par, nsteps-1)


    def gen_nlp(self, maxdisk=500, maxoverhead=0.5):
        strategies = run_nlp(Stats.instance(), self, maxdisk, maxoverhead)
        for op, strats in strategies.items():
            wlog.info("%s\t%s", op, ', '.join(map(str, strats)))
        wlog.info("\n")
        return strategies


    def get_optimizable_ops(self):
        ops = set()
        def collect(w):
            #if Mode.FULL_MAPFUNC not in w.op.supported_modes():
            #if  'CreateM' in str(w.op):
            # ops.add(w.op)
            ops.add(w.op)
        self.visit(collect)
        return list(ops)

        
    def get_matstrats(self):
        strats = [Strat.single(Mode.QUERY, Spec.default()),
                  Strat.single(Mode.FULL_MAPFUNC, Spec.default())]

        bstrats = [Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.KEY), True),
                   Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.KEY), True),
                   Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
                   Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True)]

        fstrats = [Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
                   Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False)]

        cstrats = []
        for bstrat in bstrats:
            for fstrat in fstrats:
                bstrat = bstrat.copy()
                fstrat = fstrat.copy()
                cstrat = Strat([bstrat.buckets[0], fstrat.buckets[0]])
                cstrats.append(cstrat)

        strats.extend(bstrats)
        strats.extend(fstrats)
        strats.extend(cstrats)                
        return strats


    def get_op_prob(self, op, iqs):
        if len(iqs):
            return float(len(filter(lambda iq: iq.op == op, iqs))) / len(iqs)
        return 0.1
                
        
    def get_indiv_qs(self):
        uniques = {}
        iqs = []
        for op, stats in self.indivq_stats.items():
            for stat in stats:
                key = (op, stat.direction)
                if key not in uniques:
                    uniques[key] = 0
                uniques[key] += 1
            iqs.extend(stats)
        return iqs, uniques
                

            
    def get_query_probs(self, q, uniqueqs, op):
        if q[0] != op or q not in uniqueqs: return 0.0
        return float(uniqueqs[q]) / sum(uniqueqs.values())
        

    def visit(self, f=lambda x: x):
        """
        visits each note in the workflow and executes f(wrapper)
        """
        queue = []
        queue.extend(map(lambda x: x[0], self.roots()))
        seen = set(queue)
        while len(queue):
            op = queue.pop()
            seen.add(op)
            f(self.wrapper(op))
            for x in self.children(op):
                x = x[0]
                if x not in seen:
                    queue.append(x)




class Op(object):
    _id = 0
    def __init__(self):
        self.deterministic = True
        self.oid = Op._id
        Op._id += 1

        # this is populated when operator is registered to a workflow
        self.workflow = None

        # this is populated when the operator is wrapped in a Wrapper instance
        # it's safe to assume this is set
        self.wrapper = None

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC_BOX]

    def output_shape(self, run_id):
        pass

    def can_box(self):
        return False

    def alltoall(self, arridx):
        return False

    def fmap(self, coord, run_id, arridx):
        pass

    def bmap(self, coord, run_id, arridx):
        pass

    def bbox(self, coord, run_id, arridx):
        return ((0,0), self.wrapper.get_input_shape(run_id, arridx))

    def bmap_obj(self, obj, run_id, arridx):
        pass

    def bbox_obj(self, obj, run_id, arridx):
        return ((0,0), self.wrapper.get_input_shape(run_id, arridx))

    def run(self, inputs, run_id):
        pass

    def pstore(self, run_id):
        return Runtime.instance().get_pstore(self, run_id)

    def __str__(self):
        return ('%s[%d]' % (type(self).__name__, self.oid)).ljust(21)

    def __eq__(self, op):
        return hash(op) == hash(self)



