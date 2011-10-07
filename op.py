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



diskseek = 0.0005#0.005
diskread = 0.000000017339533 * 1048576



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
        wlog.info("%s\taddinput\t%d", self.op, idx)
        self.slots[idx] = arrid

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
        wpstore = self.op.pstore(run_id)
        
        output, stats = self.op.run(newinputs, run_id)
        newinputs = None

        # calculate runtime and provenance overheads
        wpstore.close()

        # store outputs
        outputid = ArrayStore.instance().add(output)
        outputshape = output.shape
        del output
        self.outputs[run_id] = outputid


        #Stats.instance().add_wrun(run_id, self.op, runtime, inputshapes, outputshape,
        #                           wpstore, rpstore)
        return outputid

    def provsize(self, s, run_id=None):
        if s.s in (Strategy.FUNCTION, Strategy.QUERY):
            return 1.0 / 1048576.0
        return Stats.instance().get_disk(self.op, s, run_id)

    def provoverhead(self, s, run_id=None):
        """
        overhead of capturing provenance during runtime
        """
        if s.s in (Strategy.FUNCTION, Strategy.QUERY):
            return 0.000001
            return 0.002 # empirical from FUNCTION operators
        return Stats.instance().get_overhead(self.op, s, run_id)
        
        

    def cost(self, strategy, run_id=None):
        if strategy.s == Strategy.FUNCTION:
            if not self.op.implements_mapfunctions():
                return 1000
            return 0.0
        
        stats = Stats.instance()
        cost =  stats.get_provq_cost(self.op, strategy, run_id)
        if cost is not None: return cost

        if strategy.s == Strategy.QUERY:
            avg_runtime = stats.get_provq_cost(self.op, STRAT_PSET)
            return stats.get_runtime(self.op, strategy) + avg_runtime #self.prov_runtime(strategy)
        else:
            avg_runtime = stats.get_provq_cost(self.op, strategy)
            if avg_runtime: return avg_runtime
            return diskseek + diskread * self.provsize(strategy, run_id) + 0.00001
        
    

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

    def default_strategy(self):
        # decide what provenance strategy to use for the nodes
        # rules based : if fmap, then use FuncPStore, else PtrPStore (which one?)
        runtime = Runtime.instance()
        def set_strategy(w):
            if Mode.FULL_MAPFUNC in w.op.supported_modes():
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

    def forward_path(self, incoords, run_id, path):
        """
        Compile query into a query plan.
        Return: query plan
        """
        start = time.time()
        child = Scan(incoords)
        for op, arridx in path:
            wop = op.wrapper
            wop = op.wrapper
            shape = wop.get_output_shape(run_id)
            pstore = Runtime.instance().get_pstore(op, run_id) 
            if pstore is None: raise RuntimeError
            child = pstore.join(child, arridx, backward=False)

        end = time.time()
        return DedupQuery(child, shape)
        return NBDedupQuery(q)



    def backward_path(self, outcoords, run_id, path):
        """
        return: query plan
        """
        child = Scan(outcoords)
        for op, arridx in path:
            wop = op.wrapper
            shape = wop.get_input_shape(run_id, arridx)
            pstore = Runtime.instance().get_pstore(op, run_id)
            if pstore is None: raise RuntimeError

            child = pstore.join(child, arridx, backward=True)

        return DedupQuery(child, shape)

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
            if not w.op.implements_mapfunctions():
                ops.add(w.op)
        self.visit(collect)
        return list(ops)

        
    def get_matstrats(self):
        return [STRAT_F, STRAT_Q,  STRAT_PSET, STRAT_BOX, STRAT_BULK, STRAT_DIFF]
        return [STRAT_F, STRAT_Q,  STRAT_PSET, STRAT_PGRID, STRAT_BOX, STRAT_SUPERBOX, STRAT_BULK]

        
    def get_op_costs(self, op, matstrat):
        total = 0.0
        nstats = 0
        for stat in self.indivq_stats.get(op, []):
            total += op.wrapper.cost(stat, matstrat)
            nstats += 1
        if nstats == 0:
            #wlog.error('opcosts %s\t%s\t%f\t%s\t%s\tdefault', op, matstrat, total , 0, 5.0)
            return 5.0
        #wlog.error('opcosts %s\t%s\t%f\t%s\t%s', op, matstrat, total , nstats, total / nstats)
        if nstats:
            return total / nstats
        return 0.0

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
        return ('%s[%d]' % (type(self).__name__, self.oid)).ljust(20)

    def __eq__(self, op):
        return hash(op) == hash(self)



