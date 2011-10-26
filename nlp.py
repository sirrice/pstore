import os, logging, numpy
from util import print_matrix
from runtime import *
from cvxopt.glpk import ilp
from cvxopt import *


nlog = logging.getLogger('nlp')
logging.basicConfig()
nlog.setLevel(logging.INFO)



def run_nlp(stats, w, mp, maxdisk, maxoverhead):
    """
    maxoverhead is percentage of average runtime cost of optimizable operators
    """
    ops = w.get_optimizable_ops()
    matstrats = w.get_matstrats()
    pairs = [(op, s) for op in ops for s in matstrats]

    avg_runtime = sum(filter(lambda x: x > 0, [mp.get_opcost(op, s) for op,s in pairs]))
    maxoverhead *= avg_runtime

    G1 = [mp.get_disk(op,s) for (op,s) in pairs]
    G2 = [mp.get_provcost(op, s) for (op,s) in pairs]
    G = matrix([G1, G2]).trans()
    h = matrix([maxdisk, maxoverhead])

    A = []
    for i in xrange(len(ops)):
        row = []
        for col in xrange(len(pairs)):
            if len(matstrats) * i <= col and col < len(matstrats) * (i+1):
                row.append(1.0)
            else:
                row.append(0.0)
        A.append(row)
    A = matrix(A).trans()
    b = matrix([1.0] * len(ops))

    c = [mp.get_pqcost(op,s) for op,s in pairs]
    #c = map(lambda cost: cost * 100.0, c)
    minc = min(cost for cost in c) / 2.0
    # normalize disk and runcost to minc
    G1p = [g / max(G1) * minc for g in G1]
    G2p = [g / max(G2) * minc for g in G2]
    c = map(sum, zip(c, G1p, G2p))

    d = dict([(p, cost) for cost, p in zip(c, pairs)])
    c = matrix(c)


    nlog.debug("Constraints: %f\t%f" , maxdisk, maxoverhead)
    for op, s in pairs:
        pqcost = mp.get_pqcost(op, s)
        nlog.debug('%s\t%s\t%f\t%f\t%f\t%f', op, str(s).ljust(25),
                   pqcost,
                   mp.get_disk(op,s),
                   mp.get_provcost(op, s),
                   d[(op, s)])
        

    strategies = nlp_exec_cvx(c, ops, matstrats, G, h, A, b)

    return strategies
    

def nlp_exec_cvx(c, ops, matstrats, G, h, A, b):
    solvers.options['show_progress'] = False
    I = set()
    B = set(range(len(c)))

    status, x = ilp(c, G, h, A, b, I, B)
    nlog.info( "resources  \t%s", (G*x).trans() )
    nlog.info( "constraints\t%s", h.trans() )
    nlog.info( "cost       \t%s", x.trans() * c )

    x = numpy.array(x).reshape((len(ops), len(matstrats))).tolist()
    strategies = assign_strategies(x, ops, matstrats)
    return strategies
    

def assign_strategies(matrix, ops, strats):
    ret = {}
    for x, op in enumerate(ops):
        ret[op] = []
        offset = x * len(strats)
        for y, s in enumerate(strats):
            if matrix[x][y]:
                ret[op].append(s)
    return ret
