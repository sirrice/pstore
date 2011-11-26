import os, logging, numpy
from util import print_matrix
from runtime import *
from stats import Stats
from cvxopt.glpk import ilp
from cvxopt import *
from cvxopt import glpk


nlog = logging.getLogger('nlp')
logging.basicConfig()
nlog.setLevel(logging.DEBUG)



def run_nlp(stats, w, mp, maxdisk, maxoverhead):
    """
    maxoverhead is percentage of average runtime cost of optimizable operators
    """
    stats = Stats.instance()
    ops = w.get_optimizable_ops()
    matstrats = w.get_matstrats()
    currun = w._runid
    pairs = [(currun, op, s) for op in ops for s in matstrats]

    trips = list(pairs)
    existingops = []
    for r,o,s in Runtime.instance().get_disk_strategies():
        trips.append((r,o,s))
        trips.append((r,o,Strat.query()))
        existingops.append((o,r))

    # for o in ops:
    #     for s in matstrats:
    #         if 'ONE_KEY' in str(s) and 'GetNames' in str(o):
    #             import pdb
    #             pdb.set_trace()
    #             mp.get_pqcost(o,s,currun)

    
    #trips = [(r,op,s) for r in xrange(1, currun+1) for op in ops for s in matstrats]

    xold = [Runtime.instance().check_strategy(op,r,s) and 1 or 0 for r,op,s in trips]

    avg_runtime = sum(filter(lambda x: x > 0, [mp.get_opcost(op, s) for r,op,s in pairs if r == currun]))
    maxoverhead *= avg_runtime

    G1 = []
    G1dict = {}
    for r,op,s in trips:
        if op not in G1dict:
            G1dict[op] = []
        if r == currun:
            disk = mp.get_disk(op, s)
        else:
            disk = stats.get_disk(r,op,s)
        G1.append(disk)
        if disk > 0:
            G1dict[op].append(disk)
    G2 = []
    G2dict = {}
    for r,op,s in trips:
        if op not in G2dict:
            G2dict[op] = []
        if r == currun:
            ov = mp.get_provcost(op, s)
        else:
            ov = 0.0
        G2.append(ov)
        if ov > 0:
            G2dict[op].append(ov)


    
    G = matrix([G1, G2]).trans()
    h = matrix([maxdisk, maxoverhead])

    A = []
    # every operator needs exactly 1 strategy constraint
    blocksize = len(matstrats) 
    for i in xrange(len(ops)):
        row = [0.0] * len(trips)
        for col in xrange(len(trips)):
            if blocksize * i <= col and col < blocksize * (i+1) and col < len(pairs):
                row[col] = 1.0
        A.append(row)
    for i, (op,r) in enumerate(existingops):
        row = [0.0] * len(trips)
        col = len(pairs) + i * 2
        row[col] = 1.0
        row[col+1] = 1.0
        A.append(row)
    A = matrix(A).trans()
    b = matrix([1.0] * (len(ops)+len(existingops)))

    c = []
    mincs = {}
    for r,op,s in trips:
        cost = mp.get_pqcost(op,s,r)
        c.append(cost)
        if cost > 0 and op not in mincs:
            mincs[op] = cost
        if cost > 0 and cost < mincs[op]:
            mincs[op] = cost
        if cost > 1000 and 'ONE_MANY_f' in str(s) and r == w._runid-1 and 'Extract' in str(op):
            import pdb
            pdb.set_trace()
            mp.get_pqcost(op,s,r)
    #c = [mp.get_pqcost(op,s,r) for r,op,s in trips]

    # normalize disk and runcost to minc
    G1p, G2p = [], []
    for (r,o,s), g1, g2 in zip(trips, G1, G2):
        G1p.append(g1 / max(G1dict[o]) * mincs[op])
        G2p.append(g2 / max(G2dict[o]) * mincs[op])
    c = map(sum, zip(c, G1p, G2p))
    cp = list(c)
    d = dict([(t, cost) for cost, t in zip(c, trips)])
    c = matrix(c)
    
    nlog.debug("Constraints: %f\t%f" , maxdisk, maxoverhead)
    for (r, op, s), pqcost in zip(trips, G1p):
        if 'Extract' not in str(op): continue
        nlog.debug('%s\t%s\t%.15f\t%f\t%f\t%.15f', op, str(s).ljust(25),
                   pqcost,
                   mp.get_disk(op,s),
                   mp.get_provcost(op, s),
                   d[(r, op, s)])
        

    strategies, torm = nlp_exec_cvx(c, ops, matstrats, trips, G, h, A, b)

    return strategies, torm
    

def nlp_exec_cvx(c, ops, matstrats, trips, G, h, A, b):
    """
    @return newassignments, prov_to_remove
    prov_to_remove is a list of (op, runid) pairs
    """
    solvers.options['show_progress'] = False
    solvers.options['LPX_K_MSGLEV'] = 0
    solvers.options['MessageLevel'] = 3
    glpk.options['LPX_K_MSGLEV'] = 0

    I = set()
    B = set(range(len(c)))



    status, x = ilp(c, G, h, A, b, I, B)

    nlog.info( "status     \t%s", status )
    nlog.info( "resources  \t%s", np.array((G*x).trans())[0] )
    nlog.info( "constraints\t%s", np.array(h.trans())[0] )
    nlog.info( "cost       \t%s", np.array(x.trans() * c)[0] )

    pairsize = len(ops)*len(matstrats)
    newassignments = x[:pairsize]
    rmassignments = x[pairsize:]

    newassignments = numpy.array(newassignments).reshape((len(ops), len(matstrats))).tolist()
    strategies = assign_strategies(newassignments, ops, matstrats)

    torm = []
    for strat, ass in zip(trips[pairsize:], rmassignments):
        if strat[2] == Strat.query() and ass:
            torm.append((strat[1], strat[0]))

    return strategies, torm
    

def assign_strategies(matrix, ops, strats):
    ret = {}
    for x, op in enumerate(ops):
        ret[op] = []
        offset = x * len(strats)
        for y, s in enumerate(strats):
            if matrix[x][y]:
                ret[op].append(s)
    return ret
