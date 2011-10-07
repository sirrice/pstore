from util import print_matrix
import os
from runtime import *
import numpy


def run_nlp(stats, w, mp, maxdisk, maxoverhead):
    """
    maxoverhead is percentage of average runtime cost of optimizable operators
    """
    ops = w.get_optimizable_ops()
    matstrats = w.get_matstrats()
    pairs = [(op, s) for op in ops for s in matstrats] 
    triples = [(r, op, s) for r in xrange(1,w._runid) for op in ops for s in matstrats]

    x = [STRAT_Q == s and 1 or 0 for (r,op,s) in triples]
    x.extend([0] * len(pairs))

    F = [mp.get_pqcost(op, s) for (r, op, s) in triples]
    F.extend([mp.get_pqcost(op,s) for op,s in pairs])
    P = [1.0 for t in triples]
    P.extend([1.0 for p in pairs])

    avg_runtime = numpy.mean(filter(lambda x: x > 0, [mp.get_opcost(op, s) for op,s in pairs]))
    maxoverhead *= avg_runtime

    A1 = [mp.get_disk(op, s) for (r,op,s) in triples]
    A1.extend([mp.get_disk(op,s) for (op,s) in pairs])
    A2 = [0.0 for t in triples]
    A2.extend([mp.get_provcost(op, s) for (op,s) in pairs])
    A = [A1, A2]
    B = [[maxdisk, maxoverhead]]



    Aeq = []
    for i in xrange(len(ops)):
        row = [0 for t in triples]
        for col in xrange(len(pairs)):
            if len(matstrats) * i <= col and col < len(matstrats) * (i+1):
                row.append(1)
            else:
                row.append(0)
        Aeq.append(row)
    Beq = [[1] * len(ops)]

    print "Max constraints Disk(%f)\tOverhead(%f)" % (maxdisk, maxoverhead)
    
    print "%s     \t" * 7 % ("operator      ","strategy","proqcost",
                             "provsize","overhead","savecost","runtime")
    for op, s in pairs:
        print "%s\t%s\t% 5f\t% 5f\t% 5f" % (op, s,
                                            mp.get_pqcost(op,s),
                                            mp.get_disk(op,s),
                                            mp.get_provcost(op, s))
                                                        #stats.get_save(op, s),
                                                        #stats.get_runtime(op, s))

    strategies = nlp_exec(F, ops, matstrats, A, B, Aeq, Beq)
    return strategies




def old_run_nlp(stats, w, maxdisk, maxoverhead):
    """
    maxoverhead is percentage of average runtime cost of optimizable operators
    XXX: worked for long runnnig workflows
    """
    ops = w.get_optimizable_ops()
    matstrats = w.get_matstrats()
    pairs = [(op, s) for op in ops for s in matstrats] 
    triples = [(r, op, s) for r in xrange(1,w._runid) for op in ops for s in matstrats]

    x = [Runtime.instance().get_strategy(op, r) == s and 1 or 0 for (r,op,s) in triples]
    x.extend([0] * len(pairs))

    F = [w.wrapper(op).cost(s, r) for (r, op, s) in triples]
    F.extend([w.wrapper(op).cost(s) for op,s in pairs])
    P = [1.0 for t in triples]
    P.extend([1.0 for p in pairs])

    avg_runtime = numpy.mean(filter(lambda x: x > 0, [stats.get_runtime(op, s) for op,s in pairs]))
    maxoverhead *= avg_runtime

    A1 = [w.wrapper(op).provsize(s,r) for (r,op,s) in triples]
    A1.extend([w.wrapper(op).provsize(s) for (op,s) in pairs])
    A2 = [0.0 for t in triples]
    A2.extend([w.wrapper(op).provoverhead(s) for (op,s) in pairs])
    A = [A1, A2]
    B = [[maxdisk, maxoverhead]]

    Aeq = []
    for i in xrange(len(ops)):
        row = [0 for t in triples]
        for col in xrange(len(pairs)):
            if len(matstrats) * i <= col and col < len(matstrats) * (i+1):
                row.append(1)
            else:
                row.append(0)
        Aeq.append(row)
    Beq = [[1] * len(ops)]

    print "Max constraints Disk(%f)\tOverhead(%f)" % (maxdisk, maxoverhead)
    
    print "%s     \t" * 7 % ("operator      ","strategy","proqcost",
                             "provsize","overhead","savecost","runtime")
    for op, s in pairs:
        print "%s\t%s\t% 5f\t% 5f\t% 5f\t% 5f\t% 5f" % (op, s,
                                                  w.wrapper(op).cost(s),
                                                  w.wrapper(op).provsize(s),
                                                  w.wrapper(op).provoverhead(s),
                                                  stats.get_save(op, s),
                                                  stats.get_runtime(op, s))


    strategies = nlp_exec(F, ops, matstrats, A, B, Aeq, Beq)
    return strategies


def gen_F(ops, matstrats, disk, opcosts, probs):
    # operator/strategy costs
    opcosts = [x*y for x,y in zip(opcosts, probs)]
    for cost, d, (op, mat) in zip(opcosts, disk, [(op,s) for op in ops for s in matstrats]):
        print "opcosts\t%s\t%s\t%f\t%f" % (str(op).ljust(15), mat, d, cost)

    F = [opcost + (d / (1000.0 * max(disk))) for d, opcost in zip(disk, opcosts)]
    return F


def nlp_exec(F, ops, matstrats, A, B, Aeq, Beq):
    ncombos = len(matstrats) * len(ops)

    X0 = [0 for op in ops for s in xrange(len(F))] 

    #A = [list(disk)]
    #B = [[maxdisk]]

    # Aeq = []
    # for i in xrange(len(ops)):
    #     row = []
    #     for col in xrange(ncombos):
    #         if len(matstrats) * i <= col and col < len(matstrats) * (i+1):
    #             row.append(1)
    #         else:
    #             row.append(0)
    #     Aeq.append(row)
    # Beq = [[1] * len(ops)]

    nruns = len(F) / (len(ops) * len(matstrats))

    # main program
    f = file('lsst_1.m', 'w')
    print >> f, "F   = [%s]" % (','.join(map(str, F)))
    print >> f, 'X0  = [%s]' % (';'.join(map(str, X0)))
    print >> f, print_matrix('Aeq', Aeq)
    print >> f, print_matrix('Beq', Beq)
    print >> f, print_matrix('A', A)
    print >> f, print_matrix('B', B)
    print >> f, "B = B'"
    print >> f, "Beq = Beq'"        
    print >> f, "[x, fval, exitflag, output] = bintprog(F,A,B,Aeq,Beq)"
    # print >> f, "options = optimset('Algorithm','active-set');"
    print >> f, """
    res = full(x)
    res = reshape(res, %d, %d, %d)
    res = res(:,:,%d)
    res = res'

    f = fopen('nlpresult.dat', 'w')
    fprintf(f, '%s', size(res))
    fprintf(f, '\\n')
    fprintf(f, '%s', round(res))
    fclose(f)
    """ % (len(matstrats), len(ops), nruns, nruns, '%d\\n', '%d\\n')
    f.close()

    f = file("runlsst.m", 'w')
    print >> f, """function assignment = runlsst()
    cd ~/mitnotes/research/provenance/src/sim
    lsst_1""" 
    f.close()


    # run matlab and solve
    #os.system('echo "runlsst" | /Applications/MATLAB_R2010b.app/bin/matlab -nodesktop 2> /dev/null > /dev/null')
    os.system('echo "runlsst" | matlab -nodesktop 2> /dev/null > /dev/null')

    # read solution
    assignmentmatrix = parse_nlpresult()
    print ( "==Assignment Matrix==")
    for row in assignmentmatrix:
        print (str(row))
    print ('\n')

    strategies = assign_strategies(assignmentmatrix, ops, matstrats)
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

def parse_nlpresult(fname='nlpresult.dat'):
    f = file(fname, 'r')
    rows = int(f.readline())
    cols = int(f.readline())
    f.readline()

    arr = [[0 for c in xrange(cols)] for r in xrange(rows)]
    for c in xrange(cols):
        for r in xrange(rows):
            arr[r][c] = float(f.readline())
    f.close()
    return arr


def extract_nlp_params(fname='out'):
    """
    opname -> strat -> (disk, opcost)
    """
    f = file(fname, 'r')
    ret = {}
    for l in f:
        if not l.startswith('opcosts'): continue
        arr = l.split()[1:]
        opname, strat, disk, opcost = arr[0], str_to_strat(arr[1]), float(arr[2]), float(arr[3])
        if opname not in ret:
            ret[opname] = {}
        ret[opname][strat] = (disk, opcost)
    return ret

def str_to_strat(s):

    if s == 's:Func':
        return STRAT_F
    elif s == 's:Query':
        return STRAT_Q
    elif s == 's:PtrList':
        return STRAT_PLIST
    elif s == 's:PtrSet':
        return STRAT_PSET
    elif s == 's:PtrGrid':
        return STRAT_PGRID
    elif s == 's:Box':
        return STRAT_BOX
    elif s == 's:Test':
        return STRAT_PTEST
    raise RuntimeError, "can't parse strategy: %s" % s
