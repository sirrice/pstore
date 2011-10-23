import os, logging, numpy
from util import print_matrix
from runtime import *

nlog = logging.getLogger('nlp')
logging.basicConfig()
nlog.setLevel(logging.ERROR)


def run_nlp(stats, w, mp, maxdisk, maxoverhead):
    """
    maxoverhead is percentage of average runtime cost of optimizable operators
    """
    ops = w.get_optimizable_ops()
    matstrats = w.get_matstrats()
    pairs = [(op, s) for op in ops for s in matstrats] 

    x = [0] * len(pairs)

    F = [mp.get_pqcost(op,s) for op,s in pairs]
    P = [1.0 for p in pairs]

    avg_runtime = numpy.mean(filter(lambda x: x > 0, [mp.get_opcost(op, s) for op,s in pairs]))
    maxoverhead *= avg_runtime

    A1 = [mp.get_disk(op,s) for (op,s) in pairs]
    A2 = [mp.get_provcost(op, s) for (op,s) in pairs]
    A = [A1, A2]
    B = [[maxdisk, maxoverhead]]



    Aeq = []
    for i in xrange(len(ops)):
        row = []
        for col in xrange(len(pairs)):
            if len(matstrats) * i <= col and col < len(matstrats) * (i+1):
                row.append(1)
            else:
                row.append(0)
        Aeq.append(row)
    Beq = [[1] * len(ops)]

    nlog.debug( "Max constraints Disk(%f)\tOverhead(%f)", maxdisk, maxoverhead )
    
    nlog.debug( "%s     \t" * 7, "operator      ","strategy","proqcost",
                "provsize","overhead","savecost","runtime")
    for op, s in pairs:
        nlog.debug( "%s\t%s\t% 5f\t% 5f\t% 5f", op, s,
                                            mp.get_pqcost(op,s),
                                            mp.get_disk(op,s),
                                            mp.get_provcost(op, s))
    strategies = nlp_exec(F, ops, matstrats, A, B, Aeq, Beq)
    return strategies



def gen_F(ops, matstrats, disk, opcosts, probs):
    # operator/strategy costs
    opcosts = [x*y for x,y in zip(opcosts, probs)]
    for cost, d, (op, mat) in zip(opcosts, disk, [(op,s) for op in ops for s in matstrats]):
        nlog.debug( "opcosts\t%s\t%s\t%f\t%f" , str(op).ljust(15), mat, d, cost )

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
    nlog.debug ( "==Assignment Matrix==")
    for row in assignmentmatrix:
        nlog.debug (str(row))
    nlog.debug ('\n')

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
