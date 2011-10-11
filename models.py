import math
from runtime import *

# STRAT_DIFF is hard to predict so I'm punting on it



def get_parse_costs(desc):
    coord_one = 8.51511955261e-6
    coord_many = 1.36479878426e-6
    box = 4.33921813965e-6
    key = 7.15684890747e-6
    none = 2.87890434265e-6
    if desc.spec.outcoords == Spec.COORD_ONE:
        outcost = coord_one
    elif desc.spec.outcoords == Spec.COORD_MANY:
        outcost = coord_many
    elif desc.spec.outcoords == Spec.BOX:
        outcost = box
    elif desc.spec.outcoords == Spec.KEY:
        outcost = coord_many
    else:
        outcost = none

    if desc.spec.payload == Spec.COORD_ONE:
        incost = coord_one
    elif desc.spec.payload == Spec.COORD_MANY:
        incost = coord_many
    elif desc.spec.payload == Spec.BOX:
        incost = box
    elif desc.spec.payload == Spec.KEY:
        incost = coord_many
    else:
        incost = none
    return outcost, incost



def disk_model(strat, fanin, fanout, density, noutput, grid=10):
    return disk_model_desc(strat.buckets[0].descs[0], fanin, fanout, density, noutput, grid=10)
    
def disk_model_desc(desc, fanin, fanout, density, noutput, grid=10):
    #desc.mode, desc.spec, desc.backward
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        disk = 0
    elif desc.mode == Mode.PT_MAPFUNC:
        disk = noutput * 8 + 4
    if desc.mode == Mode.PTR:

        #if desc.spec.outcoords == Spec.COORD_ONE:
        #ncalls = noutput
        #else:
        ncalls = int(math.ceil(noutput / float(fanout)))

        if not desc.backward:
            fanin, fanout = fanout, fanin



        outsize = 0
        mult = 1
        if desc.spec.outcoords == Spec.COORD_ONE:
            outsize = 4
            mult *= fanout
        if desc.spec.outcoords == Spec.COORD_MANY:
            outsize = fanout * 4 + 4
        if desc.spec.outcoords == Spec.KEY:
            outsize = fanout * 4 + 4 + 20
        if desc.spec.outcoords == Spec.BOX:
            outsize = 16


        insize = 0
        if desc.spec.payload == Spec.COORD_ONE:
            insize = 4
            mult *= fanin
        if desc.spec.payload == Spec.COORD_MANY:
            insize = fanin * 4 + 4
        if desc.spec.payload == Spec.KEY:
            insize = fanin * 4 + 4 + 20
        if desc.spec.payload == Spec.BOX:
            insize = 16


        disk = (outsize + insize) * mult * ncalls
    return disk

def write_model(strat, fanin, fanout, density, noutput):
    return write_model_desc(strat.buckets[0].descs[0], fanin, fanout, density, noutput)

def write_model_desc(desc, fanin, fanout, density, noutput, grid=10):
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        return 0.001
    diskcost = (disk_model_desc(desc, fanin, fanout, density, noutput) / 1048576.0) / 13.0
    
    # cost to serialize per cell
    coord_one = 7.102e-6
    box = 3.17e-6
    other = 6.92e-7

    if desc.mode == Mode.PT_MAPFUNC:
        return diskcost+ noutput * other 
    

    if desc.spec.outcoords == Spec.COORD_ONE:
        outcost = coord_one
    elif desc.spec.outcoords == Spec.BOX:
        outcost = box
    else:
        outcost = other

    if desc.spec.payload == Spec.COORD_ONE:
        incost = coord_one
    elif desc.spec.payload == Spec.BOX:
        incost = box
    else:
        incost = other
    

    ncalls = int(math.ceil(noutput / float(fanout)))

    if not desc.backward:
        fanin, fanout = fanout, fanin

    mult = 1
    if not desc.backward and desc.spec.outcoords == Spec.COORD_ONE:
        # overhead of collisions!  How many collisions are there in inputs?
        mult = max(1, fanout / (160 / 10.0))
    
    cost_per_call = fanout * outcost + fanin * incost * mult

    return diskcost + ( ncalls * cost_per_call)

def forward_model(strat, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None, grid=10):
    return forward_model_desc(strat.buckets[0].descs[0], fanin, fanout, density, noutput,
                              runtime, nqs, sel, inputarea=inputarea)

def forward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel,
                  inputarea=None, grid=10):

    if desc.mode == Mode.QUERY:
        return runtime 
    outcost, incost = get_parse_costs(desc)
    diskcost = (disk_model_desc(desc, fanin, fanout, density, noutput) / 1048576.0) / 55.0
    ncalls = int(math.ceil(noutput / float(fanout)))
    #if not desc.backward:
    #    fanin, fanout = fanout, fanin

    if desc.mode == Mode.PT_MAPFUNC:
        return diskcost + nqs * sel * 6e-7

    # == if hash join ==
    if not desc.backward and desc.spec.outcoords == Spec.COORD_ONE:
        coord_one = 7.102e-6
        parsecost = nqs * coord_one

        extractcost = fanout * outcost
        extractcost *= nqs * sel

        qcost = parsecost + extractcost
        return qcost
    elif not desc.backward:
        # if optimized for forward queries
        parsecost = 0

        # X --> Y
        # parse all of X
        # parse some of Y

        parsecost = fanin * (noutput / fanout) * incost
        extractcost = fanout * outcost * sel
        qcost = parsecost + extractcost
        qcost *= nqs
        
    else:
        # if optimized for backward queries
        # iterate through everything, read the inputs

        # X --> Y
        # X * Y --> Y
        # parse X * Y
        # parse some of Y

        nptrs = desc.spec.outcoords == Spec.COORD_ONE and noutput or noutput / fanout
        parsecost = fanin * incost * nqs * nptrs

        if desc.spec.payload == Spec.BOX:
            extractcost = box_model(density, fanin, inputarea, runtime, nptrs)
            # perc = min((fanin / density), inputarea) / float(inputarea)
            # nreexec = (noutput / fanout) * perc * nqs
            # extractcost = runtime * perc * nreexec
        else:
            extractcost = fanout * outcost * nqs * sel

        qcost = (parsecost + extractcost) 
    
    return diskcost + qcost
    
def box_model(density, fanin, inputarea, runtime, nptrs):
    a,b,c = 0.08147001934235978, 0.00852998065764022, 1.4085629809324707
    baseline = a + b / (density ** c)
    a,b,c = -0.12499999999999964, 1.1249999999999998, 0.5642714304385626
    nintersect = a + b / (density ** c)
    nintersect = min(nptrs, nintersect)

    cost = baseline + nintersect * (min(fanin / density, inputarea) / inputarea) * runtime
    return cost

def backward_model(strat, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None, grid=10):
    return backward_model_desc(strat.buckets[0].descs[0], fanin, fanout, density, noutput,
                              runtime, nqs, sel, inputarea=inputarea)

def backward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel,
                   inputarea=None, grid=10):

    if desc.mode == Mode.QUERY:
        return runtime

    outcost, incost = get_parse_costs(desc)
    diskcost = (disk_model_desc(desc, fanin, fanout, density, noutput) / 1048576.0) / 55.0
    ncalls = int(math.ceil(noutput / float(fanout)))

    if desc.mode == Mode.PT_MAPFUNC:
        return diskcost + nqs * sel * 6e-7

    # == if hash join ==
    if desc.backward and desc.spec.outcoords == Spec.COORD_ONE:
        coord_one = 7.102e-6
        parsecost = nqs * coord_one

        if desc.spec.payload == Spec.BOX:
            extractcost = runtime * min((fanin / density), inputarea)  / float(inputarea)
        else:
            extractcost = fanin * incost
        extractcost *= nqs * sel

        qcost = parsecost + extractcost
        return qcost

    elif desc.backward:
        parsecost = fanout * (noutput / fanout) * outcost

        if desc.spec.payload == Spec.BOX:
            a,b,c = 0.10013271527592049, 0.00022428472407952143, 1.5845707477154047
            baseline = a + b / (density ** c)
            perc = min((fanin / density), inputarea) / float(inputarea)
            extractcost = (baseline + runtime * perc) * sel 
        else:
            extractcost = fanin * incost * sel

        qcost = parsecost + extractcost
        qcost *= nqs
        
    else:
        nptrs = desc.spec.outcoords == Spec.COORD_ONE and noutput or noutput / fanout
        parsecost = fanin * incost * nqs * nptrs
        extractcost = fanout * outcost * nqs * sel
        qcost = parsecost + extractcost
    
    return diskcost + qcost

