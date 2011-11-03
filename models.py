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



def disk_model(strat, fanin, fanout, density, noutput):
    cost = 0.0
    for bucket in strat.buckets:
        for desc in bucket.descs:
            cost += disk_model_desc(desc, fanin, fanout, density, noutput)
    return cost
    

    
def disk_model_desc(desc, fanin, fanout, density, noutput):
    #desc.mode, desc.spec, desc.backward
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        return 0
    elif desc.mode == Mode.PT_MAPFUNC:
        if desc.spec.outcoords == Spec.COORD_MANY:
            return noutput * 4.25 + (noutput / fanin) * 8
        return noutput * (8 + 4.25)

    ncalls = math.ceil(noutput / float(fanout))
    spec = (desc.spec.outcoords, desc.spec.payload)
    disk = 0

    if desc.backward:
        if spec == (Spec.COORD_ONE, Spec.COORD_MANY):
            disk = fanout * 4.25 + fanout * fanin * 4.25
        elif spec == (Spec.COORD_ONE, Spec.KEY):
            disk = fanout * 4.25 + fanin * 4.25 + fanout * 28 + 28
        elif spec == (Spec.COORD_ONE, Spec.BOX):
            disk = fanout * (4.25 + 16)
        elif spec == (Spec.COORD_MANY, Spec.COORD_MANY):
            disk = fanout * 4.25 + fanin * 4.25
        elif spec == (Spec.COORD_MANY, Spec.KEY):
            disk = fanout * 4.25 + fanin * 4.25 + 28 * 2
        elif spec == (Spec.COORD_MANY, Spec.BOX):
            disk = fanout * 4.25 + 16
    else:
        if spec == (Spec.COORD_ONE, Spec.COORD_MANY):
            disk = fanin * 4.25 + fanout * fanin * 4.25
        elif spec == (Spec.COORD_ONE, Spec.KEY):
            disk = fanin * 4.25 + fanout * 4.25 + fanin * 28 + 28
        elif spec == (Spec.COORD_ONE, Spec.BOX):
            disk = fanin * (4.25 + 16)
        elif spec == (Spec.COORD_MANY, Spec.COORD_MANY):
            disk = fanin * 4.25 + fanout * 4.25
        elif spec == (Spec.COORD_MANY, Spec.KEY):
            disk = fanin * 4.25 + fanout * 4.25 + 28 * 2
        elif spec == (Spec.COORD_MANY, Spec.BOX):
            disk = fanin * 4.25 + 16

    return int( math.ceil(disk * ncalls / 4096.0) * 4096 ) #+ ncalls * 312 + 12288

def write_model(strat, fanin, fanout, density, noutput, opcost):
    cost = 0.0
    for bucket in strat.buckets:
        for desc in bucket.descs:
            cost += write_model_desc(desc, fanin, fanout, density, noutput, opcost)
    return cost

def write_model_desc(desc, fanin, fanout, density, noutput, opcost):
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        return 0.00001
    if desc.mode in (Mode.STAT, Mode.NOOP):
        return 0.0
    
    diskcost = (disk_model_desc(desc, fanin, fanout, density, noutput) / 1048576.0) / 13.0

    # cost to just generate the coordinates
    #  - measured to be 1e-5 per pointer * (noutput * fanin)
    #    for the general case: transpose etc
    genprovcost = noutput * fanin * 1e-5

    # cost to serialize per cell
    coord_one  = (5.6e-6, 5.6e-6)
    coord_many = (9.8e-6, 5.9e-7)
    box        = (4.6e-6, 2.3e-6)
    key        = (2.5e-5, 6.2e-7)
    none       = (4.1e-6, 5.5e-7)


    if desc.mode == Mode.PT_MAPFUNC:
        if desc.spec.outcoords == Spec.COORD_MANY:
            if desc.backward:
                if fanout <= 10:
                    return diskcost + noutput * coord_many[0]
                return diskcost + noutput * coord_many[1]
        if desc.spec.outcoords == Spec.COORD_ONE:
            return diskcost + noutput * coord_one[0]

    # update stats costs
    statscost = (noutput / fanout) * (fanout + fanin) * 0.000003

    sercost = ser_model_desc(desc, fanin, fanout, density, noutput) * 10 # seem to be 10x off reality
    return sercost + diskcost + opcost * 0.4 + statscost

def ser_model_desc(desc, fanin, fanout, density, noutput):
    if desc.mode in (Mode.STAT, Mode.NOOP):
        return 0.0
    # cost to serialize per cell
    coord_one  = (5.6e-6, 5.6e-6)
    coord_many = (9.8e-6, 5.9e-7)
    box        = (4.6e-6, 7e-7)# 2.3e-6)
    key        = (2.5e-5, 6.2e-7)
    none       = (4.1e-6, 5.5e-7)
    coord_one = coord_many = box = key = none = (9.8e-6, 5.9e-7)

    ncalls = int(math.ceil(noutput / float(fanout)))
    spec = (desc.spec.outcoords, desc.spec.payload)
    if desc.backward:
        if spec == (Spec.COORD_ONE, Spec.COORD_MANY):
            cost = fanout * coord_one[0] + fanout * fanin * coord_many[1]
        elif spec == (Spec.COORD_ONE, Spec.KEY):
            cost = fanout * coord_one[0] + fanin * key[1]
        elif spec == (Spec.COORD_ONE, Spec.BOX):
            cost = fanout * (coord_one[0] + fanin * box)
        elif spec == (Spec.COORD_MANY, Spec.COORD_MANY):
            cost = fanout * coord_many[1] + fanin * coord_many[1]
        elif spec == (Spec.COORD_MANY, Spec.KEY):
            cost = fanout * coord_many[1] + fanin * key[1]
        elif spec == (Spec.COORD_MANY, Spec.BOX):
            cost = fanout * coord_many[1] + fanin * box[1]
        else:
            raise RuntimeError, str(desc)
    else:
        if spec == (Spec.COORD_ONE, Spec.COORD_MANY):
            cost = fanin * coord_one[0] + fanout * fanin * coord_many[1] * ncalls * 0.3
        elif spec == (Spec.COORD_ONE, Spec.KEY):
            cost = fanin * coord_one[0] + fanout * key[1] * ncalls * 0.3
        elif spec == (Spec.COORD_ONE, Spec.BOX):
            cost = fanin * (coord_one[0] + fanout * box * ncalls * 0.3 )
        elif spec == (Spec.COORD_MANY, Spec.COORD_MANY):
            cost = fanin * coord_many[1] + fanout * coord_many[1]
        elif spec == (Spec.COORD_MANY, Spec.KEY):
            cost = fanin * coord_many[1] + fanout * key[1]
        elif spec == (Spec.COORD_MANY, Spec.BOX):
            cost = fanin * coord_many[1] + fanout * box[1]
        else:
            raise
                
    return cost * ncalls


def forward_model(strat, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None):
    scost = 10
    for bucket in strat.buckets:
        bcost = 0.0
        for desc in bucket.descs:
            cost = forward_model_desc(desc, fanin, fanout, density, noutput,
                                      runtime, nqs, sel, inputarea=inputarea)
            bcost += cost
        scost = min(scost, bcost)
    return scost

def forward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None):
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        return backward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea)
    if desc.mode in (Mode.STAT, Mode.NOOP):
        return 0.0
    
    
    outcost, incost = get_parse_costs(desc)
    diskcost = (disk_model_desc(desc, fanin, fanout, density, noutput) / 1048576.0) / 55.0
    ncalls = int(math.ceil(noutput / float(fanout)))
    dedupcost = nqs * sel * 1e-6 # cost of load() in dedup.  includes overlap of results

    if desc.mode == Mode.PT_MAPFUNC:
        itercost = (noutput / fanout) * nqs * 0.000001
        return diskcost + itercost + noutput / fanout * 2e-6 + dedupcost
        
    coord_one = 8.51511955261e-6
    coord_many = 1.36479878426e-6
    box = 4.33921813965e-6
    key = 7.15684890747e-6
    none = 2.87890434265e-6

    # == if hash join ==
    if not desc.backward and desc.spec.outcoords == Spec.COORD_ONE:
        coord_one = 7.102e-6
        parsecost = nqs * coord_one

        extractcost = fanout * outcost
        extractcost *= nqs * sel
        datacost = nqs * sel *  0.000007 # cost to get data from bdb        

        qcost = parsecost + extractcost + datacost + dedupcost
        return qcost
    elif not desc.backward:
        # if optimized for forward queries
        parsecost = 0

        # X --> Y
        # parse all of X
        # parse some of Y

        parsecost = fanin * (noutput / fanout) * incost
        extractcost = fanout * outcost * sel
        matchcost = (noutput / fanout) * 0.000001   # cost to execute join predicate        
        qcost = parsecost + extractcost + matchcost
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
        matchcost = nptrs * nqs * 0.000001   # cost to execute join predicate                

        qcost = parsecost + extractcost + matchcost 
    
    return diskcost + qcost + dedupcost
    
def box_model(density, fanin, inputarea, runtime, nptrs):
    a,b,c = 0.08147001934235978, 0.00852998065764022, 1.4085629809324707
    baseline = a + b / (density ** c)
    a,b,c = -0.12499999999999964, 1.1249999999999998, 0.5642714304385626
    nintersect = a + b / (density ** c)
    nintersect = min(nptrs, nintersect)

    cost = baseline + nintersect * (min(fanin / density, inputarea) / inputarea) * runtime
    return cost

def backward_model(strat, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None):
    scost = 10
    for bucket in strat.buckets:
        bcost = 0.0
        for desc in bucket.descs:
            cost = backward_model_desc(desc, fanin, fanout, density, noutput,
                                       runtime, nqs, sel, inputarea=inputarea)
            bcost += cost
        scost = min(scost, bcost)
    return scost


def backward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None):
    if desc.mode in (Mode.STAT, Mode.NOOP):
        return 0.0
    if desc.mode == Mode.QUERY:
        return runtime + nqs * sel * 1e-6
    if desc.mode == Mode.FULL_MAPFUNC:
        return nqs * 1e-6
    
    outcost, incost = get_parse_costs(desc)
    diskcost = (disk_model_desc(desc, fanin, fanout, density, noutput) / 1048576.0) / 55.0
    ncalls = int(math.ceil(noutput / float(fanout)))
    dedupcost = nqs * fanin * 1e-6 # cost of load() in dedup.  includes overlap of results

    if desc.mode == Mode.PT_MAPFUNC:
        if desc.backward and desc.spec.outcoords== Spec.COORD_ONE:
            sercost = nqs * 0.000008  # cost to serialize coords into bdb key
            extractcost = nqs * sel * 0.000016 # cost to extract data from pointers
            datacost = nqs * sel *  0.000007 # cost to get data from bdb
            return diskcost + sercost + extractcost + datacost + dedupcost

        # n iterations * nqs in each iteration
        itercost = (noutput / fanout) * nqs * 0.000001
        #print '\titercost', itercost, (noutput / fanout * nqs), noutput, nqs, fanout
        return diskcost + itercost + noutput / fanout * 2e-6  + dedupcost


    coord_one = 8.51511955261e-6
    coord_many = 1.36479878426e-6
    box = 4.33921813965e-6
    key = 7.15684890747e-6
    none = 2.87890434265e-6

    # == if hash join ==
    if desc.backward and desc.spec.outcoords == Spec.COORD_ONE:
        coord_one = 7.102e-6
        parsecost = nqs * coord_one

        if desc.spec.payload == Spec.BOX:
            extractcost = runtime * min((fanin / density), inputarea)  / float(inputarea)
        else:
            extractcost = fanin * incost
        extractcost *= nqs * sel
        datacost = nqs * sel *  0.000007 # cost to get data from bdb

        qcost = parsecost + extractcost + extractcost + datacost + dedupcost
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
        matchcost = (noutput / fanout) * 0.000001   # cost to execute join predicate

        qcost = parsecost + extractcost + matchcost
        qcost *= nqs
        
    else:
        nptrs = desc.spec.outcoords == Spec.COORD_ONE and noutput or noutput / fanout
        parsecost = fanin * incost * nqs * nptrs
        extractcost = fanout * outcost * nqs * sel
        matchcost = nptrs * nqs * 0.000001   # cost to execute join predicate        
        qcost = parsecost + extractcost + matchcost
    
    return diskcost + qcost + dedupcost

if __name__ == '__main__':
    strats = [Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
              Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False)]
    for desc in strats:
        for fanout in [1, 25, 64]:
            print disk_model_desc(desc, 169, fanout, 1.0, 1000), desc
