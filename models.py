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


def index_model(strat, fanin, fanout, density, noutput):
    cost = 0.0
    for bucket in strat.buckets:
        for desc in bucket.descs:
            cost += index_model_desc(desc, fanin, fanout, density, noutput)
    return cost

def index_model_desc(desc, fanin, fanout, density, noutput):
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        return 0
    if desc.spec.outcoords == Spec.COORD_ONE:
        return 0
    return ( noutput / fanout ) + (24 + 60.18) + 7340


def disk_model(strat, fanin, fanout, density, noutput):
    cost = 0.0
    for bucket in strat.buckets:
        for desc in bucket.descs:
            cost += disk_model_desc(desc, fanin, fanout, density, noutput)
    return cost
    


def get_disk_params(spec, fanin, fanout, density, noutput, payload):
    a,b = 1,1
    okey = ikey = ptr = idx = 0
    nentries = 0
    if spec[0] == Spec.COORD_ONE:
        nptrs = noutput
        nentries = 1
        if spec[1] == Spec.KEY:
            a, b = 1.5,  15.9
            ikey = ( 18 + 4 * (1 + fanin) ) / fanout
            ptr = ( 18 + 4 )
            nentries = 1 + 1.0 / fanout
        elif spec[1] == Spec.COORD_MANY:
            a,b = 2.31,  9.986
            ptr = 4 + 4 * (fanin + 1)
        elif spec[1] == Spec.GRID:
            a,b = 2.31701812,  9.98665472
            negs = 4 * (2 + 1 + (math.ceil((fanin/density) ** 0.5) ** 2) - fanin)
            ptr = 4 + negs
        elif spec[1] == Spec.BOX:
            a,b = 1.75343538, 41.79651412
            ptr = 12
        elif spec[1] == Spec.BINARY:
            a,b = 1.75343538, 41.79651412
            ptr = 4 + 4 + payload
    elif spec[0] == Spec.COORD_MANY:
        nptrs = noutput / fanout
        outsize = ( 4 + 4 * (fanout + 1) )
        idx = 4 + outsize
        nentries = 2
        if spec[1] == Spec.KEY:
            a,b = 1.3, 50
            ikey = ( 18 + 4 * (1 + fanin) ) 
            ptr = 18 + outsize
            nentries = 3
        elif spec[1] == Spec.COORD_MANY:
            a,b =  2.02,  36.1
            ptr = 4 * (1 + fanin) + outsize
        elif spec[1] == Spec.GRID:
            a,b =  2.02,  36.1
            negs = 4 * ( 2 + 1 + (math.ceil((fanin/density) ** 0.5) ** 2) - fanin)
            ptr = negs + outsize
        elif spec[1] == Spec.BOX:
            a,b = 1.75343538, 41.79651412
            ptr = 8 + outsize
        elif spec[1] == Spec.BINARY:
            a,b = 1.75343538, 41.79651412
            ptr = outsize + 4 + payload

    elif spec[0] == Spec.KEY:
        nptrs = noutput / fanout
        idx = (4 + 18)
        okey = ( 18 + 4 * (1 + fanout) )
        nentries = 3
        if spec[1] == Spec.KEY:
            a,b = 2.99000236, -9.01436915
            ikey = ( 18 + 4 * (1 + fanin) )
            ptr = 18 + 18
            nentries += 1
        elif spec[1] == Spec.BOX:
            a,b = 2.68105634,  4.94986864
            ptr = 18 + 8
        elif spec[1] == Spec.COORD_MANY:
            a,b = 1.35 , 46.3
            ptr = ( 18 + (4 + 4*(fanin + 1)) )
        elif spec[1] == Spec.GRID:
            a,b = 1.35 , 46.3
            negs = 4 * ( 2 + 1 + (math.ceil((fanin/density) ** 0.5) ** 2) - fanin)
            ptr = negs + 18
        elif spec[1] == Spec.BINARY:
            a,b = 1.75343538, 41.79651412
            ptr = 18 + 4 + payload

    return a,b,okey,ikey,ptr,idx,nentries    
    
def disk_model_desc(desc, fanin, fanout, density, noutput, payload=2):
    #desc.mode, desc.spec, desc.backward
    if desc.mode in (Mode.QUERY, Mode.FULL_MAPFUNC):
        return 0

    spec = (desc.spec.outcoords, desc.spec.payload)
    idxsize = index_model_desc(desc, fanin, fanout, density, noutput)    
    disk = 0

    if spec[0] == Spec.COORD_ONE:
        nptrs = noutput
    else:
        nptrs = int(math.ceil(noutput / float(fanout)))


    if desc.backward:
        a,b,okey,ikey,ptr,idx,nentries = get_disk_params(spec, fanin, fanout, density, noutput, payload)
    else:
        a,b,okey,ikey,ptr,idx,nentries = get_disk_params(spec, fanout, fanin, density, noutput, payload)
        # fix things
        if spec == (Spec.COORD_ONE, Spec.KEY):
            ikey = ( 18 + 4 * (1 + fanin) ) / fanout
            nentries = 1 + 1.0 / fanout
        if spec[0] != Spec.COORD_ONE:
            nptrs = noutput / fanout
        else:
            nptrs = noutput / fanout * fanin
        
    disk = nptrs * ( ( okey + ikey + ptr + idx ) * a + nentries * b )    
    return idxsize + int( math.ceil(disk/ 4096.0) * 4096 ) 

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


        

    nptrs = noutput / fanout
    spec = (desc.spec.outcoords, desc.spec.payload)
    
    disk = disk_model_desc(desc, fanin, fanout, density, noutput) 
    index =  index_model_desc(desc, fanin, fanout, density, noutput)
    diskcost = (disk - index) / 1048576.0 / 13.0
    idxcost = 0
    if spec[0] != Spec.COORD_ONE:
        a,b = 3.26775170e-06,   1.39564528e-04
        idxcost = a * index + b * nptrs

    if desc.mode == Mode.PT_MAPFUNC:
        # serialization costs
        sercost = 0.0
        nser = 0.0
        if spec[0] == Spec.BOX:
            nser += 2
        elif spec[0] == Spec.GRID:
            nser += 3 + (math.ceil((fanout/density) ** 0.5) ** 2) - fanout
        else:
            if spec[0] in ( Spec.COORD_MANY, Spec.KEY ):
                nser += 1
            nser += fanin
        nser *= nptrs
        sercost += (1.07e-06 / nser + 4.87e-08) * nser

        
        accost = fanout == 1 and  5.17e-7 or 1.07e-6
        addtocachecosts = nptrs * fanin * accost

        # foo() costs
        if spec[0] == Spec.COORD_ONE:
            foocost = 2.3e-5 * (nptrs + noutput)
        else:
            foocost = 2.3e-5 * (nptrs * 2)
        ovcost = sercost + addtocachecosts + foocost
    else:
        if desc.backward:
            ovcost = overhead_costs(desc, fanin, fanout, density, noutput, nptrs)
        else:
            ovcost = overhead_costs(desc, fanout, fanin, density, noutput / fanout * fanin, nptrs)

    return idxcost + diskcost + ovcost

def overhead_costs(desc, fanin, fanout, density, noutput, nptrs):
    spec = (desc.spec.outcoords, desc.spec.payload)    
    # serialization costs
    sercost = 0.0
    nser = 0.0
    if spec[0] == Spec.BOX:
        nser += 2
    elif spec[0] == Spec.GRID:
        nser += 3 + (math.ceil((fanout/density) ** 0.5) ** 2) - fanout
    else:
        if spec[0] in ( Spec.COORD_MANY, Spec.KEY ):
            nser += 1
        nser += fanin

    if spec[1] == Spec.BOX:
        nser += 2
    elif spec[1] == Spec.GRID:
        nser += 3 + (math.ceil((fanout/density) ** 0.5) ** 2) - fanout
    elif spec[1] == Spec.BINARY:
        nser += 1
    else:
        if spec[1] in ( Spec.COORD_MANY, Spec.KEY ):
            nser += 1
        nser += fanout
    nser *= nptrs
    sercost += (1.07e-06 / nser + 4.87e-08) * nser

    # key overheads
    keycost = 0.0
    # nkeys = 0.0
    # if spec[0] == Spec.KEY:
    #     nkeys += noutput / fanout
    # if spec[1] == Spec.KEY:
    #     nkeys += noutput / fanout
    # keycost = (nkeys) * 0.00012


    # input/output addtocache costs
    addtocachecosts = 0.0
    accost = 0.0
    if fanout == 1:
        accost += 5.17e-7
    else:
        accost += 1.07e-6

    
    if spec[1] == Spec.BOX:
        if fanin == 1:
            accost += 5.17e-7
        else:
            accost += 1.07e-6
    elif spec[1] == Spec.GRID:
        a,b,c = 1.33e-9,   1.96e-06,  -7.19e-10
        if fanin == 1:
            accost += 4.57239151001e-05
        else:
            accost += a * fanin + b / density + c * fanin / density
    addtocachecosts = nptrs * fanin * accost


    # foo() costs
    if spec[0] == Spec.COORD_ONE:
        foocost = 2.3e-5 * (nptrs + noutput)
    else:
        foocost = 2.3e-5 * (nptrs * 2)

    # collision costs

    return sercost + addtocachecosts + keycost  + foocost



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
        return sum(backward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea))
    if desc.mode in (Mode.STAT, Mode.NOOP):
        return 0.0
    if desc.mode == Mode.PT_MAPFUNC:
        return 1000000000

    nentries = noutput / fanout
    if not desc.backward:
        desc = Desc(desc.mode, desc.spec, True)
        return sum(backward_model_desc(desc, fanout, fanin, density, noutput, runtime, nqs, sel, inputarea, nentries))
    return 1000000000
    
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
            cost = sum(backward_model_desc(desc, fanin, fanout, density, noutput,
                                           runtime, nqs, sel, inputarea=inputarea))
            bcost += cost
        scost = min(scost, bcost)
    return scost


def backward_model_desc(desc, fanin, fanout, density, noutput, runtime, nqs, sel, inputarea=None, nentries=None):
    if desc.mode in (Mode.STAT, Mode.NOOP):
        return 0,0,0,0
    if desc.mode == Mode.QUERY:
        return runtime + nqs * sel * 1e-6,0,0,0
    if desc.mode == Mode.FULL_MAPFUNC:
        return nqs * 1e-6,0,0,0

    coord_one = 8.51511955261e-6
    coord_many = 1.36479878426e-6
    box = 4.33921813965e-6
    key = 7.15684890747e-6
    none = 2.87890434265e-6


    if nentries == None:
        nentries = noutput / fanout

        
    boxcost = parsecost = bdbcost = 0
    if not desc.backward:
        return 1000000,0,0,0

    if desc.spec.outcoords == Spec.COORD_ONE:
        idxcost = 0

        a,b =    1.95947313e-10,   3.94434154e-08
        entrysize = 30
        bdbcost = nqs * (1.0 - sel) * 4e-06
        nlookups = nqs * sel
        if desc.spec.payload == Spec.KEY:
            nlookups *= 2
        bdbcost += nlookups * ( nentries * a + entrysize * b / 2 )

        parsecost = nqs * coord_one
        if desc.spec.payload != Spec.BINARY:
            if desc.spec.payload == Spec.GRID:
                negs = (math.ceil((fanin/density) ** 0.5) ** 2) - fanin
                parsecost += nqs * sel * negs * coord_many
            else:
                parsecost += nqs * sel * fanin * coord_many

            if desc.spec.payload == Spec.BOX:
                boxcost = runtime * min((fanin / density), inputarea)  / float(inputarea)
        else:
            parsecost = nqs * sel * 0.000016    
    else:
        #
        a,b,c = 7.4e-09,5.68e-05,1.07e-04
        nmatches = max(nqs * fanout / inputarea, nqs * sel)
        idxcost = a * nentries  + b * nmatches + c


        # bdb lookup costs
        a,b =    1.95947313e-10,   3.94434154e-08
        entrysize = 30
        nlookups = nmatches
        if desc.spec.payload == Spec.KEY:
            nlookups += nmatches
        bdbcost = nlookups * ( nentries * a + entrysize * b / 2 )

        if desc.spec.payload != Spec.BINARY:
            parsecost = 0.0
            parsecost += nmatches * fanout * coord_many
            if desc.spec.payload == Spec.GRID:
                a,b,c = 4.54e-05,   3.78e-07,   2.67e-06
                parsecost += (a / fanin + b / density + c) * fanin
            else:
                parsecost += nqs * sel * fanin * coord_many

            if desc.spec.payload == Spec.BOX:
                boxcost = runtime * min((fanin / density), inputarea)  / float(inputarea)
            # if desc.spec.payload == Spec.BOX:
            #     a,b,c = 0.10013271527592049, 0.00022428472407952143, 1.5845707477154047
            #     baseline = a + b / (density ** c)
            #     perc = min((fanin / density), inputarea) / float(inputarea)
            #     extractcost = (baseline + runtime * perc) * sel
        else:
            parsecost = nmatches * 0.000016
        
    return idxcost, bdbcost, parsecost, boxcost
            

if __name__ == '__main__':
    strats = [Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
              Desc(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False)]
    for desc in strats:
        for fanout in [1, 25, 64]:
            print disk_model_desc(desc, 169, fanout, 1.0, 1000), desc
