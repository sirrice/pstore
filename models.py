import math
from runtime import *

# STRAT_DIFF is hard to predict so I'm punting on it

def disk_model(strat, fanin, fanout, density, noutput, grid=10):
    if strat in [STRAT_Q, STRAT_F, STRAT_NOOP]:
        disk = 0
    elif strat == STRAT_DIFF:
        disk = 100
    elif strat  == STRAT_PSET:
        disk = (7791 * fanin) 
    elif strat ==  STRAT_BOX:
        disk = 45056
    elif strat == STRAT_BULK:
        disk = 140 + (((2 * 6 * 4) + (4 + fanin * 8) + (4 + fanout * 8)))  
        disk *= math.ceil(1000.0 / fanout)
        #disk = 1000 * (8 + fanin * 4) / fanout 
    return disk / 1000.0 * noutput

def write_model(strat, fanin, fanout, density, noutput, grid=10):
    disk = disk_model(strat, fanin, fanout, density, noutput, grid=grid) / 1048576.0
    if strat in [STRAT_Q, STRAT_F]:
        return 0
    if strat == STRAT_DIFF:
        return write_model(STRAT_PSET, fanin, fanout, density, noutput) / 100.0
    elif strat == STRAT_NOOP:
        return 0
        f = float('3.73325355e-8') * fanin + float('3.6e-05')
        total = f * (noutput / fanout)
    elif strat  == STRAT_PSET:
        write = disk 
        ser = disk * (4 + 0.28 * fanin)
        cpu = ( 0.001364 + (-0.0000060517 * fanout)) * fanin + 0.1529
        total = write + ser + cpu + write_model(STRAT_NOOP, fanin, fanout, density, noutput)

        total = 0.00101 * fanin + 0.219
        total += write_model(STRAT_NOOP, fanin, fanout, density, noutput)

    elif strat ==  STRAT_BOX:
        return write_model(STRAT_PSET, fanin, fanout, density, noutput)
        write = disk * 0.85
        ser = disk * 4.2236
        cpu = 0.2 + 0.0012 * fanin
        total = write + ser + cpu + 0
    elif strat == STRAT_BULK:
        a,b,c = 0.00012113359323305934, 0.005033453279969142, 1.0275079057658358
        slope = a + b / (fanout ** c)
        a,b,c = 0.2421166069903642, 0.14210971262225633, 0.43770876881419485
        fixed = a + b / (fanout ** c)
        total = fixed + slope * fanin
        
        # a,b,c = 0.0011887, -0.004104377, 0.5156220
        # total = a * fanin + b * fanout + c
        #total += write_model(STRAT_NOOP, fanin, fanout, density, noutput)
    return total / 1000.0 * noutput + float('1.1e-6') * noutput

def forward_model(strat, fanin, fanout, density, noutput, runtime, nqs, sel,
                  inputarea=None, grid=10):

    if strat in [STRAT_F, STRAT_NOOP]:
        return 0
    if strat == STRAT_Q:
        return runtime
    if strat == STRAT_DIFF:
        return forward_model(STRAT_PSET, fanin, fanout, density, noutput,
                              runtime, nqs, sel, inputarea=inputarea) / 100.0
    elif strat == STRAT_PSET:
        fixed = float('6.4664e-05') * fanout + 0.031577765
        total = fixed + 0.000162287269319807 * fanin        
    elif strat == STRAT_BULK:
        a,b,c = 0.00174984669, 0.03049307, 1.3506390
        fixed = a + b / fanout ** c
        a,b,c = float('1.90661e-06'), float('6.73957315e-05'), 0.722754
        slope = a + b / fanout ** c
        total = fixed + slope * fanin        
    elif strat == STRAT_BOX:
        a,b,c = 17.02816440, 1.45883027, 1.30902734
        bbox = ((math.ceil(fanin**0.5) ** 2) / density)
        bperq = (1 + a * (bbox / (1000*1000)))
        mincost = b * bbox / (1000*1000) * runtime
        total = bperq * nqs * mincost + c

    return total / (1000*1000) * noutput * nqs


def backward_model(strat, fanin, fanout, density, noutput, runtime, nqs, sel,
                   inputarea=None, grid=10):
    if strat in [STRAT_F, STRAT_NOOP]:
        return 0
    if strat == STRAT_Q:
        return runtime
    if strat == STRAT_DIFF:
        return backward_model(STRAT_PSET, fanin, fanout, density, noutput,
                               runtime, nqs, sel, inputarea=inputarea) / 100.0
    elif strat == STRAT_PSET:
        total = 0.000433349 + float('2.27715287889753e-05') * fanin
    elif strat == STRAT_BULK:
        a,b,c = 0.000712780, 0.0327779, 1.042021
        fixed = a + b / fanout ** c
        total = fixed + float('1.72568928627741e-05') * fanin
    elif strat == STRAT_BOX:
        a,b = 0.8101029336838939, 0.5230569678037271
        f = a / (density ** b)
        mincost = (int(math.ceil((float(fanin)**0.5)))**2 / density) * runtime
        total = f * mincost

    return total / (1000*1000) * noutput * nqs
