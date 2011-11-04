import sys, random
sys.path.insert(0, "../")
sys.path.insert(0, "../../")

from bench_util import *
from models import *
from provstore import *
from op import *
from runtime import *
from arraystore import ArrayStore
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.font_manager import FontProperties


import numpy as np
from scipy.optimize import curve_fit
from cStringIO import StringIO
import timeit



markers = ['-', '--o']#, '--', '-.', ',', 'o', 'v', '^', '>', '1', '*', ':']
colors = ['b', 'g', 'r', 'c', 'm', 'y']
markers = ['%s%s' % (color, m) for color in colors for m in markers]


def run_exp(db, strats, fanins, fanouts, noutput):
    class BenchOp(Op):
        class Wrapper(object):
            def __init__(self, shape):
                self.nargs = 1
                self.shape = shape
                self._arr = None

            def get_input_shapes(self, run_id):
                return [self.shape]

            def get_input_shape(self, run_id, arridx):
                return self.shape

            def get_inputs(self, run_id):
                return [self._arr]


        def __init__(self, shape):
            super(BenchOp, self).__init__()
            self.wrapper = BenchOp.Wrapper(shape)
            self.workflow = None

        def run(self, inputs, run_id):
            pass

        def output_shape(self, run_id):
            return (100,100)

    op = BenchOp((100,100))
    
    runid = 1

    cur = db.cursor()
    create = """create table stats(strat varchar(128), fanin int, fanout int, ser float, wcost float,
    updatecost float, bdbcost float, serin float, serout float, disk int)"""
    try:
        cur.execute(create)
        db.commit()
    except:
        db.rollback()
    cur.close()
    cur = db.cursor()
    
    for strat in strats:
        Runtime.instance().set_strategy(op, strat)
        for fanin in fanins:
            for fanout in fanouts:
                prov = []
                incoords = []
                outcoords = []
                fanin_n = 0
                for n in xrange(noutput):
                    outcoords.append((n / 100, n % 100))
                    for i in xrange(fanin): 
                        incoords.append((fanin_n / 100, fanin_n % 100))
                    if n % fanout == 0:
                        prov.append((outcoords, incoords))
                        outcoords = []
                        incoords = []
                if len(outcoords) > 0:
                    prov.append((outcoords, incoords))
                    outcoords = []
                

                costs = []
                for iterid in xrange(2):
                    pstore = op.pstore(runid)
                    runid += 1
                    for outcoords, incoords in prov:
                        pstore.write(outcoords, incoords)
                        updatecost = pstore.get_stat('update_stats', 0)
                        bdbcost = pstore.get_stat('bdb', 0)
                        serin = pstore.get_stat('serin', 0)
                        serout = pstore.get_stat('serout', 0)
                        ser = pstore.get_stat('_serialize', 0)
                        wcost = pstore.get_stat('write', 0)
                        disk = pstore.disk()
                        runcosts = (ser, wcost, updatecost, bdbcost, serin, serout, disk)
                        costs.append(runcosts)
                    pstore.close()
                costs = costs[-2:]
                costs = map(np.mean, zip(*costs))

                params = [str(strat), fanin, fanout]
                params.extend(costs)
                sql = "insert into stats values(%s)" % (','.join(["?"]*(len(params)))) 
                cur.execute(sql, tuple(params))
                db.commit()
    cur.close()


def draw(db, ylabel, fanin):
    cur = db.cursor()
#    if fanin > 1:
    where = "strat != 'ONE_MANY_f' and strat != 'ONE_KEY_f' "
    #else:
     #   where = '1=1'
    cur.execute("select strat, fanout, %s from stats where fanin = ? and %s order by strat" % (ylabel, where), (fanin,))
    xs = set()
    ys = {}
    for row in cur.fetchall():
        strat, fanout, y = row
        xs.add(int(fanout))
        vals = ys.get(strat, {})
        vals[fanout] = float(y)
        ys[strat] = vals

    xs = sorted(xs)
    lines = []
    labels = []
    ymax = 0
    for strat in ys:
        labels.append(strat)
        lines.append((strat, [ys[strat].get(x, 0) for x in xs]))
    ymax = max(map(max, map(lambda x: x[1], lines))) * 1.2

    # draw the graph
    fig = plt.figure()
    ax = fig.add_subplot(111, ylim=[0,ymax])
    for idx, (label, ys) in enumerate(lines):
        ax.plot(xs, ys, markers[idx], label=label)

    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(loc='upper center',# bbox_to_anchor=(0.5, 1.05),
              ncol=3, fancybox=True, shadow=True, prop=fontP)        
#    ax.legend(loc=2)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('fanout')
    ax.set_title("%s      fanin = %s" % (ylabel, fanin))
    plt.savefig('_figs/microexp/%s.png' % '%d_%s' % (fanin, ylabel), format='png')
    plt.cla()
    plt.clf()
    cur.close()



def stacked(db, strat, labels, fanin):
    cur = db.cursor()
    xs = set()
    ys = {}
    
    for label in labels:
        cur.execute("select fanout, %s from stats where fanin = ? and strat = ? order by fanout" % label, (fanin, strat))
        ys[label] = {}
        for row in cur.fetchall():
            fanout, y = row
            xs.add(int(fanout))
            ys[label][fanout] = float(y)

    xs = sorted(xs)
    lines = []
    labels = []
    ymax = 0
    for label in ys:
        labels.append(label)
        lines.append((label, [ys[label].get(x, 0) for x in xs]))
    ymax = max(map(max, map(lambda x: x[1], lines))) * 1.2

    # draw the graph
    fig = plt.figure()
    ax = fig.add_subplot(111, ylim=[0,ymax])
    for idx, (label, ys) in enumerate(lines):
        ax.plot(xs, ys, markers[idx], label=label)

    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(loc='upper center',# bbox_to_anchor=(0.5, 1.05),
              ncol=3, fancybox=True, shadow=True, prop=fontP)        
#    ax.legend(loc=2)
    ax.set_ylabel(strat)
    ax.set_xlabel('fanout')
    ax.set_title("%s      fanin = %s" % (strat, fanin))
    plt.savefig('_figs/microperstrat/%s.png' % '%d_%s' % (fanin, strat), format='png')
    plt.cla()
    plt.clf()
    cur.close()


def fit(db, attr, f):
    cur = db.cursor()
    xs = set()
    ys = {}

    cur.execute('select strat, fanout, fanin, %s from stats where fanout != 200 order by strat, fanout' % attr)
    for row in cur.fetchall():
        strat, fanout, fanin, y = row
        if strat not in ys:
            ys[strat] = {}
        ys[strat][(fanout, fanin)] = y

    # fit for each strategy
    for strat in ys:
        xs = []
        vals = []
        for (fanout, fanin), y in ys[strat].items():
            xs.append((fanout, fanin))
            vals.append(y)

        xdata = np.array(xs)
        ydata = np.array(vals)
        popt, pcov = curve_fit(f, xdata, ydata)
        print
        print strat
        print popt
        print pcov



def viz(db, fanins):
    labels = ('ser', 'wcost', 'updatecost', 'bdbcost', 'serin', 'serout',
              'serin+serout-bdbcost-ser', 'serout-bdbcost', 'disk')
    for label in labels:
        for fanin in fanins:
            draw(db, label, fanin)

def stackviz(db, strats, fanins):
    labels = ('ser', 'wcost', 'updatecost', 'bdbcost', 'serin','serout',
              'serin+serout-bdbcost-ser', 'serout-bdbcost')
    for strat in strats:
        for fanin in fanins:
            stacked(db, str(strat), labels, fanin)



def bbox(coords):
    if not len(coords):
        return ((-1, -1), (-1,-1))
    minc, maxc = None, None
    for coord in coords:
        if minc == None:
            minc = coord
            maxc = coord
        else:
            minc = map(min, zip(minc, coord))
            maxc = map(max, zip(maxc, coord))

    return (minc, maxc)



if __name__ == '__main__':
    import sqlite3

    
    db = sqlite3.connect('./_output/pstore_microbench.db')
    
    def fgen(f):
        def _f(xs, a, b):
            res = []
            for row in xs:
                fanout, fanin = row[0], row[1]
                res.append(f(fanout, fanin, a, b))
            return res
        return _f
    def f1(fanout, fanin, a, b):
        return a * (1000 / fanout) * fanin + b * 1000
    def f2(fanout, fanin, a, b):
        return a * fanout * fanin + b * 1000

    
    #fit(db, 'serin', fgen(f1))
    fit(db, 'bdbcost', fgen(f1))
    fit(db, 'bdbcost', fgen(f2))
    #fit(db, 'updatecost', fgen(f1))        

    exit()

    strats = [
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
        
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),        
        
        ]

    noutput = 1000
    fanins = [1,10,25,50,100,200]
    fanouts = [1,10,25,50,100,150,200,250,1000]
    fanins = [50, 100, 200]
    if len(sys.argv) <= 1:
        print "python pstore_micro.py [run | viz | all]"
        exit()
    mode = sys.argv[1]
    if mode in ( 'run', 'all'):    
        run_exp(db, strats, fanins, fanouts, noutput)
    if mode in ( 'viz', 'all'):            
        viz(db, fanins)
    if mode in ( 'stack', 'all'):            
        stackviz(db, strats, fanins)
    


    db.close()