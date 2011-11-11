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
from random import randint


import numpy as np
from scipy.optimize import curve_fit
from cStringIO import StringIO
import timeit



markers = ['-', '--o']#, '--', '-.', ',', 'o', 'v', '^', '>', '1', '*', ':']
colors = ['b', 'g', 'r', 'c', 'm', 'y']
markers = ['%s%s' % (color, m) for color in colors for m in markers]


def setup_table(db):
    cur = db.cursor()
    create = """create table stats(id integer primary key autoincrement, strat varchar(128),
                fanin int, fanout int, noutput int,
                ser float, wcost float, updatecost float, bdbcost float, serin float, serout float, disk int)"""
    try:
        cur.execute(create)
        db.commit()
    except Exception, e:
        print e
        db.rollback()

    create = "create table qcosts(sid int references stats(id), qsize int, backward bool, cost float, nres int)"
    try:
        cur.execute(create)
        db.commit()
    except  Exception, e:
        print e
        db.rollback()
        
    cur.close()



def get_prov(config):
    strat, fanin, fanout, noutput = config    
    side = int(math.ceil(math.pow(float(fanin), 0.5))) 
    prov = []
    incoords = []
    outcoords = []
    fanin_n = 0
    for n in xrange(noutput):
        outcoords.append((n / 100, n % 100))
        if outcoords >= fanout:
            ptx, pty = randint(0, 99-side), randint(0, 99-side)
            for i in xrange(fanin):
                incoords.append( (randint(0,side) + ptx,randint(0,side) + pty) )
                #incoords.append((fanin_n / 100, fanin_n % 100))
                fanin_n += 1
            prov.append((outcoords, incoords))
            outcoords = []
            incoords = []
    if len(outcoords) > 0:
        ptx, pty = randint(0, 99-side), randint(0, 99-side)
        for i in xrange(fanin):
            incoords.append( (randint(0,side) + ptx,randint(0,side) + pty) )                        
            #incoords.append( (random.randint(0, 99), random.randint(0, 99)) )
            #incoords.append((fanin_n / 100, fanin_n % 100))
            fanin_n += 1
        prov.append((outcoords, incoords))
        outcoords = []

    return prov


def run_exp(db, gen_config):
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
            self.shape = shape

        def run(self, inputs, run_id):
            pass

        def output_shape(self, run_id):
            return self.shape

    op = BenchOp((100,100))
    runid = 1
    setup_table(db)
    cur = db.cursor()    

    for config in gen_configs():
        strat, fanin, fanout, noutput = config
        Runtime.instance().set_strategy(op, strat)
        prov = get_prov(config)

        # data loading costs
        costs = []
        for iterid in xrange(1):
            pstore = op.pstore(runid)
            runid += 1
            for outcoords, incoords in prov:
                if pstore.uses_mode(Mode.PTR):
                    pstore.write(outcoords, incoords)
                else:
                    pstore.write(outcoords, 's' * 10)
            pstore.close()


            updatecost = pstore.get_stat('update_stats', 0)
            bdbcost = pstore.get_stat('bdb', 0)
            serin = pstore.get_stat('serin', 0)
            serout = pstore.get_stat('serout', 0)
            ser = pstore.get_stat('_serialize', 0)
            wcost = pstore.get_stat('write', 0)
            disk = pstore.disk()
            runcosts = (ser, wcost, updatecost, bdbcost, serin, serout, disk)
            costs.append(runcosts)

        #costs = costs[-2:]
        costs = map(np.mean, zip(*costs))

        params = [str(strat), fanin, fanout, noutput]
        params.extend(costs)
        sql = """insert into stats(strat, fanin, fanout, noutput, ser,
                                   wcost, updatecost, bdbcost, serin, serout, disk)
                  values(%s)""" % (','.join(["?"]*(len(params)))) 
        cur.execute(sql, tuple(params))
        sid = cur.lastrowid
        
        

        # query the provenance store
        # vary query size
        qsizes = [1, 10, 100, 1000, 100000]
        qcosts = []
        for qsize in qsizes:
            q = []
            for i in xrange(qsize):
                q.append((random.randint(0, 99), random.randint(0,99)))
            q = Scan(q)
            qcosts.append([])
            for backward in (True, False):
                nres = 0
                start = time.time()
                for res in pstore.join(q, 0, backward):
                    nres += 1
                cost = time.time() - start

                sql = "insert into qcosts values (?, ?, ?, ?, ?)"
                cur.execute(sql, (sid, qsize, backward, cost, nres))

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
#        cur.execute("select fanout, %s from stats where fanin = ? and strat = ? order by fanout" % label, (fanin, strat))
        cur.execute("select fanin, %s from stats where fanout = ? and strat = ? order by fanin" % label, (fanin, strat))
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
    plt.savefig('_figs/microperstrat/%s.png' % '%s_%s' % (strat, fanin), format='png')
    plt.cla()
    plt.clf()
    cur.close()


def fit(db, attr, f):
    print "====%s====" % attr
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
              'serout-bdbcost', 'disk/1048576.0')
    for strat in strats:
        for fanin in fanins:
            try:
                stacked(db, str(strat), labels, fanin)
            except:
                pass



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
    

    def gen_configs():

        strats = [
            Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),
            Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
            Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.GRID), True),

            Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.GRID), True),            
            Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.KEY), True),
            Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.BOX), True),
            Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.KEY), True),
            
            Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.GRID), True),            
            Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
            Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
            Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),

            # Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.BINARY), True),
            # Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True),

            # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
            # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
            # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),
            # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
        ]

        noutput = 100000
        fanins = [1,10,25,50,100]
        fanouts = [1,100,1000]#10,25,50,100,150,200,250,1000]
        fanins = [1, 10, 25]#, 100]#, 200]

        for strat in strats:
            for fanin in fanins:
                for fanout in fanouts:
                    yield (strat, fanin, fanout, noutput)
        
    if len(sys.argv) <= 1:
        print "python pstore_micro.py [run | viz | all]"
        exit()
    mode = sys.argv[1]
    if mode in ( 'run', 'all'):    
        run_exp(db, gen_configs)
    if mode in ( 'viz', 'all'):            
        viz(db, fanins)
    if mode in ( 'stack', 'all'):            
        stackviz(db, strats, fanouts)
    if mode in ( 'fit', 'all' ):
        def fgen(f):
            def _f(xs, a,b,c):
                res = []
                for row in xs:
                    fanout, fanin = row[0], row[1]
                    res.append(f(fanout, fanin, a,b,c))
                return res
            return _f
        def f1(fanout, fanin, a, b, c, d):
            return a * (1000 / fanout) * fanin + b * 1000
        def f2(fanout, fanin, a, b, c, d):
            return a * fanout * fanin + b * 1000
        def f3(fanout, fanin, a, b, c):
            return a * 1000 * fanin / math.pow(fanout, b) + c
            return (1000 / fanout) * (a * fanin + b * fanout ) + c + d
            return a * fanout * fanin + b * fanout + c * fanin + d
        def f4(fanout, fanin, a,b,c):
            # per entry overhead
#            return 1000 / fanout * (a + b * fanin)
            return (1000 / fanout) * (a + fanin * b) + c * 1000
            return (1000 / fanout) * (a ) * fanin
            disk = (fanout * 4.25 + fanout * fanin * 4.25) * (1000 / fanout)
            return ( 1000 / fanout ) * ( fanout * a + fanout * fanin
    +b) + c * disk

    
        fit(db, 'bdbcost', fgen(f4))
        # fit(db, 'serin', fgen(f1))
        # fit(db, 'bdbcost', fgen(f1))
        #fit(db, 'bdbcost', fgen(f2))
        #fit(db, 'updatecost', fgen(f3))        


    db.close()
