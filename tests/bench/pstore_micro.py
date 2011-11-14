import sys, random
sys.path.insert(0, "../")
sys.path.insert(0, "../../")

#from bench_util import *
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


linestyles = ['-', '--']
markers = [None, 'o', '*']#, '--', '-.', ',', 'o', 'v', '^', '>', '1', '*', ':']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
markers = [(color, m, l) for color in colors for m in markers for l in linestyles]
#markers = [(color, m, l) for m in markers for l in linestyles for color in colors ]


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
        if len(outcoords) >= fanout:
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


def run_exp(db, configs):
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

    for config in configs:
        strat, fanin, fanout, noutput = config
        Runtime.instance().set_strategy(op, strat)
        prov = get_prov(config)
        print 'running %s\t%d\t%d\t%d' % tuple(config)

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
        qsizes = [1, 10, 100, 1000, 10000]
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


def draw(db, ylabel, fanin, noutput):
    cur = db.cursor()

    where = "strat != 'ONE_MANY_f' and strat != 'ONE_KEY_f' "
    if ylabel == 'disk':
        newylabel = 'disk / 1048576.0'
    else:
        newylabel = ylabel

    sql = "select strat, fanout, %s from stats where fanin = ? and noutput = ? and %s order by strat"
    sql = sql % (newylabel, where)
    cur.execute(sql, (fanin, noutput))
    xs = set()
    ys = {}
    for row in cur.fetchall():
        strat, fanout, y = row
        xs.add(int(fanout))
        vals = ys.get(strat, {})
        vals[fanout] = float(y)
        ys[strat] = vals

    title = "%s      fanin = %s noutput = %d" % (ylabel, fanin, noutput)
    fname = '%s_%s_%s' % (ylabel, fanin, noutput)
    plot(title, xs, ys, fname, path='_figs/microexp')
    cur.close()
    return




def stacked(db, strat, labels, fanin, noutput=10000):
    cur = db.cursor()
    xs = set()
    ys = {}

    for label in labels:
        sql = """select fanout, avg(%s) from stats
                 where fanin = ? and strat = ? and noutput = ?
                group by fanout order by fanout""" % label
        cur.execute(sql, (fanin, strat, noutput))
#        cur.execute("select fanin, %s from stats where fanout = ? and strat = ? order by fanin" % label, (fanin, strat))
        ys[label] = {}
        for row in cur.fetchall():
            fanout, y = row
            xs.add(int(fanout))
            ys[label][fanout] = float(y)
        print ys[label], strat, label, fanin

    title = "%s     fanin = %s   noutput = %d" % (strat, fanin, noutput)
    fname = "%s_%s_%s" % ( noutput, strat, fanin)
    plot(title, xs, ys, fname)
    cur.close()

    
def queries(db, fanin, noutput, qsize=10000):
    cur = db.cursor()
    xs = set()
    ys = {}
    print "executing", fanin, noutput
    sql = """select fanout, strat, backward, cost from stats, qcosts
             where qcosts.sid = stats.id and fanin = ? and noutput = ? and qsize = ? order by fanout""" 

    cur.execute(sql, (fanin, noutput, qsize))

    for row in cur.fetchall():

        fanout, strat, backward, cost = row
        #if qsize >= 1000: qsize = '%dk' % (qsize / 1000)
        #label = 'qsize%s_%s' % (qsize,backward and 'b' or 'f')
        label = '%s_%s' % (strat,backward and 'b' or 'f')            
        y = cost
        xs.add(int(fanout))
        if label not in ys: ys[label] = {}
        ys[label][fanout] = float(y)

    title = "noutput=%d     fanin = %s" % (noutput, fanin)
    fname = "noutput%s_fanin%s" % (noutput, fanin)
    plot(title, xs, ys, fname)
    cur.close()

    

def plot(title, xs, ys, fname, path = '_figs/microperstrat'):
    xs = sorted(xs)
    lines = []
    labels = []
    for label in sorted(ys.keys()):
        labels.append(label)
        lines.append((label, [ys[label].get(x, 0) for x in xs]))
    ymax = max(map(max, map(lambda x: x[1], lines)))
    ymin = min(map(min, map(lambda x: x[1], lines)))

    # if ymin == 0:
    #     if ymax > 100:
    #         yscale = 'log'
    #         mult = 100
    #     else:
    #         yscale = 'linear'
    #         mult = 1.5
    # elif ymax / ymin > 50:
    #     yscale = 'log'
    #     ymin = ymin / 10
    #     mult = 100
    # else:
    yscale = 'linear'
    mult = 1.5
    ymin = 0



    # draw the graph
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, ylim=[ymin,ymax * mult], yscale = yscale)

    for idx, (label, ys) in enumerate(lines):
        c,m,l = markers[idx]
        ax.plot(xs, ys, color = c, marker=m, linestyle=l, label=label, linewidth=1.5)

    fontP = FontProperties()
    fontP.set_size('small')
    ax.legend(loc='upper center',# bbox_to_anchor=(0.5, 1.05),
              ncol=3, fancybox=True, shadow=True, prop=fontP)        
#    ax.legend(loc=2)
    ax.set_ylabel('time (sec)')
    ax.set_xlabel('fanout')
    ax.set_title(title)
    plt.savefig('%s/%s.png' % (path,fname), format='png')
    plt.cla()
    plt.clf()


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
    labels = ('wcost', 'bdbcost', 'disk')
    for label in labels:
        for fanin in fanins:
            for noutput in (100, 1000, 10000):
                draw(db, label, fanin, noutput)

def stackviz(db, strats, fanins):
    
    labels = ('ser', 'wcost', 'updatecost', 'bdbcost', 'serin','serout',
              'serout-bdbcost', 'disk/1048576.0')
    labels = ('wcost', 'bdbcost', 'serin','serout', 'ser')
#              'disk/1048576.0')
    for strat in strats:
        for fanin in fanins:
            for noutput in (100, 1000, 10000):
                try:
                    stacked(db, str(strat), labels, fanin, noutput)
                except Exception ,e:
                    #raise
                    print e, strat, fanin


def queriesviz(db, fanins, noutputs):
    for fanin in fanins:
        for noutput in noutputs:
            try:
                queries(db, fanin, noutput)
            except Exception, e:
                print e, fanin, noutput

    


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

    strats = [

        # Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.GRID), True),            
        # Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.COORD_MANY), True),
        # Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.BOX), True),
        # Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.KEY), True),
        
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.GRID), True),            
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),

        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.GRID), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),            


        # Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.BINARY), True),
        # Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True),

        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
    ]

    noutputs = (100,1000,10000)
    fanins = [1,10,25,50,100]
    fanouts = [1, 10, 100,1000]#10,25,50,100,150,200,250,1000]

    noutputs = (10000,)
    fanins = [1, 10, 25]#, 100]#, 200]
    fanins = [10]


    def gen_configs(strats, noutputs, fanins, fanouts):

        for noutput in noutputs:
            for strat in strats:
                for fanin in fanins:
                    for fanout in fanouts:
                        if fanout > noutput:
                            continue
                        yield (strat, fanin, fanout, noutput)
        
    if len(sys.argv) <= 1:
        print "python pstore_micro.py [run | viz | stack | query | fit | all]"
        exit()
    mode = sys.argv[1]
    if mode in ( 'run', 'all'):    
        run_exp(db, gen_configs(strats, noutputs, fanins, fanouts))
    if mode in ( 'viz', 'all'):            
        viz(db, fanins)
    if mode in ( 'stack', 'all'):            
        stackviz(db, strats, fanins)
    if mode in ( 'query', 'all'):
        queriesviz(db, fanins, noutputs)
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
