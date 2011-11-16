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
markers = [None, 'o',]#, '--', '-.', ',', 'o', 'v', '^', '>', '1', '*', ':']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
#markers = [(color, m, l) for color in colors for m in markers for l in linestyles]
markers = [(color, m, l) for m in markers for l in linestyles for color in colors[:4] ]
#markers = [(color, m, l) for m in markers  for color in colors for l in linestyles]


def setup_table(db):
    cur = db.cursor()
    create = """create table stats(id integer primary key autoincrement, strat varchar(128),
                fanin int, fanout int, noutput int,
                wcost float, incache float, outcache float, flush float, serout float, serin float,
                mergecost float, bdbcost float, disk int, idx int)"""
    try:
        cur.execute(create)
        db.commit()
    except Exception, e:
        print e
        db.rollback()

    create = """create table qcosts(sid int references stats(id), qsize int, backward bool,
                cost float, nres int,
                parsecost float, keycost float, datacost float, extractcost float, nhits int)"""
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
    oside = int(math.ceil(fanout ** 0.5))
    prov = []
    incoords = []
    outcoords = []
    fanin_n = 0

    optx, opty = randint(0, 99 - oside), randint(0, 99-oside)  # centroid of output coords
    for n in xrange(noutput):
        #outcoords.append((n / 100, n % 100))
        outcoords.append( (randint(0,oside) + optx, randint(0,side) + opty ) )
        if len(outcoords) >= fanout:
            ptx, pty = randint(0, 99-side), randint(0, 99-side)
            for i in xrange(fanin):
                incoords.append( (randint(0,side) + ptx,randint(0,side) + pty) )
                #incoords.append((fanin_n / 100, fanin_n % 100))
                fanin_n += 1
            prov.append((outcoords, incoords))
            outcoords = []
            incoords = []
            optx, opty = randint(0, 99 - oside), randint(0, 99-oside) 
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


def run_model(db, configs, qsizes):
    from models import *
    setup_table(db)
    cur = db.cursor()
    for config in configs:
        strat, fanin, fanout, noutput = config
        density = 1.0
        wcost = write_model(strat, fanin, fanout, density, noutput, 0)
        disk = disk_model(strat, fanin, fanout, density, noutput)
        idx = index_model(strat, fanin, fanout, density, noutput)
        params = (str(strat), fanin, fanout, noutput, wcost, disk, idx)
        sql = """insert into stats(strat, fanin, fanout, noutput,
                 wcost, disk, idx) values (%s)""" % (','.join(["?"]*(len(params))))
        cur.execute(sql, tuple(params))
        sid = cur.lastrowid

        # query the provenance store
        # vary query size
        for qsize in qsizes:
            fcost = forward_model(strat, fanin, fanout, 1.0, noutput, 0.001, qsize, 1.0, 100*100)
            bcost = backward_model(strat, fanin, fanout, 1.0, noutput, 0.001, qsize, 1.0, 100*100)
            desc = list(strat.descs())[0]
            idxcost, keycost, parsecost, extractcost = backward_model_desc(desc, fanin, fanout, 1.0, noutput,
                                                                           0.001, qsize, 1.0, 100*100)
            params = [sid, qsize, True, bcost, 0]
            # parse, key, data, extract
            params.extend([parsecost, keycost, idxcost, extractcost, 1000])
            sql = "insert into qcosts values (%s)" % ','.join(['?']*len(params))
            cur.execute(sql, tuple(params))

            params = [sid, qsize, False, bcost, 0]
            params.extend([0,0,0,0,0])
            cur.execute(sql, tuple(params))

    db.commit()
    cur.close()
        

def run_exp(db, configs, qsizes):
    class BenchOp(Op):
        class Wrapper(object):
            def __init__(self, arr):
                shape = arr.shape
                self.nargs = 1
                self.shape = shape
                self._arr = arr

            def get_input_shapes(self, run_id):
                return [self.shape]

            def get_input_shape(self, run_id, arridx):
                return self.shape

            def get_inputs(self, run_id):
                return [self._arr]


        def __init__(self, arr):
            super(BenchOp, self).__init__()
            shape = arr.shape
            self.wrapper = BenchOp.Wrapper(arr)
            self.workflow = None
            self.shape = shape

        def run(self, inputs, run_id):
            pass

        def output_shape(self, run_id):
            return self.shape

    arr = np.zeros((100,100))
    op = BenchOp(arr)
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


            incache = pstore.get_stat('incache', 0)
            outcache = pstore.get_stat('outcache', 0)
            flush = pstore.get_stat('flush', 0)
            serin = pstore.get_stat('serin', 0)
            serout = pstore.get_stat('serout', 0)
            mergecost = pstore.get_stat('mergecost', 0)
            bdbcost = pstore.get_stat('bdbcost', 0)
            wcost = pstore.get_stat('write', 0)
            disk = pstore.disk()
            idx = pstore.indexsize()
            runcosts = (wcost, incache, outcache, flush, serout, serin, mergecost, bdbcost, disk, idx)
            costs.append(runcosts)

        #costs = costs[-2:]
        costs = map(np.mean, zip(*costs))
        
        params = [str(strat), fanin, fanout, noutput]
        params.extend(costs)
        sql = """insert into stats(strat, fanin, fanout, noutput,
                 wcost, incache, outcache, flush, serout, serin,
                 mergecost, bdbcost, disk, idx ) values(%s)""" % (','.join(["?"]*(len(params)))) 
        cur.execute(sql, tuple(params))
        sid = cur.lastrowid
        
        

        # query the provenance store
        # vary query size
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

                params = [sid, qsize, backward, cost, nres]
                params.extend([pstore.get_stat('parsecost', 0),
                               pstore.get_stat('keycost', 0),
                               pstore.get_stat('datacost', 0),
                               pstore.get_stat('extractcost', 0),
                               pstore.get_stat('nhits', 0) ])
                
                sql = """insert into qcosts values (%s)""" % (','.join(['?']*len(params)))
                cur.execute(sql, tuple(params))

        db.commit()
    cur.close()


def draw(db, ylabel, fanin, noutput):
    cur = db.cursor()

    where = '1=1'#"strat != 'ONE_MANY_f' and strat != 'ONE_KEY_f' "
    if ylabel == 'disk':
        newylabel = '(disk - ( noutput / fanout ) + (24 + 60.18) + 7340) / 1048576.0'
        newylabel = 'disk / 1048576.0'
    elif ylabel == 'idx':
        newylabel = '( noutput / fanout ) + (24 + 60.18) + 7340'
        newlabel = 'idx / 1048576.0'
    elif ylabel == 'wcost':
        newylabel = 'wcost'
    #newylabel = 'wcost - .0000404 * ( ( noutput / fanout ) + (24 + 60.18) + 7340 )'
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
    print fname
    plot(title, xs, ys, fname, path='_figs/microexp')
    cur.close()
    return




def stacked(db, strat, labels, fanin, noutput=10000):
    cur = db.cursor()
    xs = set()
    ys = {}

    for label in labels:
        sql = """select fanout, avg(%s) from stats, qcosts
                 where stats.id = qcosts.sid and fanin = ? and strat = ? and noutput = ?
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

    
def queries(db, fanin, noutput, strats, qsize=2000, backward=True):
    cur = db.cursor()
    xs = set()
    ys = {}
    print "executing", fanin, noutput, qsize, map(str,strats)
    strats = ','.join(map(lambda s: "'%s'" % s, map(str, strats)))
    sql = """select fanout, strat, backward, cost from stats, qcosts
             where qcosts.sid = stats.id and fanin = ? and noutput = ? and qsize = ? and backward = ?
             and strat in (%s)
             order by fanout"""  % strats
    
    cur.execute(sql, (fanin, noutput, qsize, backward))

    for row in cur.fetchall():
        fanout, strat, backward, cost = row
        label = '%s_%s' % (strat,backward and 'b' or 'f')            
        y = cost
        xs.add(int(fanout))
        if label not in ys: ys[label] = {}
        ys[label][fanout] = float(y)

    title = "noutput=%d     fanin = %s  qsize = %s   backward = %s" % (noutput, fanin, qsize, backward)
    fname = "queries/%s_noutput%s_fanin%s_qsize%s" % (backward and 'b' or 'f', noutput, fanin, qsize)
    plot(title, xs, ys, fname)
    cur.close()



def qstacked(db, strat, labels, fanin, noutput, qsize, backward=True):
    cur = db.cursor()
    xs = set()
    ys = {}

    for label in labels:
        sql = """select fanout, avg(%s) from stats, qcosts
                 where stats.id = qcosts.sid and fanin = ? and strat = ? and noutput = ? and
                 qsize = ? and backward = ?
                group by fanout order by fanout""" % label
        cur.execute(sql, (fanin, strat, noutput, qsize, backward))
#        cur.execute("select fanin, %s from stats where fanout = ? and strat = ? order by fanin" % label, (fanin, strat))
        ys[label] = {}
        for row in cur.fetchall():
            fanout, y = row
            xs.add(int(fanout))
            ys[label][fanout] = float(y)
        print ys[label], strat, label, fanin

    title = "%s     fanin = %s   noutput = %d" % (strat, fanin, noutput)
    fname = "qstack/%s_%s_%s_%s_%s" % ( backward and 'b' or 'f', strat, noutput, fanin, qsize)
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
    ax = fig.add_subplot(111, ylim=[ymin,ymax * mult], yscale = yscale, xscale='log')

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


def fit(db, attr, f, where='1 = 1'):
    print "====%s====" % attr
    cur = db.cursor()
    xs = set()
    ys = {}

    cur.execute('''select strat, fanout, fanin, noutput, %s from
                   stats,qcosts where stats.id = qcosts.sid and %s
                   order by strat, fanout''' % (attr, where))
    for row in cur.fetchall():
        strat, fanout, fanin, noutput, y = row
        if strat not in ys:
            ys[strat] = {}
        ys[strat][(fanout, fanin, noutput)] = y

    # fit for each strategy
    for strat in sorted(ys.keys()):
        xs = []
        vals = []
        for (fanout, fanin, noutput), y in ys[strat].items():
            xs.append((fanout, fanin, noutput))
            vals.append(y)

        xdata = np.array(xs)
        ydata = np.array(vals)
        popt, pcov = curve_fit(f, xdata, ydata)


        print
        print strat
        print np.mean(f(xdata,*popt)), np.std(f(xdata,*popt)), min(f(xdata,*popt)), max(f(xdata,*popt))
        print popt
        print pcov
        print '['
        for a,b,c in zip(xdata, ydata, f(xdata, *popt)):
            print '[', ',\t'.join([',\t'.join(map(str,a)), str(b), str(c)]), ']'
        print ']'



def viz(db, fanins, noutputs):
    labels = ('ser', 'wcost', 'updatecost', 'bdbcost', 'serin', 'serout',
              'serin+serout-bdbcost-ser', 'serout-bdbcost', 'disk')
    labels = ('wcost', 'disk')
    for label in labels:
        for fanin in fanins:
            for noutput in noutputs:
                draw(db, label, fanin, noutput)

def stackviz(db, strats, fanins, noutputs):
    
    labels = ('ser', 'wcost', 'updatecost', 'bdbcost', 'serin','serout',
              'serout-bdbcost', 'disk/1048576.0')
    labels = ('wcost', 'bdbcost', 'serin','serout', 'ser')
    labels = ('wcost', 'incache', 'outcache', 'flush', 'serout', 'serin',
              'mergecost', 'bdbcost', 'flush-bdbcost-serout-serin-mergecost')
    #labels = ('bdbcost', 'disk/1048576.0', 'idx/1048576.0')

    for strat in strats:
        for fanin in fanins:
            for noutput in noutputs:
                try:
                    stacked(db, str(strat), labels, fanin, noutput)
                except Exception ,e:
                    #raise
                    print e, strat, fanin
                    #raise


def queriesstacked(db, fanins, noutputs, strats, qsizes):
    labels = ('cost', 'parsecost', 'keycost', 'datacost', 'extractcost')#, 'nhits')
    for strat in strats:
        for fanin in fanins:
            for noutput in noutputs:
                for qsize in qsizes:
                    try:
                        qstacked(db, str(strat), labels, fanin, noutput, qsize)
                    except Exception ,e:
                        print e, strat, fanin
    

def queriesviz(db, fanins, noutputs, strats, qsizes):
    for fanin in fanins:
        for noutput in noutputs:
            for backward in (True, False):
                for qsize in qsizes:
                    try:
                        queries(db, fanin, noutput, strats, qsize, backward=backward)
                    except Exception, e:
                        print e, fanin, noutput
                        raise



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

    


    strats = [
        #Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.GRID), True),            
        #Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.COORD_MANY), True),
        #Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.BOX), True),
        #Strat.single(Mode.PTR, Spec(Spec.KEY, Spec.KEY), True),
        
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.GRID), True),            
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.BOX), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), True),

        
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.GRID), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), True),            
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.BOX), True),
        Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), True),
        
        #Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_MANY, Spec.BINARY), True),
        #Strat.single(Mode.PT_MAPFUNC, Spec(Spec.COORD_ONE, Spec.BINARY), True),

        # Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_MANY, Spec.KEY), False),
        # Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.KEY), False),
        #Strat.single(Mode.PTR, Spec(Spec.COORD_ONE, Spec.COORD_MANY), False),
    ]

    noutputs = (10000)
    fanins = [1,10,25,50,100]
    fanouts = [1, 10, 100,1000]#10,25,50,100,150,200,250,1000]
    noutputs = (100000,)
    fanins = [1, 10, 100, 500 ]
    fanins = [10]
    fanouts = [1, 50, 100,1000]
    qsizes = [1000]
    #fanins = [1, 10, 100, 1000, 9500, 10000]



    def gen_configs(strats, noutputs, fanins, fanouts):

        for noutput in noutputs:
            for strat in strats:
                for fanin in fanins:
                    for fanout in fanouts:
                        if fanout > noutput:
                            continue
                        yield (strat, fanin, fanout, noutput)
        
    if len(sys.argv) <= 2:
        print "python pstore_micro.py [run | viz | stack | query | fit | all] [db path]"
        exit()
    mode = sys.argv[1]
    dbname = sys.argv[2]
#    db = sqlite3.connect('./_output/pstore_microbench.db')
#    db = sqlite3.connect('./results/pstore_microbench.db.nov.14.2011')
    db = sqlite3.connect(dbname)
    
    if mode in ( 'run', 'all'):    
        run_exp(db, gen_configs(strats, noutputs, fanins, fanouts), qsizes)
    if mode in ( 'model', 'all'):
        run_model(db, gen_configs(strats, noutputs, fanins, fanouts), qsizes)
        
    if mode in ( 'viz', 'all'):            
        viz(db, fanins, noutputs)
    if mode in ( 'stack', 'all'):            
        stackviz(db, strats, fanins, noutputs)
    if mode in ( 'query', 'all'):
        queriesviz(db, fanins, noutputs, strats, qsizes)
    if mode in ( 'qstack', 'all'):
        queriesstacked(db, fanins, noutputs, strats, qsizes)
    if mode in ( 'fit', 'all' ):
        def fgen(f):
            def _f(xs, a,b):
                res = []
                for row in xs:
                    fanout, fanin, noutput = row[0], row[1], row[2]
                    res.append(f(fanout, fanin, noutput, a,b))
                return res
            return _f

        def one_disk(fanout, fanin, noutput, a,b):#, c):#, d):
            nptrs = noutput
            # key:
            key = (18 + 4 * (1 + fanin)) * a
            ptr = (18 + 4) * a
            #return nptrs * ptr + key * noutput / fanout + (nptrs + noutput/fanout) * b

            # many
            #return nptrs * ((4 + 4* (fanin + 1)) * a + b )
        
            #grid
            negs = (int(math.ceil(fanin ** 0.5) ** 2) - fanin)
            negs = 4 * (2 + 1 + negs) 
            return nptrs * ( (4 + negs) * a + b)
            return nptrs * (8 * (a+b))

        def many_disk(fanout, fanin, noutput, a,b):#, d):
            nptrs = noutput / fanout

            # key:
            key = ( 18 + 4 * (1 + fanin) ) 
            insize = 18 
            outsize = (4 + 4*(fanout + 1))            
            ptr = insize + outsize
            idx = 4 + outsize
            #return nptrs * ( ( key + idx + ptr ) * a ) + (nptrs * 3  * b )

            # many
            insize = 4 * (1 + fanin)
            idx = 4 + outsize
            ptr = insize + outsize
            #return nptrs * ( (idx + ptr) * a ) + (nptrs * 2 * b)

            # box:
            insize = 8
            idx = 4 + outsize
            ptr = insize + outsize
            return nptrs * ( (idx + ptr) * a ) + (nptrs * 2 * b)


        def key_disk(fanout, fanin, noutput, a,b):#, d):
            nptrs = noutput / fanout

            # key:
            okey = ( 18 + 4 * (1 + fanout) )
            ikey = ( 18 + 4 * (1 + fanin) )
            ptr = ( 18 + 18 )
            idx = ( 4 + 18 )
            #return nptrs * ( (okey + ikey + ptr + idx) * a ) + (nptrs * 4 * b)

            # box
            ptr = ( 18 + 8 )
            #return nptrs * ( (okey + ptr + idx) * a ) + (nptrs * 3 * b)

            # many:            
            ptr = ( 18 + (4 + 4*(fanin + 1)) )
            return nptrs * ( (okey + ptr + idx) * a ) + (nptrs * 3 * b)


        def idx(fanout, fanin, noutput, a,b):
            return (noutput / fanout) * (24 + a ) + b #* math.log(noutput/fanout)#+ b * fanin + c)

        def bdb(fanout, fanin, noutput, a,b):
            return a * (( noutput / fanout ) + (24 + 60.18) + 7340) + (noutput / fanout) * b
                
        def bcost(fanout, fanin, noutput, a,b):
            return fanin * a + b

    
        #fit(db, 'idx', fgen(idx))
        #fit(db, 'disk - (( noutput / fanout ) + (24 + 60.18) + 7340)', fgen(key_disk))
        # fit(db, 'serin', fgen(f1))
        #fit(db, 'bdbcost', fgen(bdb))
        #fit(db, 'bdbcost', fgen(f2))
        #fit(db, 'updatecost', fgen(f3))
        fit(db, 'cost', fgen(bcost), ' backward = 1 and noutput = 10000 and qsize = 2000')


    db.close()
