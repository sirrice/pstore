import sys
sys.path.append('..')
from stats import Stats
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib
import numpy as np


OVERHEADPEROP = """SELECT op, sum(save), sum(overhead) / sum(wr.opcost),
    sum(serialize) / sum(wr.opcost), 
    sum(disk) / 32000080.0, sum(wr.opcost)
FROM pstore_overhead as po, workflow_run as wr, exec 
WHERE po.wid = wr.rowid and wr.eid = exec.rowid and exec.rowid = ? and strat != 's:Func'
GROUP BY op;"""
OPCOST = """SELECT sum(wr.opcost)
FROM workflow_run as wr, exec as e1, exec as e2
WHERE e1.rowid = wr.eid and e1.runtype = 'noop' and
      strat != 's:Func' and e2.rowid = ? and e1.runmode = e2.runmode;"""
NORMDISK = """select sum(wi.area) * 8.00002
FROM workflow_inputs as wi, workflow_run as wr
WHERE wi.wid = wr.rowid and wr.eid = ?"""

PQCOST = """SELECT avg(pq.cost)
     FROM pq
     WHERE rowid in %s"""
STRATS = """SELECT op, strat FROM workflow_run where eid = ? and strat != 's:Func';"""
ALLPATHS = """SELECT pqid, opid, 0, 0 FROM iq, pq WHERE pq.eid = 15 and iq.pqid = pq.rowid
              ORDER BY pqid, iq.rowid, opid, ninputs"""
ALLPATHS = "SELECT rowid, (rowid-1) % 7, 0, 0 FROM pq where pq.eid = ? order by rowid, (rowid-1)%7"
ALLPATHS = """SELECT pq.rowid, pq.forward, pp.opid
              FROM pq, pq_path as pp
              WHERE pp.pqid = pq.rowid and pq.eid = ?
              ORDER BY pq.forward, pq.rowid, pp.idx;"""

stats = Stats.instance('./results/lsstfull.db.nov.10.2011')
#stats = Stats.instance('./_output/pablostats.db')
#stats = Stats.instance('./_output/stats.db')
#stats = Stats.instance('./_output/lsstfull.db')
conn = stats.conn
cur = conn.cursor()


def get_exps(runmode):
    # get all the experiments
    cur.execute("""select rowid, runmode, runtype, width, height, diskconstraint, runconstraint
                from exec where runmode = ? and runtype not in ('noop', 'noop_model', 'stats', 'stats_model', 'opt')
                order by rowid, diskconstraint""", (runmode,))
    return cur.fetchall()

def get_labels(exps):
    replaces = (('KEY', 'REF'), ('PTR_', ''), ('_B', '_b'), ('_F', '_f'))
    fullreps = (('q_opt', 'Query Log'), ('F_b', 'ONE_REF_b,f'))
    
    
    labels = []
    for rowid, runmode, notes, width, height, dcon, rcon in exps:
        label = "disk(%d)\nrun(%d)" % (dcon, rcon)
        label = 'config-%s' % notes[-1:]
        label = notes.upper()
        # for k,v in replaces:
        #     label = label.replace(k,v)
        # for k,v in fullreps:
        #     if k == label:
        #         label = v
        print label
        labels.append(label)
    return labels

def get_baseline(rowid):
    cur.execute(OPCOST, (rowid,))
    opcost = cur.fetchone()[0]
    cur.execute(NORMDISK, (rowid, ))
    normdisk = cur.fetchone()[0]
    return opcost, normdisk

def get_overhead(exps):
    OVERHEAD = """SELECT sum(opcost), sum(disk)
    FROM pstore_overhead as po, workflow_run as wr, exec
    WHERE po.wid = wr.rowid and wr.eid = exec.rowid and exec.rowid = ? and strat != 's:Func' """

    table = {'overhead':[], 'disk':[]}
    for rowid, runmode, notes, width, height, dcon, rcon in exps:
        basecost, basedisk = get_baseline(rowid)
        
        # get overhead
        cur.execute(OVERHEAD, (rowid,))
        opcost, disk = cur.fetchone()

        print "opcost", notes, opcost, basecost, disk, basedisk
        
        table['overhead'].append((opcost - basecost) / basecost)
        table['disk'].append(float(disk) / basedisk)
    return table

def get_allpaths(exps):
    allpaths = set()
    table = dict()
    for rowid, runmode, notes, width, height, dcon, rcon in exps:
        queries = get_queries(rowid)
        paths = map(tuple, queries.values())
        allpaths.update(paths)
    allpaths = dict([(path, idx) for idx, path in enumerate(allpaths)])
    return allpaths

def get_queries(rowid):
    cur.execute(ALLPATHS, (rowid,))
    queries = {}
    for pqid, f, opid in cur:
        if pqid not in queries: queries[pqid] = []
        queries[pqid].append(opid)
    return queries

def get_qcosts():
    pass

def get_plot(runmode):
    exps = get_exps(runmode)
    labels = get_labels(exps)
    table = get_overhead(exps)
    ymax = max(map(max, table.values()))


    draw(ymax * 1.2, ['overhead','disk'], table, labels, 'overhead%d' % runmode, 
         'X times baseline', plotargs={'yscale':'linear'})


    # collect all the paths
    allpaths = get_allpaths(exps)
    for path, idx in allpaths.items():
        table[idx] = []


    ymax = 0
    for rowid, runmode, notes, width, height, dcon, rcon in exps:
        queries = get_queries(rowid)
        paths = map(tuple, queries.values())

        # find all pqids with the same path
        d = {}
        for path in paths:
            d[path] = []
            for k,v in queries.items():
                if path == tuple(v):
                    d[path].append(k)
        
        for path in allpaths:
            if path not in d:
                table[allpaths[path]].append(0)
                continue
            q = PQCOST % ('(%s)' % ','.join(map(str, d[path])))

            cur.execute(q)
            cost = cur.fetchone()[0]
            table[allpaths[path]].append( cost)

            print rowid, d[path], int(cost), path

            
            ymax = max(ymax, cost)

    table = dict([('Q%s' % k, v) for k,v in table.items()])
    features = sorted(['Q%s' % v for v in allpaths.values()])
    #features = sorted(allpaths.values())
    draw(ymax*1.2, features, table, labels, 'costs%d' % runmode,#'%s_qcosts' % title,
         'querycost', plotargs={'yscale':'linear'})
    return
    # get strategies
    cur.execute(STRATS, (rowid,))
    strats = dict(cur.fetchall())
 
    print '%d\t%d\t%d\t%d' % (rowid, runmode, dcon, rcon)
    print '\t%f' * 4 % ( overhead, disk, fcost, bcost)
    for op in strats:
        print '\t%s\t%s' % (op, strats[op])

def draw(ymax, features, table, labels, title, ylabel, plotargs={}):
    # draw the graph
    figparams = matplotlib.figure.SubplotParams(top=0.9)
    fig = plt.figure(figsize=(20, 5), subplotpars=figparams)
    ax = fig.add_subplot(111, ylim=[0.0, ymax*1.2], **plotargs)
    ind = np.arange(len(table[table.keys()[0]]))#3)
    width = 0.05#0.037
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#eeeeee']
    n = 0
    rects = []
    for feature in features:
        row = table[feature]
        print row
        rect = ax.bar(ind + (width * n * 2), row, width, color=colors[n % len(colors)])
        rects.append(rect)
        n += 1

    fontP = FontProperties()
    fontP.set_size('small')
    # ax.legend(loc='upper center',# 
    #           ncol=3, fancybox=True, shadow=True, prop=fontP)        
    plt.figlegend([rect[0] for rect in rects], ['%s' % f for f in features],
                  loc='upper center', ncol=5, fancybox=True, shadow=True, prop=fontP,
                  bbox_to_anchor=(0.5, .93))
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Storage Strategies')
    ax.set_xticks(ind+(width * 5))
    ax.set_xticklabels(labels)
    ax.grid(True)
    plt.savefig('./figs/%s.png' % (title), format='png')
    plt.cla()
    plt.clf()


if len(sys.argv) > 1:
    runmode = int(sys.argv[1])
    get_plot(runmode)

#get_plot(1)
#get_plot(3)
#get_plot(4)
#get_plot(1)
#get_plot(-1)
#for runmode in range(3, 8):
#    get_plot(runmode)
