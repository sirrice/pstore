import sys
import collections
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
NORMDISK = """select sum(wi.area) * 8.00002 * 2
FROM workflow_inputs as wi, workflow_run as wr
WHERE wi.wid = wr.rowid and wr.eid = ?"""
NORMDISK = """select width * height * 8.00002
FROM exec where rowid = ?"""

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

#stats = Stats.instance('./results/lsstfull.db.nov.14.2011')
#stats = Stats.instance('./results/pablostats.db.nov.14.2011')
#stats = Stats.instance('./_output/pablostats.db')
#stats = Stats.instance('./_output/stats.db')
#stats = Stats.instance('./_output/lsstfull.db')
#stats = Stats.instance('./results/lsstfull.db.nov.10.2011')

if len(sys.argv) <= 1:
    print """python gen_summary.py [ runmode | dbfilename ]*

             dbfilename: default -- ./_output/pablostats.db
          """
    exit()

dbname = './_output/pablostats.db'
runmodes = []
for arg in sys.argv:
    if '.db' in arg:
        dbname = arg
    else:
        try:
            runmode = int(arg)
            runmodes.append(runmode)
        except:
            pass

stats = Stats.instance(dbname)
conn = stats.conn
cur = conn.cursor()


def get_exps(runmode):
    # get all the experiments
    cur.execute("""select rowid, runmode, runtype, width, height, diskconstraint, runconstraint
                from exec where runmode = ? and
                runtype not in ('noop', 'noop_model', 'stats', 'stats_model', 'opt', 'noop_m', 'stats_m')
                order by rowid, diskconstraint""", (runmode,))
    return cur.fetchall()

def get_labels(exps):
    replaces = (('KEY', 'REF'), ('PTR_', '3_'), ('_B', '_b'), ('_F', '_f'), ('PTMAP_', '2_')
                )
    fullreps = (('q_opt', 'Query Log'), ('F_b', 'ONE_REF_b,f'), ('OPT_MANUAL', 'SubZero'),
                ('QUERY', 'BlackBox'), ('ONE_REF', 'FullOne'), ('many_ref', 'FullMany'), ('many_many', 'FullMany'),
        ('queryopt', 'BlackBoxOpt'))
    
    
    
    labels = []
    for rowid, runmode, notes, width, height, dcon, rcon in exps:
        label = "disk(%d)\nrun(%d)" % (dcon, rcon)
        label = 'config-%s' % notes[-1:]
        label = notes.upper()
        for k,v in replaces:
            label = label.replace(k,v)
        for k,v in fullreps:
            if k.lower() == label.lower():
                label = v
        if label.endswith('_M'):
            label = label[:-2]
        if label.startswith('OPT_'):
            arr = label.split('_')
            label = '_'.join(['SubZero', str(int(float(arr[1])))])
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
    OVERHEAD = """SELECT sum(opcost), sum(disk), 0
    FROM pstore_overhead as po, workflow_run as wr, exec
    WHERE po.wid = wr.rowid and wr.eid = exec.rowid and exec.rowid = ? and strat != 's:Func' """

    table = collections.defaultdict(list)
    for rowid, runmode, notes, width, height, dcon, rcon in exps:
        basecost, basedisk = get_baseline(rowid)
        
        # get overhead
        cur.execute(OVERHEAD, (rowid,))
        opcost, disk, idx = cur.fetchone()

        print "opcost", notes, opcost, basecost, opcost/basecost, disk, basedisk, disk/basedisk

        table['Runtime'].append(opcost)
        table['Disk'].append(disk / 1048576.0 + 7.9)
        table['Idx'].append(float(idx) / 1048576.0)


        # table['overhead'].append((opcost - basecost) / basecost)
        # table['disk'].append(float(disk) / basedisk)
        # table['idx'].append(float(idx) / basedisk)
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
        queries[pqid].append((opid, f))
    return queries

def get_qcosts():
    pass

def get_plot(runmode):
    exps = get_exps(runmode)
    labels = get_labels(exps)
    table = get_overhead(exps)
    ymax = max(map(max, table.values()))


    title = "Disk and Runtime Overhead (scaled %dx)" % runmode
    title = "Disk and Runtime Overhead"
    ylabel = "Disk (MB) & Runtime (sec) (log)"
    fname = "overhead%d" % runmode
    ymax = ymax * 10
    ymin = 1
    draw(ymax, ymin, ['Runtime', 'Disk'], table, labels, title,
         ylabel, fname, plotargs={'yscale':'log'})


    # collect all the paths
    allpaths = get_allpaths(exps)
    for path, idx in allpaths.items():
        qname = '%sQ %d' % (path[0][1] and 'F' or 'B', allpaths[path])  
        table[qname] = []


    ymax = 0
    ymin = 10000000000000000
    allcosts = []
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
            qname = '%sQ %d' % (path[0][1] and 'F' or 'B', allpaths[path])            
            if path not in d:
                table[qname].append(0)
                continue
            q = PQCOST % ('(%s)' % ','.join(map(str, d[path])))

            cur.execute(q)
            cost = cur.fetchone()[0]
            #table[allpaths[path]].append( cost)
            table[qname].append(cost)

            print rowid, d[path], float(cost), path
            
            allcosts.append(cost)

    table['FQ 0 Slow'] = table['FQ 5']
    table['FQ 0'] = table['FQ 6']
    features = ['BQ 0', 'BQ 1', 'BQ 2', 'BQ 3', 'BQ 4', 'FQ 0', 'FQ 0 Slow', ]
    #features = sorted(['%sQ %s' % (path[0][1] and 'F' or 'B', pid) for path, pid in allpaths.items()])
    ymax, ymin = max(allcosts), min(allcosts)

    if ymax / ymin > 100:
#    if max(allcosts) > np.std(allcosts) * 5 + np.mean(allcosts):
        yscale = 'log'
        ymax *= 50
        ymin = 0.01
    else:
        yscale = 'linear'
        ymin = 0
        ymax *= 1.2
        #ymax = min(ymax, 30)
    plotargs = { 'yscale' : yscale }
    
    title = 'Query Cost vs Strategies (scaled %dx)' % runmode
    title = 'Query Cost vs Strategies'
    ylabel = 'Query Cost (sec) (log)'
    fname = 'cost%d' % runmode
    draw(ymax, ymin, features, table, labels, title, ylabel, fname, plotargs)
    return
    # get strategies
    cur.execute(STRATS, (rowid,))
    strats = dict(cur.fetchall())
 
    print '%d\t%d\t%d\t%d' % (rowid, runmode, dcon, rcon)
    print '\t%f' * 4 % ( overhead, disk, fcost, bcost)
    for op in strats:
        print '\t%s\t%s' % (op, strats[op])

def draw(ymax, ymin, features, table, labels, title, ylabel, fname, plotargs={}):
    # draw the graph
    fontP = FontProperties()
    fontP.set_size(17)
    figparams = matplotlib.figure.SubplotParams(top=0.9, bottom=0.15)
    fig = plt.figure(figsize=(15, 6), subplotpars=figparams)
    ax = fig.add_subplot(111, ylim=[ymin, ymax*1.2], **plotargs)

    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(19)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(19)



    ind = np.arange(len(table[table.keys()[0]]))#3)
    width = 0.07#0.037
    colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
    n = 0
    rects = []
    for feature in features:
        row = table[feature]
        print row
        rect = ax.bar(ind + (width * n * 1.3) + 0.05, row, width+0.01,
                      color=colors[n % len(colors)], linewidth = 0)
        rects.append(rect)
        n += 1

    # ax.legend(loc='upper center',# 
    #           ncol=3, fancybox=True, shadow=True, prop=fontP)        
    plt.figlegend([rect[0] for rect in rects], ['%s' % f for f in features],
                  loc='upper center', ncol=4, fancybox=True, shadow=True, prop=fontP,
                  bbox_to_anchor=(0.5, .89))
    ax.set_ylabel(ylabel, size=20)
    ax.set_xlabel('Storage Strategies', size=20)
    ax.set_xticks(ind+(width * 5))
    ax.set_xticklabels(labels)
    ax.set_title(title, size=20)
    ax.grid(True)
    plt.savefig('./figs/%s.png' % (fname), format='png')
    plt.cla()
    plt.clf()


for runmode in runmodes:
    get_plot(runmode)

# if len(sys.argv) > 1:
#     for runmode in sys.argv[1:]:
#         runmode = int(runmode)


#get_plot(1)
#get_plot(3)
#get_plot(4)
#get_plot(1)
#get_plot(-1)
#for runmode in range(3, 8):
#    get_plot(runmode)
