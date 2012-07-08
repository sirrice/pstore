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
                'opt' = substr(runtype, 0, 4)
                order by rowid, diskconstraint""", (runmode,))

    experiments = {}  # (disk, runconstraint) -> sorted list of eids
    for rowid, _,_,_,_, disk, run in cur.fetchall():

        key = (disk, run)
        if key not in experiments:
            experiments[key] = []
        experiments[key].append(rowid)
    experiments = dict([ (key, sorted(val)) for key, val in experiments.iteritems() ])
    return experiments

def get_qlogcosts(runmode):
    cur.execute("select rowid from exec where runmode = ? and runtype = 'q_opt'", (runmode,))
    eids = [row[0] for row in cur.fetchall()]
    fcosts, bcosts = get_qcosts(eids)
    print eids
    cur.execute("select sum(opcost) from workflow_run where eid in (%s)" % ','.join(map(str, eids)))
    basecost = cur.fetchone()[0] / float(len(eids))
    return np.mean(fcosts), np.mean(bcosts),  basecost
    


def get_qcosts(eids):
    fcosts, bcosts = [], []
    for eid in eids:
        cur.execute("select avg(cost) from pq where eid = ? and forward = 1", (eid,))
        fcosts.append( cur.fetchone()[0] or 0 )
    for eid in eids:
        cur.execute("select avg(cost) from pq where eid = ? and forward = 0", (eid,))
        bcosts.append( cur.fetchone()[0] or 0 )

    return fcosts, bcosts

def get_overheads(eids):
    q = """SELECT sum(opcost), sum(po.disk), avg(es.disk)
    FROM pstore_overhead as po, workflow_run as wr, exec, exec_stats as es
    WHERE po.wid = wr.rowid and wr.eid = exec.rowid and es.eid = exec.rowid and exec.rowid = ?"""
    disks = []
    opcosts = []
    for eid in eids:
        cur.execute(q, (eid,))
        opcost, disk, totaldisk = tuple(cur.fetchone())
        disks.append(disk / 1048576.0)
        opcosts.append(opcost)
    return disks, opcosts

def get_plot(runmode):
    experiments = get_exps(runmode)

    exps = get_exps(runmode)
    fqlog, bqlog, baseopcost = get_qlogcosts(runmode)
    for (disk, runconstraint), eids in exps.iteritems():
        print '\t', disk, runconstraint, eids
        plot_experiment(disk, runconstraint, eids, fqlog, bqlog, baseopcost)

def plot_experiment(disk, runc, eids, fqlog, bqlog, baseopcost):
    fqcosts, bqcosts = get_qcosts(eids)
    disks, opcosts = get_overheads(eids)
    runs = map(lambda eid: eid-min(eids)+1, eids)
    labels = map(lambda run: 'iter %d' % run, runs)

    fqcosts.insert(0, fqlog)
    bqcosts.insert(0, bqlog)
    disks.insert(0, 0)
    opcosts.insert(0,baseopcost)
    labels.insert(0, 'Blackbox')
    maxlabellen = max(map(len, labels))
    labels = map(lambda l: l.ljust(maxlabellen), labels)


    table = {'FQ' : fqcosts,
             'BQ' : bqcosts,
             'Disk' : disks,
             'Runtime' : opcosts}
    ymax = max(max(disks), max(opcosts))
    #ymax = 40

    title = "Disk and Runtime Overhead (scaled %dx)" % runmode
    ylabel = "Disk (MB) and Runtime (sec)"
    fname = "overhead_%d" % (disk)
    plotargs = {'yscale':'linear'}
    draw(ymax * 1.2, 0, ['Runtime', 'Disk'], table, labels, title,
         ylabel, fname, plotargs=plotargs)


    ymax = max(max(fqcosts), max(bqcosts))
    #ymax = 12
    title = 'Query Cost (scaled %dx)' % runmode
    ylabel = 'Query Cost (sec)'
    fname = 'cost_%d' % (disk)
    draw(ymax, 0, ['FQ', 'BQ'], table, labels, title, ylabel, fname, plotargs)

def draw(ymax, ymin, features, table, labels, title, ylabel, fname, plotargs={}):
    # draw the graph
    fontP = FontProperties()
    fontP.set_size(17)
    figparams = matplotlib.figure.SubplotParams(top=0.9, bottom=0.15)
    fig = plt.figure(figsize=(15, 6), subplotpars=figparams)
    ax = fig.add_subplot(111, ylim=[0.0, ymax*1.2], **plotargs)
    ax.titlesize = 15


    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(19)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(19)



    ind = np.arange(len(table[table.keys()[0]]))#3)
    width = 0.2#0.037
    colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a', '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94', '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d', '#17becf', '#9edae5']
    n = 0
    rects = []
    for feature in features:
        row = table[feature]
        print feature, row
        rect = ax.bar(ind + (width * n * 1.4) + .45, row, width+0.01,
                      color=colors[n % len(colors)], linewidth = 0)
        rects.append(rect)
        n += 1

    # ax.legend(loc='upper center',# 
    #           ncol=3, fancybox=True, shadow=True, prop=fontP)        
    plt.figlegend([rect[0] for rect in rects], ['%s' % f for f in features],
                  loc='upper center', ncol=5, fancybox=True, shadow=True, prop=fontP,
                  bbox_to_anchor=(0.5, .89))
    ax.set_ylabel(ylabel, size=20)
    ax.set_xlabel('Storage Strategies', size=20)
    ax.set_xticks(ind+(width * 5 - 0.3))
    ax.set_xticklabels(labels)
    ax.set_title(title, size=20)
    ax.grid(True)
    plt.savefig('./figs/%s.png' % (fname), format='png')
    plt.cla()
    plt.clf()


for runmode in runmodes:
    get_plot(runmode)
