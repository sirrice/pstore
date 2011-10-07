import matplotlib.pyplot as plt
import psycopg2 as pg
import Gnuplot as gplot













dims = ['fanin', 'oclustsize', 'density', 'qsize']
labels = 'bulkF,bulkB,setF,setB,gridF,gridB,runF,runB'.split(',')
markers = ['-', '--o']#, '--', '-.', ',', 'o', 'v', '^', '>', '1', '*', ':']
colors = ['b', 'g', 'r', 'c', 'm', 'y']
markers = ['%s%s' % (color, m) for color in colors for m in markers]




def get_strats(conn, run, table, where='1=1'):
    c = conn.cursor()
    c.execute('select distinct strat from %s where rid = %s and %s' %  (table, run,  where))
    strats = map(lambda row: row[0], c.fetchall())
    c.close()
    return strats

def plot(conn, xlabel, ylabel, run, where='1=1', table='stats', plotargs={}):
    global cols
    c = conn.cursor()

    try:
        prefix = './_figs/%s/%s/%s' % (xlabel, table, run)
        os.makedirs(prefix)        
    except:
        pass

    cols = 'fanin,fanout,density,nqs,noutput,opcost'.split(',')
    cols.remove(xlabel)
    strats = get_strats(conn, run, table, where)

    # get all the different charts to make
    q = 'select distinct %s from %s where rid = %s and %s' % (','.join(cols), table, run, where)
    c.execute(q)
    params_all = c.fetchall()

    WHERE = [col + ' = %s' for col in cols]
    WHERE.append('strat = %s')
    if where:
        WHERE.append(where)
    WHERE = ' and '.join(WHERE)


    ymax = 0.0
    c.execute('select max(%s) from %s where rid = %s and %s' % (ylabel, table, run, where))
    print 'select max(%s) from %s where rid = %s and %s' % (ylabel, table, run, where)
    ymax = c.fetchone()[0] * 1.2

    # for each param, generate the graph
    for params in params_all:
        
        # get the points for each strategy (line)
        lines = []
        xs = set()
        for strat in strats:
            _params = list(params)
            _params.append(strat)

            q = 'select %s, %s from %s where rid = %d and %s' % (xlabel,ylabel,table,run,WHERE)
            c.execute(q, tuple(_params))
            pairs = c.fetchall()
            print pairs
            xs.update([pair[0] for pair in pairs])
            lines.append(dict(pairs))

        # special case: if not enough points, probably isn't worth plotting
        if len(xs) < 3:
            continue
        postfix = '_'.join(['='.join(map(str,pair)) for pair in zip(cols, params)])
        fname = "%s_table_%s_run%d_where_%s" % (ylabel, table, run, postfix)

        ymax = max([max(line.values()) for line in lines]) * 1.2

        # draw the graph
        fig = plt.figure()
        ax = fig.add_subplot(111, ylim=[0,ymax], **plotargs)
        xs = sorted(xs)
        for idx, line in enumerate(lines):
            ys = [line.get(pt,0) for pt in xs]
            ax.plot(xs, ys, markers[idx], label=strats[idx])
            #ax.semilogy(xs, ys, markers[idx], label=strats[idx]) 
        #ax.ylim(0, ymax)
        ax.legend()
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.set_title(postfix)
        plt.savefig('%s/%s.png' % (prefix, fname), format='png')
        plt.cla()
        plt.clf()

    c.close()
        
from bench_util import *
conn = get_db()

ptrrun = 68#34 # 16
boxrun = 33 #17
table = 'ptr_plot_abs'
where = "fanin = 100"#"density = 1.0"

for xaxis in ['density']:#['noutput', 'fanin', 'fanout', 'density', 'nqs']:
    #plot(conn, xaxis, 'overhead-(sercost+writecost)', ptrrun,  table=table, where=where)    
    plot(conn, xaxis, 'overhead', ptrrun,  table=table, where=where)
    plot(conn, xaxis, 'disk', ptrrun,  table=table, where=where)
    plot(conn, xaxis, 'fcost', ptrrun,  table=table, where=where)
    plot(conn, xaxis, 'bcost', ptrrun,  table=table, where=where)

# where = "strat in ('s:Box', 's:PtrSet', 's:Bulk')  and opcost = 10"
# for xaxis in [ 'nqs']:
#     plot(conn, xaxis, 'overhead', boxrun,  table='box_plot_abs', where=where)
#     plot(conn, xaxis, 'disk', boxrun,  table='box_plot_abs', where=where)
#     plot(conn, xaxis, 'fcost', boxrun,  table='box_plot_abs', where=where)
#     plot(conn, xaxis, 'bcost', boxrun,  table='box_plot_abs', where=where)


conn.close()
