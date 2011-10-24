import matplotlib.pyplot as plt
import psycopg2 as pg











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

    if table == 'tmp':
        cols = 'fanin,fanout,density,nqs,noutput,writecost'.split(',')
    else:
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
            print q % tuple(_params)
            c.execute(q, tuple(_params))
            pairs = c.fetchall()
            xs.update([pair[0] for pair in pairs])
            lines.append(dict(pairs))
            print pairs

        # special case: if not enough points, probably isn't worth plotting
        if len(xs) < 3:
            print xs
            continue
        postfix = '_'.join(['='.join(map(str,pair)) for pair in zip(cols, params)])
        fname = "%s_table_%s_run%d_where_%s" % (ylabel, table, run, postfix)

        try:
            ymax = max([max(line.values()) for line in lines if len(line)]) * 1.2
        except:
            print lines
            raise
        # draw the graph
        fig = plt.figure()
        ax = fig.add_subplot(111, ylim=[0,ymax], **plotargs)
        xs = sorted(xs)
        for idx, line in enumerate(lines):
            ys = [line.get(pt,0) for pt in xs]
            ax.plot(xs, ys, markers[idx], label=strats[idx])
            #ax.semilogy(xs, ys, markers[idx], label=strats[idx]) 
        #ax.ylim(0, ymax)
        ax.legend(loc=2)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.set_title(postfix)
        plt.savefig('%s/%s.png' % (prefix, fname), format='png')
        plt.cla()
        plt.clf()

    c.close()
        
from bench_util import *
conn = get_db()

#=== attributes ===
# x axis
# y axis
# lines
# files
# where clause (constraints)
# table

import sys
if len(sys.argv) < 2:
    print 'python gen_read_graphs.py run_id table_name'
    exit()
run = int(sys.argv[1])
table = sys.argv[2]

#table = 'box_model_plot_abs'
where = "strat != 'ONE_MANY_f' and strat != 'ONE_MANY_b'"#"1=1"#"fanin = 100"#"density = 1.0"

for xaxis in ['fanin', 'fanout', 'density', 'nqs']:
    #plot(conn, xaxis, 'overhead-(sercost+writecost)', run,  table=table, where=where)    
    plot(conn, xaxis, 'overhead', run,  table=table, where=where)
    plot(conn, xaxis, 'disk', run,  table=table, where=where)
    plot(conn, xaxis, 'fcost', run,  table=table, where=where)
    plot(conn, xaxis, 'bcost', run,  table=table, where=where)
conn.close()
exit()

where = "1=1"#opcost = 6"
for xaxis in [ 'density']:#, 'nqs', 'writecost']:
    plot(conn, xaxis, 'overhead', boxrun,  table=table, where=where)
    plot(conn, xaxis, 'disk', boxrun,  table=table, where=where)
    plot(conn, xaxis, 'fcost', boxrun,  table=table, where=where)
    plot(conn, xaxis, 'bcost', boxrun,  table=table, where=where)


conn.close()
