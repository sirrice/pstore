from runtime import *
from datetime import datetime
from operator import mul
import numpy, sqlite3, logging

slog = logging.getLogger('stats')
logging.basicConfig()
slog.setLevel(logging.ERROR)



class Stats(object):
    """
    prov query performance
     - query (forward, starting, steps, runid)
     - cost

    per op prov query
     - provqid, op, runid, strat, cost

    runtime performance
     - op cost
     - strat
     - prov overhead
     - prov save cost
     - disk space cost
    """

    def __init__(self, fname):
        self.timeseries = False
        self.conn = sqlite3.connect(fname)
        self.setup()


    _instance = None
    @staticmethod
    def instance(fname=None):
        if Stats._instance is None:
            Stats._instance = Stats(fname)
        return Stats._instance



    def setup(self):
        queries = ["""create table exec(id integer primary key autoincrement, 
              time text, width int, height int,  
              runmode int, runtype varchar(20), logdir varchar(128), 
              finished int, notes text, 
              diskconstraint float, runconstraint float, 
              modelparamseid int references exec(id))""",
              """create table exec_stats(eid int references exec(id),
                 disk int)""",
        """create table opt_mappings (eid int references exec(id), 
                       opid int, op varchar(20), strat varchar(20))""",

                   # there are three variables
                   # oclustsize: the size of an output cluster (for pwrite_bulk)
                   # noutcells:  the number of output cells that have provenance attached
                   # outputsize: the total # cells in output array
                   # outputdisk: the total on disk size of output array
        """create table workflow_run (eid int references exec(id), 
                 runid int, opid int, op varchar(20),
                 strat varchar(20), nptrs int,
                 oclustsize int, noutcells int, outputsize int, opcost float, outputdisk int)""",
        """create table workflow_inputs (wid int references workflow_run(rowid),
                 arridx int,  area int)""",

        """create table pstore_overhead (wid int references workflow_run(rowid), 
                    save float, overhead float, serialize float, disk float, idx float)""",
        """create table pstore_stats (wid int references workflow_run(rowid),
                 arridx int, fanin float, area float, density float)""",

        """create table pq (eid int references exec(id), runid int,
                 forward int, cost float, maxres int)""",
        """create table pq_inputs (pqid int references pq(rowid),
                 arridx, ncells)""",
        """create table pq_path (pqid int references pq(rowid), idx int,
                  opid int, op varchar(128))""",
        """create table iq (pqid int references pq(rowid),
                            opid int, op varchar(128), arridx int, strat varchar(20),
                            noutputs int, ninputs int, forward int, cost float)"""
                   ]

        cur = self.conn.cursor()
        for q in queries:
            try:
                cur.execute(q)            
                self.conn.commit()
            except Exception, err:
                slog.error( err )
                self.conn.rollback()
        cur.close()
        

    
    def close(self):
        self.conn.close()

    def finish_exec(self):
        cur = self.conn.cursor()
        cur.execute("update exec set finished = 1 where id = ?", (self.eid,))
        self.conn.commit()
        cur.close


    def add_exec(self, width, height, runmode, runtype, logdir, notes, disk=0, run=0, noopeid=-1):
        cur = self.conn.cursor()
        cur.execute("""insert into exec(time,width,height,runmode,runtype,logdir,notes,
                                        diskconstraint,runconstraint,modelparamseid)
                    values(?,?,?,?,?,?,?,?,?,?)""",
                    (str(datetime.now()), width, height, runmode, runtype, logdir,notes,
                     disk, run, noopeid))
        self.eid = cur.lastrowid
        self.conn.commit()
        cur.close()

    

    def add_opt_mappings(self, strategies):
        cur = self.conn.cursor()
        for op, s in strategies.items():
            cur.execute("insert into opt_mappings values (?,?,?,?)",
                        (self.eid, op.oid, str(op).strip(), str(s)))
        self.conn.commit()
        cur.close()



    def clear(self):
        pass


    def save(self, fname):
        try:
            f = open(fname, 'wb')
            pickle.dump(self.__dict__, f, -1)
            f.close()
        except Exception as e:
            slog.error( "save_stats error\t%s", e )
            raise e
        
    def load(self, fname):
        try:
            f = open(fname, 'rb')
            self.__dict__.update(pickle.load(f))
            f.close()
        except Exception as e:
            slog.error( "load_stats error\t%s", e )
            raise e


    def add_modelrun(self, runid, op, strat, disk, overhead, eids):
        cur = self.conn.cursor()        

        for arridx in xrange(op.wrapper.nargs):
            fanin, area, density, oclustsize, nptrs, noutcells, outputsize, wiarea, opcost = self.get_op_stats(eids, op.oid, arridx)

            if arridx == 0:
                cur.execute("insert into workflow_run values (?,?,?,?,?,?,?,?,?,?,?)",
                            (self.eid, runid, op.oid, str(op).strip(), str(strat),
                             nptrs, oclustsize, noutcells, outputsize, opcost + overhead, 0))
                wid = cur.lastrowid
                cur.execute("insert into pstore_overhead values (?,?,?,?,?,?)",
                            (wid, 0, overhead, 0, disk, 0))

            cur.execute("insert into pstore_stats values (%d,%d,%f,%d,%f)" %
                        (wid, arridx, fanin, area, density))
            cur.execute("insert into workflow_inputs values (?, ?, ?)",
                        (wid, arridx, wiarea))
        self.conn.commit()
        cur.close()
        return wid

    def add_exec_stats(self, disk):
        cur = self.conn.cursor()
        cur.execute("insert into exec_stats values (?,?) ", (self.eid, disk))
        self.conn.commit()
        cur.close()
        
    
    def add_wrun(self, runid, op, cost, inputs, output, outputdisk, pstore):
        strat = str(pstore.strat)
        overhead = pstore.stats.get('write', 0)
        save = pstore.stats.get('write', 0)
        serialize = pstore.stats.get('_serialize', 0)
        nptrs = pstore.get_nptrs()
        oclustsize = pstore.get_oclustsize()    # cluster size
        noutcells = pstore.noutcells # number of outcells that have provenance
        outputsize = reduce(mul,output)       # size of the output array
        disk = pstore.disk()
        fanins = pstore.get_fanins()
        areas = pstore.get_inareas()
        densities = pstore.get_densities()
        idx = pstore.indexsize()

        cur = self.conn.cursor()

        cur.execute("insert into workflow_run values (?,?,?,?,?,?,?,?,?,?,?)",
                    (self.eid, runid, op.oid, str(op).strip(), strat,
                     nptrs, oclustsize, noutcells, outputsize, cost, outputdisk))
        wid = cur.lastrowid
        cur.execute("insert into pstore_overhead values (?,?,?,?,?,?)",
                    (wid, save, overhead, serialize, disk, idx))


        slog.info( "addwrun", op, fanins, areas, densities, inputs, noutcells )
        for arridx, (fanin, outarea, density, input) in enumerate(zip(fanins,
                                                                   areas,
                                                                   densities,
                                                                   inputs)):
            
            cur.execute("insert into pstore_stats values (%d,%d,%f,%d,%f)" %
                        (wid, arridx, fanin, outarea, density))
            inarea = reduce(mul,input)
            cur.execute("insert into workflow_inputs values (?, ?, ?)",
                        (wid, arridx, inarea))


        self.conn.commit()
        cur.close()
        return wid
        

    def add_pq(self, run_id, path, direction, inputs, cost=-1, maxres=-1):
        if not isinstance(path, list): path = [path]
        bforward = direction == 'forward' and 1 or 0

        cur = self.conn.cursor()
        cur.execute("insert into pq values (?,?,?,?,?)", (self.eid, run_id, bforward, cost, maxres))
        pqid = cur.lastrowid

        for pathid, op in enumerate(path):
            cur.execute("insert into pq_path values (?,?,?,?)", (pqid, pathid,
                                                                 op.oid, str(op).strip()))

        for arridx, ncells in enumerate(inputs):
            cur.execute("insert into pq_inputs values (?, ?, ?)", (pqid, arridx, ncells))
        
        self.conn.commit()
        cur.close()
        return pqid

    def add_pq_cost(self, pqid, cost):
        cur = self.conn.cursor()
        cur.execute("update pq set cost = ? where rowid = ?", (cost, pqid))
        self.conn.commit()
        cur.close()

    def add_iq(self, pqid, op, arridx, strat, direction, ninputs, noutputs, cost):
        if direction == 'forward':
            bforward = 1
        else:
            bforward = 0
        cur = self.conn.cursor()
        cur.execute("insert into iq values (?,?,?,?,?,?,?,?,?)",
                    (pqid, op.oid, str(op).strip(), arridx, str(strat).strip(),
                     noutputs,ninputs,bforward,cost))
        iqid = cur.lastrowid
        self.conn.commit()
        cur.close()
        return iqid


    def get_matching_noops(self, runmode, shape, disk=None, runc=None):
        cur = self.conn.cursor()

        if disk is None :
            res = cur.execute("""select rowid from exec where runmode = ? and
            width = ? and height = ? and finished = 1 and runtype = 'stats'""",
            (runmode, shape[0], shape[1]))
        else:
            res = cur.execute("""select rowid from exec where runmode = ? and
            width = ? and height = ? and finished = 1 and runtype = 'stats' and diskconstraint = ?
            and runconstraint = ?""",
            (runmode, shape[0], shape[1], disk, runc))


        eids = [int(row[0]) for row in res]
        slog.info ('get_noops\t%s\t%s\t%s', runmode, shape, eids )
        cur.close()
        return eids

    def get_similar_eids(self, noopeid, diskc, runc):
        cur = self.conn.cursor()
        if not self.timeseries:
            q = """select e2.rowid from exec as e1, exec as e2
            where e1.rowid = ? and e2.width = e1.width and e2.height = e2.height and e2.runmode = e1.runmode
            and e2.finished = 1 and e2.runtype != 'noop' and e2.runtype != 'stats' and e2.finished = 1; """
            res = cur.execute(q, (noopeid,))
        elif diskc == None or runc == None:
            q = """select e2.rowid from exec as e1, exec as e2
            where e1.rowid = ? and e2.width = e1.width and e2.height = e2.height and e2.runmode = e1.runmode
            and e2.finished = 1 and 'opt' = substr(e2.runtype, 0, 4) and e2.finished = 1 ; """
            res = cur.execute(q, (noopeid,))
        else:
            q = """select e2.rowid from exec as e1, exec as e2
            where e1.rowid = ? and e2.width = e1.width and e2.height = e2.height and e2.runmode = e1.runmode
            and e2.finished = 1 and 'opt' = substr(e2.runtype, 0, 4) and e2.finished = 1
            and e2.diskconstraint = ? and e2.runconstraint = ?; """
            slog.debug( "SimilarEIDS\t%d\t%f\t%f", noopeid, diskc, runc)
            res = cur.execute(q, (noopeid, diskc,  runc))
        eids = [int(row[0]) for row in res]

        cur.close()
        return eids
        

        
    def get_pstore_stats(self, eids):
        eids = '(%s)' % ','.join(map(str,eids))
        query = """select wr.eid, wr.opid, wr.op, arridx, avg(fanin), avg(area),
        avg(density), avg(nptrs), avg(oclustsize), avg(noutcells), avg(opcost)
        from pstore_stats as ps, workflow_run as wr
        where ps.wid = wr.oid and wr.eid in %s and strat != 's:Func'
        group by wr.opid, arridx;""" % eids

        cur = self.conn.cursor()
        res = cur.execute(query)
        ret = [tuple(list(row)) for row in res]
        cur.close()
        return ret

    def est_avg_input_b(self, eid, opid, qsize):
        qsizes = []
        for pqid, path in getpaths(eid):
            if opid not in path: continue            
            fanins = [getfanin(op) for op in path]

            for op, fanin in zip(path, fanins):
                if opid == op: break
                qsize *= fanin
            qsizes.append(qsize)
        return qsizes
                
    def get_op_stats(self, eids, opid, arridx):
        if not isinstance(eids, list): eids = [eids]
        eids = '(%s)' % (','.join(map(str,eids)))
        q = """select avg(fanin), avg(ps.area), avg(density), avg(oclustsize),
               avg(nptrs), avg(noutcells), avg(outputsize), avg(wi.area), avg(opcost)
        from workflow_run as wr, workflow_inputs as wi, pstore_stats as ps
        where wr.rowid = wi.wid and ps.wid = wr.rowid and ps.arridx = wi.arridx and
        wr.eid in %s and wr.opid = %s and wi.arridx = %s """ % (eids, opid, arridx)

        cur = self.conn.cursor()
        cur.execute(q)
        ret = cur.fetchone()
        cur.close()

        if not ret or  ret[0] == None: return None
        return tuple(ret)
        
        

    def get_array_sizes(self, eid, opid, arridx, forward):
        q = """select avg(wi.area), avg(outputsize)
        from workflow_run as wr, workflow_inputs as wi
        where wr.rowid = wi.wid and wr.eid = ? and wr.opid = ? and wi.arridx = ?"""
        cur = self.conn.cursor()
        cur.execute(q, (eid, opid, arridx))
        ret = cur.fetchone()
        if ret is None: return None
        incount, outcount = ret[0], ret[1]
        cur.close()
        return incount, outcount


    def get_pq_bstats(self, eid, opid, arridx):

        q = """select iq.arridx, avg(iq.ncells), avg(noutputs)
        from iq, pq
        where iq.pqid = pq.rowid and iq.forward = 0 and
        pq.eid = ? and  iq.opid = ? and iq.arridx = ?"""
        cur = self.conn.cursor()
        res = cur.execute(q, (eid, opid, arridx))
        ret = [tuple(row) for row in res]
        cur.close()
        return ret

    def get_iq_stat(self, eids, opid, arridx):
        eids = sorted(eids, reverse=True)
        cur = self.conn.cursor()
        fqsizes, bqsizes = [], []
        for eid in eids:
            q = """select avg(cast(noutputs as real) / (ninputs+1)), count(*), avg(ninputs)
            from iq, pq, exec where iq.pqid = pq.rowid and exec.rowid = pq.eid and iq.opid = ?
            and iq.arridx = ? and pq.eid = ? and iq.forward = 1  and exec.finished = 1
            """
            cur.execute(q, (opid, arridx, eid))
            row = cur.fetchone()
            fqsizes.append(row)

            q = """select avg(cast(noutputs as real) / ninputs), count(*), avg(ninputs)
            from iq, pq, exec where iq.pqid = pq.rowid and exec.rowid = pq.eid and iq.opid = ?
            and iq.arridx = ? and pq.eid = ? and iq.forward = 0 and exec.finished = 1
            """
            cur.execute(q, (opid, arridx, eid))
            row = cur.fetchone()
            bqsizes.append(row)

        cur.close()
        return fqsizes, bqsizes
        

    def get_disk(self, runid, op, strat, default=0.0):
        q = """SELECT disk / 1048576.0
        FROM workflow_run as wr, pstore_overhead as po
        WHERE po.wid = wr.rowid and wr.runid = ? and wr.op = ? and wr.strat = ?"""
        cur = self.conn.cursor()
        res = cur.execute(q, (runid, str(op).strip(), str(strat)))
        for row in res:
            return row[0]
        return default


    def get_overhead(self, runid, op, strat, default=0.0):
        q = """SELECT overhead
        FROM workflow_run as wr, pstore_overhead as po
        WHERE po.wid = wr.rowid and wr.runid = ? and wr.op = ? and wr.strat = ?"""
        cur = self.conn.cursor()
        res = cur.execute(q, (runid, str(op).strip(), str(strat)))
        for row in res:
            return row[0]
        return default





if __name__ == '__main__':
    import sys
    from models import *
    from runtime import *
    Stats.instance('outputs/test//stats.db')

    def f(name, eids, opid, arridx):
        stats = Stats.instance().get_op_stats(eids, opid, arridx)
        fanin, area, density, oclustsize, nptrs, noutcells, outputsize, inputsize, opcost = stats
        strats = [STRAT_BOX, STRAT_SUPERBOX, STRAT_Q, STRAT_PSET, STRAT_BULK, STRAT_RUN]
        fqsize = bqsize = 10
        for strat in strats:
            prov = write_model(strat, fanin, oclustsize, density, noutcells) 
            disk = disk_model(strat, fanin, oclustsize, density, noutcells)
            fcost = forward_model(strat, fanin, oclustsize, density, noutcells,
                                  opcost, fqsize, 1.0, inputarea=inputsize)
            bcost = backward_model(strat, fanin, oclustsize, density, noutcells,
                                   opcost, bqsize, 1.0, inputarea=inputsize)

            print "%s\t%s\t%f\t%f\t%10d\t%f\t%f" % (name, strat, density, prov, disk, fcost, bcost)
        
        
    eid = int(sys.argv[1])
    f("removecrs(0)", [eid], 26, 0)
    f("removecrs(1)", [eid], 26, 1)
    f("cluster     ", [eid], 27, 0)
    f("detectA     ", [eid], 11, 0)
    f("detectB     ", [eid], 23, 0)
    

