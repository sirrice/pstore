import sqlite3, time, logging, sys, pickle
from provstore import *
from strat import *


log = logging.getLogger('runtime')
logging.basicConfig()
log.setLevel(logging.DEBUG)

PSTOREDIR = './_pstore'

class Runtime(object):
    """
    provstores separated into writestores and read stores.
    Stores each type separately.  Writes stores are not allowed
    to have corresponding readstores unless closed.

    All open writestores are in-memory
    Manages an in-memory cache of read stores and uses round robin replacement (LRU maybe?)
    """
    _instance = None

    def __init__(self):
        self.cur_strats = {} # op -> strategy
        self.old_strats = {} # op -> {run_id -> strat}
        self.pstores = {}    # (op,run_id) -> pstore
        self.override = {}   # (op,run_id) -> pstore

        # don't _actually_ delete old provstores
        # in experiments, want to restore them later on
        self.rm_strats = {}  # (op,run_id) -> strat
        self.rm_pstores = {}  # (op,run_id) -> pstore


    @staticmethod
    def instance():
        if Runtime._instance is None:
            Runtime._instance = Runtime()
        return Runtime._instance

    def clear(self):
        self.cur_strats = {} # op -> strategy
        self.old_strats = {} # op -> {run_id -> strat}
        self.pstores = {}    # (op,run_id) -> pstore
        self.override = {}   # (op,run_id) -> pstore
        

    def set_strategy(self, op, strat):
        self.cur_strats[op] = strat
        log.debug( "set\t%s\t%s" % (op, strat) )

    def available_strats(self, op, run_id):
        ret = []
        strat = self.get_strategy(op, run_id)
        if Mode.QUERY not in strat.modes():
            ret.append(Strat.query())
        ret.append(strat)
        return ret

    def get_strategy(self, op, run_id = None):
        if run_id == None or op not in self.old_strats or run_id not in self.old_strats[op]:
            return self.cur_strats[op]
        return self.old_strats[op][run_id]

    def check_strategy(self, op, run_id, strat):
        if op not in self.old_strats:
            return False
        if run_id not in self.old_strats[op]:
            return False
        if self.old_strats[op][run_id] != strat:
            return False
        return True

    def get_disk_strategies(self):
        print len(self.old_strats)
        ret = []
        for op, strats in self.old_strats.items():
            for runid, s in strats.items():
                if set(s.modes()).intersection([Mode.PT_MAPFUNC, Mode.PT_MAPFUNC_BOX, Mode.PTR]):
                    ret.append( (runid, op, s) )
        return ret
    

    def get_filename(self, op, run_id):
        return '%s/%s_%d_%s' % (PSTOREDIR,
                                str(op).strip(),
                                run_id,
                                self.get_strategy(op, run_id).fnamestring())

    def create_pstore(self, op, run_id, fname, strat):
        """
        Given the strategy, composes the proper provenance storage backend
        Considers
        1) Forward vs Backward (strat.backward)
        2) Storage class 1, 2 vs 3 (strat.mode)
        3) Single vs Composite (strat.mode)
        4) Specific storage implementation
        """

        ret = []
        for buckid, bucket in enumerate(strat.buckets):
            pstores = []
            for desc in bucket.descs:
                mode = desc.mode
                id = '%d_%d' % (buckid, mode)
                newfname = '%s_%s' % (fname, id)

                if mode == Mode.NOOP:
                    ptype = NoopPStore
                elif mode == Mode.STAT:
                    ptype = StatPStore
                elif mode == Mode.QUERY:
                    ptype = PStoreQuery
                elif mode == Mode.FULL_MAPFUNC:
                    ptype = PStore1
                elif mode == Mode.PT_MAPFUNC:
                    ptype = PStore2
                elif mode == Mode.FULL_MAPFUNC_BOX:
                    ptype = PStore2Box
                elif mode == Mode.PT_MAPFUNC_BOX:
                    ptype = PStore2Box
                elif mode == Mode.PTR:
                    if desc.spec.payload == Spec.BOX:
                        ptype = PStore3Box
                    else:
                        ptype = PStore3
                else:
                    raise RuntimeError, "run modes at the descriptor level cannot be composite: %d" % mode
                if mode == Mode.PTR and not desc.backward and (Spec.BOX in (desc.spec.outcoords,
                                                                            desc.spec.payload)):
                    raise RuntimeError, "Don't support Box implementation for forward prov store"

                if mode == Mode.PTR and not desc.backward:
                    ptype = create_ftype(ptype)

                pstore = ptype(op, run_id, newfname, desc)
                pstores.append(pstore)
                #print "created base store", pstore
            
            if len(pstores) == 1:
                pstore = pstores[0]
            else:
                pstore = CompositePStore(op, run_id, None, desc, pstores)
            #print "ultimately, became", pstore

            pstore.open(new=True)
            ret.append(pstore)
            
        if len(ret) == 1:
            ret = ret[0]
        else:
            # make sure we have a forward and backward store
            # and don't have two backward stores, or something
            # ridiculous like that
            fps, bps = None, None
            for ps in ret:
                if ps.strat.backward:
                    bps = ps
                else:
                    fps = ps
            if not fps or not bps:
                raise RuntimeError
            ret = MultiPStore(op, run_id, None, desc, bps, fps)
                    
                      
        return ret


        

        
    def get_query_pstore(self, op, run_id, strat):
        if (op, run_id) in self.pstores:
            pstore = self.pstores[(op, run_id)]
            for desc in strat.descs():
                if pstore.strat == desc:
                    return pstore

        # we should just be creating a QueryStore
        if Mode.QUERY not in strat.modes():
            raise RuntimeError, "trying to create query pstore :%s" % str(strat)

        return PStoreQuery(op, run_id, '', Desc(Mode.QUERY, Spec.default(), True))
            

    def get_pstore(self, op, run_id):

        if (op, run_id) in self.override:
            return self.override[(op, run_id)]

        if (op, run_id) in self.pstores:
            return self.pstores[(op, run_id)]

        fname = self.get_filename(op,run_id)
        strat = self.get_strategy(op, run_id).copy()
        pstore = self.create_pstore(op, run_id, fname, strat)
        if op not in self.old_strats:
            self.old_strats[op] = {}
        self.old_strats[op][run_id] = strat
        self.pstores[(op,run_id)] = pstore
        
        return pstore


    def delete_pstore(self, op, run_id):
        if op not in self.old_strats or run_id not in self.old_strats[op]:
            return False
        log.debug( "delete\t%s\t%s" % (op, run_id) )        
        strat = self.old_strats[op][run_id]
        self.old_strats[op][run_id] = Strat.query()
        self.rm_strats[(op,run_id)] = strat

        key = (op, run_id)
        if key in self.pstores:
            pstore = self.pstores[key]
            # pstore.delete()
            del self.pstores[key]
            self.rm_pstores[key] = pstore

        if key in self.override:
            del self.override[key]
        return True
    
    def restore_pstores(self):
        """
        """
        log.debug( "restore_pstores" )
        for (op, r), s in self.rm_strats.items():
            log.debug( "restored strat\t%s\t%s" % (op, r) )
            self.old_strats[op][r] = s
        for (op, r), p in self.rm_pstores.items():
            log.debug( "restored pstore\t%s\t%s\t%s" % (op, r, p) )
            self.pstores[(op,r)] = p
        self.rm_strats = {}
        self.rm_pstores = {}

    def set_reexec(self, op, run_id, pstore):
        self.override[(op, run_id)] = pstore
        
    def clear_reexec(self, op, run_id):
        del self.override[(op, run_id)]
    
