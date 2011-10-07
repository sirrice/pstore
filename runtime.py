import sqlite3, time, logging, sys, pickle
from provstore import *
from strat import *


log = logging.getLogger('runtime')
logging.basicConfig()
log.setLevel(logging.ERROR)

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
        

    @staticmethod
    def instance():
        if Runtime._instance is None:
            Runtime._instance = Runtime()
        return Runtime._instance

    def set_strategy(self, op, strat):
        self.cur_strats[op] = strat

    def get_strategy(self, op, run_id = None):
        if run_id == None or op not in self.old_strats or run_id not in self.old_strats[op]:
            return self.cur_strats[op]
        return self.old_strats[op][run_id]

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

                if mode == Mode.QUERY:
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
                if ps.spec.backward:
                    bps = ps
                else:
                    fps = ps
            if not fps or not bps:
                raise RuntimeError
            ret = MultiPStore(op, run_id, None, desc, bps, fps)
                    
                      
        return ret


    def get_pstore(self, op, run_id):

        if (op, run_id) in self.override:
            return self.override[(op, run_id)]

        if (op, run_id) in self.pstores:
            return self.pstores[(op, run_id)]

        fname = self.get_filename(op,run_id)
        strat = self.get_strategy(op, run_id)
        pstore = self.create_pstore(op, run_id, fname, strat)
        if op not in self.old_strats:
            self.old_strats[op] = {}
        self.old_strats[op][run_id] = strat
        self.pstores[(op,run_id)] = pstore
        
        return pstore

    def set_reexec(self, op, run_id, pstore):
        self.override[(op, run_id)] = pstore
        
    def clear_reexec(self, op, run_id):
        del self.override[(op, run_id)]
    
