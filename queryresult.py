import time, logging, os, math
from struct import *
import numpy as np
from operator import mul

qrlog = logging.getLogger('queryres')
logging.basicConfig()
qrlog.setLevel(logging.ERROR)



def check_r(f):
    def ff(*args, **kwargs):
        self = args[0]
        if self.closed:
            self.check_open()
        return f(*args, **kwargs)
    return ff

def check_w(f):
    def ff(*args, **kwargs):
        self = args[0]
        if self.saved:
            raise RuntimeError, "trying to modify read-only resultset"
        if self.closed:
            self.check_open()
        return f(*args, **kwargs)
    return ff


class PQResult(object):
    """
    If saved, it is immutable.
    States:
    saved    closed   
     0         0       write-able
     0         1       invalid state
     1         0       read only
     1         1       must be opened, read only
    On first close() or save(), resultset becomes saved
    subsequent modifications throw exceptions
    """
    _id = 0
    
    def __init__(self, shape, id=None):
        self.shape = shape
        if id == None:
            self.id = PQResult._id
            PQResult._id += 1
        else:
            self.id = id

        self.saved = False  # has it been saved before.  if so, should be read only
        self.closed = False  # is it currently closed?

    def fprefix(self):
        return "pqres/%dx%d_id%d" % (self.shape[0], self.shape[1], self.id)

    @check_r
    def __iter__(self):
        pass

    @check_r
    def __in__(self, obj):
        return self.get(obj)

    @check_w
    def union(self, pqres):
        for coord in pqres:
            self.add(coord)

    @check_w
    def add(self, coord):
        pass

    @check_r
    def get(self, coord):
        pass

    def clear(self, coord):
        pass

    def save(self):
        if self.saved:
            raise RuntimeError, "ResultSet already saved"
        self._save()
        self.saved = True


    def load(self):
        self._load()

    def check_open(self):
        if not self.closed: return
        self.load()
        self.closed = False
    
    def close(self):
        if self.closed: return
        if not self.saved:
            self.save()
        self._close()
        self.closed = True
    

class PQSparseResult(PQResult):
    def __init__(self, shape, id=None):
        super(PQSparseResult,self).__init__(shape, id)
        self.coords = set()

    @check_w
    def add(self, coord):
        self.coords.add(tuple(coord))

    @check_w
    def add_set(self, coords):
        self.coords.update(coords)

    @check_r
    def get(self, coord):
        return tuple(coord) in self.coords

    def clear(self):
        self.coords = set()

    @check_w
    def intersect(self, pqres):
        self.coords = [coord for coord in self.coords if coord in pqres]

    @check_r
    def __iter__(self):
        return self.coords.__iter__()

    @check_r
    def __len__(self):
        return len(self.coords)

    def _save(self):
        f = file('%s.arr' % self.fprefix(), 'wb')
        arr = [0] * (len(self.coords) * 2 + 1)
        arr[0] = len(self.coords)
        idx = 1
        for coord in self.coords:
            arr[idx] = coord[0]
            arr[idx+1] = coord[1]
            idx += 2
        f.write(bytearray(pack("I"*len(arr), *arr)))
        f.close()


    def _load(self):
        fname = '%s.arr' % self.fprefix()
        if not os.path.exists(fname):
            return False
        self.coords = set()        
        f = file(fname, 'rb')
        n, = unpack("I", f.read(4))
        for i in xrange(n):
            coord = unpack("II", f.read(8))
            self.coords.add(coord)
        f.close()
        return True

    def _close(self):
        self.coords = None

        

class PQDenseResult(PQResult):
    def __init__(self, shape, id=None):
        super(PQDenseResult,self).__init__(shape, id)
        self.matrix = np.zeros(shape, bool)
        self.count = 0

    @check_w
    def add(self, coord):
        if not self.get(coord):
            self.count += 1
            self.matrix[tuple(coord)] = True

    @check_w
    def add_set(self, coords):
        idxs = zip(*coords)
        self.matrix[idxs[0], idxs[1]] = True
        self.count = self.matrix.sum()


    @check_r
    def get(self, coord):
        return self.matrix[tuple(coord)] == True

    def clear(self):
        self.matrix[:] = False
        self.count = 0

    @check_w
    def intersect(self, pqres):
        for coord in self:
            if coord not in pqres:
                self.matrix[tuple(coord)] = False
                self.count -= 1

    def _save(self):
        fname = '%s.npy' % self.fprefix()
        np.save(fname, self.matrix)

    def _load(self):
        fname = '%s.npy' % self.fprefix()
        if not os.path.exists(fname):
            return False
        self.matrix = np.load(fname)
        self.count = self.matrix.sum()
        return True

    @check_r
    def __len__(self):
        return self.count

    @check_r        
    def __iter__(self):
        return map(tuple,np.argwhere(self.matrix)).__iter__()

    def _close(self):
        self.matrix = None




            

class Query(object):
    """
    Embodies a database iterator.
    Instead of open() and close(), returns a new iterator when
    __iter__ is called
    """
    
    def __init__(self, pstore, child, arridx, shape,  **kwargs):
        """
        child simply must be an iterable
        """
        self.pstore = pstore
        self.child = child
        self.arridx = arridx
        self.shape = shape
        self.maxcount = reduce(mul, self.shape)
        self.closed = False
        self.__mylen__ = None
        
    def close(self):
        if self.pstore != None:
            self.pstore.close()
        if self.child != None:
            self.child.close()
        self.pstore = None
        self.child = None
        self.closed = True


    def __len__(self):
        if self.__mylen__ == None:
            self.__mylen__ = 0
            for x in self:
                self.__mylen__ += 1
        return self.__mylen__
    
    def __iter__(self):
        if self.closed: raise RuntimeError
        return QIter(self)
    

class Scan(Query):
    """
    Query wrapper for native iterables (lists, sets, etc)
    """
    def __init__(self, child):
        super(Scan, self).__init__(None, child, -1, [1])

    def close(self):
        self.child = None

    def __len__(self):
        return len(self.child)


    def __iter__(self):
        return iter(self.child)

class DedupQuery(Query):
    """
    De-duplication operator.  Sucks data from the child iterator into
    an in-memory binary matrix to perform deduplication
    """
    def __init__(self, child, shape=None):
        if not shape:
            if not hasattr(child, 'shape'): raise RuntimeError
            shape = child.shape
        super(DedupQuery, self).__init__(None, child, -1, shape)
        if reduce(mul, self.shape) > 1000000:
            self.pqres = PQDenseResult(self.shape)
        else:
            self.pqres = PQSparseResult(self.shape)
        self.load()

    def load(self):
        nseen = 1
        loadcost = 0.0
        for coord in self.child:
            start = time.time()
            self.pqres.add(coord)
            loadcost += time.time() - start
            nseen += 1
            if nseen >= self.maxcount and  len(self.pqres) >= self.maxcount:
                qrlog.debug( "dedup: found maxcount cells!\t%d\t%d", self.maxcount, nseen )
                break

        qrlog.debug( "loadcost: \t%f\t%d", (loadcost / nseen, nseen) )
        self.child.close()

    def __len__(self):
        return len(self.pqres)

    def close(self):
        super(DedupQuery,self).close()
        self.pqres = None

    def __iter__(self):
        return iter(self.pqres)


class NBDedupQuery(Query):
    def __init__(self, child, shape=None):
        if not shape:
            if not hasattr(child, 'shape'): raise RuntimeError
            shape = child.shape
        super(NBDedupQuery, self).__init__(None, child, -1, shape)

    def __len__(self):
        return len(self.child)

    def __iter__(self):
        class NBDedupIter(object):
            def __init__(self, par):
                self.par = par
                self.iter = iter(par.child)
                self.pqres = PQSparseResult(par.shape)

            def next(self):
                for coord in self.iter:
                    if coord not in self.pqres:
                        self.pqres.add(coord)
                        return coord
                self.pqres = None
                raise StopIteration
        return NBDedupIter(self)


    

            

if __name__ == '__main__':
    import random
    random.seed(0)
    
    shape = (1000, 500)
    maxcount = reduce(mul, shape)
    n = 1000000
    coords = [ (random.randint(0, shape[0]-1),
               random.randint(0, shape[1]-1))
               for i in xrange(n) ]

    for klass in [PQSparseResult, PQDenseResult]:
        r = klass(shape)
        random.seed(0)

        start = time.time()
        nseen = 0
        for coord in coords:
            r.add(coord)
            nseen += 1
            if nseen >= maxcount and  len(r) >= maxcount:
                break
        cost = time.time() - start


        print klass, cost
