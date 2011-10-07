import math
from struct import *
import numpy as np
import os
from operator import mul
import time

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

    @check_w
    def add(self, coord):
        self.matrix[tuple(coord)] = True

    @check_w
    def add_set(self, coords):
        idxs = zip(*coords)
        self.matrix[idxs[0], idxs[1]] = True


    @check_r
    def get(self, coord):
        return self.matrix[tuple(coord)] == True

    def clear(self):
        self.matrix[:] = False

    @check_w
    def intersect(self, pqres):
        for coord in self:
            if coord not in pqres:
                self.matrix[tuple(coord)] = False

    def _save(self):
        fname = '%s.npy' % self.fprefix()
        np.save(fname, self.matrix)

    def _load(self):
        fname = '%s.npy' % self.fprefix()
        if not os.path.exists(fname):
            return False
        self.matrix = np.load(fname)
        return True

    @check_r
    def __len__(self):
        return self.matrix.sum()

    @check_r        
    def __iter__(self):
        return map(tuple,np.argwhere(self.matrix)).__iter__()

    def _close(self):
        self.matrix = None





class QIter(object):
    def __init__(self, container):
        self.c = container
    def next(self):
        pass
    def __len__(self):
        return len(self.c.child)
    def __iter__(self):
        return self
    

class PtrPstoreIter(QIter):
    def __init__(self, container, get_coords):
        """
        get_coords: takes pstore, returns coordinates
        """
        super(PtrPstoreIter, self).__init__(container)
        self.pqres = None
        self.generator = get_coords(self.c.pstore)
        self.iter = iter(self.generator)
    def next(self):
        return self.iter.next()


class ReexecPstoreIter(QIter):
    def __init__(self, container, get_coords):
        """
        get_coords: takes (pstore, pqres), returns coordinates
        """
        super(ReexecPstoreIter, self).__init__(container)
        print "reexec", self.c.pstore
        self.pqres = PQDenseResult(self.c.shape)
        self.c.pstore.start_q()
        get_coords(self.c.pstore, self.pqres)
        self.c.pstore.stop_q()
        self.iter = iter(self.pqres)

    def next(self):
        return self.iter.next()
    

class Query(object):
    def __init__(self, pstore, child, arridx, shape,  **kwargs):
        """
        child simply must be an iterablei,
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


class BPstoreQuery(Query):
    def __init__(self, *args, **kwargs):
        super(BPstoreQuery,self).__init__(*args, **kwargs)

    def __iter__(self):
        if self.closed: raise RuntimeError
        f = lambda pstore: pstore.backward_set(self.arridx, self.child, None,
                                               maxcount=self.maxcount)            
        return PtrPstoreIter(self, f)

class FPstoreQuery(Query):
    def __init__(self, *args, **kwargs):
        super(FPstoreQuery,self).__init__(*args, **kwargs)

    def __iter__(self):
        if self.closed: raise RuntimeError
        f = lambda pstore: pstore.forward_set(self.arridx, self.child, None,
                                              maxcount=self.maxcount)
        return PtrPstoreIter(self, f)


class BReexecQuery(Query):
    def __init__(self, *args, **kwargs):
        super(BReexecQuery,self).__init__(*args, **kwargs)
        self.pqres = PQDenseResult(self.shape)
        self.pstore.start_q()
        start = time.time()
        self.pstore.backward_set(self.arridx, self.child, self.pqres, maxcount=self.maxcount)
        print "reexec took", (time.time()-start), self.pstore
        self.pstore.stop_q()
        self.child.close()
        
    def close(self):
        super(BReexecQuery,self).close()
        self.pqres = None

    def __len__(self):
        return len(self.pqres)

    def __iter__(self):
        return iter(self.pqres)


class FReexecQuery(Query):
    def __init__(self, *args, **kwargs):
        super(FReexecQuery,self).__init__(*args, **kwargs)
        self.pqres = PQDenseResult(self.shape)
        self.pstore.start_q()
        start = time.time()
        self.pstore.forward_set(self.arridx, self.child, self.pqres, maxcount=self.maxcount)
        print "reexec took", (time.time()-start), self.pstore
        self.pstore.stop_q()
        self.child.close()
        
    def close(self):
        super(FReexecQuery,self).close()
        self.pqres = None


    def __len__(self):
        return len(self.pqres)


    def __iter__(self):
        return iter(self.pqres)

class FuncQuery(Query):
    def __init__(self, *args, **kwargs):
        super(FuncQuery, self).__init__(*args, **kwargs)
        if 'f' in kwargs:
            self.f = kwargs['f']
        else:
            self.f = lambda x: (x,)

    def __iter__(self):
        class FuncIter(QIter):
            def __init__(self, *args, **kwargs):
                super(FuncIter,self).__init__(*args, **kwargs)
                self.citer = iter(self.c.child)
                self.iter = None

            def next(self):
                if self.iter:
                    for outcoord in self.iter:
                        return outcoord

                nseen = 0
                for coord in self.citer:
                    self.iter = iter(self.c.f(coord))
                    nseen += 1
                    for outcoord in self.iter:
                        return outcoord
                raise StopIteration
        return FuncIter(self)
                
            
    

class Scan(Query):
    def __init__(self, child):
        super(Scan, self).__init__(None, child, -1, [1])

    def close(self):
        self.child = None

    def __len__(self):
        return len(self.child)


    def __iter__(self):
        return iter(self.child)

class DedupQuery(Query):
    def __init__(self, child, shape=None):
        if not shape:
            if not hasattr(child, 'shape'): raise RuntimeError
            shape = child.shape
        super(DedupQuery, self).__init__(None, child, -1, shape)
        self.pqres = PQDenseResult(self.shape)

        self.load()

    def load(self):
        nseen = 0
        for coord in self.child:
            self.pqres.add(coord)
            nseen += 1
            if nseen >= self.maxcount and  len(self.pqres) >= self.maxcount:
                print "dedup: found maxcount cells!", self.maxcount, nseen
                break
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


    

            
    

class ArrayPointers(object):
    """
    Data structure to store coordinate set for an array
    """
    def __init__(self, shape, grid_size=0.1):
        self.shape = shape

        blocksize = []  # dimension size of a block 
        nblocks = []    # number of blocks per dimension
        for l in shape:
            d = float(max(1, int(math.ceil(l * grid_size))))
            nblocks.append(int(math.ceil(l / d)))
            blocksize.append(int(d))
        blockmult = [reduce(mul, nblocks[i+1:], 1) for i in xrange(len(nblocks))]
        blockmult[-1] = 1
        cellmult = [reduce(mul, blocksize[i+1:], 1) for i in xrange(len(blocksize))]
        cellmult[-1] = 1
        nperblock = reduce(mul, blocksize)
        nblocks = reduce(mul, nblocks, 1)

        self.blocksize = blocksize
        self.nblocks = int(nblocks)
        self.blockmult = blockmult
        self.cellmult = cellmult
        self.nperblock = int(nperblock)

        self.enccost = 0.0
        self.deccost = 0.0
        self.addcost = 0.0

        self.ptrs = [None for i in xrange(self.nblocks)] 
        
    def enc(self, coord):
        start = time.time()
        blockid, cellid = 0,0
        
        for c, bsize, bmult, cmult in zip(coord,
                                          self.blocksize,
                                          self.blockmult,
                                          self.cellmult):
            blockid += math.floor(c/bsize) * bmult
            cellid += (c % int(bsize)) * cmult

        self.enccost += (time.time() - start)
        return int(blockid), int(cellid)


    def dec(self, blockid, cellid):
        start = time.time()
        blockcoord = dec(blockid, self.blockmult)
        cellcoord = dec(cellid, self.cellmult)
        ret = tuple((g * size + c for g, c, size in zip(blockcoord,
                                                        cellcoord,
                                                        self.blocksize)))
        self.deccost += (time.time() - start)
        return ret

    def add(self, coord):
        blockid, cellid = self.enc(coord)

        start = time.time()
        grid = self.ptrs
        if grid == True: return

        block = grid[blockid]
        if block is None:
            #block = BitVector(size=self.nperblock)
            block = set()
            grid[blockid] = block
        if block is not True:
            #block[cellid] = 1
            block.add(cellid)
        #if block is not True and block.count_bits_sparse() >= self.nperblock:
        if block is not True and len(block) >= self.nperblock:        
            grid[blockid] = True
        self.addcost += (time.time() - start)


    def check(self, blockid, cellid):
        return self.ptrs[blockid] is True or self.ptrs[blockid][cellid] is 1

    def coords(self):
        for blockid, block in enumerate(self.ptrs):
            if not block: continue
            if block is True:
                block = xrange(self.nperblock)
                for cellid in block:
                    yield self.dec(blockid, cellid)
            else:
                for cellid in block:
                    yield self.dec(blockid, cellid)
                # idx = block.next_set_bit()
                # while idx is not -1:
                #     yield self.dec(blockid, idx)                    
                #     idx = block.next_set_bit(idx+1)



if __name__ == '__main__':
    import random
    random.seed(0)
    
    shape = (100, 100)
    for klass in [PQSparseResult, PQDenseResult]:
        r = klass(shape)

        random.seed(0)
        for i in xrange(1000):
            x,y = random.randint(0, 99), random.randint(0,99)
            r.add((x,y))

        print r.closed, r.saved

        r.close()
        count = 0
        for coord in r:
            for coord in r:
                count += 1
        print "nfound^2", count
        
        r.close()        

        print r.closed, r.saved

        print "nfound", len(r)

        print r.closed, r.saved
        r.close()
        print r.closed, r.saved        
        print

        
        
