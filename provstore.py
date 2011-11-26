import struct, math, time, random
import bsddb3 as bsddb
from StringIO import StringIO
from operator import mul, or_, and_
from ctypes import create_string_buffer

from strat import *
from queryresult import *
from util import subarray, get_db
from rtree import index


plog = logging.getLogger('provstore')
logging.basicConfig()
plog.setLevel(logging.ERROR)



#
# obj is defined as unparsed inputs and outputs
# (serialized outputs, serialized inputs)
# basically, don't follow pointers 
# if KEY: the key
# if coords: list of encoded coords
# if coord:  list of single encoded coord
# 
#

def enc(coord, shape):
    return coord[0] * shape[1] + coord[1]

def dec(val, shape):
    return (int(math.floor(float(val)/shape[1])), val % shape[1])

def bbox(coords):
    if not len(coords):
        return ((-1, -1), (-1,-1))
    if len(coords) == 1:
        return (coords[0], coords[0])
    minx, maxx, miny, maxy = None, None, None, None
    minc, maxc = None, None
    for coord in coords:
        if minx == None:
            minx = maxx = coord[0]
            miny = maxy = coord[1]
        else:
            minx = min(minx, coord[0])
            maxx = max(maxx, coord[0])
            miny = min(miny, coord[1])
            maxy = max(maxy, coord[1])

    return ((minx, miny), (maxx, maxy))


__grid__ = None
__gcells__ = 0
def gengrid(coords):
    """
    @return (box, negs)
    box: bounding box containing the coordinates
    negs: encoded list of coords not in bounding box
    """
    global __grid__, __gcells__
    if len(coords) == 0:
        return ((0,0), (0,0)), ((0,0),)
    box = bbox(coords)
    shape = (1+box[1][0]-box[0][0], 1+box[1][1]-box[0][1])
    ncells = reduce(mul, shape)
    if __gcells__ < ncells:
        __grid__ = np.ones((1, ncells), dtype=bool)
        __gcells__ = reduce(mul, __grid__.shape)
    arr = __grid__[:,:ncells].reshape((shape))
    arr[:,:] = True
    newcoords = [ (coord[0]-box[0][0], coord[1]-box[0][1]) for coord in coords ]
    arr[zip(*newcoords)] = False
    #encs = map(lambda coord: enc(coord, shape), np.argwhere(arr))
    return (box, map(tuple, np.argwhere(arr)))


def decgrid(box, negs):
    global __grid__, __gcells__    
    shape = ( 1+box[1][0]-box[0][0], 1+box[1][1]-box[0][1] )
    ncells = reduce(mul, shape)
    if __gcells__ < ncells:
        __grid__ = np.ones((1, ncells), dtype=bool)
        __gcells__ = reduce(mul, __grid__.shape)
    arr = __grid__[:,:ncells].reshape((shape))
    arr[:,:] = True
    arr[zip(*negs)] = False
    return map(lambda coord: (coord[0]+box[0][0], coord[1]+box[0][1]), np.argwhere(arr))





        
class BinaryRTree(index.Index):
    def dumps(self, obj):
        return obj
    def loads(self, s):
        return s

class SpatialIndex(object):
    def __init__(self, fname):
        self.fname = '%s_rtree' % fname
        self.rtree = None
        #self.bdb = None
        self.id = 0
        p = index.Property()
        p.dimension = 2
        p.filename = self.fname
        self.p = p
        

    def set(self, box, val):
        #self.rtree.add(self.id, box)#, val)
        self.rtree.add(val, box)#, val)
        self.id += 1
        #if val is not None:
         #   self.bdb[str(self.id)] = val

    def get_pt(self, pt):
        """
        @return generator of bdb keys
        """
        for id in self.rtree.intersection(pt):
            yield 'b:%d' % id
        #return self.rtree.intersection(pt)#, objects = True)
        # for item in self.rtree.intersection(pt, objects = True):
        #     item.object = self.bdb.get(str(item.id), None)
        #     yield item

    def get_box(self, box):
        """
        @return generator of bounding boxes
        """
        for item in self.rtree.intersection(box, objects=True):
            yield (item.bbox[:2], item.bbox[2:])
        #return self.rtree.intersection(box)#, objects = True)
        # for item in self.rtree.intersection(box, objects = True):
        #     item.object = self.bdb.get(str(item.id), None)
        #     yield item

    def disk(self):
        try:
            idxsize = os.path.getsize('%s.%s' % (self.p.get_filename(), self.p.get_idx_extension()))
            datsize = os.path.getsize('%s.%s' % (self.p.get_filename(), self.p.get_dat_extension()))
            #bdbsize = os.path.getsize('%s.%s' % (self.p.get_filename(), 'bdb'))
            return idxsize + datsize
        except:
            raise
            return 0


    def open(self, new=False):
        p = self.p
        if new:
            try:
                self.delete()
            except:
                pass
        self.rtree = BinaryRTree(p.get_filename(), properties=p)
        #self.bdb = bsddb.hashopen('%s.bdb'% p.get_filename(), 'n')

    def close(self):
        self.rtree.close()
        #self.bdb.close()

    def delete(self):
        p = self.p
        try:
            os.unlink('%s.%s' % (p.get_filename(), p.get_idx_extension()))
            os.unlink('%s.%s' % (p.get_filename(), p.get_dat_extension()))
            #os.unlink('%s.%s' % (p.get_filename(), 'bdb'))
        except:
            raise

    

def in_box(incoord, box):
    if not box[0]: return False
    return reduce(and_, [box[0][dim] <= incoord[dim] and
                         incoord[dim] <= box[1][dim]
                         for dim in xrange(len(incoord))])

def instrument(fn):
    func_name = fn.__name__
    def w(self, *args, **kwargs):
        start = time.time()
        ret = fn(self, *args, **kwargs)

        if func_name not in self.stats:
            self.stats[func_name] = 0
        self.stats[func_name] += (time.time() - start)
        return ret
    return w

def alltoall(fn):
    def w(self, left, arridx, backward=True):
        
        if backward:
            shape = self.outshape
        else:
            shape = self.inshapes[arridx]
        maxcells = reduce(mul, shape)

        if not self.op.alltoall(arridx) or len(left) < maxcells:
            return fn(self, left, arridx, backward=backward)
        if backward:
            return AllToAllScan(self.inshapes[arridx])
        return AllToAllScan(self.outshape)
    return w













class IPstore(object):
    def __init__(self, op, run_id, fname, strat):
        self.op = op
        self.run_id = run_id
        self.fname = fname
        self.strat = strat
        self.spec = strat.spec
        self.inshapes = op.wrapper.get_input_shapes(run_id)
        self.outshape = self.op.output_shape(run_id)
        self.nargs = op.wrapper.nargs

        self.stats = {}

        # whether stats are operator given, or calculated by sampling
        self.sampled = True 
        self.outarea = 0
        self.inareas = [0] * self.nargs
        self.fanins = [0] * self.nargs
        self.noutcells  = 0
        self.ncalls = 0

    def clear_stats(self):
        self.stats = {}

    def get_stat(self, attr, default=0):
        return self.stats.get(attr, default)

    def inc_stat(self, attr, val):
        if attr not in self.stats:
            self.stats[attr] = 0.0
        self.stats[attr] += val            

    def uses_mode(self, mode):
        for m in self.strat.modes():
            if mode & m: return True
        return False

    def get_iter(self):
        return ()

    def disk(self):
        return 0

    def indexsize(self):
        return 0

    def extract_outcells(self, obj):
        return (obj[0],)

    def extract_outbox(self, obj):
        raise RuntimeError

    def extract_incells(self, obj, arridx):
        return ()

    def extract_inbox(self, obj):
        raise RuntimeError

    def join(self, left, arridx, backward=True):
        raise RuntimeError


    def get_fanins(self):
        return self.fanins
    def get_inareas(self):
        return self.inareas
    def get_oclustsize(self):
        return self.outarea
    def get_densities(self):
        return map(lambda p: p[1] > 0 and float(p[0]) / p[1] or 1.0, zip(self.get_fanins(), self.get_inareas()))
    def get_nptrs(self):
        return sum([self.noutcells * fanin for fanin in self.get_fanins()])
    

    def set_fanins(self, fanins):
        self.fanins = fanins
        self.sampled = False
    def set_inareas(self, areas):
        self.inareas = areas
        self.sampled = False
    def set_outarea(self, area):
        self.outarea = area
        self.sampled = False
    def set_ncalls(self, n):
        self.ncalls = n
        self.sampled = False
    def set_noutcells(self, n):
        self.noutcells = n
        self.sampled = False
    def _in_densities(self):
        return map(lambda p: p[0] / p[1], zip(self.fanins, self.inareas))
    in_densities = property(_in_densities)
    def _out_density(self):
        return self.noutcells / self.outarea
    out_density = property(_out_density)

    def open(self, new=False):
        pass
    
    def close(self):
        pass

    def delete(self):
        pass


class NoopPStore(IPstore):
    def uses_mode(self, mode):
        return False
    
    def write(self, outcoords, *incoords_arr):
        pass

class StatPStore(IPstore):
    def __init__(self, op, run_id, fname, strat):
        super(StatPStore, self).__init__(op, run_id, fname, strat)
        self.nsampled = 0
        self.nskipped = 0

    def uses_mode(self, mode):
        for m in ( Mode.PTR, Mode.FULL_MAPFUNC ):
            if mode & m: return True
        return False
    
    @instrument
    def update_stats(self, outcoords, *incoords_arr):
        outbox = bbox(outcoords)
        inboxes = map(bbox, incoords_arr)

        diff = lambda box: map(lambda p: p[1]-p[0]+1, zip(box[0], box[1]))
        outarea = reduce(mul, diff(outbox))
        inareas = map(lambda box: reduce(mul, diff(box)), inboxes)

        fanins = map(len, incoords_arr)
        noutcells = len(outcoords)

        if inareas[0] < fanins[0]:
            raise "we have a problem", inareas, fanins

        self.outarea += outarea
        self.inareas = map(sum, zip(self.inareas, inareas))
        self.fanins = map(sum, zip(self.fanins, map(lambda x: noutcells * x, fanins)))
        self.noutcells += noutcells
        self.ncalls += 1

        # update the histograms
        

    @instrument
    def write(self, outcoords, *incoords_arr):
        if self.nsampled > 10 and random.random() < 0.9:
            self.nskipped += 1
            return
        self.update_stats(outcoords, *incoords_arr)
        self.nsampled += 1

    def close(self):
        if self.nsampled == 0: return
        if self.ncalls == 0:
            self.fanins = [1] * self.nargs
            self.inareas = [1] * self.nargs
            self.outarea = 1
            
            self.ncalls = 1
            self.noutcells = 1
            return
        
        self.fanins = map(lambda fanin: float(fanin) / self.noutcells, self.fanins)
        self.inareas = map(lambda inarea: float(inarea) / self.ncalls, self.inareas)
        self.outarea = self.noutcells / float(self.ncalls)

        self.noutcells = (self.noutcells / self.nsampled) * (self.nsampled + self.nskipped)
        self.ncalls = self.nsampled + self.nskipped
        
        
    

class PStore1(IPstore):

    def __init__(self, op, run_id, fname, strat):
        super(PStore1, self).__init__(op, run_id, fname, strat)
        self.spec = strat.spec

    def uses_mode(self, obj):
        return obj & Mode.FULL_MAPFUNC != 0

    def extract_outcells(self, obj):
        return (obj[0],)


    def extract_incells(self, obj, arridx):
        return self.op.bmap(obj[0][0], self.run_id, arridx)

    def extract_box(self, obj, arridx):
        return self.op.bbox(obj[0][0], self.run_id, arridx)

    @alltoall
    def join(self, left, arridx, backward=True):
        """
        returns a generator that yields coordinates
        """
        if backward:
            extract = lambda coord: self.extract_incells(((coord,), None), arridx)
        else:
            extract = lambda coord: self.op.fmap(coord, self.run_id, arridx)

        try:
            for l in left:
                for coord in extract(l):
                    yield coord
        except Exception, e:
            plog.error('%s:\t%s\tbackward(%s)', str(e), str(self.op), backward )
            raise


        
class DiskStore(IPstore):

    def __init__(self, op, run_id, fname, strat):
        super(DiskStore, self).__init__(op, run_id, fname, strat)
        self.bdb = None
        self.n = 0
        self.outidx = SpatialIndex(fname)

    def get_iter(self):
        badprefixes = ('key:', 'b:', 'ref:')
        for key in self.bdb.keys():
            if key.startswith("key:")  or key.startswith('b:') or key.startswith('ref:'):
                continue
            yield (key, self.bdb[key])

    def disk(self):
        try:
            return os.path.getsize(self.fname) + self.outidx.disk()
        except:
            return 1

    def indexsize(self):
        try:
            return self.outidx.disk()
        except:
            raise
            return 1

    def extract_outcells_enc(self, obj):
        """
        NOTE: returns integer encoded coords
        """
        if self.spec.outcoords == Spec.KEY:
            ser_key = self._parse(StringIO(obj[0]), Spec.KEY)
            ser_vals = self.bdb[ser_key]
            vals = self._parse(StringIO(ser_vals), Spec.COORD_MANY)
            return vals
        elif self.spec.outcoords == Spec.GRID:
            coords = map(self.dec_out, obj[0])
            box, negs = coords[:2], coords[3:]
            return map(self.enc_out, decgrid(box, negs))
            
        ret = self._parse(StringIO(obj[0]), self.spec.outcoords)
        return ret

    def extract_outcells(self, obj):
        if self.spec.outcoords == Spec.GRID:
            coords = map(self.dec_out, obj[0])
            box, negs = coords[:2], coords[3:]
            return decgrid(box, negs)
        return map(self.dec_out, self.extract_outcells_enc(obj))


    def extract_outboxes(self, obj):
        if self.spec.outcoords == Spec.BOX:
            return map(lambda c: self.dec_out, self._parse(StringIO(obj[0]), self.spec.outcoords))
        raise RuntimeError

    def extract_allincells_enc(self, obj):
        buf = StringIO(obj[1])
        parsed = []
        for idx in xrange(self.nargs):
            parsed.append(self._parse(buf, self.spec.payload))

        if Spec.KEY == self.spec.payload:
            ret = []
            for ser_key in parsed:
                encs = self._parse(StringIO(self.bdb[ser_key]), Spec.COORD_MANY)
                ret.append(encs)
            return ret
        elif Spec.GRID == self.spec.payload:
            ret = []
            for arridx, encs in enumerate(parsed):
                coords = map(self.dec_out, encs)
                box, negs = coords[:2], coords[3:]
                encs = map(lambda coord: self.enc_in(coord, len(ret)), decgrid(box, negs))
                ret.append(encs)
            return ret
        elif Spec.is_coord(self.spec.payload):
            return parsed
        raise RuntimeError


    def extract_incells_enc(self, obj, arridx):
        """
        Read incells from an obj object
        NOTE: returns integer encoded coords
        """
        buf = StringIO(obj[1])
        for idx in xrange(arridx):
            self._parse(buf, self.spec.payload)
        val = self._parse(buf, self.spec.payload)

        if Spec.KEY == self.spec.payload:
            ser_key = val
            ser_encs = self.bdb[ser_key]
            encs = self._parse(StringIO(ser_encs), Spec.COORD_MANY)
            return encs
        elif Spec.GRID == self.spec.payload:
            coords = map(self.dec_out, val)
            box, negs = coords[:2], coords[3:]
            return map(lambda coord: self.enc_in(coord, arridx), decgrid(box, negs))
        elif Spec.is_coord(self.spec.payload):
            return val
        raise RuntimeError

    def extract_incells(self, obj, arridx):
        if Spec.GRID == self.spec.payload:
            buf = StringIO(obj[1])
            for idx in xrange(arridx):
                self._parse(buf, self.spec.payload)
            val = self._parse(buf, self.spec.payload)
            coords = map(self.dec_out, val)
            box, negs = coords[:2], coords[3:]
            return decgrid(box, negs)
        return map(lambda val: self.dec_in(val, arridx),
                   self.extract_incells_enc(obj, arridx))

    def extract_inboxes(self, obj):
        if self.spec.payload == Spec.BOX:
            buf = StringIO(obj[1])
            boxes = []
            for arridx in xrange(self.nargs):
                box = self._parse(buf, self.spec.payload)
                boxes.append( (self.dec_in(box[0], arridx), self.dec_in(box[1], arridx)) )
            return boxes
        raise RuntimeError

    def enc_out(self, coord):
        return enc(coord, self.outshape)

    def dec_out(self, val):
        return dec(val, self.outshape)

    def enc_in(self, coord, arridx):
        return enc(coord, self.inshapes[arridx])

    def dec_in(self, val, arridx):
        return dec(val, self.inshapes[arridx])

    def get_key(self, outcoords):
        enc_outcoords = map(self.enc_out, outcoords)
        outbuf = StringIO()
        self._serialize(enc_outcoords, outbuf, self.spec.outcoords)
        return outbuf.getvalue()


    @instrument
    def _serialize_format(self, counts, mode):
        if mode == Spec.COORD_ONE:
            return '%dI' % sum(counts)
        elif mode == Spec.COORD_MANY:
            return '%dI' % (len(counts) + sum(counts))
        elif mode == Spec.BOX:
            return '%dI' % (2 * len(counts))
        elif mode == Spec.GRID:
            return '%dI' % sum(counts)
        elif mode == Spec.KEY:
            return '%dI' % (len(counts) + sum(counts))
        elif mode == Spec.BINARY:
            raise RuntimeError
        elif mode == Spec.NONE:
            return ''
        



    @instrument
    def _serialize(self, data, buf, mode):
        """
        NOTE: expects coordinates to be in encoded into integer format already
        """
        
        if mode == Spec.COORD_ONE:
            coords = data
            buf.write(struct.pack("I", coords[0]))
        elif mode == Spec.COORD_MANY:
            coords = data
            n = len(coords)
            s = struct.pack("I%dI" % n, n, *coords)
            buf.write(s)
        elif mode == Spec.BOX:
            box = data
            buf.write(struct.pack('2I', *box))
        elif mode == Spec.GRID:
            box, encs = data
            self._serialize(box, buf, Spec.BOX)
            self._serialize(encs, buf, Spec.COORD_MANY)
        elif mode == Spec.KEY:
            ser_coords = StringIO()
            self._serialize(data, ser_coords, Spec.COORD_MANY)
            ser_coords = ser_coords.getvalue()
            h = str(hash(ser_coords) % 4294967296).rjust(10, '_')
            ser_key = 'key:%s' % str(h)
            self.bdb[ser_key] = ser_coords
            s = struct.pack("I%ds" % len(ser_key), len(ser_key), ser_key)
            buf.write(s)
        elif mode == Spec.BINARY:
            buf.write(struct.pack("I", len(data)))
            buf.write(data)
        elif mode == Spec.NONE:
            return

    @instrument
    def _parse(self, buf, mode):
        """
        NOTE: returns coordinates in integer encoded format
        XXX: No support for multiple arridxs.  must do externally
        """
        if mode == Spec.COORD_ONE:
            return struct.unpack("I", buf.read(4))
        elif mode == Spec.COORD_MANY:
            n, = struct.unpack("I", buf.read(4))
            if n > 1000000:  raise RuntimeError, "trying to parse really large string %d" % n
            return struct.unpack("%dI" % n, buf.read(4*n))
        elif mode == Spec.BOX:
            minpt, maxpt = struct.unpack('2I', buf.read(8))
            return (minpt, maxpt)
        elif mode == Spec.GRID:
            pos = buf.tell()
            buf.seek(pos+8)
            n, = struct.unpack('I', buf.read(4))
            buf.seek(pos)
            return struct.unpack('%dI' % (3+n), buf.read(8 + 4 + 4*n))
        elif mode == Spec.KEY:
            n, ser_key = struct.unpack("I14s", buf.read(4 + 14))
            return ser_key
        elif mode == Spec.BINARY:
            n, = struct.unpack("I", buf.read(4))
            return buf.read(n)
        elif mode == Spec.NONE:
            return None

    @alltoall
    def join(self, left, arridx, backward=True):
        self.open()
        if backward and self.spec.outcoords == Spec.COORD_ONE:
            for coord in self.hash_join(left, arridx):
                yield coord
            self.close()
            return


        

        extractcost, parsecost, nhits = 0.0, 0.0, 0.0
        keycost = 0.0
        datacost = 0.0
        itemcost = 0.0

        if backward:  # backward query
            if Spec.BOX == self.spec.outcoords:
                pred = lambda obj: self._parse(StringIO(obj[0]), Spec.BOX)
                matches = in_box
                def getkey(box):
                    #box = item.bbox
                    #box = (box[:2],box[2:])
                    box = self.enc_out(box)
                    key = StringIO()
                    self._serialize(box, key, Spec.BOX)
                    return key.getvalue()
            else:
                pred = self.extract_outcells_enc
                matches = lambda coord, outcoords: self.enc_out(coord) in outcoords
                getkey = lambda key: self.bdb[key]#'b:%d' % item]
            extract = lambda obj: self.extract_incells(obj, arridx)


            keysize = 0
            valsize = 0

            niter = 0
            for l in left:
                niter += 1
                l = tuple(l)
                
                if self.spec.outcoords == Spec.BOX:
                    items = self.outidx.get_box(l)
                else:
                    items = self.outidx.get_pt(l)
                for item in items:
                    start = time.time()
                    key = getkey(item)
                    keysize += len(key)
                    keycost += time.time() - start

                    start = time.time()
                    r = (key, self.bdb[key])
                    valsize += len(r[1])
                    datacost += time.time() - start 

                    
                    start = time.time()
                    coords = pred(r)
                    b = matches(l, coords)
                    parsecost += time.time() - start                    

                    nhits += 1
                    if b:
                        start = time.time()
                        e = extract(r)
                        extractcost += time.time() - start
                        for coord in e:
                            yield coord
            self.inc_stat('parsecost', parsecost)
            self.inc_stat('keycost', keycost)
            self.inc_stat('datacost', datacost)
            self.inc_stat('extractcost', extractcost)
            self.inc_stat('nhits', nhits)


            plog.debug( "keycost  \t%f", keycost)
            plog.debug( "parsecost\t%f", parsecost)
            plog.debug( "extract  \t%f", extractcost)
            plog.debug( "nhits    \t%f\t%f", nhits, niter)
            if nhits > 0:
                plog.debug( "key/vallen\t%f\t%f", keysize/nhits, valsize/nhits) 

        else:  # forward query
            if Spec.BOX == self.spec.payload:
                def pred(coord, obj):
                    buf = StringIO(obj[1])
                    buf.seek(8*arridx)
                    box = self._parse(buf, Spec.BOX)
                    return box
                matches = in_box
            else:
                pred = lambda obj: self.extract_incells_enc(obj, arridx)
                matches = lambda coord, coords: coord in coords
                left = map(lambda l: self.enc_in(l, arridx), left)
            extract = lambda obj: self.extract_outcells(obj)



            for r in self.get_iter(): # r is a pair of unparsed strings
                start = time.time()
                coords = pred(r) #
                parsecost += time.time() - start

                for l in left:
                    b = matches(l, coords)
                    if b:
                        start = time.time()
                        e = extract(r)
                        extractcost += time.time() - start
                        for coord in e:
                            yield coord
                        break
            self.inc_stat('parsecost', parsecost)
            self.inc_stat('extractcost', extractcost)
        self.close()
        

    def hash_join(self, left, arridx):
        datacost = 0.0
        extractcost = 0.0
        keycost = 0.0
        niter = 1.0
        nhits = 1.0
        total = 0.0
        tstart = time.time()
        
        for l in left:
            niter += 1
            start = time.time()
            enc = self.enc_out(l)
            key = struct.pack('I', enc)
            keycost += time.time() - start
            if key not in self.bdb: continue

            nhits += 1
            start = time.time()
            payload = self.bdb[key]
            datacost += time.time() - start
            
            if payload is not None:
                start = time.time()
                e = self.extract_incells((key, payload), arridx)
                extractcost += time.time() - start
                for coord in e:
                    yield coord
                    
        self.inc_stat('parsecost', 0)
        self.inc_stat('keycost', keycost)
        self.inc_stat('datacost', datacost)
        self.inc_stat('extractcost', extractcost)
        self.inc_stat('nhits', nhits)
                    
        plog.debug( "datacost:    %f\t%f", datacost, datacost / nhits )
        plog.debug( "extractcost: %f\t%f", extractcost, extractcost / nhits )
        plog.debug( "keycost:     %f\t%f", keycost, keycost / niter )
        plog.debug( "nhits:     %f", nhits)        


    def open(self, new=False):
        if self.bdb != None: return
        if new:
            self.bdb = get_db(self.fname, new=True, fsync=False)
        else:
            self.bdb = get_db(self.fname, new=False, fsync=False)
        self.outidx.open(new=new)

    @instrument
    def close(self):
        global __grid__, __gcells__        
        try:
            self.bdb.close()
            self.bdb = None
            self.outidx.close()
            __grid__ = None
            __gcells__ = 0
        except Exception, e:
            pass

    def delete(self):
        """
        be careful with this guy!
        """
        try:
            os.unlink(self.fname)
            self.outidx.delete()
        except:
            raise


        


class PStore2(DiskStore):
    def __init__(self, op, run_id, fname, strat):
        super(PStore2, self).__init__(op, run_id, fname, strat)
        self.spec.payload = Spec.BINARY # ignore payload spec
        self.nbdbitems = 0
        self.reset_cache()
        self.outbuf = None
        self.inbuf = None

        self.memsize = 0
        

    def reset_cache(self):
        self.outcache = []
        self.outcounts = []
        self.outboxes = []
        self.incache = []
        self.incounts = []
        
        

    def uses_mode(self, mode):
        return mode & Mode.PT_MAPFUNC != 0

    @instrument
    def update_stats(self, outcoords, payload):
        diff = lambda box: map(lambda p: p[1]-p[0], zip(box[0], box[1]))
        
        outbox = bbox(outcoords)
        outarea = reduce(mul, diff(outbox))
        noutcells = len(outcoords)

        fanins = [1] * self.nargs
        inareas = [1] * self.nargs
        # for arridx in xrange(self.nargs):
        #     incells = map(tuple, self.op.bmap_obj((outcoords, payload), self.run_id, arridx))
        #     box = bbox(incells)
        #     inareas.append(reduce(mul, diff(box)))
        #     fanins.append(len(incells))

        self.outarea += outarea
        self.inareas = map(sum, zip(self.inareas, inareas))
        self.fanins = map(sum, zip(self.fanins, map(lambda x: noutcells * x, fanins)))
        self.noutcells += noutcells
        self.ncalls += 1
        

    def extract_incells_enc(self, obj, arridx):
        enc = lambda c: self.enc_in(c, arridx)
        return map( enc, self.extract_incells(obj, arridx) )

    def extract_incells(self, obj, arridx):
        outcoords = self.extract_outcells(obj)
        return self.op.bmap_obj((outcoords, self._parse(StringIO(obj[1]), Spec.BINARY)),
                                self.run_id,
                                arridx)

    def extract_inboxes(self, obj):
        boxes = []
        outcoords = self.extract_outcells(obj)
        newobj = (outcoords, obj[1])
        for arridx in xrange(self.nargs):
            boxes.append(self.op.bbox_obj(newobj, self.run_id, arridx))
        return boxes


    def get_val(self, payload):
        paybuf = StringIO()
        self._serialize(payload, paybuf, self.spec.payload)
        return paybuf.getvalue()

    def get_key(self, outcoords):
        if self.spec.outcoords == Spec.GRID:
            raise RuntimeError, 'not supported'
        enc_outcoords = map(self.enc_out, outcoords)
        outbuf = StringIO()
        self._serialize(enc_outcoords, outbuf, self.spec.outcoords)
        return outbuf.getvalue()

    def add_to_cache(self, cache, counts, coords, mode ):
        if Spec.BOX == mode:
            box = bbox(coords)
            cache.extend(box)
            counts.append(2)
        elif Spec.GRID == mode:
            grid = gengrid(coords)
            cache.extend( grid[0] )
            cache.append( (0, len(grid[1])) )
            cache.extend( grid[1] )
            counts.append(2 + 1 + len(grid[1]))
        else:
            if mode in (Spec.COORD_MANY, Spec.KEY):
                cache.append((0, len(coords)))
            cache.extend(coords)
            counts.append(len(coords))
        

    @instrument
    def write(self, outcoords, payload):
        #self.update_stats(outcoords, payload)
        start = time.time()
        if self.spec.outcoords != Spec.COORD_ONE:
            self.outboxes.append(bbox(outcoords))
        else:
            self.outboxes.append(None)
        self.add_to_cache(self.outcache, self.outcounts, outcoords, self.spec.outcoords)
        self.inc_stat('outcache', time.time()-start)

        start = time.time()
        self.incache.append(payload)
        self.incounts.append(len(payload))
        self.inc_stat('incache', time.time()-start)

        self.memsize += len(outcoords) * 2 + len(payload)

        if self.memsize > 1000000:
            self.flush()
            self.memsize = 0

    @instrument
    def flush(self):
        if len(self.outcounts) == 0: return

        # serialize key
        start = time.time()
        oencs = [self.enc_out(outcoord) for outcoord in self.outcache]
        fmt = self._serialize_format(self.outcounts, self.spec.outcoords)
        keysize = struct.calcsize(fmt)
        if self.outbuf is None or keysize > len(self.outbuf):
            self.outbuf = create_string_buffer(keysize)
        struct.pack_into(fmt, self.outbuf, 0, *oencs)
        self.inc_stat('serout', time.time() - start)

        # serialize val
        start = time.time()
        fmt = '%dI' % len(self.incounts)
        valsize = struct.calcsize(fmt)
        if self.inbuf is None or valsize > len(self.inbuf):
            self.inbuf = create_string_buffer(valsize)
        struct.pack_into(fmt, self.inbuf, 0, *self.incounts)
        self.inc_stat('serin', time.time() - start)

        def foo(buf, count, offset, mode):
            if mode in (Spec.COORD_MANY, Spec.KEY):
                end = offset + 4 * (1 + count)
            elif mode == Spec.GRID:
                end = offset + 4 * count
            elif mode == Spec.BOX:
                end = offset + 8
            elif mode == Spec.COORD_ONE:
                end = offset + 4
            elif mode == Spec.BINARY:
                end = offset + 4
            ret = buf[offset:end]

            if mode == Spec.KEY:
                h = str(hash(ret) % 4294967296)
                key = 'key:%s' % h.rjust(10, '_')
                self.bdb[key] = ret
                ret = struct.pack('I%ds' % len(key), len(key), key)
                return ret, end
            return ret, end

        bdbcost = 0.0
        koff, voff = 0,0
        for obox, ocount, icount, payload in zip(self.outboxes, self.outcounts, self.incounts, self.incache):
            oldvoff = voff
            val, voff = foo(self.inbuf, icount, voff, self.spec.payload)
            val = '%s%s' % (val, payload)

            if self.spec.outcoords == Spec.COORD_ONE:
                for i in xrange(ocount):
                    key, koff = foo(self.outbuf, 1, koff, self.spec.outcoords)
                    start = time.time()
                    self.bdb[key] = val
                    bdbcost += time.time() - start

            else:
                obox = (obox[0][0], obox[0][1], obox[1][0], obox[1][1])
                key, koff = foo(self.outbuf, ocount, koff, self.spec.outcoords)
                idxkey = 'b:%d' % self.nbdbitems
                start = time.time()
                self.bdb[idxkey] = key
                self.bdb[key] = val
                self.outidx.set(obox, self.nbdbitems)
                bdbcost += time.time() - start
                self.nbdbitems += 1
        self.inc_stat('bdbcost', bdbcost)
        self.reset_cache()

    def close(self):
        self.flush()
        self.outbuf = self.inbuf = None
        super(PStore2, self).close()
            
class PStore3(DiskStore):
    def __init__(self, op, run_id, f, strat):
        super(PStore3, self).__init__(op, run_id, f, strat)
        self.outcache = None
        self.nbdbitems = 0
        self.outbuf = None
        self.inbufs = [None] * self.nargs
        self.mergebuf = None


        self.memdb = {}
        self.memsize = 0
        

    def reset_cache(self):
        self.outcache = []
        self.outcounts = []
        self.outboxes = []
        self.incache = [[] for n in xrange(self.nargs)]
        self.incounts = [[] for n in xrange(self.nargs)]
        if self.outbuf == None:
            self.outbuf = None
            self.inbufs = [None] * self.nargs

    def merge_strings(self, *strs):
        size = sum(map(len, strs))
        if self.mergebuf == None or size > len(self.mergebuf):
            #print "allocating new merge buffer", (size * 2)
            self.mergebuf = create_string_buffer(size * 2)
        fmt = ''.join(('%ds' % len(s) for s in strs))
        struct.pack_into(fmt, self.mergebuf, 0, *strs)
        return self.mergebuf[:size]


    def uses_mode(self, mode):
        return mode & Mode.PTR != 0

    @instrument
    def update_stats(self, outcoords, *incoords_arr):
        outbox = bbox(outcoords)
        inboxes = map(bbox, incoords_arr)

        diff = lambda box: map(lambda p: p[1]-p[0], zip(box[0], box[1]))
        outarea = reduce(mul, diff(outbox))
        inareas = map(lambda box: reduce(mul, diff(box)), inboxes)

        fanins = map(len, incoords_arr)
        noutcells = len(outcoords)

        self.outarea += outarea
        self.inareas = map(sum, zip(self.inareas, inareas))
        self.fanins = map(sum, zip(self.fanins, map(lambda x: noutcells * x, fanins)))
        self.noutcells += noutcells
        self.ncalls += 1


    def add_to_cache(self, cache, counts, coords, mode ):
        if Spec.BOX == mode:
            box = bbox(coords)
            cache.extend(box)
            counts.append(2)
        elif Spec.GRID == mode:
            grid = gengrid(coords)
            cache.extend( grid[0] )
            cache.append( (0, len(grid[1])) )
            cache.extend( grid[1] )
            counts.append(2 + 1 + len(grid[1]))
        else:
            if mode in (Spec.COORD_MANY, Spec.KEY):
                cache.append((0, len(coords)))
            cache.extend(coords)
            counts.append(len(coords))
        
    @instrument
    def write(self, outcoords, *incoords_arr):
        minus, plus = 0,0        
        if self.spec.outcoords == Spec.COORD_ONE:

            for outcoord in outcoords:
                outcoord = (outcoord,)
                if outcoord not in self.memdb:
                    self.memdb[outcoord] = [set() for i in xrange(self.nargs)]
                    plus += 1
                for s, incoords in zip(self.memdb[outcoord], incoords_arr):
                    minus += len(s)
                    s.update(incoords)
                    plus += len(s)
        else:
            outcoords = tuple(outcoords)
            if outcoords not in self.memdb:
                self.memdb[outcoords] = [set() for i in xrange(self.nargs)]
                plus += len(outcoords)
            for s, incoords in zip(self.memdb[outcoords], incoords_arr):
                minus += len(s)
                s.update(incoords)
                plus += len(s)
        self.memsize = self.memsize - minus + plus

        if self.memsize > 500000:
            for key, val in self.memdb.iteritems():
                self.dump(key, *val)
            self.memdb = {}
            self.memsize = 0

    def dump(self, outcoords, *incoords_arr):
        
        #self.write_old(outcoords, *incoords_arr)
        #self.update_stats(outcoords, *incoords_arr)
        if self.outcache is None:
            self.reset_cache()
        start = time.time()
        box = bbox(outcoords)
        self.outboxes.append(box)
        self.add_to_cache(self.outcache, self.outcounts, outcoords, self.spec.outcoords)
        self.inc_stat('outcache', time.time() - start)

        start = time.time()
        for cache, counts, incoords in zip(self.incache, self.incounts, incoords_arr):
            self.add_to_cache(cache, counts, list(incoords), self.spec.payload)
        self.inc_stat('incache', time.time() - start)

        if sum(map(len,self.incache)) + len(self.outcache) > 1000000:
            self.flush()

    @instrument
    def flush(self):
        if not len(self.outcache): return
        KEYLEN = struct.pack("I", 14)
    
        start = time.time()
        oencs = [self.enc_out(outcoord) for outcoord in self.outcache]
        # serialize outputs into preallocated buffers
        fmt = self._serialize_format(self.outcounts, self.spec.outcoords)
        # lookup or create key buffer        
        keysize = struct.calcsize(fmt)
        if self.outbuf is None or keysize > len(self.outbuf):
            self.outbuf = create_string_buffer(keysize)
        struct.pack_into(fmt, self.outbuf, 0, *oencs)
        self.inc_stat('serout', time.time() - start)

        
        start = time.time()
        iencs = [ [ self.enc_in(incoord, arridx) for incoord in incoords ]
                  for arridx, incoords in enumerate(self.incache) ]
        # serialize the inputs into preallocated buffers
        for arridx, (counts, valbuf, incoords) in  enumerate(zip(self.incounts, self.inbufs, self.incache)):
            iencs = [ self.enc_in(incoord, arridx) for incoord in incoords ]
            fmt = self._serialize_format(counts, self.spec.payload)
            # lookup or create value buffer
            valsize = struct.calcsize(fmt)
            if valbuf is None or valsize > len(valbuf):
                self.inbufs[arridx] = create_string_buffer(valsize)
            struct.pack_into(fmt, self.inbufs[arridx], 0, *iencs)
        self.inc_stat('serin', time.time() - start)
        

        def foo(buf, count, offset, mode):
            if mode in (Spec.COORD_MANY, Spec.KEY):
                end = offset + 4 * (1 + count)
            elif mode == Spec.GRID:
                end = offset + 4 * count
            elif mode == Spec.BOX:
                end = offset + 8
            elif mode == Spec.COORD_ONE:
                end = offset + 4
            ret = buf[offset:end]
            
            
            if mode == Spec.KEY:
                start = time.time()
                h = str(hash(ret) % 4294967296).rjust(10, '_')
                key = 'key:%s' % h
                refkey = 'ref:%s' % h

                if key not in self.bdb:
                    self.bdb[key] = ret
                    self.bdb[refkey] = '0'
                ret = '%s%s' % (KEYLEN, key)
                self.inc_stat('keycost', time.time()-start)
                return ret, end
            
            return ret, end

        def merge_serialized(old, new, arridx, mode = None):
            if mode == None:
                mode = self.spec.payload
            if mode == Spec.KEY:
                if old == new: return old
                oldkey = old[4:]
                newkey = new[4:]
                oldval = self.bdb[oldkey]
                newval = self.bdb[newkey]

                newval = merge_serialized(oldval, newval, arridx, Spec.COORD_MANY)

                # insert the new key
                h = str(hash(newval) % 4294967296).rjust(10, '_')
                key = 'key:%s' % h
                refkey = 'ref:%s' % h
                if key not in self.bdb:
                    self.bdb[key] = newval
                    self.bdb[refkey] = '0'

                # we know we are removing references now
                # update ref counts and garbage collect
                oldref = 'ref:%s' % oldkey[4:]
                oldrefcount = int(self.bdb[oldref])            
                if oldrefcount <= 1:
                    del self.bdb[oldkey]
                    del self.bdb[oldref]
                else:
                    self.bdb[oldref] = str(oldrefcount - 1)

                newref = 'ref:%s' % newkey[4:]
                newrefcount = int(self.bdb[newref])
                if newrefcount <= 0:
                    del self.bdb[newkey]
                    del self.bdb[newref]
                else:
                    self.bdb[newref] = str(newrefcount - 1)


                return '%s%s' % (KEYLEN, key)

            elif mode == Spec.COORD_MANY:
                oldn, = struct.unpack('I', old[:4])
                newn, = struct.unpack('I', new[:4])
                n = struct.pack('I', oldn + newn)
                #return ''.join([n, old[4:], new[4:]])
                return self.merge_strings(n, old[4:], new[4:])
            elif mode == Spec.GRID:
                oldtmp = self._parse(StringIO(old), mode)
                newtmp = self._parse(StringIO(new), mode)
                oldbox, oldencs = oldtmp[:2], oldtmp[3:]
                newbox, newencs = newtmp[:2], newtmp[3:]
                oldbox = [self.dec_in(o, arridx) for o in oldbox]
                newbox = [self.dec_in(o, arridx) for o in newbox]
                coords = []
                coords.extend(decgrid(oldbox, map(lambda c: self.dec_in(c, arridx), oldencs)))
                coords.extend(decgrid(newbox, map(lambda c: self.dec_in(c, arridx), newencs)))
                grid = gengrid(coords)
                l = []
                l.extend( grid[0] )
                l.append( (0, len(grid[1]) ) )
                l.extend( grid[1] )
                return struct.pack('%dI' % len(l), *[self.enc_in(c, arridx) for c in l])
                          
            elif mode == Spec.BOX:
                raise RuntimeError
            elif mode == Spec.COORD_ONE:
                raise RuntimeError

        def segment_serinputs(val):
            mode = self.spec.payload
            if mode == Spec.KEY:
                off = 0
                ret = []
                for i in xrange(self.nargs):
                    ret.append(val[off:off+4+14])
                    off += 4 + 14
                return ret
            elif mode == Spec.COORD_MANY:
                off = 0
                ret = []
                for i in xrange(self.nargs):
                    n, = struct.unpack("I", val[off:off+4])
                    ret.append(val[off:off+4+n*4])
                    off += 4 + 4 * n
                return ret
            elif mode == Spec.GRID:
                off = 0
                ret = []
                for i in xrange(self.nargs):
                    n, = struct.unpack("I", val[off+8:off+8+4])
                    ret.append(val[off:off+8+4*(n+1)])
                    off += 8 + 4 * (n+1)
                return ret
            elif mode == Spec.BOX:
                off = 0
                ret = []
                for i in xrange(self.nargs):
                    ret.append(val[off:off+8])
                    off += 8
                return ret
            elif mode == Spec.COORD_ONE:
                raise RuntimeError

        def add_ptr(key, vals):
            """
            @param vals: array of bytes
            """
            start = time.time()
            if self.spec.payload == Spec.KEY:
                for val in vals:
                    # update the reference
                    refkey = 'ref:%s' % val[8:]
                    self.bdb[refkey] = str(int(self.bdb[refkey])+1)

            self.bdb[key] = ''.join(vals)
            return time.time()-start                        

                
        bdbcost = 0.0
        idxcost = 0.0
        mergecost = 0.0
        keyoffset = 0
        valoffsets = [0] * len(self.inbufs)
        #print len(self.inbufs), len(self.outboxes), len(self.outcounts), map(len, self.incounts)

        for obox, ocount,  icounts in zip(self.outboxes, self.outcounts, zip(*self.incounts)):
            
            obox = (obox[0][0], obox[0][1], obox[1][0], obox[1][1])
            vals = []
            for idx in xrange(len(valoffsets)):
                count, buf, offset = icounts[idx], self.inbufs[idx], valoffsets[idx]
                ba, offset = foo(buf, count, offset, self.spec.payload)
                valoffsets[idx] = offset
                vals.append(ba)
                

            # presume that we are creating references to each of the values
            # we will remove referencs later
            if self.spec.payload == Spec.KEY:
                nrefs = self.spec.outcoords == Spec.COORD_ONE and ocount or 1
                for val in vals:
                    refkey = 'ref:%s' % val[8:]
                    self.bdb[refkey] = str(int(self.bdb[refkey])+nrefs)

            nkeysremoved = 0
            if self.spec.outcoords == Spec.COORD_ONE:
                
                for i in xrange(ocount):
                    key, keyoffset = foo(self.outbuf, 1, keyoffset, self.spec.outcoords)
                    dec, = struct.unpack('I', key)
                    if key not in self.bdb:
                        bdbcost += add_ptr(key, vals)
                    else:
                        start = time.time()
                        oldvals = segment_serinputs(self.bdb[key])
                        curvals = vals
                        newvals = []
                        
                        for arridx, (oldval, curval) in enumerate(zip(oldvals, curvals)):
                            newvals.append(merge_serialized(oldval, curval, arridx))
                        mergecost += time.time() - start

                        bdbcost += add_ptr(key, newvals)

                        if self.spec.payload == Spec.KEY:
                            nkeysremoved += 1


            else:
                key, keyoffset = foo(self.outbuf, ocount, keyoffset, self.spec.outcoords)
                if key not in self.bdb:
                    idxkey = 'b:%d' % self.nbdbitems
                    start = time.time()
                    self.bdb[idxkey] = key                
                    bdbcost += time.time()-start

                    bdbcost += add_ptr(key, vals)

                    start = time.time()
                    self.outidx.set(obox, self.nbdbitems)
                    idxcost += time.time()-start
                    self.nbdbitems += 1
                else:
                    start = time.time()
                    oldvals = segment_serinputs(self.bdb[key])
                    curvals = vals
                    newvals = []
                    for arridx, (oldval, curval) in enumerate(zip(oldvals, curvals)):
                        curvals.append(merge_serialized(oldval, curval, arridx))
                    mergecost += time.time() - start

                    bdbcost += add_ptr(key, curvals)
                    nkeysremoved += 1



        
        self.inc_stat('bdbcost', bdbcost)
        self.inc_stat('idxcost', idxcost)
        self.inc_stat('mergecost', mergecost)        
        # reset the cache
        self.reset_cache()                


    def close(self):
        if self.memdb:
            for key, val in self.memdb.iteritems():
                self.dump(key, *val)
        self.memdb = None
        self.flush()
        self.outbuf = None
        self.inbufs = [None] * self.nargs
        self.mergebuf = None
        print "Closing", self.op, '\t', len(self.bdb), " entries"
        super(PStore3, self).close()
        



class IBox(object):
    def extract_incells(self, boxes, arridx, outcoords, inputs):

        newinputs = []
        for curarridx, curbox in enumerate(boxes):
            tmparr, tmptomap, tmpfrommap = subarray(inputs[curarridx], curbox)
            if arridx == curarridx:
                newarr, tonewcoords, fromnewcoords = tmparr, tmptomap, tmpfrommap
            newinputs.append(tmparr)
        newoutcoords = map(tonewcoords, outcoords)

        results = self.re_execute(newinputs, newoutcoords, arridx, newarr.shape, True)
        for coord in map(fromnewcoords, results):
            yield coord

        
    def extract_outcells(self, boxes, arridx, incoords, inputs):

        newinputs = []
        for curarridx, curbox in enumerate(boxes):
            tmparr, tmptomap, tmpfrommap = subarray(inputs[curarridx], curbox)    
            if arridx == curarridx:
                newarr, tonewcoords, fromnewcoords = tmparr, tmptomap, tmpfrommap
            newinputs.append(tmparr)
        newincoords = map(tonewcoords, incoords)

        results = self.re_execute(newinputs, newincoords, arridx, newarr.shape, False)
        for coord in map(fromnewcoords, results):
            yield coord


    def re_execute(self, inputs, coords, arridx, outshape, backward=True):
        from runtime import Runtime
        ptype = backward and BReexec or FReexec

        pqres = PQSparseResult(outshape)
        reexec = ptype(coords, arridx, pqres)

        Runtime.instance().set_reexec(self.op, self.run_id, reexec)
        self.op.run(inputs, self.run_id)
        Runtime.instance().clear_reexec(self.op, self.run_id)
        
        for coord in pqres:
            yield coord

    @alltoall
    def join(self, left, arridx, backward=True):
        self.open()
        if backward and self.spec.outcoords == Spec.COORD_ONE:
            for coord in self.hash_join(left, arridx):
                yield coord
            self.close()
            return
        
        inputs = self.op.wrapper.get_inputs(self.run_id)
        nfound = 0

        if backward:
            if Spec.BOX == self.spec.outcoords:
                pred = lambda obj: self._parse(StringIO(obj[0]), Spec.BOX)
                getkey = lambda item: struct.pack('4I', *item.bbox)
                matches = in_box
            else:
                pred = lambda obj: super(IBox,self).extract_outcells(obj)
                getkey = lambda item: self.bdb[item.object]
                matches = lambda coord, outcoords: coord in outcoords

            bigboxes = [[] for i in xrange(self.nargs)]

            for l in left:
                l = tuple(l)

                for item in self.outidx.get_pt(l):
                    key =  getkey(item)
                    r = (key, self.bdb[key])
                    coords = pred(r)
                    b = matches(l, coords)

                    if b:
                        boxes = self.extract_inboxes(r)
                        for i in xrange(self.nargs):
                            bigboxes[i].extend(boxes[i])

            boxes = [bbox(bigboxes[i]) for i in xrange(self.nargs)]
            for coord in self.extract_incells(boxes, arridx, left, inputs):
                yield coord

            
            
        else:
            def pred(coord, obj):
                if Spec.BOX == self.spec.payload:
                    buf = StringIO(obj[1])
                    buf.seek(4*4*arridx)
                    box = map(lambda b: self.dec_in(b, arridx), self._parse(buf, Spec.BOX))
                    return in_box(coord, box)
                else:
                    return tuple(coord) in super(IBox,self).extract_incells(obj, arridx)

            bigboxes = [[] for i in xrange(self.nargs)]
            for r in self.get_iter():
                for l in left:
                    if pred(l,r):
                        nfound += 1
                        boxes = self.extract_inboxes(r)
                        for i in xrange(self.nargs):
                            bigboxes[i].extend(boxes[i])

            boxes = [bbox(bigboxes[i]) for i in xrange(self.nargs) ]
            for coord in self.extract_outcells(boxes, arridx, left, inputs):
                yield coord

        self.close()


    def hash_join(self, left, arridx):
        bigboxes = [[] for i in xrange(self.nargs)]        
        buf = []
        serbuf = create_string_buffer(1000*4)
        for l in left:
            buf.append(l)

            if len(buf) >= 1000:
                encs = map(self.enc_out, buf)
                struct.pack_into('%dI' % len(encs), serbuf, 0, *encs)
                for i in xrange(1000):
                    key = serbuf[i*4:i*4+4]
                    if key not in self.bdb: continue
                    payload = self.bdb[key]
                    if payload != None:
                        boxes = self.extract_inboxes((key, payload))
                        for i in xrange(self.nargs):
                            bigboxes[i].extend(boxes[i])
                buf = []

        if len(buf) > 0:
            encs = map(self.enc_out, buf)
            struct.pack_into('%dI' % len(encs), serbuf, 0, *encs)
            for i in xrange(len(buf)):
                key = serbuf[i*4:i*4+4]
                if key not in self.bdb: continue
                payload = self.bdb[key]
                if payload != None:
                    boxes = self.extract_inboxes((key, payload))
                    for i in xrange(self.nargs):
                        bigboxes[i].extend(boxes[i])
            buf = []

        inputs = self.op.wrapper.get_inputs(self.run_id)
        boxes = [bbox(bigboxes[i]) for i in xrange(self.nargs)]
        for coord in self.extract_incells(boxes, arridx, left, inputs):
            yield coord


class PStore3Box(IBox,PStore3):
    def __init__(self, *args, **kwargs):
        super(PStore3Box, self).__init__(*args, **kwargs)

    def uses_mode(self, mode):
        return mode & Mode.PTR  != 0




class PStore2Box(IBox, PStore2):
    def __init__(self, *args, **kwargs):
        super(PStore2Box, self).__init__(*args, **kwargs)

    def uses_mode(self, mode):
        return mode & Mode.PT_MAPFUNC_BOX != 0

class PStoreQuery(IBox,PStore1):
    def __init__(self, *args, **kwargs):
        super(PStoreQuery, self).__init__(*args)

    def uses_mode(self, mode):
        return mode & Mode.QUERY != 0

    @instrument
    def write(self, outcoords, *incoords_arr):
        pass

    @alltoall
    def join(self, left, arridx, backward=True):
        self.open()
        inputs = self.op.wrapper.get_inputs(self.run_id)
        for coord in self.re_execute(inputs, left, arridx, self.outshape, backward=backward):
            yield coord
        self.close()

def create_ftype(ptype):
    class FStore2(IPstore):
        def __init__(self, op, run_id, fname, strat):
            super(FStore2, self).__init__(op, run_id, None, strat)

            self.pstores = []
            for arridx in xrange(self.nargs):
                newfname = '%s_%d' % (fname, arridx)
                pstore = ptype(op, run_id, newfname, strat)
                
                pstore.inshapes = [self.outshape]
                pstore.outshape = self.inshapes[arridx]
                pstore.nargs = 1
                if isinstance(pstore, PStore3):
                    pstore.reset_cache()
                self.pstores.append(pstore)

        def uses_mode(self, mode):
            for pstore in self.pstores:
                if pstore.uses_mode(mode):
                    return True
            return False

        def clear_stats(self):
            self.stats = {}
            for pstore in self.pstores:
                pstore.clear_stats()
            
        def get_stat(self, attr, default=0):
            if not len(self.pstores):
                return default
            return sum([pstore.get_stat(attr, 0) for pstore in self.pstores])

        def disk(self):
            return sum([pstore.disk() for pstore in self.pstores])

        def indexsize(self):
            return sum([pstore.indexsize() for pstore in self.pstores])

        @instrument
        def write(self, outcoords, *incoords_arr):
            for pstore, incoords in zip(self.pstores, incoords_arr):
                pstore.write(incoords, outcoords)

        @alltoall
        def join(self, left, arridx, backward=True):
            return self.pstores[arridx].join(left, 0, not backward)

        def open(self, new=False):
            for pstore in self.pstores:
                pstore.open(new=new)

        def close(self):
            for pstore in self.pstores:
                pstore.close()




    class FStore(IPstore):
        def __init__(self, op, run_id, fname, strat):
            super(FStore, self).__init__(op, run_id, None, strat)

            # newstrat = Desc(strat.mode, strat.spec, strat.backward)
            # newstrat.spec.outcoords = strat.spec.payload
            # newstrat.spec.payload = strat.spec.outcoords
            
            self.pstore = ptype(op, run_id, fname, strat)
            # fix some internal data structures
            # the internal pstore's output shape must encode all possible
            # values of the outer's input shapes
            #  - (max(ncells in each input array), narrays)
            self.pstore.inshapes = [self.pstore.outshape]
            self.pstore.outshape = (max(map(lambda shape: reduce(mul, shape), self.inshapes)),
                                    len(self.inshapes))
            self.pstore.nargs = 1


        def disk(self):
            return self.pstore.disk()

        def indexsize(self):
            return self.pstore.indexsize()
        
        def compress_incoords_arr(self, incoords_arr):
            ret = []
            for arridx, incoords in enumerate(incoords_arr):
                shape = self.inshapes[arridx]
                newincoords = map(lambda coord: (arridx, enc(coord, shape)), incoords)
                ret.extend(newincoords)
            return ret

        @instrument
        def write(self, outcoords, *incoords_arr):
            newincoords = self.compress_incoords_arr(incoords_arr)
            self.pstore.write(newincoords, outcoords)

        @alltoall
        def join(self, left, arridx, backward=True):
            shape = self.inshapes[arridx]

            if backward:
                for coord in self.pstore.join(left, arridx, not backward):
                    if coord[1] != arridx: continue
                    yield dec(coord[0], shape)
            else:

                left = [(arridx, enc(coord, shape)) for coord in left]
                for coord in self.pstore.join(left, arridx, not backward):
                    yield coord

        def open(self, new=False):
            self.pstore.open(new=new)

        def close(self):
            self.pstore.close()

    return FStore2

class MultiPStore(IPstore):
    def __init__(self, op, run_id, fname, strat, bpstore, fpstore):
        super(MultiPStore,self).__init__(op, run_id, fname, strat)
        self.bpstore = bpstore
        self.fpstore = fpstore
        self.pt_stores = []
        self.ptr_stores = []
        self.pstores = []

        for pstore in [bpstore, fpstore]:
            if pstore.uses_mode(Mode.PT_MAPFUNC):
                self.pt_stores.append(pstore)
            if pstore.uses_mode(Mode.PTR):
                self.ptr_stores.append(pstore)
            self.pstores.append(pstore)

    def uses_mode(self, mode):
        return self.bpstore.uses_mode(mode) or self.fpstore.uses_mode(mode)

    def clear_stats(self):
        self.stats = {}
        self.bpstore.clear_stats()
        self.fpstore.clear_stats()

    def get_stat(self, attr, default=0):
        if not len(self.pstores):
            return default
        return sum([pstore.get_stat(attr, 0) for pstore in (self.bpstore, self.fpstore)])


    def disk(self):
        return self.bpstore.disk() + self.fpstore.disk()

    def indexsize(self):
        return self.bpstore.indexsize() + self.fpstore.indexsize()

    @instrument
    def write(self, outcoords, *incoords_arr):
        if len(incoords_arr) == 1 and isinstance(incoords_arr[0], str):
            for pstore in self.pt_stores:
                pstore.write(outcoords, *incoords_arr)
        else:
            for pstore in self.ptr_stores:
                pstore.write(outcoords, *incoords_arr)

    @alltoall
    def join(self, left, arridx, backward=True):
        if backward:
            return self.bpstore.join(left, arridx, backward=backward)
        return self.fpstore.join(left, arridx, backward=backward)

    def close(self):
        self.bpstore.close()
        self.fpstore.close()
    
        


class FReexec(IPstore):

    def __init__(self, incoords, arridx, pqres):
        self.incoords = set(map(tuple,incoords))
        self.arridx = arridx
        self.pqres = pqres
        self.strat = Desc(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False)
        self.spec = self.strat.spec
        self.stats = {}

    def uses_mode(self, mode):
        return mode & Mode.PTR != 0

    @instrument
    def write(self, outcoords, *incoords_arr):
        if len(self.incoords.intersection(incoords_arr[self.arridx])) > 0:
            for outcoord in outcoords:
                self.pqres.add(outcoord)


class BReexec(IPstore):
    def __init__(self, outcoords, arridx, pqres):
        self.outcoords = set(map(tuple, outcoords))
        self.arridx = arridx
        self.pqres = pqres
        self.strat = Desc(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), True)
        self.spec = self.strat.spec
        self.stats = {}
        
    def uses_mode(self, mode):
        return mode & Mode.PTR != 0

    @instrument
    def write(self, outcoords, *incoords_arr):
        if len(self.outcoords.intersection(outcoords)) > 0:
            for incoord in incoords_arr[self.arridx]:
                self.pqres.add(incoord)
        



class CompositePStore(IPstore):
    def __init__(self, op, run_id, fname, strat, pstores):
        super(CompositePStore, self).__init__(op, run_id, fname, strat)
        self.spec = strat.spec
        self.pstore1 = None
        self.pstore2 = None
        self.pstore3 = None        
        self.pstores = pstores

        self.strat.mode = reduce(or_, map(lambda ps: ps.strat.mode, pstores))
        
        for pstore in pstores:
            if isinstance(pstore, PStore1):
                self.pstore1 = pstore
            elif isinstance(pstore, PStore2):
                self.pstore2 = pstore
            elif isinstance(pstore, PStore3):
                self.pstore3 = pstore

    def uses_mode(self, mode):
        for pstore in self.pstores:
            if pstore.uses_mode(mode):
                return True
        return False

    def clear_stats(self):
        for pstore in pstores:
            pstore.clear_stats()
        self.stats = {}


    def get_stat(self, attr, default=0):
        if not len(self.pstores):
            return default
        return sum([pstore.get_stat(attr, 0) for pstore in self.pstores])


    def disk(self):
        return sum([pstore.disk() for pstore in self.pstores])

    def indexsize(self):
        return sum([pstore.indexsize() for pstore in self.pstores])

    @instrument
    def update_stats(self, *args):
        super(CompositePStore, self).update_stats(*args)

    @instrument
    def write(self, outcoords, *args):
        if len(args) == 1 and type(args[0]) == str:
            self.pstore2.write(outcoords, args[0])
        else:
            self.pstore3.write(outcoords, *args)

    @alltoall
    def join(self, left, arridx, backward=True):
        self.open()
        for pstore in self.pstores:
            for coord in pstore.join(left, arridx, backward=backward):
                yield coord
        self.close()

    def open(self, new=False):
        for pstore in self.pstores:
            pstore.open(new=new)

    def close(self):
        for pstore in self.pstores:
            pstore.close()
        




if __name__ == '__main__':
    # coord = 'aaaaaaaaaaaa'

    # def run(npoints):
    #     idx = SpatialIndex('/tmp/test')
    #     idx.open(True)
    #     start = time.time()
    #     for i in xrange(npoints):
    #         idx.set((i,i,i+1,i+1), coord)
    #     end1 = time.time()
    #     idx.close()
    #     end2 = time.time()
    #     return end1-start, end2-start

    # def bdb(npoints):
    #     db = bsddb.hashopen('/tmp/test.db', 'n')
    #     start = time.time()
    #     for i in xrange(npoints):
    #         db[str(i)] = coord
    #     end1 = time.time()
    #     db.close()
    #     end2 = time.time()
    #     return end1-start, end2-start

    # def gen(bsize):
    #     for i in xrange(10000):
    #         x,y = random.randint(0, 99), random.randint(0, 99)
    #         yield (i, (x,y,x+bsize, y+bsize), 'b:%d' % i)

    # from rtree import Rtree
    # for nboxes in (1, 100, 1000):
    #     for bsize in (1, 5, 6, 10):
    #         idx = Rtree(gen(bsize))
    #         nres = 0
    #         start = time.time()
    #         for i in xrange(1000):
    #             pt = (random.randint(0, 10),random.randint(0, 10))
    #             for x in idx.intersection(pt, objects=True):
    #                 nres += 1
    #         cost = time.time() - start
    #         print '%d\t%d\t%d\t' % (nboxes, bsize, nres/1000.0), cost / 1000.0
    # exit()


    # for i in (10, 100, 1000, 10000, 100000):
    #     cost1, cost2 = run(i)
    #     datsize = os.path.getsize('/tmp/test_rtree.dat')
    #     idxsize = os.path.getsize('/tmp/test_rtree.idx')
    #     bdbsize = os.path.getsize('/tmp/test_rtree.bdb')
    #     print '%d\t%f\t%f\t%f\t%f\t%f\t%f' % (i, cost2 / i, cost1, cost2,
    #                                           datsize/1048576.0, idxsize / 1048576.0, bdbsize/1048576.0)
    # exit()

    # idx.set((0,0,1,1), '0')
    # idx.close()

    
    # idx.open(False)
    # print map(str, idx.get_pt((0,0)))
    # idx.close()
    # idx.open(False)
    # print map(str, idx.get_pt((0,0)))
    # idx.close()

    # idx.open(True)
    # print map(str, idx.get_pt((0,0)))
    # idx.close()
    # exit()


    all_coords = [ (i,i) for i in xrange(1000) ]
    def ser(pstore, coords, spec):
        buf = StringIO()        
        if spec == Spec.BOX:
            encs = bbox(coords)
        else:
            encs = map(pstore.enc_out, coords)

        if spec == Spec.COORD_ONE:
            for enc in encs:
                buf.seek(0)
                pstore._serialize((enc,), buf, spec)
        else:
            pstore._serialize(encs, buf, spec)
        return buf.getvalue()
        
    def par(pstore, s, spec):
        if spec == Spec.BOX:
            pstore._parse(StringIO(s), spec)
        elif spec == Spec.KEY:
            pstore._parse(StringIO(s), spec)
        else:
            encs = pstore._parse(StringIO(s), spec)
            if encs:
                map(pstore.dec_out, encs)


    for spec in Spec.all():
        for coords in all_coords:
            if spec == Spec.BINARY: continue
            s = ser(pstore, coords, spec)

            f = lambda: ser(pstore, coords, spec)
            h = lambda: par(pstore, s, spec)

            t = timeit.Timer(f)
            t.timeit(100)
            cost = t.timeit(1000) / 1000.0
            print '%d\t%d\t' % (len(coords), spec), float(cost) / len(coords)
            
            
        

    
    
