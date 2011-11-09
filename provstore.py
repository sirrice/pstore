import struct, math, time, random
import bsddb3 as bsddb
from StringIO import StringIO
from operator import mul, or_, and_

from strat import *
from queryresult import *
from util import subarray
from rtree import index


plog = logging.getLogger('provstore')
logging.basicConfig()
plog.setLevel(logging.INFO)


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
    minc, maxc = None, None
    for coord in coords:
        if minc == None:
            minc = coord
            maxc = coord
        else:
            minc = map(min, zip(minc, coord))
            maxc = map(max, zip(maxc, coord))

    return (minc, maxc)


__grid__ = None
__gcells__ = 0
def gengrid(coords):
    """
    @return (box, negs)
    box: bounding box containing the coordinates
    negs: encoded list of coords not in bounding box
    """
    global __grid__, __gcells__
    box = bbox(coords)
    shape = (1+box[1][0]-box[0][0], 1+box[1][1]-box[0][1])
    ncells = reduce(mul, shape)
    if __gcells__ < ncells:
        __grid__ = np.ones((1, ncells), dtype=bool)
        __gcells__ = reduce(mul, __grid__.shape)
    arr = __grid__[:,:ncells].reshape((shape))
    arr[:,:] = True
    newcoords = map(lambda coord: (coord[0]-box[0][0], coord[1]-box[0][1]), coords)
    arr[zip(*newcoords)] = False
    encs = map(lambda coord: enc(coord, shape), np.argwhere(arr))
    return (box, encs)

def decgrid(box, negs):
    global __grid__, __gcells__    
    shape = ( 1+box[1][0]-box[0][0], 1+box[1][1]-box[0][1] )
    ncells = reduce(mul, shape)
    if __gcells__ < ncells:
        __grid__ = np.ones((1, ncells), dtype=bool)
        __gcells__ = reduce(mul, __grid__.shape)
    arr = __grid__[:,:ncells].reshape((shape))
    arr[:,:] = True
    arr[zip(*map(lambda enc: dec(enc, shape), negs))] = False
    return map(lambda coord: (coord[0]+box[0][0], coord[1]+box[0][1]), np.argwhere(arr))

    

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
        return mode in self.strat.modes()

    def get_iter(self):
        return ()

    def disk(self):
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
        return mode in ( Mode.PTR, Mode.FULL_MAPFUNC )
    
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
        return obj == Mode.FULL_MAPFUNC

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
        return (x for x in self.bdb.iteritems() if not x[0].startswith("key:"))

    def disk(self):
        try:
            return os.path.getsize(self.fname)
        except:
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
            box, negs = obj[0]
            return map(self.enc_out, decgrid(box, negs))
            
        ret = self._parse(StringIO(obj[0]), self.spec.outcoords)
        return ret

    def extract_outcells(self, obj):
        if self.spec.outcoords == Spec.GRID:
            return decgrid(obj[0][0], obj[0][1])
        return map(self.dec_out, self.extract_outcells_enc(obj))


    def extract_outboxes(self, obj):
        if self.spec.outcoords == Spec.BOX:
            return self._parse(StringIO(obj[0]), self.spec.outcoords)
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
            for box, negs in parsed:
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
            box, negs = val
            return map(lambda coord: self.enc_in(coord, arridx), decgrid(box, negs))
        elif Spec.is_coord(self.spec.payload):
            return val
        raise RuntimeError

    def extract_incells(self, obj, arridx):
        if Spec.GRID == self.spec.payload:
            buf = StringIO(obj[1])
            grid, negs = self._parse(buf, Spec.GRID)
            return decgrid(grid, negs)
        return map(lambda val: self.dec_in(val, arridx),
                   self.extract_incells_enc(obj, arridx))

    def extract_inboxes(self, obj):
        if self.spec.payload == Spec.BOX:
            buf = StringIO(obj[1])
            return [self._parse(buf, self.spec.payload) for arridx in  xrange(self.nargs)]
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
            buf.write(struct.pack("4I", box[0][0], box[0][1], box[1][0], box[1][1]))
        elif mode == Spec.GRID:
            box, encs = data
            self._serialize(box, buf, Spec.BOX)
            self._serialize(encs, buf, Spec.COORD_MANY)
        elif mode == Spec.KEY:
            ser_coords = StringIO()
            self._serialize(data, ser_coords, Spec.COORD_MANY)
            ser_coords = ser_coords.getvalue()
            ser_key = 'key:%s' % str(hash(ser_coords))     # poor man's key gen
            self.bdb[ser_key] = ser_coords
            s = struct.pack("I%ds" % (len(ser_key)), len(ser_key), ser_key)
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
            return (struct.unpack("2I", buf.read(8)), struct.unpack("2I", buf.read(8)))
            vals = struct.unpack("4I", buf.read(4*4))
            return (vals[:2],vals[2:])
        elif mode == Spec.GRID:
            return self._parse(buf, Spec.BOX), self._parse(buf, Spec.COORD_MANY)
        elif mode == Spec.KEY:
            n, = struct.unpack("I", buf.read(4))
            ser_key, = struct.unpack("%ds" % n, buf.read(n))
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


        
        predcost = 0.0
        matchcost = 0.0
        extractcost = 0.0
        niter = 0.0
        npred = 0.0
        nmatches = 0.0


        # pred
        # matches
        # extract
        if backward:  # backward query
            if Spec.BOX == self.spec.outcoords:
                pred = lambda obj: self._parse(StringIO(obj[0]), Spec.BOX)
                matches = in_box
                def getkey(item):
                    box = item.bbox
                    box = (box[:2],box[2:])
                    key = StringIO()
                    self._serialize(box, key, Spec.BOX)
                    return key.getvalue()
            else:
                pred = self.extract_outcells_enc
                matches = lambda coord, outcoords: self.enc_out(coord) in outcoords
                #left = map(self.enc_out, left)
                getkey = lambda item: item.object
            extract = lambda obj: self.extract_incells(obj, arridx)

            for l in left:
                l = tuple(l)
                for item in self.outidx.get_pt(l):
                    key = getkey(item)
                    r = (key, self.bdb[key])
                    coords = pred(r)
                    b = matches(l, coords)
                    if b:
                        start = time.time()
                        e = extract(r)
                        extractcost += time.time() - start
                        for coord in e:
                            yield coord
            

        else:  # forward query
            if Spec.BOX == self.spec.payload:
                def pred(coord, obj):
                    buf = StringIO(obj[1])
                    buf.seek(4*4*arridx)
                    box = self._parse(buf, Spec.BOX)
                    return box
                matches = in_box
            else:
                pred = lambda obj: self.extract_incells_enc(obj, arridx)
                matches = lambda coord, coords: coord in coords
                left = map(lambda l: self.enc_in(l, arridx), left)
            extract = lambda obj: self.extract_outcells(obj)



            for r in self.get_iter(): # r is a pair of unparsed strings
                niter += 1
                start = time.time()
                coords = pred(r) #
                predcost += time.time() - start
                npred += len(coords)

                for l in left:
                    start = time.time()
                    b = matches(l, coords)
                    matchcost += time.time() - start
                    nmatches += 1
                    if b:
                        start = time.time()
                        e = extract(r)
                        extractcost += time.time() - start
                        for coord in e:
                            yield coord
                        break
        self.close()
        
        plog.debug( "n/left:   %d\t%d", niter, len(left) )
        plog.debug( "predcost: %f\t%d", predcost, npred )
        plog.debug( "match:    %f\t%f\t%d", matchcost, matchcost / (1.0+nmatches), nmatches )
        plog.debug( "extract:  %f", extractcost )

    def hash_join(self, left, arridx):
        datacost = 0.0
        extractcost = 0.0
        sercost = 0.0
        niter = 1.0
        nfound = 1.0
        total = 0.0
        tstart = time.time()
        for l in left:
            niter += 1
            start = time.time()
            enc = (self.enc_out(l),)
            key = StringIO()
            self._serialize(enc, key, self.spec.outcoords)
            key = key.getvalue()
            sercost += time.time() - start
            if key not in self.bdb: continue
            nfound += 1
            start = time.time()
            payload = self.bdb[key]
            datacost += time.time() - start
            if payload:
                start = time.time()
                e = self.extract_incells((key, payload), arridx)
                extractcost += time.time() - start
                for coord in e:
                    yield coord
                    
        plog.debug( "datacost:    %f\t%f", datacost, datacost / nfound )
        plog.debug( "extractcost: %f\t%f", extractcost, extractcost / nfound )
        plog.debug( "sercost:     %f\t%f", sercost, sercost / niter )


    def open(self, new=False):
        if self.bdb != None: return
        if new:
            self.bdb = bsddb.hashopen(self.fname, 'n')
        else:
            self.bdb = bsddb.hashopen(self.fname)
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


        
class BinaryRTree(index.Index):
    def dumps(self, obj):
        return obj
    def loads(self, s):
        return s

class SpatialIndex(object):
    def __init__(self, fname):
        self.fname = '%s_rtree' % fname
        self.rtree = None
        self.id = 0
        p = index.Property()
        p.dimension = 2
        p.filename = self.fname
        self.p = p
        

    def set(self, box, val):
        self.rtree.add(self.id, box, val)
        self.id += 1

    def get_pt(self, pt):
        return self.rtree.intersection(pt, objects=True)

    def get_box(self, box):
        return self.rtree.intersection(box, objects=True)

    def open(self, new=False):
        p = self.p
        if new:
            try:
                os.unlink('%s.%s' % (p.get_filename(), p.get_idx_extension()))
                os.unlink('%s.%s' % (p.get_filename(), p.get_dat_extension()))
            except:
                pass
        self.rtree = BinaryRTree(p.get_filename(), properties=p)

    def close(self):
        self.rtree.close()

    def delete(self):
        p = self.p
        try:
            os.unlink('%s.%s' % (p.get_filename(), p.get_idx_extension()))
            os.unlink('%s.%s' % (p.get_filename(), p.get_dat_extension()))
        except:
            raise

        


class PStore2(DiskStore):
    def __init__(self, op, run_id, fname, strat):
        super(PStore2, self).__init__(op, run_id, fname, strat)
        self.spec.payload = Spec.BINARY # ignore payload spec

    def uses_mode(self, mode):
        return mode == Mode.PT_MAPFUNC

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

    @instrument
    def write(self, outcoords, payload):
        #self.update_stats(outcoords, payload)
        start = time.time()
        key = self.get_key(outcoords)
        box = []
        map(lambda x: box.extend(x), bbox(outcoords))
        self.outidx.set(box, key)
        self.inc_stat('serout', time.time() - start)


        
        start = time.time()
        val = self.get_val(payload)
        self.inc_stat('serin', time.time() - start)

        
        start = time.time()
        self.bdb[key] = val
        self.inc_stat('bdb', time.time() - start)

            
class PStore3(DiskStore):
    def __init__(self, op, run_id, f, strat):
        super(PStore3, self).__init__(op, run_id, f, strat)

    def uses_mode(self, mode):
        return mode == Mode.PTR

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
        
    @instrument
    def write(self, outcoords, *incoords_arr):
        #self.update_stats(outcoords, *incoords_arr)

        start = time.time()
        val = StringIO()
        if Spec.BOX == self.spec.payload:
            boxes = map(bbox, incoords_arr)
            for box in boxes:
                self._serialize(box, val, Spec.BOX)
        elif Spec.GRID == self.spec.payload:
            grids = map(gengrid, incoords_arr)
            for grid in grids:
                self._serialize(grid, val, Spec.GRID)
        else:
            for arridx, incoords in enumerate(incoords_arr):
                encs = map(lambda coord: self.enc_in(coord, arridx), incoords)
                self._serialize(encs, val, self.spec.payload)
        val = val.getvalue()
        self.inc_stat('serin', time.time() - start)


        ostart = time.time()
        if Spec.COORD_ONE == self.spec.outcoords:
            try:
                enc_outcoords = map(self.enc_out, outcoords)
            except:
                plog.error('could not encode outcoord\t%s', outcoords )
                raise
            key = StringIO()
            for enc in enc_outcoords:
                key.seek(0)
                self._serialize((enc,), key, self.spec.outcoords)
                keystr = key.getvalue()
                if keystr not in self.bdb:
                    start = time.time()
                    self.bdb[keystr] = val
                    self.inc_stat('bdb', time.time() - start)
                else:
                    # fuck, need to add provenance to existing data in the bdb
                    if Spec.BOX == self.spec.payload:
                        raise RuntimeError, "Don't support appending to existing provenance for BOX"

                    old_encs = self.extract_allincells_enc((None, self.bdb[keystr]))

                    # buf = StringIO(self.bdb[keystr])
                    # old_encs = []
                    # for idx in xrange(len(incoords_arr)):
                    #     old_encs.append(self._parse(buf, self.spec.payload))

                    newval = StringIO()
                    for incoords, old_enc in zip(incoords_arr, old_encs):
                        encs = map(lambda coord: self.enc_in(coord, arridx), incoords)
                        encs.extend(old_enc)
                        self._serialize(encs, newval, self.spec.payload)
                    newval = newval.getvalue()

                    start = time.time()
                    self.bdb[keystr] = newval
                    self.inc_stat('bdb', time.time() - start)
        elif Spec.BOX == self.spec.outcoords:
            key = StringIO()
            box = bbox(outcoords)
            self._serialize(box, key, Spec.BOX)
            key = key.getvalue()
            
            start = time.time()
            self.bdb[key] = val
            self.outidx.set((box[0][0], box[0][1], box[1][0], box[1][1]), None)
            self.inc_stat('bdb', time.time() - start)
        else:
            enc_outcoords = map(self.enc_out, outcoords)
            box = bbox(outcoords)
            key = StringIO()
            self._serialize(enc_outcoords, key, self.spec.outcoords)
            key = key.getvalue()
            
            start = time.time()
            self.bdb[key] = val
            self.outidx.set((box[0][0], box[0][1], box[1][0], box[1][1]), key)
            self.inc_stat('bdb', time.time() - start)
            
        self.inc_stat('serout', time.time() - ostart)



class IBox(object):
    def extract_incells(self, obj, arridx, outcoords, inputs):
        boxes = self.extract_inboxes(obj)

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

        
    def extract_outcells(self, obj, arridx, incoords, inputs):
        boxes = self.extract_inboxes(obj)

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
        

        if backward:
            def pred(coord, obj):
                if Spec.BOX == self.spec.outcoords:
                    box = self._parse(StringIO(obj[0]), Spec.BOX)
                    return in_box(coord, box)
                else:
                    return tuple(coord) in super(IBox,self).extract_outcells(obj)
            def extract(obj, coords, inputs):
                return self.extract_incells(obj, arridx, coords, inputs)
        else:

            def pred(coord, obj):
                if Spec.BOX == self.spec.payload:
                    buf = StringIO(obj[1])
                    buf.seek(4*4*arridx)
                    box = self._parse(buf, Spec.BOX)
                    return in_box(coord, box)
                else:
                    return tuple(coord) in super(IBox,self).extract_incells(obj, arridx)

            def extract(obj, coords, inputs):
                return self.extract_outcells(obj, arridx, coords, inputs)


        inputs = self.op.wrapper.get_inputs(self.run_id)
        nfound = 0
        for r in self.get_iter():
            for l in left:
                if pred(l,r):
                    nfound += 1
                    for coord in extract(r, left, inputs):
                        yield coord
                    break
        self.close()


    def hash_join(self, left, arridx):
        for l in left:
            enc = (self.enc_out(l),)
            key = StringIO()
            self._serialize(enc, key, self.spec.outcoords)
            key = key.getvalue()
            if key not in self.bdb: continue
            payload = self.bdb[key]
            if payload:
                for coord in self.extract_incells((key, payload), arridx):
                    yield coord

class PStore3Box(IBox,PStore3):
    def __init__(self, *args, **kwargs):
        super(PStore3Box, self).__init__(*args, **kwargs)

    def uses_mode(self, mode):
        return mode == Mode.PTR 




class PStore2Box(IBox, PStore2):
    def __init__(self, *args, **kwargs):
        super(PStore2Box, self).__init__(*args, **kwargs)

    def uses_mode(self, mode):
        return mode == Mode.PT_MAPFUNC_BOX

class PStoreQuery(IBox,PStore1):
    def __init__(self, *args, **kwargs):
        super(PStoreQuery, self).__init__(*args)

    def uses_mode(self, mode):
        return mode == Mode.QUERY

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

        for pstore in [bpstore, fpstore]:
            if pstore.uses_mode(Mode.PT_MAPFUNC):
                self.pt_stores.append(pstore)
            if pstore.uses_mode(Mode.PTR):
                self.ptr_stores.append(pstore)

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
    
        


class FReexec(IPstore):

    def __init__(self, incoords, arridx, pqres):
        self.incoords = set(map(tuple,incoords))
        self.arridx = arridx
        self.pqres = pqres
        self.strat = Desc(Mode.PTR, Spec(Spec.COORD_MANY, Spec.COORD_MANY), False)
        self.spec = self.strat.spec
        self.stats = {}

    def uses_mode(self, mode):
        return mode == Mode.PTR

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
        return mode == Mode.PTR

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
    coord = None

    def run(npoints):
        idx = SpatialIndex('/tmp/test')
        idx.open(True)
        start = time.time()
        for i in xrange(npoints):
            idx.set((i,i,i+1,i+1), coord)
        end1 = time.time()
        idx.close()
        end2 = time.time()
        return end1-start, end2-start

    def bdb(npoints):
        db = bsddb.hashopen('/tmp/test.db', 'n')
        start = time.time()
        for i in xrange(npoints):
            db[str(i)] = coord
        end1 = time.time()
        db.close()
        end2 = time.time()
        return end1-start, end2-start
        

    for i in (10, 100, 1000, 10000, 100000):
        cost1, cost2 = run(i)
        datsize = os.path.getsize('/tmp/test_rtree.dat')
        idxsize = os.path.getsize('/tmp/test_rtree.idx')
        print '%d\t%f\t%f\t%f\t%f\t%f' % (i, cost2 / i, cost1, cost2, datsize/1048576.0, idxsize / 1048576.0)
    exit()

    idx.set((0,0,1,1), '0')
    idx.close()

    
    idx.open(False)
    print map(str, idx.get_pt((0,0)))
    idx.close()
    idx.open(False)
    print map(str, idx.get_pt((0,0)))
    idx.close()

    idx.open(True)
    print map(str, idx.get_pt((0,0)))
    idx.close()
    exit()

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
            
            
        

    
    
