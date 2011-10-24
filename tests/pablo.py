import sys
sys.path.append("../")
from op import *
from common import *

# get names
# subsample just data
# transpose
# extract nonsubclass
# extract subclass
# learn model
# predict



class Table(object):
    def __init__(self, dims, vals):
        self.dims = dims
        self.vals = vals

    @staticmethod
    def new(*arrs):
        counts = {}
        dims = [set() for arr in arrs]
        for arr, dim in zip(arrs, dims):
            dim.update(filter(lambda x: x != 'NA', arr))
        dims = map(list, dims)
        matrix = np.zeros(map(len, dims))

        for vals in zip(*arrs):
            if 'NA' in vals: continue
            idxs = tuple([dim.index(v) for v, dim in zip(vals, dims)])
            matrix[idxs] += 1
        return Table(dims, matrix)

    def get(self, vals):
        idxs = tuple([dim.index(v) for v, dim in zip(vals, self.dims)])
        return self.vals[idxs]

    def dim(self, dim):
        if dim == 0:
            a = [ np.sum(self.vals[idx,:]) for idx in xrange(len(self.dims[dim]))]
        else:
            a = [ np.sum(self.vals[:,idx]) for idx in xrange(len(self.dims[dim]))]
        a, d = np.array(a), [self.dims[dim]]
        return Table(d, a)            


def get_ev(t):
    t.vals += 2
    t.vals /= np.sum(t.vals)
    marg = t.dim(1)
    pk = t.dim(0)
    for i in xrange(len(t.dims[0])):
        t.vals[i,:] /= marg.vals


    ev = {}
    for j in xrange(len(t.dims[1])):
        top = t.get(['Y', t.dims[1][j]]) / t.get(['N', t.dims[1][j]])
        bot = (pk.get(['Y']) / pk.get(['N']))
        ev[t.dims[1][j]] = math.log(top / bot, math.e)
    return ev








class GetNames(Op):
    def run(self, inputs, run_id):
        #self.wrapper.set_inputs(inputs)
        pstore = self.pstore(run_id)

        arr = inputs[0]
        output = arr[:,0]

        if pstore.uses_mode(Mode.PTR):
            for x in xrange(arr.shape[0]):
                pstore.write(((0, x),), ((x, 0), ))
        return output.reshape((1, arr.shape[0])), {}
            
    def output_shape(self, run_id):
        return (1, self.wrapper.get_input_shape(run_id, 0)[0])
        
    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]

    def fmap(self, coord, run_id, arridx):
        if coord[1] != 0:
            return ()
        return ((0, coord[0]), )
    def bmap(self, coord, run_id, arridx):
        if coord[0] != 0:
            return ()
        return ( (coord[1], 0), )

    
class ExtractNonsubclass(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        
        m2 = inputs[0]
        names = inputs[1]
        features = inputs[2]

        names = names[0]
        features = features[0]

        namesl = names.tolist()
        locs = [namesl.index(feat) for feat in features if feat in namesl]
        output = m2[:,locs]

        if pstore.uses_mode(Mode.PTR):
            for idx, loc in enumerate(locs):
                outcoords = ( (x, idx) for x in xrange(m2.shape[0]) )
                incoords1 = ((0, loc), )
                incoords2 = ((0, idx), )
                for outcoord in outcoords:
                    pstore.write( (outcoord,), ((outcoord[0], loc),), incoords1, incoords2 )
        if pstore.uses_mode(Mode.PT_MAPFUNC):
            for idx, loc in enumerate(locs):
                outcoords = ( (x, idx) for x in xrange(m2.shape[0]) )
                s = struct.pack( 'II', idx, loc )
                for outcoord in outcoords:
                    pstore.write( (outcoord,), s )
        return output, {}

    def output_shape(self, run_id):
        return (self.wrapper.get_input_shape(run_id, 0)[0],
                self.wrapper.get_input_shape(run_id, 2)[1])

    def supported_modes(self):
        return [Mode.PTR, Mode.PT_MAPFUNC]

    def fmap_obj(self, obj, run_id, arridx):
        return obj[0]

    def bmap_obj(self, obj, run_id, arridx):
        shape = self.wrapper.get_input_shape(run_id, 0)
        idx, loc = struct.unpack('II', obj[1])
        
        if arridx == 0:
            return [ (coord[0], loc) for coord in obj[0] ]
        if arridx == 1:
            return ((0, loc), )
        if arridx == 2:
            return ((1, idx), )


class CreateModel(Op):
        
    def get_ev(self, t):
        t.vals += 2
        t.vals /= np.sum(t.vals)
        marg = t.dim(1)
        pk = t.dim(0)
        for i in xrange(len(t.dims[0])):
            t.vals[i,:] /= marg.vals

        ev = {}
        for j in xrange(len(t.dims[1])):
            top = t.get(['Y', t.dims[1][j]]) / t.get(['N', t.dims[1][j]])
            bot = (pk.get(['Y']) / pk.get(['N']))
            ev[t.dims[1][j]] = math.log(top / bot, math.e)
        return ev

        
    
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        md = inputs[0]
        m2 = inputs[1]
        names = inputs[2]
        ftypes = inputs[3]

        names = names[0].tolist()
        ftypes = ftypes[0].tolist()
        subclasses = list(set(m2[:, names.index('C6_CLUSTER')].tolist())) 

        ev_list = []
        for k in xrange(len(ftypes)):

            if ftypes[k] == 'target.feature':
                t = Table.new(md[:,0])
                t.vals += 2
                t.vals /= np.sum(t.vals)
                post_OR = t.get(['Y']) / t.get(['N'])
                logval = math.log(post_OR, math.e)

                ev_list.append( (k, '', '', logval) )

                if pstore.uses_mode(Mode.PTR):
                    outcoords = ((0, len(ev_list)-1), )
                    incoord0 = [ (x, 0) for x in xrange(md.shape[0]) ]
                    incoord3 = ((0, k), )
                    pstore.write(outcoords, incoord0, (), (), incoord3)
                if pstore.uses_mode(Mode.PT_MAPFUNC):
                    outcoords = ((0, len(ev_list)-1), )
                    pstore.write(outcoords, '1_%d' % k)

            elif ftypes[k] in ('clinical', 'pathway', 'subclass', 'genom'):
                ev = self.get_ev(Table.new(md[:, 0], md[:, k]))
                ev = [ (k, key, '', val) for key, val in ev.items()]
                ev_list.extend(ev)

                if pstore.uses_mode(Mode.PTR):
                    evidxs = range(len(ev_list) - len(ev), len(ev_list))
                    outcoords = [ (0, evidx) for evidx in evidxs ]
                    incoord0 = [ (x,y) for x in xrange(md.shape[0]) for y in (0, k) ]
                    incoord3 = ((0,k),)
                    pstore.write(outcoords, incoord0, (), (), incoord3)
                if pstore.uses_mode(Mode.PT_MAPFUNC):
                    evidxs = range(len(ev_list) - len(ev), len(ev_list))
                    outcoords = [ (0, evidx) for evidx in evidxs ]
                    pstore.write(outcoords, '2_%d' % k)

                
            else:
                c6idx = names.index('C6_CLUSTER')
                cpa = Table.new(md[:,0], md[:,k], m2[:,c6idx])

                eva = {}
                for s in xrange(len(subclasses)):
                    sub = subclasses[s]
                    t = Table(cpa.dims[:2], cpa.vals[:,:,s])
                    ev = self.get_ev(t)
                    ev = [ (k, sub, key, val) for key, val in ev.items() ]
                    ev_list.extend(ev)

                    if pstore.uses_mode(Mode.PTR):
                        evidxs = range(len(ev_list) - len(ev), len(ev_list))
                        outcoords = [ (0, evidx) for evidx in evidxs ]
                        incoord0 = [ (x,y) for x in xrange(md.shape[0]) for y in (0, k) ]
                        incoord1 = [ (x,c6idx) for x in xrange(m2.shape[0]) ]
                        incoord2 = ( (0, c6idx), )
                        incoord3 = ( (0,k), )
                        pstore.write( outcoords, incoord0, incoord1, incoord2, incoord3)

                    if pstore.uses_mode(Mode.PT_MAPFUNC):
                        evidxs = range(len(ev_list) - len(ev), len(ev_list))
                        outcoords = [ (0, evidx) for evidx in evidxs ]
                        pstore.write(outcoords, '3_%d_%d' % (k, c6idx))

        dtype = ','.join(map(lambda i: 'a%d' % i, [max(len(x[1]) for x in ev_list),
                                                   max(len(x[2]) for x in ev_list)]))
        output = np.array([ev_list], dtype='i4,%s,f4' %
                                                   dtype).reshape((1,
                                                   len(ev_list)))

        return output, {}



    def output_shape(self, run_id):
        md = self.wrapper.get_input(run_id, 0)
        ftypes = self.wrapper.get_input(run_id, 3)[0].tolist()

        l = 0
        for k in xrange(len(ftypes)):
            if ftypes[k] == 'target.feature':
                n = 1
            elif ftypes[k] in ('clinical', 'pathway', 'subclass', 'genom'):
                n = len( filter(lambda x: x != 'NA', set(md[:,k].tolist())) )
            else:
                n = len( filter(lambda x: x != 'NA', set(md[:,k].tolist())) ) * 6
            l += n
        return (1, l)

    def alltoall(self, arridx):
        return arridx == (0, 3)

    def supported_modes(self):
        return [Mode.PTR, Mode.PT_MAPFUNC]

    def fmap_obj(self, obj, run_id, arridx):
        return obj[0]

    def bmap_obj(self, obj, run_id, arridx):
        encs = map(int, obj[1].split('_'))
        mdshape = self.wrapper.get_input_shape(run_id, 0)
        m2shape = self.wrapper.get_input_shape(run_id, 1)

        if encs[0] == 1:
            if arridx == 0:
                return [ (x, 0) for x in xrange(mdshape[0])]
            if arridx == 3:
                return ((0, encs[1]), )
        if encs[0] == 2:
            if arridx == 0:
                return [ (x, y) for x in xrange(mdshape[0]) for y in (0, encs[1]) ]
            if arridx == 3:
                return ((0, encs[1]), )
        if encs[0] == 3:
            if arridx == 0:
                return [ (x, y) for x in xrange(mdshape[0]) for y in (0, encs[1]) ]
            if arridx == 1:
                return [ (x, encs[2]) for x in xrange(m2shape[0]) ]
            if arridx == 2:
                return ((0, encs[2]), )
            if arridx == 3:
                return ((0, encs[1]), )
        return ()
            
        

class Predict(Op):
    def create_index(self, ev_list):
        d = {}
        for idx, rec in enumerate(ev_list):
            k1, k2, k3, v = tuple(rec)
            d[(k1,k2,k3)] = (idx, v)
        return d



    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        ev_list = inputs[0]
        md = inputs[1]
        m2 = inputs[2]
        ftypes = inputs[3]
        names = inputs[4]

        ev_list = ev_list[0]
        ftypes = ftypes[0]
        names = names[0].tolist()
        ev_idx = self.create_index(ev_list)
        logodds = np.zeros((md.shape[0], len(ftypes)))

        for i in xrange(md.shape[0]):
            for k in xrange(len(ftypes)):

                if md[i,k] != 'NA':
                    
                    if ftypes[k] == 'target.feature':
                        idx, log_or = ev_idx[(k, '', '')]

                        if pstore.uses_mode(Mode.PTR):
                            outcoords = ((i, k), )
                            incoords0 = ((0, idx), )
                            incoords3 = ((0, k), )
                            pstore.write(outcoords, incoords0, (), (), incoords3, ())
                        if pstore.uses_mode(Mode.PT_MAPFUNC):
                            outcoords = ((i, k), )
                            try:
                                pstore.write(outcoords, '1_%d' % idx)
                            except:
                                print pstore
                                print pstore.strat
                                raise
                                    
                        
                    elif ftypes[k] in ('clinical', 'subclass', 'genom', 'pathway'):
                        idx, log_or = ev_idx[(k, md[i,k], '')]

                        if pstore.uses_mode(Mode.PTR):
                            outcoords = ((i, k), )
                            incoords0 = ((0, idx), )
                            incoords1 = ((i, k), )
                            incoords3 = ((0, k), )
                            pstore.write(outcoords, incoords0, incoords1, (), incoords3, ())
                        if pstore.uses_mode(Mode.PT_MAPFUNC):
                            outcoords = ((i, k), )
                            pstore.write(outcoords, '2_%d' % idx)
                        
                    else:
                        sc =  m2[i, names.index('C6_CLUSTER')]
                        idx, log_or = ev_idx[(k, sc, md[i,k])]

                        if pstore.uses_mode(Mode.PTR):
                            outcoords = ((i, k), )
                            incoords0 = ((0, idx), )
                            incoords1 = ((i, k), )
                            incoords2 = ((i, names.index('C6_CLUSTER')), )
                            incoords3 = ((0, k), )
                            incoords4 = ((0, names.index('C6_CLUSTER')), )
                            pstore.write(outcoords, incoords0, incoords1, incoords2, incoords3, incoords4)
                        if pstore.uses_mode(Mode.PT_MAPFUNC):
                            outcoords = ((i, k), )
                            pstore.write(outcoords, '3_%d_%d' % (idx, names.index('C6_CLUSTER')))
                        
                    logodds[i,k] = log_or

        return logodds, {}
                        
    def alltoall(self, arridx):
        return arridx in (0, 1, 3)

    def output_shape(self, run_id):
        x = self.wrapper.get_input_shape(run_id, 1)[0]
        y = self.wrapper.get_input_shape(run_id, 3)[1]
        return (x,y)

    def supported_modes(self):
        return [Mode.PT_MAPFUNC, Mode.PTR]

    def bmap_obj(self, obj, run_id, arridx):
        enc = map(int, obj[1].split('_'))
        
        if arridx == 0:
            return ((0, enc[1]), )
        if arridx == 1:
            if enc[0] > 0:
                return obj[0]
        if arridx == 2:
            if enc[0] == 3:
                return ((obj[0][0][0], enc[2]), )
        if arridx == 3:
            return ((0, obj[0][0][0]), )
        if arridx == 4:
            if enc[0] == 3:
                return ((0, enc[2]), )
        return ()
            
    def fmap_obj(self, obj, run_id, arridx):
        return obj[0]
        

class CumOdds(Op):

    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        logodds = inputs[0]
        cumodds = np.array([np.sum(logodds[idx, :]) for idx in xrange(logodds.shape[0])])

        if pstore.uses_mode(Mode.PTR):
            for idx in xrange(logodds.shape[0]):
                outcoord = (0, idx)
                incoords = [ (idx, y) for y in xrange(logodds.shape[1]) ]
                pstore.write((outcoord, ), incoords)
        
        return cumodds.reshape(1, len(cumodds)), {}

    def output_shape(self, run_id):
        return (1, self.wrapper.get_input_shape(run_id, 0)[0])

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]

    def fmap(self, coord, run_id, arridx):
        return ((0, coord[0]), )
    
    def bmap(self, coord, run_id, arridx):
        return [ (coord[1], y) for y in xrange(self.wrapper.get_input_shape(run_id, 0)[1])]




class GetClasses(OneToOne) :
    def __init__(self):
        f = lambda prob: prob > 0.5 and 'Y' or 'N' 
        super(GetClasses, self).__init__(f)

    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        inputs = inputs[0] # 2d array
        output = np.empty(inputs.shape, str)
        output = self.f(inputs, out=output).astype(str)
        
        start = time.time()
        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([1])
            pstore.set_inareas([1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))
        elif pstore.uses_mode(Mode.PTR):
            nrow, ncol = inputs.shape
            for ridx in xrange(nrow):
                for cidx in xrange(ncol):
                    pstore.write(((ridx, cidx),), ((ridx, cidx),))
        end = time.time()            

        return output, {'provoverhead' : end-start}



class Validate(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        klasses = inputs[0]
        m2 = inputs[1]
        
        klasses = klasses[0].tolist()
        actual = m2[:, 0].tolist()

        d = {}
        prov = {}
        for idx, key in enumerate(zip(klasses, actual)):
            if key not in d:
                d[key] = 0
                if pstore.uses_mode(Mode.PTR):
                    prov[key] = []
            d[key] += 1
            if pstore.uses_mode(Mode.PTR):
                prov[key].append(idx)



        if pstore.uses_mode(Mode.PTR):
            for idx, key in enumerate(d.keys()):
                outcoord = (0, idx)
                incoords0 = [ (0, inidx) for inidx in prov[key] ]
                incoords1 = [ (inidx, 0) for inidx in prov[key] ]
                pstore.write( (outcoord,), incoords0, incoords1 )
        return np.array([d.items()]), {}

    def output_shape(self, run_id):
        return (1, 4)

    def supported_modes(self):
        return [Mode.PTR]

        
