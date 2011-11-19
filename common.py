import sys, os, math, numpy
from op import *
from operator import mul
import numpy as np
from strat import *


class DummyOp(Op):
    """
    extends operator with dummy wrapper for debugging purposes
    so the operator can be run independent of runtime and workflow
    """
    class Wrapper(object):
        def __init__(self):
            self.inputs = []
        def set_inputs(self, inputs):
            self.inputs = inputs
        def get_inputs(self):
            return self.inputs
        def get_input_shapes(self, run_id):
            return map(lambda i: i.shape, self.inputs)
        def get_input_shape(self, run_id, arridx):
            return self.inputs[arridx].shape
        def get_output_shape(self, run_id):
            return (100, 100)

    def __init__(self):
        super(DummyOp, self).__init__()
        self.wrapper = DummyOp.Wrapper()

    def pstore(self, run_id):
        return IPstore(self, run_id, None, Desc(Mode.NOOP, Spec.default(), True))
    
        



class OneToOne(Op):

    def __init__(self, f=lambda x:x):
        super(OneToOne, self).__init__()
        self.f = np.frompyfunc(f, 1, 1)

    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        inputs = inputs[0] # 2d array
        output = np.empty(inputs.shape, float)
        output = self.f(inputs, out=output).astype(float)
        #output = numpy.array(map(self.f, inputs.flat), dtype=float).reshape(inputs.shape)
        
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

    def alltoall(self, arridx):
        return True

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id, 0)

    def can_box(self):
        return True

    def fmap(self, coord, run_id, arridx):
        return (tuple(coord),)

    def bmap(self, coord, run_id, arridx):
        return (tuple(coord),)

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]


class Mean(Op):
    "returns array of same size, where every cell is the mean of the source array"
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        inputs = inputs[0]  # 2d array

        total = np.sum(inputs)
        ncells = reduce(mul, inputs.shape)

        start = time.time()
        x = 0
        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([reduce(mul,inputs.shape)])
            pstore.set_inareas([reduce(mul,inputs.shape)])
            pstore.set_outarea(reduce(mul,inputs.shape))
            pstore.set_ncalls(reduce(mul, inputs.shape))
            pstore.set_noutcells(reduce(mul, inputs.shape))
        elif pstore.uses_mode(Mode.PTR):
            nrow, ncol = inputs.shape
            coords = [(ridx, cidx) for ridx in xrange(nrow) for cidx in xrange(ncol)]
            pstore.write(coords, coords)
        end = time.time()

        mean = float(total / ncells)
        output = np.ones_like(inputs) * mean
        return output, {'provoverhead' : end-start}

    def alltoall(self, arridx):
        return True

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id,0)
    

    def fmap(self, coord, run_id, arridx):
        shape = self.wrapper.get_input_shape(run_id, 0)
        nrows, ncols = shape
        for i in xrange(nrows):
            for j in xrange(ncols):
                yield (i,j)

    def bmap(self, coord, run_id, arridx):
        return self.fmap(coord, run_id, arridx)

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]
    

class MeanSingleVal(Op):
    "returns array of same size, where every cell is the mean of the source array"
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        inputs = inputs[0]  # 2d array

        total = numpy.sum(inputs)
        ncells = reduce(mul, inputs.shape)
        
        start = time.time()
        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([reduce(mul,inputs.shape)])
            pstore.set_inareas([reduce(mul,inputs.shape)])
            pstore.set_outarea(1)
            pstore.set_ncalls(1)
            pstore.set_noutcells(1)
        elif pstore.uses_mode(Mode.PTR):
            nrow, ncol = inputs.shape
            pstore.write(((0,0),),  [(x,y) for x in xrange(nrow) for y in xrange(ncol)])

        end = time.time()

        mean = total / ncells
        output = np.array([[mean]])
        return output, {'provoverhead' : end-start}

    def alltoall(self, arridx):
        return True

    def output_shape(self, run_id):
        return (1,1)
    

    def fmap(self, coord, run_id, arridx):
        return ((0,0),)

    def bmap(self, coord, run_id, arridx):
        shape = self.wrapper.get_input_shape(run_id, 0)
        return ((x,y) for x in xrange(shape[0]) for y in xrange(shape[1]))

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]


class Merge(Op):
    
    def __init__(self, f=lambda a,b:a+b):
        super(Merge, self).__init__()
        self.f = numpy.frompyfunc(f,2,1)
    
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        left = inputs[0]
        right = inputs[1]

        if left.shape != right.shape:
            raise RuntimeError, "input shapes not the same\t%s\t%s\t%s" % (self, left.shape, right.shape)
        
        nrow, ncol = left.shape

        output = np.empty(left.shape, float)
        output = self.f(left, right, out=output).astype(float)

        start = time.time()
        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([1,1])
            pstore.set_inareas([1,1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))
        elif pstore.uses_mode(Mode.PTR):
            for x in xrange(nrow):
                for y in xrange(ncol):
                    coords = ((x,y),)
                    pstore.write(coords, coords, coords)
        end = time.time()
        return output, {'provoverhead' : (end-start)}

    def alltoall(self, arridx):
        return True

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id,0)

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]

    def fmap(self, coord, run_id, arridx):
        return (tuple(coord),)

    def bmap(self, coord, run_id, arridx):
        return (tuple(coord),)

class Diff(Merge):
    def __init__(self):
        super(Diff, self).__init__(lambda a,b: a-b)



class Subsample(Op):
    def __init__(self, box):
        """
        @param box: [(min x, max x), (min y, max y)]
        """
        super(Subsample, self).__init__()
        self.box = box  #box should be a numpy array

    
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        arr = inputs[0]
        #box = inputs[1] # ((minx, max x), (miny, max y))
        box = self.box

        output = arr[box[0][0]:box[0][1], box[1][0]:box[1][1]]

        start = time.time()

        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([1])
            pstore.set_inareas([1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))
        elif pstore.uses_mode(Mode.PTR):
            mins = (box[0][0], box[1][0])
            outcoords = ((x,y) for x in xrange(output.shape[0]) for y in xrange(output.shape[1]))
            #boxcoords = [(x,y) for x in xrange(2) for y in xrange(2)]
            for outcoord in outcoords:
                #print outcoord[0] + mins[0], outcoord[1]+mins[1]
                pstore.write( (outcoord,), ((outcoord[0]+mins[0], outcoord[1]+mins[1]), ))
        end = time.time()

        return output, {'provoverhead' : (end-start)}

    def takes(self, x, I):
        """
        takes(x, I): Takes a subgrid from array x.
        I is a list of list of subindices.
        """
        for i in xrange(len(I)):
            ii = I[i]
            if ii is not None:
                x = np.take(x, ii, i)
        return x

    def alltoall(self, arridx):
        return False
    
    def output_shape(self, run_id):
        return (self.box[0][1]-self.box[0][0], self.box[1][1]-self.box[1][0])

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]

    def fmap(self, coord, run_id, arridx):
        # if coord in box, return coord
        shape = self.output_shape(run_id)
        if coord[0] >= shape[0] and coord[1] >= shape[1]:
            return ()
        incoord = (coord[0]-self.box[0][0], coord[1]-self.box[1][0])
        if incoord[0] >= 0 and incoord[1] >= 0:
            return (incoord, )
        return ()

    def bmap(self, coord, run_id, arridx):
        ocoord = (coord[0]-self.box[0][0], coord[1]-self.box[1][0])
        return (ocoord, )




class Transpose(Op):
    
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        arr = inputs[0]

        output = arr.transpose().copy()

        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([1])
            pstore.set_inareas([1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))
        elif pstore.uses_mode(Mode.PTR):
            for x in xrange(output.shape[0]):
                for y in xrange(output.shape[1]):
                    pstore.write(((x,y), ), ((y,x), ))
        return output, {}

    def alltoall(self, arridx):
        return True
    
    def output_shape(self, run_id):
        shape = self.wrapper.get_input_shape(run_id, 0)
        return (shape[1], shape[0])

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]

    def fmap(self, coord, run_id, arridx):
        return ((coord[1], coord[0]), )

    def bmap(self, coord, run_id, arridx):
        return ((coord[1], coord[0]), )
