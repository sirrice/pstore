import sys, os, math, numpy
from op import *
from operator import mul
import numpy as np
from strat import *


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
        if pstore.strat.mode == Mode.PTR:
            nrow, ncol = inputs.shape
            for ridx in xrange(nrow):
                for cidx in xrange(ncol):
                    pstore.write(((ridx, cidx),), ((ridx, cidx),))
        else:
            pstore.set_fanins([1])
            pstore.set_inareas([1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))
        end = time.time()            

        return output, {'provoverhead' : end-start}

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id, 0)

    def can_box(self):
        return True

    def fmap(self, coord, run_id, arridx):
        return (tuple(coord),)

    def bmap(self, coord, run_id, arridx):
        return (tuple(coord),)

    def supported_model(self):
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
        if pstore.strat.mode == Mode.PTR:
            nrow, ncol = inputs.shape
            coords = [(ridx, cidx) for ridx in xrange(nrow) for cidx in xrange(ncol)]
            pstore.write(coords, coords)
        else:
            pstore.set_fanins([reduce(mul,inputs.shape)])
            pstore.set_inareas([reduce(mul,inputs.shape)])
            pstore.set_outarea(reduce(mul,inputs.shape))
            pstore.set_ncalls(reduce(mul, inputs.shape))
            pstore.set_noutcells(reduce(mul, inputs.shape))
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

    def supported_model(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]
    

class MeanSingleVal(Op):
    "returns array of same size, where every cell is the mean of the source array"
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        inputs = inputs[0]  # 2d array

        total = numpy.sum(inputs)
        ncells = reduce(mul, inputs.shape)
        
        start = time.time()
        if pstore.strat.mode == Mode.PTR:
            nrow, ncol = inputs.shape
            pstore.write(((0,0),),  [(x,y) for x in xrange(nrow) for y in xrange(ncol)])
        else:
            pstore.set_fanins([reduce(mul,inputs.shape)])
            pstore.set_inareas([reduce(mul,inputs.shape)])
            pstore.set_outarea(1)
            pstore.set_ncalls(1)
            pstore.set_noutcells(1)

        end = time.time()

        mean = total / ncells
        output = np.array([[mean]])
        return output, {'provoverhead' : end-start}

    def alltoall(self, arridx):
        return True

    def output_shape(self, run_id):
        return (1,1)
    

    def fmap(self, coord, run_id, arridx):
        return ((0,),)

    def bmap(self, coord, run_id, arridx):
        shape = self.wrapper.get_input_shape(run_id, 0)
        return ((x,y) for x in xrange(shape[0]) for y in xrange(shape[1]))

    def supported_model(self):
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
        if pstore.strat.mode == Mode.PTR:
            for x in xrange(nrow):
                for y in xrange(ncol):
                    coords = ((x,y),)
                    pstore.write(coords, coords, coords)
        else:
            pstore.set_fanins([1,1])
            pstore.set_inareas([1,1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))
        end = time.time()
        return output, {'provoverhead' : (end-start)}

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id,0)

    def supported_model(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]

    def fmap(self, coord, run_id, arridx):
        return (tuple(coord),)

    def bmap(self, coord, run_id, arridx):
        return (tuple(coord),)

