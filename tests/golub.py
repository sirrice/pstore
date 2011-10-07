#
# This analysis assumes that all gene expression arrays are of the same
# genes in the exact same order!!!
#
# if a gene expression row is empty, the values are 0
#

import sys, os, math, numpy
from op import *
from common import *

def load_gene_data(fname):
    names = []
    training = []
    test = []

    f = file(fname, 'r')
    for l in f:
        l = l.strip()
        if len(l) == 0:
            break
        arr = l.split()
        name = arr[0]
        row = map(int, arr[1:])
        names.append(name)
        training.append(row)
    for l in f:
        test.append(map(int, l.split()[1:]))
    return names, numpy.array(training, dtype=float), numpy.array(test, dtype=float)



class Clean(OneToOne):

    def __init__(self, minmax=(100.0, 16000.0)):
        f = lambda x: min(max(x, minmax[0]), minmax[1])
        super(Clean, self).__init__(f)

class Filter(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        data = inputs[0]
        f = lambda row: max(row)/min(row) <= 5 or max(row)-min(row) <= 500
        output = [f(row) and True or False for row in data]

        start = time.time()
        rowsize = len(data[0])
        if pstore.takes_pointers():
            for rowidx, b in enumerate(output):
                fd = pstore.popen((rowidx,))
                if fd:
                    for colidx in xrange(rowsize):
                        pstore.pwrite(fd, 0, (rowidx, colidx))
                pstore.pclose(fd)
        end = time.time()

        return numpy.array(output), {'provoverhead' : end-start}

    def can_box(self):
        return True

    def fmap_0(self, coord, run_id):
        return [(coord[0],)]

    def bmap_0(self, coord, run_id):
        data = self.wrapper.get_inputs(run_id)[0]
        row, ncols = coord[0], len(data[0])
        return [(row, col) for col in xrange(ncols)]

    def fgrp_0(self, coord, run_id):
        return self.bmap_0((coord[0],), run_id)

class NormalizeTrain(OneToOne):
    def __init__(self):
        super(NormalizeTrain, self).__init__(lambda x: math.log(x,10))


class NormalizeTest(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        data, stats = tuple(inputs)
        rowsize = len(data[0])
        output = []

        overhead = 0.0
        for row, stat in zip(data, stats):
            if sum(row) > 0:
                mean, std, b = tuple(stat)
                row = map(lambda v: (math.log(v,10) - mean / std)-b, row)
            output.append(row)


            if pstore.takes_pointers():
                start = time.time()                
                for col in xrange(rowsize):
                    fd = pstore.popen((row,col))
                    if fd:
                        [pstore.pwrite(fd, 0, (row, provcol)) for provcol in xrange(rowsize)]
                        [pstore.pwrite(fd, 1, (0, provcol)) for provcol in xrange(3)]
                    pstore.pclose(fd)
                overhead += time.time() - start

        return numpy.array(output), {'provoverhead' : overhead}


    def can_box(self):
        return True


    def fmap_0(self, coord, run_id):
        data = self.wrapper.get_inputs(run_id)[0]
        row, ncols = coord[0], len(data[0])
        return [(row, col) for col in xrange(ncols)]

    def fmap_1(self, coord, run_id):
        return self.fmap_0(coord, run_id)

    def bmap_0(self, coord, run_id):
        return self.fmap_0(coord, run_id)

    def bmap_1(self, coord, run_id):
        return [(coord[0], col) for col in xrange(3)]

    def fgrp_0(self, coord, run_id):
        return self.fmap_0(coord, run_id)

    def fgrp_1(self, coord, run_id):
        return self.bmap_1(coord, run_id)



class Correlate(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        data, npatients = inputs[0], inputs[1][0]
        corrs = []
        for row in data:
            allrow = row[:npatients]
            amlrow = row[npatients:]

            mean, std = numpy.mean(row), numpy.std(row)
            amlrow = [(val - mean)/std for val in amlrow]
            allrow = [(val - mean)/std for val in allrow]

            amlmean, allmean = numpy.mean(amlrow), numpy.mean(allrow)
            amlstd,  allstd  = numpy.std(amlrow),  numpy.std(allrow)
            corr = (allmean - amlmean) / (amlstd + allstd)
            corrs.append( corr )



        start = time.time()
        rowsize = len(data[0])
        if pstore.takes_pointers():
            for rowidx, b in enumerate(corrs):
                fd = pstore.popen((rowidx,))
                if fd:
                    for colidx in xrange(rowsize):
                        pstore.pwrite(fd, 0, (rowidx, colidx))
                        pstore.pwrite(fd, 1, (0,))
                pstore.pclose(fd)
        end = time.time()
            

        return numpy.array(corrs), {'provoverhead' : end-start}

    def can_box(self):
        return True
        
    def fmap_0(self, coord, run_id):
        return [(coord[0],)]

    def bmap_0(self, coord, run_id):
        data = self.wrapper.get_inputs(run_id)[0]
        row, ncols = coord[0], len(data[0])
        return [(row, col) for col in xrange(ncols)]

    def fgrp_0(self, coord, run_id):
        return self.bmap_0((coord[0],), run_id)


class TTest(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        data, npatients = inputs[0], inputs[1][0]
        corrs = []
        for row in data:
            allrow = row[:npatients]
            amlrow = row[npatients:]
            nall, naml = len(allrow), len(amlrow)

            allmean, amlmean = numpy.mean(allrow), numpy.mean(amlrow)
            allsq, amlsq = sum(map(lambda x: x*x, allrow)), sum(map(lambda x: x*x, amlrow))
            
            sigjoint = (allsq + amlsq - nall*allmean*allmean - naml*amlmean*amlmean) / \
                       (len(allrow) + len(amlrow) - 2)
            top = (allmean - amlmean) / sigjoint
            divfactor = math.pow((1.0/len(allrow) + 1.0/len(amlrow)), 0.5)
            T = top / divfactor

            corrs.append(T)


        start = time.time()
        rowsize = len(data[0])
        if pstore.takes_pointers():
            for rowidx, b in enumerate(corrs):
                fd = pstore.popen((rowidx,))
                if fd:
                    for colidx in xrange(rowsize):
                        pstore.pwrite(fd, 0, (rowidx, colidx))
                        pstore.pwrite(fd, 1, (0,))
                pstore.pclose(fd)
        end = time.time()

        return numpy.array(corrs), {'provoverhead' : end-start}

        
    def fmap_0(self, coord, run_id):
        return [(coord[0],)]

    def bmap_0(self, coord, run_id):
        data = self.wrapper.get_inputs(run_id)[0]
        row, ncols = coord[0], len(data[0])
        return [(row, col) for col in xrange(ncols)]

    def fgrp_0(self, coord, run_id):
        return self.bmap_0((coord[0],), run_id)




class GenStats(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        data, npatients = inputs[0], inputs[1][0]
        stats = []
        for row in data:
            mean = numpy.mean(row)
            std = numpy.std(row)

            if std == 0:
                b = 0
            else:
                all_norm = map(lambda val: (val - mean) / std, row[:npatients])
                aml_norm = map(lambda val: (val - mean) / std, row[npatients:])
                b = (numpy.mean(aml_norm) + numpy.mean(all_norm)) / 2
            stats.append((mean, std, b))

        start = time.time()
        rowsize = len(data[0])
        if pstore.takes_pointers():
            for rowidx, b in enumerate(stats):
                for colidx in xrange(3):
                    fd = pstore.popen((rowidx,colidx))
                    if fd:
                        for provcol in xrange(rowsize):
                            pstore.pwrite(fd, 0, (rowidx, provcol))
                            pstore.pwrite(fd, 1, (0,))
                    pstore.pclose(fd)
        end = time.time()
            

        return numpy.array(stats), {'provoverhead' : end-start}

    def can_box(self):
        return True
        
    def fmap_0(self, coord, run_id):
        return [(coord[0],col) for col in xrange(3)]

    def bmap_0(self, coord, run_id):
        data = self.wrapper.get_inputs(run_id)[0]
        row, ncols = coord[0], len(data[0])
        return [(row, col) for col in xrange(ncols)]

    def fgrp_0(self, coord, run_id):
        return self.bmap_0((coord[0],), run_id)



class PredictorMask(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        filtermask, corr = tuple(inputs)
        print sum([1 for b in filtermask if not b])

        corr = enumerate(corr)
        corr = filter(lambda p: not filtermask[p[0]], corr)
        ordered = sorted(corr, key=lambda p: p[1], reverse=True)
        print min(corr), max(corr)
        predictors = []
        predictors.extend(ordered[:25])
        predictors.extend(ordered[-25:])

        overhead = 0.0
        output = [False] * len(filtermask)
        for idx, v in predictors:
            output[idx] = True

            start = time.time()
            if pstore.takes_pointers():
                fd = pstore.popen((idx,))
                if fd:
                    pstore.pwrite(fd, 0, (idx,))
                    pstore.pwrite(fd, 1, (idx,))
                pstore.pclose(fd)
            overhead += time.time() - start
            
        return numpy.array(output), {'provoverhead' : overhead}

class Predict(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        data, corrs, mask = tuple(inputs)
        npatients = len(data[0])

        output = []
        for pid in xrange(npatients):
            patient = self.get_patient(data, pid)
            vs = [x*w for x, w, use in zip(patient, corrs, mask) if use and x is not None]
            ALL = abs(sum(filter(lambda v: v > 0, vs)))
            AML = abs(sum(filter(lambda v: v < 0, vs)))
            conf = abs(ALL-AML) / (AML + ALL)

            output.append((ALL > AML, conf))
            
            start = time.time()
            if pstore.takes_pointers():
                for outidx in xrange(2):
                    fd = pstore.popen((pid, outidx))
                    if fd:
                        for idx in xrange(len(data)):
                            pstore.pwrite(fd, 0, (idx,))
                            pstore.pwrite(fd, 1, (idx,))
                            pstore.pwrite(fd, 2, (idx,))
                    pstore.pclose(fd)
            end = time.time()

        return numpy.array(output), {'provoverhead' : end - start}


    def get_patient(self, test, col):
        return numpy.array([row[col] for row in test])

    def fmap_0(self, coord, run_id):
        return [coord]
    def fmap_1(self, coord, run_id):
        return self.fmap_0(coord, run_id)
    def fmap_2(self, coord, run_id):
        return self.fmap_0(coord, run_id)
    def bmap_0(self, coord, run_id):
        return self.fmap_0(coord, run_id)        
    def bmap_1(self, coord, run_id):
        return self.fmap_0(coord, run_id)        
    def bmap_2(self, coord, run_id):
        return self.fmap_0(coord, run_id)        
