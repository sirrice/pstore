# Simple 2D image processing pipeline that simulates a simplified
# LSST pipeline.

# Intakes FITS files, and produces arrays of stars
import time, math, sys, logging, os, hashlib
import numpy as np
from common import *
from scipy import ndimage, signal
from op import *
from runtime import *
from arraystore import *

log = logging.getLogger('lsst')
logging.basicConfig()
log.setLevel(logging.ERROR)



class Clip(OneToOne):
    def __init__(self, f=lambda x: True):
        def lamb(x):
            if f(x):
                return x
            return 0.0
        super(Clip, self).__init__(lamb)


    


class MergeSingleVal(Op):

    def __init__(self, f=lambda a,b: a+b):
        super(MergeSingleVal, self).__init__()
        def tmp(singleval):
            def h(el):
                return f(el, singleval)
            return numpy.frompyfunc(h,1,1)
        self.gen_f = tmp
        self.f = f

    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)

        left = inputs[0]
        right = inputs[1][0][0] # the single value

        nrow, ncol = left.shape
        f = self.gen_f(right)
        output = np.empty(left.shape, float)
        output = f(left, out=output).astype(float)
        #output = numpy.array(map(lambda x: self.f((x,right)), left.flat)).reshape(left.shape)

        start = time.time()
        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([1,1])
            pstore.set_inareas([1,1])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, output.shape))
            pstore.set_noutcells(reduce(mul, output.shape))

        if pstore.uses_mode(Mode.PTR):
            for x in xrange(nrow):
                for y in xrange(ncol):
                    pstore.write(((x,y),),  ((x,y),), ((0,0),))
        end = time.time()
        return output, {'provoverhead' : (end-start)}

    def alltoall(self, arridx):
        return True
        return arridx == 1
        
    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id,0)


    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PTR]


    def fmap(self, coord, run_id, arridx):
        if arridx == 0:
            return (coord,)
        elif arridx == 1:
            shape = self.wrapper.get_output_shape(run_id)
            return ((x,y) for x in xrange(shape[0]) for y in xrange(shape[1]))

    def bmap(self, coord, run_id, arridx):
        if arridx == 0:
            return (tuple(coord),)
        elif arridx == 1:
            return ((0,0),)

class DiffSingleVal(MergeSingleVal):
    def __init__(self):
        super(DiffSingleVal, self).__init__(lambda a,b: a-b)

class BadPixRemoval(Merge):
    def __init__(self):
        def f(a,b):
            return not b and a or 90.0
        super(BadPixRemoval, self).__init__(f)

class Add(Merge):
    def __init__(self):
        super(Add, self).__init__(lambda a,b: a+b)

class LogicalAnd(Merge):
    def __init__(self):
        def f(a,b):
            return (a and b) and 1 or 0
        super(LogicalAnd, self).__init__(f)#lambda pair: (pair[0] and pair[1] and 1) or 0)



class Div(Merge):
    def __init__(self):
        super(Div, self).__init__(lambda a,b: a/b)
        




class Convolve(Op):
    """
    array, kernel -> convolved array
    """
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        
        arr = inputs[0]
        kernel = inputs[1]

        ar, ac = arr.shape
        kr, kc = kernel.shape[0]/2, kernel.shape[1]/2

        start = time.time()
        if pstore.uses_mode(Mode.FULL_MAPFUNC):
            pstore.set_fanins([1,reduce(mul, kernel.shape)])
            pstore.set_inareas([1,reduce(mul, kernel.shape)])
            pstore.set_outarea(1)
            pstore.set_ncalls(reduce(mul, arr.shape))
            pstore.set_noutcells(reduce(mul, arr.shape))

        if pstore.uses_mode(Mode.PTR):
            for rid in xrange(ar):
                for cid in xrange(ac):
                    minr, maxr = (max(0,rid - kr), min(ar, rid + kr+1))
                    minc, maxc = (max(0,cid - kc), min(ac, cid + kc+1))
                    prov0 = [(px, py) for px in xrange(minr, maxr) for py in xrange(minc, maxc)]
                    prov1 = [(kx, ky) for kx in range(maxr-minr) for ky in xrange(maxc-minc)]
                    pstore.write(((rid, cid),), prov0, prov1)

        if pstore.uses_mode(Mode.PT_MAPFUNC):
            for x in xrange(ar):
                for y in xrange(ac):
                    pstore.write(((x,y),), '')
        end = time.time()

        output = np.empty(arr.shape, float)
        ndimage.convolve(arr, kernel, output=output, mode='constant', cval=0.0)
        return output, {'provoverhead' : end - start}

    def alltoall(self, arridx):
        return True

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC, Mode.PT_MAPFUNC, Mode.PTR]
                
    def fmap(self, coord, run_id, arridx):
        if arridx == 0:
            x,y = tuple(coord)
            arrshape = self.wrapper.get_input_shape(run_id, 0)
            kernshape = self.wrapper.get_input_shape(run_id, 1)

            kr,kc = int(kernshape[0] / 2.0), int(kernshape[1] / 2.0)
            minx, maxx = max(0, x - kr), min(arrshape[0], x + kr+1)
            miny, maxy = max(0, y - kc), min(arrshape[1], y + kc+1)

            return ((x,y) for x in xrange(minx, maxx) for y in xrange(miny, maxy))
        elif arridx == 1:
            shape = self.wrapper.get_input_shape(run_id, 0)
            return ((x,y) for x in xrange(shape[0]) for y in xrange(shape[1]))

    def bmap(self, coord, run_id, arridx):
        if arridx == 0:
            return self.fmap(coord, run_id, arridx)
        elif arridx == 1:
            shape = self.wrapper.get_input_shape(run_id, 1)
            return ((x,y) for x in xrange(shape[0]) for y in xrange(shape[1]))

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id, 0)

    def bmap_obj(self, obj, run_id, arridx):
        return self.bmap(obj[0][0], run_id, arridx)

    def fmap_obj(self, obj, run_id, arridx):
        return self.fmap(obj[0][0], run_id, arridx)



    
class Cluster(Op):

    def run(self, inputs, run_id):
        """
        inputs = [array, [threshold]]
        output = array with cells masked to cluster number value
        """
        pstore = self.pstore(run_id)


        arr = inputs[0]
        threshold = inputs[1][0][0]
        ids = np.zeros(arr.shape)
        nrow, ncol = arr.shape        
        fixups = {}
        nextid = 1
        clusters = {}


        for x in xrange(nrow):
            for y in xrange(ncol):
                if arr[x,y] > threshold:
                    if y > 0 and ids[x,y-1]:
                        oid = ids[x,y-1]
                    elif y > 0 and x > 0 and ids[x-1,y-1]:
                        oid = ids[x-1,y-1]
                    elif x > 0 and ids[x-1,y]:
                        oid = ids[x-1,y]
                    elif x > 0 and y < len(arr[0])-1 and ids[x-1,y+1]:
                        oid = ids[x-1,y+1]
                    else:
                        oid = nextid
                        nextid += 1

                    ids[x,y] = oid
                    
                    if x > 0 and y < len(arr[0])-1:# and oldid and oid != oldid:
                        oldid = ids[x-1,y+1]                        
                        if oldid and oid > oldid:
                            fixups[oid] = oldid
                            self.fix(fixups, oid)

                    if oid not in clusters:
                        clusters[oid] = []
                    clusters[oid].append((x,y))


        #
        # START PROVENANCE STORAGE CODE
        #
        
        start = time.time()
        if pstore.uses_mode(Mode.PTR):
            # stores all pointers explicitly
            for cid, coords in enumerate(clusters.values()):
                pstore.write(coords, coords, ((0,0),))
                
        if pstore.uses_mode(Mode.PT_MAPFUNC):
            # stores just enough information to recalculate provenance
            for cid, coords in enumerate(clusters.values()):
                pstore.write(coords, '')
        end = time.time()

        #
        # END PROVENANCE CODE
        #

        output = ids

        log.info("Cluster.run\toutput: %s", output)
        return np.array(output), {'provoverhead' : (end - start)}

    def supported_modes(self):
        return [Mode.PT_MAPFUNC, Mode.PTR]

    def fix(self, fixups, oid):
        if oid not in fixups: return oid
        if not fixups[oid]: return oid
        
        fixups[oid] = self.fix(fixups, fixups[oid])
        return fixups[oid]

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id, 0)
        
    def can_box(self):
        return True

    def bmap_obj(self, obj, run_id, arridx):
        if arridx == 0:
            return obj[0]
        elif arridx == 1:
            return ((0,0),)

    def fmap_obj(self, obj, run_id, arridx):
        if arridx == 0:
            return obj[0]
        elif arridx == 1:
            return []






class CRDetect(Op):
    def run(self, inputs, run_id):
        pstore = self.pstore(run_id)
        verbose = False

        if verbose:
			print "Convolving image with Laplacian kernel ..."

        
        cleanarray = inputs[0]
        laplkernel = inputs[1]

        growkernel = np.ones((3,3))
        gain = 2.2
        readnoise = 10.0
        sigclip = 5.0
        sigfrac = 0.3
        sigcliplow = sigclip * sigfrac        
        objlim = 5.0
        satlevel = 50000.0
        pssl = 0.0
            
		
        # We subsample, convolve, clip negative values, and rebin to original size
        subsam = self.subsample(cleanarray)
        conved = ndimage.convolve(subsam, laplkernel, mode='reflect')
        #conved = signal.convolve2d(subsam, laplkernel, mode="same", boundary="symm")
        cliped = conved.clip(min=0.0)
        #cliped = np.abs(conved) # unfortunately this does not work to find holes as well ...
        lplus = self.rebin2x2(cliped)
        conved = subsam = cliped = None

		
        if verbose:
			print "Creating noise model ..."
			
		# We build a custom noise map, so to compare the laplacian to
        m5 = ndimage.filters.median_filter(cleanarray, size=5, mode='mirror')
        # We keep this m5, as I will use it later for the interpolation.
        m5clipped = m5.clip(min=0.00001) # As we will take the sqrt
        noise = (1.0/gain) * np.sqrt(gain*m5clipped + readnoise*readnoise)
        m5 = m5clipped = None
 
        if verbose:
            print "Calculating Laplacian signal to noise ratio ..."
 
        # Laplacian signal to noise ratio :
        s = lplus / (2.0 * noise) # the 2.0 is from the 2x2 subsampling
        # This s is called sigmap in the original lacosmic.cl
        
        # We remove the large structures (s prime) :
        sp = s - ndimage.filters.median_filter(s, size=5, mode='mirror')
        
        if verbose:
            print "Selecting candidate cosmic rays ..."
            
        # Candidate cosmic rays (this will include stars + HII regions)
        candidates = sp > sigclip    
        nbcandidates = np.sum(candidates)
        
        if verbose:
            print "  %5i candidate pixels" % nbcandidates
            print "Building fine structure image ..."
            
        # We build the fine structure image :
        m3 = ndimage.filters.median_filter(cleanarray, size=3, mode='mirror')
        m37 = ndimage.filters.median_filter(m3, size=7, mode='mirror')
        f = m3 - m37
        m3 = m37 = None

        f = f / noise
        f = f.clip(min=0.01) # as we will divide by f. like in the iraf version.
        
        if verbose:
            print "Removing suspected compact bright objects ..."
            
        # Now we have our better selection of cosmics :
        cosmics = np.logical_and(candidates, sp/f > objlim)
        # Note the sp/f and not lplus/f ... due to the f = f/noise above.
        
        nbcosmics = np.sum(cosmics)
        
        if verbose:
            print "  %5i remaining candidate pixels" % nbcosmics
        
        # What follows is a special treatment for neighbors, with more relaxed constains.
        
        if verbose:
            print "Finding neighboring pixels affected by cosmic rays ..."
            
        # We grow these cosmics a first time to determine the immediate neighborhod  :
        growcosmics = np.cast['bool'](ndimage.convolve(np.cast['float32'](cosmics), growkernel, mode="mirror"))
        
        # From this grown set, we keep those that have sp > sigmalim
        # so obviously not requiring sp/f > objlim, otherwise it would be pointless
        growcosmics = np.logical_and(sp > sigclip, growcosmics)
        
        # Now we repeat this procedure, but lower the detection limit to sigmalimlow :
        finalsel = np.cast['bool'](ndimage.convolve(np.cast['float32'](growcosmics),
                                                    growkernel, mode="mirror"))
        finalsel = np.logical_and(sp > sigcliplow, finalsel)
        nbfinal = np.sum(finalsel)

        # dealloc everything that's not: cleanarray, finalsel, laplkernel
        lplus = cosmics = growcosmics = f = noise = None
        
        if verbose:
            print "  %5i pixels detected as cosmics" % nbfinal
        
        #
        # START PROVENANCE STORAGE CODE
        # each pixel in CR mask is either
        # a) not a CR, so 1-to-1 provenance or 
        # b) a CR, and dependent on max(kernelsize, 7) pixels surrounding it
        #

        start = time.time()
        if pstore.uses_mode(Mode.PT_MAPFUNC | Mode.FULL_MAPFUNC):
            # more efficient storage.  Storage class 2
            for (x,y) in np.argwhere(finalsel):
                pstore.write(((x,y),), '')

        if pstore.uses_mode(Mode.PTR):
            # store every single pointer explicitly
            # for comparisons, would never do this in real life
            # storage class 3
            prov1 = [(px,py)
                     for px in xrange(laplkernel.shape[0])
                     for py in xrange(laplkernel.shape[1])]
            bound = 3 #max(max(laplkernel.shape), 7) / 2
            for x in xrange(cleanarray.shape[0]):
                for y in xrange(cleanarray.shape[1]):
                    if finalsel[x][y]:
                        mins = (max(0,x-bound),max(0,y-bound))
                        maxs = (min(x+bound+1, cleanarray.shape[0]),
                                min(y+bound+1, cleanarray.shape[1]))
                        prov0 = [(px,py)
                                 for px in xrange(mins[0], maxs[0])
                                 for py in xrange(mins[1], maxs[1])]
                        pstore.write(((x,y),), prov0, prov1)
                    else:
                        coords = ((x,y),)
                        fd = pstore.write(coords, coords, ())

        end = time.time()

        #
        # END PROVENANCE STORAGE CODE
        #

        return finalsel, {'provoverhead' : end - start}


    def subsample(self, a):
        """
        iterative implementation.  only works for 2D images
        """
        newshape = (2*a.shape[0], 2*a.shape[1])
        ret = np.zeros(newshape).astype(int)
        block = 10
        xslice = slice(0,a.shape[0],float(a.shape[0])/newshape[0])
        for i in xrange(int(math.ceil(a.shape[1] / float(block)))):
            offset = i*block
            upper = min(offset+block,a.shape[1])
            slices = [xslice, slice(offset, upper, float(a.shape[1])/newshape[1])]
            indices = np.mgrid[slices].astype(int)
            ret[:,i*block*2:min((i+1)*block*2, newshape[1])] = a[tuple(indices)]
        return ret

    def rebin(self, a, newshape):
        shape = a.shape
        lenShape = len(shape)
        factor = np.asarray(shape)/np.asarray(newshape)
        #print factor
        evList = ['a.reshape('] + \
                 ['newshape[%d],factor[%d],'%(i,i) for i in xrange(lenShape)] + \
                 [')'] + ['.sum(%d)'%(i+1) for i in xrange(lenShape)] + \
                 ['/factor[%d]'%i for i in xrange(lenShape)]
        return eval(''.join(evList))


    def rebin2x2(self, a):
        """
        Wrapper around rebin that actually rebins 2 by 2
        """
        inshape = np.asarray(a.shape)
        #if not (inshape % 2 == np.zeros(2)).all(): # Modulo check to see if size is even
        #    raise RuntimeError, "I want even image shapes !"
        return self.rebin(a, inshape/2)

    def alltoall(self, arridx):
        return arridx == 0

    def can_box(self):
        return True

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id,0)

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC | Mode.PT_MAPFUNC, Mode.PTR]

    def bmap_obj(self, obj, run_id, arridx):
        coord = obj[0][0]
        
        if arridx == 0:
            shape = self.wrapper.get_input_shape(run_id, 0)
            x,y = tuple(coord)
            bound = 3
            mins = ( max(0,x-bound),max(0,y-bound) )
            maxs = ( min(x+bound+1, shape[0]), min(y+bound+1, shape[1]))
            return [(px,py)
                    for px in xrange(mins[0], maxs[0])
                    for py in xrange(mins[1], maxs[1])]
        if arridx == 1:
            shape = self.wrapper.get_input_shape(run_id, 1)
            return [(x,y) for x in xrange(shape[0]) for y in xrange(shape[1])]

    def fmap_obj(self, obj, run_id, arridx):
        if arridx == 0:
            return self.bmap_obj(obj, run_id, arridx)
        raise RuntimeError, "CRDetect: forward function for arridx=1 not implemented"

    def bmap(self, coord, run_id, arridx):
        if arridx == 0:
            return (coord,)
        return ()

    def fmap(self, coord ,run_id, arridx):
        if arridx == 0:
            return (coord,)
        return ()




class RemoveCRs(Op):
    def run(self, inputs, run_id):
        if inputs[0].shape != inputs[1].shape:
            raise RuntimeError
        pstore = self.pstore(run_id)
        verbose = False

        cleanarray = inputs[0]
        mask = inputs[1].astype(bool)

        if verbose:
            print "Cleaning cosmic affected pixels ..."
        
        # So... mask is a 2D array containing False and True, where True means "here is a cosmic"
        # We want to loop through these cosmics one by one.
        cosmicindices = np.argwhere(mask)
        # This is a list of the indices of cosmic affected pixels.

        # We put cosmic ray pixels to np.Inf to flag them :
        cleanarray[mask] = np.Inf
        
        # Now we want to have a 2 pixel frame of Inf padding around our image.
        w,h = tuple(cleanarray.shape)
        padarray = np.zeros((w+4,h+4))+np.Inf
        padarray[2:w+2,2:h+2] = cleanarray.copy() # that copy is important, we need 2 independent arrays
        
        provoverhead = 0.0

        # every cell is at least one to one with itself
        start = time.time()
        #if pstore.strat.mode == Mode.PT_MAPFUNC | Mode.FULL_MAPFUNC:
        if pstore.uses_mode(Mode.PTR):
            for x in xrange(cleanarray.shape[0]):
                for y in xrange(cleanarray.shape[1]):
                    coords = ((x,y),)
                    pstore.write(coords, coords, coords)
        provoverhead += time.time() - start


        
        # A loop through every cosmic pixel :
        for cosmicpos in cosmicindices:
            x = int(cosmicpos[0])
            y = int(cosmicpos[1])
            cutout = padarray[x:x+5, y:y+5].ravel() # remember the shift due to the padding !
            #print cutout
            # Now we have our 25 pixels, some of them are np.Inf, and we want to take the median
            goodcutout = cutout[cutout != np.Inf]
            #print np.alen(goodcutout)
            
            if np.alen(goodcutout) >= 25 :
                # This never happened, but you never know ...
                raise RuntimeError, "Mega error in clean !"
            elif np.alen(goodcutout) > 0 :
                replacementvalue = np.median(goodcutout)
            else :    
                # i.e. no good pixels : Shit, a huge cosmic, we will have to improvise ...
                raise RuntimeError, "OH NO, I HAVE A HUUUUUUUGE COSMIC !!!!!"
                replacementvalue = self.guessbackgroundlevel()
            
            # We update the cleanarray,
            # but measure the medians in the padarray, so to not mix things up...
            cleanarray[x, y] = replacementvalue


            
            start = time.time()
            if pstore.uses_mode(Mode.PTR) or pstore.uses_mode(Mode.FULL_MAPFUNC | Mode.PT_MAPFUNC):
                
                coord = ((x,y),)
                box = ((x,y), (min(x+6,cleanarray.shape[0]), min(y+6,cleanarray.shape[1])))
                if pstore.uses_mode(Mode.PTR):
                    prov0 = [(px,py)
                             for px in xrange(box[0][0],box[1][0])
                             for py in xrange(box[1][0],box[1][1])]
                    pstore.write(coord, prov0, coord)
                else:
                    pstore.write(coord, '')

            provoverhead += time.time() - start


        # That's it.
        if verbose:
            print "Cleaning done"

        return cleanarray, {'provoverhead' : provoverhead}


    def alltoall(self, arridx):
        return True

    def output_shape(self, run_id):
        return self.wrapper.get_input_shape(run_id,0)

    def supported_modes(self):
        return [Mode.FULL_MAPFUNC | Mode.PT_MAPFUNC, Mode.PTR]

    def bmap_obj(self, obj, run_id, arridx):
        coord = obj[0][0]
        if arridx == 0:
            shape = self.wrapper.get_input_shape(run_id, 0)
            x,y = tuple(coord)
            box = ((x,y), (min(x+6,shape[0]), min(y+6,shape[1])))
            print box
            return [(px,py)
                    for px in xrange(box[0][0], box[1][0])
                    for py in xrange(box[0][1],box[1][1])]
        return (coord,)



    def fmap_obj(self, obj, run_id, arridx):
        coord = obj[0][0]
        if arridx == 0:
            x,y = tuple(coord)
            box = ((max(0, x-6),max(0, y-6)), (x,y))
            return ((px,py)
                    for px in xrange(box[0][0], box[1][0])
                    for py in xrange(box[0][1],box[1][1]))
        return (coord,)

    def bmap(self, coord, run_id, arridx):
        shape = self.wrapper.get_input_shape(run_id, arridx)
        if coord[0] >= 0 and coord[0] < shape[0] and coord[1] >= 0 and coord[1] < shape[1]:
            return (coord,)
        return ()

    def fmap(self, coord, run_id, arridx):
        return self.bmap(coord, run_id, arridx)


# class Subsample(Op):
#     def run(self, inputs, run_id):
#         pstore = self.pstore(run_id)

#         a = inputs[0]
#         newshape = (2*a.shape[0], 2*a.shape[1])
#         slices = [slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
#         coordinates = np.mgrid[slices]
#         indices = coordinates.astype('i')   #choose the biggest smaller integer index

#         start = time.time()
#         if pstore.takes_pointers():
#             for x in xrange(newshape[0]):
#                 for y in xrange(newshape[1]):
#                     fd = pstore.popen([x,y])
#                     if fd:
#                         pstore.pwrite(fd, 0, (indices[0][x][y], indices[1][x][y]))
#                         pstore.pclose(fd)
#         output = a[tuple(indices)]
#         end = time.time()
#         return output, {'provoverhead' : (end - start)}

#     def fmap_0(self, coord, run_id=None):
#         return map(lambda pair: (coord[0] * 2 + pair[0], coord[1] * 2 + pair[1]),
#             [(x,y) for y in xrange(2) for x in xrange(2)])


#     def bmap_0(self, coord, run_id=None):
#         return [(coord[0] / 2, coord[1] / 2)]
    

# class Regrid(Op):
#     def run(self, inputs, run_id):
#         pstore = self.pstore(run_id)

#         arr = inputs[0]
#         newshape = np.asarray(arr.shape) / 2
#         factor = np.asarray(arr.shape)/np.asarray(newshape)
#         mult = int(1 / factor)

#         evList = ['a.reshape('] + \
#                  ['newshape[%d],factor[%d],'%(i,i) for i in xrange(lenShape)] + \
#                  [')'] + ['.sum(%d)'%(i+1) for i in xrange(lenShape)] + \
#                  ['/factor[%d]'%i for i in xrange(lenShape)]
#         output = eval(''.join(evList))
        
#         start = time.time()
#         if pstore.takes_pointers():
#             for x in xrange(0, newshape[0], mult[0]):
#                 for y in xrange(0, newshape[1], mult[1]):
#                     fd = pstore.popen([x,y])
#                     if fd:
#                         for mx in xrange(mult[0]):
#                             for my in xrange(mult[1]):
#                                 pstore.pwrite(fd, 0, (x+mx, y+my))
#                     pstore.pclose(fd)
#         end = time.time()
#         return output, {'provoverhead' : (end - start)}



