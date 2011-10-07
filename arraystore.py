import os, numpy, hashlib, struct

class ArrayStore(object):
    """
    Simulates the SciDB array storage
    """
    
    _instance = None
    def __init__(self):
        self.dirname = '_arraystore'
        if not os.path.exists('./%s'%self.dirname):
            os.mkdir("./%s"%self.dirname)

        # in memory list of keys of arrays i own        
        self.store = set()
        for name in os.listdir('./%s'%self.dirname):
            if name.endswith(".npy"):
                self.store.add(name[:-4])
        # in memory map of arrays to their shape
        self.shapes = {}
        

    @staticmethod
    def instance():
        if ArrayStore._instance is None:
            ArrayStore._instance = ArrayStore()
        return ArrayStore._instance

    def size(self, key):
        if key not in self.store:
            if not os.path.exists('./%s/%s.npy' % (self.dirname, key)):
                raise Exception('array not found.  key: %s' % key)
        return os.path.getsize("./%s/%s.npy" % (self.dirname, key))

    def get(self, key):
        if key not in self.store:
            if not os.path.exists('./%s/%s.npy' % (self.dirname, key)):
                raise Exception('array not found.  key: %s' % key)
        return numpy.load("./%s/%s.npy" % (self.dirname, key))
        f = file('./%s/%s.npy' % (self.dirname, key), 'rb')
        typeid, ndims = struct.unpack('II', f.read(struct.calcsize('II')))
        shape = struct.unpack('%dI' % ndims, f.read(struct.calcsize('%dI' % ndims)))
        rowsize = shape[1]
        if typeid == 0:
            dtype,fmt = float, '%df'
        elif typeid == 1:
            dtype,fmt = int, '%dI'
        elif typeid == 2:
            dtype,fmt = bool, '%d?'
        else:
            raise RuntimeError
        fmt = fmt%rowsize
        
        arr = numpy.empty(shape, dtype)
        for i in xrange(shape[0]):
            arr[i] = struct.unpack(fmt, f.read(struct.calcsize(fmt)))
        f.close()
        return arr

        return self.store.get(key, [])

    def add(self, arr):
        # if doen't exist, create file, add to self.store
        
        chksum = self.checksum(arr)
        if chksum not in self.store:
            numpy.save("./%s/%s.npy"%(self.dirname, chksum), arr)
            return chksum
            if arr.dtype == float:
                typeid,fmt = 0,"%df"
            elif arr.dtype == int:
                typeid,fmt = 1,"%dI"
            elif arr.dtype == bool:
                typeid,fmt = 2,"%d?"
            else:
                raise RuntimeError
            
            f = file('./%s/%s.npy' % (self.dirname, chksum), 'wb')
            shape = arr.shape
            f.write(struct.pack("%dI" % (len(shape)+2), typeid, len(shape), *shape))

            
            for row in arr:
                s = struct.pack(fmt % len(row), *row)
                f.write(s)
                del s
            f.close()
            

            self.store.add(chksum)
        self.shapes[chksum] = arr.shape
        return chksum

    def shape(self, key):
        if key not in self.shapes:
            shape = self.get(key).shape
            self.shapes[key] = shape
        return self.shapes[key]

    def checksum(self, arr):
        # md5 faster than sha1
        try:
            return hashlib.md5(arr.data).hexdigest()
        except:
            return hashlib.md5(arr.copy().data).hexdigest()



if __name__ == '__main__':
    import sys
    store = ArrayStore.instance()
    for size in [100, 500, 1000, 2000]:
        arr = numpy.zeros((size,size), float)
        key = store.add(arr)
        print size, store.size(key)
    
    # 8.00002 bytes per cell
