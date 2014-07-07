import bsddb
import sys, os, math
from cStringIO import StringIO
sys.path.append('../')
from stats import *
import matplotlib.pyplot as plt
import matplotlib




for j in [1000]:# [1, 10, 100, 1000]:
    for i in [1, 10, 100, 1000]:
        for n in (1, 10, 100, 1000, 2000, 5000, 10000):
            osize = 0
            isize = 0
            cost = 0.0
            db = bsddb.hashopen('/tmp/test.db', 'n')
            iencs = [j] * j
            ibuf = StringIO()
            ibuf.write(struct.pack('%di' % len(iencs), *iencs))
            iser = ibuf.getvalue()
            isize = len(iser)
            for x in xrange(n):
                oencs = [x] * i
                obuf = StringIO()
                obuf.write(struct.pack('%di' % len(oencs), *oencs))
                oser = obuf.getvalue()
                if not osize: osize = len(oser)
                start = time.time()
                db[oser] = iser
                end = time.time()
                cost += (end - start)

            db.close()

            dbsize = os.path.getsize('/tmp/test.db')
            estsize = math.ceil(float((osize + isize) * n) / 4096.0) * 4096
            data = (n, i, j, cost * 1000, osize, isize, dbsize, estsize, dbsize - estsize, dbsize / estsize)
            print '\t'.join(['%d']*len(data)) % data



exit()
Stats.instance('_output/pablostats.db')
path = sys.argv[1]
for fname in os.listdir(path):
    db = bsddb.hashopen('%s/%s' % (path, fname), 'r')
    keysize, valsize = 0, 0
    for key, val in db.iteritems():
        keysize += len(key)
        valsize += len(val)
    print "%d\t%d\t%s" % (keysize, valsize, fname)

    

    
    
    db.close()
