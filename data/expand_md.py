
import numpy as np
import sys

expand = int(sys.argv[1])


with file('MD_train_set.txt', 'r') as f:
    ds = np.array([l.strip().split('\t') for l in f])



header = ds[0,:].tolist()
data = ds[1:,4:].tolist()
names = ds[1:,:4].tolist()
npatients = len(data[0]) * expand

f = open('MD_train_set_%dx.txt' % expand, 'w')

header = header[:4]
header.extend(['patient_%d' % i for i in xrange(npatients)])

print >>f, '\t'.join(header)

for namerow, datarow in zip(names,data):
    row = []
    row.extend(namerow)
    for i in xrange(expand):
        row.extend(datarow)
    print >>f, '\t'.join(row)

f.close()
