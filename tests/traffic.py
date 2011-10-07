from op import *
import numpy


class SegDrives(Op):
    
    def run(self, inputs, run_id):
        """
        Generates a 2d jagged array of
        output[segmentid][pointidx] = (lat, lon, tstamp, segmentid)
        """
        pstore = Runtime.instance().get_provstore(self, run_id)
        
        points = []
        previd = None
        output = []
        cid = 0
        gps = inputs[0]


        for idx, (id,ts,lat,lon) in enumerate(gps):

            if id != previd and len(points) > 0:
                chunks = self.create_chunks(points)

                for c in chunks:
                    for ptidx, pt in enumerate(c):
                        fd = pstore.popen([cid, ptidx])
                        pstore.pwrite(fd, 0, [pt[-1]])
                        pt[-1] = cid
                        pstore.pclose(fd)                        
                    output.append(c)
                    cid += 1

                points = []
            points.append((lat, lon, ts, idx))
            previd = id
        return output, {}

    def create_chunks(self, points):
        TIME_THRESH = 60 * 2 # maximum time between positions to cause a new trace to be generated
        MIN_SEGMENT_TIME = 30
        prevt = 0
        curidx = 0
        chunks = []
        chunk = []
        for idx, (lat,lon,t, inidx) in enumerate(points):
            if t > prevt:
                if t > prevt + TIME_THRESH and prevt != 0:
                    chunks.append(chunk)
                    chunk = []
                chunk.append([lat,lon,t, inidx])
                prevt = t
        if len(chunk) > 1 and chunk[-1][2] - chunk[0][2] > MIN_SEGMENT_TIME:
            chunks.append(chunk)
        return chunks

class FakeSeg(Op):
    def run(self, inputs, run_id):
        """
        return jagged array of segments.  each segment is up to 5 pts long
        output[segmentid][ptidx] = (lat, lon, t, segmentid)
        """
        pstore = Runtime.instance().get_provstore(self, run_id)

        inputs = inputs[0]
        output = [[]]

        fd = None
        segid = 0
        for idx, (id, t, lat, lon) in enumerate(inputs):
            if idx and idx % 5 == 0:
                output.append([])

            fd = pstore.popen([len(output)-1, len(output[-1])])
            pstore.pwrite(fd, 0, [idx])
            pstore.pclose(fd)
            output[-1].append([lat, lon, t, segid])

        if not output[-1]: output.pop(len(output)-1)
        return output, {}
            
        

class SumOp(Op):
    def run(self, inputs, run_id):
        pstore = Runtime.instance().get_provstore(self, run_id)

        inputs = inputs[0]
        output = []

        for cid, inputs in enumerate(inputs):
            fd = pstore.popen([cid])
            for ptid, (lat, lon, t, x) in enumerate(inputs):
                pstore.pwrite(fd, 0, [cid, ptid])
            pstore.pclose(fd)
            output.append(sum(map(lambda x:x[0], inputs)))
        return output, {}


class AvgOp(Op):
    def run(self, inputs, run_id):
        pstore = Runtime.instance().get_provstore(self, run_id)

        inputs = inputs[0]
        output = []

        for cid, inputs in enumerate(inputs):
            fd = pstore.popen([cid])
            for ptid, (lat, lon, t, x) in enumerate(inputs):
                pstore.pwrite(fd, 0, [cid, ptid])
            pstore.pclose(fd)
            output.append(numpy.mean(map(lambda x:x[0], inputs)))
        return output, {}



    


if __name__ == '__main__':
    op = FakeSeg()#SegDrives()
    sumop = SumOp()
    avgop = AvgOp()
    w = Workflow()
    w.register(op, 1)
    w.register(sumop, 1)
    w.register(avgop, 1)    
    w.connect(op, sumop, 0)
    w.connect(op, avgop, 0)    


    from datetime import datetime
    from time import mktime
    print "Processing"
    f = open("/Users/sirrice/mitnotes/research/provenance/src/traffic/tmp.data","r")

    gpspoints = []
    N = 1000
    for idx, line in enumerate(f):
        try:
            data = line.split(",")
            if (len(data) != 10):
                    print "BAD LINE: ", data
                    continue
            (id,ts,lat,lon,spd,hdg,servertime,error,devid,plat) = data
            timeStamp = datetime.strptime(ts,"%Y-%m-%dT%H:%M:%SZ")
            ts = int(mktime(timeStamp.timetuple()))
            gpspoints.append((id,ts,float(lat),float(lon)))#,spd,hdg,servertime,error,devid,plat))
        except Exception as e:
            #print e
            continue
        if idx >= N: break
    f.close()

    # gpspoints = []
    # for i in xrange(20):
    #     gpspoints.append((i,i,i,i))

    print "RUN"
    w.run([(op, [gpspoints[:10]])])
    #w.run([(op, [gpspoints[10:]])])    
    print  w.backward([[0]], sumop, 0, 1)
    print  w.backward([[0]], avgop, 0, 1)    
    print w.forward([[[i] for i in xrange(4,7,1)]], op, 0, 1)
