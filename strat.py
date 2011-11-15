from cStringIO import StringIO

class Mode(object):
    QUERY = 1
    FULL_MAPFUNC = 1 << 1
    FULL_MAPFUNC_BOX = 1 << 2
    PT_MAPFUNC = 1 << 3
    PT_MAPFUNC_BOX = 1 << 4
    PTR = 1 << 5
    NOOP = 1 << 6
    STAT = 1 << 7

    @staticmethod
    def all():
        return [Mode.QUERY,
                Mode.FULL_MAPFUNC_BOX,
                Mode.FULL_MAPFUNC,
                Mode.PT_MAPFUNC_BOX,
                Mode.PT_MAPFUNC,
                Mode.PTR,
                Mode.NOOP,
                Mode.STAT]

    @staticmethod
    def ptype(mode):
        from provstore import PStore1, PStore2, PStore3
        if mode in (Mode.FULL_MAPFUNC_BOX, Mode.FULL_MAPFUNC):
            return PStore1
        if mode in (Mode.PT_MAPFUNC_BOX, Mode.PT_MAPFUNC):
            return PStore2
        return PStore3




# this is a tuple descriptor
class Spec(object):
    COORD_ONE  = 0
    COORD_MANY = 1
    #OUTCOORD_ONE  = 2
    #OUTCOORD_MANY = 3
    BOX        = 4
    KEY        = 5
    PRED       = 6
    BINARY     = 7  # serialized to byte string
    GRID       = 8
    NONE       = -1

    def __init__(self, outcoords, payload):
        self.outcoords = outcoords
        self.payload = payload

    @staticmethod
    def all():
        return (Spec.COORD_ONE, Spec.COORD_MANY, Spec.BOX, Spec.KEY, Spec.BINARY, Spec.NONE)

    @staticmethod
    def default():
        return Spec(Spec.NONE, Spec.NONE)

    @staticmethod
    def is_coord(spec):
        """
        Does the spec define an encoding that can be parsed into a set of
        coordinates?
        """
        return spec in (Spec.COORD_ONE, Spec.COORD_MANY, Spec.KEY)

    def copy(self):
        return Spec(self.outcoords, self.payload)


class Desc(object):
    def __init__(self, mode, spec, backward):
        self.mode = mode
        self.spec = spec
        self.backward = backward

    def copy(self):
        return Desc(self.mode, self.spec.copy(), self.backward)

    def __str__(self):
        s = ''
        if self.mode == Mode.NOOP:
            s = 'NOOP'
        elif self.mode == Mode.STAT:
            s = 'STAT'
        elif self.mode == Mode.FULL_MAPFUNC:
            s = "MAP"
        elif self.mode == Mode.PT_MAPFUNC:
            s = "PTMAP"
            x = []
            for foo in (self.spec.outcoords, self.spec.payload):
                if foo == Spec.KEY:
                    x.append('KEY')
                elif foo == Spec.COORD_MANY:
                    x.append("MANY")
                elif foo == Spec.COORD_ONE:
                    x.append("ONE")
                elif foo == Spec.BOX:
                    x.append("BOX")
                elif foo == Spec.GRID:
                    x.append("GRID")
                else:
                    x.append(str(foo))
            s = '%s_%s' % (s, "_".join(x))
        elif self.mode == Mode.QUERY:
            s = "Q"
        elif self.mode == Mode.PTR:
            x = []
            for foo in (self.spec.outcoords, self.spec.payload):
                if foo == Spec.KEY:
                    x.append('KEY')
                elif foo == Spec.COORD_MANY:
                    x.append("MANY")
                elif foo == Spec.COORD_ONE:
                    x.append("ONE")
                elif foo == Spec.BOX:
                    x.append("BOX")
                elif foo == Spec.GRID:
                    x.append("GRID")
                else:
                    x.append(str(foo))
            x.append(self.backward and 'b' or 'f')
            s = "_".join(x)

        else:                
            s = "%d_%d_%d_%d" % (self.mode,
                                 self.spec.outcoords,
                                 self.spec.payload,
                                 self.backward)
        return s
    
class Bucket(object):
    """
    Bucket describes a composite pstore
    """
    def __init__(self, descs):
        self.descs = descs

    def copy(self):
        return Bucket([desc.copy() for desc in self.descs])

class Strat(object):
    """
    Strat can contain multiple buckets
    """    
    def __init__(self, buckets):
        self.buckets = buckets

    @staticmethod
    def single(compmode, spec, backward=True):
        descs = []
        for mode in Mode.all():
            if mode & compmode:
                descs.append(Desc(mode, spec, backward))
        return Strat([Bucket(descs)])

    @staticmethod
    def noop():
        return Strat.single(Mode.NOOP, Spec.default(), True)

    @staticmethod
    def stat():
        return Strat.single(Mode.STAT, Spec.default(), True)

    @staticmethod
    def query():
        return Strat.single(Mode.QUERY, Spec.default(), True)

    @staticmethod
    def full():
        return Strat.single(Mode.FULL_MAPFUNC, Spec.default(), True)


    def modes(self):
        ret = set()
        for bucket in self.buckets:
            for desc in bucket.descs:
                ret.add(desc.mode)
        return list(ret)

    def copy(self):
        return Strat([bucket.copy() for bucket in self.buckets])

    def fnamestring(self):
        ll = []
        for buck in self.buckets:
            l = []
            for desc in buck.descs:
                l.append(str(desc))
            
            ll.append("__".join(l))
        return "___".join(ll)


    def __str__(self):
        return self.fnamestring()

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, o):
        return hash(o) == hash(self)

    def __cmp__(self, o):
        if self == o: return 0
        if str(self) > str(o):
            return 1
        return -1




# Strategy:
#     Bucket:
#         desc:
#             mode:
#                 query
#                 full map
#                 obj map
#                 full box
#                 obj box
#                 ptr
#             direction
#         composite = list(desc)
#     multi = list(bucket)
