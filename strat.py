from cStringIO import StringIO

class Mode(object):
    QUERY = 1
    FULL_MAPFUNC = 1 << 1
    FULL_MAPFUNC_BOX = 1 << 2
    PT_MAPFUNC = 1 << 3
    PT_MAPFUNC_BOX = 1 << 4
    PTR = 1 << 5

    @staticmethod
    def all():
        return [Mode.QUERY,
                Mode.FULL_MAPFUNC_BOX,
                Mode.FULL_MAPFUNC,
                Mode.PT_MAPFUNC_BOX,
                Mode.PT_MAPFUNC,
                Mode.PTR]

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
    OUTCOORD_ONE  = 2
    OUTCOORD_MANY = 3
    BOX        = 4
    KEY        = 5
    PRED       = 6
    BINARY     = 7  # serialized to byte string
    NONE       = -1

    def __init__(self, outcoords, payload):
        self.outcoords = outcoords
        self.payload = payload

    @staticmethod
    def is_coord(spec):
        """
        Does the spec define an encoding that can be parsed into a set of
        coordinates?
        """
        return spec in (Spec.COORD_ONE, Spec.COORD_MANY, Spec.KEY)


class Desc(object):
    def __init__(self, mode, spec, backward):
        self.mode = mode
        self.spec = spec
        self.backward = backward
    
class Bucket(object):
    """
    Bucket describes a composite pstore
    """
    def __init__(self, descs):
        self.descs = descs

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

    def fnamestring(self):
        ll = []
        for buck in self.buckets:
            l = []
            for desc in buck.descs:
                if desc.mode == Mode.FULL_MAPFUNC:
                    l.append("MAP")
                elif desc.mode == Mode.PT_MAPFUNC:
                    l.append("PTMAP")
                elif desc.mode == Mode.QUERY:
                    l.append("Q")
                elif desc.mode == Mode.PTR:
                    if desc.spec.payload == Spec.BOX:
                        l.append("BOX")
                    else:
                        x = []
                        for foo in (desc.spec.outcoords, desc.spec.payload):
                            if foo == Spec.KEY:
                                x.append('KEY')
                            elif foo == Spec.COORD_MANY:
                                x.append("MANY")
                            elif foo == Spec.COORD_ONE:
                                x.append("ONE")
                            elif foo == Spec.BOX:
                                x.append("BOX")
                            else:
                                x.append(str(foo))
                        l.append("_".join(x))
                else:                
                    l.append("%d_%d_%d_%d" % (desc.mode,
                                              desc.spec.outcoords,
                                              desc.spec.payload,
                                              desc.backward))
            
            ll.append("__".join(l))
        return "___".join(ll)




    def __str__(self):
        return self.fnamestring()
        buf = StringIO()
        
        for buck in self.buckets:
            buf.write("(")
            for desc in buck.descs:
                buf.write("[%d, %d / %d, %d]" % (desc.mode,
                                                 desc.spec.outcoords,
                                                 desc.spec.payload,
                                                 desc.backward))
            buf.write(")")
        return buf.getvalue()

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
