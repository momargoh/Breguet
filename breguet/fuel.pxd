cdef class Fuel():
    cpdef public double density
    cpdef public str name
    cpdef public double specificVolume(self)
    cpdef Fuel copy(self)
