cdef class Wing():
    cdef public double span
    cdef public double area
    cdef public double _oswaldEfficiency
    #
    cpdef Wing copy(self)
    cpdef public double chord(self)
    cpdef public double aspectRatio(self)
    cpdef public double oswaldEfficiency(self)
    cpdef public double k(self)
