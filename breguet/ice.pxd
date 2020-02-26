from .fuel cimport Fuel
cdef class ICE():
    cdef public double massDry
    cdef public double powerRated
    cdef public double fuelConsumptionVolumetric
    cdef public Fuel fuel
    #
    cpdef ICE copy(self)
    cpdef public double fuelConsumptionMass(self)
