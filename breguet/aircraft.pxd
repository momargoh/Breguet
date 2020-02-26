from .ice cimport ICE
from .wing cimport Wing

cdef class Aircraft():
    cdef public double massEmpty
    cpdef public double massAlternator
    cpdef public double massStarterMotor
    cpdef public double massBatteries
    cdef public double payload
    cdef public double fuelVolume # [L] not m^3
    cdef public double Cl
    cdef public double Cd0
    cdef public double rhoSoc # [kg/m^3]
    cdef public double VTAS # [m/s]
    cdef public double VTASClimb # [m/s]
    cdef public double cruiseSpeed # [m]
    cdef public double rangeMax # [m]
    cdef public ICE engine
    cdef public Wing wing
    cdef public double efficiencyPropulsive
    cdef public str model
    cdef public str modelId
    cdef public int mode
    #
    cpdef Aircraft copy(self)
    #
    cpdef public double q(self)
    cpdef public double qS(self)
    #    
    cpdef public double Cd(self)
    cpdef public void auto_Cl(self)
    #
    cpdef public double Vmp(self)
    cpdef public double Vmd(self)
    cpdef public double Vcc(self)
    #
    cpdef public double Cp(self)
    cpdef public double powerICE(self)
    cpdef public double powerShaft(self)
    cpdef public double powerAvailable(self)
    cpdef public double powerRequired(self)
    cpdef public double powerExcess(self)
    #
    cpdef public double massICE(self)
    cpdef public double massPowertrain(self)
    cpdef public double massStructure(self)
    cpdef public double massFuelSoc(self)
    cpdef public double massTotalSoc(self)
    #
    cpdef public double weightICE(self)
    cpdef public double weightPowertrain(self)
    cpdef public double weightStructure(self)
    cpdef public double weightFuelSoc(self)
    cpdef public double weightTotalSoc(self)
    #
    cpdef double cruiseRange_1(self, double massFuelBurnt)
    cpdef double massFuelBurnt_1(self, double cruiseRange)
    cpdef double volumeFuelBurnt_1(self, double cruiseRange)
    cpdef double massTotalEoc_1(self, double cruiseRange)
    #
    cpdef double cruiseRange_2(self, double massFuelBurnt)
    cpdef double massFuelBurnt_2(self, double cruiseRange)
    cpdef double volumeFuelBurnt_2(self, double cruiseRange)
    cpdef double massTotalEoc_2(self, double cruiseRange)
    #
    cpdef double cruiseRange_3(self, double massFuelBurnt)
    cpdef double massFuelBurnt_3(self, double cruiseRange)
    cpdef double volumeFuelBurnt_3(self, double cruiseRange)
    cpdef double massTotalEoc_3(self, double cruiseRange)
    #
    cpdef double cruiseRange(self, double massFuelBurnt)
    cpdef double massFuelBurnt(self, double cruiseRange)
    cpdef double volumeFuelBurnt(self, double cruiseRange)
    cpdef double massTotalEoc(self, double cruiseRange)
 
