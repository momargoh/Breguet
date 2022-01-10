from .aircraft cimport Aircraft
from .ice cimport ICE

cdef class HybridAircraft():
    cpdef public Aircraft originalAircraft
    cpdef public double downsizingFactor
    cpdef public double massHarvester
    cpdef public double efficiencyHarvester
    cpdef public double _specificPowerHarvester
    cpdef public double massGenerator
    cpdef public double efficiencyGenerator
    cpdef public double _specificPowerGenerator
    cpdef public double massRectifier
    cpdef public double efficiencyRectifier
    cpdef public double _specificPowerRectifier
    cpdef public double massInverter
    cpdef public double efficiencyInverter
    cpdef public double _specificPowerInverter
    cpdef public double massMotor
    cpdef public double efficiencyMotor
    cpdef public double _specificPowerMotor
    cpdef public double massBatteries
    cpdef public double efficiencyBatteries
    cpdef public double _specificPowerBatteries
    cpdef public double efficiencyICE
    cpdef public double efficiencyExhaust
    # Volume variables
    cpdef public double densityGenerator
    cpdef public double densityMotor
    cpdef public double densityRectifier
    cpdef public double densityInverter
    cpdef public double densityBatteries
    #
    cpdef public str modelIdPrefix
    cpdef public void update(self, kwargs)
    #
    cdef public ICE engine
    cpdef public double massEmpty(self)
    cdef public double fuelVolume
    #
    cpdef public void rebuild(self)
    cpdef public void size_Cl(self)
    #
    cpdef public double efficiencyHarvestingSystem(self)
    #
    cpdef public double powerExcessClimb(self)
    cpdef public double powerClimb(self)
    cpdef public double powerShaft(self)
    cpdef public double powerTotal(self)
    cpdef public double powerICE(self)
    cpdef public double powerExhaust(self)
    cpdef public double powerHarvester(self)
    cpdef public double powerGenerator(self)
    cpdef public double powerRectifier(self)
    cpdef public double powerInverter(self)
    cpdef public double powerInverterClimb(self)
    cpdef public double powerMotor(self)
    cpdef public double powerMotorClimb(self)
    cpdef public double powerBatteries(self)
    cpdef public double powerEGHS(self)
    cpdef public double powerAvailable(self)
    cpdef public double powerRequired(self)
    cpdef public double Cp(self)
    #
    cpdef public double volumeEGHS(self)
    cpdef public double volumeHarvester(self)
    cpdef public double volumeGenerator(self)
    cpdef public double volumeRectifier(self)
    cpdef public double volumeInverter(self)
    cpdef public double volumeBatteries(self)
    cpdef public double volumeMotor(self)
    #
    cpdef public double specificPowerHarvester(self)
    cpdef public double specificPowerGenerator(self)
    cpdef public double specificPowerRectifier(self)
    cpdef public double specificPowerInverter(self)
    cpdef public double specificPowerInverterClimb(self)
    cpdef public double specificPowerBatteries(self)
    cpdef public double specificPowerMotor(self)
    cpdef public double specificPowerMotorClimb(self)
    cpdef public double specificPowerEGHS(self)
    #
    cpdef public double massICE(self)
    cpdef public double massEGHS(self)
    cpdef public double massPowertrain(self)
    cpdef public double massEmpty(self)
    cpdef public double massStructure(self)
    cpdef public double massFuelSoc(self)
    cpdef public double massTotalSoc(self)
    #
    cpdef public double weightEGHS(self)
    #cpdef public double weightMotor(self)
    cpdef public double weightICE(self)
    cpdef public double weightPowertrain(self)
    cpdef public double weightOther(self)
    cpdef public double weightFuelSoc(self)
    cpdef public double weightTotalSoc(self)
    #
    cpdef public double size_generator(self)
    cpdef public double size_rectifier(self)
    cpdef public double size_batteries(self)
    cpdef public double size_inverter(self)
    cpdef public double size_motor(self)
    cpdef public void size_EGHS(self)
    #
    cpdef double downsizingFactorMax_1(self)
    cpdef double rangeX_1(self, double massFuelBurnt)
    cpdef double massFuelBurnt_1(self, double rangeX)
    cpdef double volumeFuelBurnt_1(self, double rangeX)
    cpdef double massTotalEoc_1(self, double rangeX)
    #
    cpdef double downsizingFactorMax_2(self)
    cpdef double rangeX_2(self, double massFuelBurnt)
    cpdef double massFuelBurnt_2(self, double rangeX)
    cpdef double volumeFuelBurnt_2(self, double rangeX)
    cpdef double massTotalEoc_2(self, double rangeX)
    #
    cpdef double downsizingFactorMax_3(self)
    cpdef double rangeX_3(self, double massFuelBurnt)
    cpdef double massFuelBurnt_3(self, double rangeX)
    cpdef double volumeFuelBurnt_3(self, double rangeX)
    cpdef double massTotalEoc_3(self, double rangeX)
    #
    cpdef double downsizingFactorMax(self)
    cpdef double rangeX(self, double massFuelBurnt)
    cpdef double massFuelBurnt(self, double rangeX)
    cpdef double volumeFuelBurnt(self, double rangeX)
    cpdef double massTotalEoc(self, double rangeX)
    #
    cpdef double downsizingFactorMaxSeries(self)
    cpdef double rangeXSeries(self, double massFuelBurnt)
    cpdef double massFuelBurntSeries(self, double rangeX)
    cpdef double volumeFuelBurntSeries(self, double rangeX)
    cpdef double massTotalEocSeries(self, double rangeX)
    cpdef double shaftToICERatio(self)
    cpdef double shaftToICERatioSeries(self)

