from .fuel cimport Fuel
from .misc import GRAVITY, isa
from . import conversions as conv
from math import pi, tan, atan, sqrt, nan, isnan, exp, log
from tabulate import tabulate

#--------------------------------------
# ICEngine
#--------------------------------------
cdef class ICE():
    def __init__(self, double massDry, double powerRated, double fuelConsumptionVolumetric, Fuel fuel=Fuel(718)):
        self.massDry = massDry #[kg]
        self.powerRated = powerRated #[W]
        self.fuelConsumptionVolumetric = fuelConsumptionVolumetric #[m^3/s]
        self.fuel = fuel

    cpdef public double fuelConsumptionMass(self):
        "float: mass-based fuel consumption [kg/s]"
        return self.fuelConsumptionVolumetric * self.fuel.density

    cpdef ICE copy(self):
        return ICE(self.massDry, self.powerRated, self.fuelConsumptionVolumetric, self.fuel)
    
    def summary(self, title="Engine summary"):
        print(title)
        table = [["massDry", self.massDry, "kg"],
                 ["powerRated", self.powerRated/1000, "kW"],
                 ["", conv.W2bhp(self.powerRated), "bhp"],
                 ["fuelConsumptionVolumetric", self.fuelConsumptionVolumetric, "m^3/s"],
                 ["fuelConsumptionVolumetric", conv.m32L(self.fuelConsumptionVolumetric), "L/s"],
                 ["fuelConsumptionMass", self.fuelConsumptionMass(), "kg/s"]]
        print(tabulate(table))
        if self.fuel:
            self.fuel.summary()
