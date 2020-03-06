from .ice cimport ICE
from .wing cimport Wing
from .misc import GRAVITY, isa
from . import conversions as conv
from math import pi, tan, atan, sqrt, nan, isnan, exp, log
from tabulate import tabulate


cdef class Aircraft():
    def __init__(self, ICE engine=None, Wing wing=None, double massEmpty=nan, double massAlternator=0, double massStarterMotor=0, double massBatteries=0, double payload=0, double fuelVolume=nan, double Cd0=nan, double Cl=nan, double efficiencyPropulsive=nan, double rhoSoc=1.25, double VTAS=nan, double VTASClimb=nan, double cruiseSpeed=nan, double rangeMax=nan,  str model="", str modelId="", int mode=1, **kwargs):
        self.engine = engine
        self.wing = wing
        assert massEmpty > engine.massDry, "Total aircraft mass must be greater than engine mass"
        self.massEmpty = massEmpty # [kg]
        self.massAlternator = massAlternator # [kg]
        self.massStarterMotor = massStarterMotor # [kg]
        self.massBatteries = massBatteries # [kg]
        self.payload = payload # [kg]
        self.fuelVolume = fuelVolume # [m^3]
        self.Cd0 = Cd0 # [-]
        self.Cl = Cl # [-]
        #self.Cd = Cd # [-]
        self.efficiencyPropulsive = efficiencyPropulsive # [-]
        self.rhoSoc = rhoSoc # [kg/m^3]
        self.VTAS = VTAS # [m/s]
        self.VTASClimb = VTASClimb
        self.cruiseSpeed = cruiseSpeed
        self.rangeMax = rangeMax # [m]
        self.model = model # full name of aircraft model
        self.modelId = modelId # short name of model, used as folder name for results
        self.mode = mode

    def __str__(self):
        return self.model

    cpdef Aircraft copy(self):
        return Aircraft(self.engine.copy(), self.wing.copy(), self.massEmpty, self.massAlternator, self.massStarterMotor, self.massBatteries, self.payload, self.fuelVolume, self.Cd0, self.Cl, self.efficiencyPropulsive, self.rhoSoc, self.VTAS, self.VTASClimb, self.cruiseSpeed, self.rangeMax, self.model, self.modelId, self.mode)            


    #------------------------------------------
    # Flight Conditions
    #------------------------------------------
    
    def set_rhoSoc(self, double altitude):
        self.rhoSoc = isa(altitude)['rho']

    cpdef public double q(self):
        return 0.5*self.rhoSoc*self.VTAS**2

    cpdef public double qS(self):
        return self.q()**self.wing.area
    
    #------------------------------------------
    # Coefficients
    #------------------------------------------
    
    cpdef public double Cd(self):
        "float: Coefficient of drag: CD = CD0 + k*CL**2"
        return self.Cd0 + self.wing.k()*self.Cl**2

    cpdef public void size_Cl(self):
        "Automatically set Cl from start-of-cruise weight, density and VTAS"
        self.Cl = self.weightTotalSoc()/(0.5*self.rhoSoc*self.VTAS**2*self.wing.area)
    
    #------------------------------------------
    # Speed
    #------------------------------------------

    cpdef public double Vmp(self):
        return (4/3*(self.weightTotalSoc()/self.wing.area)**2/self.rhoSoc**2*self.wing.k()/self.Cd0)**0.25

    cpdef public double Vmd(self):
        return (4*(self.weightTotalSoc()/self.wing.area)**2/self.rhoSoc**2*self.wing.k()/self.Cd0)**0.25

    cpdef public double Vcc(self):
        return self.Vmd()*2 - self.Vmp()

    def calc_VTASCruise(self, altitude, set_VTAS=True, set_cruiseSpeed=True):
        cdef double PavCruise = 0.75 * self.engine.powerRated * self.efficiencyPropulsive
        cdef vtas = (PavCruise/(0.5*self.rhoSoc*self.wing.area*self.Cd()))**(1./3)
        if set_VTAS:
            self.VTAS = vtas
        if set_cruiseSpeed:
            self.cruiseSpeed = vtas
        return vtas

    #------------------------------------------
    # Power
    #------------------------------------------
    
    cpdef public double Cp(self):
        "float: Power-based specific fuel consumption [N/s.W, N/J]"
        return self.engine.fuelConsumptionMass() * GRAVITY / self.powerShaft()
        
    cpdef public double powerICE(self):
        "float: Power from ICE [W]."
        return self.engine.powerRated
        
    cpdef public double powerShaft(self):
        "float: Power through shaft [W]."
        return self.engine.powerRated
    
    cpdef public double powerAvailable(self):
        "float: Power available [W]."
        return self.efficiencyPropulsive * self.powerShaft()
    
    cpdef public double powerRequired(self):
        "float: Power required [W]."
        return 0.5*self.rhoSoc*self.VTAS**3*self.wing.area*self.Cd()
    
    cpdef public double powerExcess(self):
        "float: Excess power [W]."
        return self.powerAvailable() - self.powerRequired()

    #------------------------------------------
    # Masses and Weights
    #------------------------------------------
    
    cpdef public double massICE(self):
        return self.engine.massDry
    
    cpdef public double weightICE(self):
        return self.massICE()*GRAVITY
    
    cpdef public double massPowertrain(self):
        "float: Mass of powertrain: ICE + alternator + starter motor."
        return self.engine.massDry + self.massAlternator + self.massStarterMotor
    
    cpdef public double weightPowertrain(self):
        return self.massPowertrain()*GRAVITY
    
    cpdef public double massStructure(self):
        "float: Dry mass of aircraft structure (excluding fuel); empty mass - powertrain - batteries."
        return self.massEmpty - self.massPowertrain()# - self.massBatteries
    
    cpdef public double weightStructure(self):
        return self.massStructure()*GRAVITY
    
    cpdef public double massFuelSoc(self):
        "float: Fuel mass of at start-of-cruise."
        return self.fuelVolume * self.engine.fuel.density
    
    cpdef public double weightFuelSoc(self):
        return self.massFuelSoc()*GRAVITY
    
    cpdef public double massTotalSoc(self):
        "float: Total mass of aircraft at start-of-cruise; empty mass + fuel mass + payload."
        return self.massEmpty + self.massFuelSoc() + self.payload
    
    cpdef public double weightTotalSoc(self):
        return self.massTotalSoc()*GRAVITY

    

    #-------------------------
    # Flight Mode 1
    #-------------------------
    cpdef double cruiseRange_1(self, double massFuelBurnt):
        cdef double a = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area*sqrt(self.Cd0/self.wing.k())
        cdef double Wsoc = self.massTotalSoc()*GRAVITY # start of cruise
        cdef double Wfb = massFuelBurnt*GRAVITY # end of cruise
        cdef double numerator = a * Wfb
        cdef double denominator = a**2 + Wsoc*(Wsoc-Wfb)
        return self.efficiencyPropulsive/self.Cp()/sqrt(self.Cd0*self.wing.k())*atan(numerator/denominator)
    
    cpdef double massFuelBurnt_1(self, double cruiseRange):
        cdef double a = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area*sqrt(self.Cd0/self.wing.k())
        cdef double Wsoc = self.massTotalSoc()*GRAVITY # start of cruise
        cdef double numerator = a**2 + Wsoc**2
        cdef double tan_ = tan(cruiseRange*self.Cp()*sqrt(self.Cd0*self.wing.k())/self.efficiencyPropulsive)
        cdef double denominator = a/tan_ + Wsoc
        return numerator/denominator/GRAVITY
        
    
    cpdef double volumeFuelBurnt_1(self, double cruiseRange):
        return self.massFuelBurnt_1(cruiseRange)/self.engine.fuel.density

    cpdef double massTotalEoc_1(self, double cruiseRange):
        cdef double dm = self.massFuelBurnt_1(cruiseRange)
        cdef double m1 = self.massTotalSoc()
        return m1 - dm
    
    #-------------------------
    # Flight Mode 2
    #-------------------------
    cpdef double cruiseRange_2(self, double massFuelBurnt):
        cdef double massTotalSoc = self.massTotalSoc()
        return (self.Cl/self.Cd())*self.efficiencyPropulsive/self.Cp()*log((massTotalSoc-massFuelBurnt)/massTotalSoc)
    
    cpdef double massFuelBurnt_2(self, double cruiseRange):
        return self.massTotalSoc()*(1-exp(-cruiseRange*self.Cp()/(self.Cl/self.Cd()*self.efficiencyPropulsive)))
        
    
    cpdef double volumeFuelBurnt_2(self, double cruiseRange):
        return self.massFuelBurnt_2(cruiseRange)/self.engine.fuel.density

    cpdef double massTotalEoc_2(self, double cruiseRange):
        cdef double dm = self.massFuelBurnt_2(cruiseRange)
        cdef double m1 = self.massTotalSoc()
        return m1 - dm
    
    #-------------------------
    # Flight Mode 3
    #-------------------------
    cpdef double cruiseRange_3(self, double massFuelBurnt):
        return 0
    
    cpdef double massFuelBurnt_3(self, double cruiseRange):
        return 0
        
    
    cpdef double volumeFuelBurnt_3(self, double cruiseRange):
        return 0

    cpdef double massTotalEoc_3(self, double cruiseRange):
        return 0
    
    #-------------------------
    # Auto Mode
    #-------------------------
    cpdef double cruiseRange(self, double massFuelBurnt):
        if self.mode == 1:
            return self.cruiseRange_1(massFuelBurnt)
        elif self.mode == 2:
            return self.cruiseRange_2(massFuelBurnt)
        elif self.mode == 3:
            return self.cruiseRange_3(massFuelBurnt)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))
    
    cpdef double massFuelBurnt(self, double cruiseRange):
        if self.mode == 1:
            return self.massFuelBurnt_1(cruiseRange)
        elif self.mode == 2:
            return self.massFuelBurnt_2(cruiseRange)
        elif self.mode == 3:
            return self.massFuelBurnt_3(cruiseRange)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))
        
    
    cpdef double volumeFuelBurnt(self, double cruiseRange):
        if self.mode == 1:
            return self.volumeFuelBurnt_1(cruiseRange)
        elif self.mode == 2:
            return self.volumeFuelBurnt_2(cruiseRange)
        elif self.mode == 3:
            return self.volumeFuelBurnt_3(cruiseRange)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))

    cpdef double massTotalEoc(self, double cruiseRange):
        if self.mode == 1:
            return self.massTotalEoc_1(cruiseRange)
        elif self.mode == 2:
            return self.massTotalEoc_2(cruiseRange)
        elif self.mode == 3:
            return self.massTotalEoc_3(cruiseRange)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))

    
    
    #-------------------------
    # Miscellaneous
    #-------------------------
    
    def summary(self, engine=True, wing=True, title=""):
        if not title:
            title = self.model + ' summary'
        print(title)
        table = [["massEmpty", self.massEmpty, "kg"],
                 ["massAlternator", self.massAlternator, "kg"],
                 ["massStarterMotor", self.massStarterMotor, "kg"],
                 ["massBatteries", self.massBatteries, "kg"],
                 ["payload", self.payload, "kg"],
                 ["massTotalSoc", self.massTotalSoc(), "kg"],
                 ["massICE", self.massICE(), "kg"],
                 ["massICE", self.massICE(), "kg"],
                 ["fuelVolume", self.fuelVolume, "m^3"],
                 ["", conv.m32L(self.fuelVolume), "L"],
                 ["massFuelSoc", self.massFuelSoc(), "kg"],
                 ["massStructure", self.massStructure(), "kg"],
                 ["", "", ""],
                 ["Cd0", self.Cd0, ""],
                 ["Cl", self.Cl, ""],
                 ["Cd", self.Cd(), ""],
                 ["", "", ""],
                 ["efficiencyPropulsive", self.efficiencyPropulsive, ""],
                 ["", "", ""],
                 ["Cp", self.Cp(), "N/s.W, N/J"],
                 ["", "", ""],
                 ["density", self.rhoSoc, "kg/m^3"],
                 ["VTAS", self.VTAS, "m/s"],
                 ["", conv.mps2knot(self.VTAS), "knot"],
                 ["VTASClimb", self.VTASClimb, "m/s"],
                 ["", conv.mps2knot(self.VTASClimb), "knot"],
                 ["Vmp", self.Vmp(), "m/s"],
                 ["", conv.mps2knot(self.Vmp()), "knot"],
                 ["Vmd", self.Vmd(), "m/s"],
                 ["", conv.mps2knot(self.Vmd()), "knot"],
                 ["Vcc", self.Vcc(), "m/s"],
                 ["", conv.mps2knot(self.Vcc()), "knot"],
                 ["cruise speed", self.cruiseSpeed, "m/s"],
                 ["", conv.mps2knot(self.cruiseSpeed
                 ), "knot"],
                 ["rangeMax", self.rangeMax/1000, "km"],
                 ["", conv.m2nmi(self.rangeMax), "nmi"]]
        print(tabulate(table))
        if engine:
            self.engine.summary()
        if wing:
            self.wing.summary()
