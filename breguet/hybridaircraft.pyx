from .aircraft cimport Aircraft
from .ice cimport ICE
from .wing cimport Wing
from .misc import GRAVITY, isa
from . import conversions as conv
from math import pi, tan, atan, sqrt, nan, isnan, exp, log
from tabulate import tabulate
import scipy.optimize as opt

cdef class HybridAircraft():

    def __init__(self, Aircraft originalAircraft, double downsizingFactor, double efficiencyICE, double efficiencyExhaust, double efficiencyHarvester=1, double massHarvester=0, double specificPowerHarvester=0, double efficiencyGenerator=0, double massGenerator=0, double specificPowerGenerator=0, double efficiencyRectifier=0, double massRectifier=0, double specificPowerRectifier=0, double efficiencyInverter=0, double massInverter=0, double specificPowerInverter=0, double efficiencyMotor=0, double massMotor=0, double specificPowerMotor=0, double efficiencyBatteries=0, double massBatteries=0, double specificPowerBatteries=0, double densityGenerator=0, double densityMotor=0, double densityBatteries=0, double densityRectifier=0, double densityInverter=0, str modelIdPrefix='HYB:'):
        """Set massX or specificPowerX. If both are set, specificPowerX will be used to calculate the required massX"""
        self.originalAircraft = originalAircraft
        self.modelIdPrefix = modelIdPrefix
        self.downsizingFactor = downsizingFactor
        self.efficiencyICE = efficiencyICE
        self.efficiencyExhaust = efficiencyExhaust
        #
        self.efficiencyHarvester = efficiencyHarvester
        self.massHarvester = massHarvester 
        self._specificPowerHarvester = specificPowerHarvester
        self.efficiencyGenerator = efficiencyGenerator
        self.massGenerator = massGenerator
        self._specificPowerGenerator = specificPowerGenerator
        self.efficiencyRectifier = efficiencyRectifier
        self.massRectifier = massRectifier
        self._specificPowerRectifier = specificPowerRectifier
        self.efficiencyInverter = efficiencyInverter
        self.massInverter = massInverter
        self._specificPowerInverter = specificPowerInverter
        self.efficiencyMotor = efficiencyMotor
        self.massMotor = massMotor
        self._specificPowerMotor = specificPowerMotor
        self.efficiencyBatteries = efficiencyBatteries
        self.massBatteries = massBatteries
        self._specificPowerBatteries = specificPowerBatteries
        #
        self.specificVolumeGenerator = densityGenerator
        self.specificVolumeMotor = densityMotor
        self.specificVolumeRectifier = densityRectifier
        self.specificVolumeInverter = densityInverter
        self.specificVolumeBatteries = densityBatteries
        #
        cdef ICE eng = originalAircraft.engine
        cdef double massDry = eng.massDry*(1-downsizingFactor)
        cdef double powerRated  = eng.powerRated*(1-downsizingFactor)
        cdef double fcv  = eng.fuelConsumptionVolumetric*(1-downsizingFactor)
        self.engine = ICE(massDry=massDry, powerRated=powerRated, fuelConsumptionVolumetric=fcv, fuel=eng.fuel)
        self.fuelVolume = originalAircraft.fuelVolume * (1-downsizingFactor)


    cpdef public double efficiencyHarvestingSystem(self):
        return self.efficiencyHarvester*self.efficiencyGenerator*self.efficiencyRectifier*self.efficiencyInverter*self.efficiencyMotor
    
    cpdef public double size_generator(self):
        self.massGenerator = self.powerHarvester() / self._specificPowerGenerator
        return self.massGenerator
    
    cpdef public double size_rectifier(self):
        self.massRectifier = self.powerGenerator() / self._specificPowerRectifier
        return self.massRectifier
    
    cpdef public double size_inverter(self):
        self.massInverter = (self.powerRectifier()+self.powerBatteries())/self._specificPowerInverter
        return self.massInverter
    
    cpdef public double size_motor(self):
        self.massMotor = self.powerInverterClimb()/self._specificPowerMotor
        return self.massMotor

    """
    # calculate batteries iteratively using brentq
    cpdef public double solver_massBatteries(self, mBat):
        self.massBatteries = mBat
        self.calc_massMotor() #set motor size to handle extra battery power
        cdef powerExcessClimb1 = self.powerExcessClimb() #from massTotalSoc & VTASClimb
        cdef powerExcessClimb2 = self.powerBatteries() #from battery power only
        return powerExcessClimb1 - powerExcessClimb2


    cpdef double calc_massBatteries(self, bounds):
        cdef double solvedMassBatteries = opt.brentq(
            self.solver_massBatteries,
                        bounds[0],
                        bounds[1])
        self.massBatteries = solvedMassBatteries
        return solvedMassBatteries
    """
    
    cpdef double size_batteries(self):
        cdef double denominator = self.efficiencyBatteries*self._specificPowerBatteries*(1/self.VTASClimb/GRAVITY-1/self._specificPowerMotor)
        cdef double numerator = self.payload + self.massFuelSoc() + self.massStructure() + self.massICE() + self.massHarvester + self.massGenerator + self.massRectifier + self.massInverter + self.powerInverter()/self._specificPowerMotor
        self.massBatteries = numerator/denominator
        return self.massBatteries
    
    cpdef public void size_EGHS(self):
        if self._specificPowerGenerator:
            self.size_generator()
        if self._specificPowerRectifier:
            self.size_rectifier()
        if self._specificPowerBatteries:
            self.size_batteries()
        if self._specificPowerInverter:
            self.size_inverter()
        if self._specificPowerMotor:
            self.size_motor()

            
    cpdef public void update(self, kwargs):
        """Update (multiple) class variables from a dictionary of keyword arguments.

Parameters
-----------
kwargs : dict
    Dictionary of attributes and their updated value."""
        cdef str key
        cdef list key_split
        for key, value in kwargs.items():
            if "." in key:
                key_split = key.split(".", 1)
                key_attr = getattr(self, key_split[0])
                key_attr.update({key_split[1]: value})
            elif '[' in key and ']' in key:
                key = key.replace('[',']')
                key_split = key.split(']')
                key_split.remove('')
                if len(key_split) > 2:
                    key_attr = getattr(self, key_split[0])[int(key_split[1])]
                    key_attr.update({key_split[2]: value})
                else:
                    key_attr = getattr(self, key_split[0])
                    key_attr[int(key_split[1])] = value
            else:
                setattr(self, key, value)

                
    cpdef public void rebuild(self):
        "void: 'Rebuilds' the hybrid aircraft: should be called if downsizingFactor or originalAircraft.engine is changed."
        self.engine.massDry = self.originalAircraft.engine.massDry*(1-self.downsizingFactor)
        self.engine.fuelConsumptionVolumetric = self.originalAircraft.engine.fuelConsumptionVolumetric*(1-self.downsizingFactor)
        self.engine.powerRated = self.originalAircraft.engine.powerRated*(1-self.downsizingFactor)
        self.fuelVolume = self.originalAircraft.fuelVolume * (1-self.downsizingFactor)
        #self.size_EGHS()
        
    cpdef public double Cp(self):
        "float: Power-based specific fuel consumption [N/s.W, N/J]"
        return self.engine.fuelConsumptionMass() * GRAVITY / self.powerShaft()
        

    cpdef public void size_Cl(self):
        "Automatically set Cl from start-of-cruise weight, density and speed"
        self.originalAircraft.size_Cl()

    cpdef public double powerExcessClimb(self):
        return self.weightTotalSoc()*self.VTASClimb
    

    cpdef public double powerClimb(self):
        return self.powerExcessClimb() + self.powerRequired()
    
    cpdef public double powerICE(self):
        "float: Power from ICE [W]."
        return self.engine.powerRated
    
    cpdef public double powerShaft(self):
        "float: Power through shaft [W]."
        # self.powerICE() is for hybrid, already takes into account (1-downsizingFactor)
        return self.powerICE() + self.powerEGHS()
    
    cpdef public double powerTotal(self):
        return self.powerICE() / self.efficiencyICE
    
    cpdef public double powerHarvester(self):
        "float: Power through exhaust gases [W]."
        return self.efficiencyHarvester * self.powerExhaust()
    
    cpdef public double powerGenerator(self):
        "float: Power provided by generator [W]."
        return self.efficiencyGenerator * self.powerHarvester()
    
    cpdef public double powerRectifier(self):
        "float: Power provided by rectifier [W]."
        return self.efficiencyRectifier * self.powerGenerator()
    
    cpdef public double powerInverter(self):
        "float: Power provided by Inverter [W]."
        return self.efficiencyInverter * self.powerRectifier()
    
    cpdef public double powerInverterClimb(self):
        "float: Power provided by Inverter [W]."
        return self.efficiencyInverter * (self.powerRectifier() + self.powerBatteries())
    
    cpdef public double powerMotor(self):
        "float: Cruise power provided by motor [W]."
        return self.efficiencyMotor * self.powerInverter()
    
    cpdef public double powerMotorClimb(self):
        "float: Maximum power provided by batteries (including excess for climb) [W]."
        return self.efficiencyMotor * self.powerInverterClimb()
    
    cpdef public double powerBatteries(self):
        "float: Power provided by batteries [W]."
        return self.efficiencyBatteries * self.massBatteries * self._specificPowerBatteries
    
    cpdef public double powerExhaust(self):
        "float: Power through exhaust gases [W]."
        return self.efficiencyExhaust * self.powerTotal()
    
    cpdef public double powerEGHS(self):
        "float: Power output of exhaust gas harvesting system [W]."
        return self.powerMotor()
    
    cpdef public double powerAvailable(self):
        return self.efficiencyPropulsive * self.powerShaft()
    
    cpdef public double powerRequired(self):
        return 0.5*self.rhoSoc*self.VTAS**3*self.wing.area*self.Cd()
    
    cpdef public double specificPowerHarvester(self):
        "float: Specific power of harvesting system [W/kg]."
        return self.powerHarvester()/self.massHarvester
    
    cpdef public double specificPowerGenerator(self):
        "float: Specific power of generator [W/kg]."
        return self.powerGenerator()/self.massGenerator
    
    cpdef public double specificPowerRectifier(self):
        "float: Specific power of rectifier [W/kg]."
        return self.powerRectifier()/self.massRectifier
    
    cpdef public double specificPowerBatteries(self):
        "float: Specific power of batteries [W/kg]."
        return self.powerBatteries()/self.massBatteries
    
    cpdef public double specificPowerInverter(self):
        "float: Specific power of converter in cruise [W/kg]."
        return self.powerInverter()/self.massInverter
    
    cpdef public double specificPowerInverterClimb(self):
        "float: Specific power of converter in climb [W/kg]."
        return self.powerInverterClimb()/self.massInverter
    
    cpdef public double specificPowerMotor(self):
        "float: Specific power of motor in cruise [W/kg]."
        return self.powerMotor()/self.massMotor
    
    cpdef public double specificPowerMotorClimb(self):
        "float: Specific power of motor in climb [W/kg]."
        return self.powerMotorClimb()/self.massMotor
    
    cpdef public double specificPowerEGHS(self):
        "float: Specific power of harvesting system [W/kg]."
        return self.powerEGHS()/self.massEGHS()
    
    cpdef public double massEmpty(self):
        return self.originalAircraft.massStructure() + self.massICE() + self.massEGHS()
    
    cpdef public double massEGHS(self):
        return self.massHarvester + self.massGenerator + self.massRectifier + self.massInverter + self.massMotor + self.massBatteries
    
    cpdef public double weightEGHS(self):
        return self.massEGHS()*GRAVITY

    cpdef public double massICE(self):
        return self.engine.massDry
    
    cpdef public double weightICE(self):
        return self.massICE*GRAVITY
    
    cpdef public double massPowertrain(self):
        "float: Dry mass of downsized engine; original engine dry mass * downsizing factor."
        return self.massICE() + self.massEGHS()
    
    cpdef public double weightPowertrain(self):
        return self.massPowertrain()*GRAVITY
    
    cpdef public double massStructure(self):
        return self.originalAircraft.massStructure() # other masses are unaffected by the hybridisation
    
    cpdef public double weightOther(self):
        return self.massOther()*GRAVITY
    
    cpdef public double massFuelSoc(self):
        "float: Fuel mass of at start of cruise."
        return self.fuelVolume * self.engine.fuel.density
    
    cpdef public double weightFuelSoc(self):
        return self.massFuelSoc()*GRAVITY
    
    cpdef public double massTotalSoc(self):
        "float: Total mass of hybrid aircraft at take-off; original take-off weight * downsizing factor + harvesting mass."
        return self.massEmpty() + self.massFuelSoc() + self.payload
    
    cpdef public double weightTotalSoc(self):
        return self.massTotalSoc()*GRAVITY

    #-------------------------
    # Volumes
    #-------------------------
    cpdef public double volumeEGHS(self):
        return self.volumeHarvester() + self.volumeGenerator() + self.volumeRectifier() + self.volumeInverter() + self.volumeMotor() + self.volumeBatteries()

    cpdef public double volumeHarvester(self):
        return 0
    cpdef public double volumeGenerator(self):
        return self.massGenerator*self.specificVolumeGenerator
    cpdef public double volumeRectifier(self):
        return self.massRectifier*self.specificVolumeRectifier
    cpdef public double volumeInverter(self):
        return self.massInverter*self.specificVolumeInverter
    cpdef public double volumeBatteries(self):
        return self.massBatteries*self.specificVolumeBatteries
    cpdef public double volumeMotor(self):
        return self.massMotor*self.specificVolumeMotor
    

            
    #-------------------------
    # Flight Mode 3
    #-------------------------
    cpdef double downsizingFactorMax_3(self):
        return 0
    
    cpdef double rangeX_3(self, double massFuelBurnt):
        return 0
    
    cpdef double massFuelBurnt_3(self, double rangeX):
        return 0
        
    
    cpdef double volumeFuelBurnt_3(self, double rangeX):
        return 0

    cpdef double massTotalEoc_3(self, double rangeX):
        return 0
    
    #-------------------------
    # Flight Mode 2
    #-------------------------
    
    cpdef double downsizingFactorMax_2(self):
        cdef double A = (self.Cl/self.Cd())*self.originalAircraft.powerShaft()/self.VTAS*self.efficiencyPropulsive*(1+self.efficiencyMotor*self.efficiencyEGHS*self.efficiencyExhaust/self.efficiencyICE)
        return (A-(self.originalAircraft.weightTotalSoc()+self.weightEGHS()+self.weightMotor()))/(A-(self.originalAircraft.weightICE()+self.originalAircraft.weightFuelSoc()))

    cpdef double rangeX_2(self, double massFuelBurnt):
        cdef double massTotalSoc = self.massTotalSoc()
        return (self.Cl/self.Cd())*self.efficiencyPropulsive/self.Cp()*log((massTotalSoc-massFuelBurnt)/massTotalSoc)
    
    cpdef double massFuelBurnt_2(self, double rangeX):
        return self.massTotalSoc()*(1-exp(-rangeX*self.Cp()/(self.Cl/self.Cd()*self.efficiencyPropulsive)))
    
    cpdef double volumeFuelBurnt_2(self, double rangeX):
        return self.massFuelBurnt_2(rangeX)/self.engine.fuel.density

    cpdef double massTotalEoc_2(self, double rangeX):
        cdef double dm = self.massFuelBurnt_2(rangeX)
        cdef double m1 = self.massTotalSoc()
        return m1 - dm

    #-------------------------
    # Flight Mode 1
    #-------------------------
    cpdef double downsizingFactorMax_1(self):
        return 1 - 0.75/(1+self.efficiencyExhaust/self.efficiencyICE*self.efficiencyHarvester*self.efficiencyGenerator*self.efficiencyRectifier*self.efficiencyInverter*self.efficiencyMotor)
    """
        cdef double qS = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area
        cdef double W1 = self.originalAircraft.weightTotalSoc() + self.weightEGHS() + self.weightMotor()
        cdef double W2 = self.originalAircraft.weightICE() + self.originalAircraft.weightFuelSoc()
        cdef double P = self.originalAircraft.powerShaft()/qS/self.Cd0/self.VTAS*self.efficiencyPropulsive*(1+self.efficiencyMotor*self.efficiencyEGHS*self.efficiencyExhaust/self.efficiencyICE)
        cdef double B = 2*W1/W2-P*(qS*sqrt(self.Cd0/self.wing.k())/W2)**2
        cdef double C = (W1/W2)**2+(1-P)*(qS*sqrt(self.Cd0/self.wing.k())/W2)**2
        #
        cdef double D = B**2-4*C
        if D < 0:
            return nan
        else:
            r1 = B/2 + sqrt(D)/2
            #r2 = B/2 - sqrt(D)/2
            #print("+ve root={}, -ve root={}".format(r1, r2))
            return r1
    """
        
    cpdef double rangeX_1(self, double massFuelBurnt):
        cdef double a = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area*sqrt(self.Cd0/self.wing.k())
        cdef double Wsoc = self.weightTotalSoc() # start of cruise
        cdef double Wfb = massFuelBurnt*GRAVITY # end of cruise
        cdef double numerator = a * Wfb
        cdef double denominator = a**2 + Wsoc*(Wsoc-Wfb)
        return self.efficiencyPropulsive/self.Cp()/sqrt(self.Cd0*self.wing.k())*atan(numerator/denominator)
    
    cpdef double massFuelBurnt_1(self, double rangeX):
        cdef double a = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area*sqrt(self.Cd0/self.wing.k())
        cdef double Wsoc = self.massTotalSoc()*GRAVITY # start of cruise
        cdef double numerator = a**2 + Wsoc**2
        cdef double tan_ = tan(rangeX*self.Cp()*sqrt(self.Cd0*self.wing.k())/self.efficiencyPropulsive)
        cdef double denominator = a/tan_ + Wsoc
        return numerator/denominator/GRAVITY        
    
    cpdef double volumeFuelBurnt_1(self, double rangeX):
        return self.massFuelBurnt_1(rangeX)/self.engine.fuel.density

    cpdef double massTotalEoc_1(self, double rangeX):
        cdef double dm = self.massFuelBurnt_1(rangeX)
        cdef double m1 = self.massTotalSoc()
        return m1 - dm
    
    #-------------------------
    # Auto Mode
    #-------------------------
    cpdef double shaftToICERatio(self):
        return 1+self.efficiencyExhaust/self.efficiencyICE*self.efficiencyHarvester*self.efficiencyGenerator*self.efficiencyRectifier*self.efficiencyInverter*self.efficiencyMotor
    
    cpdef double downsizingFactorMax(self):
        if self.mode == 1:
            return self.downsizingFactorMax_1()
        elif self.mode == 2:
            return self.downsizingFactorMax_2()
        elif self.mode == 3:
            return self.downsizingFactorMax_3()
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))
    
    cpdef double rangeX(self, double massFuelBurnt):
        if self.mode == 1:
            return self.rangeX_1(massFuelBurnt)
        elif self.mode == 2:
            return self.rangeX_2(massFuelBurnt)
        elif self.mode == 3:
            return self.rangeX_3(massFuelBurnt)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))
    
    cpdef double massFuelBurnt(self, double rangeX):
        if self.mode == 1:
            return self.massFuelBurnt_1(rangeX)
        elif self.mode == 2:
            return self.massFuelBurnt_2(rangeX)
        elif self.mode == 3:
            return self.massFuelBurnt_3(rangeX)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))
        
    
    cpdef double volumeFuelBurnt(self, double rangeX):
        if self.mode == 1:
            return self.volumeFuelBurnt_1(rangeX)
        elif self.mode == 2:
            return self.volumeFuelBurnt_2(rangeX)
        elif self.mode == 3:
            return self.volumeFuelBurnt_3(rangeX)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))

    cpdef double massTotalEoc(self, double rangeX):
        if self.mode == 1:
            return self.massTotalEoc_1(rangeX)
        elif self.mode == 2:
            return self.massTotalEoc_2(rangeX)
        elif self.mode == 3:
            return self.massTotalEoc_3(rangeX)
        else:
            raise ValueError("mode must be 1,2 or 3. Given: {}".format(self.mode))


    #-------------------------
    # Series topology
    #-------------------------
    cpdef double shaftToICERatioSeries(self):
        return self.efficiencyInverter*self.efficiencyMotor*(self.efficiencyGenerator*self.efficiencyRectifier+self.efficiencyExhaust/self.efficiencyICE*self.efficiencyHarvester*self.efficiencyGenerator*self.efficiencyRectifier)

    cpdef double downsizingFactorMaxSeries(self):
        return 1 - 0.75/self.shaftToICERatioSeries()
        
    cpdef double rangeXSeries(self, double massFuelBurnt):
        cdef double a = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area*sqrt(self.Cd0/self.wing.k())
        cdef double Wsoc = self.weightTotalSocSeries() # start of cruise
        cdef double Wfb = massFuelBurnt*GRAVITY # end of cruise
        cdef double numerator = a * Wfb
        cdef double denominator = a**2 + Wsoc*(Wsoc-Wfb)
        return self.efficiencyPropulsive/self.Cp()/sqrt(self.Cd0*self.wing.k())*atan(numerator/denominator)
    
    cpdef double massFuelBurntSeries(self, double rangeX):
        cdef double a = 0.5*self.rhoSoc*self.VTAS**2*self.wing.area*sqrt(self.Cd0/self.wing.k())
        cdef double Wsoc = self.massTotalSocSeries()*GRAVITY # start of cruise
        cdef double numerator = a**2 + Wsoc**2
        cdef double tan_ = tan(rangeX*self.Cp()*sqrt(self.Cd0*self.wing.k())/self.efficiencyPropulsive)
        cdef double denominator = a/tan_ + Wsoc
        return numerator/denominator/GRAVITY        
    
    cpdef double volumeFuelBurntSeries(self, double rangeX):
        return self.massFuelBurntSeries(rangeX)/self.engine.fuel.density

    cpdef double massTotalEocSeries(self, double rangeX):
        cdef double dm = self.massFuelBurntSeries(rangeX)
        cdef double m1 = self.massTotalSocSeries()
        return m1 - dm
    
    #-------------------------
    # Python stuff
    #-------------------------
        
    def __str__(self):
        return self.model

    def summary(self, engine=True, wing=True, title="Hybrid Aircraft summary", original=False, original_kwargs={'engine': True, 'wing': True}):
        print(title)
        table = [["model", self.model, ""],
                 ["downsizingFactor", self.downsizingFactor, ""],
                 ["efficiencyHarvester", self.efficiencyHarvester, ""],
                 ["powerHarvester", self.powerHarvester()/1000, "kW"],
                 ["specificPowerHarvester", self.specificPowerHarvester()/1000, "kW/kg"],
                 ["massHarvester", self.massHarvester, "kg"],
                 ["massEmpty", self.massEmpty(), "kg"],
                 ["massTotalSoc", self.massTotalSoc(), "kg"],
                 ["massICE", self.massICE(), "kg"],
                 ["fuelVolume", self.fuelVolume, "m^3"],
                 ["", conv.m32L(self.fuelVolume), "L"],
                 ["massFuelSoc", self.massFuelSoc(), "kg"],
                 ["Cl", self.Cl, ""],
                 ["Cd", self.Cd(), ""],
                 ["Cl/Cd", self.Cl/self.Cd(), ""],
                 ["Cd0", self.Cd0, ""],
                 ["efficiencyPropulsive", self.efficiencyPropulsive, ""],
                 ["Cp", self.Cp(), "N/s.W, N/J"],
                 ["density", self.rhoSoc, "kg/m^3"],
                 ["cruise speed", self.VTAS, "m/s"],
                 ["", conv.mps2knot(self.VTAS), "knot"]]
        print(tabulate(table))
        if engine:
            self.engine.summary()
        if wing:
            self.wing.summary()
        if original:
            self.originalAircraft.summary(**original_kwargs, title="Original Aircraft summary")

    @property
    def mode(self):
        return self.originalAircraft.mode

    @mode.setter
    def mode(self, value):
        self.originalAircraft.mode = value

    @property
    def model(self):
        return "EGHS-Hybrid " + self.originalAircraft.model

    @property
    def modelId(self):
        return self.modelIdPrefix + self.originalAircraft.modelId

    @property
    def payload(self):
        return self.originalAircraft.payload

    @payload.setter
    def payload(self, value):
        self.originalAircraft.payload = value

    @property
    def Cl(self):
        return self.originalAircraft.Cl

    @Cl.setter
    def Cl(self, value):
        self.originalAircraft.Cl = value

    def Cd(self):
        return self.originalAircraft.Cd()

    @property
    def Cd0(self):
        return self.originalAircraft.Cd0

    @Cd0.setter
    def Cd0(self, value):
        self.originalAircraft.Cd0 = value

    @property
    def efficiencyPropulsive(self):
        return self.originalAircraft.efficiencyPropulsive

    @efficiencyPropulsive.setter
    def efficiencyPropulsive(self, value):
        self.originalAircraft.efficiencyPropulsive = value

    @property
    def rhoSoc(self):
        return self.originalAircraft.rhoSoc

    @rhoSoc.setter
    def rhoSoc(self, value):
        self.originalAircraft.rhoSoc = value

    def set_rhoSoc(self, double altitude):
        self.originalAircraft.rhoSoc = isa(altitude)['rho']
        
    @property
    def VTAS(self):
        return self.originalAircraft.VTAS

    @VTAS.setter
    def VTAS(self, value):
        self.originalAircraft.VTAS = value
        
    @property
    def VTASClimb(self):
        return self.originalAircraft.VTASClimb

    @VTASClimb.setter
    def VTASClimb(self, value):
        self.originalAircraft.VTASClimb = value
        
    def Vmp(self):
        return self.originalAircraft.Vmp()
        
    def Vmd(self):
        return self.originalAircraft.Vmd()
        
    def Vcc(self):
        return self.originalAircraft.Vcc()

    @property
    def wing(self):
        return self.originalAircraft.wing

    @wing.setter
    def wing(self, value):
        self.originalAircraft.wing = value
        
