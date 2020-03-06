from .fuel import Fuel
from .ice import ICE
from .wing import Wing
from .aircraft import Aircraft
from .conversions import galUS2m3, lb2kg, bhp2W, hour2s, knot2mps
from math import nan
"""
Approximating AvGas 100/100LL density as 718 kg/m^3.
"""

#--------------------------------------
# Fuel database
#--------------------------------------
avgas100LL = Fuel(density=718, name='Avgas 100LL')


#--------------------------------------
# ICEngine database
#--------------------------------------
def lycomingIO720A1B(fuel=avgas100LL):
    fcv = galUS2m3(33.9) / hour2s(1)
    return ICE(
        massDry=lb2kg(601),
        powerRated=bhp2W(400),
        fuelConsumptionVolumetric=fcv,
        fuel=fuel)


def lycomingIO720D1CD(fuel=avgas100LL):
    fcv = galUS2m3(33.9) / hour2s(1) * (375 / 400)
    return ICE(
        massDry=lb2kg(607),
        powerRated=bhp2W(375),
        fuelConsumptionVolumetric=fcv,
        fuel=fuel)


def continentalIOF240(fuel=avgas100LL):
    "Based on max. recommended cruise power of 94bhp. Maximum rated power is 125bhp."
    fcv = lb2kg(34.5) / hour2s(1) / 718
    return ICE(
        massDry=lb2kg(255),
        powerRated=bhp2W(94),
        fuelConsumptionVolumetric=fcv,
        fuel=fuel)


def continentalTSIO520C(fuel=avgas100LL):
    fcv = galUS2m3(28.15) / hour2s(1)
    return ICE(
        massDry=lb2kg(454),
        powerRated=bhp2W(285),
        fuelConsumptionVolumetric=fcv,
        fuel=fuel)


#--------------------------------------
# Aircraft database
#--------------------------------------


def piperPawneeBrave(mode=1):
    engine = lycomingIO720D1CD()
    wing = Wing(span=11.82, area=21.0)
    return Aircraft(
        engine=engine,
        wing=wing,
        massEmpty=1162.,
        massAlternator=5,
        massStarterMotor=4,
        massBatteries=10.3,
        payload=150.,
        fuelVolume=galUS2m3(36),
        Cd0=0.075,
        Cl=nan,
        efficiencyPropulsive=0.85,
        rhoSoc=nan,
        VTAS=knot2mps(113),
        VTASClimb=2.794,  #550 ft/min
        cruiseSpeed=knot2mps(113),
        rangeMax=772e3,
        model="Piper PA-36-375, Pawnee Brave",
        modelId="PA-36-375",
        mode=mode)


def libertyXL2(mode=1):
    "Based on 55% power cruise"
    engine = continentalIOF240()
    wing = Wing(span=8.53, area=10.4)
    return Aircraft(
        engine=engine,
        wing=wing,
        massEmpty=lb2kg(1050),
        massAlternator=5,
        massStarterMotor=4,
        massBatteries=10.7,
        payload=150.,
        fuelVolume=galUS2m3(28),
        Cd0=0.032,
        Cl=nan,
        efficiencyPropulsive=0.85,
        rhoSoc=nan,
        VTAS=knot2mps(120),
        VTASClimb=3.46456,  #682 ft/min
        cruiseSpeed=knot2mps(120),
        rangeMax=926e3,
        model="Liberty XL-2",
        modelId="XL-2",
        mode=mode)


def cessnaTU206F(mode=1):
    "Based on 55% power cruise"
    engine = continentalTSIO520C()
    wing = Wing(span=10.92, area=16.18)
    return Aircraft(
        engine=engine,
        wing=wing,
        massEmpty=lb2kg(1050),
        massAlternator=5,
        massStarterMotor=4,
        massBatteries=13.4,
        payload=150.,
        fuelVolume=galUS2m3(65),
        Cd0=0.035,
        Cl=nan,
        efficiencyPropulsive=0.85,
        rhoSoc=nan,
        VTAS=knot2mps(148),
        VTASClimb=5.2324,  #1030 ft/min
        cruiseSpeed=knot2mps(148),
        rangeMax=1327e3,
        model="Cessna TU206F, Super Skywagon",
        modelId="TU206F",
        mode=mode)
