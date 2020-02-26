from math import exp, nan
#--------------------------------------
# Atmosphere
#--------------------------------------
GRAVITY = 9.80665 # [m/s^2]
R = 287.058 #[J/kg.K]

cdef (double, double) _isaPropsFromBase(double altitude, double lapseRate, double altitudeBase, double pStagBase, double TStagBase):
    cdef double T=nan
    cdef double p=nan
    if lapseRate != 0:
        T = TStagBase + lapseRate * (altitude - altitudeBase)
        p = pStagBase * (T / TStagBase)**(-GRAVITY / lapseRate / R)
    else:
        T = TStagBase
        p = pStagBase * exp(-GRAVITY / R / TStagBase * (altitude - altitudeBase))
    return (p, T)

cpdef dict isa(double altitude):
    """FlowState: Returns stagnation density, pressure & temperature at desired altitude in the International Standard Atmosphere.


Parameters
-----------
altitude : float
    Geopotential altitude [m].
"""
    assert altitude >= -610, "Altitude must be above -610 [m] (given: {})".format(
        altitude)
    assert altitude <= 86000, "Altitude must be below 86,000 [m] (given: {})".format(
        altitude)
    cdef double[7] refLapseRate = [-0.0065, 0, 0.001, 0.0028, 0, -0.0028, -0.002]
    cdef double[8] refAltitude = [0, 11000, 20000, 32000, 47000, 51000, 71000, 86000]
    cdef double _pStag = 101325.
    cdef double _TStag = 288.15
    cdef size_t i = 0
    while i < len(refAltitude) - 1:
        if altitude <= refAltitude[i + 1]:
            _pStag, _TStag = _isaPropsFromBase(altitude, refLapseRate[i], refAltitude[i], _pStag, _TStag)
            break
        else:
            _pStag, _TStag = _isaPropsFromBase(refAltitude[i + 1], refLapseRate[i], refAltitude[i], _pStag, _TStag)
            i += 1

    cdef rho = _pStag/(_TStag*R)
    return {'rho': rho, 'p': _pStag, 'T': _TStag}

