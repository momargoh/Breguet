from .misc import GRAVITY, isa
from . import conversions as conv
from math import pi, tan, atan, sqrt, nan, isnan, exp, log
from tabulate import tabulate

cdef class Wing():
    def __init__(self, double span, double area, double _oswaldEfficiency=0):
        self.span = span
        self.area = area
        self._oswaldEfficiency = _oswaldEfficiency
    
    cpdef public double chord(self):
        return self.area / self.span
    
    cpdef public double aspectRatio(self):
        return self.span**2 / self.area
    
    cpdef public double oswaldEfficiency(self):
        "float: Equation from Obert 2009 Aerodynamic Design of Transport Aircraft, Delft: IOS Press"
        if self._oswaldEfficiency == 0:
            return 1/(1.05+0.007*pi*self.aspectRatio())
        else:
            return self._oswaldEfficiency

    cpdef public double k(self):
        return 1./(pi * self.aspectRatio() * self.oswaldEfficiency())
    
    cpdef Wing copy(self):
        return Wing(self.span, self.area, self._oswaldEfficiency)
    
    def summary(self, title="Wing summary"):
        print(title)
        table = [["span", self.span, "m"],
                 ["area", self.area, "m^2"],
                 ["chord", self.chord(), "m"],
                 ["aspectRatio", self.aspectRatio(), ""],
                 ["oswaldEfficiency", self.oswaldEfficiency(), ""],
                 ["k", self.k(), ""]]
        print(tabulate(table))
