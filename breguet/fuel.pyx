from tabulate import tabulate
cdef class Fuel():

    def __init__(self, double density, str name='Fuel Instance', **kwargs):
        self.density = density
        self.name = name
        
    def __str__(self):
        return self.name
    
    cpdef Fuel copy(self):
        return Fuel(self.density, self.name)
    
    cpdef public double specificVolume(self):
        return 1/self.density

     
    def summary(self, title=""):
        if not title:
            title = self.name + ' summary'
        print(title)
        table = [["density", self.density, "kg/m^3"],
                 ["specificVolume", self.specificVolume(), "m^3/kg"]]
        ret = tabulate(table)
        print(ret)
        return ret


