import matplotlib.pyplot as plt
import math
import numpy as np
import Landscapes as ls
import Minimize as mn
import Utilities as ut

class lattice:
    def __init__(self):
        self.axis= None
        self.surf = None
        self.coords = None
        self.energy = 0
        self.force = None
        self.surf = None
        self.normF = None
        
        #Minimization Settings
        self.stepSize = 0.1
        self.maxIter = 1000
        self.forceTol = 0.001
        self.minIterations = 0
    
    def pullConfig(self):
        data = ut.ORlC('config.cfg',1)
        self.coords = [ float(data[0].split()[i]) for i in range(2) ]
        self.surf = self.getSurface(data[1].strip())
    
    def getSurface(self,dat):
        if dat == 'Muller Brown':
            surface = ls.Muller_Brown()
        elif dat == 'Styblinski Tang':
            surface = ls.Styblinski_Tang()
        elif dat == 'Egg Holder':
            surface = ls.Egg_Holder()
        return surface
    
if __name__=='__main__':
    pass
    
