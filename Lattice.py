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
    
    def pullConfig(self):
        data = ut.ORlC('config.cfg',1)
        self.coords = [ float(data[0].split()[i]) for i in range(2) ]
        self.surf = ls.Muller_Brown()
        #data[1]
    
if __name__=='__main__':
    pass
    
