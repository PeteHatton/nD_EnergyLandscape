import matplotlib.pyplot as plt
import math
import numpy as np
import copy

import Landscapes as ls
import Minimize as mn
import Utilities as ut

class lattice:

    def __init__(self,params):
    
        #Surface initialization
        self.coords = copy.deepcopy(params.initialCoords)
        
        #plotting axis initialization
        self.axis= None
        
        #energy and force initialization
        self.energy = 0.0
        self.force = [np.inf,np.inf]
        self.normF = None
        
        self.surf = getSurface(params.Surface)

        self.saddlePoints = []

    
def getSurface(dat):

    if dat == 'Muller Brown':
        surface = ls.Muller_Brown()

    elif dat == 'Styblinski Tang':
        surface = ls.Styblinski_Tang()

    elif dat == 'Egg Holder':
        surface = ls.Egg_Holder()
        
    elif dat == 'Schwefel':
        surface = ls.Schwefel()

    return surface
    
if __name__=='__main__':
    pass
    
