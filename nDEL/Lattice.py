import matplotlib.pyplot as plt
import math
import numpy as np
import copy

import nDEL.Landscapes as ls
import nDEL.Minimize as mn
import nDEL.Utilities as ut

class lattice:

    def __init__(self,params):
        self.params = copy.deepcopy(params)
        self.coords = copy.deepcopy(params.initialCoords)
        
        #plotting axis initialization
        self.axis= None
        
        #energy and force initialization
        self.energy = 0.0
        self.force = [np.inf,np.inf]
        self.normF = None

        #
        self.dimerCount = 0
        
        self.surf = getSurface(params.Surface,self.params)

        self.saddlePoints = []

    
def getSurface(dat,params):

    if dat == 'Muller Brown':
        surface = ls.Muller_Brown(params)

    elif dat == 'Styblinski Tang':
        surface = ls.Styblinski_Tang(params)

    elif dat == 'Egg Holder':
        surface = ls.Egg_Holder(params)
        
    elif dat == 'Schwefel':
        surface = ls.Schwefel(params)

    return surface



    
if __name__=='__main__':
    pass
    
