import matplotlib.pyplot as plt
import math
import numpy as np
import copy

import Landscapes as ls
import Minimize as mn
import Utilities as ut

class lattice:

    def __init__(self,params):
        
        #plotting axis initialization
        self.axis= None
        
        #Surface initialization
        self.coords = copy.deepcopy(params.initialCoords)
        
        #energy and force initialization
        self.energy = 0
        self.force = [np.inf,np.inf]
        self.normF = None
        
        #Dimension
        self.dimension = 0
        
        #Minimization Settings
        self.minAlgorithm = ''
        self.maxIter = 0
        self.forceTol = 0
        self.minIterations = 0
        
        self.surf = getSurface(params.Surface)
        self.minAlgorithm = params.minAlgorithm
#        self.
        self.saddlePoints = []
        
    
#    def pullConfig(self):
#
#        #Read data from config file
#        data = ut.ORlC('config.cfg',1)
#
#        #parse the config lines
#        self.coords = [ float(data[0].split()[i]) for i in range(2) ]
#        self.surf = getSurface(data[1].strip())
#        self.dimenstion = int(data[2])
#        self.minAlgorithm = data[3].strip()
#        self.forceTol = float(data[4])
#        self.maxIter = float(data[5])
#        self.stepSize = float(data[6])
    
    
def getSurface(dat):

    if dat == 'Muller Brown':
        surface = ls.Muller_Brown()

    elif dat == 'Styblinski Tang':
        surface = ls.Styblinski_Tang()

    elif dat == 'Egg Holder':
        surface = ls.Egg_Holder()

    return surface
    
if __name__=='__main__':
    pass
    
