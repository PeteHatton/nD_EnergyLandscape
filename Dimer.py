import numpy as np
import matplotlib.pyplot as plt
import copy
import sys

import Landscapes as ls
import Utilities as ut
import Lattice as lt
        
class Dimer:
    
    def __init__(self):
    
        self.minEnergy = None
        self.dimerOffset = 0.1
        self.dimerStepSize = 0.001
        self.direction = None # length 1 direction, initially random.
        self.dimerCoords_1 = None
        self.dimerCoords_2 = None
        self.dimerForce_1 = None
        self.dimerForce_2 = None
        
        
    def run(self,obj):
        
        #setup and plot
        ut.log(__name__ , 'Initializing Dimer',1)
        self.initializeDimer(obj)
        obj.axis.scatter(obj.coords[0],obj.coords[1])
        
        #some loop over this
        ut.log(__name__ , 'Walking the dimer',1)
        for i in range(100):
            self.minDimer_translate(obj)
            self.minDimer_rot(obj)
            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
            ut.log(__name__ , 'Energy: '
                            + str(round(obj.energy,5))
                            + ', Rel. En.: '
                            + str(round(obj.energy - self.minEnergy,5))
                            ,2)
        
        ut.log(__name__ , 'Dimer walker finished',1)
        
    def initializeDimer(self,obj):
        
        #initialize with random direction
        ut.log(__name__ , 'Choosing random inital direction',1)
        self.dirRand()
        
        #create the dimer structure
        ut.log(__name__ , 'Setting up dimer structure',1)
        self.dimerCoords_1 = self.dimerOffset * self.direction / 2
        self.dimerCoords_2 = - self.dimerOffset * self.direction / 2
        
        #calculate initial energy
        obj.energy = obj.surf.func_eval(obj.coords)
        self.minEnergy = copy.deepcopy(obj.energy)
        
    def dirRand(self):
        vecRand = np.random.rand(1,2)[0] # i think this only return positive numbers so need to think about how to do this.
        self.direction = vecRand / np.linalg.norm(vecRand)
        
    def minDimer_translate(self,obj):
        
        #step all components in the new direction with some step size
        self.dimerCoords_1 += self.dimerStepSize * self.direction
        self.dimerCoords_2 += self.dimerStepSize * self.direction
        obj.coords += self.dimerStepSize * self.direction
        
        #calc energy of new point
        obj.energy = obj.surf.func_eval(obj.coords)

    def minDimer_rot(self,obj):
        pass

def main():
    ut.log(__name__ , 'Begin',0)
    
    lattice = lt.lattice()
    lattice.pullConfig()
    
    lattice.axis = ls.surfPlot(lattice)
    dimer = Dimer()
    
    dimer.run(lattice)

    # show the generated plot
    ut.log(__name__ , 'Plotting...',0)
    plt.show()
    
    ut.log(__name__ , 'Fin.',0)
    
if __name__ == '__main__':
    pass

