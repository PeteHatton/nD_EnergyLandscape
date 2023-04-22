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
        self.energy = 0
        self.dimerOffset = 0.1
        self.dimerStepSize = 0.001
        self.direction = None # length 1 direction, initially random.
        self.dimerCoords_1 = None
        self.dimerCoords_2 = None
        self.dimerForce_1 = None
        self.dimerForce_2 = None
        self.transForce = None
        self.dimerEnergy = None
        self.dTheta = None
        self.dimerForce_perp = None
        self.Theta = None
    
    def run(self,obj):
        
        #setup and plot
        ut.log(__name__ , 'Initializing Dimer',1)
        self.initializeDimer(obj)
        obj.axis.scatter(obj.coords[0],obj.coords[1])
        
        #some loop over this
        ut.log(__name__ , 'Walking the dimer',1)
        for i in range(10):
            self.calcTransForce(obj)
            self.minDimer_translate(obj)
            self.calcRotForce(obj)
            self.minDimer_rot(obj)
            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
            ut.log(__name__ , 'Energy: '
                            + str(round(obj.energy,5))
                            + ', Rel. En.: '
                            + str(round(self.minEnergy - obj.energy,5))
                            ,2)
        
        ut.log(__name__ , 'Dimer walker finished',1)
        
    def initializeDimer(self,obj):
        
        #initialize with random direction
        ut.log(__name__ , 'Choosing random inital direction',1)
        self.dirRand()
        
        #create the dimer structure
        ut.log(__name__ , 'Setting up dimer structure',1)
        self.dimerCoords_1 = self.dimerOffset * self.direction
        self.dimerCoords_2 = - self.dimerOffset * self.direction
        
        #calculate initial energy on images
        self.dimerEnergy_1 = obj.surf.func_eval(self.dimerCoords_1)
        self.dimerEnergy_2 = obj.surf.func_eval(self.dimerCoords_2)
        
        #calculate initial force on images
        self.dimerForce_1 = obj.surf.func_prime_eval(self.dimerCoords_1)
        self.dimerForce_2 = obj.surf.func_prime_eval(self.dimerCoords_2)
        
        #infer energy of dimer
        self.dimerEnergy = obj.surf.func_eval(self.dimerCoords_1) + obj.surf.func_eval(self.dimerCoords_2)
        
        #infer energy of midpoint and save it
        obj.energy = (self.energy/2
                    + (self.dimerOffset/4) * np.dot(self.dimerForce_1
                    - self.dimerForce_2,self.direction ))
        self.minEnergy = copy.deepcopy(obj.energy)
        
        #infer force on midpoint
        obj.force = (self.dimerForce_1 + self.dimerForce_2)/2
        
        
    def dirRand(self):

        vecRand = np.random.uniform(low = -1.0, high = 1.0, size = (1,2))[0]
        self.direction = vecRand / np.linalg.norm(vecRand)
        
    def minDimer_translate(self,obj):
        
        #step all components in the new direction with some step size
        self.dimerCoords_1 -= self.dimerStepSize * self.transForce
        self.dimerCoords_2 -= self.dimerStepSize * self.transForce
        obj.coords -= self.dimerStepSize * self.transForce
        
        #calc energy of new point
        obj.energy = obj.surf.func_eval(obj.coords)

    def minDimer_rot(self,obj):
#        self.dimerCoords_1 += ( self.direction*np.cos(self.dTheta) + self.Theta * np.sin(self.dTheta)) * self.dimerOffset
        pass
    
    def calcTransForce(self,obj):
        
        obj.force = obj.surf.func_prime_eval(obj.coords)
        
        cos_angle = np.dot(obj.force,self.direction) / ( np.dot(obj.force,obj.force) * np.dot(self.direction,self.direction) )
        self.transForce = obj.force - 2 * obj.force * cos_angle
#        self.transForce = - self.direction #np.dot(obj.force,self.direction)
    
    def calcRotForce(self,obj):
        
        #calc force on dimer perpendicular to the direction of dimer.
        dimerForce_perp_1 = self.dimerForce_1 - (np.dot(self.dimerForce_1,self.direction)*self.direction)
        dimerForce_perp_2 = self.dimerForce_2 - (np.dot(self.dimerForce_2,self.direction)*self.direction)
        self.dimerForce_perp = dimerForce_perp_1 - dimerForce_perp_2
        
        self.Theta = self.dimerForce_perp/np.linalg.norm(self.dimerForce_perp)
    

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

