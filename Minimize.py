import numpy as np
import matplotlib.pyplot as plt
import copy
import Landscapes as ls
from Utilities import *
import Lattice as lt
    
class Minimizer:
    
    def __init__(self):
        pass
        
    def run(self,obj):
        self.runMain(obj)
        pass
        
    def runMain(self):
        pass
        
class Steepest_Descent_fixed_step(Minimizer):
    
    def __init__(self):
        self.stepSize = 0.1
        self.maxIter = 1000
        self.forceTol = 0.001
        self.iterations = 0
        
    def runMain(self,obj):
        obj.axis.scatter(obj.coords[0],
                        obj.coords[1],
                        color = 'r' ,
                        marker='s',
                        s=50) # initial coord plot
                        
        obj.energy = obj.surf.func_eval(obj.coords)
        obj.force = obj.surf.func_prime_eval(obj.coords)
        dir = obj.surf.norm_func_prime_eval(obj.coords)
        
        obj.coords = obj.coords + self.stepSize*dir
        
        i = self.iterations
        
        while np.max(np.abs(obj.force)) > self.forceTol and i < self.maxIter:

            obj.energy = obj.surf.func_eval(obj.coords)
            obj.force = obj.surf.func_prime_eval(obj.coords)
            
            dir = obj.surf.norm_func_prime_eval(obj.coords)
            
            obj.coords += self.stepSize*dir
            if obj.surf.checkBounds(obj.coords):
                print('outside of range')
                return 0
            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
            
            i=i+1
            
        self.iterations = i
        
        obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' , marker='*',s=50)
        return 0

        
class Steepest_Descent_adaptive_step(Minimizer):
    
    def __init__(self):
        self.stepSize = 0.1
        self.maxIter = 1000
        self.forceTol = 0.001
        self.iterations = 0
        
    def runMain(self,obj):
        obj.axis.scatter(obj.coords[0],
                        obj.coords[1],
                        color = 'r' ,
                        marker='s',
                        s=50) # initial coord plot
                        
        obj.energy = obj.surf.func_eval(obj.coords)
        obj.force = obj.surf.func_prime_eval(obj.coords)
        dir = obj.surf.norm_func_prime_eval(obj.coords)
        
        obj.coords = obj.coords + self.stepSize*dir
        
        i = self.iterations
        
        while np.max(np.abs(obj.force)) > self.forceTol and i < self.maxIter:

            obj.energy = obj.surf.func_eval(obj.coords)
            obj.force = obj.surf.func_prime_eval(obj.coords)
            
            dir_old = copy.deepcopy(dir)
            dir = obj.surf.norm_func_prime_eval(obj.coords)
            
            if np.dot(dir,dir_old) > 0:
                self.stepSize *= 1.2
            else:
                self.stepSize *= 0.5
            
            obj.coords += self.stepSize*dir
            if obj.surf.checkBounds(obj.coords):
                print('outside of range')
                return 0
            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
            
            i=i+1
            
        self.iterations = i
        
        obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' , marker='*',s=50)
        return 0

def main():
    
    lattice = lt.lattice()
    lattice.pullConfig()
    
    min = Steepest_Descent_adaptive_step()
    lattice.axis = lattice.surf.initialPlot()
    
    _ = min.run(lattice)
    
    print(min.iterations)
    log(__name__,'test')
    plt.show()
    
if __name__=='__main__':
    pass

