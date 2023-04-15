import numpy as np
import matplotlib.pyplot as plt
import copy
import Landscapes as ls
from Utilities import *
    
class Minimizer:
    
    def __init__(self):
        pass
        
    def main(self):
        pass

class Steepest_Descent_fixed_step(Minimizer):
    
    def __init__(self,ax,Surf,coords):
        self.stepSize = 0.02
        self.ax = ax
        self.Surf = Surf
        self.coords = coords
        self.energy = self.Surf.func_eval(self.coords)
        self.force = 0
        self.maxIter = 1000
        self.forceTol = 0.1
        self.iterations = 0
        
    def main(self):
        self.ax.scatter(self.coords[0],
                        self.coords[1],
                        color = 'r' ,
                        marker='s',
                        s=50) # initial coord plot
        print(self.coords)
        self.energy = self.Surf.func_eval(self.coords)
        self.force = self.Surf.func_prime_eval(self.coords)
        
        dir = self.Surf.norm_func_prime_eval(self.coords)
        
        self.coords = self.coords + self.stepSize*dir
        
        i = self.iterations
        
        while np.max(np.abs(self.force)) > self.forceTol and i < self.maxIter:

            self.energy = self.Surf.func_eval(self.coords)
            self.force = self.Surf.func_prime_eval(self.coords)
            
            dir = self.Surf.norm_func_prime_eval(self.coords)
            
            self.coords += self.stepSize*dir
            
            if self.Surf.checkBounds(self.coords):
                print('outside of range')
                return 0
                
            self.ax.scatter(self.coords[0],self.coords[1], color = 'r' ,alpha=0.2)
            i=i+1
        self.iterations = i
        
        self.ax.scatter(self.coords[0],self.coords[1], color = 'r' , marker='*',s=50)

        
class Steepest_Descent_adaptive_step(Minimizer):
    
    def __init__(self,ax,Surf,coords):
        self.stepSize = 0.1
        self.ax = ax
        self.Surf = Surf
        self.coords = coords
        self.energy = 0
        self.force = 0
        self.maxIter = 1000
        self.forceTol = 0.001
        self.iterations = 0
        
    def main(self):
        self.ax.scatter(self.coords[0],
                        self.coords[1],
                        color = 'r' ,
                        marker='s',
                        s=50) # initial coord plot
                        
        self.energy = self.Surf.func_eval(self.coords)
        self.force = self.Surf.func_prime_eval(self.coords)
        dir = self.Surf.norm_func_prime_eval(self.coords)
        
        self.coords = self.coords + self.stepSize*dir
        
        i = self.iterations
        
        while np.max(np.abs(self.force)) > self.forceTol and i < self.maxIter:

            self.energy = self.Surf.func_eval(self.coords)
            self.force = self.Surf.func_prime_eval(self.coords)
            
            dir_old = copy.deepcopy(dir)
            dir = self.Surf.norm_func_prime_eval(self.coords)
            
            if np.dot(dir,dir_old) > 0:
                self.stepSize *= 1.2
            else:
                self.stepSize *= 0.5
            
            self.coords += self.stepSize*dir
            if self.Surf.checkBounds(self.coords):
                print('outside of range')
                return 0
            self.ax.scatter(self.coords[0],self.coords[1], color = 'r' ,alpha=0.2)
            
            i=i+1
            
        self.iterations = i
        
        self.ax.scatter(self.coords[0],self.coords[1], color = 'r' , marker='*',s=50)

def main():
    
#    enSurf = ls.Styblinski_Tang()
    enSurf = ls.Muller_Brown()
    ax = enSurf.initialPlot()
    initialCoords = [0,1.5]
    if enSurf.checkBounds(initialCoords):
        print('outside of range')
        return 0
    min = Steepest_Descent_adaptive_step(ax,enSurf,initialCoords)
    
    min.main()
    print(min.iterations)
    log(__name__,'test')
    plt.show()
    
if __name__=='__main__':
    main()

