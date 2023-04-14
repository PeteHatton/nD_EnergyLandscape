import numpy as np
import matplotlib.pyplot as plt
import Landscapes as ls
from Utilities import *

def initialPlot():
    
    xlist = np.linspace(-5.0, 5.0, 500)
    ylist = np.linspace(-5.0, 5.0, 500)

    X, Y = np.meshgrid(xlist, ylist)
    Z = ls.Styblinski_Tang().func_eval([X,Y])
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, Z)
    fig.colorbar(cp)
    
    return ax
    
class Minimizer:
    
    def __init__(self,stepSize):
        pass
        
    def step(self,coords,direction):
        pass

class Steepest_Descent_fixed_step(Minimizer):
    
    def __init__(self,stepSize):
        self.stepSize = stepSize
        pass
        
    def step(self,coords,direction):
        return coords + self.stepSize*direction
        




    
if __name__=='__main__':
    
    ax = initialPlot()
    min = Steepest_Descent_fixed_step(0.01)
    maxIter = 1000
    
    coords = [1,0] # initial coords
    ax.scatter(coords[0],coords[1], color = 'r' , marker='s',s=50) # initial coord plot
    
    #initial
    F = ls.Styblinski_Tang().func_prime_eval(coords)
    dir = ls.Styblinski_Tang().norm_func_prime_eval(coords)
    coords = min.step(coords,dir) #initial step
    
    i=0
    while np.max(np.abs(F)) > 0.1 and i < maxIter:
        F = ls.Styblinski_Tang().func_prime_eval(coords)
        dir = ls.Styblinski_Tang().norm_func_prime_eval(coords)
        coords = min.step(coords,dir)
        ax.scatter(coords[0],coords[1], color = 'r' ,alpha=0.2)
        i=i+1
    print(i)
    ax.scatter(coords[0],coords[1], color = 'r' , marker='*',s=50)
    log('test','test')
    plt.show()

