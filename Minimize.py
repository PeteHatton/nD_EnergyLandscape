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
        
    def main(self):
        pass

class Steepest_Descent_fixed_step(Minimizer):
    
    def __init__(self,ax,Surf,stepSize,coords):
        self.stepSize = stepSize
        self.ax = ax
        self.Surf = Surf
        self.coords = coords
        self.energy = self.Surf.func_eval(self.coords)
        
    def main(self):
        self.ax.scatter(self.coords[0],
                        self.coords [1],
                        color = 'r' ,
                        marker='s',
                        s=50) # initial coord plot
                        
        E = self.Surf.func_eval(self.coords)
        F = self.Surf.func_prime_eval(self.coords)
        dir = self.Surf.norm_func_prime_eval(self.coords)
        self.coords = self.coords + self.stepSize*dir
        
        i=0
        while np.max(np.abs(F)) > 0.1 and i < 1000:
            self.energy = self.Surf.func_eval(self.coords)
            F = self.Surf.func_prime_eval(self.coords)
            dir = self.Surf.norm_func_prime_eval(self.coords)
            self.coords += self.stepSize*dir
            self.ax.scatter(self.coords[0],self.coords[1], color = 'r' ,alpha=0.2)
            i=i+1
        
        self.ax.scatter(self.coords[0],self.coords[1], color = 'r' , marker='*',s=50)

def main():
    ax = initialPlot()
    enSurf = ls.Styblinski_Tang()
    
    initialCoords = [1,0]
    
    min = Steepest_Descent_fixed_step(ax,enSurf,0.02,initialCoords)
    min.main()
    
    log('test','test')
    plt.show()
    
if __name__=='__main__':
    main()

