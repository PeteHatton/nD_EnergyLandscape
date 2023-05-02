import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import copy
from mpl_toolkits import mplot3d


class Surface:

    '''
    Generic surface object
    
    TO-DO:
        - Make generic for different dimensions depending on surface...
        - option for fully 3d plotting or contour
        - Add support for plotting the effect of deflation

    '''

    def __init__(self,params):
        self.params = copy.deepcopy(params)
        
    def checkBounds(self,coords):
        for i,_ in enumerate(coords):
            if self.boundary[i][0] > coords[i] or self.boundary[i][1] < coords[i]:
                return 1
        return 0
    
    def func_prime_eval(self,coords):

        dif = 0.00001
        x1 = [coords[0]-dif, coords[0], coords[0]+dif]
        x2 = [coords[1]-dif, coords[1], coords[1]+dif]
        
        return -1*np.asarray([ (self.func_eval([x1[2],x2[1]]) - self.func_eval([x1[0],x2[1]])) / (2*dif), (self.func_eval([x1[1],x2[2]]) - self.func_eval([x1[1],x2[0]])) / (2*dif)])
        
    def norm_func_prime_eval(self,coords):
        eval = self.func_prime_eval(coords)
        return eval/np.linalg.norm(eval)
    
    def surfPlot(self,obj):

        """
        initial plot of the energy surface on the region.
        
        """
        ax = plt.axes(projection='3d')

        # Make data.
        xlist = np.linspace(obj.surf.boundary[0][0], obj.surf.boundary[0][1], obj.surf.plotPoints)
        ylist = np.linspace(obj.surf.boundary[1][0], obj.surf.boundary[1][1], obj.surf.plotPoints)

        X, Y = np.meshgrid(xlist, ylist)
        Z = obj.surf.func_eval([X,Y])

        # Plot the surface.
        ax.contour3D(X, Y, Z,70)

        return ax

class Styblinski_Tang(Surface):

    def __init__(self,params):
        self.params = copy.deepcopy(params)
        '''
        Suggested values for the Styblinski Tang surface [https://www.sfu.ca/~ssurjano/stybtang.html].
        
        '''
        self.boundary = [[-5,5],[-5,5]]
        self.vmax = 100
        self.vmin = -100
        self.levels = 50
        self.plotPoints = 500
        
    def func_eval(self,coords):
        [x1,x2] = coords
        
        return 0.5 * (( x1**4 - 16*x1**2 + 5*x1 ) + ( x2**4 - 16*x2**2 + 5*x2 ) )
        
class Schwefel(Surface):

    def __init__(self):
    
        self.boundary = [[-100,100],[-100,100]]
        self.vmax = 10
        self.vmin = -10
        self.levels = 10
        self.plotPoints = 100
        
    def func_eval(self,coords):
        [x1,x2] = coords
        return -418.9829*2 + x1*np.sin(np.sqrt(np.abs(x1))) + x2*np.sin(np.sqrt(np.abs(x2)))


class Muller_Brown(Surface):

    def __init__(self):
        '''
        Suggested values for the Muller Brown surface [https://www.wolframcloud.com/objects/demonstrations/TrajectoriesOnTheMullerBrownPotentialEnergySurface-source.nb].
        
        '''
        self.boundary = [[-2,1],[-0.5,2]]
        self.vmax = 100
        self.vmin = -100
        self.levels = 100
        self.plotPoints = 400
        
    def func_eval(self,coords):
        [x1,x2] = coords
        A = (-200,-100,-170,15)
        xo = (1,0,-0.5,-1)
        yo = (0,0.5,1.5,1)
        a = (-1,-1,-6.5,0.7)
        b = (0,0,11,0.6)
        c = (-10,-10,-6.5,0.7)
        
        res = 0
        for i in range(4):
            res += A[i] * np.exp( a[i]*(x1 - xo[i])**2 + b[i]*(x1 - xo[i])*(x2 - yo[i]) + c[i]*(x2 - yo[i])**2 )
        return res

class Egg_Holder(Surface):

    def __init__(self):
    
        '''
        Suggested values for the Egg Holder surface [https://www.sfu.ca/~ssurjano/egg.html].
        
        '''
        
        self.boundary = [[-512,512],[-512,512]]
        self.vmax = 1000
        self.vmin = -1000
        self.levels = 100
        self.plotPoints = 400
        
    def func_eval(self,coords):
        [x1,x2] = coords
        return -(x2+47)*np.sin(np.sqrt(np.abs(x2+x1/2 + 47)))-x1*np.sin(np.sqrt(np.abs(x1-(x2+47))))
        

    
if __name__ == '__main__':
    pass
