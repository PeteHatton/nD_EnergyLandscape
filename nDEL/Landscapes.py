import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import copy

from mpl_toolkits import mplot3d
from scipy.misc import derivative

class Surface:

    '''
    Generic surface object
    
    TO-DO:
        - Add support for plotting the effect of deflation

    '''

    def __init__(self,params):
        self.params = copy.deepcopy(params)

    def checkDimCompat(self,dim,dimRange):
        
        if dim != dimRange and dimRange != None:
            print('ERROR: Requested dimension of '+str(dim[0])+' is not possible on the '+self.params.Surface+' surface. Please choose a dimension in the range: '+str(dimRange))
            sys.exit()
        
    def checkBounds(self,coords):
        """
            need to make this N dimensional
        """
        flag = 0
        for i,coord in enumerate(coords):
            if coord < self.boundary[i][0] or coord >self.boundary[i][1]:
                flag = 1
        if flag == 1:
                return 1
        return 0
    
    def func_prime_eval(self,coords):

        dif = 0.00001

        XH = [ [coords[j]+dif if i==j else coords[j] for j in range(self.params.Dimension) ] for i in range(self.params.Dimension) ]
        XL = [ [coords[j]-dif if i==j else coords[j] for j in range(self.params.Dimension) ] for i in range(self.params.Dimension) ]

        return -1*np.asarray([ (self.func_eval(XH[i]) - self.func_eval(XL[i])) / (2*dif) for i in range(self.params.Dimension) ])
    
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
        ax.contour3D(X, Y, Z,70,alpha=0.7)

        return ax

class Styblinski_Tang(Surface):

    def __init__(self,params):
        self.params = copy.deepcopy(params)
        self.implementedDims = None
        
        self.checkDimCompat([self.params.Dimension],self.implementedDims)

        '''
        Suggested values for the Styblinski Tang surface [https://www.sfu.ca/~ssurjano/stybtang.html].
        
        '''

        self.boundary = [ [-5,5] for _ in range(self.params.Dimension)]
        
        #plotting params
        self.vmax = 100
        self.vmin = -100
        self.levels = 50
        self.plotPoints = 500
        
    def func_eval(self,coords):

        res = 0

        for D in range(self.params.Dimension):
            x = coords[D]
            res += ( x**4 - 16*x**2 + 5*x )

        return 0.5 * res
        
class Schwefel(Surface):

    def __init__(self,params):
        
        self.params = copy.deepcopy(params)
        self.implementedDims = None

        self.checkDimCompat([self.params.Dimension],self.implementedDims)

        self.boundary = [[-100,100] for _ in range(self.params.Dimension)]
        
        #plotting params
        self.vmax = 10
        self.vmin = -10
        self.levels = 10
        self.plotPoints = 100
        
    def func_eval(self,coords):
        
        res = 0

        for D in range(self.params.Dimension):
            x = coords[D]
            res += x*np.sin(np.sqrt(np.abs(x)))

        return -418.9829*self.params.Dimension + res



class Muller_Brown(Surface):

    def __init__(self,params):
        
        self.params = copy.deepcopy(params)
        
        

        '''
        Suggested values for the Muller Brown surface [https://www.wolframcloud.com/objects/demonstrations/TrajectoriesOnTheMullerBrownPotentialEnergySurface-source.nb].
        
        '''
        self.implementedDims = [2]
        self.checkDimCompat([self.params.Dimension],self.implementedDims)
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
        self.params = copy.deepcopy(params)

        '''
        Suggested values for the Egg Holder surface [https://www.sfu.ca/~ssurjano/egg.html].
        
        '''
        self.implementedDims = [2]
        
        self.checkDimCompat([self.params.Dimension],self.implementedDims)

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
