import matplotlib.pyplot as plt
import math
import numpy as np
"""
make a class structure for defining each function and its deriv and eval...
"""
class Surface:
    def __init__(self):
        pass
    
    def func_prime_eval(self,coords):
        dif = 0.00001
        x1 = [coords[0]-dif, coords[0], coords[0]+dif]
        x2 = [coords[1]-dif, coords[1], coords[1]+dif]
        
        return -1*np.asarray([ (self.func_eval([x1[2],x2[1]]) - self.func_eval([x1[0],x2[1]])) / (2*dif), (self.func_eval([x1[1],x2[2]]) - self.func_eval([x1[1],x2[0]])) / (2*dif)])
        
    def norm_func_prime_eval(self,coords):
        eval = self.func_prime_eval(coords)
        return eval/np.linalg.norm(eval)
        
    def checkBounds(self,initialCoords):
        for i,_ in enumerate(initialCoords):
            if self.boundary[i][0] > initialCoords[i] or self.boundary[i][1] < initialCoords[i]:
                return 1
        return 0

class Styblinski_Tang(Surface):

    def __init__(self):
        self.boundary = [[-5,5],[-5,5]]
        
    def func_eval(self,coords):
        x1 = coords[0]
        x2 = coords[1]
        return 0.5 * (( x1**4 - 16*x1**2 + 5*x1 ) + ( x2**4 - 16*x2**2 + 5*x2 ) )

    def func_prime_eval(self,coords):
        dif = 0.1
        x1 = [coords[0]-dif, coords[0], coords[0]+dif]
        x2 = [coords[1]-dif, coords[1], coords[1]+dif]
        
        return -1*np.asarray([ (self.func_eval([x1[2],x2[1]]) - self.func_eval([x1[0],x2[1]])) / 2*dif, (self.func_eval([x1[1],x2[2]]) - self.func_eval([x1[1],x2[0]])) / 2*dif])

    def initialPlot(self):
    
        xlist = np.linspace(self.boundary[0][0], self.boundary[0][1], 500)
        ylist = np.linspace(self.boundary[1][0], self.boundary[1][1], 500)

        X, Y = np.meshgrid(xlist, ylist)
        Z = self.func_eval([X,Y])
        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(X, Y, Z,levels=50,vmin=-100,vmax=100)
        fig.colorbar(cp)
    
        return ax
        
class Muller_Brown(Surface):

    def __init__(self):
        self.boundary = [[-2,1],[-0.5,2]]
        
    def func_eval(self,coords):
        x1 = coords[0]
        x2 = coords[1]
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
        
    def initialPlot(self):
    
        xlist = np.linspace(self.boundary[0][0], self.boundary[0][1], 400)
        ylist = np.linspace(self.boundary[1][0], self.boundary[1][1], 400)

        X, Y = np.meshgrid(xlist, ylist)
        Z = self.func_eval([X,Y])
        fig,ax=plt.subplots(1,1)
        cp = ax.contourf(X, Y, Z,levels=100,vmin=-100,vmax=100)
        fig.colorbar(cp)
    
        return ax
        
if __name__ == '__main__':
    pass
    


    
