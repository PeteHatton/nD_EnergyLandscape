import numpy as np
import matplotlib.pyplot as plt
import copy
import sys

import nDEL.Landscapes as ls
import nDEL.Utilities as ut
import nDEL.Lattice as lt
import nDEL.Input as ip

def getMinimizer(params):
    if params.minAlgorithm == 'Steepest Descent':
        return Steepest_Descent_adaptive_step(params)

class Minimizer:
    
    def __init__(self,params):
        self.params = copy.deepcopy(params)
        
    def run(self,obj):
        pass

        
class Steepest_Descent_adaptive_step(Minimizer):
    
    def __init__(self,params):
        super(Steepest_Descent_adaptive_step, self).__init__(params)
        
    def run(self,obj):
        
        #pararms
        self.stepSize = self.params.minStepSize
        self.minForceTol = self.params.minForceTol
        self.MaxIt = self.params.minMaxIterations
        

                        
        #initialize step counter
        iter = 0
        
        #Energy and force calculations
        obj.energy = obj.surf.func_eval(obj.coords)
        obj.force = obj.surf.func_prime_eval(obj.coords)
        obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
        
        ut.log(__name__ , 'STEP: '+str(iter)+' E = ' +str(round(obj.energy,3))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),3)),2)

        # plot initial coordinate
        if self.params.Dimension==2:

            obj.axis.scatter(obj.coords[0],obj.coords[1],obj.energy,
                            color = 'r' ,
                            marker='s',
                            s=50)
        
        #Make steepest descent step
        obj.coords = obj.coords + self.stepSize*obj.normF
        
        #Check if we've left the bounds of the surface.
        status = obj.surf.checkBounds(obj.coords)

        if status:
            ut.log(__name__, 'MINIMIZATION FAILED! OoB! Steps: '+ str(iter)
                            +  '. E = ' + str(round(obj.energy,5))
                            + '. Final coords: ' + str(obj.coords)
                            ,1)
            return 1
        
        while np.max(np.abs(obj.force)) > self.minForceTol and iter < self.MaxIt:
            
            #increase step counter
            iter+=1

            #plot current position
            if self.params.Dimension == 2 and self.params.plotSurface:
                obj.axis.scatter(obj.coords[0],obj.coords[1],obj.energy ,color = 'r' ,alpha=0.2)
            
            #Energy and force calculations
            obj.energy = obj.surf.func_eval(obj.coords)
            obj.force = obj.surf.func_prime_eval(obj.coords)
            normF_old = copy.deepcopy(obj.normF)
            obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
            
            ut.log(__name__ , 'STEP: '+str(iter)+'. E = ' +str(round(obj.energy,3))+ ', max(F) = '+str(round(np.max(np.abs(obj.force)),5)),2)
            
            #Change the step size
            if np.dot(obj.normF,normF_old) > 0:
                self.stepSize *= 1.2
            else:
                self.stepSize *= 0.4
            
            #Make steepest descent step
            obj.coords += self.stepSize*obj.normF
            
            #Check if we've left the bounds of the surface.
            status = obj.surf.checkBounds(obj.coords)
            if status:
                ut.log(__name__, 'FAILED! OoB! Steps: '+ str(iter)
                                +  '. E = ' + str(round(obj.energy,5))
                                + '. Final coords: ' + str(obj.coords)
                                ,1)
                return 1
            
        # record number of steps we took
        obj.minIterations = iter
        
        #plot minimized position
        if self.params.Dimension==2:
            obj.axis.scatter(obj.coords[0],obj.coords[1],obj.energy, color = 'r' , marker='*',s=50)
        
        ut.log(__name__ , 'Minimization Complete, Final Coordinates: '+ str(obj.coords) +', Energy: ' +str(round(obj.energy,5))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),5)),1)
        
        return 0

def main():
    
    #read Params
    minParams = ip.getParams()
    
    #Set-up configuration
    lattice = lt.lattice(minParams)
    
    ut.log(__name__, 'Initializing Minimizer')
    min = getMinimizer(minParams)
    
    if minParams.plotSurface and minParams.Dimension == 2:
        lattice.axis = lattice.surf.surfPlot(lattice)
    
    ut.log(__name__, 'Running '+ minParams.minAlgorithm +' Minimization',1)
    if min.run(lattice):
        sys.exit()
    
    ut.log(__name__ , 'Plotting...',0)
        
    # show the generated plot
    plt.show()
    
    ut.log(__name__ , 'Fin.',0)
    
if __name__ == '__main__':
    pass

