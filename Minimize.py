import numpy as np
import matplotlib.pyplot as plt
import copy
import sys

import Landscapes as ls
import Utilities as ut
import Lattice as lt
import Input as ip

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
        
        # plot initial coordinate
        obj.axis.scatter(obj.coords[0],
                        obj.coords[1],
                        color = 'r' ,
                        marker='s',
                        s=50)
                        
        #initialize step counter
        i = 0
        
        #Energy and force calculations
        obj.energy = obj.surf.func_eval(obj.coords)
        obj.force = obj.surf.func_prime_eval(obj.coords)
        obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
        
        ut.log(__name__ , 'STEP: '+str(i)+' E = ' +str(round(obj.energy,3))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),3)),2)
        
        #Make steepest descent step
        obj.coords = obj.coords + self.stepSize*obj.normF
        
        while np.max(np.abs(obj.force)) > self.minForceTol and i < self.MaxIt:
            
            #increase step counter
            i+=1
            
            #Energy and force calculations
            obj.energy = obj.surf.func_eval(obj.coords)
            obj.force = obj.surf.func_prime_eval(obj.coords)
            normF_old = copy.deepcopy(obj.normF)
            obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
            
            ut.log(__name__ , 'STEP: '+str(i)+'. E = ' +str(round(obj.energy,3))+ ', max(F) = '+str(round(np.max(np.abs(obj.force)),5)),2)
            
            #Change the step size
            if np.dot(obj.normF,normF_old) > 0:
                self.stepSize *= 1.2
            else:
                self.stepSize *= 0.4
            
            #Make steepest descent step
            obj.coords += self.stepSize*obj.normF
            
            #plot current position
            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
        
        # record number of steps we took
        obj.minIterations = i
        
        #plot minimized position
        obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' , marker='*',s=50)
        
        ut.log(__name__ , 'Minimization Complete, Final Coordinates: [' + str(obj.coords[0])+','+str(obj.coords[1]) + '], Energy: ' +str(round(obj.energy,5))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),5)),1)
        
        return 0

def main():
    
    #read Params
    minParams = ip.getParams()
    
    #Set-up configuration
    lattice = lt.lattice(minParams)
    
    ut.log(__name__, 'Initializing Minimizer')
    min = getMinimizer(minParams)
    lattice.axis = ls.surfPlot(lattice)
    
    ut.log(__name__, 'Running '+ minParams.minAlgorithm +' Minimization',1)
    if min.run(lattice):
        sys.exit()
    
    ut.log(__name__ , 'Plotting...',0)
        
    # show the generated plot
    plt.show()
    
    ut.log(__name__ , 'Fin.',0)
    
if __name__ == '__main__':
    pass

