import numpy as np
import matplotlib.pyplot as plt
import copy
import Landscapes as ls
import Utilities as ut
import Lattice as lt
    
class Minimizer:
    
    def __init__(self):
        pass
        
    def run(self,obj):
        self.runMain(obj)
        pass
        
    def runMain(self):
        pass
        
#class Steepest_Descent_fixed_step(Minimizer):
#    
#    def __init__(self):
#        pass
#        
#    def runMain(self,obj):
#        obj.axis.scatter(obj.coords[0],
#                        obj.coords[1],
#                        color = 'r' ,
#                        marker='s',
#                        s=50) # initial coord plot
#                        
#        obj.energy = obj.surf.func_eval(obj.coords)
#        obj.force = obj.surf.func_prime_eval(obj.coords)
#        obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
#        
#        obj.coords = obj.coords + obj.stepSize*obj.normF
#        
#        i = 0
#        
#        while np.max(np.abs(obj.force)) > obj.forceTol and i < obj.maxIter:
#
#            obj.energy = obj.surf.func_eval(obj.coords)
#            obj.force = obj.surf.func_prime_eval(obj.coords)
#            
#            obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
#            
#            obj.coords += obj.stepSize*obj.normF
#            
#            if obj.surf.checkBounds(obj.coords):
#                print('outside of range')
#                return 0
#                
#            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
#            
#            i=i+1
#            
#        obj.minIterations = i
#        
#        obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' , marker='*',s=50)
#        
#        return 0

        
class Steepest_Descent_adaptive_step(Minimizer):
    
    def __init__(self):
        self.name = 'Steepest Descent'
        pass
        
    def runMain(self,obj):
        obj.axis.scatter(obj.coords[0],
                        obj.coords[1],
                        color = 'r' ,
                        marker='s',
                        s=50)
        i = 0
        #Energy and force calculations
        obj.energy = obj.surf.func_eval(obj.coords)
        obj.force = obj.surf.func_prime_eval(obj.coords)
        obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
        ut.log(__name__ , 'STEP: '+str(i)+' Energy: ' +str(round(obj.energy,3))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),3)),2)
        
        obj.coords = obj.coords + obj.stepSize*obj.normF
        
        
        
        while np.max(np.abs(obj.force)) > obj.forceTol and i < obj.maxIter:
            
            i+=1
            
            obj.energy = obj.surf.func_eval(obj.coords)
            obj.force = obj.surf.func_prime_eval(obj.coords)
            
            normF_old = copy.deepcopy(obj.normF)
            obj.normF = obj.surf.norm_func_prime_eval(obj.coords)
            
            ut.log(__name__ , 'STEP: '+str(i)+' Energy: ' +str(round(obj.energy,3))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),5)),2)
            
            if np.dot(obj.normF,normF_old) > 0:
                obj.stepSize *= 1.2
            else:
                obj.stepSize *= 0.4
            
            obj.coords += obj.stepSize*obj.normF
            
            if obj.surf.checkBounds(obj.coords):
                print('outside of range')
                return 0
                
            obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' ,alpha=0.2)
            

            
        obj.minIterations = i
        
        obj.axis.scatter(obj.coords[0],obj.coords[1], color = 'r' , marker='*',s=50)
        
        ut.log(__name__ , 'Minimization Complete, Energy: ' +str(round(obj.energy,5))+ ', Max. Force: '+str(round(np.max(np.abs(obj.force)),5)),1)
        
        return 0

def main():
    
    lattice = lt.lattice()
    lattice.pullConfig()
    
    ut.log(__name__, 'Initializing Minimizer')
    min = Steepest_Descent_adaptive_step()
    lattice.axis = lattice.surf.surfPlot()
    
    ut.log(__name__, 'Running '+min.name+' minimization',1)
    _ = min.run(lattice)
    ut.log(__name__ , 'Fin.',0)
    plt.show()
    
if __name__=='__main__':
    pass

