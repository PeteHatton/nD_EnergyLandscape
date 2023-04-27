import numpy as np
import matplotlib.pyplot as plt
import copy
import sys
import math

import Landscapes as ls
import Utilities as ut
import Lattice as lt
        
class Dimer:
    
    def __init__(self):
    
        self.minEnergy = None
        self.energy = 0
        self.dimerOffset = 0.1
        self.dimerStepSize = 0.005
        self.direction = None # length 1 direction, initially random.
        self.dimerCoords_1 = None
        self.dimerCoords_2 = None
        self.dimerForce_1 = None
        self.dimerForce_2 = None
        self.transForce = [np.inf,np.inf]
        self.dimerEnergy = None
        self.dTheta = None
        self.dimerForce_perp = None
        self.Theta = None
        self.dimerForce_par = None
        self.F_scal = 0
        self.flag = 0
        self.dimerForce_perp_1=None
        self.dimerForce_perp_2=None
        self.maxIter=2000
        
    def run(self,obj):

        #setup and plot
        ut.log(__name__ , 'Initializing Dimer',1)
        self.initializeDimer(obj)
        
        obj.axis.plot([obj.coords[0],self.dimerCoords_1[0]]
                    ,[obj.coords[1],self.dimerCoords_1[1] ]
                    ,color='orange'
                    ,alpha=0.5)
        obj.axis.plot([obj.coords[0],self.dimerCoords_2[0]]
                    ,[obj.coords[1],self.dimerCoords_2[1] ]
                    ,color='orange'
                    ,alpha=0.5)
        iter = 0
        obj.axis.scatter(obj.coords[0],obj.coords[1] ,alpha=1,marker='s',s=30,color='r')
        
        # calculate and save inital energy
        obj.energy = obj.surf.func_eval(obj.coords)
        self.minEnergy = copy.deepcopy(obj.energy)
        
        ut.log(__name__ , 'Walking the dimer',1)
        
        while np.max(np.abs(obj.force))>0.1 and self.flag == 0 and iter<self.maxIter:

            self.calcForce(obj)
            self.calcForceRot(obj)
            self.minDimer_rot(obj)
            
            self.calcForceTrans(obj)
            
            obj.axis.plot([self.dimerCoords_1[0], self.dimerCoords_2[0]]
                        ,[self.dimerCoords_1[1], self.dimerCoords_2[1]]
                        ,color='orange'
                        ,alpha=0.1)

            self.minDimer_translate(obj)
            obj.energy = obj.surf.func_eval(obj.coords)
            obj.axis.scatter(obj.coords[0],obj.coords[1] ,alpha=0.1,color='r')

            ut.log(__name__ , 'Energy: '
                            + str(round(obj.energy,5))
                            + ', Rel. En.: '
                            + str(round(obj.energy - self.minEnergy,5))
                            ,2)
            iter += 1
        ut.log(__name__ , 'Dimer walker finished',1)
        print(iter)
        
        
    def initializeDimer(self,obj):
        
        #initialize with random direction
        ut.log(__name__ , 'Choosing random inital direction',1)
        self.dirRand()
        
        ut.log(__name__ , 'Setting up dimer structure',1)
        
        # step first before start the dimer to push it out of minima first
        #obj.coords += self.dimerOffset * self.direction
        
        # create dimer
        self.dimerCoords_1 = obj.coords + self.dimerOffset * self.direction
        self.dimerCoords_2 = obj.coords - self.dimerOffset * self.direction

    def dirRand(self):

        vecRand = np.random.uniform(low = -1.0, high = 1.0, size = (1,2))[0]
#        vecRand = np.asarray([-0.06669997 , 0.1235402])
#        vecRand = np.asarray([-0.84252301,0.9203875])
#        vecRand = np.asarray([-1,-2])
        print(vecRand)
        self.direction = vecRand / np.linalg.norm(vecRand)
        
    def minDimer_translate(self,obj):
        
        #step center in the new direction with some step size and reinitialize dimer
        obj.coords += self.dimerStepSize * self.transForce
        self.dimerCoords_1 = obj.coords + self.dimerOffset * self.direction
        self.dimerCoords_2 = obj.coords - self.dimerOffset * self.direction

    def minDimer_rot(self,obj):

        self.direction = ( self.direction*math.cos(self.dTheta)
                        + self.Theta * math.sin(self.dTheta))
        self.direction /= np.linalg.norm(self.direction)
        

        self.dimerCoords_1 = obj.coords + self.dimerOffset * self.direction
        self.dimerCoords_2 = obj.coords - self.dimerOffset * self.direction
        

    def calcForce(self,obj):
        
        #force on state and 2 images:
        self.dimerForce_1 = obj.surf.func_prime_eval(self.dimerCoords_1)
        self.dimerForce_2 = obj.surf.func_prime_eval(self.dimerCoords_2)

        obj.force=(self.dimerForce_1 + self.dimerForce_2)/2
        
        
    def calcForceTrans(self,obj):
    
        #force parallel to dimer direction
        self.dimerPar_1 = np.dot(self.dimerForce_1,self.direction)*self.direction
        self.dimerPar_2 = np.dot(self.dimerForce_2,self.direction)*self.direction
        
        self.dimerForce_par = (self.dimerPar_1 + self.dimerPar_2)/2
                
        C = np.dot(self.dimerForce_2 - self.dimerForce_1,self.direction)/(2*self.dimerOffset)
        
        if C > 0:
            self.transForce = - self.dimerForce_par
        else:
            self.transForce = obj.force - 2 * self.dimerForce_par

        self.transForce /= np.linalg.norm(self.transForce)
    
    def calcForceRot(self,obj):
        #calc force on dimer perpendicular to the direction of dimer.
        self.dimerForce_perp_1 = self.dimerForce_1 - (np.dot(self.dimerForce_1,self.direction)*self.direction)
        self.dimerForce_perp_2 = self.dimerForce_2 - (np.dot(self.dimerForce_2,self.direction)*self.direction)
        
        self.dimerForce_perp = self.dimerForce_perp_1 - self.dimerForce_perp_2
        
        self.Theta = self.dimerForce_perp / np.linalg.norm(self.dimerForce_perp)
        
        self.F_scal = np.dot(self.dimerForce_perp ,self.Theta) / self.dimerOffset
        
    #######################################################
    #THIS IS FOR THE DERIV IN EQN 5 IN DIMER PAPER
        #Rot for deriv
        dTheta_deriv = 0.001
        dimerCoords_1_deriv = self.dimerCoords_1 + ( self.direction*np.cos(dTheta_deriv) + self.Theta * np.sin(dTheta_deriv)) * self.dimerOffset
        
        direction_star = ( - dimerCoords_1_deriv + obj.coords ) / self.dimerOffset
        direction_star /= np.linalg.norm(direction_star)
        
        dimerCoords_2_deriv = obj.coords - self.dimerOffset * direction_star
        
        F_dimer_1_star = obj.surf.func_prime_eval(dimerCoords_1_deriv)
        F_dimer_2_star = obj.surf.func_prime_eval(dimerCoords_2_deriv)
        F_dimer_star = (F_dimer_1_star - F_dimer_2_star)
        
        
        dimerForce_perp_1_star = F_dimer_1_star - (np.dot(F_dimer_1_star,direction_star)*direction_star)
        
        dimerForce_perp_2_star = F_dimer_2_star - (np.dot(F_dimer_2_star,direction_star)*direction_star)
        
        dimerForce_perp_star = dimerForce_perp_1_star - dimerForce_perp_2_star
        Theta_star = dimerForce_perp_star / np.linalg.norm(dimerForce_perp_star)
        
        F_prime_scaler = np.abs( (np.dot(F_dimer_star,Theta_star) - np.dot(obj.force,self.Theta)) / (2*dTheta_deriv) )
        
        self.dTheta = 0.5 * math.atan(2*self.F_scal/F_prime_scaler)
        
    #######################################################
        
def main():
    ut.log(__name__ , 'Begin',0)
    
    lattice = lt.lattice()
    
    lattice.pullConfig()
    
    lattice.axis = ls.surfPlot(lattice)
    dimer = Dimer()
    
    dimer.run(lattice)

    # show the generated plot
    ut.log(__name__ , 'Plotting...',0)
    plt.show()
    
    ut.log(__name__ , 'Fin.',0)
    
if __name__ == '__main__':
    pass

