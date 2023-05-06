import numpy as np
import matplotlib.pyplot as plt
import copy
import sys
import math

import nDEL.Landscapes as ls
import nDEL.Utilities as ut
import nDEL.Lattice as lt
import nDEL.Input as ip

class saddleSearcher:

    def __init__(self,params):
        pass

class Dimer:
    
    """
    To-Do:
        - option to just have arrow between start and stop.
        - better adaptive step size calc

    """

    def __init__(self,params):
        self.params = copy.deepcopy(params)
        
        
    def initializeDimer(self,obj):
    
        #energy info
        self.minEnergy = 0
        self.energy = 0
        self.dimerEnergy = 0
        
        #dimer params
        self.dimerOffset = self.params.dimerOffset
        self.dimerStepSize = self.params.dimerStepSize
        self.maxIter = self.params.dimerMaxIter
        self.dimerTol = self.params.dimerForceTol
        
        #dimer coord init
        self.dimerCoords_1 = None
        self.dimerCoords_2 = None
        
        #curvature init
        self.C = None
        
        #force init
        self.direction = None
        self.dimerForce_1 = None
        self.dimerForce_2 = None
        self.transForce = [ np.inf for _ in range(self.params.Dimension) ]
        self.dimerForce_perp = None
        self.dimerForce_par = None
        self.dimerForce_perp_1=None
        self.dimerForce_perp_2=None
        
        # rot angle init
        self.dTheta = None
        self.Theta = None
        self.dTheta_deriv = self.params.dTheta
        self.F_scal = 0
        
        obj.coords = self.params.initialCoords
        
        ut.log(__name__ , 'Initializing random Direction',1)
        #normalized random direction vector
        self.dirRand()
        
        ut.log(__name__ , 'Setting up Dimer structure',1)
        
        # step first before start the dimer to push it out of minima first
        obj.coords += self.dimerOffset * self.direction
        
        # create dimer
        self.dimerCoords_1 = obj.coords + self.dimerOffset * self.direction
        self.dimerCoords_2 = obj.coords - self.dimerOffset * self.direction

    def run(self,obj):

        #setup and plot
        ut.log(__name__ , 'Initializing Dimer ' + str(obj.dimerCount),1)
        self.initializeDimer(obj)
        
        # #plot the initial dimer
        # obj.axis.plot([self.dimerCoords_1[0],self.dimerCoords_2[0]]
        #             ,[self.dimerCoords_1[1],self.dimerCoords_2[1]]
        #             ,color='orange'
        #             ,alpha=0.5)
        
        
        # calculate and save inital energy
        obj.energy = obj.surf.func_eval(obj.coords)
        self.minEnergy = copy.deepcopy(obj.energy)
        
        if self.params.plotSurface and self.params.Dimension == 2:
            obj.axis.scatter(obj.coords[0],obj.coords[1],obj.energy ,alpha=1,marker='s',s=50,color='r')

        ut.log(__name__ , 'Walking the dimer',1)
        
        #initialize step counter
        iter = 0
        
        while np.max(np.abs(obj.force))>self.dimerTol and iter<self.maxIter: 
            
            #calc forces on dimer images
            self.calcForce(obj)
            
            #calc forces for rotation and rotate dimer
            self.calcForceRot(obj)
            self.minDimer_rot(obj)
            
            #calc forces for translation and translate
            self.calcForceTrans(obj)
            self.minDimer_translate(obj)
            
            #Check if we've left the bounds of the surface.
            status = obj.surf.checkBounds(obj.coords)
            if status:
                if self.params.plotSurface and self.params.Dimension == 2:
                    obj.axis.scatter(obj.coords[0],obj.coords[1],obj.energy ,alpha=1,color='r',s=20)

                ut.log(__name__, 'FAILED! OoB! Steps: '+ str(iter)
                                +  '. E = ' + str(round(obj.energy,5))
                                + ', Rel. En.: ' + str(round(obj.energy - self.minEnergy,5))
                                + '. Final coords: ' + str(obj.coords)
                                ,1)
                return 1
            
            
            self.calcStepSize(obj)
            
            #calc new energy
            obj.energy = obj.surf.func_eval(obj.coords)
            
            #plot new position and dimer.
#            obj.axis.plot([self.dimerCoords_1[0], self.dimerCoords_2[0]]
#                        ,[self.dimerCoords_1[1], self.dimerCoords_2[1]]
#                        ,color='orange'
#                        ,alpha=0.01)
            if self.params.plotSteps and self.params.plotSurface:
                obj.axis.scatter(obj.coords[0],obj.coords[1],obj.energy ,alpha=0.01,color='r',s=10)
            
            #increment step counter
            iter += 1
            
            #step log
            if self.params.printSteps:
                ut.log(__name__ ,'Step: ' + str(iter) + ' E = '
                                + str(round(obj.energy,5))
                                + ', Rel. En.: '
                                + str(round(obj.energy - self.minEnergy,5))
                                ,2)
            
        #final log
        if iter < self.maxIter and status==0:
            ut.log(__name__, 'CONVERGED! '+ str(iter) + ' steps.'
                            + ' E = ' + str(round(obj.energy,5))
                            + ', Rel. En.: ' + str(round(obj.energy - self.minEnergy,5))
                            + '. Saddle point: ' + str(obj.coords)
                            ,1)
            if self.params.plotSurface and self.params.Dimension==2:
                obj.axis.scatter(obj.coords[0],obj.coords[1] ,obj.energy ,alpha=1,color='r',s=100,marker='*')

            return 0
        else:
            ut.log(__name__, 'FAILED! Max Steps: '+ str(iter)
                            +  '. E = ' + str(round(obj.energy,5))
                            + ', Rel. En.: ' + str(round(obj.energy - self.minEnergy,5))
                            + '. Final coords: ' + str(obj.coords)
                            ,1)
            if self.params.plotSurface and self.params.Dimension==2:
                obj.axis.scatter(obj.coords[0],obj.coords[1] ,obj.energy ,alpha=1,color='b',s=100,marker='.')
            return 1

    def dirRand(self):
        
        #random vector
        vecRand = np.random.uniform(low = -1.0, high = 1.0, size = (1,self.params.Dimension))[0]
        
        #normalize vector
        self.direction = vecRand / np.linalg.norm(vecRand)
        
        ut.log(__name__, 'Dimer initialized',2)# with direction: ' + str(self.direction) ,2)
        
        return 0
        
    def minDimer_translate(self,obj):
        
        #step center point
        obj.coords += self.dimerStepSize * self.transForce
        
        #reinitialize dimer structure
        self.dimerCoords_1 = obj.coords + self.dimerOffset * self.direction
        self.dimerCoords_2 = obj.coords - self.dimerOffset * self.direction
        
        return 0

    def minDimer_rot(self,obj):
        
        self.direction = ( self.direction*math.cos(self.dTheta)
                        + self.Theta * math.sin(self.dTheta))
        self.direction /= np.linalg.norm(self.direction)
        

        self.dimerCoords_1 = obj.coords + self.dimerOffset * self.direction
        self.dimerCoords_2 = obj.coords - self.dimerOffset * self.direction
        
        return 0

    def calcForce(self,obj):
        mod_1=1
        mod_2=1
        if self.params.useDeflation and not self.params.useLocalDeflationOperator:
        
            mod_1 = self.globalDeflationOp(obj,self.dimerCoords_1)
            mod_2 = self.globalDeflationOp(obj,self.dimerCoords_2)
            
        elif self.params.useDeflation and self.params.useLocalDeflationOperator:

            mod_1 = self.localDeflationOp(obj,self.dimerCoords_1)
            mod_2 = self.localDeflationOp(obj,self.dimerCoords_2)
            
        #force on 2 images:
        self.dimerForce_1 =  mod_1 * obj.surf.func_prime_eval(self.dimerCoords_1)
        self.dimerForce_2 =  mod_2 * obj.surf.func_prime_eval(self.dimerCoords_2)

        #infer the force on the midpoint
        obj.force=(self.dimerForce_1 + self.dimerForce_2)/2
        
        return 0
        
    def globalDeflationOp(self,obj,coords):
        
        mod=1
        if len(obj.saddlePoints)>0:
            for i in range(len(obj.saddlePoints)):
                mod *= 1/np.linalg.norm(coords - obj.saddlePoints[i])**1
            mod += +1
        
        return mod
        
    def localDeflationOp(self,obj,coords):
    
        alpha = self.params.LocalDeflationOperator_alpha
        r = self.params.LocalDeflationOperator_r
        p = self.params.LocalDeflationOperator_power
    
        mod=1
        if len(obj.saddlePoints)>0:
            for i in range(len(obj.saddlePoints)):
                b = 1
                for j in range(self.params.Dimension):
                    if coords[j] < obj.saddlePoints[i][j] + r and coords[j] > obj.saddlePoints[i][j] - r:
                        b *= np.exp( -alpha / ( r**p - ( coords[j] - obj.saddlePoints[i][j])**p )) / np.exp( -alpha / ( r**p ) )
                    else:
                        b *= 0
                mod *= 1/(1-b)
        
        return mod
        
    def calcForceTrans(self,obj):
    
        #force on images parallel to dimer direction
        self.dimerPar_1 = np.dot(self.dimerForce_1,self.direction)*self.direction
        self.dimerPar_2 = np.dot(self.dimerForce_2,self.direction)*self.direction
        
        #force on midpoint parallel to dimer
        self.dimerForce_par = (self.dimerPar_1 + self.dimerPar_2)/2
                
        #local curvature calc.
        self.C = np.dot(self.dimerForce_2 - self.dimerForce_1,self.direction)/(2*self.dimerOffset)
        
        # choose translation direction
        if self.C > 0:
            self.transForce = - self.dimerForce_par
        else:
            self.transForce = obj.force - 2 * self.dimerForce_par
        
        #normalize translation direction
        self.transForce /= np.linalg.norm(self.transForce)
        
        return 0
    
    def calcForceRot(self,obj):
    
        #calc force on dimer perpendicular to the direction of dimer.
        self.dimerForce_perp_1 = self.dimerForce_1 - (np.dot(self.dimerForce_1,self.direction)*self.direction)
        self.dimerForce_perp_2 = self.dimerForce_2 - (np.dot(self.dimerForce_2,self.direction)*self.direction)
        
        #force perp to direction of the dimer on
        self.dimerForce_perp = self.dimerForce_perp_1 - self.dimerForce_perp_2
        
        #normalized vector in direction perp to dimer direction
        self.Theta = self.dimerForce_perp / np.linalg.norm(self.dimerForce_perp)
        
        #some force scaler for rotation angle
        self.F_scal = np.dot(self.dimerForce_perp ,self.Theta) / self.dimerOffset
        
    #######################################################
    #THIS IS FOR THE DERIV IN EQN 5 IN DIMER PAPER
        #Rot for deriv
        
        dimerCoords_1_deriv = self.dimerCoords_1 + ( self.direction*np.cos(self.dTheta_deriv) + self.Theta * np.sin(self.dTheta_deriv)) * self.dimerOffset
        
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
        
        F_prime_scaler = np.abs( (np.dot(F_dimer_star,Theta_star) - np.dot(obj.force,self.Theta)) / (2*self.dTheta_deriv) )
        
        self.dTheta = 0.5 * math.atan(2*self.F_scal/F_prime_scaler)
        
        return 0
        
    #######################################################
    
    def calcStepSize(self,obj):
        if self.params.dimerVarStepSize == 1:
            """
            Linear scaling according to the curvature
            
            """

            step = self.dTheta_deriv
            max = 0.001
            min = 0.00001
            
            self.dimerStepSize = step * self.C / 2.0
            
            if self.dimerStepSize > max:
                self.dimerStepSize = max
            elif self.dimerStepSize < min:
                self.dimerStepSize = min
                
        elif self.params.dimerVarStepSize == 2:
            """
            exp scaling according to the curvature
            
            """

            step = self.dTheta_deriv
            max = 0.001
            min = 0.00001
            
            self.dimerStepSize = step * np.exp(self.C) / 2.0
            
            if self.dimerStepSize > max:
                self.dimerStepSize = max
            elif self.dimerStepSize < min:
                self.dimerStepSize = min
            
        
        return 0

def dimerCounts(obj,params):

    tot_params = params.noDimers
    success = len(obj.saddlePoints)
    failed = tot_params - success
    uniqueSaddles = []

    for s1,sad_1 in enumerate(obj.saddlePoints):
        flag=0
        if s1 == 0:
            uniqueSaddles.append(sad_1)
        else:
            for s2,sad_2 in enumerate(uniqueSaddles):
                if np.all(np.abs(np.asarray(sad_2) - np.asarray(sad_1)) < np.asarray([ params.uniqueSaddleCutoff for i in range(params.Dimension) ])):
                    flag = 1
                    break
            if flag==0:
                uniqueSaddles.append(sad_1)

    return tot_params,success,failed, len(uniqueSaddles)
        
def main():

    ut.log(__name__ , 'Begin',0)
    
    #read Params
    dimerParams = ip.getParams()

    #set up 'lattice'
    lattice = lt.lattice(dimerParams)
    lattice.initialCoords = copy.deepcopy(lattice.coords)

    #plot if we want to and can do
    if dimerParams.plotSurface and dimerParams.Dimension == 2:
        lattice.axis = lattice.surf.surfPlot(lattice)
    
    #run the dimers
    for _ in range(dimerParams.noDimers):
        lattice.dimerCount += 1
        
        #reset the forces on the dimer
        lattice.force = [ np.inf for _ in range(dimerParams.Dimension) ]
        #set up a new dimer
        dimer = Dimer(dimerParams)
        #run the dimer
        status = dimer.run(lattice)
        if not status:
            #if we converged then add that to the list of saddle points
            lattice.saddlePoints.append(lattice.coords.tolist())

    #count the results of the dimer runs and print them
    tot,success,failed, uniques = dimerCounts(lattice,dimerParams)
    ut.log(__name__ , 'DIMER Results. Total: '+ str(tot) +'. Saddles: '+str(success)+'. Uniques: '+str(uniques)+'. Failed: '+str(failed),0)

    # show the generated plot if we can/want
    if dimerParams.plotSurface and dimerParams.Dimension:
        ut.log(__name__ , 'Plotting...',0)
        plt.show()

    ut.log(__name__ , 'Fin.',0)

if __name__ == '__main__':
    pass

