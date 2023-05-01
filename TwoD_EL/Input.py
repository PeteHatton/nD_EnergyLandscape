import sys
import os
import numpy as np
import TwoD_EL.Utilities as ut

class InputParams():

    def __init__(self):
        
        self.InputFile = ''
        
        #Main
        self.initialxCoord = 0.0
        self.initialyCoord = 0.0
        self.Surface = ''
        self.Dimension = 0
        
        #Minimization
        self.minAlgorithm = ''
        self.minForceTol = 0.0
        self.minMaxIterations = 0
        self.minStepSize = 0.0
        
        #Dimer
        self.printSteps = 0
        self.dimerStepSize = 0.0
        self.dimerOffset = 0.0
        self.dimerVarStepSize = 0
        self.dimerMaxIter = 0
        self.dTheta = 0.0
        self.noDimers = 0
        self.useDeflation = 0
        self.useLocalDeflationOperator = 0
        self.LocalDeflationOperator_alpha = 0.0
        self.LocalDeflationOperator_r = 0.0
        self.LocalDeflationOperator_power=0
        
        self.dimerForceTol = 0.0
        
        
def createInputObject():
    return InputParams()

def getParams(sectionName="", inputParamFile="config.cfg"):
    
    try:
        f = open(inputParamFile, "r")
    except:
        sys.exit(__name__+": ERROR: could not open file: " + inputParamFile)
    
    inputParamObject = createInputObject()
    
    lookForName = 1
    lookForValue = 0
    
    sectionStartMark = '!>'
    sectionEndMark = '!<'
    
    lineIter = 0
    
    for line in f:
        lineIter = lineIter + 1
        line.rstrip()
        
        if line[0] == '#':
            continue
        
        if ((line[0:len(sectionStartMark)] == sectionStartMark or
                              line[0:len(sectionEndMark)] == sectionEndMark)):
            continue
        
            
        paramValueStr = None
        
        if (lookForName == 1):
            if line[0] != "%":
                print(line)
                sys.exit(__name__+": ERROR: something is wrong with the input file: " + inputParamFile)
            else:
                lineLen = len(line)
                paramName = line[1:lineLen-1]
                
                if hasattr(inputParamObject, paramName):
                    paramType = type(getattr(inputParamObject, paramName)).__name__
                else:
                    sys.exit(__name__+": ERROR: unexpected parameter:" + paramName + ": " + inputParamFile)
                
                lookForName = 0
                lookForValue = 1
                
        elif (lookForValue == 1):
            lineLen = len(line)
            paramValueStr = line[0:lineLen-1]
                            
            paramValue = ut.convertStrToType(paramValueStr, paramType)
            
            setattr(inputParamObject, paramName, paramValue)
                            
            lookForName = 1
            lookForValue = 0
        else:
            sys.exit(__name__+": ERROR: undefined action while reading the input file: " + inputParamFile)
    
    f.close()
    
    inputParamObject.initialCoords =[ inputParamObject.initialxCoord,  inputParamObject.initialyCoord ]
    
    InputParams.InputFile = inputParamFile
    
    return inputParamObject
    
if __name__=='__main__':
    minParam = getParams()
