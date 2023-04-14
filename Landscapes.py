import numpy as np
"""
make a class structure for defining each function and its deriv and eval...
"""
class Surface:
    def __init__(self):
        pass
    
    def func_eval(self):
        pass

class Styblinski_Tang(Surface):
    def func_eval(self,coords):
        x1 = coords[0]
        x2 = coords[1]
        return 0.5 * (( x1**4 - 16*x1**2 + 5*x1 ) + ( x2**4 - 16*x2**2 + 5*x2 ) )
    
    def func_prime_eval(self,coords):
        x1 = coords[0]
        x2 = coords[1]
        return np.asarray([ -0.5 * ( 4*x1**3 - 32*x1 + 5 ),
                 -0.5 * ( 4*x2**3 - 32*x2 + 5 ) ])
    
    def norm_func_prime_eval(self,coords):
        eval = self.func_prime_eval(coords)
        return eval/np.linalg.norm(eval)
        
if __name__ == '__main__':
    pass
    


    
