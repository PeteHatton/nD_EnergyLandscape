
"""
make a class structure for defining each function and its deriv and eval...
"""
class Surface:
    def __init__(self):
        pass

class Styblinski_Tang(Surface):
    def func_eval(self,coords):
        x1 = coords[0]
        x2 = coords[1]
        return 0.5 * (( x1**4 - 16*x1**2 + 5*x1 ) + ( x2**4 - 16*x2**2 + 5*x2 ) )

#print(Styblinski_Tang().func_eval([0,1]))
if __name__ == '__main__':
    pass
    


    
