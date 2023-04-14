import numpy as np
import matplotlib.pyplot as plt
import Landscapes as ls


def initialPlot():
    
    xlist = np.linspace(-5.0, 5.0, 500)
    ylist = np.linspace(-5.0, 5.0, 500)

    X, Y = np.meshgrid(xlist, ylist)
    Z = ls.Styblinski_Tang().func_eval([X,Y])
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(X, Y, Z)
    fig.colorbar(cp)
    
    return ax
    
if __name__=='__main__':
    initialPlot()
    plt.show()

