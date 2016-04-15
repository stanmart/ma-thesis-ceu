import numpy as np
from numpy import exp, log
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

# Utility functions

def u(c, gamma):
    """CRRA utility function"""
    if gamma == 1:
        return log(c)
    else:
        return c ** (1-gamma) / (1-gamma)

def mu(c, gamma):
    """Marginal utility function"""
    return c ** (-gamma)

def excess_asset_demand(r, params, shocks):
    beta = params['beta']
    betahat = params['betahat']
    gamma = params['gamma']
    delta = params['delta']
    e = params['e']
    Pi = params['Pi']
    aBar = params['aBar']

    # Discretize grid for resources
    aMax = 10
    numA = 30
    aGrid = create_grid(aMax, numA, aBar)
    
    # Initialization
    mNext = (1+r) * np.expand_dims(aGrid, 1) + np.expand_dims(e, 0) - aBar
    mGrid = mNext
    c     = mGrid
    
    iter  = 0
    diffC = 1
    
    # Solve model using EGM
    
    while diffC > 1e-8 and iter < 2000:
        # print(iter)
        # Find optimal policy function for current V
        cNext = np.zeros(np.shape(c))
        for i in range(0, len(e)):
            cNext[:,i] = np.interp(mNext[:,i], mGrid[:,i], c[:,i])
            cNext[mNext[:,i]<mGrid[0,i], i] = mNext[mNext[:,i]<mGrid[0,i], i]

        #polyOrd = 3
        #dCdM = np.zeros(np.shape(c))
        #cPol = np.zeros((polyOrd+1, len(e)))
        #for i in range(0, len(e)):
        #    cPol[:,i] = np.polyfit(mNext[:,i], cNext[:,i], polyOrd)
        #    dCdM[:,i] = np.polyval(np.polyder(cPol[:, i]), mNext[:,i])

        dCdM = np.zeros(np.shape(c))
        cSpline = np.zeros(np.shape(c))
        for i in range(0, len(e)):
            spl = UnivariateSpline(mNext[:,i], cNext[:,i], k=3)
            dCdM[:,i] = spl(mNext[:, i], 1)
            cSpline[:,i] = spl(mNext[:, i], 0)

        muNext = mu(cNext, gamma)
        cNew = euler(muNext, r, beta, betahat, delta, dCdM, Pi) ** (-1/gamma)
        mGridNew = np.expand_dims(aGrid, 1) + cNew - aBar
        
        cNewComp = np.zeros(np.shape(c))
        for i in range(0, len(e)):
            cNewComp[:,i] = np.interp(mGrid[:,i], mGridNew[:,i], cNew[:,i])
            cNewComp[mGrid[:,i]<mGridNew[0,i], i] = mGrid[mGrid[:,i]<mGridNew[0,i], i]
        
        d = (cNewComp-c)/c
        diffC = np.max(np.abs((cNewComp-c)/c))
        iter += 1
        if iter == 2000: print('Maximum iteration count reached!')

        c = cNew
        mGrid = mGridNew

    plt.plot(mGrid, c, 'b-', mNext, dCdM, 'r--', mNext, cSpline, 'ko')
    plt.show()
    return np.hstack([mGrid, c])



def create_grid(aMax, numA, aBar):
    """Creates initial grid for a"""
    return exp(exp(exp(np.linspace(0, log(log(log(aMax+1)+1)+1), numA))-1)-1)-1 + aBar

def euler(muNext, r, beta, betahat, delta, dCdM, Pi):
    """Hyperbolic Euler-equation"""
    #dCdM = np.diff(cNext, axis=0) / np.diff(mGrid, axis=0)
    #dCdM = np.vstack([np.ones((1, 2)), dCdM]) # temporary solution!
    return ((1+r) * beta * delta * (dCdM + (1-dCdM)/betahat) * muNext) @ Pi

params = {'beta': 0.7,
          'delta': 0.99,
          'gamma': 3,
          'Pi': np.array([[0.5, 0.075], [0.5, 0.925]]),
          'aBar': -2,
          'e': np.array([0.1, 1]),
          'betahat': 0.7
          }
print(excess_asset_demand(0, params, ""))
