import numpy as np
from numpy import exp, log
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import broyden1, fsolve
from concurrent.futures import ProcessPoolExecutor

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

def get_policy_functions(r, params):
    beta = params['beta']
    betahat = params['betahat']
    gamma = params['gamma']
    delta = params['delta']
    e = params['e']
    Pi = params['Pi']
    aBar = params['aBar']

    # Natural borrowing constraint
    if r > 0: aBar = max(aBar * r, -np.min(e) + 1e-4) / r

    es = range(0, len(e))

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

        # Find optimal policy function for current V
        cSplines = []
        cNext = np.zeros(np.shape(c))
        dCdM = np.zeros(np.shape(c))
        for i in es:
            cSplines.append(InterpolatedUnivariateSpline(mGrid[:,i], c[:,i]))
            cNext[:,i] = cSplines[i](mNext[:,i])
            #cSplines[i].set_smoothing_factor(1)
            cNext[mNext[:,i] < mGrid[0,i], i] = mNext[mNext[:,i] < mGrid[0,i], i]
            dCdM[:,i] = cSplines[i](mNext[:,i], 1)
        
        # Calculate new policy function
        muNext = mu(cNext, gamma)
        cNew = euler(muNext, dCdM, r, beta, betahat, delta, Pi) ** (-1/gamma)
        mGridNew = np.expand_dims(aGrid, 1) + cNew - aBar

        # Compare policy functions
        cNewComp = np.zeros(np.shape(c))
        for i in es:
            spl = InterpolatedUnivariateSpline(mGridNew[:,i], cNew[:,i])
            cNewComp[:,i] = spl(mGrid[:,i])
            cNewComp[mGrid[:,i] < mGridNew[0,i], i] = mGrid[mGrid[:,i] < mGridNew[0,i], i]
        
        d = (cNewComp-c)/c
        diffC = np.max(np.abs((cNewComp-c)/c))
        iter += 1
        if iter == 2000: print('Maximum iteration count reached!')

        c = cNew
        mGrid = mGridNew
        
    sGrid = (mGrid - np.expand_dims(e, 0) + aBar) / (1+r)
    grids = {'a': aGrid,
             'm': mGrid,
             's': sGrid,
             'c': c,
             'mN': mNext
             }

    return grids

def excess_asset_demand(r, params, grids, shocks):
    sGrid = grids['s']
    mGrid = grids['m']
    aGrid = grids['a']

    Pi   = params['Pi']
    e    = params['e']
    aBar = params['aBar']
    es   = range(0, len(e))

    N, T = np.shape(shocks)
    simA = np.zeros((N, T))
    simE = np.ones((N, T), dtype=np.int)
    PSD  = np.linalg.matrix_power(Pi, 100)[:, 0]

    sSplines = []
    for i in es:
        sSplines.append(InterpolatedUnivariateSpline(sGrid[:,i], aGrid))

    for t in range(1, T):

        # Endowment
        for i in es:
            PiColCS = np.cumsum(Pi[:, i])
            idx = simE[:, t-1] == i
            simE[idx, t] = np.searchsorted(PiColCS, shocks[idx, t])

        # Enforce LLN
        excess = np.zeros(len(e), dtype=np.int)
        for i in range(0, len(e) - 1):
            targets = np.where(simE[:, t] == i)[0]
            misses  = np.where(simE[:, t] > i)[0]
            excess[i] = np.round(len(targets) - round(PSD[i]*N))
            if excess[i] < 0:
                x = misses[0:2]
                simE[x, t] = i
            if excess[i] > 0:
                simE[targets[0 : excess[i]], t] = i + 1

        # Savings
        aNow = np.zeros(N)
        for i in es:
            idx = simE[:, t] == i
            aNow[idx] = sSplines[i](simA[idx, t-1])
        simA[:, t] = np.maximum(aNow, aBar)

    return np.sum(simA[:, -1])

def get_eq_r(params, shocks):
    #r = broyden1(lambda r: get_excess_asset_demand(r, params, shocks), 0)
    r = fsolve(lambda r: get_excess_asset_demand(r, params, shocks), 0)
    return r 

def get_excess_asset_demand(r, params, shocks):
    return excess_asset_demand(r, params, get_policy_functions(r, params), shocks)

def create_grid(aMax, numA, aBar):
    """Creates initial grid for a"""
    return exp(exp(exp(np.linspace(0, log(log(log(aMax+1)+1)+1), numA))-1)-1)-1 + aBar

def euler(muNext, dCdM, r, beta, betahat, delta, Pi):
    """Hyperbolic Euler-equation"""
    return ((1+r) * beta * delta * (dCdM + (1-dCdM)/betahat) * muNext) @ Pi

PARAMS = [[{'beta': beta,
            'betahat': betahat,
            'delta': 0.99,
            'gamma': 3,
            'Pi': np.array([[0.5, 0.075], [0.5, 0.925]]),
            'aBar': -2,
            'e': np.array([0.1, 1])
            }
           for betahat in np.arange(beta, 1.0001, 0.025)
           ]
          for beta in [0.6, 0.7, 0.8, 0.9, 1]
          ]

SHOCKS = np.random.rand(5000, 1000)
MARKERS = ['bx-', 'go-', 'r^-', 'mD-', 'ks'] 

def paramfun(params): return get_eq_r(params, SHOCKS)

def main():
    for i in range(len(PARAMS)):
        with ProcessPoolExecutor(max_workers=2) as executor:
            Gz = []
            for (p, r) in zip(PARAMS[i], executor.map(paramfun, PARAMS[i])):
                Gz.append((p['betahat'], r))
        executor.shutdown()
        G = list(zip(*Gz))
        plt.plot(G[0], G[1], MARKERS[i])
    plt.legend(['beta = 0.6',
                'beta = 0.7',
                'beta = 0.8',
                'beta = 0.9',
                'beta = 1'
                ])
    plt.show()

if __name__ == '__main__':
    main()
