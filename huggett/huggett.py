"""A module for calculating value functions and simulating
   asset demand in an endowment economy with hyperbolic
   discounters"""

import numpy as np
from numpy import exp, log
from math import inf
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

def get_policy_functions(rd, rl, params):
    beta = params['beta']
    betahat = params['betahat']
    gamma = params['gamma']
    delta = params['delta']
    e = params['e']
    Pi = params['Pi']
    aBar = params['aBar']

    # Natural borrowing constraint
    if rl > 0: aBar = max(aBar * rl, -np.min(e) + 1e-4) / rl

    es = range(0, len(e))

    # Discretize grid for resources
    aMax = 10
    numA = 30
    aGrid = create_grid(aMax, numA, aBar)
    idx_l = aGrid < 0
    idx_d = aGrid >= 0
    # if rl > rd:
    #     aGrid = np.hstack([aGrid[idx_l],
    #                        np.array([-inf]),
    #                        np.array([inf]),
    #                        aGrid[idx_d]
    #                       ])
    #     idx_l = aGrid < 0
    #     idx_d = aGrid >= 0
    #     aGrid[aGrid == -inf] = 0
    #     aGrid[aGrid == inf] = 0

    # Initialization
    mNext = np.vstack(((1+rl) * np.expand_dims(aGrid[idx_l], 1) + np.expand_dims(e, 0) - aBar,
                       (1+rd) * np.expand_dims(aGrid[idx_d], 1) + np.expand_dims(e, 0) - aBar
                      ))
    mGrid = mNext
    c = mGrid

    itr = 0
    diffC = 1

    # Solve model using EGM
    while diffC > 1e-8 and itr < 2000:

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
        cNewD = euler(muNext, dCdM, rd, beta, betahat, delta, Pi) ** (-1/gamma)
        cNewL = euler(muNext, dCdM, rl, beta, betahat, delta, Pi) ** (-1/gamma)
        cNew = np.zeros_like(cNewD)
        cNew[idx_d, :] = cNewD[idx_d, :]
        cNew[idx_l, :] = cNewL[idx_l, :]
        mGridNew = np.expand_dims(aGrid, 1) + cNew - aBar

        # Compare policy functions
        cNewComp = np.zeros(np.shape(c))
        for i in es:
            spl = InterpolatedUnivariateSpline(mGridNew[:,i], cNew[:,i])
            cNewComp[:,i] = spl(mGrid[:,i])
            cNewComp[mGrid[:,i] < mGridNew[0,i], i] = mGrid[mGrid[:,i] < mGridNew[0,i], i]

        diffC = np.max(np.abs((cNewComp-c)/c))
        itr += 1
        if itr == 2000: print('Maximum itration count reached!')

        c = cNew
        mGrid = mGridNew

    # sGrid = (mGrid - np.expand_dims(e, 0) + aBar) / (1+rl)
    sGrid = c - np.expand_dims(e, 0) + np.expand_dims(aGrid, 1)
    grids = {'a': aGrid,
             'm': mGrid,
             's': sGrid,
             'c': c,
             'mN': mNext
            }

    return grids

def excess_asset_demand(params, grids, shocks):
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
    r = fsolve(lambda r: get_excess_asset_demand(r, r, params, shocks)  / np.shape(SHOCKS)[0], 0)
    return r

def get_eq_rl(params, rd, shocks):
    #r = broyden1(lambda r: get_excess_asset_demand(r, params, shocks), 0)
    r = fsolve(lambda rl: get_excess_asset_demand(rd, rl, params, shocks)  / np.shape(SHOCKS)[0], rd+0.01)
    return r

def get_excess_asset_demand(rd, rl, params, shocks):
    if rd > rl: return inf 
    return excess_asset_demand(params, get_policy_functions(rd, rl, params), shocks)

def paramfun(params):
    try:
        return get_eq_r(params, SHOCKS)
    except:
        return np.array([np.nan])

def paramfun2(params):
    try:
        return get_eq_rl(params, 0.01, SHOCKS)
    except:
        return np.array([np.nan])

def create_grid(aMax, numA, aBar):
    """Creates initial grid for a"""
    grid = exp(exp(exp(np.linspace(0, log(log(log(aMax+1)+1)+1), numA))-1)-1)-1 + aBar
    if 0 not in grid:
        grid = np.sort(np.append(grid, 0))
    return grid

def euler(muNext, dCdM, r, beta, betahat, delta, Pi):
    """Hyperbolic Euler-equation"""
    return ((1+r) * beta * delta * (dCdM + (1-dCdM)/betahat) * muNext) @ Pi

np.random.seed(19911110)
SHOCKS = np.random.rand(5000, 1000)