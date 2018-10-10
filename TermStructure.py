# -*- coding: utf-8 -*-
"""
@author: WenqiAngieWu

@para:

n: short-rate periods
bondMtry: bond maturity
N: derivatives expiration periods

u： short-rate up ratio
d : short-rate down ratio
r00: short-rate at t = 0

face: bond's face value
K: strike price of options

for swaps and swaptions only:
  rf: fixed rate of interest rate swaps
  pcpl: notional principal
  n: last payment time
  rp: if pay fixed rate, rp =1; if receive fixed rate: rp = -1
  start: forward start time of swap (if just a swap with no forward, start = 0)
  N: expiration time of swaption

"""

from __future__ import division
import numpy as np

# =============================================================================
def GenerateShortRateTree(n, r00, u, d):
    """Generate short-rate tree"""  
    shortRateTree = np.zeros((n+1, n+1))
    
    # compute the short-rate tree
    shortRateTree[0,0] = r00
    for i in range(1,n+1):
        shortRateTree[0,i] = shortRateTree[0, i-1]*u
        for j in range(1,n+1):
            shortRateTree[j,i] = shortRateTree[j-1, i-1]*d
    
    return shortRateTree

# =============================================================================
def GenerateBondTree(n, bondMtry, face, r00, u, d):
    """ Generate Bond Tree(zero coupon)"""
    q1 = q2 = 0.5
    shortRateTree = GenerateShortRateTree(n, r00, u, d)
    bondTree = np.zeros((bondMtry+1, bondMtry+1))

    # compute the bond tree (zero coupon)
    for j in range(bondMtry+1):
        bondTree[j, bondMtry] = face
    for i in range(bondMtry-1,-1,-1):
        for j in range(i+1):
            bondTree[j, i] = (q1 * bondTree[j, i+1] + q2 * bondTree[j+1, i+1])\
                            /(1 + shortRateTree[j, i]) 
                
    return bondTree

# =============================================================================
def Bond(n, bondMtry, face, r00, u, d):
    """Bond (zero coupon) pricing"""
    tree = GenerateBondTree(n, bondMtry, face, r00, u, d)
    return tree[0,0]


# =============================================================================

def BondOptionAM(n, bondMtry, face, N, r00, K, u, d, cp):
    """first return: American options on bond pricing
       second return: when is the earliest time to exercise"""
    q1 = q2 = 0.5
    shortRateTree = GenerateShortRateTree(n, r00, u, d)
    bondTree = GenerateBondTree(n, bondMtry, face, r00, u, d)
    optionTree = np.zeros((N+1, N+1))
    
    # compute the option tree
    flag = 0 
    list = []
    for j in range(N+1):
        optionTree[j, N] = max(0, cp * (bondTree[j, N]-K))
    for i in range(N-1,-1,-1):
        for j in range(i+1):
            optionTree[j, i] = max((q1 * optionTree[j, i+1] + q2 * optionTree[j+1, i+1])\
                      /(1 + shortRateTree[j, i]), cp * (bondTree[j, i] - K))                        
            if (optionTree[j, i] - cp * (bondTree[j, i] - K)) < 1e-10: # early exercise
                flag += 1
                list.append(i)
    
    when = N
    if(flag):  when = list[-1]
    
                
    return (optionTree[0,0], when)


# =============================================================================
def BondOptionEU(n, bondMtry, face, N, r00, K, u, d, cp):
    """European options on bond pricing"""
    q1 = q2 = 0.5
    shortRateTree = GenerateShortRateTree(n, r00, u, d)
    bondTree = GenerateBondTree(n, bondMtry, face, r00, u, d)
    optionTree = np.zeros((N+1, N+1))
     
    # compute the option tree
    for j in range(N+1):
        optionTree[j, N] = max(0, cp * (bondTree[j, N]-K))
    for i in range(N-1,-1,-1):
        for j in range(i+1):
            optionTree[j, i] = (q1 * optionTree[j, i+1] + q2 * optionTree[j+1, i+1])\
                                /(1 + shortRateTree[j, i])
                        
                
    return optionTree[0,0]

# =============================================================================
def BondForward(n, bondMtry, face, N, r00, u, d):
    """Bond forwards pricing"""
    down = Bond(n, N, 1, r00, u, d)
    up = Bond(n, bondMtry, face, r00, u, d)
    return up/down
    
# =============================================================================
def BondFutures(n, bondMtry, face, N, r00, u, d):  
    """bond futures pricing"""
    q1 = q2 = 0.5
    bondTree = GenerateBondTree(n, bondMtry, face, r00, u, d)
    futuresTree = np.zeros((N+1, N+1))
    
    # compute the futures tree
    for j in range(N+1):
        futuresTree[j, N] = bondTree[j, N]
    
    for i in range(N-1, -1,-1):  #列
        for j in range(i+1):
            futuresTree[j,i] = q1 * futuresTree[j, i+1] + q2 * futuresTree[j+1, i+1]
    
    return futuresTree[0,0]

# =============================================================================
def GenerateElementary(n, N, r00, u, d):
    """elemantary prices (The Forward Equation)"""
    shortRateTree = GenerateShortRateTree(n, r00, u, d)
    elementaryTree = np.zeros((N+1, N+1))
    
    elementaryTree[0,0] = 1
    for j in range(1, N+1):
        elementaryTree[j, j] = 0.5/(1+shortRateTree[j-1, j-1])*elementaryTree[j-1, j-1]
    for i in range(1, N+1):
        elementaryTree[0, i] = 0.5/(1+shortRateTree[0, i-1])*elementaryTree[0, i-1]
    for i in range(2, N+1):
        for j in range(1, i):
            elementaryTree[j, i] = 0.5/(1+shortRateTree[j, i-1])*elementaryTree[j, i-1]+\
                                    0.5/(1+shortRateTree[j-1, i-1])*elementaryTree[j-1, i-1]

    return elementaryTree
    
# =============================================================================
def ForwardSwap(n, pcpl, r00, rf, u, d, rp, start):
    """ swap / forward-start swap pricing using elementary prices"""

    # swap last payment time: n
    # pay fixed rate: rp =1, receive fixed rate: rp = -1
    # forward start time: start, (if just swap, start = 0)
    eleTree = GenerateElementary(n, n, r00, u, d)
    shortRateTree = GenerateShortRateTree(n, r00, u, d)
    
    sum = 0
    for i in range(start, n):
        for j in range(i+1):
            sum += rp*(shortRateTree[j, i] - rf)/(1 + shortRateTree[j, i])*eleTree[j, i]
    
    return sum*pcpl
    
# =============================================================================
def Swaption(n, N, pcpl, K, r00, rf, u, d, rp, start):
    """ swaption pricing"""

    # N: expiration time of swaption
    # n: swap last payment time
    q1 = q2 = 0.5
    shortRateTree = GenerateShortRateTree(n, r00, u, d)
    swapTree = np.zeros((n, n))
    
    # compute swap tree
    for j in range(n):
        swapTree[j, n-1] = rp*(shortRateTree[j, n-1] - rf)/(1 + shortRateTree[j, n-1])
        
    for i in range(n-2,start-1,-1):
        for j in range(i+1):
            swapTree[j, i] = (q1*swapTree[j, i+1] + q2*swapTree[j+1, i+1] \
                    + rp*(shortRateTree[j, i] - rf))/(1 + shortRateTree[j, i])
    if(start):
        for i in range(start-1,-1,-1):
            for j in range(i+1):
                swapTree[j,i] = (q1*swapTree[j, i+1] + q2*swapTree[j+1, i+1])\
                        /(1 + shortRateTree[j, i])
    
    # compute swaption tree
    swaptionTree = np.zeros((N+1, N+1))
    for j in range(N+1):
        swaptionTree[j, N] = (swapTree[j, N] - K)
        
    for i in range(N-1,-1,-1):
        for j in range(i+1):
            swaptionTree[j, i] = (q1*max(swaptionTree[j, i+1] - K, 0) + \
                        q2*max(swaptionTree[j+1, i+1] - K, 0))/(1 + shortRateTree[j, i])
                        
    return swaptionTree[0, 0]*pcpl
 
# =============================================================================
