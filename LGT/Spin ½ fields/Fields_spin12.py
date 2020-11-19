from Phase import *
from qutip import *
from Fields import *
import numpy as np 
import os
from Storage import *
from Operator_projection import *

# =============================================================================
# Build gauge fields (Spin-1)
# =============================================================================

def Ux(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [sigmam(), identity(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field

    return field

def Uy(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [identity(2), sigmam()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field

    return field

def Ux_dag(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [sigmam().dag(), identity(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field

    return field

def Uy_dag(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [identity(2), sigmam().dag()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field

    return field

def void(x, y):
    
    N = x*y
    
    return ([identity(2)] * N * 2)
    
def Ex(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    U = [sigmaz(), identity(2)]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field

    return field
    
def Ey(x, y, i, j):    
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
        
    N = x*y
    target = (j*y + i) % N
    
    U = [identity(2), sigmaz()]
    field = ([identity(2)] * (target) * 2 + U + [identity(2)] * (N - target - 1) * 2)
    field = void(x, y) + field

    return field
    
    
    
    
    
    
    
    
    
    
    
