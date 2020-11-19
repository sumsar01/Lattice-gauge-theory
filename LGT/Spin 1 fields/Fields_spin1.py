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

def vac(x, y):
    
    N = x*y
    
    return ([identity(2)] * N * 2)

def Ux(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [spin_Jm(1), identity(3)]
    field = ([identity(3)] * (target) * 2 + U + [identity(3)] * (N - target - 1) * 2)
    field = vac(x, y) + field

    return field

def Uy(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [identity(3), spin_Jm(1)]
    field = ([identity(3)] * (target) * 2 + U + [identity(3)] * (N - target - 1) * 2)
    field = vac(x, y) + field

    return field

def Ux_dag(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [spin_Jm(1).dag(), identity(3)]
    field = ([identity(3)] * (target) * 2 + U + [identity(3)] * (N - target - 1) * 2)
    field = vac(x, y) + field

    return field

def Uy_dag(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [identity(3), spin_Jm(1).dag()]
    field = ([identity(3)] * (target) * 2 + U + [identity(3)] * (N - target - 1) * 2)
    field = vac(x, y) + field

    return field

def void(x, y):
    
    N = x*y
    
    return ([identity(3)] * N * 2)
    
def Ex(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    U = [spin_Jz(1), identity(3)]
    field = ([identity(3)] * (target) * 2 + U + [identity(3)] * (N - target - 1) * 2)
    field = vac(x, y) + field

    return field
    
def Ey(x, y, i, j):    
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
        
    N = x*y
    target = (j*y + i) % N
    
    U = [identity(3), spin_Jz(1)]
    field = ([identity(3)] * (target) * 2 + U + [identity(3)] * (N - target - 1) * 2)
    field = vac(x, y) + field

    return field
    

    
    
    
    
    
    
    
    
