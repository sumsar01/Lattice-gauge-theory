from Phase import *
from qutip import *
from Fields import *
import numpy as np 


#husk at tjekke rep

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
    
    U = [spin_Jm(1), identity(3)]
    
    return tensor([identity(3)] * (target) * 2 + U +
                  [identity(3)] * (N - target - 1) * 2)

def Uy(x, y, i, j):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    target = j*y + i
    N = x*y
    
    U = [identity(3), spin_Jm(1)]
    
    return tensor([identity(3)] * (target) * 2 + U +
                  [identity(3)] * (N - target - 1) * 2)

def void(x, y):
    
    N = x*y
    
    return tensor([identity(3)] * N * 2)
    
def Ex(x, y, i, j, a, g):
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
    
    N = x*y
    target = (j*y + i) % N
    
    U = [-a*g/2*spin_Jz(1), identity(3)]
    
    return tensor([identity(3)] * (target) * 2 + U +
                  [identity(3)] * (N - target - 1) * 2)
    
def Ey(x, y, i, j, a, g):    
    
    if i >= x:
        i = i % x
    
    if j >= y:
        j = j % y
        
    N = x*y
    target = (j*y + i) % N
    
    U = [identity(3), -a*g/2*spin_Jz(1)]
    
    return tensor([identity(3)] * (target) * 2 + U +
                  [identity(3)] * (N - target - 1) * 2)    
    
    
    
    
    
    
    
    
    
    
    
