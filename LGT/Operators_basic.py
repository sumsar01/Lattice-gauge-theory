from Phase import *
from qutip import *
from Fields import *
from Fields_spin12 import *
from Momentum_spin import *
from Potential import *
import numpy as np 

# =============================================================================
# Number operator
# =============================================================================

def num_ope(x, y):
    
    num = 0
    
    j = 0
    while j < y:
        i = 0
        while i < x:
            
            num += ferm_1(x, y, i, j).dag()*ferm_1(x, y, i, j) + ferm_2(x, y, i, j)*ferm_2(x, y, i, j).dag()
            
            i += 1
        j += 1
    
    num = tensor(num, void(x, y))

    return num

def num_ope_site(x, y, i, j):

    num = ferm_1(x, y, i, j).dag()*ferm_1(x, y, i, j) + ferm_2(x, y, i, j)*ferm_2(x, y, i, j).dag()

    num = tensor(num, void(x, y))

    return num



# Loops

def Wilson_loop(x, y, i, j):
    N = x*y
    ind = tensor([identity(2)] * N * 2)
    
    plaquette = Ux(x, y, i, j)*Uy(x, y, i+1, j)*Ux(x, y, i, j+1).dag()*Uy(x, y, i, j).dag()
    loop = plaquette + plaquette.dag()
    loop = tensor(ind, loop)
    
    return loop

def Wilson_loop2(x, y, i, j):
    N = x*y
    ind = tensor([identity(2)] * N * 2)
    
    plaquette = Ux(x, y, i, j)*Ux(x, y, i+1, j)
    loop = plaquette + plaquette.dag()
    loop = tensor(ind, loop)
    
    return loop

def Wilson_loop3(x, y, i, j):
    N = x*y
    ind = tensor([identity(2)] * N * 2)
    
    plaquette = Uy(x, y, i, j)*Uy(x, y, i, j+1)
    loop = plaquette + plaquette.dag()
    loop = tensor(ind, loop)
    
    return loop


























