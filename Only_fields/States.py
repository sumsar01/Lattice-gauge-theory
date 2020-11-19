from Phase import *
from qutip import *
from Fields import *
from Potential import *
import numpy as np 

# =============================================================================
# Create Vacuum without gauge DoF
# =============================================================================

def vac(x, y):
    
    N = x*y
    
    state = tensor([basis(2, 0), basis(2, 1)] * N)
    
    return state

# =============================================================================
# Create Vacuum with gauge DoF
# =============================================================================

def vac_spin(x, y):
    
    N = x*y
    
    state = tensor([basis(2, 0), basis(2, 1)] * N + [basis(3,1), basis(3,1)] * N)
    
    return state

