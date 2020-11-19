from Phase import *
from qutip import *
from Fields import *
from Fields_spin1 import *
from Momentum_spin import *
from Potential import *
import numpy as np 

    
# =============================================================================
# MISSING: evt. time dependance
# =============================================================================


# =============================================================================
# We  the Hamiltonian for a linear, periodic chain of fermionic matter sites,
# with local U(1) symmetry, i.e. interaction with an EM-like gauge field.
# Both matter field are represented by spin-1/2.
# =============================================================================

def Hamiltonian(x, y, m, r, a, g, C):
        
# =============================================================================
# electric term term
# =============================================================================

    H_E = electric_pot(x, y, a, g)
#    print('\r 2/5 done, electric term finished' ,end='\n')
    
# =============================================================================
# electric term term
# =============================================================================   
    
    H_B = magnetic_pot(x, y, a, g)
#    print('\r 3/5 done, magnetic term finished' ,end='\n')
    

# =============================================================================
# Build Hamiltonian
# =============================================================================
    
    H = a**2*(H_E + H_B)
    
    return H

