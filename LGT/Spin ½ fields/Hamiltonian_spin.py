from Phase import *
from qutip import *
from Fields import *
from Fields_spin12 import *
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
# mass term
# =============================================================================
       
    H_m = mass_spin(x, y, r, m, a, C)
#    print('\r 1/5 done, mass term finished' ,end='\n')
        
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
# x momentum
# =============================================================================
    
    H_Tx = momentum_x(C, x, y, a, r)
#    print('\r 4/5 done, momentum x term finished' ,end='\n')    
    
# =============================================================================
# y momentum
# =============================================================================
    
    H_Ty = momentum_y(C, x, y, a, r)
#    print('\r 5/5 done, momentum y term finished' ,end='\n')    

# =============================================================================
# Build Hamiltonian
# =============================================================================
    
    H = a**2*(H_m + H_Tx + H_Ty + H_E + H_B)
    
    return H

