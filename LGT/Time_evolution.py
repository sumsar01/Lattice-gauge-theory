from Phase import *
from qutip import *
from Fields import *
from Momentum import *
from Potential import *
from States import *
import numpy as np 


# =============================================================================
# Function for time evolution
# =============================================================================

def time_evol(H, psi, times, coll_ops, e_ops, args):

    #    Options for solver
    opts = Options(atol=1e-10, 
                   rtol=1e-8, 
                   order=14, 
                   nsteps=10000, 
                   store_states=True)
    
#    Return result object from mesolve
    return mesolve(H, psi, times, coll_ops, e_ops, args, options=opts)


#   Build H for time evolution

def build_H(H):
    
    H0 = H[0] + \
    H[1][0][0] +\
    H[1][1][0] + \
    H[2][0][0] + H[2][0][1] + H[2][0][2] + \
    H[2][1][0] + H[2][1][1] + H[2][1][2]
    
    
    H101 = H[1][0][1]
    H102 = H[1][0][2]
    H111 = H[1][1][1]
    H112 = H[1][1][2]
    
    H = [H0, [H101, U], [H102, U2], [H111, U_dag], [H112, U2_dag]]

    return H

def build_H_uniformU(H):

    
    H0 = H[0] + \
    H[1][0][0] + \
    H[1][1][0] + \
    H[2][0][0] + \
    H[2][1][0]
    
    H101 = H[1][0][1]
    H102 = H[1][0][2]
    H111 = H[1][1][1]
    H112 = H[1][1][2]
    
    H201 = H[2][0][1]
    H202 = H[2][0][2]
    H211 = H[2][1][1]
    H212 = H[2][1][2]
    
    
    
    H = [H0, [H101, U], [H102, U2], [H111, U_dag], [H112, U2_dag], [H201, U], [H202, U2], [H211, U_dag], [H212, U2_dag]]

    return H

def build_H_spin(H):

    H = H[0] + \
    H[1][0][0] + H[1][0][1] + H[1][0][2] + \
    H[1][1][0] + H[1][1][1] + H[1][1][2] + \
    H[2][0][0] + H[2][0][1] + H[2][0][2] + \
    H[2][1][0] + H[2][1][1] + H[2][1][2]

    return H











