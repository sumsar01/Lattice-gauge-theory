from Phase import *
from Hamiltonian import *
from qutip import *
from Fields import *
from Momentum import *
from Potential import *
from States import *
from Time_evolution import *
from Current_x import *
import numpy as np

# https://arxiv.org/pdf/1606.00342.pdf
# =============================================================================
# Building current
# =============================================================================


# =============================================================================
# Vector current
# =============================================================================

def j_vector_x(x, y, i, j, r, t, args, state):
    
    C = [1.5, -0.3, -1/30]
    
    J_v = 0
    #k = 0
    
    j_0 = 1.j*C[0]/2*(ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j) - 
               ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j))
    
    j_1 = 1.j*C[1]/2*((ferm_bar_1(x, y, i, j)*x_phase(x, y, i)*ferm_2(x, y, i+1, j) - 
               ferm_bar_2(x, y, i, j)*x_phase(x, y, i)*ferm_1(x, y, i+1, j)) - 
    r*(ferm_bar_1(x, y, i, j)*x_phase(x, y, i)*ferm_1(x, y, i+1, j) - 
       ferm_bar_2(x, y, i, j)*x_phase(x, y, i)*ferm_2(x, y, i+1, j)))*U(t, args)
    
    j_2 = 1.j*C[2]/2*((ferm_bar_1(x, y, i, j)*x_phase2(x, y, i)*ferm_2(x, y, i+2, j) - 
               ferm_bar_2(x, y, i, j)*x_phase2(x, y, i)*ferm_1(x, y, i+2, j)) - 
    2*r*(ferm_bar_1(x, y, i, j)*x_phase2(x, y, i)*ferm_1(x, y, i+2, j) - 
       ferm_bar_2(x, y, i, j)*x_phase2(x, y, i)*ferm_2(x, y, i+2, j)))*U2(t, args)
    
    J_v += j_0 + j_1 + j_2
    #k = 1

    j_0 = 1.j*C[0]/2*(ferm_bar_1(x, y, i-1, j)*ferm_2(x, y, i-1, j) - 
               ferm_bar_2(x, y, i-1, j)*ferm_1(x, y, i-1, j))
    
    j_1 = 1.j*C[1]/2*((ferm_bar_1(x, y, i-1, j)*x_phase(x, y, i-1)*ferm_2(x, y, i, j) - 
               ferm_bar_2(x, y, i-1, j)*x_phase(x, y, i-1)*ferm_1(x, y, i, j)) - 
    r*(ferm_bar_1(x, y, i-1, j)*x_phase(x, y, i-1)*ferm_1(x, y, i, j) - 
       ferm_bar_2(x, y, i-1, j)*x_phase(x, y, i-1)*ferm_2(x, y, i, j)))*U(t, args)
    
    j_2 = 1.j*C[2]/2*((ferm_bar_1(x, y, i-1, j)*x_phase2(x, y, i-1)*ferm_2(x, y, i+1, j) - 
               ferm_bar_2(x, y, i-1, j)*x_phase2(x, y, i-1)*ferm_1(x, y, i+1, j)) - 
    2*r*(ferm_bar_1(x, y, i-1, j)*x_phase2(x, y, i-1)*ferm_1(x, y, i+1, j) - 
       ferm_bar_2(x, y, i-1, j)*x_phase2(x, y, i-1)*ferm_2(x, y, i+1, j)))*U2(t, args)

    J_v += j_0 + j_1 + j_2

    # creating H.c
    
    #k = 0
    
    j_0_dag = -1.j*C[0]/2*(ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j) - 
                    ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j))
    
    j_1_dag = -1.j*C[1]/2*((ferm_bar_1(x, y, i+1, j)*x_phase(x, y, i)*ferm_2(x, y, i, j) - 
               ferm_bar_2(x, y, i+1, j)*x_phase(x, y, i)*ferm_1(x, y, i, j)) - 
    r*(ferm_bar_1(x, y, i+1, j)*x_phase(x, y, i)*ferm_1(x, y, i, j) - 
       ferm_bar_2(x, y, i+1, j)*x_phase(x, y, i)*ferm_2(x, y, i, j)))*U_dag(t, args)
    
    j_2_dag = -1.j*C[2]/2*((ferm_bar_1(x, y, i+2, j)*x_phase2(x, y, i)*ferm_2(x, y, i, j) - 
               ferm_bar_2(x, y, i+2, j)*x_phase2(x, y, i)*ferm_1(x, y, i, j)) - 
    2*r*(ferm_bar_1(x, y, i+2, j)*x_phase2(x, y, i)*ferm_1(x, y, i, j) - 
       ferm_bar_2(x, y, i+2, j)*x_phase2(x, y, i)*ferm_2(x, y, i, j)))*U2_dag(t, args)
    
    J_v += j_0_dag + j_1_dag + j_2_dag
    
    #k = 1

    j_0_dag = -1.j*C[0]/2*(ferm_bar_1(x, y, i-1, j)*ferm_2(x, y, i-1, j) - 
               ferm_bar_2(x, y, i-1, j)*ferm_1(x, y, i-1, j))
    
    j_1_dag = -1.j*C[1]/2*((ferm_bar_1(x, y, i, j)*x_phase(x, y, i-1)*ferm_2(x, y, i-1, j) - 
               ferm_bar_2(x, y, i, j)*x_phase(x, y, i-1)*ferm_1(x, y, i-1, j)) - 
    r*(ferm_bar_1(x, y, i, j)*x_phase(x, y, i-1)*ferm_1(x, y, i-1, j) - 
       ferm_bar_2(x, y, i, j)*x_phase(x, y, i-1)*ferm_2(x, y, i-1, j)))*U_dag(t, args)
    
    j_2_dag = -1.j*C[2]/2*((ferm_bar_1(x, y, i+1, j)*x_phase2(x, y, i-1)*ferm_2(x, y, i-1, j) - 
               ferm_bar_2(x, y, i+1, j)*x_phase2(x, y, i-1)*ferm_1(x, y, i-1, j)) - 
    2*r*(ferm_bar_1(x, y, i+1, j)*x_phase2(x, y, i-1)*ferm_1(x, y, i-1, j) - 
       ferm_bar_2(x, y, i+1, j)*x_phase2(x, y, i-1)*ferm_2(x, y, i-1, j)))*U2_dag(t, args)
    
    J_v += j_0_dag + j_1_dag + j_2_dag
    
    j_vector = expect(J_v, state)
    
    return j_vector


# =============================================================================
# Axial current
# =============================================================================

def j_axial_x(x, y, i, j, r, t, args, state):

    C = [1.5, -0.3, -1/30]
    
    J_a = 0
    #k = 0
    
    j_0 = 1.j*C[0]/2*(ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j) - 
               ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j))
    
    j_1 = 1.j*C[1]/2*(ferm_bar_2(x, y, i, j)*x_phase(x, y, i)*ferm_1(x, y, i+1, j) - 
               ferm_bar_1(x, y, i, j)*x_phase(x, y, i)*ferm_2(x, y, i+1, j))*U(t, args)


    j_2 = 1.j*C[2]/2*(ferm_bar_2(x, y, i, j)*x_phase2(x, y, i)*ferm_1(x, y, i+2, j) - 
               ferm_bar_1(x, y, i, j)*x_phase2(x, y, i)*ferm_2(x, y, i+2, j))*U2(t, args)

    J_a += j_0 + j_1 + j_2
    #k = 1
    
    j_0 = 1.j*C[0]/2*(ferm_bar_2(x, y, i-1, j)*ferm_1(x, y, i-1, j) - 
               ferm_bar_1(x, y, i-1, j)*ferm_2(x, y, i-1, j))
    
    j_1 = 1.j*C[1]/2*(ferm_bar_2(x, y, i-1, j)*x_phase(x, y, i-1)*ferm_1(x, y, i, j) - 
               ferm_bar_1(x, y, i-1, j)*x_phase(x, y, i-1)*ferm_2(x, y, i, j))*U(t, args)


    j_2 = 1.j*C[2]/2*(ferm_bar_2(x, y, i-1, j)*x_phase2(x, y, i-1)*ferm_1(x, y, i+1, j) - 
               ferm_bar_1(x, y, i-1, j)*x_phase2(x, y, i-1)*ferm_2(x, y, i+1, j))*U2(t, args)

    J_a += j_0 + j_1 + j_2

    # creating H.c
    
    #k = 0
    
    j_0_dag = - 1.j*C[0]/2*(ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j) - 
               ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j))
    
    j_1_dag = - 1.j*C[1]/2*(ferm_bar_2(x, y, i+1, j)*x_phase(x, y, i)*ferm_1(x, y, i, j) - 
               ferm_bar_1(x, y, i+1, j)*x_phase(x, y, i)*ferm_2(x, y, i, j))*U_dag(t, args)


    j_2_dag = - 1.j*C[2]/2*(ferm_bar_2(x, y, i+2, j)*x_phase2(x, y, i)*ferm_1(x, y, i, j) - 
               ferm_bar_1(x, y, i+2, j)*x_phase2(x, y, i)*ferm_2(x, y, i, j))*U2_dag(t, args)

    J_a += j_0_dag + j_1_dag + j_2_dag
    #k = 1
    
    j_0_dag = - 1.j*C[0]/2*(ferm_bar_2(x, y, i-1, j)*ferm_1(x, y, i-1, j) - 
               ferm_bar_1(x, y, i-1, j)*ferm_2(x, y, i-1, j))
    
    j_1_dag = - 1.j*C[1]/2*(ferm_bar_2(x, y, i, j)*x_phase(x, y, i-1)*ferm_1(x, y, i-1, j) - 
               ferm_bar_1(x, y, i, j)*x_phase(x, y, i-1)*ferm_2(x, y, i-1, j))*U_dag(t, args)


    j_2_dag = - 1.j*C[2]/2*(ferm_bar_2(x, y, i+1, j)*x_phase2(x, y, i-1)*ferm_1(x, y, i-1, j) - 
               ferm_bar_1(x, y, i+1, j)*x_phase2(x, y, i-1)*ferm_2(x, y, i-1, j))*U2_dag(t, args)

    J_a += j_0_dag + j_1_dag + j_2_dag
    
    j_axial = expect(J_a, state)

    return j_axial

















