from Phase import *
from Hamiltonian import *
from qutip import *
from Fields import *
from Momentum import *
from Potential import *
from States import *
from Time_evolution import *
from Current_x import *
from Current_y import *
import numpy as np

# https://arxiv.org/pdf/1606.00342.pdf
# =============================================================================
# Building current
# =============================================================================


# =============================================================================
# Vector current
# =============================================================================

def j_vector_y(x, y, i, j, r, t, args, state):
    
    C = [1.5, -0.3, -1/30]
    
    J_v = 0
    #k = 0
    
    j_0 = - C[0]/2*(ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j))
    
    j_1 = C[1]/2*( - (ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1) + 
               ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1)) - 
    1.j*r*(ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1) - 
       ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1)))*U(t, args)
    
    j_2 = C[2]/2*( - (ferm_bar_1(x, y, i, j)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2) + 
               ferm_bar_2(x, y, i, j)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2)) - 
    1.j*2*r*(ferm_bar_1(x, y, i, j)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2) - 
       ferm_bar_2(x, y, i, j)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2)))*U2(t, args)
    
    J_v += j_0 + j_1 + j_2
    #k = 1

    j_0 = - C[0]/2*(ferm_bar_1(x, y, i, j-1)*ferm_2(x, y, i, j-1) + 
               ferm_bar_2(x, y, i, j-1)*ferm_1(x, y, i, j-1))
    
    j_1 = C[1]/2*( - (ferm_bar_1(x, y, i, j-1)*y_phase(x, y, i, j-1)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j-1)*y_phase(x, y, i, j-1)*ferm_1(x, y, i, j)) - 
    1.j*r*(ferm_bar_1(x, y, i, j-1)*y_phase(x, y, i, j-1)*ferm_1(x, y, i, j) - 
       ferm_bar_2(x, y, i, j-1)*y_phase(x, y, i, j-1)*ferm_2(x, y, i, j)))*U(t, args)
    
    j_2 = 1.j*C[2]/2*( - (ferm_bar_1(x, y, i, j-1)*y_phase2(x, y, i, j-1)*ferm_2(x, y, i, j+1) + 
               ferm_bar_2(x, y, i, j-1)*y_phase2(x, y, i, j-1)*ferm_1(x, y, i, j+1)) - 
    1.j*2*r*(ferm_bar_1(x, y, i, j-1)*y_phase2(x, y, i, j-1)*ferm_1(x, y, i, j+1) - 
       ferm_bar_2(x, y, i, j-1)*y_phase2(x, y, i, j-1)*ferm_2(x, y, i, j+1)))*U2(t, args)

    J_v += j_0 + j_1 + j_2

    # creating H.c
    
    #k = 0
    
    j_0_dag = - C[0]/2*(ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j) + 
                    ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j))
    
    j_1_dag = C[1]/2*( - (ferm_bar_1(x, y, i, j+1)*y_phase(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j+1)*y_phase(x, y, i, j)*ferm_1(x, y, i, j)) + 
    1.j*r*(ferm_bar_1(x, y, i, j+1)*y_phase(x, y, i, j)*ferm_1(x, y, i, j) - 
       ferm_bar_2(x, y, i, j+1)*y_phase(x, y, i, j)*ferm_2(x, y, i, j)))*U_dag(t, args)
    
    j_2_dag = C[2]/2*( - (ferm_bar_1(x, y, i, j+2)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j+2)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j)) + 
    1.j*2*r*(ferm_bar_1(x, y, i, j+2)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j) - 
       ferm_bar_2(x, y, i, j+2)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j)))*U2_dag(t, args)
    
    J_v += j_0_dag + j_1_dag + j_2_dag
    
    #k = 1

    j_0_dag = - C[0]/2*(ferm_bar_1(x, y, i, j-1)*ferm_2(x, y, i, j-1) + 
               ferm_bar_2(x, y, i, j-1)*ferm_1(x, y, i, j-1))
    
    j_1_dag = C[1]/2*( - (ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j-1)*ferm_2(x, y, i, j-1) - 
               ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j-1)*ferm_1(x, y, i, j-1)) + 
    1.j*r*(ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j-1)*ferm_1(x, y, i, j-1) - 
       ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j-1)*ferm_2(x, y, i, j-1)))*U_dag(t, args)
    
    j_2_dag = C[2]/2*( - (ferm_bar_1(x, y, i, j+1)*y_phase2(x, y, i, j-1)*ferm_2(x, y, i, j-1) - 
               ferm_bar_2(x, y, i, j+1)*y_phase2(x, y, i, j-1)*ferm_1(x, y, i, j-1)) + 
    1.j*2*r*(ferm_bar_1(x, y, i, j+1)*y_phase2(x, y, i, j-1)*ferm_1(x, y, i, j-1) - 
       ferm_bar_2(x, y, i, j+1)*y_phase2(x, y, i, j-1)*ferm_2(x, y, i, j-1)))*U2_dag(t, args)
    
    J_v += j_0_dag + j_1_dag + j_2_dag
    
    j_vector = expect(J_v, state)
    
    return j_vector


# =============================================================================
# Axial current
# =============================================================================

def j_axial_y(x, y, i, j, r, t, args, state):

    C = [1.5, -0.3, -1/30]
    
    J_a = 0
    #k = 0
    
    j_0 = C[0]/2*(ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j))
    
    j_1 = C[1]/2*(ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1) + 
               ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1))*U(t, args)


    j_2 = C[2]/2*(ferm_bar_1(x, y, i, j)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2) + 
               ferm_bar_2(x, y, i, j)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2))*U2(t, args)

    J_a += j_0 + j_1 + j_2
    #k = 1
    
    j_0 = C[0]/2*(ferm_bar_1(x, y, i, j-1)*ferm_2(x, y, i, j-1) + 
               ferm_bar_2(x, y, i, j-1)*ferm_1(x, y, i, j-1))
    
    j_1 = C[1]/2*(ferm_bar_1(x, y, i, j-1)*y_phase(x, y, i, j-1)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j-1)*y_phase(x, y, i, j-1)*ferm_1(x, y, i, j))*U(t, args)


    j_2 = C[2]/2*(ferm_bar_1(x, y, i, j-1)*y_phase2(x, y, i, j-1)*ferm_2(x, y, i, j+1) + 
               ferm_bar_2(x, y, i, j-1)*y_phase2(x, y, i, j-1)*ferm_1(x, y, i, j+1))*U2(t, args)

    J_a += j_0 + j_1 + j_2

    # creating H.c
    
    #k = 0
    
    j_0_dag = C[0]/2*(ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j))
    
    j_1_dag = C[1]/2*(ferm_bar_1(x, y, i, j+1)*y_phase(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j+1)*y_phase(x, y, i, j)*ferm_1(x, y, i, j))*U_dag(t, args)


    j_2_dag = C[2]/2*(ferm_bar_1(x, y, i, j+2)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j) + 
               ferm_bar_2(x, y, i, j+2)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j))*U2_dag(t, args)

    J_a += j_0_dag + j_1_dag + j_2_dag
    #k = 1
    
    j_0_dag = C[0]/2*(ferm_bar_1(x, y, i, j-1)*ferm_2(x, y, i, j-1) + 
               ferm_bar_2(x, y, i, j-1)*ferm_1(x, y, i, j-1))
    
    j_1_dag = C[1]/2*(ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j-1)*ferm_2(x, y, i, j-1) + 
               ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j-1)*ferm_1(x, y, i, j-1))*U_dag(t, args)


    j_2_dag = C[2]/2*(ferm_bar_1(x, y, i, j+1)*y_phase2(x, y, i, j-1)*ferm_2(x, y, i, j-1) + 
               ferm_bar_2(x, y, i, j+1)*y_phase2(x, y, i, j-1)*ferm_1(x, y, i, j-1))*U2_dag(t, args)

    J_a += j_0_dag + j_1_dag + j_2_dag
    
    j_axial = expect(J_a, state)

    return j_axial





















