from Phase import *
from qutip import *
from Fields_spin12 import *
from Fields import *
import numpy as np 
import os
from Storage import *

# =============================================================================
# Building x momentum
# =============================================================================


def Make_momentum_x(C, x, y, a, r):
    
    
    save_dir = './Data/Operators/Momentum/'
    filename = 'momentum_x_a=' + str(a) + '_x=' + str(x) + '_y=' +str(y)
    
    
    if os.path.isfile(save_dir + filename + '.p'):
        H_Tx = load_data(save_dir, filename)
    else:
        H_Tx = momentum_x_spin(C, x, y, a, r)
        save_data(H_Tx, save_dir, filename)

        
    return H_Tx


"only first order correction"
def momentum_x_spin(C, x, y, a, r):
    
    moment0 = 0
    moment1 = 0
    moment2 = 0

    moment0_dag = 0
    moment1_dag = 0
    moment2_dag = 0

    j = 0
    while j < y:
        i = 0
        while i < x:
            
            moment0 += C[0]/(a*2)*(ferm_dag_1(x, y, i, j)*ferm_2(x, y, i, j) - 
                        ferm_dag_2(x, y, i, j)*ferm_1(x, y, i, j))
            
            
            
            moment1 += C[1]/(a*2)*((ferm_dag_1(x, y, i, j)*1*r*x_phase(x, y, i, j)*ferm_1(x, y, i+1, j) - 
                        ferm_dag_2(x, y, i, j)*1*r*x_phase(x, y, i, j)*ferm_2(x, y, i+1, j)) + 
                       (ferm_dag_2(x, y, i, j)*x_phase(x, y, i, j)*ferm_1(x, y, i+1, j) - 
                            ferm_dag_1(x, y, i, j)*x_phase(x, y, i, j)*ferm_2(x, y, i+1, j)))*Ux(x, y, i, j)
            
            moment2 += C[2]/(a*2)*((ferm_1(x, y, i, j).dag()*2*r*x_phase2(x, y, i, j)*ferm_1(x, y, i+2, j) -
                        ferm_2(x, y, i, j).dag()*2*r*x_phase2(x, y, i, j)*ferm_2(x, y, i+2, j)) + 
            (ferm_2(x, y, i, j).dag()*x_phase2(x, y, i, j)*ferm_1(x, y, i+2, j) - 
             ferm_1(x, y, i, j).dag()*x_phase2(x, y, i, j)*ferm_2(x, y, i+2, j)))*Ux(x, y, i, j)*Ux(x, y, i+1, j)
            
            i += 1
        j += 1

    P = moment0 + moment1 + moment2
    P_dag = moment0.dag() + moment1.dag() + moment2.dag()
    H_Tx = P + P_dag

    return H_Tx

# =============================================================================
# Building y momentum
# =============================================================================

def Make_momentum_y(C, x, y, a, r):
    
    save_dir = './Data/Operators/Momentum/'
    filename = 'momentum_y_a=' + str(a) + '_x=' + str(x) + '_y=' +str(y)
    
    if os.path.isfile(save_dir + filename + '.p'):
        H_Ty = load_data(save_dir, filename)
    else:
        H_Ty = momentum_y_spin(C, x, y, a, r)
        save_data(H_Ty, save_dir, filename)
        
    return H_Ty


def momentum_y_spin(C, x, y, a, r):

    moment0 = 0
    moment1 = 0
    moment2 = 0

    moment0_dag = 0
    moment1_dag = 0
    moment2_dag = 0
        
    j = 0
    while j < y:
        i = 0
        while i < x:
            
            moment0 += -1.j*C[0]/(a*2)*(ferm_dag_1(x, y, i, j)*ferm_2(x, y, i, j) + 
                           ferm_dag_2(x, y, i, j)*ferm_1(x, y, i, j))
            
            moment1 += C[1]/(a*2)*((ferm_dag_1(x, y, i, j)*(1*r)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1) - 
                        ferm_dag_2(x, y, i, j)*(1*r)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1)) -
                       1.j*(ferm_dag_1(x, y, i, j)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1) + 
                            ferm_dag_2(x, y, i, j)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1)))*Uy(x, y, i, j)
            
            moment2 += C[2]/(a*2)*((ferm_dag_1(x, y, i, j)*(2*r)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2) - 
                        ferm_dag_2(x, y, i, j)*(2*r)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2)) -
                       1.j*(ferm_dag_1(x, y, i, j)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2) + 
                            ferm_dag_2(x, y, i, j)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2)))*Uy(x, y, i, j)*Uy(x, y, i, j+1)
            
            i += 1
        j += 1
    
    P = moment0 + moment1 + moment2
    P_dag = moment0.dag() + moment1.dag() + moment2.dag()

    H_Ty = P + P_dag
    
    return H_Ty