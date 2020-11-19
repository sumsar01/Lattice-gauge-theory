from Phase import *
from qutip import *
from Fields import *
import numpy as np 

# =============================================================================
# Building x momentum
# =============================================================================

def momentum_x(C, x, y, a, r):
    
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
            
            moment0 += C[0]*a/(2)*(ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j) - ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j))
        
#            moment1 += C[1]*a/(2)*((ferm_bar_1(x, y, i, j)*1*r*x_phase(x, y, i)*ferm_1(x, y, i+1, j) - 
 #                       ferm_bar_1(x, y, i, j)*1*r*x_phase(x, y, i)*ferm_1(x, y, i+1, j)) - 
   #                    (ferm_bar_1(x, y, i, j)*x_phase(x, y, i)*ferm_2(x, y, i+1, j) - 
    #                        ferm_bar_2(x, y, i, j)*x_phase(x, y, i)*ferm_1(x, y, i+1, j)))
            
#            moment2 += C[2]*a/(2)*((ferm_bar_1(x, y, i, j)*2*r*x_phase2(x, y, i)*ferm_1(x, y, i+2, j) -
 #                       ferm_bar_1(x, y, i, j)*2*r*x_phase2(x, y, i)*ferm_1(x, y, i+2, j)) - 
  #                     (ferm_bar_1(x, y, i, j)*x_phase2(x, y, i)*ferm_2(x, y, i+2, j) - 
   #                         ferm_bar_2(x, y, i, j)*x_phase2(x, y, i)*ferm_1(x, y, i+2, j)))
          
            
            moment0_dag += C[0]*a/(2)*(ferm_2(x, y, i, j)*ferm_bar_1(x, y, i, j) - ferm_1(x, y, i, j)*ferm_bar_2(x, y, i, j))
         
#            moment1_dag += C[1]*a/(2)*((ferm_1(x, y, i, j)*1*r*x_phase(x, y, i)*ferm_bar_1(x, y, i+1, j) -
 #                           ferm_1(x, y, i, j)*1*r*x_phase(x, y, i)*ferm_bar_1(x, y, i+1, j)) - 
  #                         (ferm_1(x, y, i, j)*x_phase(x, y, i)*ferm_bar_2(x, y, i+1, j) - 
   #                             ferm_2(x, y, i, j)*x_phase(x, y, i)*ferm_bar_1(x, y, i+1, j)))
            
#            moment2_dag += C[2]*a/(2)*((ferm_1(x, y, i, j)*2*r*x_phase2(x, y, i)*ferm_bar_1(x, y, i+2, j) - 
 #                           ferm_1(x, y, i, j)*2*r*x_phase2(x, y, i)*ferm_bar_1(x, y, i+2, j)) - 
  #                         (ferm_1(x, y, i, j)*x_phase2(x, y, i)*ferm_bar_2(x, y, i+2, j) - 
   #                             ferm_2(x, y, i, j)*x_phase2(x, y, i)*ferm_bar_1(x, y, i+2, j)))
        
            i += 1
        j += 1

    P = [moment0, moment1, moment2]
    P_dag = [moment0_dag, moment1_dag, moment2_dag]

    H_Tx = [P, P_dag]

    return H_Tx

# =============================================================================
# Building y momentum
# =============================================================================

def momentum_y(C, x, y, a, r):

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
            
            moment0 += -1.j*C[0]*a/(2)*(ferm_bar_2(x, y, i, j)*ferm_1(x, y, i, j) + 
                           ferm_bar_1(x, y, i, j)*ferm_2(x, y, i, j))
            
            moment1 += C[1]*a/(2)*((ferm_bar_1(x, y, i, j)*(1*r)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1) - 
                        ferm_bar_2(x, y, i, j)*(1*r)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1)) -
                       1.j*(ferm_bar_2(x, y, i, j)*y_phase(x, y, i, j)*ferm_1(x, y, i, j+1) - 
                            ferm_bar_1(x, y, i, j)*y_phase(x, y, i, j)*ferm_2(x, y, i, j+1)))
            
            moment2 += C[2]*a/(2)*((ferm_bar_1(x, y, i, j)*(2*r)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2) - 
                        ferm_bar_2(x, y, i, j)*(2*r)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2)) -
                       1.j*(ferm_bar_2(x, y, i, j)*y_phase2(x, y, i, j)*ferm_1(x, y, i, j+2) - 
                            ferm_bar_1(x, y, i, j)*y_phase2(x, y, i, j)*ferm_2(x, y, i, j+2)))
            
            
            moment0_dag += 1.j*C[0]*a/(2)*(ferm_2(x, y, i, j)*ferm_bar_1(x, y, i, j) + 
                           ferm_1(x, y, i, j)*ferm_bar_2(x, y, i, j))
            
            moment1_dag += C[1]*a/(2)*((ferm_1(x, y, i, j)*(1*r)*y_phase(x, y, i, j)*ferm_bar_1(x, y, i, j+1) - 
                            ferm_1(x, y, i, j)*(1*r)*y_phase(x, y, i, j)*ferm_bar_1(x, y, i, j+1)) +
                       1.j*(ferm_2(x, y, i, j)*y_phase(x, y, i, j)*ferm_bar_1(x, y, i, j+1) + 
                            ferm_1(x, y, i, j)*y_phase(x, y, i, j)*ferm_bar_2(x, y, i, j+1)))
            
            moment2_dag += C[2]*a/(2)*((ferm_1(x, y, i, j)*(2*r)*y_phase2(x, y, i, j)*ferm_bar_1(x, y, i, j+2) - 
                            ferm_1(x, y, i, j)*(2*r)*y_phase2(x, y, i, j)*ferm_bar_1(x, y, i, j+2)) +
                       1.j*(ferm_2(x, y, i, j)*y_phase2(x, y, i, j)*ferm_bar_1(x, y, i, j+2) + 
                            ferm_1(x, y, i, j)*y_phase2(x, y, i, j)*ferm_bar_2(x, y, i, j+2)))
            
            i += 1
        j += 1
    
    P = [moment0, moment1, moment2]
    P_dag = [moment0_dag, moment1_dag, moment2_dag]

    H_Ty = [P, P_dag]
    
    return H_Ty































