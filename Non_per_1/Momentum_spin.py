from Phase import *
from qutip import *
from Fields_spin1 import *
from Fields import *
import numpy as np 
import os
from Storage import *
from Utility import *
from Projector_advanced import *

# =============================================================================
# Building x momentum
# =============================================================================

"only first order correction"
def momentum_x(C, x, y, a, r):
    
    save_dir = './Data/Operators/Momentum/'
    filename = 'momentum_x_a=' + str(a) + '_x=' + str(x) + '_y=' +str(y)

    if os.path.isfile(save_dir + filename + '.p'):
        H_Tx = load_data(save_dir, filename)
    else:
        P = [] 
        j = 0
        while j < y:
            i = 0
            while i < x:
                save_dir2 = './Data/Operators/Simple/'
                filename2 = 'momentum_site_x_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)

                if os.path.isfile(save_dir2 + filename2 + '.p'):
                    moment = load_data(save_dir2, filename2)
                else:
                    jump_P = ope_prod(ferm_dag_1(x, y, i, j), x_phase(x, y, i, j), ferm_1(x, y, i+1, j), Ux(x, y, i, j))
                    jump_A = ope_prod(ferm_dag_2(x, y, i, j), x_phase(x, y, i, j), ferm_2(x, y, i+1, j), Ux(x, y, i, j))
                    annihilate = ope_prod(ferm_dag_2(x, y, i, j), x_phase(x, y, i, j), ferm_1(x, y, i+1, j), Ux(x, y, i, j))
                    create = ope_prod(ferm_dag_1(x, y, i, j), x_phase(x, y, i, j), ferm_2(x, y, i+1, j), Ux(x, y, i, j))

                    jump_P = project_op(x, y, jump_P)
                    jump_A = project_op(x, y, jump_A)
                    annihilate = project_op(x, y, annihilate)
                    create = project_op(x, y, create)
            
                    moment = r*(jump_P - jump_A) + (annihilate - create)
                    moment = moment + moment.dag() 
                    save_data(moment, save_dir2, filename2)
                
                P.append(moment)
          
                i += 1
            j += 1
        H_Tx = sum(P)
        H_Tx = H_Tx*C[1]/(a*2)
        save_data(H_Tx, save_dir, filename)

    return H_Tx

# =============================================================================
# Building y momentum
# =============================================================================

def momentum_y(C, x, y, a, r):

    save_dir = './Data/Operators/Momentum/'
    filename = 'momentum_y_a=' + str(a) + '_x=' + str(x) + '_y=' +str(y)

        
    if os.path.isfile(save_dir + filename + '.p'):
        H_Ty = load_data(save_dir, filename)
    else:
        P = [] 
        j = 0
        while j < y:
            i = 0
            while i < x:
                save_dir2 = './Data/Operators/Simple/'
                filename2 = 'momentum_site_y_' + '_x=' + str(x) + '_y=' + str(y) + '_i=' + str(i) + '_j=' + str(j)

                if os.path.isfile(save_dir2 + filename2 + '.p'):
                    moment = load_data(save_dir2, filename2)
                else:
                    jump_P = ope_prod(ferm_dag_1(x, y, i, j), y_phase(x, y, i, j), ferm_1(x, y, i, j+1), Uy(x, y, i, j))
                    jump_A = ope_prod(ferm_dag_2(x, y, i, j), y_phase(x, y, i, j), ferm_2(x, y, i, j+1), Uy(x, y, i, j))
                    annihilate = ope_prod(ferm_dag_1(x, y, i, j), y_phase(x, y, i, j), ferm_2(x, y, i, j+1), Uy(x, y, i, j))
                    create = ope_prod(ferm_dag_2(x, y, i, j), y_phase(x, y, i, j), ferm_1(x, y, i, j+1), Uy(x, y, i, j))

                    jump_P = project_op(x, y, jump_P)
                    jump_A = project_op(x, y, jump_A)
                    annihilate = project_op(x, y, annihilate)
                    create = project_op(x, y, create)

                    moment = r*(jump_P - jump_A) - 1.j*(annihilate - create)
                    moment = moment + moment.dag() 
                    save_data(moment, save_dir2, filename2)
                
                P.append(moment)
          
                i += 1
            j += 1
        H_Ty = sum(P)
        H_Ty = H_Ty*C[1]/(a*2)
        save_data(H_Ty, save_dir, filename)

    return H_Ty






















