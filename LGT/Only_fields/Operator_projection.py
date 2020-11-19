from Gauss_law_advanced import *
from qutip import *
import numpy as np 
import os
from Storage import *

def project_op(x, y, op_list):
#    Der modtages en liste, op_list, med to-niveau operatorer, hvis tensorprodukt skal projiceres
#    Først lige et par praktiske detaljer omkring hvor vi gemmer og at vi har G_n = 0 for alle sites
    save_dir = './Data/Lattice_Projection/'
    Nm = (x*y)*2
    G_n = [0]*x*y
    filename = 'good_stateReps_Spin_1_N=' + str(Nm) + '_x=' + str(x) + '_y=' +str(y)
    
#    Jeg tjekker om jeg tidligere har beregnet repræsentationerne af de gode tilstande
#    og hvis jeg ikke har så beregner jeg dem nu
    if os.path.isfile(save_dir + filename + '.p'):
        good_stateReps = load_data(save_dir, filename)
    else:
        good_stateReps = advanced_projector(x, y)
        save_data(good_stateReps, save_dir, filename)
    num_good_stateReps = len(good_stateReps)
    
#    Nu laver jeg qutip operatorerne om til np.array's
    op_list = [op.full() for op in op_list]
#    Her laver jeg et tomt array, som i sidste ende skal blive til den projicerede operator
    projected_op = np.empty((num_good_stateReps,num_good_stateReps),dtype=np.complex_)
#    Her beregnes matrixelementerne og de sættes ind i det ellers tomme array
    
    l = 0
    starttime = time.time()
    for i, row_state in enumerate(good_stateReps):        
        for j, col_state in enumerate(good_stateReps):
            projected_op[i][j] = np.prod([op[row_state[k]][col_state[k]] for k, op in enumerate(op_list)])
            
            l += 1
            
            progress = round(l/num_good_stateReps**2,4)
            speed = round(progress/(time.time() - starttime),6)
            print('\rProject ope: N = ' + str(Nm) + ', progress = ' + str(progress) + r'%, Speed = ' + str(speed) + r'%/s',end='')    
    
#    Det nu fyldte array laves om til en qutip operator
    projected_op = Qobj(projected_op,dims=[[num_good_stateReps],[num_good_stateReps]])
    
    return projected_op

















