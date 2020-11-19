from qutip import *
from itertools import product
import itertools as itertools
from pathlib import Path
import numpy as np 
import time
from Utility import *
import os
from Storage import *




def string_projector(x, y):
    
    Nm = (x*y)*2
    G_n = [0]*x*y
    good_stateReps = [] 
    
    def Generator_on_site(Rep):
        charge = []
        if Rep[16] == 0 and Rep[19] == 0:
            charge.extend([Rep[0] + Rep[1] - 1 + (Rep[16] + Rep[18] - Rep[12] - Rep[13])])
        else:
            charge.extend([1])

        if Rep[21] == 0:
            charge.extend([Rep[2] + Rep[3] - 1 + (Rep[12] + Rep[21] - Rep[14] - Rep[15])])
        else:
            charge.extend([1])       

        if Rep[23] == 0 and Rep[16] == 0:
            charge.extend([Rep[4] + Rep[5] - 1 + (Rep[14] + Rep[23] - Rep[16] - Rep[17])])
        else:
            charge.extend([1]) 

        if Rep[22] == 0 and Rep[19] == 0:
            charge.extend([Rep[6] + Rep[7] - 1 + (Rep[22] + Rep[13] - Rep[19] - Rep[18])])
        else:
            charge.extend([1]) 

        if Rep[21] == 0:
            charge.extend([Rep[8] + Rep[9] - 1 + (Rep[19] + Rep[15] - Rep[20] - Rep[21])])
        else:
            charge.extend([1]) 
 
        if Rep[22] == 0 and Rep[23] == 0:
            charge.extend([Rep[10] + Rep[11] - 1 + (Rep[20] + Rep[17] - Rep[22] - Rep[23])])
        else:
            charge.extend([1])        

        return charge

    
    i = 0
    starttime = time.time()
    while i < 2**(Nm):
        stateRep_fermions = [int(x) for x in '{0:0{1}b}'.format(i,Nm)]
        j = 0
        while j < 3**(Nm):
            num = int(find_ternary(j))
            stateRep_gauge = [int(x)-1 for x in '{0:0{1}}'.format(num,Nm)]
            
            stateRep = stateRep_fermions + stateRep_gauge
              
            if Generator_on_site(stateRep) == [0]*x*y:
                good_stateReps.append(stateRep)
            
            j += 1
        progress = round(i/2**(Nm)*100,4)
        speed = round(progress/(time.time() - starttime),6)
        print('\rN = ' + str(Nm*2) + ', progress = ' + str(progress) + r'%, Speed = ' + str(speed) + r'%/s',end='')   
        
        i += 1
    
    return good_stateReps 
      
   


def string_op(x, y, op_list):
#    Der modtages en liste, op_list, med to-niveau operatorer, hvis tensorprodukt skal projiceres
#    Først lige et par praktiske detaljer omkring hvor vi gemmer og at vi har G_n = 0 for alle sites
    save_dir = './Data/Lattice_Projection/'
    Nm = (x*y)*4
    G_n = [0]*x*y
    filename = 'string_Spin_1_N=' + str(Nm) + '_x=' + str(x) + '_y=' +str(y)
    
#    Jeg tjekker om jeg tidligere har beregnet repræsentationerne af de gode tilstande
#    og hvis jeg ikke har så beregner jeg dem nu
    if os.path.isfile(save_dir + filename + '.p'):
        good_stateReps = load_data(save_dir, filename)
    else:
        good_stateReps = string_projector(x, y)
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



def Generator_on_site(Rep,x ,y):
    charge = []
    if Rep[0] == 1 and Rep[1] == 1 and Rep[16] == 0 and Rep[18] == 0:
        charge.extend([Rep[0] + Rep[1] - 1 + (Rep[16] + Rep[18] - Rep[12] - Rep[13])])
    else:
        charge.extend([1])

    if Rep[21] == 0:
        charge.extend([Rep[2] + Rep[3] - 1 + (Rep[12] + Rep[21] - Rep[14] - Rep[15])])
    else:
        charge.extend([1])       

    if Rep[23] == 0 and Rep[16] == 0:
        charge.extend([Rep[4] + Rep[5] - 1 + (Rep[14] + Rep[23] - Rep[16] - Rep[17])])
    else:
        charge.extend([1]) 

    if Rep[22] == 0 and Rep[18] == 0:
        charge.extend([Rep[6] + Rep[7] - 1 + (Rep[22] + Rep[13] - Rep[19] - Rep[18])])
    else:
        charge.extend([1]) 

    if Rep[21] == 0:
        charge.extend([Rep[8] + Rep[9] - 1 + (Rep[19] + Rep[15] - Rep[20] - Rep[21])])
    else:
        charge.extend([1]) 
 
    if Rep[10] == 0 and Rep[11] == 0 and Rep[22] == 0 and Rep[23] == 0:
        charge.extend([Rep[10] + Rep[11] - 1 + (Rep[20] + Rep[17] - Rep[22] - Rep[23])])
    else:
        charge.extend([1])        

    return charge