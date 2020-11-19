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
        j = 0
        i = 0
        P = (j*y + i)*2
        A = P+1
        N_tot = x*y*2
        i1 = (i-1) % x
        j1 = (j-1) % y
        E_x = P + N_tot
        E_y = A + N_tot
        E_x_back = (j*y + i1)*2 + N_tot
        E_y_back = (j1*y + i)*2 + N_tot + 1
        if Rep[P] == 1 and Rep[A] == 1:
            charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
        else:
            charge.extend([1])
        
        
        j = 0
        i = 1
        P = (j*y + i)*2
        A = P+1
        N_tot = x*y*2
        i1 = (i-1) % x
        j1 = (j-1) % y
        E_x = P + N_tot
        E_y = A + N_tot
        E_x_back = (j*y + i1)*2 + N_tot
        E_y_back = (j1*y + i)*2 + N_tot + 1
        charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
        
        
        j = 0
        i = 2
        P = (j*y + i)*2
        A = P+1
        N_tot = x*y*2
        i1 = (i-1) % x
        j1 = (j-1) % y
        E_x = P + N_tot
        E_y = A + N_tot
        E_x_back = (j*y + i1)*2 + N_tot
        E_y_back = (j1*y + i)*2 + N_tot + 1
        charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])        
        
        j = 1
        i = 0
        P = (j*y + i)*2
        A = P+1
        N_tot = x*y*2
        i1 = (i-1) % x
        j1 = (j-1) % y
        E_x = P + N_tot
        E_y = A + N_tot
        E_x_back = (j*y + i1)*2 + N_tot
        E_y_back = (j1*y + i)*2 + N_tot + 1
        charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
        
        j = 1
        i = 1
        P = (j*y + i)*2
        A = P+1
        N_tot = x*y*2
        i1 = (i-1) % x
        j1 = (j-1) % y
        E_x = P + N_tot
        E_y = A + N_tot
        E_x_back = (j*y + i1)*2 + N_tot
        E_y_back = (j1*y + i)*2 + N_tot + 1
        charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
    
    
        j = 1
        i = 2
        P = (j*y + i)*2
        A = P+1
        N_tot = x*y*2
        i1 = (i-1) % x
        j1 = (j-1) % y
        E_x = P + N_tot
        E_y = A + N_tot
        E_x_back = (j*y + i1)*2 + N_tot
        E_y_back = (j1*y + i)*2 + N_tot + 1
        if Rep[P] == 0 and Rep[A] == 0:
            charge.extend([Rep[P] + Rep[A] - 1 + (Rep[E_x_back] + Rep[E_y_back] - Rep[E_x] - Rep[E_y])])
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

